import sys
import datetime
import numpy as np
import xarray as xr
import scipy.constants
import scipy.stats
import xopr.radar_util
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

# Snakemake inputs/outputs
input_surfbed = snakemake.input.surfbed
output_file = snakemake.output.nc if hasattr(snakemake.output, 'nc') else snakemake.output[0]
status_file = snakemake.output.status
log_file = snakemake.log[0]

# Parameters
chunk_path = snakemake.params.chunk_path
ice_permittivity = snakemake.params.ice_permittivity
min_window_size_km = snakemake.params.min_window_size_km
max_window_size_km = snakemake.params.max_window_size_km
window_step_km = snakemake.params.window_step_km
alpha_range = snakemake.params.alpha_range_db_km
alpha_step = snakemake.params.alpha_step_db_km
correlation_threshold = snakemake.params.correlation_threshold

sys.stderr = open(log_file, 'w')


def estimate_attenuation_adaptive(distance_km, ice_thickness_m, P_geom,
                                   alpha_range, alpha_step,
                                   min_window_km, max_window_km, window_step_km,
                                   correlation_threshold):
    """Estimate attenuation rate using adaptive windowing and decorrelation."""

    n_points = len(distance_km)
    attenuation_rate = np.full(n_points, np.nan)
    correlation_length = np.full(n_points, np.nan)
    correlation_coef = np.full(n_points, np.nan)

    # Test attenuation values
    alpha_test = np.arange(alpha_range[0], alpha_range[1] + alpha_step, alpha_step)

    for i in range(n_points):
        center_dist = distance_km[i]

        # Try increasing window sizes
        for window_size in np.arange(min_window_km, max_window_km + window_step_km, window_step_km):
            # Select data within window
            half_window = window_size / 2
            in_window = (distance_km >= center_dist - half_window) & (distance_km <= center_dist + half_window)

            if np.sum(in_window) < 10:  # Need at least 10 points
                continue

            z_window = ice_thickness_m[in_window]
            P_window = P_geom[in_window]

            # Remove NaN values
            valid = ~(np.isnan(z_window) | np.isnan(P_window))
            if np.sum(valid) < 10:
                continue

            z_window = z_window[valid]
            P_window = P_window[valid]

            # Test different attenuation rates
            correlations = []
            for alpha in alpha_test:
                # Apply attenuation correction: P_corrected = P_geom * 10^(2*alpha*z/10)
                # Convert alpha from dB/km to dB/m and apply
                P_corrected = P_window * np.power(10, 2 * alpha * z_window / 1000 / 10)

                # Calculate correlation with ice thickness
                if len(P_corrected) > 2:
                    corr, _ = scipy.stats.pearsonr(z_window, P_corrected)
                    correlations.append(np.abs(corr))
                else:
                    correlations.append(np.inf)

            # Find alpha that minimizes correlation
            correlations = np.array(correlations)
            best_idx = np.argmin(correlations)
            best_corr = correlations[best_idx]

            # Check if decorrelation is sufficient
            if best_corr < correlation_threshold:
                attenuation_rate[i] = alpha_test[best_idx]
                correlation_length[i] = window_size
                correlation_coef[i] = best_corr
                break

        # If no window converged, use result from largest window
        if np.isnan(attenuation_rate[i]) and len(correlations) > 0:
            attenuation_rate[i] = alpha_test[best_idx]
            correlation_length[i] = window_size
            correlation_coef[i] = best_corr
            print(f"Warning: Point at distance {center_dist:.2f} km did not converge within max window size, using largest window result", file=sys.stderr)

    return attenuation_rate, correlation_length, correlation_coef


try:
    print(f"[{datetime.datetime.now()}] Loading surface/bed data from {input_surfbed}", file=sys.stderr)

    # Load surface/bed data
    ds = xr.open_dataset(input_surfbed)

    # Check if processing was successful
    if 'processing_error' in ds.attrs:
        raise ValueError(f"Input data has processing error: {ds.attrs['processing_error']}")

    # Check required variables
    required_vars = ['surface_twtt', 'bed_twtt', 'surface_power_dB', 'bed_power_dB', 'Latitude', 'Longitude']
    missing_vars = [v for v in required_vars if v not in ds]
    if missing_vars:
        raise ValueError(f"Missing required variables: {missing_vars}")

    # Check if we have valid data
    if len(ds.slow_time) < 10:
        raise ValueError(f"Insufficient data points: {len(ds.slow_time)}")

    print(f"[{datetime.datetime.now()}] Calculating ice thickness and geometric correction", file=sys.stderr)

    # Calculate ice velocity and thickness
    c = scipy.constants.c  # Speed of light in m/s
    v_ice = c / np.sqrt(ice_permittivity)  # Speed in ice

    # Ice thickness from two-way travel time
    ice_thickness = (ds['bed_twtt'] - ds['surface_twtt']) * v_ice / 2  # meters

    # Add distance along track using xopr utility
    ds = xopr.radar_util.add_along_track(ds)
    distance_km = ds['along_track'].values / 1000  # Convert meters to km

    # Calculate aircraft altitude above ice surface from surface TWTT
    # Surface TWTT is two-way travel time through air to ice surface
    # One-way distance = surface_twtt * c / 2
    h = ds['surface_twtt'].values * c / 2  # Aircraft height above ice surface (meters)

    print(f"[{datetime.datetime.now()}] Aircraft altitude above surface: mean={np.nanmean(h):.1f}m, range={np.nanmin(h):.1f}-{np.nanmax(h):.1f}m", file=sys.stderr)

    # Geometric spreading correction (Matsuoka et al. 2010, eq. 2.3)
    # G = (h + z/sqrt(epsilon))^2
    # P_geom_corrected = P_bed * G
    epsilon = ice_permittivity
    z = ice_thickness.values
    G = np.power(h + z / np.sqrt(epsilon), 2)

    # Convert bed power from dB to linear
    P_bed_dB = ds['bed_power_dB'].values
    P_bed_linear = np.power(10, P_bed_dB / 10)

    # Apply geometric correction
    P_geom = P_bed_linear * G

    print(f"[{datetime.datetime.now()}] Estimating attenuation rate using adaptive windowing", file=sys.stderr)
    print(f"[{datetime.datetime.now()}]   Window range: {min_window_size_km} - {max_window_size_km} km", file=sys.stderr)
    print(f"[{datetime.datetime.now()}]   Alpha range: {alpha_range[0]} - {alpha_range[1]} dB/km", file=sys.stderr)

    # Estimate attenuation rate
    attenuation_rate, correlation_length_km, correlation_coef = estimate_attenuation_adaptive(
        distance_km, z, P_geom,
        alpha_range, alpha_step,
        min_window_size_km, max_window_size_km, window_step_km,
        correlation_threshold
    )

    # Calculate corrected power using estimated attenuation
    P_atten_corrected = P_geom * np.power(10, 2 * attenuation_rate * z / 1000 / 10)

    # Calculate basal reflectivity (relative to surface)
    P_surface_dB = ds['surface_power_dB'].values
    P_surface_linear = np.power(10, P_surface_dB / 10)

    # Relative basal reflectivity
    basal_reflectivity = P_atten_corrected / P_surface_linear
    basal_reflectivity_dB = 10 * np.log10(basal_reflectivity)

    print(f"[{datetime.datetime.now()}] Creating output dataset", file=sys.stderr)

    # Create output dataset
    result_ds = xr.Dataset({
        'ice_thickness': (['slow_time'], z),
        'aircraft_altitude': (['slow_time'], h),
        'distance_along_track': (['slow_time'], distance_km),
        'geometric_spreading_factor': (['slow_time'], G),
        'bed_power_geom_corrected': (['slow_time'], P_geom),
        'attenuation_rate': (['slow_time'], attenuation_rate),
        'correlation_length': (['slow_time'], correlation_length_km),
        'correlation_coef': (['slow_time'], correlation_coef),
        'basal_reflectivity': (['slow_time'], basal_reflectivity),
        'basal_reflectivity_dB': (['slow_time'], basal_reflectivity_dB),
    }, coords={
        'slow_time': ds.slow_time,
        'Latitude': ds.Latitude,
        'Longitude': ds.Longitude,
    })

    # Copy Elevation if present
    if 'Elevation' in ds.coords:
        result_ds.coords['Elevation'] = ds.Elevation

    # Add metadata
    result_ds.attrs['chunk_path'] = chunk_path
    result_ds.attrs['ice_permittivity'] = ice_permittivity
    result_ds.attrs['aircraft_altitude_agl_method'] = 'calculated from surface_twtt * c / 2'
    result_ds.attrs['mean_aircraft_altitude_agl_m'] = float(np.nanmean(h))
    result_ds.attrs['min_window_size_km'] = min_window_size_km
    result_ds.attrs['max_window_size_km'] = max_window_size_km
    result_ds.attrs['alpha_range_db_km'] = alpha_range
    result_ds.attrs['correlation_threshold'] = correlation_threshold

    # Copy relevant attributes from input
    for attr in ['num_frames', 'frame_ids', 'data_product', 'season', 'segment']:
        if attr in ds.attrs:
            result_ds.attrs[attr] = ds.attrs[attr]

    # Add variable attributes
    result_ds['ice_thickness'].attrs = {'units': 'm', 'long_name': 'Ice thickness'}
    result_ds['aircraft_altitude'].attrs = {'units': 'm', 'long_name': 'Aircraft altitude above ice surface', 'description': 'Calculated from surface_twtt * c / 2'}
    result_ds['distance_along_track'].attrs = {'units': 'km', 'long_name': 'Distance along track'}
    result_ds['attenuation_rate'].attrs = {'units': 'dB/km', 'long_name': 'Englacial attenuation rate'}
    result_ds['correlation_length'].attrs = {'units': 'km', 'long_name': 'Spatial correlation length'}
    result_ds['basal_reflectivity_dB'].attrs = {'units': 'dB', 'long_name': 'Basal reflectivity (relative to surface)'}

    # Save
    result_ds.to_netcdf(output_file)
    print(f"[{datetime.datetime.now()}] Saved basal properties to {output_file}", file=sys.stderr)

    # Statistics
    valid_atten = attenuation_rate[~np.isnan(attenuation_rate)]
    if len(valid_atten) > 0:
        print(f"[{datetime.datetime.now()}] Attenuation rate statistics:", file=sys.stderr)
        print(f"[{datetime.datetime.now()}]   Mean: {np.mean(valid_atten):.2f} dB/km", file=sys.stderr)
        print(f"[{datetime.datetime.now()}]   Median: {np.median(valid_atten):.2f} dB/km", file=sys.stderr)
        print(f"[{datetime.datetime.now()}]   Range: {np.min(valid_atten):.2f} - {np.max(valid_atten):.2f} dB/km", file=sys.stderr)

    # Write success status
    with open(status_file, 'w') as f:
        f.write("success")
    print(f"[{datetime.datetime.now()}] Basal properties estimated successfully", file=sys.stderr)

except Exception as e:
    error_msg = f"{type(e).__name__}: {str(e)}"
    print(f"[{datetime.datetime.now()}] Error estimating basal properties: {error_msg}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)

    # Create minimal output file
    empty_ds = xr.Dataset()
    empty_ds.attrs['chunk_path'] = chunk_path
    empty_ds.attrs['processing_error'] = error_msg
    empty_ds.to_netcdf(output_file)

    # Write failure status
    with open(status_file, 'w') as f:
        f.write(f"failed: {error_msg}")

    print(f"[{datetime.datetime.now()}] Chunk marked as failed, workflow will continue", file=sys.stderr)

sys.stderr.close()
