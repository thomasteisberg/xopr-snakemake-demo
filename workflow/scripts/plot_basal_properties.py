import sys
import datetime
import numpy as np
import xarray as xr
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec

# Snakemake inputs/outputs
input_basal_props = snakemake.input.basal_props
input_surfbed = snakemake.input.surfbed
output_plot = snakemake.output.png
log_file = snakemake.log[0]

# Parameters
chunk_path = snakemake.params.chunk_path

sys.stderr = open(log_file, 'w')

try:
    print(f"[{datetime.datetime.now()}] Loading basal properties from {input_basal_props}", file=sys.stderr)

    # Load basal properties
    props = xr.open_dataset(input_basal_props)

    # Check if processing was successful
    if 'processing_error' in props.attrs:
        raise ValueError(f"Input data has processing error: {props.attrs['processing_error']}")

    # Load original surface/bed data for comparison
    surfbed = xr.open_dataset(input_surfbed)

    print(f"[{datetime.datetime.now()}] Creating multi-panel plot", file=sys.stderr)

    # Create figure with 6 panels on left and map on right
    fig = plt.figure(figsize=(18, 14))
    gs = GridSpec(6, 2, figure=fig, hspace=0.3, wspace=0.3, width_ratios=[2, 1])

    # Get distance along track
    distance = props['distance_along_track'].values

    # Get lat/lon for map (handle both uppercase and lowercase coordinate names)
    lat_coord = 'Latitude' if 'Latitude' in props else 'latitude'
    lon_coord = 'Longitude' if 'Longitude' in props else 'longitude'
    lat = props[lat_coord].values
    lon = props[lon_coord].values

    # Panel 1: Ice Thickness Profile
    ax1 = fig.add_subplot(gs[0, 0])
    ax1.plot(distance, props['ice_thickness'].values, linewidth=1.5, color='steelblue')
    ax1.fill_between(distance, 0, props['ice_thickness'].values, alpha=0.3, color='steelblue')
    ax1.set_ylabel('Ice Thickness (m)', fontsize=11, fontweight='bold')
    ax1.set_title(f'Chunk {chunk_path}: Basal Properties Analysis', fontsize=13, fontweight='bold')
    ax1.grid(True, alpha=0.3)
    ax1.set_xlim(distance[0], distance[-1])

    # Add mean thickness annotation
    mean_thickness = np.nanmean(props['ice_thickness'].values)
    ax1.text(0.02, 0.95, f'Mean: {mean_thickness:.0f} m',
             transform=ax1.transAxes, verticalalignment='top',
             bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Panel 2: Bed Power (original vs corrected)
    ax2 = fig.add_subplot(gs[1, 0], sharex=ax1)

    # Original bed power
    ax2.plot(distance, surfbed['bed_power_dB'].values,
             linewidth=1.2, color='gray', alpha=0.6, label='Original')

    # Geometric corrected power (convert to dB)
    P_geom_dB = 10 * np.log10(props['bed_power_geom_corrected'].values)
    ax2.plot(distance, P_geom_dB,
             linewidth=1.2, color='orange', label='Geom. Corrected')

    ax2.set_ylabel('Bed Power (dB)', fontsize=11, fontweight='bold')
    ax2.legend(loc='upper right', framealpha=0.9)
    ax2.grid(True, alpha=0.3)

    # Panel 3: Basal Reflectivity
    ax3 = fig.add_subplot(gs[2, 0], sharex=ax1)

    ax3.plot(distance, props['basal_reflectivity_dB'].values,
             linewidth=1.5, color='darkblue', marker='o', markersize=2)
    ax3.set_ylabel('Basal Reflectivity (dB)', fontsize=11, fontweight='bold')
    ax3.grid(True, alpha=0.3)

    # Add mean reflectivity annotation
    valid_refl = props['basal_reflectivity_dB'].values[~np.isnan(props['basal_reflectivity_dB'].values)]
    if len(valid_refl) > 0:
        mean_refl = np.mean(valid_refl)
        ax3.text(0.02, 0.95, f'Mean: {mean_refl:.1f} dB',
                 transform=ax3.transAxes, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Panel 4: Correlation Coefficient
    ax4 = fig.add_subplot(gs[3, 0], sharex=ax1)

    corr_coef = props['correlation_coef'].values
    ax4.plot(distance, corr_coef, linewidth=1.5, color='purple', marker='o', markersize=2)
    ax4.set_ylabel('|Correlation|\nCoefficient', fontsize=11, fontweight='bold')
    ax4.grid(True, alpha=0.3)

    # Add mean correlation annotation
    valid_corr = corr_coef[~np.isnan(corr_coef)]
    if len(valid_corr) > 0:
        mean_corr = np.mean(valid_corr)
        ax4.text(0.02, 0.95, f'Mean: {mean_corr:.3f}',
                 transform=ax4.transAxes, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Panel 5: Attenuation Rate
    ax5 = fig.add_subplot(gs[4, 0], sharex=ax1)

    # Plot attenuation rate
    valid = ~np.isnan(props['attenuation_rate'].values)
    ax5.plot(distance[valid], props['attenuation_rate'].values[valid],
             linewidth=1.5, color='darkred', marker='o', markersize=3)

    ax5.set_ylabel('Attenuation Rate\n(dB/km)', fontsize=11, fontweight='bold')
    ax5.grid(True, alpha=0.3)

    # Add statistics
    if np.sum(valid) > 0:
        atten_values = props['attenuation_rate'].values[valid]
        mean_atten = np.mean(atten_values)
        median_atten = np.median(atten_values)
        ax5.axhline(mean_atten, color='darkred', linestyle='--',
                   linewidth=1, alpha=0.5, label=f'Mean: {mean_atten:.1f}')
        ax5.text(0.02, 0.95,
                 f'Mean: {mean_atten:.1f} dB/km\nMedian: {median_atten:.1f} dB/km',
                 transform=ax5.transAxes, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Panel 6: Correlation Length
    ax6 = fig.add_subplot(gs[5, 0], sharex=ax1)

    valid = ~np.isnan(props['correlation_length'].values)
    ax6.plot(distance[valid], props['correlation_length'].values[valid],
             linewidth=1.5, color='darkgreen', marker='s', markersize=3)

    ax6.set_ylabel('Correlation Length\n(km)', fontsize=11, fontweight='bold')
    ax6.set_xlabel('Distance Along Track (km)', fontsize=11, fontweight='bold')
    ax6.grid(True, alpha=0.3)

    # Add statistics
    if np.sum(valid) > 0:
        corr_length_values = props['correlation_length'].values[valid]
        mean_length = np.mean(corr_length_values)
        ax6.text(0.02, 0.95, f'Mean: {mean_length:.1f} km',
                 transform=ax6.transAxes, verticalalignment='top',
                 bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Right panel: Map showing segment location
    ax_map = fig.add_subplot(gs[:, 1])

    # Plot the segment path
    ax_map.plot(lon, lat, linewidth=2, color='red', marker='o', markersize=3,
                markeredgecolor='darkred', alpha=0.8)

    # Mark start and end points
    ax_map.plot(lon[0], lat[0], 'go', markersize=10, label='Start', markeredgecolor='darkgreen', markeredgewidth=2)
    ax_map.plot(lon[-1], lat[-1], 'rs', markersize=10, label='End', markeredgecolor='darkred', markeredgewidth=2)

    ax_map.set_xlabel('Longitude', fontsize=11, fontweight='bold')
    ax_map.set_ylabel('Latitude', fontsize=11, fontweight='bold')
    ax_map.set_title('Segment Location', fontsize=12, fontweight='bold')
    ax_map.grid(True, alpha=0.3)
    ax_map.legend(loc='best', framealpha=0.9)
    ax_map.set_aspect('equal', adjustable='box')

    # Add lat/lon range annotation
    lat_range = f"{lat.min():.2f}° to {lat.max():.2f}°"
    lon_range = f"{lon.min():.2f}° to {lon.max():.2f}°"
    ax_map.text(0.02, 0.98, f'Lat: {lat_range}\nLon: {lon_range}',
                transform=ax_map.transAxes, verticalalignment='top',
                fontsize=9, bbox=dict(boxstyle='round', facecolor='white', alpha=0.8))

    # Remove x-tick labels for all but bottom panel on left
    for ax in [ax1, ax2, ax3, ax4, ax5]:
        ax.set_xticklabels([])

    plt.tight_layout()
    plt.savefig(output_plot, dpi=150, bbox_inches='tight')
    plt.close()

    print(f"[{datetime.datetime.now()}] Saved plot to {output_plot}", file=sys.stderr)

except Exception as e:
    error_msg = f"{type(e).__name__}: {str(e)}"
    print(f"[{datetime.datetime.now()}] Error creating plot: {error_msg}", file=sys.stderr)
    import traceback
    traceback.print_exc(file=sys.stderr)

    # Create error plot
    fig, ax = plt.subplots(figsize=(12, 6))
    ax.text(0.5, 0.5, f'Plot Creation Failed:\n{error_msg}',
            horizontalalignment='center', verticalalignment='center',
            transform=ax.transAxes, fontsize=14, color='red',
            bbox=dict(boxstyle='round', facecolor='white', edgecolor='red', linewidth=2))
    ax.set_title(f"Chunk {chunk_path}")
    ax.axis('off')
    plt.tight_layout()
    plt.savefig(output_plot, dpi=100)
    plt.close()

sys.stderr.close()
