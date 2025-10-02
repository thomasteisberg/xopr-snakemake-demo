import sys
import datetime
import numpy as np
import xarray as xr
import pandas as pd
import geopandas as gpd
import scipy.constants
import xopr.opr_access
import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt

# Snakemake inputs/outputs
input_frames = snakemake.input.frames
output_file = snakemake.output.nc if hasattr(snakemake.output, 'nc') else snakemake.output[0]
log_file = snakemake.log[0]

# Parameters
cache_dir = snakemake.params.cache_dir
data_product = snakemake.params.data_product
resample_interval = snakemake.params.resample_interval
layer_margin_m = snakemake.params.layer_margin_m
frame_id = snakemake.params.frame_id
create_debug_plot = snakemake.params.get('create_debug_plot', False)
debug_plot_path = snakemake.params.get('debug_plot_path', None)

sys.stderr = open(log_file, 'w')


def extract_layer_peak_power(radar_ds, layer_twtt, margin_twtt):
    """Extract peak power within margin around layer."""
    # Align time coordinates
    t_start = np.minimum(radar_ds.slow_time.min(), layer_twtt.slow_time.min())
    t_end = np.maximum(radar_ds.slow_time.max(), layer_twtt.slow_time.max())
    layer_twtt = layer_twtt.sel(slow_time=slice(t_start, t_end))
    radar_ds = radar_ds.sel(slow_time=slice(t_start, t_end))
    layer_twtt = layer_twtt.reindex(
        slow_time=radar_ds.slow_time,
        method='nearest',
        tolerance=pd.Timedelta(seconds=1),
        fill_value=np.nan
    )

    # Extract data within margin
    start_twtt = layer_twtt - margin_twtt
    end_twtt = layer_twtt + margin_twtt
    data_within_margin = radar_ds.where(
        (radar_ds.twtt >= start_twtt) & (radar_ds.twtt <= end_twtt),
        drop=True
    )

    # Calculate power in dB
    power_dB = 10 * np.log10(np.abs(data_within_margin.Data))

    # Find peak
    peak_twtt_index = power_dB.argmax(dim='twtt')
    peak_twtt = power_dB.twtt[peak_twtt_index]
    peak_power = power_dB.isel(twtt=peak_twtt_index)

    # Clean up dimensions
    peak_twtt = peak_twtt.drop_vars('twtt')
    peak_power = peak_power.drop_vars('twtt')

    return peak_twtt, peak_power


try:
    # Load frame list from Parquet
    stac_items = gpd.read_parquet(input_frames)

    stac_item = stac_items.loc[frame_id]

    print(f"[{datetime.datetime.now()}] Processing frame {frame_id}: {stac_item['id']}", file=sys.stderr)

    # Create OPR connection (each parallel process gets its own connection)
    opr = xopr.opr_access.OPRConnection(cache_dir=cache_dir)

    # Load frame
    frame = opr.load_frame(stac_item, data_product=data_product)

    # Resample
    frame = frame.sortby('slow_time')
    frame = frame.resample(slow_time=resample_interval).mean()

    # Get layers
    try:
        layers = opr.get_layers_db(frame, include_geometry=False)
    except Exception as e:
        print(f"Warning: Could not retrieve layers: {e}", file=sys.stderr)
        # Create empty dataset with metadata
        reflectivity_dataset = xr.Dataset({
            'surface_twtt': (('slow_time',), np.full(len(frame.slow_time), np.nan)),
            'bed_twtt': (('slow_time',), np.full(len(frame.slow_time), np.nan)),
            'surface_power_dB': (('slow_time',), np.full(len(frame.slow_time), np.nan)),
            'bed_power_dB': (('slow_time',), np.full(len(frame.slow_time), np.nan)),
        })
        reflectivity_dataset.coords['slow_time'] = frame.slow_time
        reflectivity_dataset.coords['Latitude'] = frame.Latitude
        reflectivity_dataset.coords['Longitude'] = frame.Longitude
        reflectivity_dataset.attrs['frame_id'] = stac_item['id']
        reflectivity_dataset.attrs['processing_error'] = 'Could not retrieve layers'
        reflectivity_dataset.to_netcdf(output_file)
        sys.stderr.close()
        sys.exit(0)

    # Calculate margin in TWTT
    speed_of_light_in_ice = scipy.constants.c / np.sqrt(3.17)
    layer_selection_margin_twtt = layer_margin_m / speed_of_light_in_ice

    # Extract surface and bed powers
    surface_twtt, surface_power = extract_layer_peak_power(
        frame, layers[1]['twtt'], layer_selection_margin_twtt
    )
    bed_twtt, bed_power = extract_layer_peak_power(
        frame, layers[2]['twtt'], layer_selection_margin_twtt
    )

    # Create output dataset
    reflectivity_dataset = xr.Dataset({
        'surface_twtt': surface_twtt,
        'bed_twtt': bed_twtt,
        'surface_power_dB': surface_power,
        'bed_power_dB': bed_power,
    })

    # Add coordinates (use capitalized names as in xOPR)
    reflectivity_dataset.coords['Latitude'] = frame.Latitude
    reflectivity_dataset.coords['Longitude'] = frame.Longitude
    if 'Elevation' in frame:
        reflectivity_dataset.coords['Elevation'] = frame.Elevation

    # Add metadata
    reflectivity_dataset.attrs['frame_id'] = stac_item['id']
    reflectivity_dataset.attrs['collection'] = stac_item['collection']
    reflectivity_dataset.attrs['data_product'] = data_product
    reflectivity_dataset.attrs['resample_interval'] = resample_interval
    reflectivity_dataset.attrs['frame_id'] = frame_id

    # Copy frame attributes
    for attr in ['season', 'segment', 'doi', 'ror', 'funder_text']:
        if attr in frame.attrs:
            reflectivity_dataset.attrs[attr] = frame.attrs[attr]

    # Save
    reflectivity_dataset.to_netcdf(output_file)
    print(f"[{datetime.datetime.now()}] Saved processed frame to {output_file}", file=sys.stderr)

    # Create debug plot if requested
    if create_debug_plot and debug_plot_path:
        fig, ax = plt.subplots(figsize=(10, 6))
        reflectivity_dataset['surface_power_dB'].plot(ax=ax, x='slow_time', label='Surface', linewidth=2)
        reflectivity_dataset['bed_power_dB'].plot(ax=ax, x='slow_time', label='Bed', linewidth=2)
        ax.set_ylabel('Power [dB]')
        ax.set_title(f"Frame {stac_item['id']}: Surface and Bed Power")
        ax.legend()
        ax.grid(True, alpha=0.3)
        plt.tight_layout()
        plt.savefig(debug_plot_path, dpi=100)
        plt.close()
        print(f"[{datetime.datetime.now()}] Saved debug plot to {debug_plot_path}", file=sys.stderr)

except Exception as e:
    print(f"[{datetime.datetime.now()}] Error processing frame {frame_id}: {e}", file=sys.stderr)
    raise

sys.stderr.close()