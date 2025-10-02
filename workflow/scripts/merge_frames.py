import sys
import xarray as xr

# Snakemake inputs/outputs
input_frames = snakemake.input.frames
output_file = snakemake.output[0]
log_file = snakemake.log[0]

# Parameters
output_format = snakemake.params.format

sys.stderr = open(log_file, 'w')

try:
    print(f"Merging {len(input_frames)} processed frames...", file=sys.stderr)

    # Load all frame datasets
    datasets = []
    for frame_file in input_frames:
        ds = xr.open_dataset(frame_file)
        # Add frame identifier as coordinate
        if 'frame_id' in ds.attrs:
            ds = ds.assign_coords(frame_id=ds.attrs['frame_id'])
        datasets.append(ds)
        print(f"Loaded {frame_file}", file=sys.stderr)

    # Merge datasets along new dimension
    # Use compat='override' to handle different coordinate values
    merged = xr.concat(datasets, dim='frame', coords='minimal', compat='override')

    # Preserve global attributes from first frame
    if datasets:
        global_attrs = {}
        # Copy common attributes
        for attr in ['data_product', 'resample_interval', 'doi', 'ror', 'funder_text']:
            if attr in datasets[0].attrs:
                global_attrs[attr] = datasets[0].attrs[attr]

        # Add merge metadata
        global_attrs['n_frames'] = len(datasets)
        global_attrs['frame_ids'] = [ds.attrs.get('frame_id', 'unknown') for ds in datasets]
        merged.attrs = global_attrs

    # Save based on format
    if output_format == 'zarr':
        merged.to_zarr(output_file)
        print(f"Saved merged dataset to {output_file} (zarr format)", file=sys.stderr)
    else:
        merged.to_netcdf(output_file)
        print(f"Saved merged dataset to {output_file} (netcdf format)", file=sys.stderr)

    # Print summary
    print(f"\nMerge summary:", file=sys.stderr)
    print(f"  Total frames: {len(datasets)}", file=sys.stderr)
    print(f"  Total points: {merged.sizes['slow_time'] * merged.sizes['frame']}", file=sys.stderr)
    print(f"  Variables: {list(merged.data_vars)}", file=sys.stderr)

except Exception as e:
    print(f"Error merging frames: {e}", file=sys.stderr)
    raise

sys.stderr.close()