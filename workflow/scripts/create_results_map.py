import json
import sys
import xarray as xr
import hvplot.xarray
import geoviews as gv
import geoviews.feature as gf
import cartopy.crs as ccrs
from shapely.geometry import shape
import xopr.geometry

# Snakemake inputs/outputs
input_region = snakemake.input.region
segment_files = snakemake.input.segments
status_files = snakemake.input.status_files
output_map = snakemake.output[0]
log_file = snakemake.log[0]

# Parameters
projection_str = snakemake.params.projection
color_field = snakemake.params.color_field
hover_fields = snakemake.params.get('hover_fields', [])
map_title = snakemake.params.get('title', 'Results Map')
color_label = snakemake.params.get('color_label', color_field)
input_pattern = snakemake.params.get('input_pattern', 'segment_surfbed_')
cmap = snakemake.params.get('cmap', 'turbo')
clim = snakemake.params.get('clim', None)

sys.stderr = open(log_file, 'w')

try:
    print(f"Creating results map...", file=sys.stderr)

    # Filter successful segments
    successful_segments = []
    failed_segments = []

    for segment_file, status_file in zip(segment_files, status_files):
        with open(status_file, 'r') as f:
            status = f.read().strip()

        if status == "success":
            successful_segments.append(segment_file)
        else:
            # Extract segment path from filename for logging
            segment_path = segment_file.split(input_pattern)[1].replace('.nc', '')
            failed_segments.append((segment_path, status))

    print(f"\nSegment processing summary:", file=sys.stderr)
    print(f"  Total segments: {len(segment_files)}", file=sys.stderr)
    print(f"  Successful: {len(successful_segments)}", file=sys.stderr)
    print(f"  Failed: {len(failed_segments)}", file=sys.stderr)

    if failed_segments:
        print(f"\nFailed segments:", file=sys.stderr)
        for segment_path, status in failed_segments:
            print(f"  - {segment_path}: {status}", file=sys.stderr)

    # Use only successful segments for visualization
    segment_files = successful_segments

    if len(segment_files) == 0:
        print(f"\nWarning: No successful segments to visualize!", file=sys.stderr)
        print(f"Creating empty map with region only...", file=sys.stderr)

    # Set up projection
    if projection_str == 'EPSG:3031':
        projection = ccrs.Stereographic(central_latitude=-90, true_scale_latitude=-71)
    elif projection_str == 'EPSG:3413':
        projection = ccrs.Stereographic(central_latitude=90, true_scale_latitude=70)
    else:
        projection = ccrs.PlateCarree()

    latlng = ccrs.PlateCarree()

    # Load region (it's already GeoJSON)
    with open(input_region, 'r') as f:
        region_geometry = json.load(f)

    # Create base map features
    features = (
        gf.ocean.options(scale='50m').opts(projection=projection) *
        gf.coastline.options(scale='50m').opts(projection=projection)
    )

    # Create region polygon - convert GeoJSON to shapely for GeoViews
    
    region_shape = shape(region_geometry)
    region_gv = gv.Polygons(region_shape, crs=latlng).opts(
        line_color='black',
        fill_alpha=0,
        projection=projection,
        scalebar=True
    )
    
    # Load segment datasets
    data_lines = []
    if len(segment_files) > 0:
        for ds_fn in segment_files:
            ds = xr.open_dataset(ds_fn)

            # Create derived fields if needed (before checking if field exists)
            if color_field == 'bed_minus_surf':
                if 'bed_power_dB' in ds and 'surface_power_dB' in ds:
                    ds['bed_minus_surf'] = ds['bed_power_dB'] - ds['surface_power_dB']
                else:
                    print(f"Warning: Cannot create bed_minus_surf in {ds_fn}, missing required fields", file=sys.stderr)
                    continue

            # Check if color field exists (after creating derived fields)
            if color_field not in ds:
                print(f"Warning: {color_field} not found in {ds_fn}, skipping", file=sys.stderr)
                continue

            ds = ds.dropna(dim='slow_time')
            ds = xopr.geometry.project_dataset(ds, target_crs='EPSG:3031')

            # Build hvplot kwargs
            plot_kwargs = {
                'x': 'x',
                'y': 'y',
                'c': color_field,
                'cmap': cmap,
                'size': 3,
            }

            # Add hover fields if specified
            if hover_fields:
                plot_kwargs['hover_cols'] = hover_fields

            # Add color limits if specified
            if clim is not None:
                plot_kwargs['clim'] = tuple(clim)

            # Add color label
            plot_kwargs['clabel'] = color_label

            sc = ds.hvplot.scatter(**plot_kwargs)
            data_lines.append(sc)

    # Combine all elements
    if data_lines:
        map_plot = features * region_gv * gv.Overlay(data_lines)
    else:
        map_plot = features * region_gv

    # Configure plot
    map_plot = map_plot.opts(
        title=map_title,
        projection=projection,
        width=800,
        height=600,
        aspect='equal'
    )

    # Save to HTML
    hvplot.save(map_plot, output_map)
    print(f"\nSaved results map to {output_map}", file=sys.stderr)

    # Print summary
    print(f"\nMap summary:", file=sys.stderr)
    print(f"  Region: Selected region", file=sys.stderr)
    print(f"  Number of visualized segments: {len(data_lines)}", file=sys.stderr)
    print(f"  Number of failed segments: {len(failed_segments)}", file=sys.stderr)
    print(f"  Projection: {projection_str}", file=sys.stderr)

except Exception as e:
    print(f"Error creating results map: {e}", file=sys.stderr)
    raise

sys.stderr.close()