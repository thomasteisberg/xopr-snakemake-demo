import json
import sys
import geopandas as gpd
import holoviews as hv
import hvplot.pandas
import geoviews as gv
import geoviews.feature as gf
import cartopy.crs as ccrs
from shapely.geometry import shape
import xopr

# Snakemake inputs/outputs
input_region = snakemake.input.region
input_frames = snakemake.input.frames
output_map = snakemake.output[0]
log_file = snakemake.log[0]

# Parameters
projection_str = snakemake.params.projection

sys.stderr = open(log_file, 'w')

try:
    print(f"Creating frame coverage map...", file=sys.stderr)

    # Set up projection
    if projection_str == 'EPSG:3031':
        projection = ccrs.Stereographic(central_latitude=-90, true_scale_latitude=-71)
    elif projection_str == 'EPSG:3413':
        projection = ccrs.Stereographic(central_latitude=90, true_scale_latitude=70)
    else:
        projection = ccrs.PlateCarree()

    latlng = ccrs.PlateCarree()

    # Create base map features
    features = (
        gf.ocean.options(scale='50m').opts(projection=projection) *
        gf.coastline.options(scale='50m').opts(projection=projection)
    )

    # Load region and create a polygon
    with open(input_region, 'r') as f:
        region_geometry = json.load(f)
    
    region_polygon = shape(region_geometry)
    region_projected = xopr.geometry.project_geojson(region_polygon, target_crs=projection_str)


    region_hv = hv.Polygons([region_projected]).opts(
        color='green',
        line_color='black',
        fill_alpha=0.5)


    # Load frames from Parquet
    stac_items = gpd.read_parquet(input_frames)
    stac_items = stac_items.to_crs(projection_str)

    # Create frame paths
    flight_lines = stac_items.hvplot(by='collection')

    # Combine all elements
    map_plot = features * region_hv * gv.Overlay(flight_lines)

    # Configure plot
    map_plot = map_plot.opts(
        title=f"Frame Coverage Map",
        projection=projection,
        width=800,
        height=600,
        aspect='equal'
    )

    # Save to HTML
    hvplot.save(map_plot, output_map)
    print(f"Saved frame coverage map to {output_map}", file=sys.stderr)

    # Print summary
    print(f"\nMap summary:", file=sys.stderr)
    print(f"  Region: Selected region", file=sys.stderr)
    print(f"  Total frames: {len(stac_items)}", file=sys.stderr)
    print(f"  Projection: {projection_str}", file=sys.stderr)

except Exception as e:
    print(f"Error creating frame map: {e}", file=sys.stderr)
    raise

sys.stderr.close()