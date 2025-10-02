import json
import sys
import xopr.geometry

# Snakemake inputs/outputs
output_geojson = snakemake.output.geojson
output_metadata = snakemake.output.metadata
log_file = snakemake.log[0]

# Parameters
region_name = snakemake.params.name
region_type = snakemake.params.type
regions = snakemake.params.regions

sys.stderr = open(log_file, 'w')

try:
    # Select region from MEaSUREs dataset - returns Shapely geometry when merge_regions=True
    region_geometry = xopr.geometry.get_antarctic_regions(
        name=region_name,
        type=region_type,
        regions=regions,
        merge_regions=True,
        simplify_tolerance=100
    )

    # Get GeoJSON representation using __geo_interface__
    region_geojson = region_geometry.__geo_interface__

    # Save region geometry as GeoJSON
    with open(output_geojson, 'w') as f:
        json.dump(region_geojson, f, indent=2)

    # Calculate area and save metadata
    projected = xopr.geometry.project_geojson(region_geometry, target_crs='EPSG:3031')
    area_km2 = projected.area / 1e6

    metadata = {
        "name": region_name,
        "type": region_type,
        "regions": regions,
        "area_km2": round(area_km2, 2)
    }

    with open(output_metadata, 'w') as f:
        json.dump(metadata, f, indent=2)

    print(f"Selected region: {metadata['name']}", file=sys.stderr)
    print(f"Area: {metadata['area_km2']} km²", file=sys.stderr)

except Exception as e:
    print(f"Error selecting region: {e}", file=sys.stderr)
    raise

sys.stderr.close()