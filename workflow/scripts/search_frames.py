import json
import sys
import xopr.opr_access
import geopandas as gpd

# Snakemake inputs/outputs
input_region = snakemake.input.region
output_frames = snakemake.output.frames
output_summary = snakemake.output.summary
log_file = snakemake.log[0]

# Parameters
cache_dir = snakemake.params.cache_dir
max_items = snakemake.params.max_items

sys.stderr = open(log_file, 'w')

try:
    # Load region geometry (it's already GeoJSON)
    with open(input_region, 'r') as f:
        region_geometry = json.load(f)

    # Create OPR connection
    opr = xopr.opr_access.OPRConnection(cache_dir=cache_dir)

    # Search for frames
    print(f"Searching for frames intersecting region...", file=sys.stderr)
    stac_items = opr.query_frames(geometry=region_geometry)

    # Limit if specified
    if max_items:
        stac_items = stac_items.head(max_items)

    # Save as Parquet (preserves GeoDataFrame structure)
    stac_items.to_parquet(output_frames)

    # Create summary
    summary_lines = [
        f"Search Results Summary",
        f"=" * 50,
        f"Total frames found: {len(stac_items)}",
    ]

    if len(stac_items) > 0:
        collections = stac_items['collection'].unique()
        summary_lines.extend([
            f"Collections: {', '.join(sorted(collections))}",
            f"Frame IDs:",
        ])
        for _, frame in stac_items.iterrows():
            summary_lines.append(f"  - {frame['id']}")

    with open(output_summary, 'w') as f:
        f.write('\n'.join(summary_lines))

    print(f"Found {len(stac_items)} frames", file=sys.stderr)

except Exception as e:
    print(f"Error searching frames: {e}", file=sys.stderr)
    raise

sys.stderr.close()