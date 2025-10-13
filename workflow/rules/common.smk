# Rules for finding radar frames in OPR and selecting region

rule select_region:
    """Select geographic region from MEaSUREs Antarctic dataset."""
    output:
        geojson="results/region/selected_region.geojson",
        metadata="results/region/region_metadata.json",
    params:
        name=config["region"]["name"],
        type=config["region"]["type"],
        regions=config["region"]["regions"],
    conda:
        "../envs/xopr.yaml"
    log:
        "results/logs/select_region.log",
    script:
        "../scripts/select_region.py"


checkpoint search_frames:
    """Search for radar frames intersecting the selected region."""
    input:
        region="results/region/selected_region.geojson",
    output:
        frames="results/search/frame_items.parquet",
        summary="results/search/search_summary.txt",
        frame_ids="results/search/frame_ids.txt",
    params:
        cache_dir=config["opr"]["cache_dir"],
        max_items=config["search"]["max_items"],
    conda:
        "../envs/xopr.yaml"
    log:
        "results/logs/search_frames.log",
    script:
        "../scripts/search_frames.py"


checkpoint divide_segments:
    """Divide segments into continuous chunks (split on frame number gaps)."""
    input:
        frame_ids="results/search/frame_ids.txt",
    output:
        chunk_ids="results/search/chunk_ids.txt",
        chunk_frames="results/search/chunk_frames.txt",
    conda:
        "../envs/xopr.yaml"
    log:
        "results/logs/divide_segments.log",
    script:
        "../scripts/divide_segments.py"

# Visualization helpers

rule create_frame_map:
    """Create interactive map of frame coverage."""
    input:
        region="results/region/selected_region.geojson",
        frames="results/search/frame_items.parquet",
    output:
        report("results/visualizations/frame_coverage_map.html",
               caption="../report/frame_coverage_map.rst",
               category="Visualizations",
               subcategory="Frame Coverage"),
    params:
        projection=config["visualization"]["projection"],
    conda:
        "../envs/xopr_viz.yaml"
    log:
        "results/logs/create_frame_map.log",
    script:
        "../scripts/create_frame_map.py"


