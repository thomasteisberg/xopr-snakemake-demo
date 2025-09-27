# Rules for radar data processing workflow

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


rule search_frames:
    """Search for radar frames intersecting the selected region."""
    input:
        region="results/region/selected_region.geojson",
    output:
        frames="results/search/frame_items.json",
        summary="results/search/search_summary.txt",
    params:
        cache_dir=config["opr"]["cache_dir"],
        stac_api_url=config["opr"]["stac_api_url"],
        max_items=config["search"]["max_items"],
    conda:
        "../envs/xopr.yaml"
    log:
        "results/logs/search_frames.log",
    script:
        "../scripts/search_frames.py"


rule process_frame:
    """Process individual radar frame to extract reflection powers."""
    input:
        frames="results/search/frame_items.json",
    output:
        "results/processed_frames/frame_{frame_idx}.nc",
    params:
        cache_dir=config["opr"]["cache_dir"],
        data_product=config["processing"]["data_product"],
        resample_interval=config["processing"]["resample_interval"],
        layer_margin_m=config["processing"]["layer_margin_m"],
        frame_idx=lambda wildcards: int(wildcards.frame_idx),
    conda:
        "../envs/xopr.yaml"
    log:
        "results/logs/process_frame_{frame_idx}.log",
    script:
        "../scripts/process_frame.py"


rule merge_frames:
    """Merge all processed frames into single dataset."""
    input:
        frames=lambda wildcards: expand(
            "results/processed_frames/frame_{frame_idx}.nc",
            frame_idx=range(get_frame_count())
        ),
    output:
        "results/merged/merged_dataset.nc",
    params:
        format=config["output"]["format"],
    conda:
        "../envs/xopr.yaml"
    log:
        "results/logs/merge_frames.log",
    script:
        "../scripts/merge_frames.py"


rule create_frame_map:
    """Create interactive map of frame coverage."""
    input:
        region="results/region/selected_region.geojson",
        frames="results/search/frame_items.json",
    output:
        "results/visualizations/frame_coverage_map.html",
    params:
        projection=config["visualization"]["projection"],
    conda:
        "../envs/xopr_viz.yaml"
    log:
        "results/logs/create_frame_map.log",
    script:
        "../scripts/create_frame_map.py"


rule create_results_map:
    """Create map of processed radar data."""
    input:
        region="results/region/selected_region.geojson",
        data="results/merged/merged_dataset.nc" if config["output"]["merge_frames"] else [],
        frames=lambda wildcards: expand(
            "results/processed_frames/frame_{frame_idx}.nc",
            frame_idx=range(get_frame_count())
        ) if not config["output"]["merge_frames"] else [],
    output:
        "results/visualizations/results_map.html",
    params:
        projection=config["visualization"]["projection"],
    conda:
        "../envs/xopr_viz.yaml"
    log:
        "results/logs/create_results_map.log",
    script:
        "../scripts/create_results_map.py"


# Helper function to get frame count
def get_frame_count():
    """Get number of frames to process from search results."""
    import json
    import os
    
    frames_file = "results/search/frame_items.json"
    if os.path.exists(frames_file):
        with open(frames_file) as f:
            frames = json.load(f)
        return len(frames)
    return 0