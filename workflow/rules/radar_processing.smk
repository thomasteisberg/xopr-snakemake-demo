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


rule process_frame:
    """Process individual radar frame to extract reflection powers."""
    input:
        frames="results/search/frame_items.parquet",
    output:
        nc="results/processed_frames/frame_{frame_id}.nc",
        plot="results/processed_frames/frame_{frame_id}.png" if config["processing"]["create_debug_plots"] else [],
    params:
        cache_dir=config["opr"]["cache_dir"],
        data_product=config["processing"]["data_product"],
        resample_interval=config["processing"]["resample_interval"],
        layer_margin_m=config["processing"]["layer_margin_m"],
        frame_id=lambda wildcards: wildcards.frame_id,
        create_debug_plot=config["processing"]["create_debug_plots"],
        debug_plot_path=lambda wildcards: f"results/processed_frames/frame_{wildcards.frame_id}.png" if config["processing"]["create_debug_plots"] else None,
    conda:
        "../envs/xopr.yaml"
    log:
        "results/logs/process_frame_{frame_id}.log",
    script:
        "../scripts/process_frame.py"

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


rule create_results_map:
    """Create map of processed radar data."""
    input:
        region="results/region/selected_region.geojson",
        frame_items="results/search/frame_items.parquet",
        frames=lambda wildcards: expand(
            "results/processed_frames/frame_{frame_id}.nc",
            frame_id=get_frame_ids()),
    output:
        report("results/visualizations/results_map.html",
               caption="../report/results_map.rst",
               category="Visualizations",
               subcategory="Results"),
    params:
        projection=config["visualization"]["projection"],
    conda:
        "../envs/xopr_viz.yaml"
    log:
        "results/logs/create_results_map.log",
    script:
        "../scripts/create_results_map.py"