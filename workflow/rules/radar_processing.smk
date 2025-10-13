# Rules for radar data processing workflow

rule process_frame:
    """Process individual radar frame to extract reflection powers."""
    input:
        frames=ancient("results/search/frame_items.parquet"),
    output:
        nc="results/processed_frames/frame_{frame_id}.nc",
        status="results/processed_frames/frame_{frame_id}.status",
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


rule extract_surface_bed_power:
    """Process continuous chunk to extract reflection powers.

    Chunks are continuous segments of frames without gaps.
    This rule merges all frames from a chunk and extracts surface/bed power.
    """
    input:
        frames=ancient("results/search/frame_items.parquet"),
        chunk_frames="results/search/chunk_frames.txt",
        chunk_ids="results/search/chunk_ids.txt",
    output:
        nc="results/processed_segments/segment_surfbed_{chunk_path}.nc",
        status="results/processed_segments/segment_surfbed_{chunk_path}.status",
        plot="results/processed_segments/segment_surfbed_{chunk_path}.png" if config["processing"]["create_debug_plots"] else [],
    params:
        cache_dir=config["opr"]["cache_dir"],
        data_product=config["processing"]["data_product"],
        resample_interval=config["processing"]["resample_interval"],
        layer_margin_m=config["processing"]["layer_margin_m"],
        chunk_path=lambda wildcards: wildcards.chunk_path,
        create_debug_plot=config["processing"]["create_debug_plots"],
        debug_plot_path=lambda wildcards: f"results/processed_segments/segment_surfbed_{wildcards.chunk_path}.png" if config["processing"]["create_debug_plots"] else None,
    conda:
        "../envs/xopr.yaml"
    log:
        "results/logs/extract_segment_surfbed_{chunk_path}.log",
    script:
        "../scripts/extract_segment.py"

rule create_results_map:
    """Create map of processed radar data."""
    input:
        region="results/region/selected_region.geojson",
        frame_items="results/search/frame_items.parquet",
        chunk_frames="results/search/chunk_frames.txt",
        segments=lambda wildcards: expand(
            "results/processed_segments/segment_surfbed_{chunk_path}.nc",
            chunk_path=get_chunk_paths()),
        status_files=lambda wildcards: expand(
            "results/processed_segments/segment_surfbed_{chunk_path}.status",
            chunk_path=get_chunk_paths()),
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