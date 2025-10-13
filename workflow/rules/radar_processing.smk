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

rule estimate_basal_properties:
    """Estimate basal reflectivity and attenuation rate from surface/bed power.

    Uses geometric spreading correction and adaptive windowing to decorrelate
    bed power from ice thickness, following Matsuoka et al. (2010) and
    Schroeder et al. (2016) methodologies.
    """
    input:
        surfbed="results/processed_segments/segment_surfbed_{chunk_path}.nc",
    output:
        nc="results/basal_properties/basal_props_{chunk_path}.nc",
        status="results/basal_properties/basal_props_{chunk_path}.status",
    params:
        chunk_path=lambda wildcards: wildcards.chunk_path,
        ice_permittivity=config["attenuation"]["ice_permittivity"],
        min_window_size_km=config["attenuation"]["min_window_size_km"],
        max_window_size_km=config["attenuation"]["max_window_size_km"],
        window_step_km=config["attenuation"]["window_step_km"],
        alpha_range_db_km=config["attenuation"]["alpha_range_db_km"],
        alpha_step_db_km=config["attenuation"]["alpha_step_db_km"],
        correlation_threshold=config["attenuation"]["correlation_threshold"],
    conda:
        "../envs/xopr.yaml"
    log:
        "results/logs/basal_props_{chunk_path}.log",
    script:
        "../scripts/estimate_basal_properties.py"


rule plot_basal_properties:
    """Create visualization of basal properties analysis."""
    input:
        basal_props="results/basal_properties/basal_props_{chunk_path}.nc",
        surfbed="results/processed_segments/segment_surfbed_{chunk_path}.nc",
    output:
        png="results/basal_properties/basal_props_{chunk_path}.png",
    params:
        chunk_path=lambda wildcards: wildcards.chunk_path,
    conda:
        "../envs/xopr.yaml"
    log:
        "results/logs/plot_basal_props_{chunk_path}.log",
    script:
        "../scripts/plot_basal_properties.py"


rule create_bed_power_map:
    """Create map of bed power results."""
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
        report("results/visualizations/bed_power_map.html",
               caption="../report/bed_power_map.rst",
               category="Visualizations",
               subcategory="Bed Power"),
    params:
        projection=config["visualization"]["projection"],
        color_field="bed_minus_surf",
        hover_fields=["surface_power_dB", "bed_power_dB"],
        title="Bed Power Results (Bed - Surface)",
        color_label="Bed - Surface Power (dB)",
        input_pattern="segment_surfbed_",
        cmap="turbo",
    conda:
        "../envs/xopr_viz.yaml"
    log:
        "results/logs/create_bed_power_map.log",
    script:
        "../scripts/create_results_map.py"


rule create_basal_reflectivity_map:
    """Create map of basal reflectivity results."""
    input:
        region="results/region/selected_region.geojson",
        frame_items="results/search/frame_items.parquet",
        chunk_frames="results/search/chunk_frames.txt",
        segments=lambda wildcards: expand(
            "results/basal_properties/basal_props_{chunk_path}.nc",
            chunk_path=get_chunk_paths()),
        status_files=lambda wildcards: expand(
            "results/basal_properties/basal_props_{chunk_path}.status",
            chunk_path=get_chunk_paths()),
    output:
        report("results/visualizations/basal_reflectivity_map.html",
               caption="../report/basal_reflectivity_map.rst",
               category="Visualizations",
               subcategory="Basal Reflectivity"),
    params:
        projection=config["visualization"]["projection"],
        color_field="basal_reflectivity_dB",
        hover_fields=["ice_thickness", "attenuation_rate", "basal_reflectivity_dB"],
        title="Basal Reflectivity",
        color_label="Basal Reflectivity (dB)",
        input_pattern="basal_props_",
        cmap="viridis",
        clim=[0, 100],
    conda:
        "../envs/xopr_viz.yaml"
    log:
        "results/logs/create_basal_reflectivity_map.log",
    script:
        "../scripts/create_results_map.py"


rule create_attenuation_map:
    """Create map of englacial attenuation rates."""
    input:
        region="results/region/selected_region.geojson",
        frame_items="results/search/frame_items.parquet",
        chunk_frames="results/search/chunk_frames.txt",
        segments=lambda wildcards: expand(
            "results/basal_properties/basal_props_{chunk_path}.nc",
            chunk_path=get_chunk_paths()),
        status_files=lambda wildcards: expand(
            "results/basal_properties/basal_props_{chunk_path}.status",
            chunk_path=get_chunk_paths()),
    output:
        report("results/visualizations/attenuation_map.html",
               caption="../report/attenuation_map.rst",
               category="Visualizations",
               subcategory="Attenuation"),
    params:
        projection=config["visualization"]["projection"],
        color_field="attenuation_rate",
        hover_fields=["ice_thickness", "attenuation_rate", "correlation_length", "correlation_coef"],
        title="Englacial Attenuation Rate",
        color_label="Attenuation Rate (dB/km)",
        input_pattern="basal_props_",
        cmap="plasma",
        clim=[5, 30],
    conda:
        "../envs/xopr_viz.yaml"
    log:
        "results/logs/create_attenuation_map.log",
    script:
        "../scripts/create_results_map.py"