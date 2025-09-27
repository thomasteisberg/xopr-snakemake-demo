rule list_collections:
    output:
        "results/opr_collections/collections.txt",
    conda:
        "../envs/xopr.yaml"
    params:
        cache_dir=config.get("cache_dir", "/home/thomasteisberg/Documents/opr/xopr/docs/notebooks/radar_cache"),
    log:
        "results/logs/list_collections.log",
    script:
        "../scripts/list_collections.py"