# Snakemake workflow: xOPR Radar Processing

[![Snakemake](https://img.shields.io/badge/snakemake-≥8.0.0-brightgreen.svg)](https://snakemake.github.io)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049.svg?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)

A Snakemake workflow demonstrating the use of xOPR (Open Polar Radar) for ice-penetrating radar data analysis. This workflow replicates the analysis from the xOPR `search_and_scaling` notebook in a reproducible, parallelizable pipeline.

## Overview

This workflow includes the following steps:
1. **Region Selection**: Select a geographic region of interest from Antarctica
2. **Frame Search**: Find radar frames intersecting the selected region
3. **Parallel Processing**: Extract surface and bed reflection powers from each frame
5. **Visualization**: Create interactive maps of frame coverage and results

## Table of Contents
- [Installation](#installation)
- [Usage](#usage)
- [Configuration](#configuration)
- [Output Structure](#output-structure)
- [Parallelization](#parallelization)
- [Authors](#authors)
- [References](#references)

## Installation

1. Clone this repository:
```bash
git clone https://github.com/thomasteisberg/xopr-snakemake-demo
cd xopr-snakemake-demo
```

2. Setup a conda/mamba environment (if not already installed):
```bash
mamba env create -n snakemake -f environment.yaml
mamba activate snakemake
```

## Usage

### Quick Start

```bash
# Activate snakemake environment
mamba activate snakemake

# Dry run to preview execution plan
snakemake --dry-run

# Run complete workflow (use multiple cores for parallel processing)
snakemake --cores 4 --sdm conda
```

### Other options:

**Search only workflow:**

This workflow runs only the geographic select and frame search parts of the workflow, 
producing a map of where data has been found without actually processing any of it.

```bash
# Only search for frames (no processing)
snakemake search_only --sdm conda
```

**Clearing results:**

This will delete all results files that would be created by the current workflow:

```bash
# Delete current results files
snakemake --delete-all-output
```

**Workflow visualization:**

Generate visual representations of the workflow:

```bash
# Generate workflow DAG
snakemake --dag | dot -Tpng > workflow_dag.png

# Generate rule graph
snakemake --rulegraph | dot -Tpng > workflow_rules.png
```

## Configuration

Edit `config/config.yaml` to customize your analysis:

```yaml
# OPR connection settings
opr:
  cache_dir: "/path/to/cache"  # Local cache for downloaded data

# Region selection (MEaSUREs Antarctic dataset)
region:
  name: "Dotson"    # Region name (e.g., "Dotson", "George_VI", "LarsenC")
  type: "FL"        # Type: FL (floating), GR (grounded), IS (island)
  regions: null     # Geographic regions: Peninsula, West, East

# Frame search parameters
search:
  max_items: 10     # Maximum frames to process (null for all)

# Processing parameters
processing:
  data_product: "CSARP_standard"  # Radar data product
  resample_interval: "5s"         # Temporal resampling
  layer_margin_m: 50               # Margin around layers (meters)

# Output settings
output:
  format: "netcdf"      # Output format: netcdf or zarr

# Visualization settings
visualization:
  create_maps: true         # Create interactive maps
  create_plots: true        # Create static plots
  projection: "EPSG:3031"   # Map projection (3031: Antarctic, 3413: Greenland)
```

## Output Structure

The workflow generates the following output structure:

```
results/
├── region/                 # Selected region
│   ├── selected_region.geojson
│   └── region_metadata.json
├── search/                 # Frame search results
│   ├── frame_items.json
│   └── search_summary.txt
├── processed_frames/       # Individual processed frames
│   ├── frame_0.nc
│   ├── frame_1.nc
│   └── ...
├── merged/                 # Combined dataset
│   └── merged_dataset.nc
└── visualizations/         # Interactive HTML maps
    ├── frame_coverage_map.html
    └── results_map.html
```

## References

- [xOPR Package](https://github.com/thomasteisberg/xopr) - Python library for Open Polar Radar data access
- [Open Polar Radar](https://data.openradardata.org/) - Ice-penetrating radar data repository
