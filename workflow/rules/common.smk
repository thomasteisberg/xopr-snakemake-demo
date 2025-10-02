# Common setup and validation for xOPR workflow

from snakemake.utils import validate

# validate config file
validate(config, schema="../../config/schemas/config.schema.yaml")


# Helper function to get frame IDs
def get_frame_ids():
    """Get frame IDs from search results."""
    import os

    frame_ids_file = "results/search/frame_ids.txt"

    if os.path.exists(frame_ids_file):
        with open(frame_ids_file, 'r') as f:
            return [line.strip() for line in f if line.strip()]
    else:
        return []