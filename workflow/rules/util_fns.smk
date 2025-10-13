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


def get_chunk_paths():
    """Get chunk paths from divide_segments checkpoint.

    Chunks are continuous segments of frames. Segments with gaps in frame numbers
    are split into multiple chunks (e.g., Data_20131120_01_chunk1, Data_20131120_01_chunk2).
    """
    import os

    chunk_ids_file = "results/search/chunk_ids.txt"

    if os.path.exists(chunk_ids_file):
        with open(chunk_ids_file, 'r') as f:
            return [line.strip() for line in f if line.strip()]
    else:
        return []