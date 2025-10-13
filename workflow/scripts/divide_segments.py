import sys
from collections import defaultdict

# Snakemake inputs/outputs
input_frame_ids = snakemake.input.frame_ids
output_chunk_ids = snakemake.output.chunk_ids
output_chunk_frames = snakemake.output.chunk_frames
log_file = snakemake.log[0]

sys.stderr = open(log_file, 'w')

print("Dividing segments into continuous chunks...", file=sys.stderr)

# Read frame IDs
with open(input_frame_ids, 'r') as f:
    frame_ids = [line.strip() for line in f if line.strip()]

print(f"Read {len(frame_ids)} frame IDs", file=sys.stderr)

# Group frames by segment path (all but last part)
segments = defaultdict(list)
for frame_id in frame_ids:
    parts = frame_id.split('_')
    if len(parts) >= 2:
        segment_path = '_'.join(parts[:-1])
        frame_num = int(parts[-1])
        segments[segment_path].append((frame_num, frame_id))

print(f"Found {len(segments)} unique segments", file=sys.stderr)

# For each segment, detect gaps and create chunks
chunks = {}  # chunk_id -> list of frame_ids

for segment_path, frames in segments.items():
    # Sort by frame number
    frames.sort(key=lambda x: x[0])

    # Detect gaps (frame numbers that differ by more than 1)
    chunk_groups = []
    current_chunk = [frames[0]]

    for i in range(1, len(frames)):
        prev_num = frames[i-1][0]
        curr_num = frames[i][0]

        # If gap detected (difference > 1)
        if curr_num - prev_num > 1:
            chunk_groups.append(current_chunk)
            current_chunk = [frames[i]]
        else:
            current_chunk.append(frames[i])

    # Add the last chunk
    chunk_groups.append(current_chunk)

    # Name chunks
    if len(chunk_groups) == 1:
        # No gaps, use original segment name
        chunk_id = segment_path
        chunks[chunk_id] = [f[1] for f in chunk_groups[0]]
        print(f"  {segment_path}: {len(frames)} frames, continuous", file=sys.stderr)
    else:
        # Has gaps, add chunk suffixes
        for idx, chunk_group in enumerate(chunk_groups, start=1):
            chunk_id = f"{segment_path}_chunk{idx}"
            chunks[chunk_id] = [f[1] for f in chunk_group]
            frame_range = f"{chunk_group[0][0]}-{chunk_group[-1][0]}"
            print(f"  {segment_path}: split into {len(chunk_groups)} chunks", file=sys.stderr)
            print(f"    {chunk_id}: frames {frame_range} ({len(chunk_group)} frames)", file=sys.stderr)

print(f"\nTotal chunks created: {len(chunks)}", file=sys.stderr)

# Write chunk IDs (one per line)
with open(output_chunk_ids, 'w') as f:
    for chunk_id in sorted(chunks.keys()):
        f.write(f"{chunk_id}\n")

# Write chunk->frames mapping (for the script to use)
with open(output_chunk_frames, 'w') as f:
    for chunk_id in sorted(chunks.keys()):
        frame_list = ','.join(chunks[chunk_id])
        f.write(f"{chunk_id}\t{frame_list}\n")

print(f"Written chunk IDs to {output_chunk_ids}", file=sys.stderr)
print(f"Written chunk frames mapping to {output_chunk_frames}", file=sys.stderr)

sys.stderr.close()
