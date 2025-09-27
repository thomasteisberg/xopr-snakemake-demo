import sys
import xopr

cache_dir = snakemake.params.cache_dir
output_file = snakemake.output[0]
log_file = snakemake.log[0]

sys.stderr = open(log_file, 'w')

try:
    opr = xopr.OPRConnection(cache_dir=cache_dir)
    collections = opr.get_collections()
    
    with open(output_file, 'w') as f:
        f.write("Available OPR Collections:\n")
        f.write("=" * 50 + "\n")
        for collection in collections:
            f.write(f"{collection}\n")
        f.write("=" * 50 + "\n")
        f.write(f"Total collections: {len(collections)}\n")
    
    print(f"Successfully wrote {len(collections)} collections to {output_file}", file=sys.stderr)
    
except Exception as e:
    print(f"Error: {e}", file=sys.stderr)
    raise

sys.stderr.close()