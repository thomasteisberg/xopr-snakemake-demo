# Common setup and validation for xOPR workflow

from snakemake.utils import validate

# validate config file
validate(config, schema="../../config/schemas/config.schema.yaml")