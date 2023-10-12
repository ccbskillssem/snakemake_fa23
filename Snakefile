# This is a sample snakefile for the workshop.
# The way that snakemake works is that it starts by reading the Snakefile, then checking the files that the Snakefile points to.
# Each listed .smk file is evaluated to identify available rules and their relationships.
# This is why the (optional) rule all directive is placed in the Snakefile rather than in a separate .smk file: snakemake will read the Snakefile first, meaning "rule all" is triggered with highest priority.

configfile: "config/config.yml"

include: "rules/common_utils.smk"
include: "rules/get_data.smk"
include: "rules/preprocessing.smk"
include: "rules/mapping.smk"
include: "rules/call_variants.smk"
include: "rules/plotting.smk"
