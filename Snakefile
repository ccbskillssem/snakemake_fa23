configfile: "config/config.yml"

include: "rules/common_utils.smk"
include: "rules/get_data.smk"
include: "rules/preprocessing.smk"
include: "rules/mapping.smk"
include: "rules/call_variants.smk"
include: "rules/plotting.smk"
