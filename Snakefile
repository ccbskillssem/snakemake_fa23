# To run this workflow through the mapping step for one sample, use the command:
# snakemake -pj{n} --use-conda output/bwa/SRR23032907.sorted.bam.bai
# You can adjust n according to your computing resources.

# To run this workflow through the final variant calling step, use the command:
# snakemake -pj{n} --use-conda output/bcftools/all_samples.vcf

# To create the rule graph ("simplified" graph) of this workflow, use one of the following:
# snakemake --rulegraph --use-conda output/bcftools/all_samples.vcf | dot -Tpng > output/workflow_rulegraph.png
# snakemake --rulegraph --use-conda output/bcftools/all_samples.vcf | dot -Tsvg > output/workflow_rulegraph.svg

# To create the full graph (showing all samples) of this workflow, use one of the following:
# snakemake --dag --use-conda output/bcftools/all_samples.vcf | dot -Tpng > output/workflow_dag.png
# snakemake --dag --use-conda output/bcftools/all_samples.vcf | dot -Tsvg > output/workflow_dag.svg

configfile: "config/config.yml"

include: "rules/common_utils.smk"
include: "rules/get_data.smk"
include: "rules/preprocessing.smk"
include: "rules/mapping.smk"
include: "rules/call_variants.smk"