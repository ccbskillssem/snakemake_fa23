# To run this workflow through the final variant calling step + create plots, use the command:
# snakemake -pj{n} --use-conda output/visuals/vcf_heatmap.pdf output/visuals/sample_coverage.pdf
# You can adjust n according to your computing resources.

configfile: "config/config.yml"

include: "rules/common_utils.smk"
include: "rules/get_data.smk"
include: "rules/preprocessing.smk"
include: "rules/mapping.smk"
include: "rules/call_variants.smk"
include: "rules/plotting.smk"

# "rule all" is a trick to ensure that you always run your full workflow through the point of the end files.
# how does this work? 
# snakemake parses files starting from the Snakefile, then the .smk files.
# because "rule all" is the first rule it encounters, it will trigger a check for all of the files listed under "input" to see if they exist.
# if not, then the tasks to create those files are spawned.
# In the case of my work, I prefer not to use "rule all" until I have a stable build that yields some summary figures or files.
# Instead, I tell snakemake to generate specific files, as we did in the tutorial. 

# rule all:
#     input:
#         coverage_html = "output/visuals/sample_coverage.html",
#         coverage_pdf = "output/visuals/sample_coverage.pdf",
#         vcf = "output/variants/all_samples.vcf.gz",
#         vcf_heatmap = "output/visuals/vcf_heatmap.pdf"

# To create the rule graph ("simplified" graph) of this workflow:
# snakemake --rulegraph --use-conda output/visuals/vcf_heatmap.pdf output/visuals/sample_coverage.pdf | dot -Tpng > output/visuals/workflow_rulegraph.png

# To create the full graph (showing all samples) of this workflow:
# snakemake --dag --use-conda output/visuals/vcf_heatmap.pdf output/visuals/sample_coverage.pdf | dot -Tpng > output/visuals/workflow_dag.png
