# Notes on rules here:

# We use wildcards to indicate a variable part of the file path. 
# Example: data/SRR2303207.fastq.gz and data/SRR2303208.fastq.gz vary by the indicated {id} component in data/{id}.fastq.gz. 
# If we specify a target output file with the path output/cutadapt/SRR2303207.trimmed.fastq.gz, snakemake understands to form-fill every instance of {id} with SRR2303207.
# You'll also see references to a variable called config in the params field of some rules: refer to the slides on config file setup to learn how/why this works.

# You'll see frequent use of the "expand" function: this is a snakemake helper utility that uses string formatting to create a list of paths with the {id} filled in, where id is a "wildcard" (placeholder).
# expand() is useful for aggregating lists of expected inputs or outputs with similar file paths.

# In all rules below, we use expand() to fill in an ID field within a path string: this is signified with {id} as a "wildcard" (placeholder).
# We specifically tell snakemake to expand {id} with the specimens variable: this variable is a list of unique specimen IDs from our sample table, and it's instantiated in the common_utils.smk file.
# By using the expand() function in each rule's input, we tell snakemake that we need the relevant files to exist for *every id* before we can create any of the summary plots: in this manner, your summary plots should *always* be up to date.
# If this isn't making complete sense, check the GitHub homepage for the "full DAG" command: this will generate a visualization of what snakemake is checking.

rule plotly_coverage:
    # Creates an ~interactive~ plot summary of coverage per sample.
    # The input is a list of paths to bam files, where snakemake expects to see one bam file per sample in the sample table.
    # See notes above for info on the expand() function.
    input: 
        bam = expand("output/bwa/{id}.sorted.bam", id = specimens),
        bai = expand("output/bwa/{id}.sorted.bam.bai", id = specimens)
    output:
        "output/visuals/sample_coverage.html"
    conda:
        "../envs/plotting.yml"
    threads: 5
    params: 
        format = "plotly",
        title = "'Example plot describing sequencing coverage'"
    shell: 
        """
        plotCoverage -p {threads} --bamfiles {input.bam} --plotFile {output} --plotFileFormat {params.format} --plotTitle {params.title} --ignoreDuplicates --minMappingQuality 10 
        """

rule pdf_coverage:
    # Creates a PDF plot summary of coverage per sample.
    # The input is a list of paths to bam files, where snakemake expects to see one bam file per sample in the sample table.
    # See notes above for info on the expand() function.
    input: 
        bam = expand("output/bwa/{id}.sorted.bam", id = specimens),
        bai = expand("output/bwa/{id}.sorted.bam.bai", id = specimens)
    output: 
        "output/visuals/sample_coverage.pdf"
    conda:
        "../envs/plotting.yml"
    threads: 5
    params: 
        format = "pdf",
        title = "'Example plot describing sequencing coverage'"
    shell: 
        """
        plotCoverage -p {threads} --bamfiles {input.bam} --plotFile {output} --plotFileFormat {params.format} --plotTitle {params.title} --ignoreDuplicates --minMappingQuality 10 
        """

rule vcf_viewer:
    # Generates a SNP heatmap of the first 1000 variants from the vcf file.
    # This rule and script is wholly derived from the BAGEP pipeline:
    # https://github.com/idolawoye/BAGEP
    # Many thanks to the BAGEP team for their work!
    input:
        "output/variants/all_samples.vcf.gz"
    output:
        "output/visuals/vcf_heatmap.pdf"
    conda:
        "../envs/plotting.yml"
    params:
        outdir = "output/visuals"
    shell:
        """
        Rscript scripts/vcf2heatmap.R {input} {output}
        """