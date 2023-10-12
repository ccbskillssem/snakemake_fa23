rule plotly_coverage:
    # Creates an interactive plot summary of coverage per sample.
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