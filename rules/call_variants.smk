rule bcftools_call:
    # Takes all specimen bamfiles for variant calling.
    # The output is a single vcf file summarizing variants from all input specimens specified in config/samples.tsv.
    # Note that these variants are not likely to be super meaningful because of the read sub-sampling performed in the workshop :)
    input:
        refgenome = config['reference']['fasta'],
        bam = expand("output/bwa/{id}.sorted.bam", id = specimens),
        index = expand("output/bwa/{id}.sorted.bam.bai", id = specimens)
    output:
        "output/bcftools/all_samples.vcf"
    conda:
        "../envs/call_variants.yml"
    threads: 10
    shell:
        """
        bcftools mpileup -f {input.refgenome} {input.bam} | bcftools call --ploidy 1 -mv - > {output}
        """