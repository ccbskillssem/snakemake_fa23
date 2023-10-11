rule bwa_index:
    # Creates index database for the reference genome.
    input:
        refgenome = config['reference']['fasta']
    output:
        expand(config['reference']['fasta'] + ".{suffix}", suffix=['amb', 'ann', 'bwt', 'pac', 'sa'])
    conda:
        "../envs/mapping.yml"
    shell:
        """
        bwa index {input.refgenome}
        """

rule bwa_map:
    # Performs mapping, varying memory requirements depending on the input
    # We pipe the output of bwa mem directly to samtools view, s.t. our mapped output is automatically saved in bam format.
    # Additionally, we mark the output *unsorted.bam as a temporary file that can be cleaned up automatically after it is no longer necessary in the workflow.
    input:
        refgenome = config['reference']['fasta'],
        index = expand(config['reference']['fasta'] + ".{suffix}", suffix=['amb', 'ann', 'bwt', 'pac', 'sa']),
       	filtered = "output/cutadapt/{id}.filt.fastq.gz"
    output:
        mapped=temp("output/bwa/{id}.unsorted.bam")
    conda:
        "../envs/mapping.yml"
    threads: 10
    shell:
        """
        bwa mem -t {threads} {input.refgenome} -p {input.filtered} | samtools view -bh - > {output.mapped}
        """

rule samtools_sort:
    # Sorts a mapped bam file in preparation for indexing.
    input:
        mapped="output/bwa/{id}.unsorted.bam"
    output:
        sorted="output/bwa/{id}.sorted.bam"
    conda:
        "../envs/mapping.yml"
    threads: 10
    shell:
        """
        samtools sort -@ {threads} {input.mapped} -o {output.sorted}
        """

rule samtools_index:
    # Creates .bai index files from a sorted bam file.
    input:
        mapped="output/bwa/{id}.sorted.bam"
    output:
        index="output/bwa/{id}.sorted.bam.bai"
    conda:
        "../envs/mapping.yml"
    threads: 10
    shell:
        """
        samtools index {input.mapped} -@ {threads}
        """