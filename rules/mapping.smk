# Notes on these rules:
# We use wildcards to indicate a variable part of the file path. 
# Example: data/SRR2303207.fastq.gz and data/SRR2303208.fastq.gz vary by the indicated {id} component in data/{id}.fastq.gz. 
# If we specify a target output file with the path output/cutadapt/SRR2303207.trimmed.fastq.gz, snakemake understands to form-fill every instance of {id} with SRR2303207.
# # You'll also see references to a variable called config in the params field of some rules: refer to the slides on config file setup to learn how/why this works.

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