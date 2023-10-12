# Notes on these rules:
# We use wildcards to indicate a variable part of the file path. 
# Example: data/SRR2303207.fastq.gz and data/SRR2303208.fastq.gz vary by the indicated {id} component in data/{id}.fastq.gz. 
# If we specify a target output file with the path output/cutadapt/SRR2303207.trimmed.fastq.gz, snakemake understands to form-fill every instance of {id} with SRR2303207.
# You'll also see references to a variable called config in the params field of some rules: refer to the slides on config file setup to learn how/why this works.

rule get_reads:
    # Given an NCBI Sequence Read Archive (SRA) ID, downloads a subset of reads associated with the ID.
    # This rule doesn't operate on a file as input, so we don't provide an input field.
    # Note: We apply the temporary flag to these files because we want to combine them into a single interleaved (read 1 + read 2) file (see rules/preprocssing.smk).
    output:
        temp("data/{id}_1.fastq"),
        temp("data/{id}_2.fastq")
    conda:
        '../envs/get_data.yml'
    params:
        outdir = 'data',
        n_reads = config['fastq-dump']['n_reads']
    shell:
        """
        fastq-dump --split-files -X {params.n_reads} -O {params.outdir} {wildcards.id} 
        """

rule get_refgenome:
    # Downloads the relevant reference genome.
    # This rule doesn't operate on a file as input, so we don't provide an input field. 
    # Instead, we provide the FTP address of the file in the parameters field, where it will be (correctly) interpreted as a string.
    output:
        config['reference']['fasta']
    conda:
        '../envs/get_data.yml'
    params:
        ftp = 'https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/009/496/975/GCF_009496975.1_ASM949697v1/GCF_009496975.1_ASM949697v1_genomic.fna.gz',
        outdir = 'data'
    shell:
        """
        wget -P {params.outdir} {params.ftp}
        gunzip {output}.gz
        """