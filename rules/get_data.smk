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