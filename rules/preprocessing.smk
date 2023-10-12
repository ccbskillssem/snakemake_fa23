# Notes on these rules:
# We use wildcards to indicate a variable part of the file path. 
# Example: data/SRR2303207.fastq.gz and data/SRR2303208.fastq.gz vary by the indicated {id} component in data/{id}.fastq.gz. 
# If we specify a target output file with the path output/cutadapt/SRR2303207.trimmed.fastq.gz, snakemake understands to form-fill every instance of {id} with SRR2303207.
# You'll also see references to a variable called config in the params field of some rules: refer to the slides on config file setup to learn how/why this works.

rule interleave_fastq:
    # Reformats the independent read 1 and read 2 fastqs into a single interleaved fastq using a script from bbmap.
    # After the interleaved file is created, the *_R1.fastq and *_R2.fastq files will be cleaned up.
    input:
        r1 = "data/{id}_1.fastq",
        r2 = "data/{id}_2.fastq",
    output:
        interleaved = "data/{id}.fastq.gz"
    conda:
        "../envs/preprocessing.yml"
    threads: 10
    shell:
        """
        reformat.sh in1={input.r1} in2={input.r2} out={output.interleaved}
        """

rule raw_fastqc:
    # Creates a report describing the length and quality of the raw (not trimmed or filtered) reads.
    input:
        fastq = "data/{id}.fastq.gz"
    output:
        "output/fastq/raw/{id}.html",
        "output/fastq/raw/{id}.zip"
    conda:
        "../envs/preprocessing.yml"
    threads: 5
    params:
        outdir = "output/fastqc/raw"
    shell:
        """
        fastqc {input.fastq} -d {params.outdir}
        """

rule cutadapt_trim:
    # Trims pesky sequencing adapters from reads.
    # Notice that we've marked *trimmed.fastq files as "temporary files", meaning that they will be removed once they're no longer necessary for downstream rules.
    input:
        fastq = "data/{id}.fastq.gz"
    output:
        trimmed = temp("output/cutadapt/{id}.trimmed.fastq.gz")
    conda:
        "../envs/preprocessing.yml"
    threads: 5
    params:
        overlap = config['cutadapt']['overlap']
    log:
        "logs/cutadapt/{id}.trimmed.json"
    shell:
        """
        cutadapt --interleaved \
        -A CTGTCTCTTATACACATCT \
        -G AGATGTGTATAAGAGACAG \
        -A CAAGCAGAAGACGGCATACGAGAT \
        -G AATGATACGGCGACCACCGA \
        -B TCTACACATATTCTCTGTC \
        --cores={threads} \
        -O {params.overlap} \
        -o {output.trimmed} \
        --json={log} \
        {input.fastq}
        """

rule cutadapt_filter:
    # Filters out reads that are too short after trimming.
    # "Too short" reads are retained in *.tooshort.fastq, just in case.
    # After this rule runs, *trimmed.fastq files should be automatically cleaned up.
    input:
        trimmed = "output/cutadapt/{id}.trimmed.fastq.gz"
    output:
        filtered = "output/cutadapt/{id}.filt.fastq.gz",
        tooshort = "output/cutadapt/{id}.tooshort.fastq.gz"
    conda:
        "../envs/preprocessing.yml"
    threads: 5
    params:
        min_length = config['cutadapt']['min_length']
    log:
        "logs/cutadapt/{id}.filtered.json"
    shell:
        """
        cutadapt --interleaved \
        -m {params.min_length} \
        --cores={threads} \
        -o {output.filtered} \
        --too-short-output={output.tooshort} \
        --json={log} \
        {input.trimmed}
        """

rule filt_fastqc:
    # Examines the trim-filtered reads with fastqc.
    input:
        filtered = "output/cutadapt/{id}.filt.fastq"
    output:
        "output/fastqc/filt/{id}.filt.html",
        "output/fastqc/filt/{id}.filt.zip"
    conda:
        "../envs/preprocessing.yml"
    threads: 5
    params:
        outdir = "output/fastqc/filt"
    shell:
        """
        fastqc {input.fastq} -d {params.outdir}
        """