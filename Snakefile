import pandas as pd
SRA_df=pd.read_csv("SRA_accessions_filtered.csv")
sample_list=list(SRA_df['Run'])

rule all:
    input:
        expand("fastq_files/{sample}_1.fastq", sample = sample_list),
        expand("fastq_files/{sample}_2.fastq", sample = sample_list)

rule download_SRA:
    output:
        "fastq_files/{sample}_1.fastq",
        "fastq_files/{sample}_2.fastq"
    shell:
        '''fastq-dump --split-files {wildcards.sample} -outdir fastq_files'''
