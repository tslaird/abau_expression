import pandas as pd
import urllib.request
import hashlib

SRA_df=pd.read_csv("SRA_accessions_filtered.csv")
sample_list=list(SRA_df['Run'])

rule all:
    input:
        expand("infofiles/{sample}.info", sample = sample_list),
        expand("fastq_files/{sample}_1.fastq.gz", sample = sample_list),
        expand("fastq_files/{sample}_2.fastq.gz", sample = sample_list),
        expand("fastq_files/{sample}_1.fastq",sample = sample_list),
        expand("fastq_files/{sample}_2.fastq.gz",sample = sample_list),
        expand("cleaned_fastq/{sample}_1.clean.fastq",sample =sample_list),
        expand("cleaned_fastq/{sample}_2.clean.fastq",sample =sample_list),
        "genome/GCF_001593425.2_ASM159342v2_genomic.gbff.gz",
        "genome/GCF_001593425.2_ASM159342v2_genomic.gbff"
        directory("salmon_transcriptome_index"),
        "genome/GCF_001593425.2_ASM159342v2_genomic_transcripts.fasta"


rule download_infofiles:
    output:"infofiles/{sample}.info"
    run:
        url= "https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&accession="+wildcards.sample
        urllib.request.urlretrieve(url,"infofiles/"+wildcards.sample+".info")

rule download_fastq:
    input:"infofiles/{sample}.info"
    output:
        "fastq_files/{sample}_1.fastq.gz",
        "fastq_files/{sample}_2.fastq.gz"
    run:
        info_df=pd.read_table(str(input), sep='\t')
        ftp_links=info_df['fastq_ftp'][0].split(";")
        md5_sums=info_df['fastq_md5'][0].split(";")
        for i in [0,1]:
            path_fastq = "http://"+ftp_links[i]
            name = 'fastq_files/'+path_fastq.split('/')[-1]
            md5= md5_sums[i]
            retries=20
            while(retries > 0):
                try:
                    urllib.request.urlretrieve(path_fastq,name)
                    with open(name, 'rb') as file_to_check:
                        data = file_to_check.read()
                        md5_returned = str(hashlib.md5(data).hexdigest())
                    if md5_returned==md5:
                        print("Fetched " + name)
                        break
                except:
                    print("Retrying download from " + path_fastq)
                    retries = retries - 1
                    continue

rule uncompress_fastq:
    input:
        read1="fastq_files/{sample}_1.fastq.gz",
        read2="fastq_files/{sample}_2.fastq.gz"
    output:
        "fastq_files/{sample}_1.fastq",
        "fastq_files/{sample}_2.fastq"
    shell:
        """gunzip -k {input.read1}
        gunzip -k {input.read2}"""

rule bbduk_trim:
    input:
        read1="fastq_files/{sample}_1.fastq",
        read2="fastq_files/{sample}_2.fastq",
    output:
        read1="cleaned_fastq/{sample}_1.clean.fastq",
        read2="cleaned_fastq/{sample}_2.clean.fastq",
    shell: """bbduk.sh in1={input.read1} in2={input.read2} out1={output.read1} out2={output.read2} ref=~/miniconda3/envs/rna_seq/opt/bbmap-38.79-0/resources/adapters.fa ktrim=r k=23 mink=11 hdist=1 trimq=20 -Xmx10g"""

rule download_genome:
    output:
        zipped="genome/GCF_001593425.2_ASM159342v2_genomic.gbff.gz"
        unzipped="genome/GCF_001593425.2_ASM159342v2_genomic.gbff"
    shell:
        '''
        wget -O {output.zipped} https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/001/593/425/GCF_001593425.2_ASM159342v2/GCF_001593425.2_ASM159342v2_genomic.gbff.gz
        gunzip -k {output.zipped}
        '''

rule generate_transcriptome:
    input: "genome/GCF_001593425.2_ASM159342v2_genomic.gbff"
    output: "genome/GCF_001593425.2_ASM159342v2_genomic_transcripts.fasta"
    shell:
        '''./src/parse_gbff_transcripts.py {input}'''

rule index_transcriptome:
    input: "genome/GCF_001593425.2_ASM159342v2_genomic_transcripts.fasta"
    output: directory("salmon_transcriptome_index")
    shell:
        '''salmon index -t {input} -i {output} -k 31'''
