import pandas as pd
import urllib.request
import hashlib

SRA_df=pd.read_csv("SRA_accessions_filtered.csv")
sample_list=list(SRA_df['Run'])

rule all:
    input:
        expand("infofiles/{sample}.info", sample = sample_list),
        expand("fastq_files/{sample}_1.fastq.gz", sample = sample_list),
        expand("fastq_files/{sample}_2.fastq.gz", sample = sample_list)

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
