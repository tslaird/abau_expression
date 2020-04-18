
import urllib.request
import pandas as pd



SRA_df=pd.read_csv('/home/tslaird/leveau_lab/abau_expression/SRA_accessions_filtered.csv')

ftp_fwd=[]
ftp_rev=[]
md5_fwd=[]
md5_rev=[]
for i in list(SRA_df['Run']): 
    url="https://www.ebi.ac.uk/ena/portal/api/filereport?result=read_run&accession="+i
    urllib.request.urlretrieve(url,str(i)+".info")
    info_df=pd.read_table(str(i)+".info", sep='\t')
    ftp_links=info_df['fastq_ftp'][0].split(";")
    md5_sums=info_df['fastq_md5'][0].split(";")   
    ftp_fwd.append(ftp_links[0])
    ftp_rev.append(ftp_links[1])
    md5_fwd.append(md5_sums[0])
    md5_rev.append(md5_sums[1])