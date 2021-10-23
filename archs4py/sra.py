import os
import urllib.request
import pandas as pd

def download():
    os.makedirs("downloads", exist_ok=True)
    urllib.request.urlretrieve("ftp://ftp.ncbi.nlm.nih.gov/sra/reports/Metadata/SRA_Accessions.tab", "downloads/SRA_Accessions.tab")
    os.system("grep ^SRR downloads/SRA_Accessions.tab | grep GSM > downloads/gsm_sra_match.tsv")
    #os.system("grep ^SRR downloads/SRA_Accessions.tab > downloads/sra_all.tsv")

def read():
    #all_srr = pd.read_csv("downloads/sra_all.tsv", sep="\t")
    gsm_srr = pd.read_csv("downloads/gsm_sra_match.tsv", sep="\t")
    gsmtemp = gsm_srr.iloc[:,9].str.replace('_r[0-9]+', '', regex=True)
    gsmtemp = gsmtemp.str.replace('_[0-9]+', '', regex=True)
    gsmtemp = gsmtemp.str.replace('\\.[0-9]+', '', regex=True)
    srr = gsm_srr.iloc[:,0]
    srr.index = gsmtemp
    return srr

def get_gsms(overwrite=False):
    if not os.path.exists("downloads/gsm_sra_match.tsv") or overwrite:
        download()
    return read()

