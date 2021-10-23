import urllib
import GEOparse
import os
import pandas as pd
from tqdm import tqdm
import multiprocessing
import contextlib
import io
import sys

def parse_platform(platform, srr, processed_gsms):
    os.makedirs("downloads/soft", exist_ok=True)
    p = platform
    p1 =  p[0:5]+"nnn"
    p2 = p[0:9]
    url = "ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/"+p1+"/"+p2+"/soft/"+p2+"_family.soft.gz"
    urllib.request.urlretrieve(url, "downloads/soft/"+p2+".soft.gz")
    
    geo = GEOparse.get_GEO(filepath="downloads/soft/"+p2+".soft.gz", silent=True)
    
    gsmids = []
    for gsmid, gsm in geo.gsms.items():
        gsmids.append(gsmid)
    matching_gsms = list(set(srr.index).intersection(set(gsmids)))
    new_gsms = list(set(matching_gsms).difference(set(list(processed_gsms[0]))))
    sll = srr.loc[new_gsms]
    #chunk_samples = chunk(matching_gsms, 1000)
    #with multiprocessing.Pool(6) as pool:
    #    parameters = [(x, geo, sll) for x in chunk_samples]
    #    res = list(tqdm(pool.imap(check_gsm_star, parameters), desc="Scanning gsms", total=len(parameters)))
    #res = pd.concat(res)
    return check_gsm(new_gsms, geo, sll)

def chunk(l, n):
	return [l[i:i+n] for i in range(0, len(l), n)]

def check_gsm_star(args):
    return check_gsm(*args)

def check_gsm(new_gsms, geo, sll):
    sralist = []
    for gsmid in new_gsms:
        gsm = geo.gsms[gsmid]
        if gsm.metadata['library_strategy'][0] == 'RNA-Seq':
            series = gsm.metadata["series_id"][0]
            species = gsm.metadata["taxid_ch1"][0]
            sralist.append([series, gsmid, species])

    mm = pd.DataFrame(sralist)
    if mm.shape[0] == 0:
        return pd.DataFrame(columns=["gse", "gsm", "sra", "species"])
    else:
        mm.index = mm.iloc[:,1]
        mm = mm.join(sll, how="inner").iloc[:,[0,1,3,2]]
        mm.columns = ["gse", "gsm", "sra", "species"]
        return mm

def scan_platforms(srr, processed_gsm):
    platforms = ["GPL24676", "GPL24247", "GPL21626", "GPL21697", "GPL21273", "GPL20795", "GPL21493", "GPL21103", "GPL19057", "GPL18480", "GPL17021", "GPL15103", "GPL13112", "GPL21290", "GPL20301", "GPL18573", "GPL18460", "GPL16791", "GPL15433", "GPL11154", "GPL23227", "GPL23479"]
    platform_results = []
    for p in tqdm(platforms):
        res = parse_platform(p, srr, processed_gsm)
        platform_results.append(res)
    return pd.concat(platform_results)
