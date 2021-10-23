import urllib.request
import os
import pandas as pd
import pymysql
from sqlalchemy import create_engine
import feather
import GEOparse
from GEOparse.logger import set_verbosity
import numpy as np
import json
import hashlib

set_verbosity("INFO")

def get_data_path() -> str:
    path = os.path.join(
        os.path.dirname(__file__),
        'data/'
    )
    return(path)

def get_db_cred():
    f = open(get_data_path()+'/config.json',)
    dbcred = json.load(f)
    f.close()
    return dbcred

def connect():
    credentials = get_db_cred()
    return pymysql.connect(credentials["dbhost"], credentials["dbuser"], credentials["dbpass"], credentials["dbname"])

def get(select):
    db = connect()
    cursor = db.cursor()
    cursor.execute(select)
    data = pd.DataFrame(cursor.fetchall())
    db.close()
    return data

def insert(df, table):
    db = connect()
    credentials = get_db_cred()
    db_data = 'mysql+mysqldb://' + credentials["dbuser"] + ':' + credentials["dbpass"] + '@' + credentials["dbhost"] + ':3306/' \
           + credentials["dbname"] + '?charset=utf8mb4'
    engine = create_engine(db_data)
    df.to_sql(table, engine, if_exists='append', index=False)
    db.close()

def processed_gsms():
    return get("SELECT DISTINCT(gsm) FROM samplemapping")

def get_max_id():
    return get("SELECT MAX(listid) FROM samplemapping")[0][0]

def add_mapping(samples):
    insert(samples, "samplemapping")

def processed_sras():
    processed_sras = get("SELECT DISTINCT id, datalinks FROM sequencing")
    processed_sras = processed_sras.iloc[:,1].str.replace('^.*\\/', '', regex=True).str.replace('.sra','')
    return processed_sras

def add_sequencing_jobs(new_samples):
    pro_sras = processed_sras()
    sequencing = new_samples
    sequencing.index = sequencing.loc[:,"sra"]
    new_sra = list(set(sequencing.index).difference(set(list(pro_sras[0]))))
    sequencing = sequencing.loc[new_sra,:]
    sequencing.loc[:,"species"] = sequencing.loc[:,"species"].str.replace('10090','organism:mouse')
    sequencing.loc[:,"species"] = sequencing.loc[:,"species"].str.replace('9606','organism:human')
    sequencing.loc[:,"uid"] = list(hash(sequencing, "sra"))
    status = ["waiting"]*sequencing.shape[0]
    sequencing["status"] = status
    bucket = ["seq-results"]*sequencing.shape[0]
    sequencing["bucket"] = bucket
    sequencing = sequencing.sort_values("listid")
    sequencing = sequencing.loc[:, ["listid", "uid", "bucket", "sra", "species", "status"]]
    sequencing["datalinks"] = ["ftp://ftp-trace.ncbi.nih.gov/sra/sra-instant/reads/ByRun/sra/SRR/"+x[0:6]+"/"+x+"/"+x+".sra" for x in sequencing.loc[:,"sra"]]
    seq = sequencing.loc[:,["listid", "uid", "bucket", "datalinks", "species", "status"]]
    seq.columns = ["id", "uid", "resultbucket", "datalinks", "parameters", "status"]
    insert(seq, "sequencing")

def hash(sourcedf, column):
    return(pd.DataFrame(sourcedf.loc[:,column].values)[0].str.encode('utf-8').apply(lambda x: (hashlib.sha512(x).hexdigest().upper()[0:20])))

    
