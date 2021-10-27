#!/usr/bin/python
import datetime, time
import subprocess
import shlex
import shutil
import os
import os.path
import sys
import tinys3
import glob
import urllib.request
import boto
from boto.s3.key import Key
import requests
import json
import random
from time import sleep
import time
from datetime import datetime

sleep(random.uniform(0, 5))

awsid = os.environ['AWSID']
awskey = os.environ['AWSKEY']
jp = "cloud18"

import subprocess, threading

class Command(object):
    def __init__(self, cmd):
        self.cmd = cmd
        self.process = None
    def run(self, timeout, errorfile):
        def target():
            f = open(errorfile, "w")
            self.process = subprocess.Popen(self.cmd, shell=True, stderr=f)
            self.process.communicate()
            f.close()
        thread = threading.Thread(target=target)
        thread.start()
        thread.join(timeout)
        if thread.is_alive():
            self.process.terminate()
            thread.join()
        return self.process.returncode

def upload_s3(file, key, bucket):
    conn = tinys3.Connection(awsid, awskey, tls=True)
    f = open(file,'rb')
    conn.upload(key, f, bucket)

def basename(p):
    temp = p.split("/")
    return temp[len(temp)-1]

r = requests.get("https://maayanlab.cloud/cloudalignment/givejobarchs4")
jj = r.json()

if jj['id'] != "empty":
    if str(jj['type']) == "sequencing":
        links = jj['datalinks'].split(";")
        
        try:
            for ll in links:
                ll = str(ll)
                fb = basename(ll)
                
                print(ll)
                
                if fb.endswith(".sra"):
                    print("do SRA dump...")
                    ffb = fb.split(".")[0]
                    os.mkdir("/alignment/data/uploads/"+ffb)
                    command = Command("wget https://sra-pub-run-odp.s3.amazonaws.com/sra/"+ffb+"/"+ffb+" -O /alignment/data/uploads/"+ffb+"/"+ffb)
                    command.run(timeout=15*60, errorfile="/alignment/data/results/wget.txt")
                    if os.path.isfile("/alignment/data/uploads/"+ffb+"/"+ffb):
                        command = Command('tools/sratools/fasterq-dump_2.11.3 -f --mem 2G --threads 2 --split-3 --skip-technical -O /alignment/data/uploads/'+ffb+' /alignment/data/uploads/'+ffb+'/'+ffb)
                        command.run(timeout=20*60, errorfile="/alignment/data/results/fasterq.txt")
                    else:
                        command = Command('tools/sratools/fasterq-dump_2.11.3 -f --mem 2G --threads 2 --split-3 --skip-technical -O /alignment/data/uploads/'+ffb+' '+ffb)
                        command.run(timeout=20*60, errorfile="/alignment/data/results/fasterq.txt")
                if fb.endswith(".gz"):
                    urllib.request.urlretrieve(ll, "/alignment/data/uploads/"+ffb+"/"+fb)
                    os.chdir("/alignment/data/uploads")
                    subprocess.call(shlex.split("gunzip "+fb))
                    os.chdir("/alignment")
            
            filenames = next(os.walk("/alignment/data/uploads/"+ffb))[2]
            print(filenames)
            organism = str(jj['parameters']).split(":")[1]
            index = "/alignment/data/index/"+organism+"_index.idx"
            indexlink = "https://s3.amazonaws.com/mssm-seq-index/"+organism+"_index.idx"
            
            print("Start Kallisto quantification")
            
            # load index if not already loaded
            if not os.path.isfile(index):
                print("Load index file: "+organism)
                urllib.request.urlretrieve(indexlink, index)
            
            run_status = 0
            if len(filenames) == 1:
                #with open("/alignment/data/results/runinfo.txt", "w") as f:
                #    subprocess.call(shlex.split("/alignment/tools/kallisto/kallisto quant -t 1 -i "+index+" --single -l 200 -s 20 -o /alignment/data/results /alignment/data/uploads/"+ffb+"/"+filenames[0]), stderr=f)
                command = Command("/alignment/tools/kallisto/kallisto quant -t 2 -i "+index+" --single -l 200 -s 20 -o /alignment/data/results /alignment/data/uploads/"+ffb+"/"+filenames[0])
                run_status = command.run(timeout=60*60, errorfile="/alignment/data/results/runinfo.txt")
            if len(filenames) > 1:
                #with open("/alignment/data/results/runinfo.txt", "w") as f:
                #subprocess.call(shlex.split("/alignment/tools/kallisto/kallisto quant -t 1 -i "+index+" -o /alignment/data/results /alignment/data/uploads/"+ffb+"/"+filenames[0]+" /alignment/data/uploads/"+ffb+"/"+filenames[1]), stderr=f)
                command = Command("/alignment/tools/kallisto/kallisto quant -t 2 -i "+index+" -o /alignment/data/results /alignment/data/uploads/"+ffb+"/"+filenames[0]+" /alignment/data/uploads/"+ffb+"/"+filenames[1])
                run_status = command.run(timeout=60*60, errorfile="/alignment/data/results/runinfo.txt")
            
            if run_status == 0: #alignment success
                print("Kallisto quantification completed")
                upload_s3("/alignment/data/results/abundance.tsv", str(jj['id'])+"-"+str(jj['uid'])+"_kallisto.tsv", "mssm-seq-results")
                print("Uploaded raw counts to S3")

                #mapping = "/alignment/data/mapping/"+organism+"_mapping.rda"

                #if not os.path.isfile(mapping):
                #    print("Load gene mapping information: "+organism)
                #    mappinglink = "https://s3.amazonaws.com/mssm-seq-genemapping/"+organism+"_mapping.rda"
                #    urllib.request.urlretrieve(mappinglink, mapping)
                
                #subprocess.call(shlex.split("Rscript --vanilla scripts/genelevel.r "+mapping))
                #upload_s3("/alignment/data/results/gene_abundance.tsv", str(jj['id'])+"-"+str(jj['uid'])+"_kallisto_gene.tsv", "mssm-seq-generesults")
                
                print("Uploaded raw gene counts to S3")
                
                numreads = 0
                numaligned = 0
                estimatedlength = 0
                with open('/alignment/data/results/runinfo.txt') as f:
                    lines = f.readlines()
                    for l in lines:
                        if "[quant] processed " in l:
                            sp = l.strip().split(" reads, ")
                            numreads = sp[0].replace("[quant] processed ","").replace(",","")
                            numaligned = sp[1].replace(" reads pseudoaligned","").replace(",","")
                        if "[quant] estimated average fragment length: " in l:
                            estimatedlength = l.replace("[quant] estimated average fragment length: ","").replace(",","").strip()
                
                sample = {
                    'id': str(jj['id']),
                    'uid': str(jj['uid']),
                    'nreads': int(numreads),
                    'naligned': int(numaligned),
                    'nlength': int(float(estimatedlength)),
                    'pass': jp
                }
                
                #send kallisto quantification to database
                r = requests.post('https://maayanlab.cloud/cloudalignment/finishjobarchs4', json=sample)
                print("Sent result to database")
            else: # alignment timed out
                sample = {
                    'id': str(jj['id']),
                    'uid': str(jj['uid']),
                    'nreads': -1,
                    'naligned': -1,
                    'nlength': -1,
                    'pass': jp
                }
                r = requests.post('https://maayanlab.cloud/cloudalignment/finishjobarchs4', json=sample)
        except Exception as e:
            sample = {
                'id': str(jj['id']),
                'uid': str(jj['uid']),
                'nreads': -2,
                'naligned': -2,
                'nlength': -2,
                'pass': jp
            }
            r = requests.post('https://maayanlab.cloud/cloudalignment/finishjobarchs4', json=sample)