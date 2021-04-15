import yaml
import pandas as pd
import gspread
from oauth2client.service_account import ServiceAccountCredentials 
#from snakemake.utils import validate
#from snakemake.utils import min_version

#min_version("5.7.1")

#report: "../report/workflow.rst"

###### Config file and sample sheets #####
configfile: "config.yaml"
#validate(config, schema="../schemas/config.schema.yaml")

#samples = pd.read_table(config["samples"]).set_index("sample", drop=False)
#validate(samples, schema="../schemas/samples.schema.yaml")

units = pd.read_table(config["units"], dtype=str).set_index(["sample","time"], drop=False)
units.index = units.index.set_levels([i.astype(str) for i in units.index.levels])  # enforce str in index
#validate(units, schema="../schemas/units.schema.yaml")

def get_fastq(wildcards):
    '''
    Get fastq files of given sample-unit.
    '''
    fastqs = units.loc[(wildcards.sample, wildcards.time), ["fq1", "fq2"]].dropna()
    
    return {"r1": fastqs.fq1, "r2": fastqs.fq2}

def pipeline_out(d,mig_size):
    '''
    Generate a list of all mixcr output for each seq run and chain.
    '''
    with open("config.yaml") as f:
        docs = yaml.full_load(f)
        for item, doc in docs.items():
            if "chains" in item:
                chains = doc

    fqs = []
    for time in d:
        for chain in chains:
            for key, val in d[time].items():
                for v in val:
                    fqs.append(f"MigecPipeline/{key}/{time}/{v}/Assembly/MIG{mig_size}/merged/mixcr/{chain}/clones_{chain}.txt")
    return fqs


def gen_migec_logs(l, s, mig_size):
    '''
    Generate a list of checkout or assemble logs depending on input s.
    '''
    import os
    
    logs = [
        os.path.join("/".join(i.split("/")[:4]),"Checkout",s) 
        if "checkout" in s 
        else 
            os.path.join("/".join(i.split("/")[:4]),f"Assemble/MIG{mig_size}",s) 
        for i in l
    ]
    return list(set(logs))

 
#def case_lookup(auth_path):
#    '''
#    Get dog ids and case status from google sheet, return case status for dog id input s.
#    '''
#    scope = ['https://spreadsheets.google.com/feeds','https://www.googleapis.com/auth/drive']   
#    creds = ServiceAccountCredentials.from_json_keyfile_name(auth_path, scope)  
#    client = gspread.authorize(creds) 
#    sheet = client.open('Copy of clean_all_clinical_and_seq').sheet1  
#
#    df = pd.DataFrame(sheet.get_all_records())
#    #d = dict(zip(df.dog_id,df.status))
#    
#    df.set_index("dog_id", drop=True, inplace=True)
#    return df.to_dict(orient="index")
#   #return dict(zip(df.dog_id,df.status))


def get_percent(s):
    '''
    Extract percent value from mixcr logs given log line s.
    '''
    val = s.split(" ")[-1]
    return ''.join(c for c in val if c not in '()%')

