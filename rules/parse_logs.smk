import os
from collections import defaultdict
import pandas as pd
import gspread
from oauth2client.service_account import ServiceAccountCredentials 

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

rule get_case_data:
    output:
        clin_data = "analysis/{imgt_release}/clin_data.csv"
    params:
        auth_path = config["auth"]
    run:
        scope = ['https://spreadsheets.google.com/feeds','https://www.googleapis.com/auth/drive']   
        creds = ServiceAccountCredentials.from_json_keyfile_name(params.auth_path, scope)  
        client = gspread.authorize(creds) 
        sheet = client.open('Copy of clean_all_clinical_and_seq').sheet1  

        df = pd.DataFrame(sheet.get_all_records())
       #df.set_index("dog_id", drop=True, inplace=True)
        df.to_csv(output.clin_data, sep=",", index=False)

rule parse_mitools_fastqs:
    input:
        merged = expand(
                "pipeline/{u.sample}/{u.time}/mitools/MIG{size}/{chain}_merged.t1.cf.fastq",
                u=units.itertuples(), size=config["mig_size"],
                chain=config["chains"]
        ),
        clin_csv = f"analysis/{config['imgt_release']}/clin_data.csv"
    output:
        umi_out = f".logs/{config['imgt_release']}/sample_umi_reads.csv"
    run:
        # read in clincial data and convert to dict
        clin_data = pd.read_csv(input.clin_csv)
        clin_data.set_index("dog_id", drop=True, inplace=True)
        clin_data.to_dict(orient="index")
        clin_dict = clin_data.to_dict(orient="index")
        
        d = {}
        for i in input.merged:
            tmp = i.split("/")
            # get dog, time, and chain from merged fastq file    
            dog,time,chain = tmp[-5],tmp[-4],tmp[-1].split("_")[0]
            key = f"{dog}_{time}_{chain}"
            if key not in d:
                d[key] = {}

            with open(i) as infile:
                for line in infile:
                    if line.startswith("@MIG"):
                        # get number of reads/umi and use as key
                        umi_count = line.strip().split(":")[-1]
                        if umi_count not in d[key]:
                            d[key][umi_count] = 1
                        # increase count if already seen
                        else:
                            d[key][umi_count] += 1

        with open(output.umi_out,"w") as outfile:
            print("dog_id,case_status,time,type,umi,count", file=outfile)
            for i in d:
                for k,v in sorted(d[i].items(), key=lambda item: int(item[0])):
                    dog,time,chain = i.split("_")
                    print(
                        dog,clin_dict[dog]["status"],
                        time,chain,k,v,
                        sep = ",",file=outfile
                    )

rule parse_migec_logs:
    input:
        assem = expand(
                "pipeline/{u.sample}/{u.time}/migec/Assemble/MIG{size}/assemble.cmd.txt",
                u=units.itertuples(), size=config["mig_size"]
        ),
        clin_csv = "analysis/{imgt_release}/clin_data.csv"
    output:
        logs_out = ".logs/{imgt_release}/all_dogs.migec_logs.csv"
    run:
        # read in clincial data and convert to dict
        clin_data = pd.read_csv(input.clin_csv)
        clin_data.set_index("dog_id", drop=True, inplace=True)
        clin_data.to_dict(orient="index")
        clin_dict = clin_data.to_dict(orient="index")

        checkout_logs = gen_migec_logs(
            input.assem,
            "checkout.log.txt",
            config["mig_size"]
        )

        assemble_logs = gen_migec_logs(
            input.assem,
            "assemble.log.txt",
            config["mig_size"]
        )

        chains = ["IGM","IGG","IGA","IGE","IGH","IGK","IGL"]
        d_migec = defaultdict(list)

        for l in checkout_logs:
            # parse the dog id and time point
            check_dog_time = ".".join(l.split("/")[1:3])
            with open(l) as f:
                for line in f:
                    line = line.strip()
                    # only check chains of interest
                    if any(c in line for c in chains):
                        tmp = line.split("\t")
                        # parse chain and number of reads where primary (master) barcode was detected
                        chain_count = [tmp[2],tmp[3]]
                        d_migec[check_dog_time].append(chain_count)
        #print(d_migec['D02293.T1.23621_1#1'])
        
        for l in assemble_logs:
            # parse the dog id and time point
            assem_dog_time = ".".join(l.split("/")[1:3])
            with open(l) as f:
                for line in f:
                    line = line.strip()
                    if any(c in line for c in chains):
                        tmp = line.split("\t")
                        # migs or reads dropped due to collisions (erroneous (1-mismatch) 
                        # variant of some UMI with higher count) by subtracting TOTAL 
                        # from GOOD_TOTAL. migs or reads dropped due overseq (MIG size
                        # threshold) no collected since MIG size set 1, nothing will be dropped
                        assem = [
                            tmp[0], # chain
                            tmp[9], # MIGS_GOOD_TOTAL number of succesfully assembled consensuses that have both R1 and R2 parts
                            tmp[10], # MIGS_TOTAL total number of input UMIs prior to coverage filtering
                            tmp[13], # READS_GOOD_TOTAL number of paired reads in succesfully assembled consensuses that have both R1 and R2 parts
                            tmp[14], # READS_TOTAL total number of input reads prior to coverage filtering
                            tmp[15], # READS_DROPPED_WITHIN_MIG_1 number of reads dropped during consensus assembly - high number of mismatches to the consensus in R1
                            tmp[16] # READS_DROPPED_WITHIN_MIG_2
                        ]
                        if assem_dog_time in d_migec:
                            d_migec[assem_dog_time].append(assem)
                        else:
                            print(f"WARNING - check {assem_dog_time}")
        #print(d_migec['D02293.T1.23621_1#1'])
        
        with open(output.logs_out, "w") as migec_out:
            print(
                "dog_id","case_status","time","type",
                "master_barcode_reads",
                "migs_good_total",
                "migs_total",
                "reads_good_total",
                "reads_total",
                "reads_drop_in_mig_r1",
                "reads_drop_in_mig_r2",
                sep=",",file=migec_out
            )
            for k,v in d_migec.items():
                dog,time = k.split(".")
                for i in range(6):
                # combine checkout and assemble metrics by chain
                    ig_data = d_migec[k][i] + d_migec[k][i+6][1:]
                    print(
                        dog,clin_dict[dog]["status"],time,
                        ",".join(ig_data),
                        sep=",",file=migec_out
                    )

        
rule parse_mixcr_logs:
    input:
        #paths_in = utils.pipeline_out(config["samples"],config["mig_size"])
        clones = expand(
                "pipeline/{u.sample}/{u.time}/mixcr/MIG{size}/{imgt_release}/{chain}/clones_{chain}.txt",
                u=units.itertuples(), size=config["mig_size"],
                chain=config["chains"], imgt_release=config["imgt_release"]
        ),
        clin_csv = "analysis/{imgt_release}/clin_data.csv"
    output:
        logs_out = ".logs/{imgt_release}/all_dogs.mixcr_logs.csv"
    #params:
        #status = utils.case_lookup(config["auth"])
    run:
        # read in clincial data and convert to dict
        clin_data = pd.read_csv(input.clin_csv)
        clin_data.set_index("dog_id", drop=True, inplace=True)
        clin_dict = clin_data.to_dict(orient="index")

        mixcr_align = [
            os.path.join(
                os.path.dirname(i),
                f"log_align_{i.split('/')[-2]}.txt"
            ) 
            for i in input.clones
        ]

        mixcr_assem = [
            os.path.join(
                os.path.dirname(i),
                f"log_assemble_{i.split('/')[-2]}.txt"
            )
            for i in input.clones
        ]

        d_mixcr = defaultdict(list)

        for l in mixcr_align:
            # parse the dog id, time point, and seq run name
            check_dog_time = ".".join(l.split("/")[1:3])
            # get chain and add to dict
            chain = os.path.basename(l).split("_")[-1].split(".")[0]
            key = f"{check_dog_time}.{chain}"
            d_mixcr[key].append(chain)
            with open(l) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("Total sequencing reads"):
                        d_mixcr[key].append(line.split(" ")[-1])
                    elif line.startswith("Successfully aligned reads"):
                        # get count and parse percent
                        d_mixcr[key].append(line.split(" ")[-2])
                        d_mixcr[key].append(get_percent(line))
                   #elif line.startswith("Alignment failed, no hits"):
                   #    d_mixcr[key].append(line.split(" ")[-2])
                   #elif line.startswith("Alignment failed because of absence of J hits"):
                   #    d_mixcr[key].append(line.split(" ")[-2])

        for l in mixcr_assem:
            # parse the dog id, time point, and seq run name
            check_dog_time = ".".join(l.split("/")[1:3])
            # get chain and add to dict
            chain = os.path.basename(l).split("_")[-1].split(".")[0]
            key = f"{check_dog_time}.{chain}"
            d_mixcr[key].append(chain)
            with open(l) as f:
                for line in f:
                    line = line.strip()
                    if line.startswith("Final clonotype count"):
                        d_mixcr[key].append(line.split(" ")[-1])
                    elif line.startswith("Average number of reads per clonotype"):
                        d_mixcr[key].append(line.split(" ")[-1])
                    elif line.startswith("Reads used in clonotypes, percent of total"):
                        # get count and parse percent
                        d_mixcr[key].append(line.split(" ")[-2])
                        d_mixcr[key].append(get_percent(line))
                    elif line.startswith("Reads clustered in PCR error correction"):
                        d_mixcr[key].append(line.split(" ")[-2])

        with open(output.logs_out, "w") as mixcr_out:
            print(
                "dog_id","case_status","time","type",
                "reads_analyzed", # total number of analysed sequencing
                "aligned_reads", # number of successful alignments
                "aligned_reads_percent", # reads_analyzed/aligned_reads
               #"align_fail_no_hits", # reads that did not align against imgt
               #"align_fail_no_j_hits", # reads without a valid j hit
                "number_of_clonotypes", # number of clonotypes after all error correction steps
                "avg_reads_per_clone", # sum_of_clontypes/number_of_clonotypes
                "total_clones", # sum of all clonotype abundances
                "reads_used_in_clone_percent", # sum_of_clonotypes/reads_analyzed (align step)
                "reads_clustered", # number of reads in clonotypes that were clustered
                sep=",",file=mixcr_out
            )
            for k,v in d_mixcr.items():
                dog,time,chain = k.split(".")
                print(
                    dog,clin_dict[dog]["status"],time,chain,
                    ",".join(v[1:4] + v[5:]),
                    sep=",",file=mixcr_out
                )

rule combine_logs:
    input:
        migec   = ".logs/{imgt_release}/all_dogs.migec_logs.csv",
        mixcr   = ".logs/{imgt_release}/all_dogs.mixcr_logs.csv",
        umi_out = ".logs/{imgt_release}/sample_umi_reads.csv"
    output:
        merge = ".logs/{imgt_release}/all_dogs.all_logs.csv"
    run:
        df_migec = pd.read_csv(input.migec)
        df_mixcr = pd.read_csv(input.mixcr)

        merged = df_migec.merge(
            df_mixcr, 
            on=["dog_id","case_status","time","type"],
            how="right"
        )

        merged.to_csv(output.merge,index=False,header=True)
