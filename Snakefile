localrules: heavy_chain_merge,get_case_data,parse_mitools_fastqs,parse_migec_logs,parse_mixcr_logs,combine_logs
include: "rules/common.smk"

# DOUBLE CHECK ALL THESE ARE OKAY
id_fixes = {
        'D0220_T1':'D02220_T1',
        'D0220_T2':'D02220_T2', 
        'CAIA6773725':'D02482_IG',
        'CAIA6773726':'D02483_IG', 
        'CAIA6773723':'D02484_IG', 
        'CAIA6773720':'D02488_IG',
        'D03782':'D03728', 
        'D02118':'D02118_'
}

SAMPLE = ['D02293', 'D02318']
TIME = ['T1']
SIZE = ['1']
CHAIN = ['IGA', 'IGE', 'IGG', 'IGH', 'IGK', 'IGL', 'IGM']

# dictionary of dog_id : [status, remission, breed, gender, etc.]
#clin_data = case_lookup(config["auth"])

rule all:
    input:
       #expand(
       #    "pipeline/{u.sample}/{u.time}/mixcr/MIG{size}/{chain}/clones_{chain}.txt",
       #    u=units.itertuples(), size=config["mig_size"],
       #    chain=config["chains"]
       #),
        ############################### 
        #     sequence processing     #
        ###############################
        # generate the log output from processing fastqs to clones - extracts
        # data from migec and mixcr log files
        #f".logs/{config['imgt_release']}/all_dogs.all_logs.csv",
        expand("pipeline/{sample}/{time}/mixcr/MIG{size}/{imgt_release}/{chain}/clones_contigs_{chain}.tsv", 
                sample=SAMPLE, time=TIME, size=SIZE, imgt_release=config["imgt_release"], chain=CHAIN)
        # expand(".logs/{imgt_release}/all_dogs.all_logs.csv",
        #     imgt_release=config["imgt_release"]
        # ),
        # and mitools merged fastqs - umi reads
       #f".logs/{config['imgt_release']}/sample_umi_reads.csv",
        ###############################
        #           analyses          #
        ###############################
       #expand("analysis/{imgt_release}/filtered_chains.tsv",
       #    imgt_release=config["imgt_release"]
       #)
        # BASIC STATS - DONT NEED THE ABOVE filtered_chains AS THIS RULE REQS THAT
        # expand("analysis/{imgt_release}/{chain}/CalcBasicStats/{chain}.basicstats.txt", 
        #         chain=["IGH","IGK","IGL"], imgt_release=config["imgt_release"]),
       ## SPECTRATYPING
       #expand("analysis/{imgt_release}/{chain}/CalcSpectratype/{molec}/{chain}.spectratype.{molec}.wt.txt", 
       #    chain=["IGH","IGK","IGL"], molec=["nt","aa"],
       #    imgt_release=config["imgt_release"]),
       ## PLOT SPECTRATYPE AND VJ USAGE
       #expand("analysis/{imgt_release}/{chain}/CalcSpectratype/plots/{u.sample}_{u.time}_{chain}.fancyspectra.txt",
       #    u=units.itertuples(), chain=["IGH","IGK","IGL"],
       #    imgt_release=config["imgt_release"]),
       ## PLOT QUANTILE STATS - SHOULD CONSIDER COMBINING SOME OF THESE FOR CLEAN SAKE!
       #expand("analysis/{imgt_release}/{chain}/PlotQuantileStats/{u.sample}_{u.time}_{chain}.qstat.txt",
       #    u=units.itertuples(), chain=["IGH","IGK","IGL"],
       #    imgt_release=config["imgt_release"]),
       ## CDR3 AA PROFILE
       #expand("analysis/{imgt_release}/{chain}/CalcCdrAAProfile/Annotate/metadata.txt",
       #    chain=["IGH","IGK","IGL"], imgt_release=config["imgt_release"]),
       ## DIVERSITY STATS AND PLOTS
       #expand("analysis/{imgt_release}/{chain}/CalcDiversityStats/plots/{chain}_reads{cutoff}.{molec}.png",
       #    chain=["IGH","IGK","IGL"], molec=["aa"],
       #    cutoff=config["filtered"], imgt_release=config["imgt_release"]),
       ## IGBLAST/MAKEDB
       ## ORIGINAL DATABASE DIRECT FROM IMGT
       #expand("analysis/{imgt_release}/{chain}/MakeDb/{u.sample}_{u.time}_{chain}_db-pass.tsv",
       #    u=units.itertuples(), chain=["IGH","IGK","IGL"],
       #    size=config["mig_size"], imgt_release=config["imgt_release"])
       ## GAPPED BASED DATABASE FROM ALBERT
       #expand("analysis/{chain}/MakeDb_gapped/{u.sample}_{u.time}_{chain}_db-pass.tsv",
       #    u=units.itertuples(), chain=["IGH","IGK","IGL"],
       #    size=config["mig_size"])



include: "rules/migec.smk"
include: "rules/mitools.smk"
include: "rules/mixcr.smk"
include: "rules/parse_logs.smk"
include: "rules/vdjtools.smk"
include: "rules/changeo.smk"

