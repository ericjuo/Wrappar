from collections import defaultdict
import pandas as pd
import os

rule convert_pool_filter:
    input:
       #clones = "pipeline/{sample}/{time}/mixcr/MIG{size}/{chain}/clones_IgH.txt"
        clones = "pipeline/{sample}/{time}/mixcr/MIG1/{imgt_release}/{chain}/clones_{chain}.txt"
       #expand(
       #    "pipeline/{u.sample}/{u.time}/mixcr/MIG{size}/{chain}/clones_{chain}.txt",
       #    u=units.itertuples(), size=config["mig_size"],
       #    chain="IgH"
       #),
    output:
        filtered = "analysis/{imgt_release}/{chain}/convert/aaVJpool/filter/{sample}_{time}_{chain}.vdjtc.pool.aaVJ.table.txt"
        #D02220_T2_IgH.vdjtc.pool.aaVJ.table.txt
    params:
        convert_dir = "analysis/{imgt_release}/{chain}/convert",
        pooled_dir  = "analysis/{imgt_release}/{chain}/convert/aaVJpool",
        filter_dir  = "analysis/{imgt_release}/{chain}/convert/aaVJpool/filter/",
        nf_dir      = "analysis/{imgt_release}/{chain}/convert/aaVJpool/filter/nonfunc/"
   #log:
   #    "logs/{imgt_release}/convert_pool_filter/{sample}_{time}_{chain}.log"
    threads: 4
    resources:
         time   = 30,
         mem_mb = 18000
    run:
        # get sample, time, and chain
        sample = wildcards.sample
        time = wildcards.time
        chain = wildcards.chain
        
        # if no clones identified at processing
        if not len(pd.read_csv(input.clones,sep="\t")): 
            # print to sample and time point to file
            with open(os.path.join(params.convert_dir,"no_results.err"),"a") as f:
                print(f"WARNING: no clones for {sample}_{time} [{chain}]",file=f)
        
        shell(f'''
            # log stdout and stderr
           #mkdir -p logs/{{rule}}
           #exec &> >(tee {{log}})
            
            VDJTOOLS={config["vdjtools"]}

            # convert to vdjtools format
            java -jar $VDJTOOLS Convert \
                -S mixcr \
                {{input.clones}} \
                {{params.convert_dir}}/{sample}_{time}_{chain}.vdjtc

            # pool each sample by CDR3aa AND V AND J
            java -jar $VDJTOOLS PoolSamples \
                -i aaVJ \
                {{params.convert_dir}}/{sample}_{time}_{chain}.vdjtc.clones_{chain}.txt \
                {{params.pooled_dir}}/{sample}_{time}_{chain}.vdjtc

            # filter clonotypes that contain a stop codon or frameshift 
            # in their receptor sequence
            java -jar $VDJTOOLS FilterNonFunctional \
                {{params.pooled_dir}}/{sample}_{time}_{chain}.vdjtc.pool.aaVJ.table.txt \
                {{params.filter_dir}}
            
            # retain non-functional clones
            java -jar $VDJTOOLS FilterNonFunctional \
                --negative \
                {{params.pooled_dir}}/{sample}_{time}_{chain}.vdjtc.pool.aaVJ.table.txt \
                {{params.nf_dir}}
            
            # cannot save filter summary when running in parallel...??
        ''')
        

rule vdjtools_meta_input:
    input:
       #clones_fltrd = expand("analysis/IgH/convert/aaVJpool/filter/{sample}_{time}_IgH.vdjtc.pool.filter.txt")
        clones_filtrd = expand(
            "analysis/{imgt_release}/{chain}/convert/aaVJpool/filter/{u.sample}_{u.time}_{chain}.vdjtc.pool.aaVJ.table.txt",
            u=units.itertuples(), size=config["mig_size"],
            chain=["IGH","IGK","IGL"], imgt_release=config["imgt_release"]
        ),
        proc_logs = ".logs/{imgt_release}/all_dogs.all_logs.csv",
        clin_csv  = "analysis/{imgt_release}/clin_data.csv"
    output:
        filtrd_chains = "analysis/{imgt_release}/filtered_chains.tsv",
        igh_meta      = "analysis/{imgt_release}/IGH/meta_all_IGH.txt",
        igk_meta      = "analysis/{imgt_release}/IGK/meta_all_IGK.txt",
        igl_meta      = "analysis/{imgt_release}/IGL/meta_all_IGL.txt",
        missing       = ".logs/{imgt_release}/vdjtools_meta_input/no_data.log"
    params:
        read_filtr   = config["read_filtr"],
        percent_used = config["percent_filtr"],
   #log:
   #    missing = ".logs/{imgt_release}/vdjtools_meta_input/no_data.log"
    threads: 4
    resources:
         time   = 20,
         mem_mb = 12000
    run:  
        # read in clincial data and convert to dict
        clin_data = pd.read_csv(input.clin_csv)
        clin_data.set_index("dog_id", drop=True, inplace=True)
        clin_dict = clin_data.to_dict(orient="index")
        
        # generate dictionary as sample_time : [file_path,key,sample,time,chain]
        # NOTE: chain needs to be changed to wildcards.chain to do all IGH, IGK, and IGL
        dogs = defaultdict(list)

        for i in input.clones_filtrd:
            # if no filtered data for given sample
            if not len(pd.read_csv(i,sep="\t")):
                with open(output.missing, "a") as f:
                    print(f"WARNING: no data - {i}",file=f)
                continue
            # extract dog id, time, and chain from input
            tmp = i.split("/")[-1].split(".")[0].split("_")
            dog,time,chain = tmp[0],tmp[1],tmp[2]
            key = f"{dog}_{time}_{chain}"

            # get clinical data for each sample (same data regardless of time)
            status = clin_dict[dog]["status"]
           #if not pd.isnull(status):
            if not status:
                dx = "NA"
            else:
                dx = clin_dict[dog]["diagnosis"]
    
            remis = clin_dict[dog]["remission"]
            if pd.isnull(remis):
                remis = "NA"

            if key not in dogs:
                dogs[key].extend(
                                 [os.path.abspath(i),key,dog,
                                  time,chain,status,dx,remis]
                                 )

        header = [
            "file_name","sample_id","dog_id","time",
            "chain","case_status","diagnosis","remission"
        ]

        # find dog_time_chain to filter out based on reads analyzed and percent
        # no results (can be adjusted in config.yaml)
        logs_df = pd.read_csv(input.proc_logs)
        filtr = logs_df.loc[(logs_df["type"].isin(["IGH","IGK","IGL"]))
                    & ((logs_df["reads_analyzed"] < params.read_filtr) 
                        | (logs_df["reads_used_in_clone_percent"] < params.percent_used))] \
            .loc[:,["dog_id","time","type",
                "reads_analyzed","reads_used_in_clone_percent"]]
        # write the filtered sample chains to file
        filtr.to_csv(output.filtrd_chains, index=None, sep="\t")
        # get list of columns as keys
        filtr_list = list(filtr[["dog_id", "time", "type"]].apply(lambda x: '_'.join(x), axis = 1))
        # remove key
        for key in filtr_list:
            if key in dogs:
                dogs.pop(key, None)
       
        all_df = pd.DataFrame.from_dict(dogs, orient="index", columns=header)
        
        # write each chain to csv and exclude LGL or CLL
        all_df.loc[(all_df["chain"] == "IGH") 
                & (~all_df["diagnosis"].isin(["LGL","CLL"]))] \
                .to_csv(output.igh_meta, index=None, sep="\t")
        
        all_df.loc[(all_df["chain"] == "IGK") 
                & (~all_df["diagnosis"].isin(["LGL","CLL"]))] \
                .to_csv(output.igk_meta, index=None, sep="\t")
        
        all_df.loc[(all_df["chain"] == "IGL")
                & (~all_df["diagnosis"].isin(["LGL","CLL"]))] \
                .to_csv(output.igl_meta, index=None, sep="\t")

rule basic_stats_and_vj_usage:
    input:
        meta = "analysis/{imgt_release}/{chain}/meta_all_{chain}.txt"
    output:
        "analysis/{imgt_release}/{chain}/CalcBasicStats/{chain}.basicstats.txt"
   #log:
   #    "logs/{imgt_release}/basic_stats_and_vj_usage/meta_{chain}.log"
    threads: 4
    resources:
         time   = 20,
         mem_mb = 12000
    run:
        # get chain and base
        chain = wildcards.chain
        base = os.path.dirname(input.meta)

        shell(f'''
            # log stdout and stderr
           #mkdir -p logs/{rule}
           #exec &> >(tee {log})

            VDJTOOLS='java -jar {config["vdjtools"]}'

            # split meta data by case status
            $VDJTOOLS SplitMetadata \
                -c case_status \
                {{input.meta}} {base}/SplitMetadata/CaseStatus/

            # split meta data by time 
            $VDJTOOLS SplitMetadata \
                -c time \
                 {{input.meta}} {base}/SplitMetadata/Time/

            # calculate basic stats on all samples
            $VDJTOOLS CalcBasicStats \
                -m {{input.meta}} {base}/CalcBasicStats/{chain}

            # STILL???? NOTE: NOT WORKING AS EXPECTED - NEED TO DISTINGUISH BETWEEN
            # CASES AND CONTROLS - ALSO SHOULD PLOT WITH FULL META DATA
            # cases - segment usage by time and remission
            $VDJTOOLS CalcSegmentUsage \
                -m {base}/SplitMetadata/CaseStatus/metadata.1.txt \
                -p -f time \
                {base}/CalcSegmentUsage/cases_{chain}_xTime

            $VDJTOOLS CalcSegmentUsage \
                -m {base}/SplitMetadata/CaseStatus/metadata.1.txt \
                -p -f remission \
                {base}/CalcSegmentUsage/cases_{chain}_xRemiss

            # t1 cases and controls
            $VDJTOOLS CalcSegmentUsage \
                -m {base}/SplitMetadata/Time/metadata.T1.txt \
                -p -f case_status \
                {base}/CalcSegmentUsage/t1_{chain}_xStatus
        ''')

rule calc_spectratype:
    input:
        meta = "analysis/{imgt_release}/{chain}/meta_all_{chain}.txt"
    output:
        "analysis/{imgt_release}/{chain}/CalcSpectratype/{molec}/{chain}.spectratype.{molec}.wt.txt"
    params:
        out_dir = "analysis/{imgt_release}/{chain}/CalcSpectratype/{molec}/{chain}"
    log:
        "logs/{imgt_release}/calc_spectratype/{chain}_{molec}.log"
    run:
        # use flag for amino acid
        flag = ""
        if "aa" in wildcards.molec:
            flag = "-a"

        shell(f'''
            # log stdout and stderr
            mkdir -p logs/{{rule}}
            exec &> >(tee {{log}})

            VDJTOOLS='java -jar {config["vdjtools"]}'

            # calculate spectratype for nt and aa
            $VDJTOOLS CalcSpectratype \
                -m {{input.meta}} {flag} \
                {{params.out_dir}}
        ''')

rule plot_spectra_and_vj_usage:
    input:
        filtered = "analysis/{imgt_release}/{chain}/convert/aaVJpool/filter/{sample}_{time}_{chain}.vdjtc.pool.aaVJ.table.txt"
    output:
        "analysis/{imgt_release}/{chain}/CalcSpectratype/plots/{sample}_{time}_{chain}.fancyspectra.txt"
    params:
        fancy_spec = "analysis/{imgt_release}/{chain}/CalcSpectratype/plots/{sample}_{time}_{chain}",
        spec_v     = "analysis/{imgt_release}/{chain}/PlotSpectratypeV/{sample}_{time}_{chain}",
        fancy_vj   = "analysis/{imgt_release}/{chain}/PlotFancyVJUsage/{sample}_{time}_{chain}"
    log:
        "logs/{imgt_release}/plot_spectra_and_vj_usage/{sample}_{time}_{chain}.log"
    shell:
        f'''
            # log stdout and stderr
            mkdir -p logs/{{rule}}
            exec &> >(tee {{log}})

            VDJTOOLS='java -jar {config["vdjtools"]}'
            
            $VDJTOOLS PlotFancySpectratype \
                {{input}} {{params.fancy_spec}}
            
            $VDJTOOLS PlotSpectratypeV \
                {{input}} {{params.spec_v}}
            
            # ISSUE WITH CIRCLIZE...STILL WORKING ON IT
            $VDJTOOLS PlotFancyVJUsage \
                {{input}} {{params.fancy_vj}}
        '''

rule plot_quantile_stats:
    input:
        filtered = "analysis/{imgt_release}/{chain}/convert/aaVJpool/filter/{sample}_{time}_{chain}.vdjtc.pool.aaVJ.table.txt"
    output:
        "analysis/{imgt_release}/{chain}/PlotQuantileStats/{sample}_{time}_{chain}.qstat.txt"
    params:
        quant_plot = "analysis/{imgt_release}/{chain}/PlotQuantileStats/{sample}_{time}_{chain}",
    log:
        "logs/{imgt_release}/plot_quantile_stats/{sample}_{time}_{chain}.log"
    run:
        # do not run if less than 5 rows (clonotypes) - errors
        if pd.read_csv(input.filtered, sep="\t").shape[0] > 4:
            # plot donut chart
            shell(f'''
                # log stdout and stderr
                mkdir -p logs/{{rule}}
                exec &> >(tee {{log}})

                VDJTOOLS='java -jar {config["vdjtools"]}'
                
                $VDJTOOLS PlotQuantileStats \
                    {{input}} {{params.quant_plot}}
            ''')
        else:
            # generate txt output with warning
            with open(str(output), "w") as out:
                print(
                    f"WARNING: {os.path.basename(input.filtered)} has insufficient reads", 
                    file=out
                )

rule cdr_aa_profile_and_annotate:
    input:
        meta = "analysis/{imgt_release}/{chain}/meta_all_{chain}.txt"
    output:
        "analysis/{imgt_release}/{chain}/CalcCdrAAProfile/Annotate/metadata.txt",
       # "analysis/{chain}/CalcCdrAAProfile/{chain}.cdr3aa.stat.wt.norm.txt",
       # "analysis/{chain}/CalcCdrAAProfile/{chain}.cdr3aa.stat.unwt.norm.txt"
    params:
        annotate = "analysis/{imgt_release}/{chain}/CalcCdrAAProfile/Annotate/{chain}",
        profile  = "analysis/{imgt_release}/{chain}/CalcCdrAAProfile/{chain}"
    log:
        "logs/{imgt_release}/calc_cdr_aa_profile/{chain}_aa_profile.log"
    run:
        shell(f'''
            # log stdout and stderr
            mkdir -p logs/{{rule}}
            exec &> >(tee {{log}})

            VDJTOOLS='java -jar {config["vdjtools"]}'
            
            # calculated weighted aa profile
            $VDJTOOLS CalcCdrAaStats \
                --weighted \
                --normalize \
                -m {{input.meta}} \
                {{params.profile}}
            
            # calculated unweighted aa profile
            $VDJTOOLS CalcCdrAaStats \
                --normalize \
                -m {{input.meta}} \
                {{params.profile}}
            
            # annotate at clonotype level
            $VDJTOOLS Annotate \
                -m {{input.meta}} \
                 {{params.annotate}}
        ''')

#rule annotate


rule calc_all_diversity:
    input:
        meta = "analysis/{imgt_release}/{chain}/meta_all_{chain}.txt"
    output:
        "analysis/{imgt_release}/{chain}/CalcDiversityStats/{chain}_readsAll.diversity.{molec}VJ.resampled.txt"
    params:
        diverse_dir = "analysis/{imgt_release}/{chain}/CalcDiversityStats/{chain}_readsAll"
    log:
        "logs/{imgt_release}/calc_diversity/{chain}_{molec}_readsAll.log"
    run:
        # 
        chain = wildcards.chain
        molec = wildcards.molec

        shell(f'''
            # log stdout and stderr
            mkdir -p logs/{rule}
            exec &> >(tee {log})

            VDJTOOLS='java -jar {config["vdjtools"]}'

            $VDJTOOLS CalcDiversityStats \
                -m {{input.meta}} \
                -i {molec}VJ \
                --resample-trials 1000 \
                {{params.diverse_dir}}
        ''')

rule calc_filtered_diversity:
    input:
        meta   = "analysis/{imgt_release}/{chain}/meta_all_{chain}.txt",
        resamp = "analysis/{imgt_release}/{chain}/CalcDiversityStats/{chain}_readsAll.diversity.{molec}VJ.resampled.txt"
    output:
        "analysis/{imgt_release}/{chain}/CalcDiversityStats/{chain}_reads{cutoff}.diversity.{molec}VJ.resampled.txt"
    params:
        base = "analysis/{imgt_release}/{chain}/CalcDiversityStats"
    log:
        "logs/{imgt_release}/calc_filtered_diversity/{chain}_{molec}_reads{cutoff}.log"
    run:
        chain = wildcards.chain
        molec = wildcards.molec
        cutoff = wildcards.cutoff
        
        meta_df = pd.read_csv(input.meta, delimiter="\t")
        
        # generate list of sample ids with reads greater than cutoff
        df = pd.read_csv(input.resamp, delimiter="\t")
        subset = list(df.loc[df["reads"] > int(cutoff)].sample_id)
       
        # subset the meta sheet to contain only above ids and write to csv
        filt_meta = os.path.join(params.base,f"meta_reads{cutoff}_{chain}.txt")
        meta_df.loc[meta_df["sample_id"].isin(subset)] \
            .replace(r'^\s*$', "NA", regex=True) \
            .to_csv(filt_meta, sep="\t", index=None, na_rep="NA")

        shell(f'''
            # log stdout and stderr
            mkdir -p logs/{rule}
            exec &> >(tee {log})

            VDJTOOLS='java -jar {config["vdjtools"]}'

            $VDJTOOLS CalcDiversityStats \
                -m {filt_meta} \
                -i {molec}VJ \
                --resample-trials 1000 \
                {{params.base}}/{chain}_reads{cutoff}
        ''')

rule plot_diversity:
    input:
        plot_script = "src/plot_diversity.R",
        resamp      = "analysis/{imgt_release}/{chain}/CalcDiversityStats/{chain}_reads{cutoff}.diversity.{molec}VJ.resampled.txt",
    output:
        "analysis/{imgt_release}/{chain}/CalcDiversityStats/plots/{chain}_reads{cutoff}.{molec}.png"
    log:
        "logs/{imgt_release}/plot_diversity/plot_{chain}_reads{cutoff}.{molec}.log"
    shell:
        '''
            # log stdout and stderr
            mkdir -p logs/{rule}
            exec &> >(tee {log})
            
            Rscript {input.plot_script} \
                -i {input.resamp} \
                -c {wildcards.chain} \
                -t {wildcards.cutoff} \
                -o {output}
        '''


#
# CalcPairwiseDistances
#
#`$VDJTOOLS CalcPairwiseDistances -i aaVJ -p -m meta_all_samples.pool.filter.igh.txt ./CalcPairwiseDistances/pdist`
#
### cases (exclude LGL and CLL)
#`$VDJTOOLS CalcPairwiseDistances -i aaVJ -p -m meta_hgl_cases_samples.pool.filter.igh.txt ./CalcPairwiseDistances/cases_pdist`
#
## ClusterSamples
#
#`$VDJTOOLS ClusterSamples -i aaVJ -f case_status -p ./CalcPairwiseDistances/pdist`
#
### cases (exclude LGL and CLL)
#`$VDJTOOLS ClusterSamples -i aaVJ -f time -p ./CalcPairwiseDistances/cases_pdist`
#
## CalcDiversityStats
#
#`$VDJTOOLS CalcDiversityStats -m meta_all_samples.pool.filter.igh.txt -i aaVJ --resample-trials 1000 CalcDiversityStats/all`
#`$VDJTOOLS CalcDiversityStats -m meta_all_samples.pool.filter.igh.txt -i ntVJ --resample-trials 1000 CalcDiversityStats/all`
#
#
#**cut out samples with less than 1k reads - D04081_01 D04390_01**
#`$VDJTOOLS CalcDiversityStats -m meta_reads1k_samples.pool.filter.igh.txt -i aaVJ --resample-trials 1000 CalcDiversityStats/Reads1k/exclude1k.rs1k.igh`
#
#**cut out samples with less than 5k reads**
#```
#D02292_02
#D02386_02
#D02976_02
#D02977_02
#D02158_01
#D02348_01
#D02354_01
#D02485_01
#D02487_01
#D03213_01
#D03345_01
#D04081_01
#D04390_01
#```
#
#`$VDJTOOLS CalcDiversityStats -m meta_reads5k_samples.pool.filter.igh.txt -i aaVJ --resample-trials 1000 CalcDiversityStats/Reads5k/exclude5k.rs1k.igh`
#
#
