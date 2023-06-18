## Download SRA files from internet
-   BioProject: PRJNA790470
-   SRA list: SRR_Acc_List.txt
The following code downloads SRA files in this BioProject  
```
prefetch --option-file SRR_Acc_List.txt --output-directory data/raw_fq
```

Parse partial SRA file to fastq  
```
fasterq-dump data/raw_fq/SRR17271990 --outdir data/raw_fq
fasterq-dump data/raw_fq/SRR17271989 --outdir data/raw_fq
```

Create sample.tsv file  
```
touch samples.tsv
vim samples.tsv
# Take 2 control samples to make trial run
# Table content
sample  time    fq1 fq2
D02293  T1  SRR17271990_1.fastq    SRR17271990_2.fastq
D02318  T1  SRR17271989_1.fastq    SRR17271989_2.fastq
```

Download required software  
```
cd ~/bin/
# migec (version:1.2.9)
wget https://github.com/mikessh/migec/releases/download/1.2.9/migec-1.2.9.zip

# mitools (version:1.5)
wget https://github.com/milaboratory/mitools/releases/download/v1.5/mitools-1.5.zip

# mixcr
wget https://github.com/milaboratory/mixcr/releases/download/v3.0.13/mixcr-3.0.13.zip
```

activate mixcr license
```
cd ~/bin/mixcr-3.0.13/
unzip 'MiXCR License_HYC.zip'
mv 'MiXCR License_HYC.zip' mi.license
```

Download imgt library
```
cd ~/bin/mixcr-3.0.13/libraries/
wget https://github.com/repseqio/library-imgt/releases/download/v8/imgt.202214-2.sv8.json.gz
gunzip imgt.202214-2.sv8.json.gz
```
