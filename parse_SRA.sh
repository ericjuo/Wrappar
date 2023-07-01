for i in $(ls -d data/raw_fq/*/)
do
fasterq-dump $i --outdir data/raw_fq/
done
