#RUN 16 03 2021 
cp -r /srv/nfs/disnap/NGS/NgsWeb/FastQ/210316_NB501480_0591_AHVJ7KBGXH_1616019003/SARSCOV2-ARTIC/ data_20210316/

data_loc="data_20210316/" #RAW FASTQ
result_loc="results_20210316/" # REP OUTPUT
pref="RUN_20210316" # PREFIX

# Create missing file
rm samplesheet/*
cp /srv/nfs/ngs-stockage/NGS_Virologie/PROTO_sars-cov2/artic_Illumina/samplesheet/coro20210316_samplesheet.csv samplesheet/
awk -F "," '{ print $1 }' samplesheet/*_samplesheet.csv |  sed -e "1d" > samplesheet/list_sample.txt

for sample in `cat samplesheet/list_sample.txt`;
do
    num_file=`ls ${data_loc}${sample}* | wc -l`
    #echo $num_file
    if [ $num_file -lt 2 ]
    then
        #echo "YAPA"
        touch ${data_loc}${sample}_R1.fastq.gz 
        touch ${data_loc}${sample}_R2.fastq.gz
    #else
        #echo "YA"
    fi
done

screen -S run_nf_210317
#Launch Nextflow
mkdir -p $result_loc
./nextflow run main.nf \
    -profile singularity \
    --ref ARTICDATA/nCoV-2019.reference.fasta \
    --bed ARTICDATA/nCoV-2019.bed \
    --illumina \
    --prefix $pref \
    --directory $data_loc \
    --outdir $result_loc
singularity exec --bind /data/HadrienR/ncov2019-artic-nf nextclade26.sif nextclade \
    -i ${result_loc}ncovIllumina_sequenceAnalysis_reportAllConsensus/all_consensus.fasta  \
    -t ${result_loc}nextclade_report.tsv
cp -r ${result_loc} /srv/nfs/ngs-stockage/NGS_Virologie/PROTO_sars-cov2/artic_Illumina/
chmod 777 -R /srv/nfs/ngs-stockage/NGS_Virologie/PROTO_sars-cov2/artic_Illumina/${result_loc}
