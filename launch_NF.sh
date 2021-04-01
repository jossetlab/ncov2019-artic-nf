#RUN 30 03 2021 
cp -r /srv/nfs/disnap/NGS/NgsWeb/FastQ/210330_NB552333_0097_AHM7CJAFX2_1617192002/ViroEst-Routine/ data_20210330/

data_loc="data_20210330/" #RAW FASTQ
result_loc="results_20210330/" # REP OUTPUT
pref="RUN_20210330" # PREFIX

# Create missing file
rm samplesheet/*
cp /srv/nfs/ngs-stockage/NGS_Virologie/PROTO_sars-cov2/artic_Illumina/samplesheet/UF34403-2021-03-30-NS-run61_samplesheet.csv samplesheet/
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
singularity shell --bind /data/HadrienR/ncov2019-artic-nf/ artic-ncov2019-illumina.sif
mkdir -p $result_loc
./nextflow run main.nf \
    -resume \
    -profile insingularity \
    --ref ARTICDATA/nCoV-2019.reference.fasta \
    --bed ARTICDATA/nCoV-2019.bed \
    --illumina \
    --script_dir /data/HadrienR/ \
    --posc /data/HadrienR/ncov2019-artic-nf/ARTICDATA/COVIDSEQ.fasta \
    --prefix $pref \
    --directory $data_loc \
    --outdir $result_loc
exit
singularity exec --bind /data/HadrienR/ncov2019-artic-nf nextclade26.sif nextclade \
    -i ${result_loc}ncovIllumina_sequenceAnalysis_reportAllConsensus/all_consensus.fasta  \
    -t ${result_loc}ncovIllumina_sequenceAnalysis_makeSummary/nextclade_report.tsv
cp -r ${result_loc} /srv/nfs/ngs-stockage/NGS_Virologie/PROTO_sars-cov2/artic_Illumina/
chmod 777 -R /srv/nfs/ngs-stockage/NGS_Virologie/PROTO_sars-cov2/artic_Illumina/${result_loc}
