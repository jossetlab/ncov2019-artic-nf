#RUN 12 03 2021 
cp -r /srv/nfs/ngs-stockage/NGS_commun/disnap/NgsWeb/FastQ/210312_NB552333_0085_AHVJ7LBGXH_1615670402/SARSCOV2-ARTIC/ data_20210312/

data_loc="data_20210312/" #RAW FASTQ
result_loc="results_20210312/" # REP OUTPUT
#Launch Nextflow
mkdir -p $result_loc
./nextflow run main.nf \
    -profile singularity \
    -resume \
    --ref ARTICDATA/nCoV-2019.reference.fasta \
    --bed ARTICDATA/nCoV-2019.bed \
    --illumina \
    --prefix "RUN_20210312" \
    --directory $data_loc \
    --outdir $result_loc

#########################################################################################################

#test Bruno:
cp -r /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/Pipeline_ncov_Bruno/ARTICDATA/ ./
cp -r /srv/nfs/ngs-stockage/NGS_Virologie/HadrienR/Pipeline_ncov_Bruno/ARTICDATA/ ./
cp /srv/nfs/ngs-stockage/NGS_Virologie/PROTO_sars-cov2/artic_Illumina/samplesheet/coro20210310_samplesheet.csv samplesheet/
curl -s https://get.nextflow.io | bash
cp /srv/scratch/HadrienRegue/ncov2019-artic-nf/nextflow /srv/scratch/HadrienRegue/BS/ncov2019-artic-nf/

export http_proxy=http://ge91097:8888/
export https_proxy=http://ge91097:8888/
data_loc="data_20210310/" #RAW FASTQ
result_loc="results_20210310/" # REP OUTPUT
#Launch Nextflow
mkdir -p $result_loc
./nextflow run ../ncov2019-artic-nf \
    -profile singularity \
    --illumina \
    --ref ARTICDATA/nCoV-2019.reference.fasta \
    --bed ARTICDATA/nCoV-2019.bed \
    --prefix "RUN20210310" \
    --directory $data_loc \
    --outdir $result_loc

./nextflow run ../ncov2019-artic-nf \
    -profile singularity \
    --illumina \
    -resume \
    --ref ARTICDATA/nCoV-2019.reference.fasta \
    --bed ARTICDATA/nCoV-2019.bed \
    --prefix "RUN20210310" \
    --directory $data_loc \
    --outdir $result_loc



awk -F "," '{ print $1 }' samplefile.txt |  sed -e "1d" > list_sample.txt

for sample in `cat list_sample.txt`;
do
    num_file=`ls data_20210301/${sample}* | wc -l`
    #echo $num_file
    if [ $num_file -lt 2 ]
    then
        echo "YAPA"
        touch data_20210301/${sample}_R1.fastq | gzip
        touch data_20210301/${sample}_R2.fastq | gzip
    #else
    #    echo "YA"
    fi
done




export http_proxy=http://ge91097:8888/
export https_proxy=http://ge91097:8888/
data_loc="DATASET_EVAL/" #RAW FASTQ
result_loc="SINGULARITY/" # REP OUTPUT
#Launch Nextflow
mkdir -p $result_loc
./nextflow run ../ncov2019-artic-nf \
    -profile singularity \
    --illumina \
    --prefix "LOL" \
    --cache /srv/scratch/HadrienRegue/ncov2019-artic-nf/conda \
    --directory $data_loc \
    --outdir $result_loc



#rework consensus
mkdir -p ${result_loc}CONS/
for file in `ls ${result_loc}ncovIllumina_sequenceAnalysis_makeConsensus/*.primertrimmed.consensus.fa`;
do
    withoutpath="${file##*/}"
    samplename="${withoutpath%%.*}"
    sed  "s/>.*/>${samplename}/" $file > ${result_loc}CONS/${samplename}.fasta
done
cat ${result_loc}CONS/*.fasta > ${result_loc}all_consensus.fasta
seqtk trimfq -b 150 -e 150 ${result_loc}all_consensus.fasta > ${result_loc}all_trimcons.fasta

#split fasta for consensus analysis
awk -v rep=$result_loc 'BEGIN {n_seq=0;} /^>/ {if(n_seq%50==0){file=sprintf(rep"consensus%d.fa",n_seq);} print >> file; n_seq++; next;} { print >> file; }' < ${result_loc}all_consensus.fasta


#compute coverage for eval summary
mkdir -p ${result_loc}COV/
for file in `ls ${result_loc}ncovIllumina_sequenceAnalysis_trimPrimerSequences/*primertrimmed.sorted.bam`;
do
    withoutpath="${file##*/}"
    samplename="${withoutpath%%.*}"
    bedtools genomecov -ibam $file -d -split | sed "s/$/\t${samplename}/" > ${result_loc}COV/${samplename}.cov
done
cat ${result_loc}COV/*.cov > ${result_loc}all_cov.cov

#index bam
for file in `ls ${result_loc}ncovIllumina_sequenceAnalysis_trimPrimerSequences/*`;
do
    samtools index $file
done

Rscript /srv/nfs/ngs-stockage/NGS_Virologie/PROTO_sars-cov2/artic_Illumina/AI_analysis.R $result_loc

#copy paste data
cp -r $result_loc /srv/nfs/ngs-stockage/NGS_Virologie/PROTO_sars-cov2/covidseq/
chmod -R 777 /srv/nfs/ngs-stockage/NGS_Virologie/PROTO_sars-cov2/covidseq/${result_loc}


#########################################################################################################

#PREPARATION DATASET EVAL
mkdir -p FASTQ_EVAL
#sample avec del D34
cp -r /srv/nfs/ngs-stockage/NGS_commun/disnap/NgsWeb/FastQ/201127_NB552333_0040_AHJTN3AFX2_1606580401/ViroEst-Routine/VRES20-411-53* FASTQ_EVAL/
for file in `ls FASTQ_EVAL`;
do
    mv FASTQ_EVAL/${file} FASTQ_EVAL/DEL_${file}
done
#sample PROTO artic illumina run 14012021
mkdir -p data_ARTIC_EVAL/  
for file in `cat list_sample.txt`;
do 
    cp -r /srv/nfs/ngs-stockage/NGS_commun/disnap/NgsWeb/FastQ/210114_NB552333_0053_AHKYYKAFX2_1610730602/ViroEst-Routine/${file}* data_ARTIC_EVAL/    
done
for file in `ls data_ARTIC_EVAL`;
do
    mv data_ARTIC_EVAL/${file} data_ARTIC_EVAL/ARTIC_${file}
done
#RUN 20 01 2021 covidseq
mkdir -p data_COVIDSEQ/
for file in `cat list_sample.txt`;
do 
    cp -r /srv/nfs/ngs-stockage/NGS_commun/disnap/NgsWeb/FastQ/210120_NB552333_0055_AHL27HAFX2_1611312601/ViroEst-Routine/${file}* data_COVIDSEQ/    
done
for file in `ls data_COVIDSEQ`;
do
    mv data_COVIDSEQ/${file} data_COVIDSEQ/COVIDSEQ_${file}
done
#regroup all FASTQ
mkdir -p DATASET_EVAL/
cp FASTQ_EVAL/* DATASET_EVAL/
cp data_COVIDSEQ/* DATASET_EVAL/
cp data_ARTIC_EVAL/* DATASET_EVAL/
#ls DATASET_EVAL/ > list_sample_eval.txt
#get 75 length for each reads
mkdir -p DATASET_EVAL_T75
for file in `cat list_sample_eval.txt`;
do
    echo $file;
    cutadapt --pair-filter=any -l 75 -o DATASET_EVAL_T75/${file}_T75_R1.fastq.gz -p DATASET_EVAL_T75/${file}_T75_R2.fastq.gz DATASET_EVAL/${file}_R1.fastq.gz DATASET_EVAL/${file}_R2.fastq.gz
done
gunzip DATASET_EVAL_T75/*
#only trimmed reads without subsampling
mkdir -p DATASET_EVAL_T75_NOSS/
for file in `ls DATASET_EVAL_T75/`;
do
    cp DATASET_EVAL_T75/${file} DATASET_EVAL_T75_NOSS/${file}
done
gzip DATASET_EVAL_T75_NOSS/*

#subsampling 1M
mkdir -p DATASET_EVAL_T75_SS1M/
for file in `ls DATASET_EVAL_T75/`;
do
    head -n 4000000 DATASET_EVAL_T75/${file} >  DATASET_EVAL_T75_SS1M/${file}
done
gzip DATASET_EVAL_T75_SS1M/*
#subsampling 500K
mkdir -p DATASET_EVAL_T75_SS500K/
for file in `ls DATASET_EVAL_T75/`;
do
    head -n 2000000 DATASET_EVAL_T75/${file} >  DATASET_EVAL_T75_SS500K/${file}
done
gzip DATASET_EVAL_T75_SS500K/*
#subsampling 200K
mkdir -p DATASET_EVAL_T75_SS200K/
for file in `ls DATASET_EVAL_T75/`;
do
    head -n 800000 DATASET_EVAL_T75/${file} >  DATASET_EVAL_T75_SS200K/${file}
done
gzip DATASET_EVAL_T75_SS200K/*

#get 36 length for each reads
mkdir -p DATASET_EVAL_T36
for file in `cat list_sample_eval.txt`;
do
    echo $file;
    cutadapt --pair-filter=any -l 36 -o DATASET_EVAL_T36/${file}_T36_R1.fastq.gz -p DATASET_EVAL_T36/${file}_T36_R2.fastq.gz DATASET_EVAL/${file}_R1.fastq.gz DATASET_EVAL/${file}_R2.fastq.gz
done
cp -r DATASET_EVAL_T36/ DATASET_EVAL_T36_NOSS/
gunzip DATASET_EVAL_T36/*
#subsampling 1M
mkdir -p DATASET_EVAL_T36_SS1M/
for file in `ls DATASET_EVAL_T36/`;
do
    head -n 4000000 DATASET_EVAL_T36/${file} >  DATASET_EVAL_T36_SS1M/${file}
done
gzip DATASET_EVAL_T36_SS1M/*
#subsampling 500K
mkdir -p DATASET_EVAL_T36_SS500K/
for file in `ls DATASET_EVAL_T36/`;
do
    head -n 2000000 DATASET_EVAL_T36/${file} >  DATASET_EVAL_T36_SS500K/${file}
done
gzip DATASET_EVAL_T36_SS500K/*
#subsampling 200K
mkdir -p DATASET_EVAL_T36_SS200K/
for file in `ls DATASET_EVAL_T36/`;
do
    head -n 800000 DATASET_EVAL_T36/${file} >  DATASET_EVAL_T36_SS200K/${file}
done
gzip DATASET_EVAL_T36_SS200K/*


#test minion
cp -r /srv/nfs/ngs-stockage/NGS_Virologie/CGINEVRA/Artic_Minion/2020-07-03_Run_5_ARTIC/2020-07-03_Run_5_ARTIC/20200703_1336_MN31515_FAO01103_e6ec20f7/* ONT_03072020


#commande de demultiplexage des données nanopore apres basecalling
/srv/nfs/ngs-stockage/NGS_Virologie/BINARIES/basecalling/ont-guppy-cpu/bin/guppy_barcoder \
    --num_barcode_threads 28 \
    --require_barcodes_both_ends \
    --input_path data_ONT/fastq_pass/ \
    --save_path data_ONT/barcoding/ \
    --arrangements_files "barcode_arrs_nb12.cfg barcode_arrs_nb24.cfg"

#Launch Nextflow nanopore
./nextflow run ../ncov2019-artic-nf -bg -profile conda --nanopolish --prefix "DATASET_TEST" \
    --cache /srv/scratch/HadrienRegue/ncov2019-artic-nf/conda \
    --basecalled_fastq data_ONT/barcoding/ \
    --fast5_pass data_ONT/fast5/ \
    --sequencing_summary data_ONT/sequencing_summary_FAN53258_cabfcaea.txt \
    --outdir data_ONT/results
