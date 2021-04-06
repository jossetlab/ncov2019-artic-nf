library(seqinr)
library(stringr)
#library(DECIPHER)

#setwd("Z:/PROTO_sars-cov2/artic_Illumina")
#rep_res<-"results_20210212/"
#samplesheet_name<-"coro20210219_samplesheet.csv"

argv <- commandArgs(TRUE)
rep_res<-as.character(argv[1])
#samplesheet_path<-as.character(argv[2])

cov_table<-read.table(paste0(rep_res,"all_cov.cov"),header=F)
colnames(cov_table)<-c("reference","position","coverage","sample")

count_table<-read.table(paste0(rep_res,"all_counts.tsv"),header=F)
colnames(count_table)<-c("reference","type","count","sample")
var_1020<-subset(count_table, type == "10-20%", select=c(sample, count))
colnames(var_1020)<-c("sample","10-20%")
var_2050<-subset(count_table, type == "20-50%", select=c(sample, count))
colnames(var_2050)<-c("sample","20-50%")

mean_cov<-aggregate(cov_table$coverage,by=list(cov_table$sample),mean)
colnames(mean_cov)<-c("sample","mean_coverage")
mean_cov$mean_coverage<-format(mean_cov$mean_coverage,scientific = FALSE)
mean_cov$mean_coverage<-as.numeric(as.character(mean_cov$mean_coverage))

## FASTA
alnfa <- read.fasta(file = paste0(rep_res,"all_consensus.fasta"))
seqnames<-names(alnfa)
lengthseq = unlist(lapply(alnfa, function(l) length(l)))
length_N_tot = unlist(lapply(alnfa, function(l) length(grep("n",l))))
alnfa_trimN = unlist(lapply(alnfa, function(l) trimSpace(c2s(l), leading = TRUE, trailing = TRUE, space ="[n]{1,}")))
length_N_int = unlist(lapply(alnfa_trimN, function(l) length(grep("n",s2c(l)))))

trimcons <- read.fasta(file = paste0(rep_res,"all_trimcons.fasta"))
seqnames<-names(trimcons)
Nseqtrim = unlist(lapply(trimcons, function(l) length(grep("n",l))))

N_check_results<-as.data.frame(cbind(seqnames,length_N_int,length_N_tot,lengthseq,Nseqtrim))
N_check_results$percCOV<-(as.numeric(as.character(N_check_results$lengthseq))-as.numeric(as.character(N_check_results$length_N_tot)))/29903*100

posc_table<-read.table(paste0(rep_res,"posc.tsv"),header=T)
posc_table$SAMPLEID <- NULL
posc <- data.frame(sample = character(0), hasposc = character(0), stringsAsFactors=FALSE)
for (i in colnames(posc_table)) {posc[nrow(posc)+1,]<-c(gsub("\\.", "-", i), any(posc_table[,i] >= 0.90))}
posc$hasposc<-gsub("FALSE", "FAILED", posc$hasposc)
posc$hasposc<-gsub("TRUE", "OK", posc$hasposc)

conta_table<-read.table(paste0(rep_res,"contamination_common_poolt.tsv"),header=T)
conta_table$SAMPLEID <- NULL
conta <- data.frame(sample = character(0), hasdp = character(0), stringsAsFactors=FALSE)
for (i in colnames(conta_table)) {conta[nrow(conta)+1,]<-c(gsub("\\.", "-", i), any(conta_table[,i] >= 5))}
conta$hasdp<-gsub("FALSE", "NO", conta$hasdp)
conta$hasdp<-gsub("TRUE", "HASDP", conta$hasdp)

alldata<-merge(mean_cov,N_check_results,by.x = "sample",by.y = "seqnames", all.x=TRUE)
alldata<-merge(alldata,var_1020,by = "sample", all.x=TRUE)
alldata<-merge(alldata,var_2050,by = "sample", all.x=TRUE)
alldata<-merge(alldata,posc,by = "sample", all.x=TRUE)
alldata<-merge(alldata,conta,by = "sample", all.x=TRUE)

#samplesheet<-read.table(paste0(samplesheet_path),header=T,sep=",")

#getsamplename<-function(fullname){
#  split<-str_split(fullname,"_")
#  sample<-as.character(split[[1]][[1]])
#  return(sample)
#}

#sample_procceed<-sapply(alldata$sample,getsamplename)
#sample_procceed<-unique(sample_procceed)

#not_here<-subset(samplesheet,!(Sample_ID %in% sample_procceed))

#if(nrow(not_here)>0){
#  missing_sample<-as.character(not_here$Sample_ID)
#  missing_sample<-cbind(missing_sample,0,0,0,0,0,0)
#  colnames(missing_sample)<-colnames(alldata)
#  alldata<-rbind(alldata,missing_sample)
#}

write.csv2(alldata,paste0(rep_res,"summary.csv"),quote=F,row.names = F)
