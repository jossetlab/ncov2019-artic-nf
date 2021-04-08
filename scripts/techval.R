library(stringr)
library(reshape)
library(ggplot2)
options(stringsAsFactors = FALSE)

argv <- commandArgs(TRUE)
BILANfile<-as.character(argv[1])
NEXTCLADEfile<-as.character(argv[2])
MATRICEMUTfile<-as.character(argv[3])
VALIDATION_REPORT<-as.character(argv[4])
EXPORT_FASTFINDER<-as.character(argv[5])

####################################
##SCRIPT
##################################
read.csv2(BILANfile)->bilan
read.delim(NEXTCLADEfile)->nextclade
merge(bilan,nextclade,by.x="sample",by.y="seqName",all=TRUE)->nextstrain_valid


nextstrain_valid$sample<-sapply(nextstrain_valid$sample,function(x) str_split(x,"_")[[1]][1])
nextstrain_valid$sample<-sapply(nextstrain_valid$sample,function(x) str_split(x,"-")[[1]][2])
nextstrain_valid$sample<-sapply(nextstrain_valid$sample,function(x) ifelse(is.na(x)==TRUE,paste("NEG"),paste(x)))
#nextstrain_valid$sample<-sapply(nextstrain_valid$sample,function(x) paste("0",x,sep=""))nextstrain_valid$sample
nextstrain_valid$Target_1_cq<-sapply(nextstrain_valid$sample,function(x) if (nextstrain_valid$percCOV[nextstrain_valid$sample==x]>90)
{paste(nextstrain_valid$clade[nextstrain_valid$sample==x])}
#else if (nextclade$percCOV[nextclade$sample==x]<99)
#if(nextclade$percCOV[nextclade$sample==x]>90)
#{paste("AVIS BIO")}
else {paste("ININT")}
)

#1st check : RBM couvert ? colonne RBMcov
library(stringr)
RBM<-as.character(c(22874:23080))
nextstrain_valid$RBMcov<-sapply(nextstrain_valid$missing,function(x) ifelse(str_detect(x,paste(RBM,collapse="|")),paste("bad"),paste("good")))


#2nd check : repertorier positions RBM + heatmap
RBMlist<-c("N439","L452","S477","E484","N501")
nextstrain_valid$N439<-sapply(nextstrain_valid$aaSubstitutions,
                              function(x) ifelse(str_detect(x,"S:N439"),str_sub(grep("S:N439",str_split(x,",")[[1]],value=T),-1),paste("N")))
nextstrain_valid$L452<-sapply(nextstrain_valid$aaSubstitutions,
                              function(x) ifelse(str_detect(x,"S:L452"),str_sub(grep("S:L452",str_split(x,",")[[1]],value=T),-1),paste("L")))
nextstrain_valid$S477<-sapply(nextstrain_valid$aaSubstitutions,
                              function(x) ifelse(str_detect(x,"S:S477"),str_sub(grep("S:S477",str_split(x,",")[[1]],value=T),-1),paste("S")))
nextstrain_valid$E484<-sapply(nextstrain_valid$aaSubstitutions,
                              function(x) ifelse(str_detect(x,"S:E484"),str_sub(grep("S:E484",str_split(x,",")[[1]],value=T),-1),paste("E")))
nextstrain_valid$N501<-sapply(nextstrain_valid$aaSubstitutions,
                              function(x) ifelse(str_detect(x,"S:N501"),str_sub(grep("S:N501",str_split(x,",")[[1]],value=T),-1),paste("N")))

nextstrain_valid$N439state<-sapply(nextstrain_valid$N439,
                                   function(x) ifelse(x=="N",paste("or"),paste(x)))
nextstrain_valid$L452state<-sapply(nextstrain_valid$L452,
                                   function(x) ifelse(x=="L",paste("or"),paste(x)))
nextstrain_valid$S477state<-sapply(nextstrain_valid$S477,
                                   function(x) ifelse(x=="S",paste("or"),paste(x)))
nextstrain_valid$E484state<-sapply(nextstrain_valid$E484,
                                   function(x) ifelse(x=="E",paste("or"),paste(x)))
nextstrain_valid$N501state<-sapply(nextstrain_valid$N501,
                                   function(x) ifelse(x=="N",paste("or"),paste(x)))

nextcladeRBM<-data.frame(nextstrain_valid$sample,nextstrain_valid$Target_1_cq,nextstrain_valid$N439state,nextstrain_valid$L452state,nextstrain_valid$S477state,nextstrain_valid$E484state,nextstrain_valid$N501state)
reshapeRBM<-melt(nextcladeRBM,id=c("nextstrain_valid.Target_1_cq","nextstrain_valid.sample"))
reshapeRBM<-reshapeRBM[reshapeRBM$nextstrain_valid.Target_1_cq!="ININT",]
reshapeRBM$variable<-gsub("nextstrain_valid.","",reshapeRBM$variable)
reshapeRBM$variable<-gsub("state","",reshapeRBM$variable)

ggplot()
heatmut<-ggplot(reshapeRBM,aes(x=variable,y=nextstrain_valid.sample))+geom_tile(fill=ifelse(reshapeRBM$value=="or","white","red"))+facet_grid(nextstrain_valid.Target_1_cq~.,scales="free_y")

#3rd check : verifier positions spike et repertorier les absents

read.delim(MATRICEMUTfile,check.names = FALSE)->matricemut
read.delim(MATRICEMUTfile,check.names = FALSE)->matricemut2
matricemut<-sapply(matricemut,function(x) paste(matricemut[,1],x,sep=""))



matricemutdel<-apply(matricemut,2,function(x) ifelse(grepl("-",x)==1,x,NA))
matricemut<-apply(matricemut,2,function(x) ifelse(grepl("[A-Z]$|-",x)==1,x,NA))
matricemut<-apply(matricemut,2,function(x) ifelse(grepl("[A-Z]$",x)==1,x,NA))

matricemut3<-matricemut
matricemut3[,1]<-matricemut2[,1]
matricemut3[,2]<-as.numeric(matricemut2[,2])
matricemut3[,3]<-as.numeric(matricemut2[,3])
matricemut3[,4]<-as.numeric(matricemut2[,4])
matricemut3[,2]<-as.numeric(matricemut3[,2])
matricemut3[,3]<-as.numeric(matricemut3[,3])
matricemut3[,4]<-as.numeric(matricemut3[,4])
matricemut3[,colnames="ININT"]<-NA
matricemut<-matricemut3


clades<-colnames(matricemut)

nextstrain_valid$spikesubcheck<-sapply(nextstrain_valid$sample,function(x) ifelse(grepl(paste0("(?=.*",matricemut[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==x]][!is.na(matricemut[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==x]])],")",collapse=""),nextstrain_valid$aaSubstitutions[nextstrain_valid$sample==x],perl=T),"ok","no"))

nextstrain_valid$commentspikesubcheck<-sapply(nextstrain_valid$sample,function(y) ifelse(nextstrain_valid$spikesubcheck[nextstrain_valid$sample==y]=="ok","",ifelse(is.na(nextstrain_valid$aaSubstitutions[nextstrain_valid$sample==y]),"",paste("Mutation(s) ",paste(matricemut[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==y]][!is.na(matricemut[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==y]])][sapply(matricemut[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==y]][!is.na(matricemut[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==y]])],function(x) str_detect(nextstrain_valid$aaSubstitutions[nextstrain_valid$sample==y],x)==FALSE)],collapse=",")," non presente(s).",sep=""))))

nextstrain_valid$spikedelcheck<-sapply(nextstrain_valid$sample,function(x) ifelse(grepl(paste0("(?=.*",matricemutdel[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==x]][!is.na(matricemutdel[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==x]])],")",collapse=""),nextstrain_valid$aaDeletions[nextstrain_valid$sample==x],perl=T),"ok","no"))

nextstrain_valid$commentspikedelcheck<-sapply(nextstrain_valid$sample,function(y) ifelse(nextstrain_valid$spikedelcheck[nextstrain_valid$sample==y]=="ok","",ifelse(is.na(nextstrain_valid$aaDeletions[nextstrain_valid$sample==y]),"",paste("Deletion(s) ",paste(matricemutdel[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==y]][!is.na(matricemutdel[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==y]])][sapply(matricemutdel[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==y]][!is.na(matricemutdel[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==y]])],function(x) str_detect(nextstrain_valid$aaDeletions[nextstrain_valid$sample==y],x)==FALSE)],collapse=",")," non presente(s).",sep=""))))


#Commentaires missings sur positions clefs de la matrice

nextstrain_valid$missingsspiketot<-sapply(nextstrain_valid$sample,function(s) paste("Position(s) d'interet suivante(s) non couvertes sur la spike : ",paste(sapply(which("no"==sapply(1:nrow(matricemut),function (x) ifelse(str_detect(paste(matricemut[x,2:4]%in%eval(parse(text=paste("c(",gsub("-",":",nextstrain_valid$missing[nextstrain_valid$sample==s]),")",sep=""))),collapse=""),"TRUE"),"no","ok"))),function(y) paste(matricemut[,1][[y]])),collapse=",")," ."))


#4th check : verifier positions du RBM specifiques
#Commentaires missings par rapport a attendu


#nextstrain_valid$commentmissingsattendus<-sapply(nextstrain_valid$sample,function(s) ifelse(nextstrain_valid$spikecheck[nextstrain_valid$sample==s]=="ok","ok",
                                                                                    #paste("Mutation(s) ",paste(sapply(which("no"==sapply(1:nrow(matricemut),function (x) ifelse(str_detect(paste(matricemut[x,2:4]%in%eval(parse(text=paste("c(",gsub("-",":",nextstrain_valid$missing[nextstrain_valid$sample==s]),")",sep=""))),collapse=""),"TRUE"),"no","ok"))),function(y) paste(matricemut[,1][[y]])),collapse=",")," sur la spike non couverte(s)")))


#Commentaire mutations atypiques dans le RBD



nextstrain_valid$commentextsub<-sapply(nextstrain_valid$sample, function(a) 
  ifelse(is.character(sapply(
    which(FALSE==str_detect(
      paste(matricemut[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==a]][!is.na(matricemut[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==a]])],collapse=""),
      str_c(grep("(3[3-9][3-9]|4[0-9][0-9]|5[0-2][0-7])",gsub(",*","",str_split(nextstrain_valid$aaSubstitutions[nextstrain_valid$sample==a],"S:")[[1]][-1]),value=T)))),
    function(z)grep("(3[3-9][3-9]|4[0-9][0-9]|5[0-2][0-7])",gsub(",*","",str_split(nextstrain_valid$aaSubstitutions[nextstrain_valid$sample==a],"S:")[[1]][-1]),value=T)[z])),paste("Substitution(s) ",paste(sapply(which(FALSE==str_detect(paste(matricemut[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==a]][!is.na(matricemut[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==a]])],collapse=""),str_c(grep("(3[3-9][3-9]|4[0-9][0-9]|5[0-2][0-7])",gsub(",*","",str_split(nextstrain_valid$aaSubstitutions[nextstrain_valid$sample==a],"S:")[[1]][-1]),value=T)))),function(g)grep("(3[3-9][3-9]|4[0-9][0-9]|5[0-2][0-7])",gsub(",*","",str_split(nextstrain_valid$aaSubstitutions[nextstrain_valid$sample==a],"S:")[[1]][-1]),value=T)[g]),collapse=",")," presente(s) sur le RBD de la spike.",sep=""),""))


nextstrain_valid$commentextdel<-sapply(nextstrain_valid$sample, function(a) 
  ifelse(is.character(sapply(
    which(FALSE==str_detect(
      paste(matricemutdel[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==a]][!is.na(matricemutdel[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==a]])],collapse=""),
      str_c(grep("([1-9]|[1-9][0-9]|[1-9][0-9][0-9]|1[0-9][0-9][0-9])",gsub(",*","",str_split(nextstrain_valid$aaDeletions[nextstrain_valid$sample==a],"S:")[[1]][-1]),value=T)))),
    function(z)grep("([1-9]|[1-9][0-9]|[1-9][0-9][0-9]|1[0-9][0-9][0-9])",gsub(",*","",str_split(nextstrain_valid$aaDeletions[nextstrain_valid$sample==a],"S:")[[1]][-1]),value=T)[z])),paste("Deletion(s) ",paste(sapply(which(FALSE==str_detect(paste(matricemutdel[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==a]][!is.na(matricemutdel[,nextstrain_valid$Target_1_cq[nextstrain_valid$sample==a]])],collapse=""),str_c(grep("([1-9]|[1-9][0-9]|[1-9][0-9][0-9]|1[0-9][0-9][0-9])",gsub(",*","",str_split(nextstrain_valid$aaDeletions[nextstrain_valid$sample==a],"S:")[[1]][-1]),value=T)))),function(g)grep("([1-9]|[1-9][0-9]|[1-9][0-9][0-9]|1[0-9][0-9][0-9])",gsub(",*","",str_split(nextstrain_valid$aaDeletions[nextstrain_valid$sample==a],"S:")[[1]][-1]),value=T)[g]),collapse=",")," presente(s) sur la spike.",sep=""),""))



#5th check : interpretation finale 

#nextstrain_valid$avis<-sapply(nextstrain_valid$sample, function(a) ifelse(nextstrain_valid$spikecheck[nextstrain_valid$sample==a]=="ok"&nextstrain_valid$commentspikecheck[nextstrain_valid$sample==a]==""&nextstrain_valid$missingsspike[nextstrain_valid$sample==a]=="Mutation(s)    non couverte(s)"&nextstrain_valid$commentmissingsattendus[nextstrain_valid$sample==a]=="ok","VALIDATION ok","AVIS BIO"))
nextstrain_valid$avis<-sapply(nextstrain_valid$sample, function(a) ifelse(nextstrain_valid$spikesubcheck[nextstrain_valid$sample==a]=="ok"&nextstrain_valid$spikedelcheck[nextstrain_valid$sample==a]=="ok"&nextstrain_valid$commentspikesubcheck[nextstrain_valid$sample==a]=="","VALIDATION ok","AVIS BIO"))

#nextstrain_valid$finalcomm<-sapply(1:nrow(nextstrain_valid),function (x) paste(nextstrain_valid$commentextsub[x],nextstrain_valid$commentextdel[x],ifelse(nextstrain_valid$missingsspiketot[x]=="Position(s) d'interet suivante(s) non couvertes sur la spike :    .","",nextstrain_valid$missingsspiketot[x]),ifelse(is.na(nextstrain_valid$totalInsertions),"",paste("Positions d'insertions nucleotidiques presentes sur la sequence consensus :",nextstrain_valid$totalInsertions,sep=""),sep="")))
nextstrain_valid$finalcomm<-sapply(1:nrow(nextstrain_valid),function (x) paste(nextstrain_valid$commentextsub[x],nextstrain_valid$commentextdel[x],ifelse(nextstrain_valid$missingsspiketot[x]=="Position(s) d'interet suivante(s) non couvertes sur la spike :    .","",nextstrain_valid$missingsspiketot[x]),ifelse(nextstrain_valid$insertions[x]=="","",paste("Positions d'insertions nucleotidiques presentes sur la sequence consensus : ",nextstrain_valid$insertions[x],sep="")),sep=""))
nextstrain_valid$finalcomm<-gsub("\\.","\\.\n",nextstrain_valid$finalcomm)


###
write.csv2(nextstrain_valid, VALIDATION_REPORT)
###

nextstrain_valid$CLADE_GLIMS = nextstrain_valid$clade

nextstrain_valid$CLADE_GLIMS[nextstrain_valid$percCOV<90] = "ININT"



###################
## avec les rates (contas + pl 29)
################
#nextstrain_valid$CLADE_GLIMS[grep("CONTA",nextstrain_valid$COM.EXT)] = "CONTA"
#nextstrain_valid$CLADE_GLIMS[is.element(nextstrain_valid$Sample,pl59)] = ""

summary(as.factor(nextstrain_valid$CLADE_GLIMS))
# 19B         20A   20E (EU1) 20H/501Y.V2 20I/501Y.V1 20J/501Y.V3       CONTA       ININT 
# 4           1           1          12         201           1          94          70 


data.frame("Sample ID"=nextstrain_valid$sample, "Instrument ID(s)"="NB552333", "Analysis authorized by"="laurence.josset@chu-lyon.fr", "Assay"="SeqArt", "Target_1_cq"=nextstrain_valid$Target_1_cq)->Export

colnames(Export)=c("Sample ID", "Instrument ID(s)", "Analysis authorized by","Assay","Target_1_cq")

write.csv(Export, EXPORT_FASTFINDER, quote=FALSE, row.names = FALSE)

