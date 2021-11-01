#!/moto/ziab/users/jr3950/.conda/envs/samtools-110/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

# need samtools conda R to have vcftools available

indiv = args[1]
indiv_id=c(indiv)
library(stringr)

library(readr)
library(stringr)
library(dplyr)
library(data.table)
setwd(args[2])
truffle=fread(paste(args[3],"-truffle.ibd.iqr", sep=""), data.table=FALSE)
genomics_sex=fread(args[4], data.table=FALSE)

n=as.numeric(args[6])
truffle$IBD0_IQR=as.numeric(truffle$IBD0_IQR)
summary(truffle$IBD0_IQR)
truffle=subset(truffle,!is.na(IBD0_IQR))

kin_cutoff=sort(truffle$IBD0_IQR[truffle$ID1 %in% indiv_id | truffle$ID2 %in% indiv_id],decreasing=FALSE)[n]

PID=c(truffle$ID1[truffle$IBD0_IQR <= kin_cutoff & (truffle$ID1 %in% indiv_id | truffle$ID2 %in% indiv_id)],truffle$ID2[truffle$IBD0_IQR <= kin_cutoff & (truffle$ID1 %in% indiv_id | truffle$ID2 %in% indiv_id)])
PID=PID[!(PID %in% indiv_id)]




df_mendel_print=expand.grid(Family="all", Offspring=indiv_id, Father=PID[PID %in% as.matrix(genomics_sex[genomics_sex$Genomics_Sex %in% c("M","Q"), c("indv")])], Mother=PID[PID %in% as.matrix(genomics_sex[genomics_sex$Genomics_Sex %in% c("F","Q"), c("indv")])], stringsAsFactors=FALSE)

nrow(df_mendel_print)
df_mendel_print=distinct(df_mendel_print)
nrow(df_mendel_print)

df_mendel_print=subset(df_mendel_print, Father != Mother)

# indiv_1=paste("SM_",str_replace(indiv, "-","_"),sep="")
# indiv_2=paste("SM_",str_replace(indiv, "_","-"),sep="")

# df_mendel_print=subset(df_mendel_print, Offspring %in% c(indiv_1,indiv_2) | Father %in% c(indiv_1,indiv_2) | Mother %in% c(indiv_1,indiv_2))

write_tsv(x = df_mendel_print, file = paste(args[2],indiv,".ped",sep=""))
remove(truffle)
remove(genomics_sex)
gc()
#system(command=paste("date && cd ",args[2]," && vcftools --gzvcf ",args[3]," --not-chr X --not-chr Y --not-chr MT --mendel ",args[2],"",indiv,".ped -c | bgzip > ",indiv,".mendel.gz && zcat ",indiv,".mendel.gz | cut -f 5 | sort | uniq -c | sort -nr > ",indiv,".mendel.summary && rm -f ",indiv,".mendel.gz",sep=""), intern=TRUE)

if(args[5] == "nopossummary")
    {
    thecommand=paste("date && cd ",args[2]," && vcftools --gzvcf ",args[3]," --not-chr X --not-chr Y --not-chr MT --mendel ",args[2],"",indiv,".ped -c | cut -f 5 | gzip > ",indiv,".mendel.gz && cat ", args[2],indiv,".ped | cut -f 2,3,4 | sed -e 's/\t/_/g' | gzip >> ",indiv,".mendel.gz && zcat ",indiv,".mendel.gz | sort -T /moto/ziab/users/jr3950/data/genomes/tmp/ --buffer-size=6g | uniq -c | sort -nr > ",indiv,".mendel.summary && rm -f ",indiv,".mendel.gz && rm -f ", args[2],indiv,".ped",sep="")
    print(thecommand)
    system(command=thecommand, intern=TRUE)
    }else{
    # this does not fix 0 error trio disappearances
    system(command=paste("date && cd ",args[2]," && vcftools --gzvcf ",args[3]," --not-chr X --not-chr Y --not-chr MT --mendel ",args[2],"",indiv,".ped -c | cut -f 1,2,5 | gzip > ",indiv,".mendel.gz && zcat ",indiv,".mendel.gz | cut -f 3 | sort -T /moto/ziab/users/jr3950/data/genomes/tmp/ --buffer-size=6g | uniq -c | sort -nr > ",indiv,".mendel.summary &&  zcat ",indiv,".mendel.gz | cut -f 1,2 | sort -T /moto/ziab/users/jr3950/data/genomes/tmp/ --buffer-size=6g | uniq -c | sort -nr | gzip > ",indiv,".mendel.pos-summary.gz && rm -f ",indiv,".mendel.gz && rm -f ", args[2],indiv,".ped",sep=""), intern=TRUE)
    
    }

system(command="date && echo done", intern=TRUE)