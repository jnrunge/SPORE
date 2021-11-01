#!/moto/ziab/users/jr3950/.conda/envs/my_conda2/bin/Rscript
args = commandArgs(trailingOnly=TRUE)

pureGT_file=args[1]
dir=args[2]
samples_file=args[3]
percent_of_pureGT_file=as.numeric(args[4])
prefix=args[5]

library(data.table)
library(MASS)
setwd(dir)



pureGT=fread(pureGT_file, data.table=FALSE)

samples_ahmm=read.table(samples_file, header = FALSE, stringsAsFactors = FALSE)

compare_two_ids <- function(df, ID1, ID2)
    {
    return(c(sum(df[,ID1]!="./."&df[,ID2]!="./."),
                      sum((df[,ID1]=="1/1"&df[,ID2]=="0/0")|(df[,ID2]=="1/1"&df[,ID1]=="0/0"))))
}

overlapping_loci = matrix(ncol=nrow(samples_ahmm),nrow=nrow(samples_ahmm))
mendel_homo_error = matrix(ncol=nrow(samples_ahmm),nrow=nrow(samples_ahmm))

for(i in 1:nrow(samples_ahmm))
    {
    
    for(j in 1:nrow(samples_ahmm))
        {
        
        if(i==j|j>i)
            {
            
            ID1=which(samples_ahmm$V1 == samples_ahmm$V1[i])
            ID2=which(samples_ahmm$V1 == samples_ahmm$V1[j])
            
            current <- compare_two_ids(pureGT, ID1, ID2)
    
            overlapping_loci[i,j]=current[1]/percent_of_pureGT_file
            mendel_homo_error[i,j]=current[2]/percent_of_pureGT_file
            
            overlapping_loci[j,i]=current[1]/percent_of_pureGT_file
            mendel_homo_error[j,i]=current[2]/percent_of_pureGT_file
            
            }
        
           # if(runif(1) < 0.001)
       # {print(paste(round(j/nrow(samples_ahmm)*100, 2), " % done in row.", sep=""))}
        
        }
    
    print(paste(round(i/nrow(samples_ahmm)*100, 2), " %", sep=""))
    if(runif(1) < 0.001)
        {
        print("writing matrix")
        fwrite(x=overlapping_loci, file=paste(dir,prefix,"-ahmm.new_overlap", sep=""), col.names=FALSE, row.names=FALSE, sep="\t")
        system(command = paste("cp -f ",dir,prefix,"-ahmm.new_overlap"," ", dir,prefix,"-ahmm.new_overlap.bak", sep=""), intern=TRUE)
        fwrite(x=mendel_homo_error, file=paste(dir,prefix,"-ahmm.mendel_homo_error", sep=""), col.names=FALSE, row.names=FALSE, sep="\t")
        system(command = paste("cp -f ",dir,prefix,"-ahmm.mendel_homo_error"," ", dir, prefix,"-ahmm.mendel_homo_error.bak", sep=""), intern=TRUE)
    }
    }
print("writing matrix")
fwrite(x=overlapping_loci, file=paste(dir,prefix,"-ahmm.new_overlap", sep=""), col.names=FALSE, row.names=FALSE, sep="\t")
fwrite(x=mendel_homo_error, file=paste(dir,prefix,"-ahmm.mendel_homo_error", sep=""), col.names=FALSE, row.names=FALSE, sep="\t")