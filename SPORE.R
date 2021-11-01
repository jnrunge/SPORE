args = commandArgs(trailingOnly=TRUE)
library(this.path)
KAISER_folder=this.dir()
homozygous_mendel_script=paste(KAISER_folder,"/scripts/Mendel1.R",sep="")
Fx_Mendel_Script=paste(KAISER_folder,"/scripts/Mendel2.R",sep="")
previously_computed_three_thresholds=FALSE
BarnFixPlates=FALSE
samtools_activation=""

readIndividuals=function()
    {
    individuals=fread(paste(vcf,".samples",sep=""),data.table=FALSE, header = FALSE)
    return(individuals)
    }
getChrsInFile=function()
    {
    chromosomes=system(command=paste("zcat ",vcf," | grep -v ^# | cut -f 1 | sort | uniq",sep=""), intern=TRUE)
    return(chromosomes)
    }
SbatchTruffle=function()
    {
    truffle_chrs_command=""
    for(chr in chromosomes)
        {
        truffle_chrs_command=paste(truffle_chrs_command,'sbatch -t 0-11:59:00 --job-name=TruffleSingleChr --mem=',max_memory,' -c ',trufflecpu,' -A ',truffle_slurm_account,' --wrap "',samtools_activation,'
    cd ',folder,';
    bcftools view -r ',chr,' -Oz -o ',vcf,'.Chr',chr,'.vcf.gz ',vcf,';
    bcftools index -f ',vcf,'.Chr',chr,'.vcf.gz;
    ',trufflepath,' --vcf ',vcf,'.Chr',chr,'.vcf.gz --cpu ',trufflecpu,' --maf ',truffle_maf,' --missing ',truffle_missing,' --out ',vcf,'.Chr',chr,'.vcf.gz-truffle; gzip -f ',vcf,'.Chr',chr,'.vcf.gz-truffle.ibd";
    sleep 1; ', sep="")
    }
    truffle_chrs_command=gsub("\r?\n|\r", " ", truffle_chrs_command)
    print(truffle_chrs_command)
    current_jobs<-system(command=truffle_chrs_command, intern=TRUE)
    
    truffle_full_command=paste('sbatch -t 0-11:59:00 --job-name=TruffleAllChr --mem=',max_memory,' -c ',trufflecpu,' -A ',truffle_slurm_account,' --wrap "',samtools_activation,'
cd ',folder,';
',trufflepath,' --vcf ',vcf,' --cpu ',trufflecpu,' --maf ',truffle_maf,' --missing ',truffle_missing,' --out ',vcf,'-truffle; gzip -f ',vcf,'-truffle.ibd";', sep="")
truffle_full_command=gsub("\r?\n|\r", " ", truffle_full_command)
    print(truffle_full_command)
    
    current_job<-system(command=truffle_full_command, intern=TRUE)
    
    current_jobs=c(current_jobs[which(grepl(pattern = "Submitted", current_jobs))], current_job)
    
    for(job in current_jobs)
    {
    current_job_done=grepl("COMPLETED", system(command=paste("sacct --format State -j ",strsplit(job, " ")[[1]][4], sep=""), intern=TRUE)[3])
    while(!current_job_done)
        {
        current_job_done=grepl("COMPLETED", system(command=paste("sacct --format State -j ",strsplit(job, " ")[[1]][4], sep=""), intern=TRUE)[3])
        Sys.sleep(60)
    }
}
   
    
    }




InlineTruffle=function()
    {
    
    for(chr in chromosomes)
        {
        print(paste("Chromosome ", chr,sep=""))
        
        print("Deleting old files...")
        system(command=paste("rm -f ", vcf,'.Chr',chr,'.vcf.gz-truffle.ibd', sep=""), intern=TRUE)
        system(command=paste("rm -f ", vcf,'.Chr',chr,'.vcf.gz-truffle.ibd.gz', sep=""), intern=TRUE)
        
        print("Extract region from VCF")
        truffle_chrs_command=paste('bcftools view -r ',chr,' -Oz -o ',vcf,'.Chr',chr,'.vcf.gz ',vcf,';',sep="")
        truffle_chrs_command=gsub("\r?\n|\r", " ", truffle_chrs_command)
        system(command=truffle_chrs_command, intern=TRUE)
        
        print("Indexing...")
        truffle_chrs_command=paste('bcftools index -f ',vcf,'.Chr',chr,'.vcf.gz;',sep="")
        truffle_chrs_command=gsub("\r?\n|\r", " ", truffle_chrs_command)
        system(command=truffle_chrs_command, intern=TRUE)
        
        print("Executing TRUFFLE")
        #truffle_chrs_command=paste('--vcf ',vcf,'.Chr',chr,'.vcf.gz --cpu ',trufflecpu,' --maf ',truffle_maf,' --missing ',truffle_missing,' --out ',vcf,'.Chr',chr,'.vcf.gz-truffle',sep="")
        #truffle_chrs_command=gsub("\r?\n|\r", " ", truffle_chrs_command)
        #print(truffle_chrs_command)
        out<-exec_wait(cmd=trufflepath, args=c("--vcf",paste(vcf,'.Chr',chr,'.vcf.gz',sep=""), "--cpu", trufflecpu, "--maf", truffle_maf, "--missing", truffle_missing, "--out", paste(vcf,'.Chr',chr,'.vcf.gz-truffle',sep="")), std_out=TRUE, std_err=TRUE)
        
        system(command=paste("rm -f ", vcf,'.Chr',chr,'.vcf.gz',sep=""), intern=TRUE)
        
        #print(out$status)
        #print(as_text(out$stdout))
        #print(as_text(out$stderr))
        
        print("Compressing TRUFFLE output...")
        truffle_chrs_command=paste('gzip -f ',vcf,'.Chr',chr,'.vcf.gz-truffle.ibd; ', sep="")
        truffle_chrs_command=gsub("\r?\n|\r", " ", truffle_chrs_command)
        system(command=truffle_chrs_command, intern=TRUE)
    }
    
    
    
    
    
    system(command=paste("rm -f ", vcf,'-truffle.ibd', sep=""), intern=TRUE)
    system(command=paste("rm -f ", vcf,'-truffle.ibd.gz', sep=""), intern=TRUE)
    
    truffle_full_command=paste(trufflepath,' --vcf ',vcf,' --cpu ',trufflecpu,' --maf ',truffle_maf,' --missing ',truffle_missing,' --out ',vcf,'-truffle; gzip -f ',vcf,'-truffle.ibd', sep="")
    
    
truffle_full_command=gsub("\r?\n|\r", " ", truffle_full_command)
    
    print(system(command=truffle_full_command, intern=TRUE))
    
    print(truffle_full_command)
    
    }





SbatchGetIndividuals=function()
    {
    get_IDs_command=paste('sbatch -t 0-11:59:00 --job-name=GetIDs --mem=2G -c 1 -A ',truffle_slurm_account,' --wrap "',samtools_activation,'
cd ',folder,';
bcftools query -l ',vcf,' > ',vcf,'.samples"', sep="")
get_IDs_command=gsub("\r?\n|\r", " ", get_IDs_command)
    current_job<-system(command=get_IDs_command, intern=TRUE)
    current_job_done=grepl("COMPLETED", system(command=paste("sacct --format State -j ",strsplit(current_job, " ")[[1]][4], sep=""), intern=TRUE)[3])
    while(!current_job_done)
        {
        current_job_done=grepl("COMPLETED", system(command=paste("sacct --format State -j ",strsplit(current_job, " ")[[1]][4], sep=""), intern=TRUE)[3])
        Sys.sleep(60)
    }
    }



InlineGetIndividuals=function()
    {
    
    
    get_IDs_command=paste('bcftools query -l ',vcf,' > ',vcf,'.samples', sep="")
    get_IDs_command=gsub("\r?\n|\r", " ", get_IDs_command)
    system(command=get_IDs_command, intern=TRUE)
    
    
    
    }



SbatchGetIndividualLociCounts=function()
    {
    for (i in 1:nrow(individuals))
    {
    cmd=paste('sbatch -t 0-1:00:00 --mem=128M -c 1 -A ',truffle_slurm_account,' --job-name=CountLociImputed --wrap "',samtools_activation,' sh ',countLociPath,' \'',folder,vcf,'\' ',i+9,' \'',individuals$V1[i],'\'"', sep="")
    cmd=gsub("\r?\n|\r", " ", cmd)
    if(i == 1){current_jobs=system(command = cmd, intern=TRUE)}
    if(i > 1){current_jobs<-c(current_jobs,system(command = cmd, intern=TRUE))}
    Sys.sleep(1)
}
    for(job in current_jobs)
    {
    current_job_done=grepl("COMPLETED", system(command=paste("sacct --format State -j ",strsplit(job, " ")[[1]][4], sep=""), intern=TRUE)[3])
    while(!current_job_done)
        {
        current_job_done=grepl("COMPLETED", system(command=paste("sacct --format State -j ",strsplit(job, " ")[[1]][4], sep=""), intern=TRUE)[3])
        Sys.sleep(60)
    }
}
    }



InlineGetIndividualLociCounts=function()
    {
    
    for (i in 1:nrow(individuals))
    {
    cmd=paste('sh ',countLociPath,' \'',folder,vcf,'\' ',i+9,' \'',individuals$V1[i],'\'', sep="")
    cmd=gsub("\r?\n|\r", " ", cmd)
    system(command = cmd, intern=TRUE)
}
    
    
    }


ReadIndividualLociCounts=function()
    {
    individuals$ahmm_loci=NA
for(i in 1:nrow(individuals))
    {

    individuals$ahmm_loci[i]=fread(cmd = paste('cat ',folder,vcf,'.',individuals$V1[i],'.loci_imputed', sep=""))
    
    }
     individuals$ahmm_loci<-as.numeric(individuals$ahmm_loci)
    return(individuals)
    }




IBD0_row_IQR=function(x)
    {
    return(as.numeric(IQR(t(df[x,which(startsWith(colnames(df), "IBD0_") & colnames(df) != "IBD0_IQR"& colnames(df) != "IBD0_sd")]))))
}
IBD0_row_sd=function(x)
    {
    return(as.numeric(sd(t(df[x,which(startsWith(colnames(df), "IBD0_") & colnames(df) != "IBD0_IQR"& colnames(df) != "IBD0_sd")]))))
}
IBD1_row_IQR=function(x)
    {
    return(as.numeric(IQR(t(df[x,which(startsWith(colnames(df), "IBD1_"))]))))
}
IBD2_row_IQR=function(x)
    {
    return(as.numeric(IQR(t(df[x,which(startsWith(colnames(df), "IBD2_"))]))))
}

IncorporateIBDandIQR=function(df)
    {
    

    df$kinship=(df$IBD1/4)+(df$IBD2/2)
     
     for(chr in chromosomes)
    {
    truffle_tmp=fread(paste(folder,vcf,".Chr",chr,".vcf.gz-truffle.ibd.gz",sep=""), data.table=FALSE)
    truffle_tmp[,paste("IBD0_",chr,sep="")]=truffle_tmp$IBD0
    #truffle_tmp[,paste("IBD1_",chr,sep="")]=truffle_tmp$IBD1
    #truffle_tmp[,paste("IBD2_",chr,sep="")]=truffle_tmp$IBD2
    truffle_tmp=truffle_tmp[,c("ID1", "ID2", paste("IBD0_",chr,sep=""))]
    df=left_join(df, truffle_tmp, by = c("ID1", "ID2"))
    df=subset(df, !duplicated(paste(ID1,ID2)))     
    system(command=paste("rm -f ",folder,vcf,".Chr",chr,".vcf.gz-truffle.ibd.gz",sep=""), intern=TRUE)
}
     
    
    
    return(df)
    }

SbatchHomozygousMendel=function()
    {
    homozygous_mendel_command=paste('zcat ',vcf,' | grep -v ^# | cut -f10- > ',vcf,'.pureGT;
    cat ',vcf,'.pureGT | awk ',"'BEGIN {srand()} !/^$/ { if (rand() <= ", str_replace(as.character(downsample_for_homozygous_mendel),"0",""),") print $0}'",' > ',vcf,'.pureGT.sample;
    sbatch -t 0-119:59:00 --mem=',max_memory,' -c 1 -A ',homozygous_mendel_slurm_account,' --job-name=ManualMendelOverlappingLoci --wrap "',samtools_activation,'
    cd ',folder,';
    Rscript ',homozygous_mendel_script,' ',vcf,'.pureGT.sample ',folder,' ',vcf,'.samples ', downsample_for_homozygous_mendel,' "', vcf,'"', sep="")
    homozygous_mendel_command=gsub("\r?\n|\r", " ", homozygous_mendel_command)

    print(homozygous_mendel_command)
    
    system(command=paste("rm -f ",vcf,"-homozygous_mendel_command.sh",sep=""), intern=TRUE)

    fileConn<-file(paste(vcf,"-homozygous_mendel_command.sh",sep=""))
    writeLines(c("#!/bin/bash",homozygous_mendel_command), fileConn)
    close(fileConn)
    
    system(command=paste("sh ",vcf,"-homozygous_mendel_command.sh",sep=""), intern=TRUE)

    current_job_done=grepl("COMPLETED", system(command=paste("sacct --format State -j ",strsplit(current_job, " ")[[1]][4], sep=""), intern=TRUE)[3])
    while(!current_job_done)
        {
        current_job_done=grepl("COMPLETED", system(command=paste("sacct --format State -j ",strsplit(current_job, " ")[[1]][4], sep=""), intern=TRUE)[3])
        Sys.sleep(60)
    }
    
    system(command=paste("rm -f ./",vcf,".pureGT",sep=""), intern = TRUE)
    system(command=paste("rm -f ./",vcf,".pureGT.sample",sep=""), intern = TRUE)
    system(command=paste("gzip -f ",vcf,"-ahmm.new_overlap",sep=""), intern = TRUE)
    system(command=paste("gzip -f ",vcf,"-ahmm.mendel_homo_error",sep=""), intern = TRUE)
    }


InlineHomozygousMendel=function()
    {
    
    
    homozygous_mendel_command=paste('zcat ',vcf,' | grep -v ^# | cut -f10- > ',vcf,'.pureGT;
    cat ',vcf,'.pureGT | awk ',"'BEGIN {srand()} !/^$/ { if (rand() <= ", str_replace(as.character(downsample_for_homozygous_mendel),"0",""),") print $0}'",' > ',vcf,'.pureGT.sample; ',homozygous_mendel_script,' ',vcf,'.pureGT.sample ',folder,' ',vcf,'.samples ', downsample_for_homozygous_mendel,' "', vcf,'"', sep="")
    homozygous_mendel_command=gsub("\r?\n|\r", " ", homozygous_mendel_command)
    
    system(command=paste("rm -f ",vcf,"-homozygous_mendel_command.sh",sep=""), intern=TRUE)
    
    fileConn<-file(paste(vcf,"-homozygous_mendel_command.sh",sep=""))
    writeLines(c("#!/bin/bash",homozygous_mendel_command), fileConn)
    close(fileConn)
    
    
    #print(system(command=paste("sh ",vcf,"-homozygous_mendel_command.sh",sep=""), intern=TRUE))
    
    exec_wait(cmd="sh", args=c(paste(vcf,"-homozygous_mendel_command.sh",sep="")), std_out=TRUE, std_err=TRUE)
    
    system(command=paste("rm -f ./",vcf,".pureGT",sep=""), intern = TRUE)
    system(command=paste("rm -f ./",vcf,".pureGT.sample",sep=""), intern = TRUE)
    system(command=paste("gzip -f ",vcf,"-ahmm.new_overlap",sep=""), intern = TRUE)
    system(command=paste("gzip -f ",vcf,"-ahmm.mendel_homo_error",sep=""), intern = TRUE)
    
    }


MergeHomozygousResults=function(x, overlapping_loci_done, homo_mendel, individuals, df)
    {
    
    i=x
    
    return(data.frame(Overl=as.numeric(overlapping_loci_done[which(individuals$V1 == df$ID1[i]), which(individuals$V1 == df$ID2[i])]), 
                                      ID1L=as.numeric(overlapping_loci_done[which(individuals$V1 == df$ID1[i]), which(individuals$V1 == df$ID1[i])]),
                                      ID2L=as.numeric(overlapping_loci_done[which(individuals$V1 == df$ID2[i]), which(individuals$V1 == df$ID2[i])]),
                                      HM=as.numeric(homo_mendel[which(individuals$V1 == df$ID1[i]), which(individuals$V1 == df$ID2[i])])))
    
    }


ReadHomozygousMendel=function(df)
    {
    print("Reading Homozygous Mendel Results...")
    overlapping_loci_done=as.matrix(fread(paste(folder,vcf,"-ahmm.new_overlap.gz",sep=""), header = FALSE))
     df$overlapping_loci=NA
    homo_mendel= as.matrix(fread(paste(folder,vcf,"-ahmm.mendel_homo_error.gz",sep=""), header = FALSE))
    df$homo_mendel=NA
     print("Merging into data frame...")
    
    merged_homo=bind_rows(mclapply(1:nrow(df), FUN=MergeHomozygousResults, overlapping_loci_done=overlapping_loci_done, homo_mendel=homo_mendel, individuals=individuals, df=df, mc.cores=max_cores))
                 
     df$overlapping_loci=merged_homo$Overl
                          df$ID1_loci=merged_homo$ID1L
                          df$ID2_loci=merged_homo$ID2L
                          df$homo_mendel=merged_homo$HM
    
#     for(i in 1:nrow(df))
#         {
#         df$overlapping_loci[i] = as.numeric(overlapping_loci_done[which(individuals$V1 == df$ID1[i]), which(individuals$V1 == df$ID2[i])])
#         df$ID1_loci[i]=as.numeric(overlapping_loci_done[which(individuals$V1 == df$ID1[i]), which(individuals$V1 == df$ID1[i])])
#         df$ID2_loci[i]=as.numeric(overlapping_loci_done[which(individuals$V1 == df$ID2[i]), which(individuals$V1 == df$ID2[i])])
#         df$homo_mendel[i] = as.numeric(homo_mendel[which(individuals$V1 == df$ID1[i]), which(individuals$V1 == df$ID2[i])])
#     }
    df$overlapping_loci[df$overlapping_loci == 0] = NA
     df$homo_mendel_rel = df$homo_mendel / df$overlapping_loci
    
    
    df$mean_loci=(df$ID1_loci+df$ID2_loci)/2
    
    before=nrow(df)

    df=subset(df, ID1_loci >= min_loci & ID2_loci >= min_loci)
    
    print(paste("Removed ", before-nrow(df), " rows because of too few loci [set in configuration file].",sep=""))
    
    
    return(df)
    }

computeROC=function(x, dataframe, column, df_to_use, countIndiv)
    {
#     dataframe=dataframe[x,]
#     df_tmp=df_to_use
#     df_tmp$threshold=FALSE
#     df_tmp$threshold[df_tmp[,column] <= dataframe[,column]]=TRUE
#     df_tmp_f=filter(df_tmp, threshold == TRUE)
#     df_tmp_sum1=summarise(group_by(df_tmp_f, ID1), n1=n(), .groups="drop_last")
#     df_tmp_sum2=summarise(group_by(df_tmp_f, ID2), n2=n(), .groups="drop_last")
#     df_tmp_sum=full_join(df_tmp_sum1,df_tmp_sum2, by=c("ID1"="ID2"))
#     df_tmp_sum$n1[is.na(df_tmp_sum$n1)]=0
#     df_tmp_sum$n2[is.na(df_tmp_sum$n2)]=0
#     df_tmp_sum$n=df_tmp_sum$n1+df_tmp_sum$n2
    
#     dataframe$AtLeastOnePOFocal=sum(df_tmp_sum$n>=1)/length(unique(c(df$ID1,df$ID2)))
    
#     dataframe$moreThanTwoPOFocal=sum(df_tmp_sum$n>when_PO_error)/length(unique(c(df$ID1,df$ID2)))
        
#     return(dataframe)
    
    dataframe=dataframe[x,]
    df_tmp=df_to_use[df_to_use[,column] <= dataframe[,column], c("ID1","ID2")]
    df_tmp_freq=as.data.frame(table(c(df_tmp$ID1,df_tmp$ID2)))
    
    dataframe$AtLeastOnePOFocal=sum(df_tmp_freq$Freq>=1)/countIndiv
    
    dataframe$moreThanTwoPOFocal=sum(df_tmp_freq$Freq>when_PO_error)/countIndiv
        
    return(dataframe)
}

computeAndPlotIBD0threshold=function()
    {
    
    sequence=unique(df$IBD0[order(df$IBD0)])
    
    if(length(sequence)>10000)
        {
        sequence=sequence[1:10000]
        sequence[9000:10000]=seq(sequence[9000], max(df$IBD0, na.rm=TRUE), length.out=1001)
        sequence=sequence[sequence<=1.0 & !is.na(sequence)]
        sequence=sequence[order(sequence)]
    }

    
    #sequence=seq(0, quantile(df$IBD0, probs = 1, na.rm=TRUE), by=(unique(df$IBD0[order(df$IBD0)])[2]-unique(df$IBD0[order(df$IBD0)])[1]))
    #sequence=sequence+(1:length(sequence))*(unique(df$IBD0[order(df$IBD0)])[2]-unique(df$IBD0[order(df$IBD0)])[1])
    #if(length(sequence)>10000)
    #    {
    #    sequence=sequence[1:10000]
    #}
    
    print(head(sequence))
    print(length(sequence))

    IBD0_threshold=data.frame(IBD0=sequence, AtLeastOnePOFocal=NA, moreThanTwoPOFocal=NA)
    countIndiv=length(unique(c(df$ID1,df$ID2)))

    IBD0_threshold=bind_rows(mclapply(1:nrow(IBD0_threshold), FUN=computeROC, dataframe=IBD0_threshold, column="IBD0", df_to_use=df, countIndiv=countIndiv, mc.cores = max_cores))

    IBD0_threshold$distance=IBD0_threshold$AtLeastOnePOFocal-IBD0_threshold$moreThanTwoPOFocal
    
     #print(summary(IBD0_threshold$distance))
        IBD0_threshold$distance_raw=IBD0_threshold$distance

    #IBD0_threshold$distance=predict(loess(distance~IBD0,IBD0_threshold,span=0.1))
    
    #print(summary(IBD0_threshold$distance))


    #ggplot(IBD0_threshold, aes(moreThanTwoPOFocal,AtLeastOnePOFocal,color=IBD0))+
    #geom_line()
    if(plots==TRUE)
        {
        
        p<-ggplot(IBD0_threshold, aes(IBD0,AtLeastOnePOFocal))+
        geom_line(color="blue")+
        geom_line(mapping=aes(IBD0, moreThanTwoPOFocal),color="red")+
        geom_line(mapping=aes(IBD0, distance),color="black")+
                                geom_line(mapping=aes(IBD0, distance_raw),color="grey")+
geom_label(x=quantile(IBD0_threshold$IBD0, 0.2), y=quantile(IBD0_threshold$AtLeastOnePOFocal, 0.75), label=IBD0_threshold$IBD0[which(IBD0_threshold$distance==max(IBD0_threshold$distance))[length(which(IBD0_threshold$distance==max(IBD0_threshold$distance)))]])+
        ylab("")+
        geom_point(mapping=aes(x=IBD0_threshold$IBD0[which(IBD0_threshold$distance==max(IBD0_threshold$distance))[length(which(IBD0_threshold$distance==max(IBD0_threshold$distance)))]],
                              y=max(IBD0_threshold$distance)), color="black")
        ggsave(filename = paste(vcf,"-",pedigree_file_add_name,"-IBD0_threshold.png",sep=""), plot = p)
        
        }
    thr=IBD0_threshold$IBD0[which(IBD0_threshold$distance==max(IBD0_threshold$distance))[length(which(IBD0_threshold$distance==max(IBD0_threshold$distance)))]]
    
    system(command=paste("echo ", thr, " > ", vcf,"-",pedigree_file_add_name,"-IBD0_threshold.txt",sep=""), intern=TRUE)
    
    return(thr)

    }


computeAndPlotIBD0IQRthreshold=function()
    {
    
    df_IBD0low=subset(df, IBD0 < 0.5)
    
    # Im doing this because unlike the other measurements, IQR can goes low on both ends of the IBD0 spectrum
    
    sequence=unique(df_IBD0low$IBD0_IQR[order(df_IBD0low$IBD0_IQR)])
    
    if(length(sequence)>10000)
        {
        sequence=sequence[1:10000]
        sequence[9000:10000]=seq(sequence[9000], max(df_IBD0low$IBD0_IQR, na.rm=TRUE), length.out=1001)
        sequence=sequence[sequence<=1.0 & !is.na(sequence)]
        sequence=sequence[order(sequence)]
    }
#     sequence=seq(0, quantile(df$IBD0_IQR, probs = 1, na.rm=TRUE), by=(unique(df$IBD0_IQR[order(df$IBD0_IQR)])[2]-unique(df$IBD0_IQR[order(df$IBD0_IQR)])[1]))
    
#     sequence=sequence+(1:length(sequence))*(unique(df$IBD0_IQR[order(df$IBD0_IQR)])[2]-unique(df$IBD0_IQR[order(df$IBD0_IQR)])[1])
    
#     if(length(sequence)>10000)
#         {
#         sequence=sequence[1:10000]
#     }
    
    print(head(sequence))
    print(length(sequence))

    IBD0_IQR_threshold=data.frame(IBD0_IQR=sequence, AtLeastOnePOFocal=NA, moreThanTwoPOFocal=NA)
    countIndiv=length(unique(c(df$ID1,df$ID2)))

    IBD0_IQR_threshold=bind_rows(mclapply(1:nrow(IBD0_IQR_threshold), FUN=computeROC, dataframe=IBD0_IQR_threshold, column="IBD0_IQR", df_to_use=df_IBD0low, countIndiv=countIndiv, mc.cores = max_cores))

    IBD0_IQR_threshold$distance=IBD0_IQR_threshold$AtLeastOnePOFocal-IBD0_IQR_threshold$moreThanTwoPOFocal
    IBD0_IQR_threshold$distance_raw=IBD0_IQR_threshold$distance

    #IBD0_IQR_threshold$distance=predict(loess(distance~IBD0_IQR,IBD0_IQR_threshold,span=0.1))


    #ggplot(IBD0_IQR_threshold, aes(moreThanTwoPOFocal,AtLeastOnePOFocal,color=IBD0_IQR))+
    #geom_line()
    if(plots==TRUE)
        {
        p<-ggplot(IBD0_IQR_threshold, aes(IBD0_IQR,AtLeastOnePOFocal))+
        geom_line(color="blue")+
        geom_line(mapping=aes(IBD0_IQR, moreThanTwoPOFocal),color="red")+
        geom_line(mapping=aes(IBD0_IQR, distance),color="black")+
                        geom_line(mapping=aes(IBD0_IQR, distance_raw),color="grey")+
geom_label(x=quantile(IBD0_IQR_threshold$IBD0_IQR, 0.2), y=quantile(IBD0_IQR_threshold$AtLeastOnePOFocal, 0.75), label=IBD0_IQR_threshold$IBD0_IQR[which(IBD0_IQR_threshold$distance==max(IBD0_IQR_threshold$distance))[length(which(IBD0_IQR_threshold$distance==max(IBD0_IQR_threshold$distance)))]])+
        ylab("")+
        geom_point(mapping=aes(x=IBD0_IQR_threshold$IBD0_IQR[which(IBD0_IQR_threshold$distance==max(IBD0_IQR_threshold$distance))[length(which(IBD0_IQR_threshold$distance==max(IBD0_IQR_threshold$distance)))]],
                              y=max(IBD0_IQR_threshold$distance)), color="black")
        ggsave(filename = paste(vcf,"-",pedigree_file_add_name,"-IBD0_IQR_threshold.png",sep=""), plot = p)
        }
    thr=IBD0_IQR_threshold$IBD0_IQR[which(IBD0_IQR_threshold$distance==max(IBD0_IQR_threshold$distance))[length(which(IBD0_IQR_threshold$distance==max(IBD0_IQR_threshold$distance)))]]
    
    system(command=paste("echo ", thr, " > ", vcf,"-",pedigree_file_add_name,"-IBD0_IQR_threshold.txt",sep=""), intern=TRUE)
    
    return(thr)

    }

# computeAndPlotIBD0sdthreshold=function()
#     {
#     sequence=seq(0, quantile(df$IBD0_sd, probs = 1, na.rm=TRUE), by=0.001*max(df$IBD0_sd, na.rm=TRUE))
#     if(length(sequence)>10000)
#         {
#         sequence=sequence[1:10000]
#     }

#     IBD0_sd_threshold=data.frame(IBD0_sd=sequence, AtLeastOnePOFocal=NA, moreThanTwoPOFocal=NA)

#     IBD0_sd_threshold=bind_rows(mclapply(1:nrow(IBD0_sd_threshold), FUN=computeROC, dataframe=IBD0_sd_threshold, column="IBD0_sd", mc.cores = max_cores))

#     IBD0_sd_threshold$distance=IBD0_sd_threshold$AtLeastOnePOFocal-IBD0_sd_threshold$moreThanTwoPOFocal

#     #IBD0_sd_threshold$distance=predict(loess(distance~IBD0_sd,IBD0_sd_threshold,span=1))


#     #ggplot(IBD0_sd_threshold, aes(moreThanTwoPOFocal,AtLeastOnePOFocal,color=IBD0_sd))+
#     #geom_line()
#     if(plots==TRUE)
#         {
#         p<-ggplot(IBD0_sd_threshold, aes(IBD0_sd,AtLeastOnePOFocal))+
#         geom_line(color="blue")+
#         geom_line(mapping=aes(IBD0_sd, moreThanTwoPOFocal),color="red")+
#         geom_line(mapping=aes(IBD0_sd, distance),color="black")+
#         geom_label(x=quantile(IBD0_sd_threshold$IBD0_sd, 0.2), y=quantile(IBD0_sd_threshold$AtLeastOnePOFocal, 0.75), label=IBD0_sd_threshold$IBD0_sd[which(IBD0_sd_threshold$distance==max(IBD0_sd_threshold$distance))[length(which(IBD0_sd_threshold$distance==max(IBD0_sd_threshold$distance)))]])+
#         ylab("")+
#         geom_point(mapping=aes(x=IBD0_sd_threshold$IBD0_sd[which(IBD0_sd_threshold$distance==max(IBD0_sd_threshold$distance))[length(which(IBD0_sd_threshold$distance==max(IBD0_sd_threshold$distance)))]],
#                               y=max(IBD0_sd_threshold$distance)), color="black")
#         ggsave(filename = paste(vcf,"-",pedigree_file_add_name,"-IBD0_sd_threshold.png",sep=""), plot = p)
#         }
    
#     thr=IBD0_sd_threshold$IBD0_sd[which(IBD0_sd_threshold$distance==max(IBD0_sd_threshold$distance))[length(which(IBD0_sd_threshold$distance==max(IBD0_sd_threshold$distance)))]]
    
#         system(command=paste("echo ", thr, " > ", vcf,"-",pedigree_file_add_name,"-IBD0_sd_threshold.txt",sep=""), intern=TRUE)

    
    
#     return(thr)

#     }

computeAndPlothomomendelrelthreshold=function()
    {
    
    
    sequence=unique(df$homo_mendel_rel[order(df$homo_mendel_rel)])
    
    if(length(sequence)>10000)
        {
        sequence=sequence[1:10000]
        sequence[9000:10000]=seq(sequence[9000], max(df$homo_mendel_rel, na.rm=TRUE), length.out=1001)
        sequence=sequence[sequence<=1.0 & !is.na(sequence)]
        sequence=sequence[order(sequence)]
    }
    
#     sequence=seq(0, quantile(df$homo_mendel_rel, probs = 1, na.rm=TRUE), by=(unique(df$homo_mendel_rel[order(df$homo_mendel_rel)])[2]-unique(df$homo_mendel_rel[order(df$homo_mendel_rel)])[1]))
    
    
#     sequence=sequence+(1:length(sequence))*(unique(df$homo_mendel_rel[order(df$homo_mendel_rel)])[2]-unique(df$homo_mendel_rel[order(df$homo_mendel_rel)])[1])
    
   
#     if(length(sequence)>10000)
#         {
#         sequence=sequence[1:10000]
#     }
    
    print(head(sequence))
    print(length(sequence))

    homo_mendel_rel_threshold=data.frame(homo_mendel_rel=sequence, AtLeastOnePOFocal=NA, moreThanTwoPOFocal=NA)
    countIndiv=length(unique(c(df$ID1,df$ID2)))

    homo_mendel_rel_threshold=bind_rows(mclapply(1:nrow(homo_mendel_rel_threshold), FUN=computeROC, dataframe=homo_mendel_rel_threshold, column="homo_mendel_rel", df_to_use=df, countIndiv=countIndiv, mc.cores = max_cores))

    homo_mendel_rel_threshold$distance=homo_mendel_rel_threshold$AtLeastOnePOFocal-homo_mendel_rel_threshold$moreThanTwoPOFocal
    
    homo_mendel_rel_threshold$distance_raw=homo_mendel_rel_threshold$distance

    #homo_mendel_rel_threshold$distance=predict(loess(distance~homo_mendel_rel,homo_mendel_rel_threshold,span=0.1))


    #ggplot(homo_mendel_rel_threshold, aes(moreThanTwoPOFocal,AtLeastOnePOFocal,color=homo_mendel_rel))+
    #geom_line()
    if(plots==TRUE)
        {
        p<-ggplot(homo_mendel_rel_threshold, aes(homo_mendel_rel,AtLeastOnePOFocal))+
        geom_line(color="blue")+
        geom_line(mapping=aes(homo_mendel_rel, moreThanTwoPOFocal),color="red")+
        geom_line(mapping=aes(homo_mendel_rel, distance),color="black")+
                geom_line(mapping=aes(homo_mendel_rel, distance_raw),color="grey")+
geom_label(x=quantile(homo_mendel_rel_threshold$homo_mendel_rel, 0.2), y=quantile(homo_mendel_rel_threshold$AtLeastOnePOFocal, 0.75), label=homo_mendel_rel_threshold$homo_mendel_rel[which(homo_mendel_rel_threshold$distance==max(homo_mendel_rel_threshold$distance))[length(which(homo_mendel_rel_threshold$distance==max(homo_mendel_rel_threshold$distance)))]])+
        ylab("")+
        geom_point(mapping=aes(x=homo_mendel_rel_threshold$homo_mendel_rel[which(homo_mendel_rel_threshold$distance==max(homo_mendel_rel_threshold$distance))[length(which(homo_mendel_rel_threshold$distance==max(homo_mendel_rel_threshold$distance)))]],
                              y=max(homo_mendel_rel_threshold$distance)), color="black")
        ggsave(filename = paste(vcf,"-",pedigree_file_add_name,"-homo_mendel_rel_threshold.png",sep=""), plot = p)
        }
    
    thr=homo_mendel_rel_threshold$homo_mendel_rel[which(homo_mendel_rel_threshold$distance==max(homo_mendel_rel_threshold$distance))[length(which(homo_mendel_rel_threshold$distance==max(homo_mendel_rel_threshold$distance)))]]
    
            system(command=paste("echo ", thr, " > ", vcf,"-",pedigree_file_add_name,"-homo_mendel_rel_threshold.txt",sep=""), intern=TRUE)

    return(thr)

    }

PlotIBD2DPThreshold=function()
    {
    
    if(plots==TRUE)
        {
        png(paste(vcf,"-",pedigree_file_add_name,"-IBD2_DP_Threshold.png"))
        hist(df$IBD2, 100)
        abline(v = IBD2_DP_Threshold, col="red")
        dev.off()
        }
    }

ClassifyRelationshipCandidates=function(df)
    {
    df$suspect=NA
    if(no_IQR==TRUE)
        {
            df$suspect[(df$IBD0 <= (IBD0_threshold) | df$homo_mendel_rel <= (homo_mendel_rel_threshold))]="PO"
        }else{
            df$suspect[(df$IBD0 <= (IBD0_threshold) | df$homo_mendel_rel <= (homo_mendel_rel_threshold)) | df$IBD0_IQR <= (IBD0_IQR_threshold)]="PO"
        }
    
    df$suspect[df$IBD2 >= IBD2_DP_Threshold]="DP"

    print(summary(as.factor(df$suspect)))
    
    return(df)
    }

duplicate_Q_comparisons=function(x, df_mendel_print)
    {
    return(paste0(naturalsort(c(df_mendel_print$Father[x], df_mendel_print$Mother[x])), collapse="-"))
}

PED_to_trios_txt=function(ped_file)
{
    trios_file=paste(ped_file,".trios.txt",sep="")
  if(previously_computed_mendel != TRUE)
    {
      tmp_ped=fread(ped_file,header=FALSE)
      tmp_ped=tmp_ped[,c(4,3,2)]
  
      fwrite(x = tmp_ped, file = paste(ped_file,".trios.txt",sep=""), sep=",", col.names = FALSE, row.names = FALSE)
      }
  return(trios_file)
}
split_trios_txt=function(trios_file)
{
  system(command=paste("split -n l/",max_cores," ",trios_file," ",trios_file,".",sep=""))
}
InlineBcftoolsMendelian=function(trios_file)
{
  if(max_cores==1)
  {
    system(command=paste("bcftools +mendelian ",vcf," -T ",trios_file," -c > ",trios_file,".out",sep=""))
  }else
  {
    split_trios_txt(trios_file)
    split_files=system(command=paste("ls ", trios_file, ".??",sep=""), intern=TRUE)
    for(spf in split_files)
    {
      if(spf == split_files[1])
      {
        bcfmendelian_cmd=paste("bcftools +mendelian ",vcf," -T ",spf," -c > ",spf,".out", sep="")
        
      }else
      {
        bcfmendelian_cmd=paste(bcfmendelian_cmd, "& bcftools +mendelian ",vcf," -T ",spf," -c > ",spf,".out",sep="")
      }
     
    }
    bcfmendelian_cmd=paste(bcfmendelian_cmd, " & wait",sep="")
    system(command=bcfmendelian_cmd, intern=TRUE)
  }
}

convertOutputToOldFormat=function(trios_file)
{
  if(max_cores==1)
  {
    all_mendel_summary=fread(paste(trios_file,".out",sep=""), header=TRUE, data.table=FALSE)
    
  }else{
    split_files=system(command=paste("ls ", trios_file, ".??.out",sep=""), intern=TRUE)
    for(spf in split_files)
    {
      if(spf == split_files[1])
      {
        all_mendel_summary=fread(spf, header=TRUE, data.table=FALSE)
      }else{
        all_mendel_summary=bind_rows(all_mendel_summary,fread(spf, header=TRUE, data.table=FALSE))
      }
    }
    all_mendel_summary=unique(all_mendel_summary)
  }
  
  all_mendel_summary$V1=all_mendel_summary[,2]/(all_mendel_summary[,1]+all_mendel_summary[,2])
  all_mendel_summary=tidyr::separate(all_mendel_summary, "[4]Trio (mother,father,child)", into=c("P1","P2","O"), sep=",")
  all_mendel_summary$V2=paste(all_mendel_summary$O,all_mendel_summary$P1,all_mendel_summary$P2,sep="_")
  all_mendel_summary=all_mendel_summary[,c("V1","V2")]
  
  return(all_mendel_summary)
}

SbatchMendelTrios=function(Genomics_Sex, focals, df)
    {
    
    # still vcftools
    
    system(command=paste("rm -f ",vcf,pedigree_file_add_name,"_all_mendel.ped",sep=""), intern=TRUE)

for(focal in focals)
    {
    #print(focal)
    indiv_id=c(focal)

    PID=c(df$ID1[df$suspect=="PO" & !is.na(df$suspect) & (df$ID1 %in% indiv_id | df$ID2 %in% indiv_id)],df$ID2[df$suspect=="PO" & !is.na(df$suspect) & (df$ID1 %in% indiv_id | df$ID2 %in% indiv_id)])

    PID=PID[!(PID %in% indiv_id)]
    PID

    df_mendel_print=expand.grid(Family="all", Offspring=indiv_id, Father=PID[PID %in% as.matrix(Genomics_Sex[Genomics_Sex$Genomics_Sex %in% c("M","Q"), c("indv")])], Mother=PID[PID %in% as.matrix(Genomics_Sex[Genomics_Sex$Genomics_Sex %in% c("F","Q"), c("indv")])], stringsAsFactors=FALSE)

    df_mendel_print=distinct(df_mendel_print)
    df_mendel_print=subset(df_mendel_print, Father != Mother)
    
    
    if(nrow(df_mendel_print)>1)
        {
        # if many Q sexes, then many trios will be computed twice. This will avoid that
        df_mendel_print$Trio = mclapply(1:nrow(df_mendel_print), FUN=duplicate_Q_comparisons, df_mendel_print=df_mendel_print, mc.cores = max_cores)

        df_mendel_print=subset(df_mendel_print, !duplicated(Trio))

        df_mendel_print=df_mendel_print[,-1*which(colnames(df_mendel_print)=="Trio")]
    }
    
    
    write_tsv(x = df_mendel_print, path = paste(vcf,pedigree_file_add_name,"_all_mendel.ped",sep=""), append = TRUE)
}
    
    thecommand=paste(" date && cd ",folder," && vcftools --gzvcf ",vcf,
                 " --not-chr X --not-chr Y --not-chr MT --mendel ",vcf,pedigree_file_add_name,"_all_mendel.ped -c | cut -f 5 | gzip > ",vcf,pedigree_file_add_name,"_all.mendel.gz && cat ",vcf,pedigree_file_add_name,"_all_mendel.ped | cut -f 2,3,4 | sed -e 's/\t/_/g' | gzip >> ",vcf,pedigree_file_add_name,"_all.mendel.gz && zcat ",vcf,pedigree_file_add_name,"_all.mendel.gz | sort -T /moto/ziab/users/jr3950/data/genomes/tmp/ --buffer-size=",as.character(as.numeric(str_replace(max_memory, "G", ""))-1),"g | uniq -c | sort -nr | gzip > ",vcf,pedigree_file_add_name,"_all.mendel.summary.gz && rm -f ",vcf,pedigree_file_add_name,"_all.mendel.gz",sep="")

    thecommand=paste('sbatch -t 0-11:59:00 --job-name=MendelFx --mem=',max_memory,' -c 1 -A ',truffle_slurm_account,' --wrap "',samtools_activation,'',thecommand,'"', sep="")

    thecommand=gsub("\r?\n|\r", " ", thecommand)

    current_job<-system(command=thecommand, intern=TRUE)
     
     
     current_job_done=grepl("COMPLETED", system(command=paste("sacct --format State -j ",strsplit(current_job, " ")[[1]][4], sep=""), intern=TRUE)[3])
    while(!current_job_done)
        {
        current_job_done=grepl("COMPLETED", system(command=paste("sacct --format State -j ",strsplit(current_job, " ")[[1]][4], sep=""), intern=TRUE)[3])
        Sys.sleep(60)
    }
    }


InlineMendelTrios=function(Genomics_Sex, focals, df)
    {
    if(previously_computed_mendel != TRUE)
    {
        system(command=paste("rm -f ",vcf,pedigree_file_add_name,"_all_mendel.ped",sep=""), intern=TRUE)

        if(length(focals)==0)
            {
            stop("No PO candidates found?")
            }

    #     if(folder=='/moto/ziab/users/jr3950/data/OtherPedigrees/Cattle/WGS_Holstein_cattle_inbreeding/converted/')
    #         {

    #         Genomics_Sex$indv=str_replace_all(Genomics_Sex$indv, "1_T", "T")
    #         Genomics_Sex$indv=str_replace_all(Genomics_Sex$indv, "1_P", "P")

    #         df$ID1=str_replace_all(df$ID1, "1_T", "T")
    #         df$ID1=str_replace_all(df$ID1, "1_P", "P")

    #         df$ID2=str_replace_all(df$ID2, "1_T", "T")
    #         df$ID2=str_replace_all(df$ID2, "1_P", "P")
    #     }


    for(focal in focals)
        {
        #print(focal)
        indiv_id=c(focal)

        PID=c(df$ID1[df$suspect=="PO" & !is.na(df$suspect) & (df$ID1 %in% indiv_id | df$ID2 %in% indiv_id)],df$ID2[df$suspect=="PO" & !is.na(df$suspect) & (df$ID1 %in% indiv_id | df$ID2 %in% indiv_id)])

        PID=PID[!(PID %in% indiv_id)]
        PID

        df_mendel_print=expand.grid(Family="all", Offspring=indiv_id, Father=PID[PID %in% as.matrix(Genomics_Sex[Genomics_Sex$Genomics_Sex %in% c("M","Q"), c("indv")])], Mother=PID[PID %in% as.matrix(Genomics_Sex[Genomics_Sex$Genomics_Sex %in% c("F","Q"), c("indv")])], stringsAsFactors=FALSE)

        df_mendel_print=distinct(df_mendel_print)
        df_mendel_print=subset(df_mendel_print, Father != Mother)


        if(nrow(df_mendel_print)>1)
            {
            # if many Q sexes, then many trios will be computed twice. This will avoid that
            df_mendel_print$Trio = mclapply(1:nrow(df_mendel_print), FUN=duplicate_Q_comparisons, df_mendel_print=df_mendel_print, mc.cores = max_cores)

            df_mendel_print=subset(df_mendel_print, !duplicated(Trio))

            df_mendel_print=df_mendel_print[,-1*which(colnames(df_mendel_print)=="Trio")]
        }


        write_tsv(x = df_mendel_print, file = paste(vcf,pedigree_file_add_name,"_all_mendel.ped",sep=""), append = TRUE)
    }
        }
    
    trios_file=PED_to_trios_txt(paste(vcf,pedigree_file_add_name,"_all_mendel.ped",sep=""))
    
#     thecommand=paste("date && cd ",folder," && vcftools --gzvcf ",vcf,
#                  " --not-chr X --not-chr Y --not-chr MT --mendel ",vcf,pedigree_file_add_name,"_all_mendel.ped -c | cut -f 5 | gzip > ",vcf,pedigree_file_add_name,"_all.mendel.gz && cat ",vcf,pedigree_file_add_name,"_all_mendel.ped | cut -f 2,3,4 | sed -e 's/\t/_/g' | gzip >> ",vcf,pedigree_file_add_name,"_all.mendel.gz && zcat ",vcf,pedigree_file_add_name,"_all.mendel.gz | sort -T /moto/ziab/users/jr3950/data/genomes/tmp/ --buffer-size=",as.character(as.numeric(str_replace(max_memory, "G", ""))-1),"g | uniq -c | sort -nr | gzip > ",vcf,pedigree_file_add_name,"_all.mendel.summary.gz && rm -f ",vcf,pedigree_file_add_name,"_all.mendel.gz",sep="")


    #thecommand=gsub("\r?\n|\r", " ", thecommand)
    
    

   #print(system(command=thecommand, intern=TRUE))
    
    if(previously_computed_mendel != TRUE)
    {
    
    InlineBcftoolsMendelian(trios_file)
        
        }
    
    return(trios_file)
    }

getMendelID=function(x)
        {
        return(as.character(strsplit(x, "_", TRUE)[[1]][1]))
    }

GetMendelTriosData=function(trios_file)
{
  # df$mean_mendel_fwd=NA
  # df$min_mendel_fwd=NA
  # 
  # df$mean_mendel_rev=NA
  # df$min_mendel_rev=NA
  # 
  # df$mean_mendel=NA
  # df$min_mendel=NA
  
  df$ID1=as.character(df$ID1)
  df$ID2=as.character(df$ID2)
  
  #all_mendel_summary=fread(paste(vcf,pedigree_file_add_name,"_all.mendel.summary.gz",sep=""), data.table=FALSE)
  all_mendel_summary=convertOutputToOldFormat(trios_file)
  
  #     if(folder=='/moto/ziab/users/jr3950/data/OtherPedigrees/Cattle/WGS_Holstein_cattle_inbreeding/converted/')
  #         {
  #         all_mendel_summary$V2=str_replace_all(all_mendel_summary$V2, "1_T", "T")
  #         all_mendel_summary$V2=str_replace_all(all_mendel_summary$V2, "1_P", "P")
  
  #         df$ID1=str_replace_all(df$ID1, "1_T", "T")
  #         df$ID1=str_replace_all(df$ID1, "1_P", "P")
  
  #         df$ID2=str_replace_all(df$ID2, "1_T", "T")
  #         df$ID2=str_replace_all(df$ID2, "1_P", "P")
  #     }
  
  
  
  all_mendel_summary$ID=unlist(mclapply(all_mendel_summary$V2, getMendelID, mc.cores=max_cores))
  
  #all_mendel_summary$V1=all_mendel_summary$V1-1
  
  # length(unique(all_mendel_summary$ID))
  # 
  # head(all_mendel_summary)
  # 
  # date()
  # for(cur in unique(all_mendel_summary$ID))
  # {
  #   mendel_summary=subset(all_mendel_summary, ID == cur)
  #   ID=cur
  #   #print(ID)
  #   mendel_summary$V2 = paste(mendel_summary$V2, "_", sep="")
  #   for(i in which(df$ID1 == ID | df$ID2 == ID))
  #   {
  #     
  #     
  #     if(length(mendel_summary$V1[str_detect(mendel_summary$V2, pattern = paste(df$ID1[i],"_",sep="")) & str_detect(mendel_summary$V2, pattern = paste(df$ID2[i],"_",sep=""))]) > 0)
  #       
  #     {
  #       
  #       if(length(mendel_summary$V1[startsWith(mendel_summary$V2, prefix=df$ID1[i]) & str_detect(mendel_summary$V2, pattern = paste(df$ID1[i],"_",sep="")) & str_detect(mendel_summary$V2, pattern = paste(df$ID2[i],"_",sep=""))]) > 0)
  #       {
  #         df$mean_mendel_fwd[i]=mean(mendel_summary$V1[startsWith(mendel_summary$V2, prefix=df$ID1[i]) & str_detect(mendel_summary$V2, pattern = paste(df$ID1[i],"_",sep="")) & str_detect(mendel_summary$V2, pattern = paste(df$ID2[i],"_",sep=""))], na.rm=TRUE)
  #         df$min_mendel_fwd[i]=min(mendel_summary$V1[startsWith(mendel_summary$V2, prefix=df$ID1[i]) & str_detect(mendel_summary$V2, pattern = paste(df$ID1[i],"_",sep="")) & str_detect(mendel_summary$V2, pattern = paste(df$ID2[i],"_",sep=""))], na.rm=TRUE)
  #         
  #       }
  #       
  #       if(length(mendel_summary$V1[startsWith(mendel_summary$V2, prefix=df$ID2[i]) & str_detect(mendel_summary$V2, pattern = paste(df$ID1[i],"_",sep="")) & str_detect(mendel_summary$V2, pattern = paste(df$ID2[i],"_",sep=""))]) > 0)
  #       {
  #         
  #         df$mean_mendel_rev[i]=mean(mendel_summary$V1[startsWith(mendel_summary$V2, prefix=df$ID2[i]) & str_detect(mendel_summary$V2, pattern = paste(df$ID1[i],"_",sep="")) & str_detect(mendel_summary$V2, pattern = paste(df$ID2[i],"_",sep=""))], na.rm=TRUE)
  #         df$min_mendel_rev[i]=min(mendel_summary$V1[startsWith(mendel_summary$V2, prefix=df$ID2[i]) & str_detect(mendel_summary$V2, pattern = paste(df$ID1[i],"_",sep="")) & str_detect(mendel_summary$V2, pattern = paste(df$ID2[i],"_",sep=""))], na.rm=TRUE)
  #         
  #       }
  #       
  #       df$min_mendel[i]=min(c(df$min_mendel_fwd[i],df$min_mendel_rev[i]), na.rm=TRUE)
  #       df$mean_mendel[i]=mean(c(df$mean_mendel_fwd[i],df$mean_mendel_rev[i]), na.rm=TRUE)
  #       
  #     }
  #     
  #     
  #     
  #   }
  # }
  # date()
  # 
  # df$min_mendel_rel[!is.na(df$min_mendel) & !is.na(df$overlapping_loci)] = df$min_mendel[!is.na(df$min_mendel) & !is.na(df$overlapping_loci)]/df$overlapping_loci[!is.na(df$min_mendel) & !is.na(df$overlapping_loci)]
  # df$min_mendel_rel[df$min_mendel_rel > 1] = NA
  # summary(df$min_mendel_rel[!is.na(df$min_mendel) & !is.na(df$overlapping_loci)])
  # 
    all_mendel_summary=subset(all_mendel_summary, V2 != "FAMILY")
    all_mendel_summary=tidyr::separate(all_mendel_summary, V2, c("O", "P1", "P2"), "_", remove=FALSE)
    

  return(list(df,all_mendel_summary))
}

SearchForDuplicateSamples=function()
    {
    pedigree=data.frame(Family="All", ID=unique(c(df$ID1,df$ID2)), Father=NA, Mother=NA, Sex=NA, stringsAsFactors = FALSE)
pedigree$Dup1=NA
    for(i in which(df$suspect=="DP"))
        {
        current_row=NA
        if(df$ID1[i] %in% pedigree$ID)
            {
            current_row=which(pedigree$ID==df$ID1[i])
        }
        if(!(df$ID1[i] %in% pedigree$ID))
            {
            for(j in which(grepl("Dup",colnames(pedigree))))
                {
                if(df$ID1[i] %in% pedigree[,j])
                    {
                        current_row=which(pedigree[,j]==df$ID1[i])
                }
            }
        }
        #print(current_row)
        if(!is.na(current_row))
            {
            if(is.na(pedigree$Dup1[current_row]))
                {
                #print(df[i,c("ID1", "ID2")])
                #print(pedigree[current_row,])
                if(!(df$ID2[i] %in% t(pedigree[current_row,])))
                    {
                    pedigree$Dup1[current_row]=df$ID2[i]
                    #print(pedigree[current_row,])
                    pedigree=subset(pedigree, ID != df$ID2[i])
                }
                # current_row is now wrong and should not be re-used!

            }
            else {
                dups=2
                while(TRUE)
                    {
                    if(!(paste("Dup",dups,sep="") %in% colnames(pedigree)))
                        {
                            pedigree[, paste("Dup",dups,sep="")] = NA
                        }
                    if(is.na(pedigree[current_row, paste("Dup",dups,sep="")]))
                        {
                        #print(i)
                        if(!(df$ID2[i] %in% t(pedigree[current_row,])))
                        {
                                            #print(!(df$ID2[i] %in% pedigree[current_row,]))

                            pedigree[current_row, paste("Dup",dups,sep="")]=df$ID2[i]
                            pedigree=subset(pedigree, ID != df$ID2[i])
                            }
                        # current_row is now wrong and should not be re-used!
                        break
                    }
                    else
                        {
                        dups=dups+1
                    }
                }
            }
        }
        if(is.na(current_row))
            {
            #print(i)
        }


    }
    return(pedigree)
    }



GatherDeepMendelianData=function(x,no_dups)
    {
    # The following code can be drastically simplified, but is pretty nasty for historic reasons.
    
    
    #print(which_to_look_at)
    #which_to_look_at=which(!is.na(df$suspect) & df$suspect == "PO" & (df$ID1 == "m25_39" | df$ID2 == "m25_39"))
        PO_test_ID=x
            #PO_test_ID=df[row,row_col]
            #PO_test_ID=as.character(PO_test_ID)
            
            if(no_dups)
                {
                cur_focal=c(PO_test_ID)
            }else{
                for(dup in 1:length(Dup_plus_ID_cols))
                {
                if(dup == 1)
                    {
                    cur_focal=t(pedigree[pedigree[,Dup_plus_ID_cols[dup]]==PO_test_ID, Dup_plus_ID_cols[dup]])
                }else
                    {
                    cur_focal=c(cur_focal, t(pedigree[pedigree[,Dup_plus_ID_cols[dup]]==PO_test_ID, Dup_plus_ID_cols[dup]]))
                }

            }
            }
            

            cur_focal=cur_focal[!is.na(cur_focal)]        
            cur_subset=subset(df, ID1 %in% cur_focal | ID2 %in% cur_focal)
            cur_subset=subset(cur_subset, suspect=="PO")
            cur_POs=unique(c(cur_subset$ID1, cur_subset$ID2))
            cur_POs=cur_POs[!(cur_POs %in% cur_focal)]
            
            
            if(no_dups)
                {
                
                # cur_POs already has the right ids because there are no dups
                
                }else{
                for(i in cur_POs)
                {
                if(i %in% pedigree$ID)
                    {
                    cur_POs=c(cur_POs, t(pedigree[pedigree$ID == i, Dup_plus_ID_cols]))
                }
                else
                    {
                    
                    for(j in colnames(pedigree)[which(startsWith(colnames(pedigree), "Dup"))])
                        {
                        if(i %in% pedigree[,j])
                            {
                            cur_POs=c(cur_POs, t(pedigree[pedigree[,j] == i, Dup_plus_ID_cols]))
                        }
                    }
                }
            }
            }
            
            
            cur_POs=cur_POs[!is.na(cur_POs)]
            cur_POs=unique(cur_POs)
            #cur_POs

            cur_IDs=c(cur_focal, cur_POs)
            #cur_IDs=c(str_replace(cur_IDs, "_", "?"), str_replace(cur_IDs, "_", "?0"))
            #cur_IDs
            #print(cur_IDs)
            #cur_mendels_files=NA
            
            cur_mendels=subset(all_mendel_summary, ID %in% cur_IDs)

            
            #for(i in cur_IDs)
            #    {
            #    if(i == cur_IDs[1])
            #         {
            #        cur_mendels=subset(all_mendel_summary, ID==i)
            #        }
            #    else
            #        {
            #        cur_mendels=bind_rows(cur_mendels,subset(all_mendel_summary, ID==i))
            #        }
            #    }
            #
            #cur_mendels_files
            #cur_mendels=NA
            #if(is.character(cur_mendels_files))
                #{
    #             for(i in cur_mendels_files)
    #             {
    #             if(i == cur_mendels_files[1])
    #                 {
    #                 cur_mendels=fread(i, data.table=FALSE)
    #                 #print(head(cur_mendels))
    #             }
    #             else
    #                 {
    #                 cur_mendels=rbind(cur_mendels,fread(i, data.table=FALSE))
    #                 #print(nrow(cur_mendels))
    #             }
    #         }
            if(exists("cur_mendels"))
            {
                #print(cur_mendels)
                if(nrow(cur_mendels) > 0)
                    {
                    #cur_mendels=subset(cur_mendels, V2 != "FAMILY")
                    #cur_mendels=tidyr::separate(cur_mendels, V2, c("O", "P1", "P2"), "_")
        #             cur_mendels$O=str_replace(cur_mendels$O, "SM_", "")
        #             cur_mendels$O=str_replace(cur_mendels$O, "-", "_")
        #             cur_mendels$O=str_replace(cur_mendels$O, "_0", "_")

        #             cur_mendels$P1=str_replace(cur_mendels$P1, "-", "_")
        #             cur_mendels$P1=str_replace(cur_mendels$P1, "_0", "_")
        #             cur_mendels$P2=str_replace(cur_mendels$P2, "-", "_")
        #             cur_mendels$P2=str_replace(cur_mendels$P2, "_0", "_")
                    #head(cur_mendels)

                    #cur_mendels$V1_percentile=ecdf(as.numeric(cur_mendels$V1))((as.numeric(cur_mendels$V1)))


                    print(cur_focal)
                    #print(cur_POs)

                    cur_mendels=subset(cur_mendels, O %in% cur_focal & P1 %in% cur_POs & P2 %in% cur_POs)
                    
                    if(nrow(cur_mendels) > 0)
                    {
                        
                    
                    cur_mendels$V1_mean=mean((as.numeric(cur_mendels$V1)))
                        if(nrow(cur_mendels)==1){
                            cur_mendels$V1_mean=NA
                            }
                        
                        }else{
                        return(NULL)
                        }


                }


               # }
                
                #mendels=bind_rows(mendels, cur_mendels[,which(colnames(cur_mendels)=="ID")*-1])
                #remove(cur_mendels)
                return(cur_mendels)
               }else{
                return(NULL)
            }

        


    
    
    #print(mendels)
    
    
    
    
    }

computeROC2=function(x, dataframe)
    {
    
#     dataframe=dataframe[x,]
    
#     df_analyze_V1_cutoffs_tmp=summarise(group_by(filter(mendels, V1_V1_mean <= dataframe$V1_V1_mean & 
#                                                        V1_percentile <= dataframe$V1_percentile),
#                                                 O), n=n(), .groups="drop_last")
#     dataframe$focals_with_one_potential_trio=sum(df_analyze_V1_cutoffs_tmp$n == 1)
#     dataframe$focals_with_more_than_one_potential_trio=sum(df_analyze_V1_cutoffs_tmp$n > 1)
        
#     return(dataframe)
    
    dataframe=dataframe[x,]
    
    df_analyze_V1_cutoffs_tmp=mendels[mendels$V1_V1_mean <= dataframe$V1_V1_mean & mendels$V1_percentile <= dataframe$V1_percentile, c("O")]
    
    df_analyze_V1_cutoffs_tmp_freq=as.data.frame(table(df_analyze_V1_cutoffs_tmp))
    
    
    dataframe$focals_with_one_potential_trio=sum(df_analyze_V1_cutoffs_tmp_freq$Freq == 1)
    dataframe$focals_with_more_than_one_potential_trio=sum(df_analyze_V1_cutoffs_tmp_freq$Freq > 1)
        
    return(dataframe)
}

CompareMendelianErrorsForThresholds=function()
    {
    
    V1_V1_means=unique(mendels$V1_V1_mean[order(mendels$V1_V1_mean)])
    V1_V1_means=V1_V1_means[!is.na(V1_V1_means)]
    
    if(length(V1_V1_means)>10000)
        {
        V1_V1_means=V1_V1_means[1:10000]
        V1_V1_means[9000:10000]=seq(V1_V1_means[9000], max(mendels$V1_V1_mean, na.rm=TRUE), length.out=1001)
        V1_V1_means=V1_V1_means[order(V1_V1_means)]
    }

    print(head(V1_V1_means))
    
    V1_percentiles=unique(mendels$V1_percentile[order(mendels$V1_percentile)])
    V1_percentiles=V1_percentiles[!is.na(V1_percentiles)]
    
    if(length(V1_percentiles)>10000)
        {
        V1_percentiles=V1_percentiles[1:10000]
        V1_percentiles[9000:10000]=seq(V1_percentiles[9000], max(mendels$V1_percentile, na.rm=TRUE), length.out=1001)
        V1_percentiles=V1_percentiles[order(V1_percentiles)]
    }
    
    print(head(V1_percentiles))
    
    
#     V1_V1_means=seq(min(mendels$V1_V1_mean), quantile(mendels$V1_V1_mean,0.5), by = 0.001)
#     V1_percentiles=seq(min(mendels$V1_percentile), quantile(mendels$V1_percentile, 0.5), by = 0.001)
#     if(length(V1_V1_means) > 10000)
#         {
#         V1_V1_means=V1_V1_means[1:10000]
#     }
#     if(length(V1_percentiles) > 10000)
#         {
#         V1_percentiles=V1_percentiles[1:10000]
#     }

    # How many possible trios per threshold?

    df_analyze_V1_cutoffs=expand.grid(V1_V1_mean=V1_V1_means, V1_percentile=V1_percentiles, focals_with_one_potential_trio=NA, focals_with_more_than_one_potential_trio=NA)

    df_analyze_V1_cutoffs=subset(df_analyze_V1_cutoffs, V1_V1_mean==max(df_analyze_V1_cutoffs$V1_V1_mean) | V1_percentile==max(df_analyze_V1_cutoffs$V1_percentile))
    
    countIndiv=length(unique(mendels$O))

    df_analyze_V1_cutoffs=bind_rows(mclapply(1:nrow(df_analyze_V1_cutoffs), dataframe=df_analyze_V1_cutoffs, FUN=computeROC2, mc.cores=max_cores))

    # for(i in 1:nrow(df_analyze_V1_cutoffs))
    #     {
    #     df_analyze_V1_cutoffs_tmp=summarise(group_by(filter(mendels, V1_V1_mean <= df_analyze_V1_cutoffs$V1_V1_mean[i] & 
    #                                                        V1_percentile <= df_analyze_V1_cutoffs$V1_percentile[i]),
    #                                                 O), n=n())
    #     df_analyze_V1_cutoffs$focals_with_one_potential_trio[i]=sum(df_analyze_V1_cutoffs_tmp$n == 1)
    #     df_analyze_V1_cutoffs$focals_with_more_than_one_potential_trio[i]=sum(df_analyze_V1_cutoffs_tmp$n > 1)
    # }
    df_analyze_V1_cutoffs$focals_with_one_potential_trio=df_analyze_V1_cutoffs$focals_with_one_potential_trio/countIndiv
    df_analyze_V1_cutoffs$focals_with_more_than_one_potential_trio=df_analyze_V1_cutoffs$focals_with_more_than_one_potential_trio/countIndiv
    
    fwrite(df_analyze_V1_cutoffs, paste(vcf,"-",pedigree_file_add_name,"-df_analyze_V1_cutoffs.csv",sep=""))
    
    #df_analyze_V1_cutoffs$distance=df_analyze_V1_cutoffs$focals_with_one_potential_trio-df_analyze_V1_cutoffs$focals_with_more_than_one_potential_trio

    plot_V1=filter(df_analyze_V1_cutoffs, V1_V1_mean==max(df_analyze_V1_cutoffs$V1_V1_mean))
    plot_V1_mean=filter(df_analyze_V1_cutoffs, V1_percentile==max(df_analyze_V1_cutoffs$V1_percentile))

    #plot_V1$focals_with_one_potential_trio=predict(loess(focals_with_one_potential_trio~V1_percentile,plot_V1, span=1))
    #plot_V1$focals_with_more_than_one_potential_trio=predict(loess(focals_with_more_than_one_potential_trio~V1_percentile,plot_V1, span=1))

    #plot_V1_mean$focals_with_one_potential_trio=predict(loess(focals_with_one_potential_trio~V1_V1_mean,plot_V1_mean, span=1))
    #plot_V1_mean$focals_with_more_than_one_potential_trio=predict(loess(focals_with_more_than_one_potential_trio~V1_V1_mean,plot_V1_mean, span=1))

    plot_V1_calc=subset(plot_V1, focals_with_more_than_one_potential_trio < focals_with_one_potential_trio)
    plot_V1_mean_calc=subset(plot_V1_mean, focals_with_more_than_one_potential_trio < focals_with_one_potential_trio)

    V1_percentile_threshold=plot_V1_calc$V1_percentile[which(plot_V1_calc$focals_with_one_potential_trio==max(plot_V1_calc$focals_with_one_potential_trio))[length(which(plot_V1_calc$focals_with_one_potential_trio==max(plot_V1_calc$focals_with_one_potential_trio)))]]
    V1_V1_mean_threshold=plot_V1_mean_calc$V1_V1_mean[which(plot_V1_mean_calc$focals_with_one_potential_trio==max(plot_V1_mean_calc$focals_with_one_potential_trio))[length(which(plot_V1_mean_calc$focals_with_one_potential_trio==max(plot_V1_mean_calc$focals_with_one_potential_trio)))]]
    
    if(plots==TRUE)
        {
        p<-ggplot(plot_V1, aes(V1_percentile,focals_with_one_potential_trio))+
        geom_line(color="blue")+
        geom_line(mapping=aes(V1_percentile, focals_with_more_than_one_potential_trio),color="red")+
        #geom_line(mapping=aes(V1_percentile, distance),color="black")+
        geom_label(x=quantile(plot_V1$V1_percentile, 0.2), y=quantile(plot_V1$focals_with_one_potential_trio, 0.5), label=V1_percentile_threshold)+
        ylab("")+
        geom_point(mapping=aes(x=V1_percentile_threshold,
                              y=max(plot_V1_calc$focals_with_one_potential_trio)), color="blue")
        ggsave(filename = paste(vcf,"-",pedigree_file_add_name,"-TrioMendelErrorsPercentileThreshold.png",sep=""), plot = p)

        p<-ggplot(plot_V1_mean, aes(V1_V1_mean,focals_with_one_potential_trio))+
        geom_line(color="blue")+
        geom_line(mapping=aes(V1_V1_mean, focals_with_more_than_one_potential_trio),color="red")+
        #geom_line(mapping=aes(V1_V1_mean, distance),color="black")+
        geom_label(x=quantile(plot_V1_mean$V1_V1_mean, 0.2), y=quantile(plot_V1_mean$focals_with_one_potential_trio, 0.5), label=V1_V1_mean_threshold)+
        ylab("")+
        geom_point(mapping=aes(x=V1_V1_mean_threshold,
                              y=max(plot_V1_mean_calc$focals_with_one_potential_trio)), color="blue")
        ggsave(filename = paste(vcf,"-",pedigree_file_add_name,"-TrioMendelErrorsComparedToMeanThreshold.png",sep=""), plot = p)
        }
    
    return(list(V1_V1_mean_threshold, V1_percentile_threshold))
    
    }

IBD0_row_mean=function(x,mendels_123123)
    {
    return(as.numeric(mean(as.numeric(t(mendels_123123[x,which(endsWith(colnames(mendels_123123), "_IBD0"))])), na.rm = TRUE)))
}

IBD0_IQR_row_mean=function(x,mendels_123123)
    {
    return(as.numeric(mean(as.numeric(t(mendels_123123[x,which(endsWith(colnames(mendels_123123), "_IBD0_IQR"))])), na.rm = TRUE)))
}

IBD0_sd_row_mean=function(x,mendels_123123)
    {
    return(as.numeric(mean(as.numeric(t(mendels_123123[x,which(endsWith(colnames(mendels_123123), "_IBD0_sd"))])), na.rm = TRUE)))
}

homo_mendel_rel_row_mean=function(x,mendels_123123)
    {
    return(as.numeric(mean(as.numeric(t(mendels_123123[x,which(endsWith(colnames(mendels_123123), "_homo_mendel_rel"))])), na.rm = TRUE)))
}


DeeperTrioComparison=function(mendels_tmp)
    {
    PO_columns=c("IBD0","IBD0_IQR","homo_mendel_rel")
    
    if(no_IQR==TRUE)
        {
        PO_columns=c("IBD0","homo_mendel_rel")
        }
    
    mendels_tmp$O=as.character(mendels_tmp$O)
    mendels_tmp$P1=as.character(mendels_tmp$P1)
    mendels_tmp$P2=as.character(mendels_tmp$P2)
    
    df$ID1=as.character(df$ID1)
    df$ID2=as.character(df$ID2)
    
    mendels_tmp=left_join(mendels_tmp, df[,c("ID1","ID2", PO_columns)], by=c("O"="ID1", "P1"="ID2"))
    colnames(mendels_tmp)[which(colnames(mendels_tmp) %in% PO_columns)]=paste("Nr1_",colnames(mendels_tmp)[which(colnames(mendels_tmp) %in% PO_columns)],sep="")
    mendels_tmp=left_join(mendels_tmp, df[,c("ID1","ID2", PO_columns)], by=c("O"="ID2", "P1"="ID1"))
    colnames(mendels_tmp)[which(colnames(mendels_tmp) %in% PO_columns)]=paste("Nr2_",colnames(mendels_tmp)[which(colnames(mendels_tmp) %in% PO_columns)],sep="")
    mendels_tmp=left_join(mendels_tmp, df[,c("ID1","ID2", PO_columns)], by=c("O"="ID1", "P2"="ID2"))
    colnames(mendels_tmp)[which(colnames(mendels_tmp) %in% PO_columns)]=paste("Nr3_",colnames(mendels_tmp)[which(colnames(mendels_tmp) %in% PO_columns)],sep="")
    mendels_tmp=left_join(mendels_tmp, df[,c("ID1","ID2", PO_columns)], by=c("O"="ID2", "P2"="ID1"))
    colnames(mendels_tmp)[which(colnames(mendels_tmp) %in% PO_columns)]=paste("Nr4_",colnames(mendels_tmp)[which(colnames(mendels_tmp) %in% PO_columns)],sep="")

    mendels_tmp$IBD0_mean=unlist(mclapply(1:nrow(mendels_tmp), IBD0_row_mean, mendels_123123=mendels_tmp, mc.cores = max_cores))
    mendels_tmp$IBD0_IQR_mean=unlist(mclapply(1:nrow(mendels_tmp), IBD0_IQR_row_mean, mendels_123123=mendels_tmp, mc.cores = max_cores))
    mendels_tmp$homo_mendel_rel_mean=unlist(mclapply(1:nrow(mendels_tmp), homo_mendel_rel_row_mean, mendels_123123=mendels_tmp, mc.cores = max_cores))

    mendels_tmp_trio_suspects=subset(mendels_tmp, V1_V1_mean <= V1_V1_mean_threshold | V1_percentile <= V1_percentile_threshold)

    return(mendels_tmp)
    }


find_pedigree_row=function(x)
    {
    for(dup_col in Dup_plus_ID_cols)
        {
        if(x[1] %in% pedigree[,dup_col])
            {
            return(which(pedigree[,dup_col]%in%x))
        }
    }
}

fill_in_ped=function(ped,data,x)
    {
    #print(head(data))
    #print(find_pedigree_row(data$O[x]))
    if(ped$Sex[find_pedigree_row(data$P1[x])]!=ped$Sex[find_pedigree_row(data$P2[x])])
                {
        
                if(ped$Sex[find_pedigree_row(data$P1[x])] == "M" | ped$Sex[find_pedigree_row(data$P2[x])] == "F")
                    {
                    ped$Father[find_pedigree_row(data$O[x])]=ped$ID[find_pedigree_row(data$P1[x])]
                    ped$Mother[find_pedigree_row(data$O[x])]=ped$ID[find_pedigree_row(data$P2[x])]
                }
                else
                    {
                    ped$Father[find_pedigree_row(data$O[x])]=ped$ID[find_pedigree_row(data$P2[x])]
                    ped$Mother[find_pedigree_row(data$O[x])]=ped$ID[find_pedigree_row(data$P1[x])]
                }
            }
            else if(ped$Sex[find_pedigree_row(data$P1[x])] == "Q" & ped$Sex[find_pedigree_row(data$P2[x])] == "Q")
                {
                ped$ParentSexQ1[find_pedigree_row(data$O[x])]=ped$ID[find_pedigree_row(data$P1[x])]
                ped$ParentSexQ2[find_pedigree_row(data$O[x])]=ped$ID[find_pedigree_row(data$P2[x])]
                #print("Sex unclear:")
                #print(x)
            }
    return(ped)
}

checkBirthdates=function(x)
    {
    if(Birthdate_File != "")
        {
        O=mendels_trio_suspects$O[x]
        P1=mendels_trio_suspects$P1[x]
        P2=mendels_trio_suspects$P2[x]
        if(O %in% Birthdates$ID)
            {
            if(P1 %in% Birthdates$ID & P2 %in% Birthdates$ID)
                {
                if(Birthdates$Birthdate[Birthdates$ID == O] > Birthdates$Birthdate[Birthdates$ID == P1] &
                   Birthdates$Birthdate[Birthdates$ID == O] > Birthdates$Birthdate[Birthdates$ID == P2])
                    {
                    return(TRUE)
                }else{
                    print(paste("Parent(s) younger than offspring. O=",O,", P1=",P1,", P2=",P2,sep=""))
                    return(FALSE)
                }
            }
            if(P1 %in% Birthdates$ID & !(P2 %in% Birthdates$ID))
                {
                if(Birthdates$Birthdate[Birthdates$ID == O] > Birthdates$Birthdate[Birthdates$ID == P1])
                    {
                    return(TRUE)
                }else{
                    print(paste("Parent(s) younger than offspring. O=",O,", P1=",P1,sep=""))
                    return(FALSE)
                }
            }
            if(P2 %in% Birthdates$ID & !(P1 %in% Birthdates$ID))
                {
                if(Birthdates$Birthdate[Birthdates$ID == O] > Birthdates$Birthdate[Birthdates$ID == P2])
                    {
                    return(TRUE)
                }else{
                    print(paste("Parent(s) younger than offspring. O=",O,", P2=",P2,sep=""))
                    return(FALSE)
                }
            
            }
            if(!(P1 %in% Birthdates$ID | P2 %in% Birthdates$ID))
                {
                return(TRUE)
            }
        }else{
            return(TRUE)
        }
    }else{
        return(TRUE)
    }  
}

FindGoodTrios=function(pedigree, Dup_plus_ID_cols)
    {
#     if(sum(is.na(pedigree$Dup1)) == nrow(pedigree))
#     {
#     pedigree=pedigree[,-1*which(colnames(pedigree)=="Dup1")]
#     Dup_plus_ID_cols=c("ID")
# }

    Genomics_Sex=fread(Genomics_Sex_File, data.table=FALSE)

    for(i in 1:nrow(pedigree))
        {
        print(i)
        cur_sex=Genomics_Sex$Genomics_Sex[Genomics_Sex$indv %in% as.list(t(pedigree[i,Dup_plus_ID_cols]))]
        if(length(unique(cur_sex))==1)
            {
            pedigree$Sex[i]=cur_sex[1]
        }
        else if(length(unique(cur_sex[cur_sex!="Q"]))==1)
            {
            pedigree$Sex[i]=cur_sex[cur_sex!="Q"][1]
        }
        else
            {
            pedigree$Sex[i]="Q"
        }

    }

    

    pedigree$Father=NA
    pedigree$Mother=NA
    pedigree$ParentSexQ1 = NA
    pedigree$ParentSexQ2 = NA
    if(nrow(mendels_trio_suspects) > 0)
        {
        
    for(i in 1:nrow(mendels_trio_suspects))
        {
        if(nrow(pedigree[c(find_pedigree_row(mendels_trio_suspects$O[i]),find_pedigree_row(mendels_trio_suspects$P1[i]),find_pedigree_row(mendels_trio_suspects$P2[i])),])==3)
            {
            if(checkBirthdates(i)) # birthdates make sense (or rather, they dont not make sense)
                {
                trio_subset=mendels_trio_suspects[mendels_trio_suspects$O %in% t(pedigree[find_pedigree_row(mendels_trio_suspects$O[i]), Dup_plus_ID_cols]),]
            if((is.na(pedigree$Father[c(find_pedigree_row(mendels_trio_suspects$O[i]))]) & is.na(pedigree$Mother[c(find_pedigree_row(mendels_trio_suspects$O[i]))])  & is.na(pedigree$ParentSexQ1[c(find_pedigree_row(mendels_trio_suspects$O[i]))])  & is.na(pedigree$ParentSexQ2[c(find_pedigree_row(mendels_trio_suspects$O[i]))])))

                {
                if(nrow(trio_subset)==1){
                pedigree=fill_in_ped(pedigree,mendels_trio_suspects,i)
                    }else{
                
                print("Multiple possible parents, choosing best trio.")
                    #criteria=c("IBD0_mean", "IBD0_IQR_mean", "homo_mendel_rel_mean", "V1_percentile")
                    #criteria=c("homo_mendel_rel_mean", "V1_percentile")
                    criteria=c("V1_percentile")
                    if(no_IQR==TRUE)
                        {
                        criteria=c("V1_percentile")
                        }
                    trio_subset=mendels_trio_suspects[mendels_trio_suspects$O %in% t(pedigree[find_pedigree_row(mendels_trio_suspects$O[i]), Dup_plus_ID_cols]),]
                    trio_subset$score=0
                    for(crit in criteria)
                        {
                        trio_subset$score[which(trio_subset[,crit] == min(trio_subset[,crit]))]=trio_subset$score[which(trio_subset[,crit] == min(trio_subset[,crit]))]+1
                    }
                    trio_subset$ChosenTrio=FALSE
                    trio_subset$ChosenTrio[trio_subset$score == max(trio_subset$score, na.rm=TRUE)]=TRUE
                    
                    if(sum(trio_subset$ChosenTrio==TRUE)>1)
                        {
                        trio_subset$ChosenTrio=FALSE
                        
                        }
                    
                    fwrite(trio_subset, file = paste(folder,vcf,"-",pedigree_file_add_name, "-trios-choices.csv", sep=""), col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE, append=TRUE)
                    
                    trio_subset=subset(trio_subset, ChosenTrio == TRUE)
                    if(nrow(trio_subset)>0)
                        {
                        pedigree=fill_in_ped(pedigree,trio_subset,1)
                        print(pedigree[find_pedigree_row(mendels_trio_suspects$O[i]),])
                        }else
                        {
                        print("Unclear trio choice.")
                        pedigree$Father[find_pedigree_row(data$O)]=NA
                        pedigree$Mother[find_pedigree_row(data$O)]=NA
                        }
                }
            }
            else
                {
                if(pedigree$ID[find_pedigree_row(mendels_trio_suspects$P1[i])] %in% c(pedigree$Mother[c(find_pedigree_row(mendels_trio_suspects$O[i]))],pedigree$Father[c(find_pedigree_row(mendels_trio_suspects$O[i]))], pedigree$ParentSexQ1[c(find_pedigree_row(mendels_trio_suspects$O[i]))], pedigree$ParentSexQ2[c(find_pedigree_row(mendels_trio_suspects$O[i]))]) & pedigree$ID[find_pedigree_row(mendels_trio_suspects$P2[i])] %in% c(pedigree$Mother[c(find_pedigree_row(mendels_trio_suspects$O[i]))],pedigree$Father[c(find_pedigree_row(mendels_trio_suspects$O[i]))], pedigree$ParentSexQ1[c(find_pedigree_row(mendels_trio_suspects$O[i]))], pedigree$ParentSexQ2[c(find_pedigree_row(mendels_trio_suspects$O[i]))]))
                    {

                }
                else
                    {
                    # formerly only here using scores
                    
                }

            }
            }
            }
        else
            {
            print("not enough rows:")
            print(i)
        }
        }
    }

    nrow(pedigree[!is.na(pedigree$Father),])
    nrow(pedigree[!is.na(pedigree$Mother),])

    nrow(pedigree[!is.na(pedigree$ParentSexQ1),])
    nrow(pedigree[!is.na(pedigree$ParentSexQ2),])
    
    return(pedigree)
    }

checkBirthdates2=function(O, P)
    {
    if(Birthdate_File != "")
        {
        

        if(O %in% Birthdates$ID)
            {
            if(P %in% Birthdates$ID)
                {
                if(Birthdates$Birthdate[Birthdates$ID == O] > Birthdates$Birthdate[Birthdates$ID == P])
                    {
                    return(TRUE)
                }else{
                    print(paste("Parent younger than offspring. O=",O,", P=",P,sep=""))
                    return(FALSE)
                }
            }else{
            return(TRUE)
                }
        }else{
            return(TRUE)
        }
    }else{
        return(TRUE)
    }  
}

checkBirthdates3=function(PO1, PO2)
    {
    if(Birthdate_File != "")
        {
        

        if(PO1 %in% Birthdates$ID)
            {
            if(PO2 %in% Birthdates$ID)
                {
                if(Birthdates$Birthdate[Birthdates$ID == PO1] != Birthdates$Birthdate[Birthdates$ID == PO2])
                    {
                    return(TRUE)
                }else{
                    print(paste("Both share the same birthdate. PO1=",PO1,", PO2=",PO2,sep=""))
                    return(FALSE)
                }
            }else{
            return(TRUE)
                }
        }else{
            return(TRUE)
        }
    }else{
        return(TRUE)
    }  
}

CheckIfUnknownSinglePoRemain=function(x)
    {
    checkIfTwoIDsHavePedigreedRelationship=df$ID2[which_to_look_at[x]] %in% t(pedigree[find_pedigree_row(df$ID1[which_to_look_at[x]]),])
    checkIfTwoIDsHavePedigreedRelationship=checkIfTwoIDsHavePedigreedRelationship | df$ID1[which_to_look_at[x]] %in% t(pedigree[find_pedigree_row(df$ID2[which_to_look_at[x]]),])
    

        
        
        
        # This still needs a check whether there are actually free boxes to fill them in.
        
        #Are there free boxes?
        
        ID1_empty=(pedigree$Sex[find_pedigree_row(df$ID2[which_to_look_at[x]])]=="M" & is.na(pedigree$Father[find_pedigree_row(df$ID1[which_to_look_at[x]])]) & (is.na(pedigree$ParentSexQ1[find_pedigree_row(df$ID1[which_to_look_at[x]])]) | is.na(pedigree$ParentSexQ2[find_pedigree_row(df$ID1[which_to_look_at[x]])]))) | (pedigree$Sex[find_pedigree_row(df$ID2[which_to_look_at[x]])]=="F" & is.na(pedigree$Mother[find_pedigree_row(df$ID1[which_to_look_at[x]])]) & (is.na(pedigree$ParentSexQ1[find_pedigree_row(df$ID1[which_to_look_at[x]])]) | is.na(pedigree$ParentSexQ2[find_pedigree_row(df$ID1[which_to_look_at[x]])]))) | (pedigree$Sex[find_pedigree_row(df$ID2[which_to_look_at[x]])]=="Q" &  (is.na(pedigree$ParentSexQ1[find_pedigree_row(df$ID1[which_to_look_at[x]])]) | is.na(pedigree$ParentSexQ2[find_pedigree_row(df$ID1[which_to_look_at[x]])])))
        
        ID2_empty=(pedigree$Sex[find_pedigree_row(df$ID1[which_to_look_at[x]])]=="M" & is.na(pedigree$Father[find_pedigree_row(df$ID2[which_to_look_at[x]])]) & (is.na(pedigree$ParentSexQ1[find_pedigree_row(df$ID2[which_to_look_at[x]])]) | is.na(pedigree$ParentSexQ2[find_pedigree_row(df$ID2[which_to_look_at[x]])]))) | (pedigree$Sex[find_pedigree_row(df$ID1[which_to_look_at[x]])]=="F" & is.na(pedigree$Mother[find_pedigree_row(df$ID2[which_to_look_at[x]])]) & (is.na(pedigree$ParentSexQ1[find_pedigree_row(df$ID2[which_to_look_at[x]])]) | is.na(pedigree$ParentSexQ2[find_pedigree_row(df$ID2[which_to_look_at[x]])]))) | (pedigree$Sex[find_pedigree_row(df$ID1[which_to_look_at[x]])]=="Q" &  (is.na(pedigree$ParentSexQ1[find_pedigree_row(df$ID2[which_to_look_at[x]])]) | is.na(pedigree$ParentSexQ2[find_pedigree_row(df$ID2[which_to_look_at[x]])])))
        

        
        
 

    
    
    return((ID1_empty | ID2_empty) & !checkIfTwoIDsHavePedigreedRelationship & checkBirthdates3(df$ID1[which_to_look_at[x]],df$ID2[which_to_look_at[x]]))
                                                                                                                                  
}

FindSingleParent=function()
{
  if(length(which_to_look_at) > 0)
  {
    for(do_again in 1:3)
    {
      if(length(which_to_look_at) > 0)
      {
        for(i in 1:length(which_to_look_at))
        {
          ID1_noP=FALSE
          ID2_noP=FALSE
          all_PO_ID1=""
          all_PO_ID2=""
          if(is.na(pedigree$Father[find_pedigree_row(df$ID1[which_to_look_at[i]])]) & is.na(pedigree$Mother[find_pedigree_row(df$ID1[which_to_look_at[i]])]))
          {
            ID1_noP=TRUE
          }
          if(is.na(pedigree$Father[find_pedigree_row(df$ID2[which_to_look_at[i]])]) & is.na(pedigree$Mother[find_pedigree_row(df$ID2[which_to_look_at[i]])]))
          {
            ID2_noP=TRUE
          }
          #if(ID1_noP == TRUE)
          #{
          all_PO_ID1=unique(c(df$ID1[!is.na(df$suspect) & df$suspect_PO_strict & df$ID2 %in% t(pedigree[find_pedigree_row(df$ID1[which_to_look_at[i]]), Dup_plus_ID_cols])], df$ID2[!is.na(df$suspect) & df$suspect_PO_strict & df$ID1 %in% t(pedigree[find_pedigree_row(df$ID1[which_to_look_at[i]]), Dup_plus_ID_cols])]))
          all_PO_ID1=all_PO_ID1[!(all_PO_ID1 %in% t(pedigree[find_pedigree_row(df$ID1[which_to_look_at[i]]), Dup_plus_ID_cols]))]
          all_PO_ID1=pedigree$ID[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row))]
          all_PO_ID1=all_PO_ID1[!duplicated(all_PO_ID1)]
          all_PO_ID1=all_PO_ID1[!is.na(all_PO_ID1)]
          
          #}
          #if(ID2_noP == TRUE)
          #{
          all_PO_ID2=unique(c(df$ID2[!is.na(df$suspect) & df$suspect_PO_strict & df$ID1 %in% t(pedigree[find_pedigree_row(df$ID2[which_to_look_at[i]]), Dup_plus_ID_cols])], df$ID1[!is.na(df$suspect) & df$suspect_PO_strict & df$ID2 %in% t(pedigree[find_pedigree_row(df$ID2[which_to_look_at[i]]), Dup_plus_ID_cols])]))
          all_PO_ID2=all_PO_ID2[!(all_PO_ID2 %in% t(pedigree[find_pedigree_row(df$ID2[which_to_look_at[i]]), Dup_plus_ID_cols]))]
          all_PO_ID2=pedigree$ID[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row))]
          all_PO_ID2=all_PO_ID2[!duplicated(all_PO_ID2)]
          all_PO_ID2=all_PO_ID2[!is.na(all_PO_ID2)]
          
          #}
          if(ID1_noP == TRUE)
          {
            # check if there is a known reverse relationship already
            all_PO_ID1=all_PO_ID1[!(pedigree$ID[find_pedigree_row(df$ID1[which_to_look_at[i]])] %in% 
                                      t(pedigree[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row)),c("Father")])) & 
                                    !(pedigree$ID[find_pedigree_row(df$ID1[which_to_look_at[i]])] %in% 
                                        t(pedigree[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row)),c("Mother")]))]
            all_PO_ID1=all_PO_ID1[!is.na(all_PO_ID1)]
            if(length(all_PO_ID1) > 0)
            {
              # there has to be at least one known parent on the other side
              all_PO_ID1=all_PO_ID1[!is.na(t(pedigree[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row)),c("Father")])) | !is.na(t(pedigree[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row)),c("Mother")]))]
            }
            
            if(length(all_PO_ID1)==1){
              #if you find one PO individual, and one of the two has a known parent of the same sex, 
              #just put the parent in the remaining slot
              if(checkBirthdates2(O=df$ID1[which_to_look_at[i]], P=all_PO_ID1))
              {
                if(pedigree$Sex[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row))] == "M" & !is.na(t(pedigree[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row)),c("Father")])))
                {
                  pedigree$Father[find_pedigree_row(df$ID1[which_to_look_at[i]])] = pedigree$ID[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row))] 
                }
                if(pedigree$Sex[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row))] == "F" & !is.na(t(pedigree[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row)),c("Mother")])))
                {
                  pedigree$Mother[find_pedigree_row(df$ID1[which_to_look_at[i]])] = pedigree$ID[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row))] 
                }
              }}
            if(length(all_PO_ID1)>1){
              
              
              
              criteria=c("IBD0", "IBD0_IQR", "homo_mendel_rel")
                if(no_IQR==TRUE)
                    {
                    criteria=c("IBD0", "homo_mendel_rel")
                    }
              df_subset_M=df[(df$ID1 == df$ID1[which_to_look_at[i]] & 
                                df$ID2 %in% all_PO_ID1[which(pedigree$Sex[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row))] == "M")]) |
                               (df$ID2 == df$ID1[which_to_look_at[i]] & 
                                  df$ID1 %in% all_PO_ID1[which(pedigree$Sex[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row))] == "M")]),]
              df_subset_F=df[(df$ID1 == df$ID1[which_to_look_at[i]] & 
                                df$ID2 %in% all_PO_ID1[which(pedigree$Sex[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row))] == "F")]) |
                               (df$ID2 == df$ID1[which_to_look_at[i]] & 
                                  df$ID1 %in% all_PO_ID1[which(pedigree$Sex[unlist(lapply(all_PO_ID1, FUN=find_pedigree_row))] == "F")]),]
              
              if(nrow(df_subset_M)>1)
              {
                print("Multiple possible fathers, finding best fit.")
                df_subset_M$score=0
                for(crit in criteria)
                {
                  df_subset_M$score[which(df_subset_M[,crit] == min(df_subset_M[,crit]))]=
                    df_subset_M$score[which(df_subset_M[,crit] == min(df_subset_M[,crit]))]+1
                }
                df_subset_M=subset(df_subset_M, score == max(df_subset_M$score, na.rm=TRUE))
              }
              if(nrow(df_subset_F)>1)
              {
                print("Multiple possible mothers, finding best fit.")
                df_subset_F$score=0
                for(crit in criteria)
                {
                  df_subset_F$score[which(df_subset_F[,crit] == min(df_subset_F[,crit]))]=
                    df_subset_F$score[which(df_subset_F[,crit] == min(df_subset_F[,crit]))]+1
                }
                df_subset_F=subset(df_subset_F, score == max(df_subset_F$score, na.rm=TRUE))
              }
              
              filtered_PO=c(df_subset_M$ID1, df_subset_F$ID1, df_subset_F$ID2, df_subset_M$ID2)
              
              all_PO_ID1=all_PO_ID1[all_PO_ID1 %in% filtered_PO]
              
              if(length(all_PO_ID2)>0)
              {
              for(po_cur in 1:length(all_PO_ID1))
              {
                if(checkBirthdates2(O=df$ID1[which_to_look_at[i]], P=all_PO_ID1[po_cur]))
                {
                  if(pedigree$Sex[unlist(lapply(all_PO_ID1[po_cur], FUN=find_pedigree_row))] == "M" & !is.na(t(pedigree[unlist(lapply(all_PO_ID1[po_cur], FUN=find_pedigree_row)),c("Father")])))
                  {
                    pedigree$Father[find_pedigree_row(df$ID1[which_to_look_at[i]])] = pedigree$ID[unlist(lapply(all_PO_ID1[po_cur], FUN=find_pedigree_row))] 
                  }
                  if(pedigree$Sex[unlist(lapply(all_PO_ID1[po_cur], FUN=find_pedigree_row))] == "F" & !is.na(t(pedigree[unlist(lapply(all_PO_ID1[po_cur], FUN=find_pedigree_row)),c("Mother")])))
                  {
                    pedigree$Mother[find_pedigree_row(df$ID1[which_to_look_at[i]])] = pedigree$ID[unlist(lapply(all_PO_ID1[po_cur], FUN=find_pedigree_row))] 
                  }
                }
              }
              }
            }
          }
          
          if(ID2_noP == TRUE)
          {
            # check if there is a known reverse relationship already
            all_PO_ID2=all_PO_ID2[!(pedigree$ID[find_pedigree_row(df$ID2[which_to_look_at[i]])] %in% 
                                      t(pedigree[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row)),c("Father")])) & 
                                    !(pedigree$ID[find_pedigree_row(df$ID2[which_to_look_at[i]])] %in% 
                                        t(pedigree[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row)),c("Mother")]))]
            all_PO_ID2=all_PO_ID2[!is.na(all_PO_ID2)]
            if(length(all_PO_ID2) > 0)
            {
              # only keep known parents on one side
              all_PO_ID2=all_PO_ID2[!is.na(t(pedigree[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row)),c("Father")])) | !is.na(t(pedigree[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row)),c("Mother")]))]
            }
            
            if(length(all_PO_ID2)==1){
              #if you find one PO individual, and one of the two has known parents, 
              #just put the parent in the remaining slot
              if(checkBirthdates2(O=df$ID2[which_to_look_at[i]], P=all_PO_ID2))
              {
                if(pedigree$Sex[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row))] == "M" & !is.na(t(pedigree[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row)),c("Father")])))
                {
                  pedigree$Father[find_pedigree_row(df$ID2[which_to_look_at[i]])] = pedigree$ID[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row))] 
                }
                if(pedigree$Sex[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row))] == "F" & !is.na(t(pedigree[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row)),c("Mother")])))
                {
                  pedigree$Mother[find_pedigree_row(df$ID2[which_to_look_at[i]])] = pedigree$ID[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row))] 
                }
              }}
            if(length(all_PO_ID2)>1){
              
              print(i)
              
              criteria=c("IBD0", "IBD0_IQR", "homo_mendel_rel")
                if(no_IQR==TRUE)
                    {
                    criteria=c("IBD0", "homo_mendel_rel")
                    }
              df_subset_M=df[(df$ID1 == df$ID2[which_to_look_at[i]] & 
                                df$ID2 %in% all_PO_ID2[which(pedigree$Sex[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row))] == "M")]) |
                               (df$ID2 == df$ID2[which_to_look_at[i]] & 
                                  df$ID1 %in% all_PO_ID2[which(pedigree$Sex[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row))] == "M")]),]
              df_subset_F=df[(df$ID1 == df$ID2[which_to_look_at[i]] & 
                                df$ID2 %in% all_PO_ID2[which(pedigree$Sex[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row))] == "F")]) |
                               (df$ID2 == df$ID2[which_to_look_at[i]] & 
                                  df$ID1 %in% all_PO_ID2[which(pedigree$Sex[unlist(lapply(all_PO_ID2, FUN=find_pedigree_row))] == "F")]),]
              
              if(nrow(df_subset_M)>1)
              {
                print("Multiple possible fathers, finding best fit.")
                df_subset_M$score=0
                for(crit in criteria)
                {
                  df_subset_M$score[which(df_subset_M[,crit] == min(df_subset_M[,crit]))]=
                    df_subset_M$score[which(df_subset_M[,crit] == min(df_subset_M[,crit]))]+1
                }
                df_subset_M=subset(df_subset_M, score == max(df_subset_M$score, na.rm=TRUE))
              }
              if(nrow(df_subset_F)>1)
              {
                print("Multiple possible mothers, finding best fit.")
                df_subset_F$score=0
                for(crit in criteria)
                {
                  df_subset_F$score[which(df_subset_F[,crit] == min(df_subset_F[,crit]))]=
                    df_subset_F$score[which(df_subset_F[,crit] == min(df_subset_F[,crit]))]+1
                }
                df_subset_F=subset(df_subset_F, score == max(df_subset_F$score, na.rm=TRUE))
              }
              
              filtered_PO=c(df_subset_M$ID1, df_subset_F$ID1, df_subset_F$ID2, df_subset_M$ID2)
              
              all_PO_ID2=all_PO_ID2[all_PO_ID2 %in% filtered_PO]
              
              if(length(all_PO_ID2) > 0)
              {
              for(po_cur in 1:length(all_PO_ID2))
              {
                if(checkBirthdates2(O=df$ID2[which_to_look_at[i]], P=all_PO_ID2[po_cur]))
                {
                  if(pedigree$Sex[unlist(lapply(all_PO_ID2[po_cur], FUN=find_pedigree_row))] == "M" & !is.na(t(pedigree[unlist(lapply(all_PO_ID2[po_cur], FUN=find_pedigree_row)),c("Father")])))
                  {
                    pedigree$Father[find_pedigree_row(df$ID2[which_to_look_at[i]])] = pedigree$ID[unlist(lapply(all_PO_ID2[po_cur], FUN=find_pedigree_row))] 
                  }
                  if(pedigree$Sex[unlist(lapply(all_PO_ID2[po_cur], FUN=find_pedigree_row))] == "F" & !is.na(t(pedigree[unlist(lapply(all_PO_ID2[po_cur], FUN=find_pedigree_row)),c("Mother")])))
                  {
                    pedigree$Mother[find_pedigree_row(df$ID2[which_to_look_at[i]])] = pedigree$ID[unlist(lapply(all_PO_ID2[po_cur], FUN=find_pedigree_row))] 
                  }
                }
              }
              }
            }
          }
          
          
          
          
          
        }
      }
    }
  }
  return(pedigree)
}

FindRemainingPOs=function()
    {
    RemainingPOs=which_to_look_at[unlist(lapply(1:length(which_to_look_at), CheckIfUnknownSinglePoRemain))]
    length(RemainingPOs)
    OutSinglePOsUndirected=df[RemainingPOs,c("ID1","ID2")]
    #OutSinglePOsUndirected
    
    

    nrow(pedigree[!is.na(pedigree$Father),])
    nrow(pedigree[!is.na(pedigree$Mother),])
    
    return(OutSinglePOsUndirected)
    }

LastEffortDirectingPOs=function(OutSinglePOsUndirected)
    {
    
    for(focal in unique(c(OutSinglePOsUndirected$ID1, OutSinglePOsUndirected$ID2)))
        {
        
        
        
        
        
        }
    
    }

source(args[1])

SbatchOrInline=1

if(!exists("no_IQR"))
    {
    no_IQR=FALSE
    }

library(readr)
library(stringr)
library(dplyr)
library(tidyr)
library(data.table)
library(MASS)
library(ggplot2)
library(parallel)
library(naturalsort)
library(sys)
options(datatable.fread.datatable=FALSE)

if(file.exists(paste(folder,vcf,"-",pedigree_file_add_name, "-trios-choices.csv", sep=""))){
    system(command=paste("rm -f ",folder,vcf,"-",pedigree_file_add_name, "-trios-choices.csv", sep=""), intern=TRUE)
   } 

options(scipen=999)

setwd(folder)

print(system(command=paste("export TMPDIR=",folder,sep=""), intern=TRUE))

print(folder)
print(vcf)
print(paste("when_PO_error=",when_PO_error,sep=""))

if(Birthdate_File != "")
        {
        Birthdates=fread(Birthdate_File, data.table=FALSE, stringsAsFactors = FALSE)
        Birthdates$Birthdate=as.Date(Birthdates$Birthdate)
        Birthdates=subset(Birthdates, !is.na(Birthdate))
    }

if(file.exists(paste(vcf,".csi",sep="")))
   {
       
    if(file.info(paste(vcf,".csi",sep=""))$mtime[1] < file.info(vcf)$mtime[1])
        {
        system(command=paste("bcftools index -f ",vcf,sep=""))
    }
       }else{
       system(command=paste("bcftools index -f ",vcf,sep=""))
       }
   
   

if(previously_computed != TRUE)
    {
    print("Finding unique chromosomes in VCF.")
    chromosomes=getChrsInFile()
    print("Running TRUFFLE.")
    print(system(command="echo Start && date", intern=TRUE))
    if(SbatchOrInline==2)
        SbatchTruffle()
    if(SbatchOrInline==1)
        InlineTruffle()
    print(system(command="echo Stop && date", intern=TRUE))
    df=fread(paste(vcf,"-truffle.ibd.gz",sep=""),data.table=FALSE)
    
    if(BarnFixPlates==TRUE)
        {
        df=subset(df, !grepl(pattern = "m14", ID1) & !grepl(pattern = "m14", ID2) & !grepl(pattern = "m23", ID1) & !grepl(pattern = "m23", ID2))
    }
    
    }

if(previously_computed != TRUE)
    {
    print("Find unique individuals in VCF.")
    if(SbatchOrInline==2)
        SbatchGetIndividuals()
    if(SbatchOrInline==1)
        InlineGetIndividuals()
        
}

individuals=readIndividuals()

if(previously_computed != TRUE)
    {
    #print("Calculating individual loci counts.")
    #if(SbatchOrInline==2)
        #SbatchGetIndividualLociCounts()
    #if(SbatchOrInline==1)
        #InlineGetIndividualLociCounts()
    #individuals=ReadIndividualLociCounts()
    }

if(previously_computed != TRUE)
    {
    print("Incorporating IBD IQR")
    print(system(command="echo Start && date", intern=TRUE))
    df=IncorporateIBDandIQR(df)
    print(system(command="echo Stop && date", intern=TRUE))
     }

if(previously_computed != TRUE)
    {
    print("Saving progress.")
    fwrite(df, file=paste(vcf,"-truffle.ibd.iqr",sep=""), quote = FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
    df=subset(df, !duplicated(paste(ID1,ID2)))
    system(command=paste("gzip -f ",vcf,"-truffle.ibd.iqr",sep=""), intern = TRUE) 
    }
    
if(previously_computed == TRUE & previously_computed_homozygous_mendel != TRUE)
    {
    print("Reading previously computed file.")
    df=fread(paste(vcf,"-truffle.ibd.iqr.gz",sep=""), data.table=FALSE, stringsAsFactors = FALSE)
    df=subset(df, !duplicated(paste(ID1,ID2)))
    }

if(previously_computed != TRUE | previously_computed_homozygous_mendel != TRUE)
    {
    print("Calculating homozygous errors.")
    print(system(command="echo Start && date", intern=TRUE))
    if(SbatchOrInline==2)
        SbatchHomozygousMendel()
    if(SbatchOrInline==1)
        InlineHomozygousMendel()
    print(system(command="echo Stop && date", intern=TRUE))
    df=ReadHomozygousMendel(df)
    
    print("Saving progress.")
    df=subset(df, !duplicated(paste(ID1,ID2)))
    fwrite(df, file=paste(vcf,"-truffle.ibd.iqr",sep=""), quote = FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
    system(command=paste("gzip -f ",vcf,"-truffle.ibd.iqr",sep=""), intern = TRUE) 

}

if(previously_computed == TRUE & previously_computed_homozygous_mendel == TRUE)
    {
    print("Reading previously computed file.")
    df=fread(paste(vcf,"-truffle.ibd.iqr.gz",sep=""), data.table=FALSE, stringsAsFactors = FALSE)
    df=subset(df, !duplicated(paste(ID1,ID2)))
    }

if(!("IBD0_IQR" %in% colnames(df)))
{
    df$IBD0_IQR=unlist(lapply(1:nrow(df), IBD0_row_IQR))
    print("Saving progress.")
    fwrite(df, file=paste(vcf,"-truffle.ibd.iqr",sep=""), quote = FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
    system(command=paste("gzip -f ",vcf,"-truffle.ibd.iqr",sep=""), intern = TRUE) 
}



#df$IBD0_sd=unlist(lapply(1:nrow(df), IBD0_row_sd))
#df$IBD1_IQR=unlist(lapply(1:nrow(df), IBD1_row_IQR))
#df$IBD2_IQR=unlist(lapply(1:nrow(df), IBD2_row_IQR))

#df=subset(df, IBD0<0.5)
print("Calculating thresholds.")
print(system(command="echo Start && date", intern=TRUE))
if(previously_computed_three_thresholds==FALSE)
    {
    IBD0_threshold=computeAndPlotIBD0threshold()
    IBD0_IQR_threshold=computeAndPlotIBD0IQRthreshold()
    homo_mendel_rel_threshold=computeAndPlothomomendelrelthreshold()
    }else{
    IBD0_threshold=as.numeric(readLines(paste(vcf,"-",pedigree_file_add_name,"-IBD0_threshold.txt",sep="")))
    IBD0_IQR_threshold=as.numeric(readLines(paste(vcf,"-",pedigree_file_add_name,"-IBD0_IQR_threshold.txt",sep="")))
    homo_mendel_rel_threshold=as.numeric(readLines(paste(vcf,"-",pedigree_file_add_name,"-homo_mendel_rel_threshold.txt",sep="")))
    print(IBD0_threshold)
    print(IBD0_IQR_threshold)
    print(homo_mendel_rel_threshold)
    }
print(system(command="echo Stop && date", intern=TRUE))
PlotIBD2DPThreshold()
print("Classiying relationships.")
df=ClassifyRelationshipCandidates(df)

if(previously_computed != TRUE)
    {
    print("Saving progress.")
    fwrite(df, file=paste(vcf,"-truffle.ibd.iqr",sep=""), quote = FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
    system(command=paste("gzip -f ",vcf,"-truffle.ibd.iqr",sep=""), intern = TRUE) 
    }
print("Reading sex file.")
Genomics_Sex=fread(Genomics_Sex_File, data.table=FALSE)


print("Finding Mendelian trio errors.")
print(system(command="echo Start && date", intern=TRUE))
focals=unique(c(df$ID1[df$suspect=="PO" & !is.na(df$suspect)], df$ID2[df$suspect=="PO" & !is.na(df$suspect)]))
#if(SbatchOrInline==2)
    #SbatchMendelTrios(Genomics_Sex, focals, df)
#if(SbatchOrInline==1)
trios_file=InlineMendelTrios(Genomics_Sex, focals, df)

print("Reading Mendelian trio data.")
out=GetMendelTriosData(trios_file)
all_mendel_summary=out[[2]]
df=out[[1]]
print(system(command="echo Stop && date", intern=TRUE))
print("Finding duplicate samples.")
pedigree=SearchForDuplicateSamples()
Dup_plus_ID_cols=c("ID", colnames(pedigree)[which(startsWith(colnames(pedigree), "Dup"))])
for(i in Dup_plus_ID_cols)
    {
    if(sum(!is.na(pedigree[, i]))==0)
        {Dup_plus_ID_cols=Dup_plus_ID_cols[-1*which(Dup_plus_ID_cols==i)]}
}
print("Making Mendelian data frame.")
print(system(command="echo Start && date", intern=TRUE))
which_to_look_at=which(!is.na(df$suspect) & df$suspect=="PO")
PO_test_IDs=unique(c(df$ID1[which_to_look_at],df$ID2[which_to_look_at]))
no_dups=all(Dup_plus_ID_cols=="ID")
#mendels=GatherDeepMendelianData()
mendels=bind_rows(mclapply(PO_test_IDs, FUN=GatherDeepMendelianData, no_dups=no_dups, mc.cores = max_cores))
mendels=subset(mendels,!is.na(V1))
mendels=unique(mendels)
mendels$V1_percentile=ecdf(as.numeric(mendels$V1))((as.numeric(mendels$V1)))
mendels$V1_V1_mean=mendels$V1/(mendels$V1_mean)
print(system(command="echo Stop && date", intern=TRUE))
fwrite(mendels, paste(vcf,"-",pedigree_file_add_name,"-mendels.csv",sep=""))
print("Finding Mendelian thresholds.")
print(system(command="echo Start && date", intern=TRUE))
out=CompareMendelianErrorsForThresholds()
V1_V1_mean_threshold=out[[1]]
V1_percentile_threshold=out[[2]]
print(system(command="echo Stop && date", intern=TRUE))
print("Comparing putative trios.")
print(system(command="echo Start && date", intern=TRUE))
out=DeeperTrioComparison(mendels)
mendels=out
mendels_trio_suspects=subset(mendels, (!is.na(V1_V1_mean)&V1_V1_mean <= V1_V1_mean_threshold) | V1_percentile <= V1_percentile_threshold)
print("Finding good trios.")
pedigree=FindGoodTrios(pedigree,Dup_plus_ID_cols)
print(system(command="echo Stop && date", intern=TRUE))
df$suspect_PO_strict=FALSE
if(no_IQR==TRUE)
    {
    df$suspect_PO_strict[df$IBD0 <= IBD0_threshold & df$homo_mendel_rel <= homo_mendel_rel_threshold]=TRUE
    }else{
    df$suspect_PO_strict[df$IBD0 <= IBD0_threshold & df$IBD0_IQR <= IBD0_IQR_threshold & df$homo_mendel_rel <= homo_mendel_rel_threshold]=TRUE
    }

which_to_look_at=which(df$suspect_PO_strict)
print("Searching for single parents to fill into pedigree.")
print(system(command="echo Start && date", intern=TRUE))
pedigree=FindSingleParent()
print(system(command="echo Stop && date", intern=TRUE))
print("Outputting pedigree and remaining undirected POs.")
OutSinglePOsUndirected=FindRemainingPOs()
    
if(previously_computed!=TRUE | previously_computed_homozygous_mendel != TRUE)
    {
    print("Saving progress.")
    fwrite(df, file=paste(vcf,"-truffle.ibd.iqr",sep=""), quote = FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
    system(command=paste("gzip -f ",vcf,"-truffle.ibd.iqr",sep=""), intern = TRUE) 
    }
    
fwrite(OutSinglePOsUndirected, file = paste(folder,vcf,"-",pedigree_file_add_name, "-PO-undirected.tsv", sep=""), col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)

fwrite(pedigree, file = paste(folder,vcf,"-",pedigree_file_add_name, "-imputed_pedigree.ped", sep=""), col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)

system(command=paste("rm -f ",vcf,pedigree_file_add_name,"_all_mendel.ped",sep=""), intern=TRUE)
    
print("Done.")