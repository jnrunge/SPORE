args = commandArgs(trailingOnly=TRUE)

homozygous_mendel_script="./scripts/Mendel1.R"
Fx_Mendel_Script="./scripts/Mendel2.R"
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
    return(as.numeric(IQR(t(df[x,which(startsWith(colnames(df), "IBD0_"))]))))
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
    truffle_tmp[,paste("IBD1_",chr,sep="")]=truffle_tmp$IBD1
    truffle_tmp[,paste("IBD2_",chr,sep="")]=truffle_tmp$IBD2
    truffle_tmp=truffle_tmp[,c("ID1", "ID2", paste("IBD0_",chr,sep=""), paste("IBD1_",chr,sep=""), paste("IBD2_",chr,sep=""))]
    df=left_join(df, truffle_tmp, by = c("ID1", "ID2"))
         
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
    ',homozygous_mendel_script,' ',vcf,'.pureGT.sample ',folder,' ',vcf,'.samples ', downsample_for_homozygous_mendel,' "', vcf,'"', sep="")
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
    
    
    system(command=paste("sh ",vcf,"-homozygous_mendel_command.sh",sep=""), intern=TRUE)
    
    system(command=paste("rm -f ./",vcf,".pureGT",sep=""), intern = TRUE)
    system(command=paste("rm -f ./",vcf,".pureGT.sample",sep=""), intern = TRUE)
    system(command=paste("gzip -f ",vcf,"-ahmm.new_overlap",sep=""), intern = TRUE)
    system(command=paste("gzip -f ",vcf,"-ahmm.mendel_homo_error",sep=""), intern = TRUE)
    
    }




ReadHomozygousMendel=function(df)
    {
    overlapping_loci_done=as.matrix(read.table(paste(folder,vcf,"-ahmm.new_overlap.gz",sep=""), as.is = TRUE, header = FALSE))
     df$overlapping_loci=NA
    for(i in 1:nrow(df))
        {
        df$overlapping_loci[i] = as.numeric(overlapping_loci_done[which(individuals$V1 == df$ID1[i]), which(individuals$V1 == df$ID2[i])])
        df$ID1_loci[i]=as.numeric(overlapping_loci_done[which(individuals$V1 == df$ID1[i]), which(individuals$V1 == df$ID1[i])])
        df$ID2_loci[i]=as.numeric(overlapping_loci_done[which(individuals$V1 == df$ID2[i]), which(individuals$V1 == df$ID2[i])])
    }
    homo_mendel= as.matrix(read.table(paste(folder,vcf,"-ahmm.mendel_homo_error.gz",sep=""), as.is = TRUE, header = FALSE))
    df$homo_mendel=NA
    for(i in 1:nrow(df))
        {
        df$homo_mendel[i] = as.numeric(homo_mendel[which(individuals$V1 == df$ID1[i]), which(individuals$V1 == df$ID2[i])])
    }
    df$overlapping_loci[df$overlapping_loci == 0] = NA
     df$homo_mendel_rel = df$homo_mendel / df$overlapping_loci
    
    
    df$mean_loci=(df$ID1_loci+df$ID2_loci)/2
    
    before=nrow(df)

    df=subset(df, ID1_loci >= min_loci & ID2_loci >= min_loci)
    
    print(paste("Removed ", before-nrow(df), " rows because of too few loci [set in configuration file].",sep=""))
    
    
    return(df)
    }

computeROC=function(x, dataframe, column)
    {
    dataframe=dataframe[x,]
    df_tmp=df_IBD0low
    df_tmp$threshold=FALSE
    df_tmp$threshold[df_tmp[,column] <= dataframe[,column]]=TRUE
    df_tmp_f=filter(df_tmp, threshold == TRUE)
    df_tmp_sum1=summarise(group_by(df_tmp_f, ID1), n1=n(), .groups="drop_last")
    df_tmp_sum2=summarise(group_by(df_tmp_f, ID2), n2=n(), .groups="drop_last")
    df_tmp_sum=full_join(df_tmp_sum1,df_tmp_sum2, by=c("ID1"="ID2"))
    df_tmp_sum$n1[is.na(df_tmp_sum$n1)]=0
    df_tmp_sum$n2[is.na(df_tmp_sum$n2)]=0
    df_tmp_sum$n=df_tmp_sum$n1+df_tmp_sum$n2
    
    dataframe$AtLeastOnePOFocal=sum(df_tmp_sum$n>=1)/length(unique(c(df$ID1,df$ID2)))
    
    dataframe$moreThanTwoPOFocal=sum(df_tmp_sum$n>when_PO_error)/length(unique(c(df$ID1,df$ID2)))
        
    return(dataframe)
}

computeAndPlotIBD0threshold=function()
    {
    sequence=seq(0, quantile(df_IBD0low$IBD0, probs = 1, na.rm=TRUE), by=0.001*max(df_IBD0low$IBD0, na.rm=TRUE))
    if(length(sequence)>10000)
        {
        sequence=sequence[1:10000]
    }

    IBD0_threshold=data.frame(IBD0=sequence, AtLeastOnePOFocal=NA, moreThanTwoPOFocal=NA)

    IBD0_threshold=bind_rows(mclapply(1:nrow(IBD0_threshold), FUN=computeROC, dataframe=IBD0_threshold, column="IBD0", mc.cores = max_cores))

    IBD0_threshold$distance=IBD0_threshold$AtLeastOnePOFocal-IBD0_threshold$moreThanTwoPOFocal

    #IBD0_threshold$distance=predict(loess(distance~IBD0,IBD0_threshold,span=1))


    #ggplot(IBD0_threshold, aes(moreThanTwoPOFocal,AtLeastOnePOFocal,color=IBD0))+
    #geom_line()
    if(plots==TRUE)
        {
        
        p<-ggplot(IBD0_threshold, aes(IBD0,AtLeastOnePOFocal))+
        geom_line(color="blue")+
        geom_line(mapping=aes(IBD0, moreThanTwoPOFocal),color="red")+
        geom_line(mapping=aes(IBD0, distance),color="black")+
        geom_label(x=quantile(IBD0_threshold$IBD0, 0.2), y=quantile(IBD0_threshold$AtLeastOnePOFocal, 0.75), label=IBD0_threshold$IBD0[which(IBD0_threshold$distance==max(IBD0_threshold$distance))[length(which(IBD0_threshold$distance==max(IBD0_threshold$distance)))]])+
        ylab("")+
        geom_point(mapping=aes(x=IBD0_threshold$IBD0[which(IBD0_threshold$distance==max(IBD0_threshold$distance))[length(which(IBD0_threshold$distance==max(IBD0_threshold$distance)))]],
                              y=max(IBD0_threshold$distance)), color="black")
        ggsave(filename = paste(vcf,"-",pedigree_file_add_name,"-IBD0_threshold.png",sep=""), plot = p)
        
        }
    return(IBD0_threshold$IBD0[which(IBD0_threshold$distance==max(IBD0_threshold$distance))[length(which(IBD0_threshold$distance==max(IBD0_threshold$distance)))]])

    }


computeAndPlotIBD0IQRthreshold=function()
    {
    sequence=seq(0, quantile(df_IBD0low$IBD0_IQR, probs = 1, na.rm=TRUE), by=0.001*max(df_IBD0low$IBD0_IQR, na.rm=TRUE))
    if(length(sequence)>10000)
        {
        sequence=sequence[1:10000]
    }

    IBD0_IQR_threshold=data.frame(IBD0_IQR=sequence, AtLeastOnePOFocal=NA, moreThanTwoPOFocal=NA)

    IBD0_IQR_threshold=bind_rows(mclapply(1:nrow(IBD0_IQR_threshold), FUN=computeROC, dataframe=IBD0_IQR_threshold, column="IBD0_IQR", mc.cores = max_cores))

    IBD0_IQR_threshold$distance=IBD0_IQR_threshold$AtLeastOnePOFocal-IBD0_IQR_threshold$moreThanTwoPOFocal

    #IBD0_IQR_threshold$distance=predict(loess(distance~IBD0_IQR,IBD0_IQR_threshold,span=1))


    #ggplot(IBD0_IQR_threshold, aes(moreThanTwoPOFocal,AtLeastOnePOFocal,color=IBD0_IQR))+
    #geom_line()
    if(plots==TRUE)
        {
        p<-ggplot(IBD0_IQR_threshold, aes(IBD0_IQR,AtLeastOnePOFocal))+
        geom_line(color="blue")+
        geom_line(mapping=aes(IBD0_IQR, moreThanTwoPOFocal),color="red")+
        geom_line(mapping=aes(IBD0_IQR, distance),color="black")+
        geom_label(x=quantile(IBD0_IQR_threshold$IBD0_IQR, 0.2), y=quantile(IBD0_IQR_threshold$AtLeastOnePOFocal, 0.75), label=IBD0_IQR_threshold$IBD0_IQR[which(IBD0_IQR_threshold$distance==max(IBD0_IQR_threshold$distance))[length(which(IBD0_IQR_threshold$distance==max(IBD0_IQR_threshold$distance)))]])+
        ylab("")+
        geom_point(mapping=aes(x=IBD0_IQR_threshold$IBD0_IQR[which(IBD0_IQR_threshold$distance==max(IBD0_IQR_threshold$distance))[length(which(IBD0_IQR_threshold$distance==max(IBD0_IQR_threshold$distance)))]],
                              y=max(IBD0_IQR_threshold$distance)), color="black")
        ggsave(filename = paste(vcf,"-",pedigree_file_add_name,"-IBD0_IQR_threshold.png",sep=""), plot = p)
        }
    return(IBD0_IQR_threshold$IBD0_IQR[which(IBD0_IQR_threshold$distance==max(IBD0_IQR_threshold$distance))[length(which(IBD0_IQR_threshold$distance==max(IBD0_IQR_threshold$distance)))]])

    }



computeAndPlothomomendelrelthreshold=function()
    {
    sequence=seq(0, quantile(df_IBD0low$homo_mendel_rel, probs = 1, na.rm=TRUE), by=0.001*max(df_IBD0low$homo_mendel_rel, na.rm=TRUE))
   
    if(length(sequence)>10000)
        {
        sequence=sequence[1:10000]
    }

    homo_mendel_rel_threshold=data.frame(homo_mendel_rel=sequence, AtLeastOnePOFocal=NA, moreThanTwoPOFocal=NA)

    homo_mendel_rel_threshold=bind_rows(mclapply(1:nrow(homo_mendel_rel_threshold), FUN=computeROC, dataframe=homo_mendel_rel_threshold, column="homo_mendel_rel", mc.cores = max_cores))

    homo_mendel_rel_threshold$distance=homo_mendel_rel_threshold$AtLeastOnePOFocal-homo_mendel_rel_threshold$moreThanTwoPOFocal

    #homo_mendel_rel_threshold$distance=predict(loess(distance~homo_mendel_rel,homo_mendel_rel_threshold,span=1))


    #ggplot(homo_mendel_rel_threshold, aes(moreThanTwoPOFocal,AtLeastOnePOFocal,color=homo_mendel_rel))+
    #geom_line()
    if(plots==TRUE)
        {
        p<-ggplot(homo_mendel_rel_threshold, aes(homo_mendel_rel,AtLeastOnePOFocal))+
        geom_line(color="blue")+
        geom_line(mapping=aes(homo_mendel_rel, moreThanTwoPOFocal),color="red")+
        geom_line(mapping=aes(homo_mendel_rel, distance),color="black")+
        geom_label(x=quantile(homo_mendel_rel_threshold$homo_mendel_rel, 0.2), y=quantile(homo_mendel_rel_threshold$AtLeastOnePOFocal, 0.75), label=homo_mendel_rel_threshold$homo_mendel_rel[which(homo_mendel_rel_threshold$distance==max(homo_mendel_rel_threshold$distance))[length(which(homo_mendel_rel_threshold$distance==max(homo_mendel_rel_threshold$distance)))]])+
        ylab("")+
        geom_point(mapping=aes(x=homo_mendel_rel_threshold$homo_mendel_rel[which(homo_mendel_rel_threshold$distance==max(homo_mendel_rel_threshold$distance))[length(which(homo_mendel_rel_threshold$distance==max(homo_mendel_rel_threshold$distance)))]],
                              y=max(homo_mendel_rel_threshold$distance)), color="black")
        ggsave(filename = paste(vcf,"-",pedigree_file_add_name,"-homo_mendel_rel_threshold.png",sep=""), plot = p)
        }
    return(homo_mendel_rel_threshold$homo_mendel_rel[which(homo_mendel_rel_threshold$distance==max(homo_mendel_rel_threshold$distance))[length(which(homo_mendel_rel_threshold$distance==max(homo_mendel_rel_threshold$distance)))]])

    }

PlotIBD2DPThreshold=function()
    {
    
    if(plots==TRUE)
        {
        png('IBD2_DP_Threshold.png')
        hist(df$IBD2, 100)
        abline(v = IBD2_DP_Threshold, col="red")
        dev.off()
        }
    }

ClassifyRelationshipCandidates=function(df)
    {
    df$suspect=NA
    df$suspect[(df$IBD0 <= (IBD0_threshold) | df$homo_mendel_rel <= (homo_mendel_rel_threshold)) & df$IBD0_IQR <= (IBD0_IQR_threshold)]="PO"
    df$suspect[df$IBD2 >= IBD2_DP_Threshold]="DP"

    print(summary(as.factor(df$suspect)))
    
    return(df)
    }

duplicate_Q_comparisons=function(x, df_mendel_print)
    {
    return(paste0(naturalsort(c(df_mendel_print$Father[x], df_mendel_print$Mother[x])), collapse="-"))
}

SbatchMendelTrios=function(Genomics_Sex, focals, df)
    {
    system(command=paste("rm -f ",vcf,when_PO_error,"_all_mendel.ped",sep=""), intern=TRUE)

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
    
    
    write_tsv(x = df_mendel_print, path = paste(vcf,when_PO_error,"_all_mendel.ped",sep=""), append = TRUE)
}
    
    thecommand=paste(" date && cd ",folder," && vcftools --gzvcf ",vcf,
                 " --not-chr X --not-chr Y --not-chr MT --mendel ",vcf,when_PO_error,"_all_mendel.ped -c | cut -f 5 | gzip > ",vcf,when_PO_error,"_all.mendel.gz && cat ",vcf,when_PO_error,"_all_mendel.ped | cut -f 2,3,4 | sed -e 's/\t/_/g' | gzip >> ",vcf,when_PO_error,"_all.mendel.gz && zcat ",vcf,when_PO_error,"_all.mendel.gz | sort -T /moto/ziab/users/jr3950/data/genomes/tmp/ --buffer-size=",as.character(as.numeric(str_replace(max_memory, "G", ""))-1),"g | uniq -c | sort -nr | gzip > ",vcf,when_PO_error,"_all.mendel.summary.gz && rm -f ",vcf,when_PO_error,"_all.mendel.gz",sep="")

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
    system(command=paste("rm -f ",vcf,when_PO_error,"_all_mendel.ped",sep=""), intern=TRUE)
    
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
    
    
    write_tsv(x = df_mendel_print, file = paste(vcf,when_PO_error,"_all_mendel.ped",sep=""), append = TRUE)
}
    
    thecommand=paste("date && cd ",folder," && vcftools --gzvcf ",vcf,
                 " --not-chr X --not-chr Y --not-chr MT --mendel ",vcf,when_PO_error,"_all_mendel.ped -c | cut -f 5 | gzip > ",vcf,when_PO_error,"_all.mendel.gz && cat ",vcf,when_PO_error,"_all_mendel.ped | cut -f 2,3,4 | sed -e 's/\t/_/g' | gzip >> ",vcf,when_PO_error,"_all.mendel.gz && zcat ",vcf,when_PO_error,"_all.mendel.gz | sort -T /moto/ziab/users/jr3950/data/genomes/tmp/ --buffer-size=",as.character(as.numeric(str_replace(max_memory, "G", ""))-1),"g | uniq -c | sort -nr | gzip > ",vcf,when_PO_error,"_all.mendel.summary.gz && rm -f ",vcf,when_PO_error,"_all.mendel.gz",sep="")


    thecommand=gsub("\r?\n|\r", " ", thecommand)

   print(system(command=thecommand, intern=TRUE))
    }

getMendelID=function(x)
        {
        return(strsplit(x, "_", TRUE)[[1]][1])
    }

GetMendelTriosData=function()
    {
    df$mean_mendel_fwd=NA
    df$min_mendel_fwd=NA

    df$mean_mendel_rev=NA
    df$min_mendel_rev=NA

    df$mean_mendel=NA
    df$min_mendel=NA

    df$ID1=as.character(df$ID1)
    df$ID2=as.character(df$ID2)

    all_mendel_summary=fread(paste(vcf,when_PO_error,"_all.mendel.summary.gz",sep=""), data.table=FALSE)

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

    all_mendel_summary$V1=all_mendel_summary$V1-1

    length(unique(all_mendel_summary$ID))

    head(all_mendel_summary)

    date()
    for(cur in unique(all_mendel_summary$ID))
        {
        mendel_summary=subset(all_mendel_summary, ID == cur)
        ID=cur
        #print(ID)
        mendel_summary$V2 = paste(mendel_summary$V2, "_", sep="")
        for(i in which(df$ID1 == ID | df$ID2 == ID))
            {


            if(length(mendel_summary$V1[str_detect(mendel_summary$V2, pattern = paste(df$ID1[i],"_",sep="")) & str_detect(mendel_summary$V2, pattern = paste(df$ID2[i],"_",sep=""))]) > 0)

               {

                    if(length(mendel_summary$V1[startsWith(mendel_summary$V2, prefix=df$ID1[i]) & str_detect(mendel_summary$V2, pattern = paste(df$ID1[i],"_",sep="")) & str_detect(mendel_summary$V2, pattern = paste(df$ID2[i],"_",sep=""))]) > 0)
                        {
                        df$mean_mendel_fwd[i]=mean(mendel_summary$V1[startsWith(mendel_summary$V2, prefix=df$ID1[i]) & str_detect(mendel_summary$V2, pattern = paste(df$ID1[i],"_",sep="")) & str_detect(mendel_summary$V2, pattern = paste(df$ID2[i],"_",sep=""))], na.rm=TRUE)
                    df$min_mendel_fwd[i]=min(mendel_summary$V1[startsWith(mendel_summary$V2, prefix=df$ID1[i]) & str_detect(mendel_summary$V2, pattern = paste(df$ID1[i],"_",sep="")) & str_detect(mendel_summary$V2, pattern = paste(df$ID2[i],"_",sep=""))], na.rm=TRUE)

                    }

                    if(length(mendel_summary$V1[startsWith(mendel_summary$V2, prefix=df$ID2[i]) & str_detect(mendel_summary$V2, pattern = paste(df$ID1[i],"_",sep="")) & str_detect(mendel_summary$V2, pattern = paste(df$ID2[i],"_",sep=""))]) > 0)
                        {

                    df$mean_mendel_rev[i]=mean(mendel_summary$V1[startsWith(mendel_summary$V2, prefix=df$ID2[i]) & str_detect(mendel_summary$V2, pattern = paste(df$ID1[i],"_",sep="")) & str_detect(mendel_summary$V2, pattern = paste(df$ID2[i],"_",sep=""))], na.rm=TRUE)
                    df$min_mendel_rev[i]=min(mendel_summary$V1[startsWith(mendel_summary$V2, prefix=df$ID2[i]) & str_detect(mendel_summary$V2, pattern = paste(df$ID1[i],"_",sep="")) & str_detect(mendel_summary$V2, pattern = paste(df$ID2[i],"_",sep=""))], na.rm=TRUE)

                        }

                    df$min_mendel[i]=min(c(df$min_mendel_fwd[i],df$min_mendel_rev[i]), na.rm=TRUE)
                    df$mean_mendel[i]=mean(c(df$mean_mendel_fwd[i],df$mean_mendel_rev[i]), na.rm=TRUE)

            }



            }
    }
    date()

    df$min_mendel_rel[!is.na(df$min_mendel) & !is.na(df$overlapping_loci)] = df$min_mendel[!is.na(df$min_mendel) & !is.na(df$overlapping_loci)]/df$overlapping_loci[!is.na(df$min_mendel) & !is.na(df$overlapping_loci)]
    df$min_mendel_rel[df$min_mendel_rel > 1] = NA
    summary(df$min_mendel_rel[!is.na(df$min_mendel) & !is.na(df$overlapping_loci)])
    
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

GatherDeepMendelianData=function()
    {
    # The following code can be drastically simplified, but is pretty nasty for historic reasons.

    which_to_look_at=which(!is.na(df$suspect) & df$suspect=="PO")
    mendels=data.frame(V1=NA,O=NA,P1=NA,P2=NA,V1_mean=NA)
    #print(which_to_look_at)
    #which_to_look_at=which(!is.na(df$suspect) & df$suspect == "PO" & (df$ID1 == "m25_39" | df$ID2 == "m25_39"))
    for(row_raw in 1:length(which_to_look_at))
        {
        row=which_to_look_at[row_raw]
        for(row_col in c("ID1", "ID2"))
            {
            PO_test_ID=df[row,row_col]
            PO_test_ID=as.character(PO_test_ID)

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

            cur_focal=cur_focal[!is.na(cur_focal)]        
            cur_subset=subset(df, ID1 %in% cur_focal | ID2 %in% cur_focal)
            cur_subset=subset(cur_subset, suspect=="PO")
            cur_POs=unique(c(cur_subset$ID1, cur_subset$ID2))
            cur_POs=cur_POs[!(cur_POs %in% cur_focal)]
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
            cur_POs=cur_POs[!is.na(cur_POs)]
            cur_POs=unique(cur_POs)
            cur_POs

            cur_IDs=c(cur_focal, cur_POs)
            #cur_IDs=c(str_replace(cur_IDs, "_", "?"), str_replace(cur_IDs, "_", "?0"))
            cur_IDs
            #print(cur_IDs)
            cur_mendels_files=NA
            for(i in cur_IDs)
                {
                if(i == cur_IDs[1])
                     {
                    cur_mendels=subset(all_mendel_summary, ID==i)
                    }
                else
                    {
                    cur_mendels=bind_rows(cur_mendels,subset(all_mendel_summary, ID==i))
                    }
                }
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
                    cur_mendels=subset(cur_mendels, V2 != "FAMILY")
                    cur_mendels=tidyr::separate(cur_mendels, V2, c("O", "P1", "P2"), "_")
        #             cur_mendels$O=str_replace(cur_mendels$O, "SM_", "")
        #             cur_mendels$O=str_replace(cur_mendels$O, "-", "_")
        #             cur_mendels$O=str_replace(cur_mendels$O, "_0", "_")

        #             cur_mendels$P1=str_replace(cur_mendels$P1, "-", "_")
        #             cur_mendels$P1=str_replace(cur_mendels$P1, "_0", "_")
        #             cur_mendels$P2=str_replace(cur_mendels$P2, "-", "_")
        #             cur_mendels$P2=str_replace(cur_mendels$P2, "_0", "_")
                    #head(cur_mendels)

                    #cur_mendels$V1_percentile=ecdf(as.numeric(cur_mendels$V1))((as.numeric(cur_mendels$V1)))
                    cur_mendels$V1_mean=mean((as.numeric(cur_mendels$V1)))


                    #print(cur_focal)
                    #print(cur_POs)

                    cur_mendels=subset(cur_mendels, O %in% cur_focal & P1 %in% cur_POs & P2 %in% cur_POs)

                }


               # }

                mendels=rbind(mendels, cur_mendels[,which(colnames(cur_mendels)=="ID")*-1])
                remove(cur_mendels)
               }

        }


    }
    
    #print(mendels)
    mendels=subset(mendels,!is.na(V1))
    mendels=unique(mendels)
    
    mendels$V1_percentile=ecdf(as.numeric(mendels$V1))((as.numeric(mendels$V1)))
    mendels$V1_V1_mean=mendels$V1/(mendels$V1_mean+1)
    return(mendels)
    
    }

computeROC2=function(x, dataframe)
    {
    
    dataframe=dataframe[x,]
    
    df_analyze_V1_cutoffs_tmp=summarise(group_by(filter(mendels, V1_V1_mean <= dataframe$V1_V1_mean & 
                                                       V1_percentile <= dataframe$V1_percentile),
                                                O), n=n(), .groups="drop_last")
    dataframe$focals_with_one_potential_trio=sum(df_analyze_V1_cutoffs_tmp$n == 1)
    dataframe$focals_with_more_than_one_potential_trio=sum(df_analyze_V1_cutoffs_tmp$n > 1)
        
    return(dataframe)
}

CompareMendelianErrorsForThresholds=function()
    {
    V1_V1_means=seq(min(mendels$V1_V1_mean), quantile(mendels$V1_V1_mean,0.5), by = 0.001)
    V1_percentiles=seq(min(mendels$V1_percentile), quantile(mendels$V1_percentile, 0.5), by = 0.001)
    if(length(V1_V1_means) > 10000)
        {
        V1_V1_means=V1_V1_means[1:10000]
    }
    if(length(V1_percentiles) > 10000)
        {
        V1_percentiles=V1_percentiles[1:10000]
    }

    # How many possible trios per threshold?

    df_analyze_V1_cutoffs=expand.grid(V1_V1_mean=V1_V1_means, V1_percentile=V1_percentiles, focals_with_one_potential_trio=NA, focals_with_more_than_one_potential_trio=NA)

    df_analyze_V1_cutoffs=subset(df_analyze_V1_cutoffs, V1_V1_mean==max(df_analyze_V1_cutoffs$V1_V1_mean) | V1_percentile==max(df_analyze_V1_cutoffs$V1_percentile))

    df_analyze_V1_cutoffs=bind_rows(mclapply(1:nrow(df_analyze_V1_cutoffs), dataframe=df_analyze_V1_cutoffs, FUN=computeROC2, mc.cores=max_cores))

    # for(i in 1:nrow(df_analyze_V1_cutoffs))
    #     {
    #     df_analyze_V1_cutoffs_tmp=summarise(group_by(filter(mendels, V1_V1_mean <= df_analyze_V1_cutoffs$V1_V1_mean[i] & 
    #                                                        V1_percentile <= df_analyze_V1_cutoffs$V1_percentile[i]),
    #                                                 O), n=n())
    #     df_analyze_V1_cutoffs$focals_with_one_potential_trio[i]=sum(df_analyze_V1_cutoffs_tmp$n == 1)
    #     df_analyze_V1_cutoffs$focals_with_more_than_one_potential_trio[i]=sum(df_analyze_V1_cutoffs_tmp$n > 1)
    # }
    df_analyze_V1_cutoffs$focals_with_one_potential_trio=df_analyze_V1_cutoffs$focals_with_one_potential_trio/length(unique(mendels$O))
    df_analyze_V1_cutoffs$focals_with_more_than_one_potential_trio=df_analyze_V1_cutoffs$focals_with_more_than_one_potential_trio/length(unique(mendels$O))
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

IBD0_row_mean=function(x)
    {
    return(as.numeric(mean(t(mendels[x,which(endsWith(colnames(mendels), "_IBD0"))]), na.rm = TRUE)))
}

IBD0_IQR_row_mean=function(x)
    {
    return(as.numeric(mean(t(mendels[x,which(endsWith(colnames(mendels), "_IBD0_IQR"))]), na.rm = TRUE)))
}

homo_mendel_rel_row_mean=function(x)
    {
    return(as.numeric(mean(t(mendels[x,which(endsWith(colnames(mendels), "_homo_mendel_rel"))]), na.rm = TRUE)))
}


DeeperTrioComparison=function()
    {
    PO_columns=c("IBD0","IBD0_IQR","homo_mendel_rel")

    mendels=left_join(mendels, df[,c("ID1","ID2", PO_columns)], by=c("O"="ID1", "P1"="ID2"))
    colnames(mendels)[which(colnames(mendels) %in% PO_columns)]=paste("Nr1_",colnames(mendels)[which(colnames(mendels) %in% PO_columns)],sep="")
    mendels=left_join(mendels, df[,c("ID1","ID2", PO_columns)], by=c("O"="ID2", "P1"="ID1"))
    colnames(mendels)[which(colnames(mendels) %in% PO_columns)]=paste("Nr2_",colnames(mendels)[which(colnames(mendels) %in% PO_columns)],sep="")
    mendels=left_join(mendels, df[,c("ID1","ID2", PO_columns)], by=c("O"="ID1", "P2"="ID2"))
    colnames(mendels)[which(colnames(mendels) %in% PO_columns)]=paste("Nr3_",colnames(mendels)[which(colnames(mendels) %in% PO_columns)],sep="")
    mendels=left_join(mendels, df[,c("ID1","ID2", PO_columns)], by=c("O"="ID2", "P2"="ID1"))
    colnames(mendels)[which(colnames(mendels) %in% PO_columns)]=paste("Nr4_",colnames(mendels)[which(colnames(mendels) %in% PO_columns)],sep="")

    mendels$IBD0_mean=unlist(lapply(1:nrow(mendels), IBD0_row_mean))
    mendels$IBD0_IQR_mean=unlist(lapply(1:nrow(mendels), IBD0_IQR_row_mean))
    mendels$homo_mendel_rel_mean=unlist(lapply(1:nrow(mendels), homo_mendel_rel_row_mean))

    mendels_trio_suspects=subset(mendels, V1_V1_mean <= V1_V1_mean_threshold | V1_percentile <= V1_percentile_threshold)

    return(list(mendels,mendels_trio_suspects))
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

FindGoodTrios=function()
    {
    if(sum(is.na(pedigree$Dup1)) == nrow(pedigree))
    {
    pedigree=pedigree[,-1*which(colnames(pedigree)=="Dup1")]
}

    Genomics_Sex=fread(Genomics_Sex_File, data.table=FALSE)

    for(i in 1:nrow(pedigree))
        {

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

            if(is.na(pedigree$Father[c(find_pedigree_row(mendels_trio_suspects$O[i]))]) & is.na(pedigree$Mother[c(find_pedigree_row(mendels_trio_suspects$O[i]))])  & is.na(pedigree$ParentSexQ1[c(find_pedigree_row(mendels_trio_suspects$O[i]))])  & is.na(pedigree$ParentSexQ2[c(find_pedigree_row(mendels_trio_suspects$O[i]))]))

                {
                pedigree=fill_in_ped(pedigree,mendels_trio_suspects,i)
            }
            else
                {
                if(pedigree$ID[find_pedigree_row(mendels_trio_suspects$P1[i])] %in% c(pedigree$Mother[c(find_pedigree_row(mendels_trio_suspects$O[i]))],pedigree$Father[c(find_pedigree_row(mendels_trio_suspects$O[i]))], pedigree$ParentSexQ1[c(find_pedigree_row(mendels_trio_suspects$O[i]))], pedigree$ParentSexQ2[c(find_pedigree_row(mendels_trio_suspects$O[i]))]) & pedigree$ID[find_pedigree_row(mendels_trio_suspects$P2[i])] %in% c(pedigree$Mother[c(find_pedigree_row(mendels_trio_suspects$O[i]))],pedigree$Father[c(find_pedigree_row(mendels_trio_suspects$O[i]))], pedigree$ParentSexQ1[c(find_pedigree_row(mendels_trio_suspects$O[i]))], pedigree$ParentSexQ2[c(find_pedigree_row(mendels_trio_suspects$O[i]))]))
                    {

                }
                else
                    {
                    print("Already has different parents, choosing best trio.")
                    criteria=c("IBD0_mean", "IBD0_IQR_mean", "homo_mendel_rel_mean", "V1_percentile")
                    trio_subset=mendels_trio_suspects[mendels_trio_suspects$O %in% t(pedigree[find_pedigree_row(mendels_trio_suspects$O[i]), Dup_plus_ID_cols]),]
                    trio_subset$score=0
                    for(crit in criteria)
                        {
                        trio_subset$score[which(trio_subset[,crit] == min(trio_subset[,crit]))]=trio_subset$score[which(trio_subset[,crit] == min(trio_subset[,crit]))]+1
                    }
                    trio_subset=subset(trio_subset, score == max(trio_subset$score, na.rm=TRUE))
                    pedigree=fill_in_ped(pedigree,trio_subset,1)
                    print(pedigree[find_pedigree_row(mendels_trio_suspects$O[i]),])
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

CheckIfUnknownSinglePoRemain=function(x)
    {
    checkIfTwoIDsHavePedigreedRelationship=df$ID2[which_to_look_at[x]] %in% t(pedigree[find_pedigree_row(df$ID1[which_to_look_at[x]]),])
    checkIfTwoIDsHavePedigreedRelationship=checkIfTwoIDsHavePedigreedRelationship | df$ID1[which_to_look_at[x]] %in% t(pedigree[find_pedigree_row(df$ID2[which_to_look_at[x]]),])
    return(checkIfTwoIDsHavePedigreedRelationship)
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
    return(pedigree)
}

FindRemainingPOs=function()
    {
    RemainingPOs=which_to_look_at[!unlist(lapply(1:length(which_to_look_at), CheckIfUnknownSinglePoRemain))]
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

options(scipen=999)

setwd(folder)

print(folder)
print(vcf)
print(paste("when_PO_error=",when_PO_error,sep=""))

if(Birthdate_File != "")
        {
        Birthdates=fread(Birthdate_File, data.table=FALSE, stringsAsFactors = FALSE)
        Birthdates$Birthdate=as.Date(Birthdates$Birthdate)
        Birthdates=subset(Birthdates, !is.na(Birthdate))
    }

#system(command=paste("bcftools index -f ",vcf,sep=""))

if(previously_computed != TRUE)
    {
    print("Finding unique chromosomes in VCF.")
    chromosomes=getChrsInFile()
    print("Running TRUFFLE.")
    if(SbatchOrInline==2)
        SbatchTruffle()
    if(SbatchOrInline==1)
        InlineTruffle()
    
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
    print("Incorporating IBD IQR.")
    df=IncorporateIBDandIQR(df)
     }

if(previously_computed != TRUE)
    {
    print("Saving progress.")
    fwrite(df, file=paste(vcf,"-truffle.ibd.iqr",sep=""), quote = FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
    system(command=paste("gzip -f ",vcf,"-truffle.ibd.iqr",sep=""), intern = TRUE) 
    }
    
if(previously_computed == TRUE & previously_computed_homozygous_mendel != TRUE)
    {
    print("Reading previously computed file.")
    df=fread(paste(vcf,"-truffle.ibd.iqr.gz",sep=""), data.table=FALSE, stringsAsFactors = FALSE)
    
    }

if(previously_computed != TRUE | previously_computed_homozygous_mendel != TRUE)
    {
    print("Calculating homozygous errors.")
    if(SbatchOrInline==2)
        SbatchHomozygousMendel()
    if(SbatchOrInline==1)
        InlineHomozygousMendel()
    df=ReadHomozygousMendel(df)
    
    print("Saving progress.")
    fwrite(df, file=paste(vcf,"-truffle.ibd.iqr",sep=""), quote = FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
    system(command=paste("gzip -f ",vcf,"-truffle.ibd.iqr",sep=""), intern = TRUE) 

}

if(previously_computed == TRUE & previously_computed_homozygous_mendel == TRUE)
    {
    print("Reading previously computed file.")
    df=fread(paste(vcf,"-truffle.ibd.iqr.gz",sep=""), data.table=FALSE, stringsAsFactors = FALSE)
    
    }  


df$IBD0_IQR=unlist(lapply(1:nrow(df), IBD0_row_IQR))
#df$IBD1_IQR=unlist(lapply(1:nrow(df), IBD1_row_IQR))
#df$IBD2_IQR=unlist(lapply(1:nrow(df), IBD2_row_IQR))

df_IBD0low=subset(df, IBD0<0.5)
print("Calculating thresholds.")
IBD0_threshold=computeAndPlotIBD0threshold()
IBD0_IQR_threshold=computeAndPlotIBD0IQRthreshold()
homo_mendel_rel_threshold=computeAndPlothomomendelrelthreshold()

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

if(previously_computed_mendel != TRUE)
    {
    print("Finding Mendelian trio errors.")
    focals=unique(c(df$ID1[df$suspect=="PO" & !is.na(df$suspect)], df$ID2[df$suspect=="PO" & !is.na(df$suspect)]))
    if(SbatchOrInline==2)
        SbatchMendelTrios(Genomics_Sex, focals, df)
    if(SbatchOrInline==1)
        InlineMendelTrios(Genomics_Sex, focals, df)
    }
print("Reading Mendelian trio data.")
out=GetMendelTriosData()
all_mendel_summary=out[[2]]
df=out[[1]]
print("Finding duplicate samples.")
pedigree=SearchForDuplicateSamples()
Dup_plus_ID_cols=c("ID", colnames(pedigree[,which(startsWith(colnames(pedigree), "Dup"))]))
print("Making Mendelian data frame.")
mendels=GatherDeepMendelianData()
print("Finding Mendelian thresholds.")
out=CompareMendelianErrorsForThresholds()
V1_V1_mean_threshold=out[[1]]
V1_percentile_threshold=out[[2]]
print("Comparing putative trios.")
out=DeeperTrioComparison()
mendels=out[[1]]
mendels_trio_suspects=out[[2]]
print("Finding good trios.")
pedigree=FindGoodTrios()

df$suspect_PO_strict=FALSE
df$suspect_PO_strict[df$IBD0 <= IBD0_threshold & df$IBD0_IQR <= IBD0_IQR_threshold & df$homo_mendel_rel <= homo_mendel_rel_threshold]=TRUE

which_to_look_at=which(df$suspect_PO_strict)
print("Searching for single parents to fill into pedigree.")
pedigree=FindSingleParent()
print("Outputting pedigree and remaining undirected POs.")
OutSinglePOsUndirected=FindRemainingPOs()
    
print("Saving progress.")
fwrite(df, file=paste(vcf,"-truffle.ibd.iqr",sep=""), quote = FALSE, col.names=TRUE, row.names=FALSE, sep="\t")
system(command=paste("gzip -f ",vcf,"-truffle.ibd.iqr",sep=""), intern = TRUE) 

fwrite(OutSinglePOsUndirected, file = paste(folder,vcf,"-",pedigree_file_add_name, "-PO-undirected.tsv", sep=""), col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)

fwrite(pedigree, file = paste(folder,vcf,"-",pedigree_file_add_name, "-imputed_pedigree.ped", sep=""), col.names = TRUE, row.names = FALSE, sep="\t", quote=FALSE)
    
print("Done.")