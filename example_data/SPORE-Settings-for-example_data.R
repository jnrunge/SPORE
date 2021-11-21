trufflepath="/moto/ziab/users/jr3950/software/TRUFFLE/truffle/truffle" # where is the TRUFFLE executable on your system?
truffle_maf=0.0001 # What is the MAF setting you want to set on TRUFFLE?
truffle_missing=0.95 # What is the missing % setting you want to set on TRUFFLE?

when_PO_error=1 # this is the assumed number of offspring that an individual has. think of it as an average number of offspring you expect individuals to have in this population. This metric essentially changes the sensitivity of the script. 

downsample_for_homozygous_mendel=0.01 # when the script looks for homozygous mendelian errors and loci overlap, what is the downsampling you want to apply to reduce memory load?

pedigree_file_add_name="APO1" # this will add a string to the output filenames, for example to 

min_loci=1 # minimum number of loci that individuals must have to be included in the analysis. For all, set 1.

IBD2_DP_Threshold=0.95 # At what IBD2 % should individuals be treated as duplicate samples? I have used 0.95 in a very inbred population successfully.

folder="./" # folder in which the vcf (below) is located

vcf="simulated_genotypes.vcf.gz" # the vcf you want to analyse

Genomics_Sex_File="" # needs columns "indv" and "Genomics_Sex" (F,M,Q)

Birthdate_File="" # can be empty if you do not have birthdates! - columns "ID" and "Birthdate" (YYYY-MM-DD format)

plots=TRUE # whether you want to produce plots that are saved in the same directory.

max_memory="2G" # Some parts of the script need to know what maximum amount of RAM you have available.
max_cores=1 # Some parts of the script need to know what maximum amount of CPU cores you have available.
trufflecpu=max_cores # Can TRUFFLE use the same amount of cores?


##### SLURM

# not fully implemented! do not use these settings

SbatchOrInline=1 # 1 = Inline; 2 = Sbatch; This allows you to use the inbuild SLURM scheduling. If you are not using SLURM, then you should leave this at 1. SLURM will allow for better paralellization.

truffle_slurm_account="" # which slurm account to schedule tasks in
homozygous_mendel_slurm_account=truffle_slurm_account # which slurm account to schedule tasks in


##### Advanced: Continue previous runs


previously_computed=FALSE # Global setting to tell the script that any part of it was run already
previously_computed_homozygous_mendel=FALSE # Define that the homozygous errors were already computed
previously_computed_mendel=FALSE # Define that the trio Mendelian errors were already computed (Do not do this if you use a new 'when_PO_errorr' value)
previously_computed_three_thresholds=FALSE # Use previously computed IBD0, IQR, and HM thresholds.
