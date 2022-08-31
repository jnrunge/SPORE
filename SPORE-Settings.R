trufflepath="truffle" # where is the TRUFFLE executable on your system?
truffle_maf=0.0001 # What is the MAF setting you want to set on TRUFFLE?
truffle_missing=0.95 # What is the missing % setting you want to set on TRUFFLE?

APO=6 # this is the assumed number of offspring that an individual has. think of it as an average number of offspring you expect individuals to have in this population. This metric essentially changes the sensitivity of the script. 

IntermediateSamplingMode=TRUE # Extensive or Intermediate Sampling Mode? If unsure, use Intermediate

downsample_for_homozygous_mendel=0.01 # when the script looks for homozygous mendelian errors and loci overlap, what is the downsampling you want to apply to reduce memory load?

pedigree_file_add_name="" # this will add a string to the output filenames, for example to distinguish different runs with different settings on the same VCF

min_loci=1 # minimum number of loci that individuals must have to be included in the analysis. For all, set 1. This also precludes relationships with less than x overlapping loci from being investigated further

IBD2_DP_Threshold=0.95 # At what IBD2 % should individuals be treated as duplicate samples? I have used 0.95 in a very inbred population successfully.

folder="./" # folder in which the vcf (below) is located

vcf="genotypes.vcf.gz" # the vcf you want to analyse

Genomics_Sex_File="Genomics_Sex.tsv" # needs columns "indv" and "Genomics_Sex" (F,M,Q)

Birthdate_File="" # can be empty if you do not have birthdates! - columns "ID" and "Birthdate" (YYYY-MM-DD format)

plots=TRUE # whether you want to produce plots that are saved in the same directory.

max_memory="12G" # Some parts of the script need to know what maximum amount of RAM you have available, but this is not guarantee to prevent out of memory situations. 

max_cores=1 # Some parts of the script need to know what maximum amount of CPU cores you have available.
trufflecpu=max_cores # Can TRUFFLE use the same amount of cores?


##### Advanced: Change the variables that are used

no_IBD0=FALSE # disable
no_IQR=FALSE # the use of 
no_HM=FALSE # the default variables that SPORE uses

# optional input file, format is ID1,ID2,value,value2,...
ExtraVariables="none" # this should be a path if you want to use this feature


##### Advanced: Continue previous runs


previously_computed=FALSE # Global setting to tell the script that any part of it was run already
previously_computed_homozygous_mendel=FALSE # Define that the homozygous errors were already computed
previously_computed_mendel=FALSE # Define that the trio Mendelian errors were already computed (Do not do this if you use a new 'APO' value)
previously_computed_three_thresholds=FALSE # Use previously computed IBD0, IQR, and HM thresholds.
