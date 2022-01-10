# SPORE (Specific Parent-Offspring Relationship Estimation)

*SPORE* will take as input a VCF file of a population and will output a pedigree based on the genotypes in the VCF file. However, only parent-offspring relationships will be inferred. The difference of SPORE to other methods is that it works well in very inbred populations. A manuscript detailing the method can be found here: <https://www.biorxiv.org/content/10.1101/2021.11.23.469656v2>.

## Prerequisites

You must be reasonably sure that there is at least one parents-offspring trio among the individuals in the VCF as that *SPORE* assumes this and will deliver false results if this assumption is not met.

You will likely need a Unix-based operating system (Linux) to fulfill the prerequisites for the script. 

### TRUFFLE

Tested with beta 1.38.

Available at <https://adimitromanolakis.github.io/truffle-website/>.

### BCFTOOLS

Tested with 1.10.

Available here: <https://github.com/samtools/bcftools>

### R

The script requires R with some packages installed. The versions with which the script was tested are listed below.

SPORE 1.0 was tested with / written for:

* R 3.6.3
* readr 1.4.0
* stringr 1.4.0
* dplyr 1.0.2
* tidyr 1.1.2
* data.table 1.13.2
* ggplot2 3.3.2
* naturalsort 0.1.3
* sys 3.4
* this.path 0.4.4

## Installation

To use the script, simply make sure you have all the prerequisites mentioned above, clone this repository, and edit the SPORE-Settings.R file to suit your needs.

## Configuration

Please edit the configuration in SPORE-Settings.R and follow the descriptions in there.

## Execution

Execute the script (after configuration) with the command

```` bash
RScript SPORE.R SPORE-Settings.R
````

### Potential issues

* Individual names should not include underscores "_".

## Output files

`*-imputed_pedigree.ped` This is the main result file and includes the pedigree as SPORE was able to infer it. If you input sexes, then "Father" and "Mother" columns are used, otherwise "ParentSexQ1" and "ParentSexQ2" will still contain trio results. If samples were very similar genetically (by default, IBD2 $\geq$ 0.95), then they will be merged into one sample and can be found in the "Dup" columns.

`*-PO-undirected.tsv` This contains a list of relationships where ID1 and ID2 were detected as parent and offspring, but could not be 'directed', meaning that it is unclear who is parent and who is offspring.

`*-IBD0_threshold.png` This shows a plot where on the x axis are different thresholds of IBD0 and the y gives the fraction for three different variables: in blue, how many individuals have at least one relationship with another individual below the IBD0 threshold; in red, how many individuals have more than *when_PO_error* relationships (see configuration file); and in black the distance between those two fractions. The black dot highlights the maximum value of the black line, which is also written on the top left, and was subsequently used as a threshold by SPORE.

`*-IBD0_IQR_threshold.png` Same as above, but for the inter-quartile range of IBD0 between the chromosomes.

`*-homo_mendel_rel_threshold.png` Same as above, but for the fraction of the genome that contains incompatible homozygous genotypes (e.g. 0/0 and 1/1) for two focals.

`*-TrioMendelErrorsComparedToMeanThreshold.png` Similar to before, but here we show on the x a threshold of how a trio's Mendelian errors compare to the mean Mendelian errors of all trios. In blue we show focals with exactly one trio that is below that threshold and in red focals with more than 1 trio below the threshold. The threshold is then picked (blue dot) based on the highest fraction of focals that have exactly one trio.

`*-TrioMendelErrorsPercentileThreshold.png` Same as above, but the x is the percentile of Mendelian errors within all trios that were looked at for a given focal.
