# TE_Ancestral_reconstruction
A repository of scripts that were used to analyze TE genotype calls of the BL6/10 mouse MA lineage in order to recreate ancestral states. 

Briefly, this repository can be used to convert the VCF output file from a MELT run and reconstructs the ancestral states of each TE genotype call with a known phylogeny, and then performs a poisson regression on each branch. The phylogeny of the BL6/10 lines is hard coded into the ancestralTEs.R, and if you want to modify the phylogeny you must copy and paste a new tree into the script.



###### SCRIPTS #######

# filterVCF.py

This script will convert the MELT output VCF file into a TSV that is formatted to be inputted into ancestralTEs.R. This script is hard coded to point to my local directory in order to have it point to your directory simply modify lines 92-96 with your directory. The script will by default only output genotype calls that have at least 1 polymorphic site and are filtered with a split read on either side. 

The output from this script is a TSV that where each column is a strain and each row is position along the genome. Every cell in this table is a 1, or 0 depending on whether there is a presence or absence of a TE, respectively.

# ancestralTEs.R

The primary workhorse of this process. This script requires APE and ggplot2. This script will take in as input a TSV in the format created by filterVCF.py and output a plot with a poisson regression as well as a TSV file with p-values of a deviance of fit test for the poisson regression. Deviance of fit is calculate via a chi-square and is a measure of whether or not we can reject the poisson model as a good explanation of the way the insertions/deletions are gained with respect to branch lengths.

Modifying the directory in the wrapper function will allow you to change which TSVs the script points to to generate the poisson regression. By default this script is set to my local directory and will produce output for B1, B2 and L1Mds. 

# Rmarkdown

The Rmd script is a fairly rough document that contains all of my notes for the creation of the R script in the ancestral reconstruction.
