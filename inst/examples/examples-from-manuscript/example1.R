# Example 1 , read a region from 1000 genomes, simulate 300 individuals, compute and display LD patterns

library(sim1000G)
library(gplots)


# Read the example file included in sim1000G

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = sprintf("%s/region.vcf.gz", examples_dir)

# Alternatively provide a vcf file here:
#vcf_file = "~/fs/tmp/sim4/pop1/region-chr4-312-GABRB1.vcf.gz"

vcf = readVCF( vcf_file, maxNumberOfVariants = 200 , min_maf = 0.15 ,max_maf = NA)


startSimulation(vcf, totalNumberOfIndividuals = 8000)
ids = generateUnrelatedIndividuals(20)

genotype = retrieveGenotypes(ids)

rownames(genotype) = sprintf("individual %d",1:nrow(genotype))

heatmap(genotype,col=c("white","orange","red"), Rowv=F)
