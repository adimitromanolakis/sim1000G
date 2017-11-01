# Example 1 , read a region from 1000 genomes, simulate 300 individuals, compute and display LD patterns

library("sim1000G")
library(gplots)

downloadGeneticMap( chromosome  = 4)
readGeneticMap( chromosome = 4)


initializeSimulation  = function(vcf_file_name = NA,
                                 min_maf = 1e-10,
                                 max_maf = 0.02,
                                 maxNumberOfVariants = 300,
                                 maxNumberOfIndividuals = 2000)
{


    vcf = readVCF(vcf_file_name,
                  maxNumberOfVariants = maxNumberOfVariants ,
                  min_maf = min_maf ,
                  max_maf = max_maf)


    startSimulation(vcf, totalNumberOfIndividuals = maxNumberOfIndividuals)
}

simulateSinglePopulation = function(num_individuals = 300)

{
    # SIM$reset()

    ids = generateUnrelatedIndividuals(100)
    genotypes = retrieveGenotypes(ids)

}





vcf_file = "~/fs/tmp/data-chr4/region-chr4-84-CD38.vcf.gz"
initializeSimulation(vcf_file,max_maf=0.3)
genotypes = ( simulateSinglePopulation(20) )




save(vcf,file="/tmp/1.rdata")

vcf = readVCF("/tmp/1.rdata")





plotRegionalGeneticMap(vcf$vcf[,2]+1)

startSimulation(vcf, totalNumberOfIndividuals = 1200)

#hist(generateRecombinationDistances(20000),n=100)


## Generatae LD plot

SIM$reset()


ids = generateUnrelatedIndividuals(100)
genotypes = retrieveGenotypes(ids)


ld_population =   cor   (  t(SIM$population_gt1 + SIM$population_gt2 )   ) ^2

ld_simulated_data = cor(genotypes)^2
ld_simulated_data[1:10,1:10]

n = nrow(ld_simulated_data)
ld_simulated_data [ lower.tri(ld_simulated_data)  ] = ld_population[ lower.tri(ld_population) ]


library(gplots)


#pdf(file="/tmp/2.pdf",w=15,h=11)
heatmap.2(ld_simulated_data , col=rev(heat.colors(899)) , trace="none",Rowv=F,Colv=F,
          rowsep = 0:10e3, colsep =  0:10e3,           sepwidth=c(0.001,0.001),sepcolor=rgb(0,0,0,0.04),
          offsetRow = 100,offsetCol = 100, breaks=seq(0,1,l=900),

          row.names = "")
#dev.off()

heatmap.2(ld_population , col=rev(heat.colors(899)) , trace="none",Rowv=F,Colv=F,
          rowsep = 0:10e3, colsep =  0:10e3,           sepwidth=c(0.001,0.001),sepcolor=rgb(0,0,0,0.04),
          offsetRow = 100,offsetCol = 100, breaks=seq(0,1,l=900),

          row.names = "")





##


