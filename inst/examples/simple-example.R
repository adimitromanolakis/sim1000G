# Example 1 , read a region from 1000 genomes, simulate 300 individuals, compute and display LD patterns

library("sim1000G")
library(gplots)


vcf = readVCF("./region.vcf.gz", maxNumberOfVariants = 300 , min_maf = 0.12 ,max_maf = NA)

downloadGeneticMap( chromosome  = 4)
readGeneticMap( chromosome = 4)

plotRegionalGeneticMap(vcf$vcf[,2]+1)

startSimulation(vcf, totalNumberOfIndividuals = 1200)

#hist(generateRecombinationDistances(20000),n=100)


## Generatae LD plot

SIM$reset()


id = c()
for(i in 1:330) id[i] = SIM$addUnrelatedIndividual()

# Show haplotype 1  of first 5 individuals
#print(SIM$gt1[1:5,1:6])

# Show haplotype 2
#print(SIM$gt1[1:5,1:6])



genotypes = SIM$gt1[1:320,] + SIM$gt2[1:320,]

print(dim(genotypes))

str(genotypes)


ld_population =   cor   (  t(SIM$population_gt1 + SIM$population_gt2 )   ) ^2
str(ld_population)

ld_simulated_data = cor(genotypes)^2

n = nrow(ld_simulated_data)
ld_simulated_data [ lower.tri(ld_simulated_data)  ] = ld_population[ lower.tri(ld_population) ]



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


