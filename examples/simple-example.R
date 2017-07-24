
library("sim1000G")

vcf = readVCF("haplosims/example1.vcf", maxNumberOfVariants = 222 , min_maf = 0.02 ,max_maf = NA)

#downloadGeneticMap( chromosome  = 4)

readGeneticMap( chromosome = 4)

str(geneticMap)
str(vcf$vcf[,2])

ls(env=vcf)

dim(vcf$vcf)


plotRegionalGeneticMap(vcf$vcf[,2]+1)

startSimulation(vcf, totalNumberOfIndividuals = 1200,randomdata = 0)

hist(generateRecombinationDistances(20000),n=100)


##

SIM$origin1[1:30,1:10]

f = function() {
    SIM$individuals_generated  = 0
    for(i in 1:10) SIM$addUnrelatedIndividual()
}

system.time(f())


SIM$mate(1,2)

