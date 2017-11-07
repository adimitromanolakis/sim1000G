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

    vcf <<- vcf

    startSimulation(vcf, totalNumberOfIndividuals = maxNumberOfIndividuals)


}

simulatePopulation = function(num_individuals = 300)

{
    ids = generateUnrelatedIndividuals(num_individuals)
    genotypes = retrieveGenotypes(ids)
}






vcf_file = "region-chr4-93-TMEM156.vcf.gz"
initializeSimulation(vcf_file,min_maf=1e-10, max_maf=0.03)
genotypes = simulatePopulation(100)









examineMAFandLD = function() {

    SIM$reset()


    ids = generateUnrelatedIndividuals(500)
    genotypes = retrieveGenotypes(ids)


    ld_population =   cor   (  t(SIM$population_gt1 + SIM$population_gt2 )   ) ^2

    ld_simulated_data = cor(genotypes)^2
    ld_simulated_data[1:10,1:10]

    s1 = ld_simulated_data[ lower.tri(ld_simulated_data) ]
    s2 = ld_population[ lower.tri(ld_simulated_data) ]

    plot(s1,s2)

    maf1 = apply(genotypes,2,mean)/2
    maf2 = apply( t(vcf$gt1+vcf$gt2)  , 2 ,mean)/2
    maf1[maf1>0.5] = 1-maf1[maf1>0.5]
    maf2[maf2>0.5] = 1-maf2[maf2>0.5]

    plot(maf1,maf2)


    n = nrow(ld_simulated_data)
    ld_simulated_data [ lower.tri(ld_simulated_data)  ] = ld_population[ lower.tri(ld_population) ]


}



