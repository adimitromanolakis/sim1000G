# Example 1 , read a region from 1000 genomes, simulate 300 individuals, compute and display LD patterns

library("sim1000G")
library(gplots)

vcf_file_name = "~/fs/tmp/CEU-TSI-ASW-region-chr4-93-TMEM156.vcf.gz"

maxNumberOfVariants = 300
min_maf = 0.02
max_maf = 0.2

vcf = readVCF(vcf_file_name,
              maxNumberOfVariants = maxNumberOfVariants ,
              min_maf = min_maf ,
              max_maf = max_maf)

startSimulation(vcf)


ids = generateUnrelatedIndividuals(150)
genotypes = retrieveGenotypes(ids)




#### Generate multiple datasets ####

datasets = list()

for(i in 1:50) {

    SIM$reset()

    ids = generateUnrelatedIndividuals(150)
    genotypes = retrieveGenotypes(ids)

    datasets[[i]] = genotypes

}








examineMAFandLD = function() {

    SIM$reset()


    ids = generateUnrelatedIndividuals(1500)
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




    # n = nrow(ld_simulated_data)
    # ld_simulated_data [ lower.tri(ld_simulated_data)  ] = ld_population[ lower.tri(ld_population) ]


}



