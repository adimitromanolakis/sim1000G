# Example 1 , read a region from 1000 genomes, simulate 300 individuals, compute and display LD patterns

library(sim1000G)
library(gplots)


# Read the example file included in sim1000G

examples_dir = system.file("examples", package = "sim1000G")

vcf_file = file.path(examples_dir,"region.vcf.gz")


# Alternatively provide a vcf file here:
#vcf_file = "~/fs/tmp/sim4/pop1/region-chr4-312-GABRB1.vcf.gz"

vcf = readVCF( vcf_file, maxNumberOfVariants = 200 , min_maf = 0.03 ,max_maf = 0.4)


startSimulation(vcf, totalNumberOfIndividuals = 8000)
ids = generateUnrelatedIndividuals(20)

genotype = retrieveGenotypes(ids)

heatmap(genotype,col=c("white","orange","red"))



str( SIM$gt1[,]  )


#q = SIM$addUnrelatedIndividualWithHaplotype(SIM$gt1[1,], -100)




examineMAFandLD = function() {

    SIM$reset()

    ids = generateUnrelatedIndividuals(2000)
    ids
    genotypes = retrieveGenotypes(ids)

    table(genotypes)
    #image(genotypes<0)

    ids = generateUnrelatedIndividuals(2)
    ids

    x= retrieveGenotypes(ids)
    str(x)
    plot(x[1,])



    ld_population =   cor   (  t(SIM$population_gt1 + SIM$population_gt2 )   ) ^2

    ld_simulated_data = cor(genotypes)^2




    s1 = ld_simulated_data[ lower.tri(ld_simulated_data) ]
    s2 = ld_population[ lower.tri(ld_simulated_data) ]

    plot(s1,s2)

    maf1 = apply(genotypes,2,mean)/2
    maf2 = apply( t(vcf$gt1+vcf$gt2)  , 2 ,mean)/2

    flip = which(maf1>0.5)

  #  maf1[flip] = 1-maf1[flip]
  #  maf2[flip] = 1-maf2[flip]

    plot(maf2,maf1,xlab="True MAF",ylab="Simulated MAF")


}


examineMAFandLD()

















LDplot = function() {

    SIM$reset()


    ids = generateUnrelatedIndividuals(1400)
    genotypes = retrieveGenotypes(ids)


    ld_population =   cor   (  t(SIM$population_gt1 + SIM$population_gt2 )   ) ^2

    ld_simulated_data = cor(genotypes)^2

    n = nrow(ld_simulated_data)
    ld_simulated_data [ lower.tri(ld_simulated_data)  ] = ld_population[ lower.tri(ld_population) ]




    heatmap.2(ld_simulated_data , col=rev(heat.colors(899)) , trace="none",
              Rowv=F,Colv=F,
              #rowsep = 0:10e3, colsep =  0:10e3,
             # sepwidth=c(0.001,0.001),sepcolor=rgb(0,0,0,0.04),
             # offsetRow = 100,offsetCol = 100,
             breaks=seq(0,1,l=900)

              )




}

LDplot()


