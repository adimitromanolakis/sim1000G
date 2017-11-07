# Example 1 , read a region from 1000 genomes, simulate 300 individuals, compute and display LD patterns

library("sim1000G")
library(gplots)

#downloadGeneticMap( chromosome  = 4)
#readGeneticMap( chromosome = 4)


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




initializeSimulationWithStratification  =
                    function(vcf_file_name = NA,
                        subset1 = NA,
                        subset2 = NA,
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

    startSimulation(vcf, totalNumberOfIndividuals = maxNumberOfIndividuals, subset = subset1)
    saveSimulation("pop1")

    startSimulation(vcf, totalNumberOfIndividuals = maxNumberOfIndividuals, subset = subset2)
    saveSimulation("pop2")


}





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




}




vcf_file = "~/fs/tmp/data-chr4/region-chr4-93-TMEM156.vcf.gz"
vcf_file = "~/fs/tmp/CEU-TSI-GBR-region-chr4-34-ABLIM2.vcf.gz"

initializeSimulation(vcf_file,min_maf=1e-10, max_maf=0.03)
genotypes = simulatePopulation(100)


# for(i in 1:100) {genotypes = simulatePopulation(100) }



z = vcf$vcf
str(z)
str(z[,5])
varid = paste(z[,1],z[,2],z[,3],z[,4],z[,5] )




save(vcf,file="/tmp/1.rdata")
vcf = readVCF("/tmp/1.rdata")





ped = read.table("20130606_g1k.ped",h=T,as=T,sep="\t")
pops = split(ped$Individual.ID,list(ped$Population))

id_firstPopulation = c( pops$CEU, pops$TSI )
id_secondPopulation = c( pops$ASW )


initializeSimulationWithStratification (
   vcf_file,
   min_maf = 1e-10,
   max_maf = 1,
   maxNumberOfVariants = 200,
   subset1 = id_firstPopulation,
   subset2 = id_secondPopulation,
   maxNumberOfIndividuals = 3000
)


loadSimulation("pop1")
gt1 = simulatePopulation(100)

loadSimulation("pop2")
gt2 = simulatePopulation(30)

x=gt1





loadSimulation("pop1")

maf1 = apply(simulatePopulation(200) , 2 , mean)/2
maf1[maf1>0.5] = 1 - maf1[maf1>0.5]

loadSimulation("pop2")

maf2 = apply(simulatePopulation(200) , 2 , mean)/2
maf2[maf2>0.5] = 1 - maf2[maf2>0.5]

plot(maf1,maf2)


