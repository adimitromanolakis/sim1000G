library(sim1000G)

a = commandArgs(T)

nmarker = as.numeric(a[1])
ngen = as.numeric(a[2])


downloadGeneticMap(4)

readGeneticMap(4)


for(nmarker in c(100,200,400,800,1600))
    for(ngen in c(100,200,500,1000,2000,4000,8000))

        {

        vcf  = readVCF("~/fs/tmp/chr4-80-filt.vcf.gz",
                       maxNumberOfVariants = nmarker,
                       min_maf = 0.1, max_maf = NA)

ped = read.table("20130606_g1k.ped",h=T,as=T,sep="\t")

id1 = ped$Individual.ID [ ped$Population %in% c("CEU","TSI") ]
id2 = ped$Individual.ID [ ped$Population == "ASW" ]
id3 = ped$Individual.ID [ ped$Population == "CHS" ]




x<- system.time( startSimulation(vcf, totalNumberOfIndividuals = ngen + 500, subset = id1) )
cat("TIME Init", nmarker, ngen, x,"\n")

saveSimulation("pop1")




SIM$reset()

f = function(N) {
    id = c()
    for(i in 1:N) id[i] = SIM$addUnrelatedIndividual()
    id
}

v = system.time ( id <- f(ngen ) )
genotypes = retrieveGenotypes(id)

cat("TIME Run",nmarker, ngen, v, "\n");

}


