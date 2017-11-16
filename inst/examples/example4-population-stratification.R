# Example 1 , read a region from 1000 genomes, simulate 300 individuals, compute and display LD patterns

library("sim1000G")
library(gplots)

simulatePopulation = function(num_individuals = 300)

{
    ids = generateUnrelatedIndividuals(num_individuals)
    genotypes = retrieveGenotypes(ids)
}




vcf1 = readVCF("~/fs/tmp/CEU-TSI-GBR/CEU-TSI-GBR-region-chr4-205-MAML3.vcf.gz",
               maxNumberOfVariants = 5000,  min_maf =  1e-6, max_maf = 0.02)


vcf2 = readVCF("~/fs/tmp/AFR/ASW-LWK-YRI-region-chr4-205-MAML3.vcf.gz",
               maxNumberOfVariants = 5000,  min_maf = 1e-6, max_maf = 0.05)

common = intersect(vcf1$varid,vcf2$varid)

vcf1 = subsetVCF(vcf1, common)
vcf2 = subsetVCF(vcf2, common)

table(vcf1$varid == vcf2$varid)


startSimulation(vcf1, totalNumberOfIndividuals = 2000)
saveSimulation("pop1")

startSimulation(vcf2, totalNumberOfIndividuals = 2000)
saveSimulation("pop2")






ped = read.table("20130606_g1k.ped",h=T,as=T,sep="\t")
pops = split(ped$Individual.ID,list(ped$Population))

loadSimulation("pop1")
gt1 = simulatePopulation(100)

loadSimulation("pop2")
gt2 = simulatePopulation(30)






loadSimulation("pop1")

maf1 = apply(simulatePopulation(200) , 2 , mean)/2
maf1[maf1>0.5] = 1 - maf1[maf1>0.5]

loadSimulation("pop2")

maf2 = apply(simulatePopulation(200) , 2 , mean)/2
maf2[maf2>0.5] = 1 - maf2[maf2>0.5]

plot(maf1,maf2)









v1 = list.files(path = "~/fs/tmp/ASW/")
id = sub(".*chr4-","",v1)
id = sub(".vcf.gz","",id)
id

v2 = list.files(path = "~/fs/tmp/CEU-TSI-GBR//")

v1

data.frame(v1,v2)


readAllVCF = function() {
    v1 = list.files(path = "~/fs/tmp/AFR/", full.names = TRUE)
    v2 = list.files(path = "~/fs/tmp/CEU-TSI-GBR/", full.names = TRUE)



    id = sub(".*chr4-","",v1)
    id = sub(".vcf.gz","",id)
    id


    vcf_pop1 = list()
    vcf_pop2 = list()


        parallel::mclapply(mc.cores=10, 1:length(v1), function(i)
        {


        geneid = id[i]

        vcf1 = readVCF( v1  [ grep(geneid, v1)  ] ,
                       maxNumberOfVariants = 5000,  min_maf =  1e-6, max_maf = 0.01)


        vcf2 = readVCF( v2  [ grep(geneid, v2)  ] ,
                       maxNumberOfVariants = 5000,  min_maf = 1e-6, max_maf = 0.05)


        if(is.na(vcf1)) return(1);
        if(is.na(vcf2)) return(1);

        common = intersect(vcf1$varid,vcf2$varid)

        if(length(common) < 10) return(1);

        vcf1 = subsetVCF(vcf1, common)
        vcf2 = subsetVCF(vcf2, common)

        cat(geneid, common," \n");


        file = sprintf("/tmp/rd-%s",geneid)

        save(vcf1,vcf2,file=file)



      # vcf_pop1[[ geneid ]] = vcf1
      # vcf_pop2[[ geneid ]] = vcf2


    })


}



sapply(vcf_pop1,function(x) length(x$maf))

vcf1 = vcf_pop1[["149-PTPN13"]]

vcf2 = vcf_pop2[["149-PTPN13"]]


