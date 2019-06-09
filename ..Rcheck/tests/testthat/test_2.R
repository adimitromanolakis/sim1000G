library(sim1000G)


vcf = readVCF("test1.vcf.rdata")


if(0) {  # Create small genetic map for testing
        a = read.table("genetic_map_GRCh37_chr4.txt.gz", as.is=TRUE, header=TRUE)
        str(a)
        colnames(a) = c("chr","bp","rate.cm","cm")
        a = a[ order(a$cm),]
        r = range( SIM$bp )
        a = a [  a$bp > r[1]-10e3 & a$bp < r[2] + 10e3 , ]
        write.table(a,file="test1.geneticmap.txt",sep=" ",quote=F,row=F)
}

readGeneticMapFromFile("test1.geneticmap.txt")



startSimulation(vcf, totalNumberOfIndividuals = 1200)



SIM$reset()


id = c()
for(i in 1:30) id[i] = SIM$addUnrelatedIndividual()

# Show haplotype 1  of first 5 individuals



genotypes = SIM$gt1[1:20,] + SIM$gt2[1:20,]

tbl = table(genotypes)




context("Genotype sanity")

test_that("Check genotypes are read correctly", {
    expect_gt( tbl[1] , 10 );
    expect_gt( tbl[2] , 10 );
    expect_gt( tbl[3] , 10 );
    expect_equal( sum(genotypes>2) , 0 );


})



