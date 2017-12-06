
library(sim1000G)
download.file( "https://adimitromanolakis.github.io/sim1000G/data/region.vcf.gz", destfile = "region.vcf.gz" )

vcf = readVCF( "~/fs/tmp/chr4-300kb-CEU-TSI-GBR.vcf.gz", maxNumberOfVariants = 400 ,
               min_maf = 0.12 , max_maf = NA )


#downloadGeneticMap( 4 )
readGeneticMap( chromosome = 4 )

startSimulation( vcf, totalNumberOfIndividuals = 500 )

SIM$reset()


id = c()
for(i in 1:355) id[i] = SIM$addUnrelatedIndividual()

# Show haplotype 1  of first 5 individuals
#print(SIM$gt1[1:5,1:6])

# Show haplotype 2
#print(SIM$gt1[1:5,1:6])



genotypes = SIM$gt1[1:220,] + SIM$gt2[1:220,]

print(dim(genotypes))

str(genotypes)

library(gplots)
X1 = cor(genotypes)^2
heatmap.2(X1, col=rev(heat.colors(100)) , trace="none",Rowv=F,Colv=F)

heatmap.2(SIM$haplodata$cor^2 , col=rev(heat.colors(100)) , trace="none",Rowv=F,Colv=F)






SIM$pool = hapsim::haplosim(SIM$npool, SIM$haplodata, summary = F)$data
SIM$pool2 = hapsim::haplosim(SIM$npool, SIM$haplodata, summary = F)$data

X1orig = cor(( SIM$pool + SIM$pool2 )) ^2
gplots::heatmap.2(X1orig, col=rev( heat.colors(100) ) ,Rowv=F,Colv=F,trace="none",breaks=seq(0,1,l=101))






X1orig = cor(t( SIM$population_gt1+SIM$population_gt2))^2
gplots::heatmap.2(X1orig, col=rev( heat.colors(100) ) ,Rowv=F,Colv=F,trace="none",breaks=seq(0,1,l=101))


X1orig = cor(t( SIM$population_gt2))^2
gplots::heatmap.2(X1orig, col=rev( heat.colors(100) ) ,Rowv=F,Colv=F,trace="none",breaks=seq(0,1,l=101))

gplots::heatmap.2(X1, col=rev( heat.colors(100) ) ,Rowv=F,Colv=F,trace="none",breaks=seq(0,1,l=101))



f=(SIM$population_gt1+SIM$population_gt2)
f = apply(f,1,mean)

plot(f,apply(genotypes,2,mean))
abline(0,1)


X1 = cor(( genotypes ))^2

z1 = X1[ lower.tri(X1[1:50,1:50]) ]
z2 = X1orig[ lower.tri(X1orig[1:50,1:50]) ]
plot(z1,z2)

