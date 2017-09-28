#### Load VCF file #####
setwd("~/fs/1000genomes/")

library(sim1000G)

ped = read.table("20130606_g1k.ped",h=T,as=T,sep="\t")

id1 = ped$Individual.ID [ ped$Population == "CEU" ]
id2 = ped$Individual.ID [ ped$Population == "ASW" ]

cat(c(id1,id2),file="/tmp/samples1.txt",sep="\n")


vcf  = readVCF("~/fs/tmp/chr4-80-filt.vcf.gz",   maxNumberOfVariants = 200 ,
               min_maf = 0.05, max_maf = 0.5)


gplots::heatmap.2(cor(t( vcf$gt1+vcf$gt2))^2,col=rev( heat.colors(100) ) ,Rowv=F,Colv=F,trace="none",breaks=seq(0,1,l=101))


#gplots::heatmap.2(cor(t(vcf2$gt1))^2,col=rev( heat.colors(100) ) ,Rowv=F,Colv=F,trace="none")



ids = vcf$individual_ids


id_pop1 = which(ids %in% id1)
id_pop2 = which(ids %in% id2)



gplots::heatmap.2(cor(t( vcf$gt1[,id_pop1]+vcf$gt2[,id_pop1]))^2,col=rev( heat.colors(100) ) ,Rowv=F,Colv=F,trace="none",breaks=seq(0,1,l=101))


gplots::heatmap.2(cor(t( vcf$gt1[,id_pop2]+vcf$gt2[,id_pop2]))^2,col=rev( heat.colors(100) ) ,Rowv=F,Colv=F,trace="none",breaks=seq(0,1,l=101))

#downloadGeneticMap(4)
readGeneticMap(4)





#### Read genetic map and start sim ####


startSimulation(vcf, totalNumberOfIndividuals = 400, subset = id1)
saveSimulation("pop1")


startSimulation(vcf, totalNumberOfIndividuals = 400, subset = id2)
saveSimulation("pop2")

SEED=1




#### Simulate Individuals ####




library(SKAT)

set.seed(SEED)
SEED = SEED+1


loadSimulation("pop1")



length(SIM$individual_ids)


N_pop1 = 90
N_pop2 = 90


N = N_pop1
id = c()
for(i in 1:N) id[i] = SIM$addUnrelatedIndividual()

genotypes = retrieveGenotypes(id)




loadSimulation("pop2")



length(SIM$individual_ids)


N = N_pop2
id = c()
for(i in 1:N) id[i] = SIM$addUnrelatedIndividual()

genotypes2 = retrieveGenotypes(id)





rownames(genotypes)= rep(  "CEU" , nrow(genotypes) )
rownames(genotypes2)= rep(  "ASW" , nrow(genotypes) )



gt = rbind(genotypes,genotypes2)

rownames(gt)




plot(apply(genotypes,2,mean)/2, apply(genotypes2,2,mean)/2)


maf = apply(gt,2,mean,na.rm=T)/2
apply(gt,2,function(x) sum(is.na(x)))
flip  = which(maf > 0.5) ; gt[,flip] = 2 - gt[,flip]
maf = apply(gt,2,mean,na.rm=T)/2

plot(maf)














effect_sizes = rep(0, ncol(gt))
nvar = length(effect_sizes)

s = sample(1:nvar, 60)
effect_sizes[s] = 0






predictor2 = function(b, geno) {
    x = b[1]
    for(i in 1:ncol(geno)) { x = x  + b[i+1] * ( geno[,i] > 0) + b[i+1] * ( geno[,i] > 1)   }
    exp(x) / (1+exp(x) )
}



S=-4
S=-0.4
s = rownames(gt) == "ASW"
p = rep(NA, nrow(gt) )
p[s] = predictor2 (  c(S+2.3,effect_sizes) ,  gt[s,])
p[!s] = predictor2 (  c(S,effect_sizes) ,  gt[!s,])


plot(apply(gt[,s],1,sum))


phenotype = rbinom( length(p) , 1 , p )


table(phenotype)

mean(phenotype[ rownames(gt)  == "CEU" ] )
mean(phenotype[ rownames(gt)  == "ASW" ] )



table(rownames(gt)  , phenotype )



if(0) {
s = which( ( rownames(gt) == "ASW" & phenotype == 0)  )[1:40]
s = setdiff(1:nrow(gt), s)
s


gt = gt[s,]
phenotype = phenotype[s]
table(rownames(gt)  , phenotype )
}



s = apply(gt[,effect_sizes>0],1,sum)
s2 = apply(gt[,],1,sum)

dim(gt)
length(effect_sizes)


mean(s[phenotype==1])
mean(s[phenotype==0])

summary(lm(s ~ phenotype))


summary(lm(s2 ~ phenotype))



#phenotype = sample(phenotype)
obj<-SKAT_Null_Model(phenotype ~ 1, out_type="D")

SKAT((gt),obj)$p.value

#SKATBinary((gt),obj)$p.value




