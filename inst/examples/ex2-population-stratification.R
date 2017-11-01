#### Load VCF file #####
library(sim1000G)


if( ! "vcf" %in% ls() ) {

setwd("~/fs/1000genomes/")


ped = read.table("20130606_g1k.ped",h=T,as=T,sep="\t")

pops = split(ped$Individual.ID,list(ped$Population))




id_firstPopulation = c( pops$CEU, pops$TSI )
id_secondPopulation = c( pops$ASW )



cat(c(id_firstPopulation,id_secondPopulation),file="/tmp/samples1.txt",sep="\n")


vcf  = readVCF("~/fs/tmp/chr4-80-filt.vcf.gz",
               maxNumberOfVariants = 300 ,
               min_maf = 1e-6, max_maf = 0.01)






ids = vcf$individual_ids


which_pop1 = which(ids %in% id_firstPopulation)
which_pop2 = which(ids %in% id_secondPopulation )

}


if(0) {
    gplots::heatmap.2(cor(t( vcf$gt1+vcf$gt2))^2,
                      col=rev( heat.colors(100) ) ,
                      Rowv=F,Colv=F,trace="none",
                      breaks=seq(0,1,l=101))


    #gplots::heatmap.2(cor(t(vcf2$gt1))^2,col=rev( heat.colors(100) ) ,Rowv=F,Colv=F,trace="none")



X1orig = cor(t( vcf$gt1[,which_pop1]+vcf$gt2[,which_pop1]))^2

gplots::heatmap.2(X1orig, col=rev( heat.colors(100) ) ,Rowv=F,Colv=F,trace="none",breaks=seq(0,1,l=101))


X2orig = cor(t( vcf$gt1[,which_pop2]+vcf$gt2[,which_pop2]))^2
gplots::heatmap.2(X2orig,  col=rev( heat.colors(100) ) ,Rowv=F,Colv=F,trace="none",breaks=seq(0,1,l=101))

}

#downloadGeneticMap(4)

readGeneticMap(4)





#### Read genetic map and start sim ####


startSimulation(vcf, totalNumberOfIndividuals = 1400, subset = id_firstPopulation)
saveSimulation("pop1")



#X1orig = cor(t( SIM$population_gt1+SIM$population_gt2))^2
#gplots::heatmap.2(X1orig, col=rev( heat.colors(100) ) ,Rowv=F,Colv=F,trace="none",breaks=seq(0,1,l=101))



startSimulation(vcf, totalNumberOfIndividuals = 1400, subset = id_secondPopulation)
saveSimulation("pop2")






s =  which( id_firstPopulation %in% colnames(vcf$gt1) )
id_firstPopulation = id_firstPopulation[ s ]

s =  which( id_secondPopulation %in% colnames(vcf$gt1) )
id_secondPopulation = id_secondPopulation[ s ]


str(vcf$gt1)

image(vcf$gt1+vcf$gt2)

maf1 = apply( (vcf$gt1+vcf$gt2)[,id_firstPopulation],1,mean)/2
maf2 = apply( (vcf$gt1+vcf$gt2)[,id_secondPopulation],1,mean)/2

plot(maf1,maf2)
plot(maf1-maf2)



common_causal_pool = which(maf1 > 0 & maf2 > 0 & maf1 < 0.03 & maf2 < 0.03)
common_causal_pool




SEED = 1

#### Simulate Individuals ####


computeMaf = function(genotypes) {

    maf = apply(genotypes,2,mean)/2
    maf
}

library(SKAT)
library(parallel)

seed=NA
numcausal=7
effect_size = 10
pop_strat = 10
N_pop1 = 400
N_pop2 = 200
N_aswcontrols = 10


#for(loop_num in 1:100) {
replicate = function(seed=NA,
                     numcausal=10,
                     effect_size = 10,
                     pop_strat = 0,
                     N_pop1 = 400,
                     N_pop2 = 200,
                     N_aswcontrols = 10

                    )
{

if(!is.na(seed)) set.seed(seed)
#SEED = SEED+1


loadSimulation("pop1")

id = generateUnrelatedIndividuals(N_pop1)

genotypes = retrieveGenotypes(id)
rownames(genotypes)= rep(  "CEU" , nrow(genotypes) )


genotypes2 = NA


if(N_pop2 > 0) {
    loadSimulation("pop2")
    id = generateUnrelatedIndividuals(N_pop2)


    genotypes2 = retrieveGenotypes(id)
    rownames(genotypes2)= rep(  "ASW" , nrow(genotypes2) )
}


if(0) {

X1 = cor(( genotypes ))^2
X2 = cor(( genotypes2 ))^2

z1 = X2[ lower.tri(X2[1:50,1:50]) ]
z2 = X2orig[ lower.tri(X2orig[1:50,1:50]) ]
plot(z1,z2)



f=(SIM$population_gt1+SIM$population_gt2)
f = apply(f,1,mean)

plot(f,apply(genotypes2,2,mean))




gplots::heatmap.2(cor(( genotypes2 ))^2,col=rev( heat.colors(300) ) ,Rowv=F,Colv=F,trace="none",breaks=seq(0,1,l=301))

X1 = cor(( genotypes ))^2
X2 = cor(( genotypes2 ))^2
X1[ lower.tri(X1) ] = X1orig[ lower.tri(X1orig)  ]

gplots::heatmap.2(X1,col=rev( heat.colors(100) ) ,Rowv=F,Colv=F,trace="none",breaks=seq(0,1,l=101))


}





if(0) {
maf1 = apply(genotypes,2,mean)/2
maf2 = apply(genotypes2,2,mean)/2

common_causal_pool = which(maf1 > 0 & maf2 > 0 & maf1 < 0.03 & maf2 < 0.03)
plot(maf1[common_causal_pool],maf2[common_causal_pool])

cat("Similar maf variants: ",length(common_causal_pool),"\n")
}


if(N_pop2 > 0) {
    gt = rbind(genotypes,genotypes2)
} else {
    gt = genotypes
}

#print( table(rownames(gt)) )
#image(gt)
#rownames(gt)




#plot(apply(genotypes,2,mean)/2, apply(genotypes2,2,mean)/2)


maf = apply(gt,2,mean,na.rm=T)/2
apply(gt,2,function(x) sum(is.na(x)))
flip  = which(maf > 0.5) ; gt[,flip] = 2 - gt[,flip]
maf = apply(gt,2,mean,na.rm=T)/2




#EFFECT_SIZE = 10


maf_gt0 = which(maf > 0)

effect_sizes = rep(0, ncol(gt))
nvar = length(effect_sizes)

s = sample(common_causal_pool, numcausal)
effect_sizes[s] = effect_size






predictor2 = function(b, geno) {
    x = b[1]
    for(i in 1:ncol(geno)) { x = x  + b[i+1] * ( geno[,i] > 0) + b[i+1] * ( geno[,i] > 1)   }
    exp(x) / (1+exp(x) )
}



S=-4
S=-0.4
s = rownames(gt) != "ASW"


p = rep(NA, nrow(gt) )
p[s] = predictor2 (  c(S,effect_sizes) ,  gt[s,])
p[!s] = predictor2 (  c(S+pop_strat,effect_sizes) ,  gt[!s,])


#plot(apply(gt[!s,],1,sum))
phenotype = rbinom( length(p) , 1 , p )


s = which(rownames(gt) == "CEU")


ss = apply(gt,1,sum)
summary(lm(phenotype[s] ~ ss[s]))



print(table(phenotype, rownames(gt)))

#table(phenotype)
#mean(phenotype[ rownames(gt)  == "CEU" ] )
#mean(phenotype[ rownames(gt)  == "ASW" ] )



# N_aswcontrols = 10


if(!is.na(N_aswcontrols) ) {



        s1 = which( ( rownames(gt) == "ASW" & phenotype == 0)  )
        s = setdiff(1:nrow(gt), s1)

    if(1) {
        s2 = which( ( rownames(gt) == "ASW" & phenotype == 1)  )
        l = length(s2)
        s2 = sample(s2, l-N_aswcontrols)
        s = setdiff(s, s2)
    }


    gt = gt[s,]
    phenotype = phenotype[s]
    table(rownames(gt)  , phenotype )


}



if(0) {
s = apply(gt[,effect_sizes>0],1,sum)
s2 = apply(gt[,],1,sum)


mean(s[phenotype==1])
mean(s[phenotype==0])
summary(lm(s ~ phenotype))
summary(lm(s2 ~ phenotype))
}



t=( table(rownames(gt)  , phenotype ) )
#for(i in 1:2) t[i,] = t[i,]/sum(t[i,])
print(t)



#phenotype = sample(phenotype)
obj<-SKAT_Null_Model(phenotype ~ 1, out_type="D")

pv = SKAT((gt),obj)$p.value

if(sum(rownames(gt) != "CEU") > 0) {
    population = factor(rownames(gt))
    obj<-SKAT_Null_Model(phenotype ~ population, out_type="D")
    pv_cov = SKAT((gt),obj)$p.value
} else {
    pv_cov = -1
}


cat("SKAT PV " , pv, "\n")


cat("EFF= ", effect_size, " POP_STRAT=" , pop_strat,
    " SKATPV=", pv,
    " SKATPV2=", pv_cov,
    " t1= ", table(rownames(gt)  , phenotype  ), "\n");


tbl = table(rownames(gt),phenotype)
data.frame(eff=effect_size,strat=pop_strat,pop1=N_pop1, pop2=N_pop2, pv1=pv,pv2=pv_cov)
}





#dfgds();


#study1 = function() {
#z = lapply(1:3,function(x) replicate(x,3,-29,300,200,0))
#z
#}

#study1()





results = list()


for(pop_strat in c(0,0.5,1,2))
for(EFF in c(0,log(1.5),log(3)))

for(p2 in c(0,100,200))
    {
        z2 = mclapply(1:50,function(x)

                  replicate(seed=NA, numcausal = 10, effect_size = EFF,
                            N_pop1 = 400-p2,
                            N_pop2 = p2,
                            pop_strat = pop_strat,
                            N_aswcontrols = NA
                  )

                    , mc.cores=10)
        print(z2)
        results = c(results,z2)

        save(results,file="sim-results.rdata")


    }


pz = plyr::ldply(results)

pz$pv2 [pz$pv2 < 0] = NA

reshape2::dcast(pz, eff + factor(pop2) ~ 1,value.var="pv1",function(x) mean(x<0.05))
reshape2::dcast(pz, eff + factor(pop2) ~ 1,value.var="pv2",function(x) mean(x<0.05))


save(results,pz,file="sim-results.rdata")



sdfdsdf();




#replicate(1,2,-5,400,200)





p2 = 220
replicate(seed=NA, numcausal = 10, effect_size = 0.5,
          N_pop1 = 400-p2,
          N_pop2 = p2,
          pop_strat = 1,
          N_aswcontrols = NA


replicate(seed=10, numcausal = 5, effect_size = 0,
          N_pop1 = 400,
          N_pop2 = 0,
          pop_strat = 0,
          N_aswcontrols = NA
)
