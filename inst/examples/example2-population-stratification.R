

library(sim1000G)
library(SKAT)
library(parallel)



par = commandArgs(T)





ped = read.table("20130606_g1k.ped",h=T,as=T,sep="\t")
pops = split(ped$Individual.ID,list(ped$Population))
id_firstPopulation = c( pops$CEU, pops$TSI )
id_secondPopulation = c( pops$ASW )


#### Initialize the simulation ####

x = list.files(path="~/fs/tmp/genes",full.names = T)
x = sample(x,1)
load(x)
genename = x


#load("~/fs/tmp/genes/rd-89-ADGRA3")

initSimulation = function() {
    startSimulation(vcf1)
    saveSimulation("pop1")

    startSimulation(vcf2)
    saveSimulation("pop2")
}

initSimulation();


if(0) {
maf1 = apply( (vcf$gt1+vcf$gt2)[,id_firstPopulation],1,mean)/2
maf2 = apply( (vcf$gt1+vcf$gt2)[,id_secondPopulation],1,mean)/2

plot(maf1,maf2)
plot(maf1-maf2)

common_causal_pool = which(maf1 > 0 & maf2 > 0 & maf1 < 0.03 & maf2 < 0.03)
}


common_causal_pool = 1:length(vcf1$maf)





SEED = 1

#### Function to generate one simulation replicate ####


computeMaf = function(genotypes) {

    maf = apply(genotypes,2,mean)/2
    maf
}




predictor_logistic = function(b, geno) {
    x = b[1]
    for(i in 1:ncol(geno)) { x = x  + b[i+1] * ( geno[,i] > 0) + b[i+1] * ( geno[,i] > 1)   }
    exp(x) / (1+exp(x) )
}




seed=NA
numcausal=7
effect_size = 0
pop_strat = 10
N_pop1 = 400
N_pop2 = 200
N_aswcontrols = 10



replicate1 = function(seed=NA,
                     numcausal=10,
                     effect_size = 10,
                     pop_strat = 0,
                     N_pop1 = 400,
                     N_pop2 = 200,
                     N_aswcontrols = 10

                    )
{

        if(!is.na(seed)) set.seed(seed)



        loadSimulation("pop1")

        id = generateUnrelatedIndividuals(N_pop1)

        genotypes = retrieveGenotypes(id)
        rownames(genotypes)= rep(  "CEU" , nrow(genotypes) )


        genotypes2 = NA
        gt = genotypes

        if(N_pop2 > 0) {
            loadSimulation("pop2")
            id = generateUnrelatedIndividuals(N_pop2)


            genotypes2 = retrieveGenotypes(id)
            rownames(genotypes2)= rep(  "ASW" , nrow(genotypes2) )

            gt = rbind(genotypes,genotypes2)
        }




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




        S=0

        while(1) {
                p = rep(NA, nrow(gt) )

                s = rownames(gt) != "ASW"

                p[s] = predictor_logistic (  c(S,effect_sizes) ,  gt[s,])
                p[!s] = predictor_logistic (  c(S+pop_strat,effect_sizes) ,  gt[!s,])

                phenotype = rbinom( length(p) , 1 , p )


                ncases = sum(phenotype==0)
                ncontrols = sum(phenotype==1)

                cat(ncases,ncontrols,"\n")
                if( abs(ncases-ncontrols) < 5 ) break();

                if(ncases > ncontrols) S = S + 0.1
                if(ncases < ncontrols) S = S - 0.1
        }




        s = which(rownames(gt) == "CEU")



        print(table(phenotype, rownames(gt)))




        t=( table(rownames(gt)  , phenotype ) )
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


        cat("EFF= ", effect_size,
            "POP1=",N_pop1,
            " POP_STRAT=" , pop_strat,
            " SKATPV=", pv,
            " SKATPV_COV=", pv_cov,
            " GENE=", genename,
            # " t1= ", table(rownames(gt)  , phenotype  ),

            "\n");


        tbl = table(rownames(gt),phenotype)
        data.frame(ID="EFFZZ", eff=effect_size,strat=pop_strat,pop1=N_pop1, pop2=N_pop2,
                   pv1=pv,pv2=pv_cov,
                   gene=genename

                   )
}





#### Generate multiple simulation replicates ####


results = list()

if(0) {
    library(parallel)
    cl <- makeCluster(getOption("cl.cores", 2))
    clusterExport(cl,c("vcf1","vcf2","initSimulation","replicate1",ls() ))
    clusterEvalQ(cl, library(sim1000G))
    clusterEvalQ(cl, initSimulation())
}



#for(pop_strat in c(0,0.5,1,2))
#for(EFF in c(0,log(1.5),log(3)))
#for(p2 in c(0,100,200))

par = as.numeric(par)
EFF = par[1]
p2 = par[2]
pop_strat = par[3]

nsim = par[4]


    {

        # clusterExport(cl,c("p2","EFF","pop_strat"))

        z2 = mclapply(mc.cores=12,mc.preschedule=T,
                      1:nsim, function(x)

                  replicate1(seed=NA, numcausal = 10, effect_size = EFF,
                            N_pop1 = 400-p2,
                            N_pop2 = p2,
                            pop_strat = pop_strat,
                            N_aswcontrols = NA
                  )

                    )
        print(z2)
        results = c(results,z2)

        save(results,file="sim-results.rdata")


    }


pz = plyr::ldply(results)

options(width = 3000)

print(pz)


rand_name = function() {

    x = sprintf("%04d",round(runif(10,1,1000)) )
    paste(x,sep="",collapse="")

}






pz$pv2 [pz$pv2 < 0] = NA

reshape2::dcast(pz, eff + factor(pop2) ~ 1,value.var="pv1",function(x) mean(x<0.05))
reshape2::dcast(pz, eff + factor(pop2) ~ 1,value.var="pv2",function(x) mean(x<0.05))


save(results,pz,file=sprintf("sim-results-%s.rdata",rand_name()) )




sdfdsdf();




#replicate(1,2,-5,400,200)





p2 = 220
replicate(seed=NA, numcausal = 10, effect_size = 4,
          N_pop1 = 400-p2,
          N_pop2 = p2,
          pop_strat = 3,
          N_aswcontrols = NA
)



replicate(seed=NA, numcausal = 5, effect_size = 0,
          N_pop1 = 400,
          N_pop2 = 200,
          pop_strat = 1,
          N_aswcontrols = NA
)




##### Other pieces of code that are not used now #####




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

    X1 = cor(( genotypes ))^2
    X2 = cor(( genotypes2 ))^2

    z1 = X2[ lower.tri(X2[1:50,1:50]) ]
    z2 = X2orig[ lower.tri(X2orig[1:50,1:50]) ]
    plot(z1,z2)



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
