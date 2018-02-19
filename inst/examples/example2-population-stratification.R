#
#for j in 1 2 3 4 5 6 7 8 9 10; do  for i in 1 2 3 4 5 6 7 8 9 10; do
#submit j3-$j-$i Rscript example2-population-stratification.R ; done ; done
#

library(sim1000G)
library(parallel)
library(SKAT)



par = commandArgs(T)

NREP = 15

run_parallel = 1

numcpus = 10

#### Initialize the simulation ####


# List all available gene region vcf files
all_vcf_files = list.files(path="pop1",  full.names = T)
selected_vcf_file = sample(all_vcf_files , 1)


# Find corresponding vcf file of second population
selected_vcf_file2 = gsub("pop1","pop2", selected_vcf_file)



cat("Reading the two vcf files: " , selected_vcf_file,selected_vcf_file2 ,"\n")




initSimulation = function() {


    ## Read two VCF file with regions from two different populations
    ## vcf1: variants from european population
    ## vcf2: variants from african populations


    gn = strsplit(selected_vcf_file,"-")
    gn = sapply(gn, function(x) x[4])

    genename <<- gn
    geneid = gn

    vcf1 <<- readVCF(selected_vcf_file ,
                    maxNumberOfVariants = 800,  min_maf =  1e-6, max_maf = 0.02)


    if(class(vcf1) != "environment") stop(1);


    vcf2 <<- readVCF( selected_vcf_file2  ,
                    maxNumberOfVariants = 800,  min_maf = 1e-6, max_maf = 0.02)


    print(class(vcf1))
    print(class(vcf2))


    if(class(vcf1) != "environment") stop(1);
    if(class(vcf2) != "environment") stop(1);




    ## We select only the common variants between the 2 VCF files
    ##

    common = intersect(vcf1$varid,vcf2$varid)
    print(length(common))

    if(length(common) < 10) stop(0);

    vcf1 <<- subsetVCF(vcf1, common)
    vcf2 <<- subsetVCF(vcf2, common)

    cat(geneid, length(common)," \n");


    common_causal_pool <<- 1:length(vcf1$maf)

    cat("start sim\n");


    startSimulation(vcf1, totalNumberOfIndividuals = 10000)
    saveSimulation("pop1")

    startSimulation(vcf2, totalNumberOfIndividuals = 10000)
    saveSimulation("pop2")


    print(length(vcf1$varid))

}





SEED = 1

#### Function to generate one simulation replicate ####




predictor_logistic = function(b, geno) {
    x = b[1]
    for(i in 1:ncol(geno)) { x = x  + b[i+1] * ( geno[,i] > 0) + b[i+1] * ( geno[,i] > 1)   }
    exp(x) / (1+exp(x) )
}





replicate1 = function(seed=NA,
                     numcausal = 10,
                     effect_size = 10,
                     pop_strat = 0,
                     N_pop1 = 400,
                     N_pop2 = 200
                    )
{

        # if(!is.na(seed)) set.seed(seed)



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



        pca = prcomp(gt)


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

                # cat(ncases,ncontrols,"\n")
                if( abs(ncases-ncontrols) < 5 ) break();

                if(ncases > ncontrols) S = S + 0.1
                if(ncases < ncontrols) S = S - 0.1
        }



        t=( table(rownames(gt)  , phenotype ) )
        print(t)




        #phenotype = sample(phenotype)

        obj<-SKAT_Null_Model(phenotype ~ 1, out_type="D")

        pv = SKAT((gt),obj)$p.value




        if(sum(rownames(gt) != "CEU") > 0) {

            population = factor(rownames(gt))

            # obj<-SKAT_Null_Model(phenotype ~ population, out_type="D")

            obj<-SKAT_Null_Model(phenotype ~
                 pca$x[,1] + pca$x[,2] + pca$x[,3] + pca$x[,4] + pca$x[,5] + pca$x[,6], out_type="D")

            pv_cov = SKAT((gt),obj)$p.value


        } else {
            pv_cov = -1
        }



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


# replicate1(NA,10,0, 2,1000,1000)


#### Generate multiple simulation replicates ####

rand_name = function() {

    x = sprintf("%04d",round(runif(10,1,1000)) )
    paste(x,sep="",collapse="")

}


seed=NA
numcausal=7
effect_size = 0
pop_strat = 10
N_pop1 = 400
N_pop2 = 200





par = as.numeric(par)
EFF = par[1]
p2 = par[2]
pop_strat = par[3]

nsim = par[4]

if(0) {
    EFF=2
    p2 = 100
    pop_strat = 2
    nsim=2
}



base_seed = round( runif(1,0,1e9) )

initSimulation()


if(run_parallel == 1) {


    library(parallel)

    cat("Cluster setup\n");
    cl <- makeCluster(getOption("cl.cores", numcpus))

    clusterExport(cl, c(
            "selected_vcf_file","selected_vcf_file2",
            "initSimulation","replicate1",
            "par","EFF","p2","pop_strat","nsim",
            "predictor_logistic",
            "base_seed"
       ) )

    cat("exportDone\n");

    clusterEvalQ(cl, library(sim1000G))
    clusterEvalQ(cl, library(SKAT))
    cat("initSim\n");

    clusterEvalQ(cl, initSimulation())

    cat("initSim Done\n");



}





clusterRep = function(x) {

    replicate1(seed=base_seed + 101*x, numcausal = 10, effect_size = EFF,
               N_pop1 = 2000-p2,
               N_pop2 = p2,
               pop_strat = pop_strat
    )
}


results = list()



pop_strat=4; EFF=0; p2=1400; clusterRep()





pop_strat = 0
EFF = log(3)
p2 = 0


effect_sizes = c(0,log(1.5), log(1.8) , log(3), log(5)  )

# pop_strat_list = c(0,0.5,1,2)
pop_strat_list = c(2)

for(pop_strat in pop_strat_list)
    for(EFF in effect_sizes)
        for(p2 in c(0,100,200,400))

{


        cat(pop_strat, EFF, p2,"\n");

        if(run_parallel == 1) {
            clusterExport(cl,c("p2","EFF","pop_strat","genename"))

             v = clusterApply(cl,1:NREP, clusterRep)
        }

        if(run_parallel == 0) v = lapply(1:NREP, clusterRep)

        results = c( results  ,  v   )
}





pz = plyr::ldply(results)
pz$MASTER_ID = "REPLICATE"

options(width = 3000)
print(pz)








pz$pv2 [pz$pv2 < 0] = NA

reshape2::dcast(pz, eff + factor(pop2) ~ 1,value.var="pv1",function(x) mean(x<0.05))
reshape2::dcast(pz, eff + factor(pop2) ~ 1,value.var="pv2",function(x) mean(x<0.05))


save(results,pz,file=sprintf("sim-results-%s.rdata",rand_name()) )

write.table(pz,  row=F,quote=F,sep=" ", col=F, file=sprintf("out-%s.txt", rand_name()) )

