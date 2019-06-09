# install.packages("/Users/Xiong/Downloads/sim1000G_1.36.tar.gz")
#
# readr hapsim


library(sim1000G)
vcf_file = "CEU-TSI-GBR-region-chr4-357-ANK2.vcf.gz" #nvariants = 442, ss=1000


vcf = readVCF( vcf_file, maxNumberOfVariants = 442  ,min_maf = 0.0005,max_maf = 0.01) #lowest MAF
dim( vcf$gt1 ) #rows represent number of variants, columns represent number of individuals



## Download and use full chromosome genetic map
downloadGeneticMap(4)
readGeneticMap(4)



sample.size=3000



data_sim = function(seed.num){
  startSimulation(vcf, totalNumberOfIndividuals = sample.size) #what does this step do, usually how can I choose number of individuals

  SIM$reset() #what does this step do

  id = generateUnrelatedIndividuals(sample.size)

  gt = retrieveGenotypes(id)


  freq = apply(gt,2,sum)/(2*nrow(gt))
  causal = sample(setdiff(1:ncol(gt),which(freq==0)),45)

  beta.sign = rep(1,45)
  c.value = 0.402
  beta.abs = c.value*abs(log10(freq[causal]))
  beta.val = beta.sign*beta.abs
  x.bar = apply(gt[,causal],2,mean)
  x.bar = as.matrix(x.bar)
  beta.val = t(as.matrix(beta.val))
  #disease prvalance = 1%
  #beta0 = -log(99)-beta.val %*% x.bar
  #disease prvalance = 1.5%
  beta0 = 0-beta.val %*% x.bar

  eta = beta.val %*% t(gt[,causal])
  eta = as.vector(eta) + rep(beta0,nrow(gt))
  prob = exp(eta)/(1+exp(eta))

  genocase = rep(NA, sample.size)

  set.seed(seed.num)
  for(i in 1:sample.size){
    genocase[i] = rbinom(1, 1, prob[i])
  }
  case.idx = sample(which(genocase==1),1000)
  control.idx = sample(which(genocase==0),1000)

  return(rbind(gt[case.idx,],gt[control.idx,]))
}


library(SKAT)
  res = NULL

  for(seed_number in 1:100){
    set.seed(seed_number)
    print(seed_number)

    Z1.skat = data_sim(seed_number)
    # 72 variants, 1000 individuals
    # write.csv(Z1.skat,paste("/Volumes/Briollais_lab/jingxiong/Project/sequencing_data/sim1000G/result/simulated_dataset/data_SKAT_2000_",seed_number,".csv",sep=""),row.names = F,quote=F)
    # 147 variants, 1000 individuals
    # write.csv(Z1.skat,paste("/Volumes/Briollais_lab/jingxiong/Project/sequencing_data/sim1000G/result/simulated_dataset/data_SKAT_2000_147_",seed_number,".csv",sep=""),row.names = F,quote=F)
    # 442 variants, 1000 individuals
    write.csv(Z1.skat,paste("/Volumes/Briollais_lab/jingxiong/Project/sequencing_data/sim1000G/result/simulated_dataset/data_SKAT_2000_442_",seed_number,".csv",sep=""),row.names = F,quote=F)

    obj = SKAT_Null_Model(c(rep(1,1000),rep(0,1000)) ~ 1,out_type="D")

    p1_skat = SKAT(as.matrix(Z1.skat),obj)$p.value
    p1_burden = SKAT(as.matrix(Z1.skat),obj,r.corr=1)$p.value
    p1_skat_O = SKAT(as.matrix(Z1.skat),obj,method="optimal.adj")$p.value

    res = rbind(res,c(p1_skat,p1_burden,p1_skat_O))
  }


  write.csv(res,"/Volumes/Briollais_lab/jingxiong/Project/sequencing_data/sim1000G/result/SKAT/skat_442_2000.csv",row.names = F,quote=F)






