# Load packages for parallel computing
library(snow, lib.loc = "/home/briollaislab/akrol/R/x86_64-pc-linux-gnu-library/3.3")
library(iterators, lib.loc = "/home/briollaislab/akrol/R/x86_64-pc-linux-gnu-library/3.3")
library(foreach, lib.loc = "/home/briollaislab/akrol/R/x86_64-pc-linux-gnu-library/3.3")
library(doSNOW, lib.loc = "/home/briollaislab/akrol/R/x86_64-pc-linux-gnu-library/3.3")
library(sim1000G, lib.loc = "/home/briollaislab/akrol/R/x86_64-pc-linux-gnu-library/3.3")
cl <- makeCluster(12, type="SOCK") #define how many jobs you want the computer to run at the same time, define number of cores to use

clusterEvalQ(cl, library(stringr, lib.loc="/home/briollaislab/akrol/R/x86_64-pc-linux-gnu-library/3.3"))
registerDoSNOW(cl) #use cluster cl
clusterEvalQ(cl, library(hapsim, lib.loc="/home/briollaislab/akrol/R/x86_64-pc-linux-gnu-library/3.3"))
registerDoSNOW(cl) #use cluster cl
clusterEvalQ(cl, library(sim1000G, lib.loc="/home/briollaislab/akrol/R/x86_64-pc-linux-gnu-library/3.3"))
registerDoSNOW(cl) #use cluster cl
clusterEvalQ(cl, library(frailtypack, lib.loc="/home/briollaislab/akrol/R/x86_64-pc-linux-gnu-library/3.3"))
registerDoSNOW(cl) #use cluster cl
clusterEvalQ(cl, library(FamEvent, lib.loc="/home/briollaislab/akrol/R/x86_64-pc-linux-gnu-library/3.3"))
registerDoSNOW(cl) #use cluster cl



set.seed(120)
# Initialize simulator of family data using 1000 genome haplotypes for each cluster
clusterEvalQ(cl, {setwd("/home/briollaislab/akrol/sim1000G/")
vcf <- readVCF("region.vcf.gz", maxNumberOfVariants = 100, min_maf = 0.02, max_maf = 0.1)
readGeneticMap(chromosome = 4)
startSimulation(vcf, totalNumberOfIndividuals = 2000)})
registerDoSNOW(cl) #use cluster cl

# Which Scenario
Scenario <- "A"
if(Scenario == "A"){
  results <- "/home/briollaislab/akrol/Parallel_simulations/Simulations_Apostolos/Results/ScenarioA_mIBD0.7.csv"
}else if(Scenario == "B"){
  results <- "/home/briollaislab/akrol/Parallel_simulations/Simulations_Apostolos/Results/ScenarioB_mIBD0.7.csv"
}
# Number of families and replications
nfam <- 50
nbou <- 500
# Define parameters of the model
Beta <- c(0.5,1.0,0.5,0.3)
ver <- length(Beta)
sig2 <- 0.7
Weibull <- c(0.007,3)
variation <- "IBD"

res.sim.bf <- foreach(ii=1:nbou, .packages=NULL,.combine=rbind,.inorder=FALSE) %dopar% {


   # Simulate families
    SIM$reset()
    ## For simulations, we set the cM that the regions spans to 4000 to diversify the IBD matrix
    ## Remove for normal use
    SIM$cm = seq( 0,4000, length = length(SIM$cm) )

    # Generate family pedigries and their genotypes
    time_x_families <- function(){
      fam <- lapply(paste(1:nfam), function(x){ noffspring2 <- sample(c(1:2), 1, replace=TRUE,prob=c(0.5, 0.5))#, 0.1482
      noffspring3 = sample(c(1,2), noffspring2, replace=TRUE,prob=c(0.8518, 0.1482))#rep(1, noffspring2)#
      newFamily3generations(x,  noffspring2 = noffspring2, noffspring3)})
      fam <- do.call(rbind, fam)

      fam
    }

    fam <- time_x_families()

    # Transform the data to have ids from 1 to n, n - number of individuals
    nmemb <-as.vector(table(fam$fid))
    cumsum_nmemb <- c(0,cumsum(nmemb))
    fam$ID <- 1:dim(fam)[1]
    fam$fatherID <- 0
    fam$fatherID <- sapply(1:dim(fam)[1], function(x) fam$ID[as.numeric(fam$father[which(fam$ID==x)])+cumsum_nmemb[fam$fid[x]]])
    fam$fatherID <- ifelse(fam$father==0, 0, fam$fatherID)
    fam$fatherID <- as.numeric(fam$fatherID)
    fam$motherID <- 0
    fam$motherID <- sapply(1:dim(fam)[1], function(x) fam$ID[as.numeric(fam$mother[which(fam$ID==x)])+cumsum_nmemb[fam$fid[x]]])
    fam$motherID <- ifelse(fam$mother==0, 0, fam$motherID)
    fam$motherID <- as.numeric(fam$motherID)
    fam <- fam[,c(1,8:10,5:7)]
    colnames(fam) <- c("famID", "ID", "fatherID", "motherID", "gender", "generation", "gtindex")



    # Create data of genotypes
    genotype <- rep(0,length(SIM$gt1[1,]))
    for(i in  1:dim(fam)[1]){
      gen <- paste(SIM$gt1[as.numeric(fam$gtindex[i]),],SIM$gt2[as.numeric(fam$gtindex[i]),], sep = "")
      gen <- ifelse(gen == "00", 3, gen)
      gen <- ifelse(gen%in%c("10","01"), 2, gen)
      gen <- ifelse(gen == "11", 1, gen)
      genotype <- rbind(genotype, as.numeric(gen))
    }

    genotype <- genotype[-1,]
    colnames(genotype) <- paste("Gen", 1:dim(genotype)[2], sep="")

    # Create mean IBD matrix
    n = SIM$individuals_generated
    vec <- rep(0, length(nmemb))
    vec2 <- vec
    vec2[2:length(vec)] <- vec[2:length(vec)] + nmemb[-length(nmemb)]
    vec3 <- cumsum(vec2)+1
    nmemb2 <- cumsum(nmemb)
    IBD2matrix <- matrix(rep(0,n^2), nrow = n)


    meanIBD = function(x) { x[2] + x[1]/2 }

    for(i in 1:length(nmemb)){
      IBD2matrix[vec3[i]:nmemb2[i],vec3[i]:nmemb2[i]] <- sapply(vec3[i]:nmemb2[i], function(y) {
        z = sapply(vec3[i]:nmemb2[i], function(x) meanIBD( computePairIBD12(x,y) ) )
        # names(z) = 1:n
        z
      })
    }
    colnames(IBD2matrix) = 1:nrow(IBD2matrix);rownames(IBD2matrix) = 1:nrow(IBD2matrix)

    # Generate time-to-event data using simfam function from the package FamEvent
    if(Scenario=="A"){
      fam_FamEvent <- simfam(N.fam = length(unique(fam$famID)), design = "pop+", variation = variation,
                           base.dist = "Weibull", frailty.dist = "lognormal", depend = sig2,
                           base.parms = Weibull, vbeta = Beta,
                           agemin = 0, data.fam = fam, genotype = NULL, IBD = IBD2matrix,
                          age1 = c(95, 2.5), age2 = c(75, 2))
    }else if(Scenario == "B"){
      fam_FamEvent <- simfam(N.fam = length(unique(fam$famID)), design = "pop+", variation = variation,
                             base.dist = "Weibull", frailty.dist = "lognormal", depend = sig2,
                             base.parms = Weibull, vbeta = Beta,
                             agemin = 0, data.fam = fam, genotype = genotype[,1:3], IBD = IBD2matrix,
                             age1 = c(95, 2.5), age2 = c(75, 2))
    }

    fam_FamEvent$gender2 <- ifelse(fam_FamEvent$gender == 0, 2, fam_FamEvent$gender)
    family <- fam_FamEvent[order(fam_FamEvent$indID),]

    family$t0 <- 0 # for frailtypack, we create dummy variable 0 for starting time of observation

    # Fit the model using function frailtyPenal from frailtypack package
    fit <- frailtyPenal(Surv(t0, time, status) ~ gender +cluster(famID), data = family, hazard="Weibull",
                        RandDist = "LogN", print.times = FALSE, init.B = Beta[1], covMatrix = IBD2matrix,
                        recurrentAG = TRUE, maxit = 35,
                        proband = family$proband, currentage = family$currentage)

    # If the model converged, we save the results
    if(all(fit$EPS<1e-3)){
      iter <- fit$n.iter
      # Data characteristics
      censored_rate <- (dim(family)[1]-sum(family$status))/dim(family)[1]
      N_events_per_family <-  sum(family$status)/length(unique(family$famID))
      if(dim(table(family$status,family$famID))[1]==2){ N_families_no_event <- length(which(table(family$status,family$famID)[2,]==0))}else{N_families_no_event<-0}
      Av_time_observed <-   mean(family$time[family$status ==1])
      N_members_per_family <- dim(family)[1]/length(unique(family$famID))
      dim_data <- dim(family)[1]
      beta <- as.vector(fit$coef)
      nvar <- length(beta)

      # Estimations of the model parameters
      if(length(fit$varH)>1){
        se <- sqrt(diag(fit$varH))
      }else{se <- sqrt(fit$varH)}


      beta_sex<-beta[1]
      beta_mgene1<-beta[2]
      beta_mgene2<-beta[3]
      beta_mgene3<-beta[4]
      var_frailty <- fit$b[3]^2

      sebeta_sex <- se[1]
      sebeta_mgene1 <- se[2]
      sebeta_mgene2 <- se[3]
      sebeta_mgene3 <- se[4]
      sevar_frailty_delta <- sqrt(fit$varTheta[1])
      sevar_frailty <-  sqrt(diag(fit$varHtotal)[3])
      shape <- fit$b[1]^2
      scale <-  fit$b[2]^2
      seshape <- sqrt(fit$varHtotal[1,1])
      seshape_delta <-  sqrt((2*fit$b[1]^2)*fit$varHtotal[1,1])
      sescale <- sqrt(fit$varHtotal[2,2])
      sescale_delta <-  sqrt((2*fit$b[2]^2)*fit$varHtotal[2,2])



    }else{
      # If the model did not converge we fix all the values to 0
    	censored_rate <-0; N_events_per_family <-0; N_families_no_event <-0;Av_time_observed <-0; N_members_per_family <-0; dim_data <-0; iter <-0
	    beta_sex <-0; beta_mgene1 <-0; beta_mgene2 <-0; beta_mgene3 <-0; var_frailty <-0; sebeta_sex <-0; sebeta_mgene1 <-0; sebeta_mgene2 <-0
      sebeta_mgene3 <-0; sevar_frailty_delta <-0; sevar_frailty <-0; shape <-0; scale <-0; seshape <-0; seshape_delta <-0; sescale <-0; sescale_delta <-0
    }
    c(iteration=ii,censored_rate = censored_rate, N_events_per_family = N_events_per_family, N_families_no_event = N_families_no_event,
	    Av_time_observed = Av_time_observed, N_members_per_family = N_members_per_family, dim_data = dim_data, n.iter = iter,
	    beta_sex = beta_sex, beta_mgene1 = beta_mgene1, beta_mgene2 = beta_mgene2, beta_mgene3 = beta_mgene3,
	    var_frailty = var_frailty, sebeta_sex =sebeta_sex, sebeta_mgene1 = sebeta_mgene1 , sebeta_mgene2 = sebeta_mgene2,
      sebeta_mgene3 = sebeta_mgene3, sevar_frailty_delta = sevar_frailty_delta, sevar_frailty = sevar_frailty,
      shape = shape, scale = scale, seshape = seshape, seshape_delta = seshape_delta, sescale = sescale,
      sescale_delta =sescale_delta)

}

# Save the results
write.csv(res.sim.bf, results)
