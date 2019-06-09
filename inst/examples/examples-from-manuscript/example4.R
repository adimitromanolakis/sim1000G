library(sim1000G)
library(FamEvent)
library(frailtypack)

set.seed(122)
# Initialize simulator of family data using 1000 genome haplotypes for each cluster
vcf <- readVCF("../region.vcf.gz", maxNumberOfVariants = 100, min_maf = 0.02, max_maf = 0.1)
readGeneticMap(chromosome = 4)
startSimulation(vcf, totalNumberOfIndividuals = 2000)


# Number of families
nfam <- 200
# Define parameters of the model
Beta <- c(0.5,1.0,0.5,0.3)
ver <- length(Beta)
sig2 <- 0.7
Weibull <- c(0.007,3)
variation <- "IBD"


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

  for(i in 1:length(nmemb)){
    IBD2matrix[vec3[i]:nmemb2[i],vec3[i]:nmemb2[i]] <- sapply(vec3[i]:nmemb2[i], function(y) {
      z = sapply(vec3[i]:nmemb2[i], function(x) computePairIBD12(x,y))
      # names(z) = 1:n
      z
    })
  }
  colnames(IBD2matrix) = 1:nrow(IBD2matrix);rownames(IBD2matrix) = 1:nrow(IBD2matrix)

  # Generate time-to-event data using simfam function from the package FamEvent
  fam_FamEvent <- simfam(N.fam = length(unique(fam$famID)), design = "pop+", variation = variation,
                           base.dist = "Weibull", frailty.dist = "lognormal", depend = sig2,
                           base.parms = Weibull, vbeta = Beta,
                           agemin = 0, data.fam = fam, genotype = NULL, IBD = IBD2matrix,
                           age1 = c(95, 2.5), age2 = c(75, 2))

  fam_FamEvent$gender2 <- ifelse(fam_FamEvent$gender == 0, 2, fam_FamEvent$gender)
  family <- fam_FamEvent[order(fam_FamEvent$indID),]

  family$t0 <- 0 # for frailtypack, we create dummy variable 0 for starting time of observation

  # Fit the model using function frailtyPenal from frailtypack package
  fit <- frailtyPenal(Surv(t0, time, status) ~ gender +cluster(famID), data = family, hazard="Weibull",
                      RandDist = "LogN", print.times = FALSE, init.B = Beta[1], covMatrix1 = IBD2matrix,
                      recurrentAG = TRUE, maxit = 35,
                      proband = family$proband, currentage = family$currentage)


  #Genotypes included in generating time-to-event data
    fam_FamEvent <- simfam(N.fam = length(unique(fam$famID)), design = "pop+", variation = variation,
                           base.dist = "Weibull", frailty.dist = "lognormal", depend = sig2,
                           base.parms = Weibull, vbeta = Beta,
                           agemin = 0, data.fam = fam, genotype = genotype[,1:3], IBD = IBD2matrix,
                           age1 = c(95, 2.5), age2 = c(75, 2))


  fam_FamEvent$gender2 <- ifelse(fam_FamEvent$gender == 0, 2, fam_FamEvent$gender)
  family <- fam_FamEvent[order(fam_FamEvent$indID),]

  family$t0 <- 0 # for frailtypack, we create dummy variable 0 for starting time of observation

  # Fit the model using function frailtyPenal from frailtypack package
  fit2 <- frailtyPenal(Surv(t0, time, status) ~ gender +cluster(famID), data = family, hazard="Weibull",
                      RandDist = "LogN", print.times = FALSE, init.B = Beta[1], covMatrix1 = IBD2matrix,
                      recurrentAG = TRUE, maxit = 35,
                      proband = family$proband, currentage = family$currentage)


