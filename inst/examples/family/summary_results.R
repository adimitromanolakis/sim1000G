summary_sim <- function(sum, true.values, genotypes = FALSE, nfam, variation, message){
  
  
  nbou <- dim(sum)[1]
  print(paste(nbou, "models converged"))
  print(paste("mean iter = ", round(mean(sum$n.iter),2)))
  print(paste("mean right-censoring rate", round(mean(sum$censored_rate)*100,1),"%"))
  print(paste("mean number of individuals", round(mean(sum$dim_data),0)))
  print(paste("number of families", nfam))
  print(paste("Average number of family members", round(mean(sum$N_members_per_family),1)))
  print(paste("average number of events per family", round(mean(sum$N_events_per_family), 1)))
  print(paste("Number of families without an event", round(mean(sum$N_families_no_event), 1)))
  print(paste("Average time to observed events", round(mean(sum$Av_time_observed), 1)))
  print(paste("Variation via", variation))
  print(paste(message))
  
  # statistics for beta_sex
  moybeta_sex <- mean(sum$beta_sex)
  
  print(paste('--> mean of beta sex = ',round(moybeta_sex, 3),'; true value = ', true.values[1]))
  
  eca <- sum(sum$beta_sex^2)
  print(paste('---> SD( beta sex) = ',round(sqrt((eca-(sum(sum$beta_sex))^2/nbou)/(nbou-1)),3)))
  
  moysebeta_sex = mean(sum$sebeta_sex)
  print(paste('--> mean of SE ( beta sex) =', round(moysebeta_sex, 3)))
  
  IC_inf=sum$beta_sex-1.96*sum$sebeta_sex
  IC_sup=sum$beta_sex+1.96*sum$sebeta_sex
  cp_ic <- length(which(IC_inf<true.values[1] & IC_sup>true.values[1] ))
  txcouv1<-cp_ic*100/nbou
  print(paste('CI 95 % for beta sex:  ',round(txcouv1,2), '%'))
  
 
   
  if(genotypes == TRUE){
    ngen <- 3
    # statistics for beta_mgene1
    moybeta_mgene1 <- mean(sum$beta_mgene1)
    print(paste('--> mean of beta mgene1 = ', round(moybeta_mgene1, 3),'; true value = ', true.values[2]))
    
    eca <- sum(sum$beta_mgene1^2) 
    print(paste('---> SD( beta mgene1) = ',round(sqrt((eca-(sum(sum$beta_mgene1))^2/nbou)/(nbou-1)),3)))
    
    moysebeta_mgene1 = mean(sum$sebeta_mgene1)
    print(paste('--> mean of SE ( beta mgene1) =',round(moysebeta_mgene1,3) ))
   
    IC_inf=sum$beta_mgene1-1.96*sum$sebeta_mgene1
    IC_sup=sum$beta_mgene1+1.96*sum$sebeta_mgene1
    cp_ic <- length(which(IC_inf<true.values[2] & IC_sup>true.values[2] ))
    txcouv1<-cp_ic*100/nbou
    print(paste('CI 95 % for beta  mgene1:  ',round(txcouv1,2), '%'))
    
  
    # statistics for beta_mgene2
    moybeta_mgene2 <- mean(sum$beta_mgene2)
    
    print(paste('--> mean of beta mgene2 = ', round(moybeta_mgene2, 3),'; true value = ', true.values[3]))
    
    eca <- sum(sum$beta_mgene2^2) 
    print(paste('---> SD( beta mgene2) = ',round(sqrt((eca-(sum(sum$beta_mgene2))^2/nbou)/(nbou-1)),3)))
    
    moysebeta_mgene2 = mean(sum$sebeta_mgene2)
    print(paste('--> mean of SE ( beta mgene2) =',round(moysebeta_mgene2,3) ))
    IC_inf=sum$beta_mgene2-1.96*sum$sebeta_mgene2
    IC_sup=sum$beta_mgene2+1.96*sum$sebeta_mgene2
    cp_ic <- length(which(IC_inf<true.values[3] & IC_sup>true.values[3] ))
    txcouv1<-cp_ic*100/nbou
    print(paste('CI 95 % for beta  mgene1:  ',round(txcouv1,2), '%'))
    
 
    # statistics for beta_mgene3
    moybeta_mgene3 <- mean(sum$beta_mgene3)
    
    print(paste('--> mean of beta mgene3 = ', round(moybeta_mgene3, 3),'; true value = ', true.values[4]))
    
    eca <- sum(sum$beta_mgene3^2) 
    print(paste('---> SD( beta mgene3) = ',round(sqrt((eca-(sum(sum$beta_mgene3))^2/nbou)/(nbou-1)),3)))
    
    moysebeta_mgene3 = mean(sum$sebeta_mgene3)
    print(paste('--> mean of SE ( beta mgene3) =',round(moysebeta_mgene3,3) ))
    IC_inf=sum$beta_mgene3-1.96*sum$sebeta_mgene3
    IC_sup=sum$beta_mgene3+1.96*sum$sebeta_mgene3
    cp_ic <- length(which(IC_inf<true.values[4] & IC_sup>true.values[4] ))
    txcouv1<-cp_ic*100/nbou
    print(paste('CI 95 % for beta  mgene3:  ',round(txcouv1,2), '%'))
    
  }else{ngen <- 0}
  
  # statistics for varfrailty
  moyvar_frailty <- mean(sum$var_frailty)
  
  print(paste('--> mean of variance frailty = ',round(moyvar_frailty, 3 ) ,'; true value = ', true.values[ngen+2]))
  
  
  eca <- sum(sum$var_frailty) 
  print(paste('---> SD(variance frailty) = ',round(sqrt((eca-(sum(sqrt(sum$var_frailty)))^2/nbou)/(nbou-1)),3)))
  
  moysevar_frailty= mean(sum$sevar_frailty)
  print(paste('--> mean of SE (variance frailty) =', round(moysevar_frailty,3 ) ))
  IC_inf=sqrt(sum$var_frailty)-1.96*sum$sevar_frailty
  IC_sup=sqrt(sum$var_frailty)+1.96*sum$sevar_frailty
  cp_ic <- length(which(IC_inf<sqrt(true.values[ngen+2]) & IC_sup>sqrt(true.values[ngen+2]) ))
  txcouv1<-cp_ic*100/nbou
  print(paste('CI 95% for variance frailty:  ',round(txcouv1,2), '%'))
  
  
  
  # statistics for shape
  moyshape <- mean(sum$shape)
  
  print(paste('--> mean of shape = ',round(moyshape, 3 ),'; true value = ', true.values[ngen+3]))
  
  
  eca <- sum(sum$shape)
  print(paste('---> SD(shape) = ',round(sqrt((eca-(sum(sqrt(sum$shape)))^2/nbou)/(nbou-1)),3)))
  
  moyseshape= mean(sum$seshape)
  print(paste('--> mean of SE (shape) =', round(moyseshape,3 ) ))

  IC_inf=sqrt(sum$shape)-1.96*sum$seshape
  IC_sup=sqrt(sum$shape)+1.96*sum$seshape
  cp_ic <- length(which(IC_inf<sqrt(true.values[ngen+3]) & IC_sup>sqrt(true.values[ngen+3]) ))
  txcouv1<-cp_ic*100/nbou
  print(paste('--> mean of SE (shape) =',round(txcouv1,2), '%'))
  
  
  
  # statistics for scale
  moyscale <- mean(sum$scale)
  
  print(paste('--> mean of scale = ',round(moyscale, 3 ),'; true value = ', 1/true.values[ngen+4]))
  
  
  
  eca <- sum(sum$scale)
  print(paste('---> SD(scale) = ',round(sqrt((eca-(sum(sqrt(sum$scale)))^2/nbou)/(nbou-1)),3)))
  
  moysescale= mean(sum$sescale)
  print(paste('--> mean of SE (scale) =', round(moysescale,3 ) ))
  IC_inf=sqrt(sum$scale)-1.96*sum$sescale
  IC_sup=sqrt(sum$scale)+1.96*sum$sescale
  cp_ic <- length(which(IC_inf<sqrt(1/true.values[ngen+4]) & IC_sup>sqrt(1/true.values[ngen+4]) ))
  txcouv1<-cp_ic*100/nbou
  print(paste('--> mean of SE (scale) =',round(txcouv1,2), '%'))
  
  
}


sum<-read.table("ScenarioB_mIBD0.7.csv",dec=".",sep=",",header=T,na.strings = "NA")

sum <- sum[which(sum$censored_rate!=0),] 
true.values <- c(0.5,0.7,3.0,0.007)
genotypes <- FALSE
message <- "Scenario B"
nfam <- 200 
variation <- "mean IBD"

sink("ScenarioB_mIBD0.7_results.txt")
summary_sim(sum, true.values = true.values, genotypes = genotypes, nfam = nfam, variation = variation, message = message)

sink()