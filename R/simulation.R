

#' Description of parameters
#'
#' @param variants_per_individual number of variants per individual (n)
#' @param non_missing_sites numer of non missing sites per individual
#' @param pheno phenotype in 0,1
#' @param method Which method to use
#'
#' @return BF
#' @seealso \code{\link{BF}}
#' @export
testfun = function(variants_per_individual, non_missing_sites, pheno, method, verbose = F) {
  
}



pkg.env <- new.env()

###' @xexport
###set = function(x) { pkg.env$x = x }

####' @xexport
###get = function() { pkg.env$x }





#' @export
SIM = new.env()





#' @export
startSimulation = function(vcf, totalNumberOfIndividuals = 250, randomdata = 0) {
        
        cat("[#####...] Creating SIM object\n");
        
        # Original haplotypes from 1000 genomes / other data
        
        SIM$population_gt1 = vcf$gt1
        SIM$population_gt2 = vcf$gt2
        
        
        
        if(randomdata) {
            
            str(SIM$population_gt1)
            # indoviduals in columns
            
            x = SIM$population_gt1*0
            
            for(i in 1:nrow(x)) x[i,] = rbinom( ncol(x), 1 , 0.25)
            
            SIM$population_gt1 = x
            SIM$population_gt2 = x
            
            # SIM$haplodata = haplodata(  t( cbind(SIM$population_gt1, SIM$population_gt2) ) ) 
            #    ldplot(SIM$haplodata$cor,ld.type="r")
        
        }
        
        
        
        
        # Generate haplodata object 
        
        dim( t( cbind(SIM$population_gt1, SIM$population_gt2) )  )
        
        
        SIM$haplodata = haplodata(  t( cbind(SIM$population_gt1, SIM$population_gt2) ) ) 
        cat("[#####...] Haplodata object created\n");
        
        
        
        # Variant information
        
        SIM$varinfo = vcf$vcf[,1:8]
        SIM$bp = vcf$vcf[,2]
        SIM$cm = approx( geneticMap$bp, geneticMap$cm, SIM$bp )$y
        
        SIM$N_markers = nrow(vcf$gt1)
        
        SIM$pool = NA
        SIM$npool = 0
        
        SIM$total_individuals = totalNumberOfIndividuals
        SIM$individuals_generated = 0
        
        SIM$gt1 = matrix(nrow = SIM$total_individuals, ncol = SIM$N_markers, -1)
        SIM$gt2 = matrix(nrow = SIM$total_individuals, ncol = SIM$N_markers, -1)
        
        
        SIM$origin1 = matrix(nrow = SIM$total_individuals, ncol = SIM$N_markers, -1)
        SIM$origin2 = matrix(nrow = SIM$total_individuals, ncol = SIM$N_markers, -1)
        
        
        #SIM$cm = seq(0, 3200, l=dim(SIM$gt1)[2] )
        SIM$npool = 0
        SIM$last_ancestral_index = 10
        

}



SIM$generateNewHaplotypes = function(n = -1) {
    
    if(SIM$npool < 2) {
        cat("Generate new individual pool n=200\n");
        
        
        SIM$npool = 500
        SIM$pool = haplosim2(SIM$npool, SIM$haplodata, summary = F)$data
        
        
    }
    
    GT = list(gt1 =  SIM$pool[SIM$npool,], gt2 = SIM$pool[SIM$npool-1,] )
    SIM$npool = SIM$npool - 2
    
    SIM$last_ancestral_index = SIM$last_ancestral_index  + 1
    
    return(GT)
    
}


SIM$addUnrelatedIndividual = function() {
    
    newGenotypes = SIM$generateNewHaplotypes()
    
    if(SIM$individuals_generated >= SIM$total_individuals) {
        
        stop("No more space for saving new individual genotypes")
        
    }
    
    
    SIM$individuals_generated = SIM$individuals_generated + 1
    j = SIM$individuals_generated
    
    cat("Adding individual ",j, " from pool\n");
    
    
    SIM$gt1[j,] = newGenotypes$gt1
    SIM$gt2[j,] = newGenotypes$gt2
    
    SIM$origin1[j,] = SIM$last_ancestral_index
    SIM$origin2[j,] = -SIM$last_ancestral_index
    
    
    return(j)
}


SIM$addIndividualFromGenotypes = function(gt1,gt2) {
    
    if(SIM$individuals_generated >= SIM$total_individuals) {
        
        stop("No more space for generating new individual genotypes")
        
    }
    
    SIM$individuals_generated = SIM$individuals_generated + 1
    j = SIM$individuals_generated
    
    cat("Adding individual ",j, " from specified genotypes\n");
    
    SIM$gt1[j,] = gt1
    SIM$gt2[j,] = gt2
    
   
    return(j)
}




debug_flag = 0


SIM$mate = function(i, j) {
    
    # For now, we hope i,j are of different sex
    
    recomb1 = generateSingleRecombinationVector(SIM$cm)
    recomb2 = generateSingleRecombinationVector(SIM$cm)
    
    
    # FATHER1 = SIM$gt1[i,]
    # FATHER2 = SIM$gt2[i,]
    # MOTHER1 = SIM$gt1[i,]
    # MOTHER2 = SIM$gt2[i,]
    
    
    
    gt1 =  recomb1 * SIM$gt1[i,]  +  (1-recomb1) * SIM$gt2[i,]
    gt2 =  recomb2 * SIM$gt1[j,]  +  (1-recomb2) * SIM$gt2[j,]
    
    
    
    #gt1 = SIM$gt1[i,]
    #gt1[recomb1==1] = SIM$gt2[i,][recomb1==1]
    
    #gt2 = SIM$gt1[j,]
    #gt2[recomb1==1] = SIM$gt2[j,][recomb1==1]
    
    #cat("Mate: ",i," " ,j,"\n");
    
    #p1 = paste(SIM$gt1[i,],collapse="")
    #p2 = paste(SIM$gt2[i,],collapse="")
    #p3 = paste(gt1,collapse="")
    
    #cat("R1:", paste(recomb1,collapse=""),"\n",sep="")
    #cat("R2:", paste(recomb2,collapse=""),"\n",sep="")
    
    #cat(p1,"\n",p2,"\n",p3,"\n",sep="")
    
    SWAP = runif(1) < 0.5
    if(SWAP) { t = gt1; gt1 = gt2; gt2 = t; }
    
    
    index = SIM$addIndividualFromGenotypes(gt1, gt2)
    

    #gt1 =  recomb1 * SIM$gt1[i,]  +  (1-recomb1) * SIM$gt2[i,]
    #gt2 =  recomb2 * SIM$gt1[j,]  +  (1-recomb2) * SIM$gt2[j,]
    
    if(!SWAP) {
        SIM$origin1[index,] = recomb1 * SIM$origin1[i,] + (1-recomb1) * SIM$origin2[i,];
        SIM$origin2[index,] = recomb2 * SIM$origin1[j,] + (1-recomb2) * SIM$origin2[j,];
    } else {
        SIM$origin2[index,] = recomb1 * SIM$origin1[i,] + (1-recomb1) * SIM$origin2[i,];
        SIM$origin1[index,] = recomb2 * SIM$origin1[j,] + (1-recomb2) * SIM$origin2[j,];
    }

    # last_gt <<- gt1 + gt2
    
    return (index)    
    
}



SIM$reset = function() {
    SIM$individuals_generated = 0
}

    
    
    
     
    

















#' @export
newNuclearFamily = function(family_id) {
    
    fam = data.frame(family_id = family_id  , 
                     id = c(1,2,3) , 
                     father = c(0,0,1) , 
                     mother = c(0,0,2) , 
                     sex = c(1,2,1)
    )
    
    
    j1 = SIM$addUnrelatedIndividual()
    j2 = SIM$addUnrelatedIndividual()
    j3 = SIM$mate(j1,j2)
    
    
    fam$gtindex = c(j1,j2,j3)
    
    fam
}


#' @export
newFamilyWithOffspring = function(family_id, noffspring = 2) {
    
    fam = data.frame(fid = family_id  , 
                     id = c(1:2) , 
                     father = c(0,0), 
                     mother = c(0,0), 
                     sex = c(1,2)
    )
    
    
    j1 = SIM$addUnrelatedIndividual()
    j2 = SIM$addUnrelatedIndividual()
    
    fam$gtindex = c(j1,j2)
    
    for(i in 1:noffspring) {
        j3 = SIM$mate(j1,j2)
        newFamLine = c(family_id, i+10, 1,2, 1 , j3)
        fam = rbind(fam, newFamLine)
    }
    
    return (fam)
}






#' @export
computePairIBD12 = function(i,j) {
    
    q1 = SIM$origin1[i,]
    q2 = SIM$origin2[i,]
    
    r1 = SIM$origin1[j,]
    r2 = SIM$origin2[j,]
    
    
    table(q1,q2)
    
    IBD1H1 = q1 == r1 | q1 == r2
    IBD1H2 = q2 == r1 | q2 == r2
    
    
    IBD12 = mean(IBD1H1 | IBD1H2)
    IBD2 = mean(   (q1 == r1 & q2 == r2 )  | (q1 == r2 & q2 == r1) )
    
    
    c(IBD1=IBD12-IBD2, IBD2=IBD2)
}







#' @export
printMatrix = function(m) {
    cat("      " , " " , sprintf("[%4d] ", 1:nrow(m) )  , "\n")   
    
    for(i in 1:nrow(m)) {
        
        cat(sprintf("[%4d]",i) , " " , sprintf(" %.3f ", m[i,]) , "\n")   
    }
    
}


















