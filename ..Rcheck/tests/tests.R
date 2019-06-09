
#hist(generateRecombinationDistances(130000),n=100)

if(0) {
    recombVector = generateSingleRecombinationVector(cm)
    paste(recombVector[seq(1,length(recombVector),by=50)],collapse="")
    
    
    
    
    n = length(cm)
    haplo3 = rep(0, n)
    
    
    
    
    cm = seq(0,331,l = length(1:1000) )
    
    ibd2 = function() {
        
        origin3father = generateSingleRecombinationVector(cm)
        
        origin3mother =  generateSingleRecombinationVector(cm)
        
        origin4father = generateSingleRecombinationVector(cm)
        
        origin4mother = generateSingleRecombinationVector(cm)
        
        mean(origin3father == origin4father | origin3mother == origin4mother) 
        
    }
    
    
    
    system.time ( t <- replicate(300,ibd2()) )
    mean(t)
}


# PLOT: Recombination points for 200 individuals across region
if(0) {
    cm = seq(0,111,l = length(1:1000) )
    
    n1=length(cm)
    n2 = 200
    m = matrix(NA,nrow=n1,ncol=n2)
    m=apply(m,2,function(x) generateSingleRecombinationVector(cm))
    
    #m[,seq(1,n2,by=2)  ] = 0
    
    image(m)
}





# PLOT: Mean number of crossover points in a single region per generation
#hist(replicate(11900, sum(diff(generateSingleRecombinationVector(seq(0,300,l=120)) ) != 0  ) ), breaks=0:12)

# SYSTEM: Timing
#system.time(replicate(1000,   generateRecombinationDistances(10)       ))
#system.time(replicate(1000,   generateSingleRecombinationVector(cm)    ))


if(0) {
    # Compute mean number of recombinations per generation
    
    cm = seq(0,25,l=200)
    v = generateSingleRecombinationVector(cm)
    paste(v,collapse="")
    
    mean( sapply(1:5000, function(i) { v = generateSingleRecombinationVector(cm); abs( max(v)-min(v) ); } ) )
}














test1_timingOfMatingFunction = function() {
    
    SIM$individuals_generated <- 0
    addUnrelatedIndividual()
    addUnrelatedIndividual()
    #addUnrelatedIndividual()
    #addUnrelatedIndividual()
    
    
    
    #mate(1,2)
    #g1 = last_gt
    
    for(i in 1:1) {
        print( system.time(x <- sapply(1:50, function(x) mate(1,2))) )
    }
    
}; 
#test1_timingOfMatingFunction()


test2_pIBS0_problem = function() {
    SIM$individuals_generated <- 0
    
    ibs0 = c()
    for(i in 1:20) {
        x = newFamilyWithOffspring(2,5)
        
        
        ind1 = x$gtindex[3]
        ind2 = x$gtindex[5]
        
        v = pIBS0(ind1,ind2)
        ibs0 = c(ibs0, v)
        
    }
    
    cat("IBS0 of full-siblings: " , ibs0, "\n");
    
}; 
#test2_pIBS0_problem()




pIBS0 = function(i,j) {
    
    x1 = SIM$gt1[i,]
    x2 = SIM$gt2[i,]
    y1 = SIM$gt1[j,]
    y2 = SIM$gt2[j,]
    
    
    v = 2- abs( (x1+x2)- (y1+y2) )
    plot(v,t="l")
    
    table(v)
    sum(v==0)
    
    
}



