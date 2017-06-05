

#' @export
generateRecombinationDistances = function ( n ) {
    
    ( rexp(n, 0.01))
}




# Genetates a recombination vector of origin of segments (0 - haplotype1 ,  1 - haplotype2 )

#' @export
generateSingleRecombinationVector = function(cm) {
    
    n = length(cm)
    maxCM = max(cm)
    
    pos = generateRecombinationDistances(14)
    
    while(sum(pos) < maxCM) pos = c(pos, generateRecombinationDistances(10) )
    if(sum(pos) < maxCM) { stop("Not enough recombination events"); }
    
    pos = cumsum(pos)
    
    recombVector = rep(TRUE, n)
    
    p = findInterval(pos,cm)
    
    
    for(i in 1:length(pos)) {
        
        if(p[i] < n) recombVector[ p[i]:n ] =  !recombVector[ p[i]:n ]  
        
        #s = cm > pos[i] 
        #if(runif(1) < 0.5) 
        #recombVector [ s ]  = recombVector [ s ]  + 1
        
    }
    
    #  recombVector = recombVector %% 2
    if(runif(1,0,1) < 0.5) recombVector = ! recombVector
    
    recombVector+0
    
}














#' @export

geneticMap <- new.env()


#' @export
readGeneticMap = function(chromosome) {
    
    filelocation = sprintf( "geneticmap/genetic_map_GRCh37_chr%s.txt.gz", as.character( chromosome) ) 
    
    if(! file.exists(filelocation) ) {
     stop(sprintf("Genetic map file %s not found",filelocation))   
    }
    
    a = read.table(filelocation, as=TRUE, h=TRUE)
    colnames(a) = c("chr","bp","rate.cm","cm")
    a = a[ order(a$cm),]
    a
    geneticMap$chr = a$chr
    geneticMap$bp = a$bp
    geneticMap$cm = a$cm
    
}



#' @export
getCMfromBP = function(bp) approx( geneticMap$bp, geneticMap$cm, bp )$y



#' @export
plotRegionalGeneticMap = function(bp) {
    
    #bp = vcf[,2]
    cm = approx( geneticMap$bp, geneticMap$cm, bp )$y
    
    plot(bp/1e6,cm, t="l",cex=0.2)
    segments( bp/1e6,cm,  bp/1e6,cm -0.2, lwd=0.3,col="blue")
    
    
}
