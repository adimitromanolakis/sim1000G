

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
downloadGeneticMap = function(chromosome) {
 
    
    fname = sprintf("genetic_map_GRCh37_chr%s.txt.gz" , chromosome )
    
    url = sprintf("https://github.com/adimitromanolakis/geneticMap-GRCh37/raw/master/genetic_map_GRCh37_chr4.txt.gz")
    
    dest_dir = system.file("data", fname, package = "sim1000G")
    dest_dir
    
     cat("Downloading genetic map\n")
        dest_dir = system.file("data", package = "sim1000G")
        
        dest_dir = "./"
        dest_path = sprintf("%s/genetic_map_GRCh37_chr4.txt.gz", dest_dir)
        cat(dest_path,"\n")
        file.exists(dest_path)
        
        download.file(url  ,  dest = dest_path)
        

}





#' @export
readGeneticMap = function(chromosome, dir=".") {
    
    filelocation = sprintf( "%s/genetic_map_GRCh37_chr%s.txt.gz", dir, as.character( chromosome) ) 
    
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
    
    
    cat("      -> Genetic map has" , length(geneticMap$bp), "entries\n");
    
    makeCDF();
    0;
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







#' @export
crossoverCDFvector = NA


#' @export
fgammamod = function(x,L, n) {
    
    (x**(n-1))*exp(-L*x)    *  (L**n)/gamma(n) 
}


#' @export
fstar = function(x,q,l) {
    
    
    sum1 = 0
    for(k in 1:5) 
        sum1 = sum1 + fgammamod(x, 2*q*(l+1) , k*(l+1))/ (2**k)
    
    sum1
    
}





#plot(function(x) fstar(x,1,4.5),xlim=c(0,3), xlab="Morgans", ylab="Density")
#title("Crossover distance")
#abline(h=1,col="gray")

#' @export
crossoverCDF = new.env()

#' @export
makeCDF = function() {
    
    dens = sapply(seq(0,4,by=0.001), function(x) fstar(x,1,1) )  # originally was 3.2
    
    
    cdf = cumsum(dens)
    cdf = cdf/max(cdf)
    
    cdf[1:10]
    crossoverCDF$vector = cdf
    crossoverCDFvector
    
}








#' @export
generateRecombinationDistances = function ( n ) {
    
    
    t = findInterval( runif(n), crossoverCDF$vector )
    t = t * 4 / length(crossoverCDF$vector)
    100 * t
}


#hist(generateRecombinationDistances(20000),n=100)













