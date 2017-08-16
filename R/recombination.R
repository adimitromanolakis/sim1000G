
#' Generate recombination distances using a no-interference model.
#'
#'
#' @param n Number of distances to generate
#'
#' @return recombination distances in centimorgan
#' @examples
#'
#' library("sim1000G")
#' mean ( generateRecombinationDistances_noInterference ( 200 ) )
#'
#' @export
generateRecombinationDistances_noInterference = function ( n ) {

    ( rexp(n, 0.01))
}




#' Genetates a recombination vector arising from one meiotic event.
#' The origin of segments is coded as (0 - haplotype1 ,  1 - haplotype2 )
#'
#' @param cm The length of the region that we want to generate recombination distances.
#'
#' @examples
#'
#' library("sim1000G")
#'
#' examples_dir = system.file("examples", package = "sim1000G")
#' vcf_file = sprintf("%s/region.vcf.gz", examples_dir)
#' vcf = readVCF( vcf_file, maxNumberOfVariants = 100 , min_maf = 0.12 ,max_maf = NA)
#'
#' # For realistic data use the functions downloadGeneticMap / readGeneticMap
#' generateFakeGeneticMap()
#' generateSingleRecombinationVector( 1:100 )
#'
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













#' Holds the genetic map information that is used for simulations.
#' @export
geneticMap <- new.env()




#' Downloads a genetic map for a particular chromosome under GRCh37 coordinates for use with sim1000G.
#'
#' @param chromosome Chromosome number to download recombination distances from.
#'
#'
#' @examples
#'
#' \dontrun{
#'
#' downloadGeneticMap(22)
#'
#' }
#' @export
downloadGeneticMap = function(chromosome) {


        fname = sprintf("genetic_map_GRCh37_chr%s.txt.gz" , chromosome )

        url = sprintf(
            "https://github.com/adimitromanolakis/geneticMap-GRCh37/raw/master/genetic_map_GRCh37_chr%s.txt.gz",
            chromosome
        )

        dest_dir = system.file("data", fname, package = "sim1000G")
        dest_dir = "."

        cat("Downloading genetic map\n")
        dest_dir = system.file("data", package = "sim1000G")

        dest_dir = "./"
        dest_path = sprintf("%s/genetic_map_GRCh37_chr4.txt.gz", dest_dir)
        cat(" -> Saving genetic map to: " , dest_path,"\n")
        #file.exists(dest_path)

        download.file(url  ,  destfile = dest_path, quiet=TRUE)


}




#' Reads a genetic map downloaded from the function downloadGeneticMap.
#'
#' The map must contains a complete chromosome or enough markers to cover the area that
#' will be simulated.
#'
#' The file can be downloaded from the provided genetic maps by using the function downloadGeneticMap.
#'
#'
#' @param chromosome Chromosome number to download recombination distances from.
#' @param dir Directory the map file is located.
#'
#'
#'
#' @examples
#'
#' \dontrun{
#'
#' downloadGeneticMap(22)
#' readGeneticMap(22)
#'
#' }
#'
#' @export
readGeneticMap = function(chromosome, dir=".") {

    filelocation = sprintf( "%s/genetic_map_GRCh37_chr%s.txt.gz", dir, as.character( chromosome) )
    readGeneticMapFromFile(filelocation)

}




#' Reads a genetic map to be used for simulations.
#'
#' The file must be contain the following columns in order: chromosome, basepaire, rate(not used), centimorgan
#'
#'
#' @param filelocation Filename containing the genetic map
#' @examples
#'
#' \dontrun{
#'
#'
#' readGeneticMapFromFile("genetic_map_GRCh37_chr4.txt.gz")
#'
#' }
#' @export
readGeneticMapFromFile = function(filelocation) {

    if(! file.exists(filelocation) ) {
        stop(sprintf("Genetic map file %s not found",filelocation))
    }

    a = read.table(filelocation, as.is=TRUE, header=TRUE)
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



#' Generates a fake genetic map.
#'
#'
#' Generates a fake genetic map by approximating 1 cm / Mbp. Only used for examples.
#'
#'
#' @examples
#'
#' library("sim1000G")
#'
#' examples_dir = system.file("examples", package = "sim1000G")
#' vcf_file = sprintf("%s/region.vcf.gz", examples_dir)
#' vcf = readVCF( vcf_file, maxNumberOfVariants = 100 , min_maf = 0.12 ,max_maf = NA)
#'
#' # For realistic data use the functions downloadGeneticMap / readGeneticMap
#' generateFakeGeneticMap()
#'
#' @export
generateFakeGeneticMap = function() {

    n = 10000
    bp = seq(0,300e6, by=5000)
    cm = bp / 1e6

    geneticMap$chr = rep(4, length(bp))
    geneticMap$bp = bp
    geneticMap$cm = cm

    makeCDF();
    0;

}




#' Converts centimorgan position to base-pair. Return a list of centimorgan positions that correspond
#' to the bp vector (in basepairs).
#'
#' @param bp vector of base-pair positions
#'
#' @examples
#'
#' library("sim1000G")
#'
#' examples_dir = system.file("examples", package = "sim1000G")
#' vcf_file = sprintf("%s/region.vcf.gz", examples_dir)
#' vcf = readVCF( vcf_file, maxNumberOfVariants = 100 , min_maf = 0.12 ,max_maf = NA)
#'
#' # For realistic data use the functions downloadGeneticMap / readGeneticMap
#' generateFakeGeneticMap()
#' getCMfromBP(seq(1e6,100e6,by=1e6))
#'
#'
#' @export
getCMfromBP = function(bp) {
    approx( geneticMap$bp, geneticMap$cm, bp )$y
}



#' Generates a plot of the genetic map for a specified region.
#'
#' The plot shows the centimorgan vs base-pair positions.
#' The position of markers that have been read is also depicted as vertical lines
#'
#'
#' @param bp Vector of base-pair positions to generate a plot for
#' library("sim1000G")
#'
#' examples_dir = system.file("examples", package = "sim1000G")
#' vcf_file = sprintf("%s/region.vcf.gz", examples_dir)
#' vcf = readVCF( vcf_file, maxNumberOfVariants = 100 , min_maf = 0.12 ,max_maf = NA)
#'
#' # For realistic data use the functions downloadGeneticMap / readGeneticMap
#' generateFakeGeneticMap()
#' plotRegionalGeneticMap(seq(1e6,100e6,by=1e6/2))
#'
#'
#' @export
plotRegionalGeneticMap = function(bp) {

    #bp = vcf[,2]
    cm = approx( geneticMap$bp, geneticMap$cm, bp )$y

    plot(bp/1e6,cm, t="l",cex=0.2, xlab="MBp", ylab = "Centimorgan")
    segments( bp/1e6,cm,  bp/1e6,cm -0.2, lwd=0.3,col="blue")


}






#' Contains recombination model information.
#'
#' This vector contains the density between two recombination events, as a cumulative density function.
##' @export
crossoverCDFvector = NA


fgammamod = function(x,L, n) {

    (x**(n-1))*exp(-L*x)    *  (L**n)/gamma(n)
}


fstar = function(x,q,l) {


    sum1 = 0
    for(k in 1:5)
        sum1 = sum1 + fgammamod(x, 2*q*(l+1) , k*(l+1))/ (2**k)

    sum1

}





#plot(function(x) fstar(x,1,4.5),xlim=c(0,3), xlab="Morgans", ylab="Density")
#title("Crossover distance")
#abline(h=1,col="gray")

crossoverCDF = new.env()

makeCDF = function() {

    dens = sapply(seq(0,4,by=0.001), function(x) fstar(x,1,1) )  # originally was 3.2


    cdf = cumsum(dens)
    cdf = cdf/max(cdf)

    cdf[1:10]
    crossoverCDF$vector = cdf
    crossoverCDFvector

}







#' Generate recombination distances using a chi-square model.
#'
#'
#' @param n Number of distances to generate
#'
#' @return vector of recombination distances in centimorgan
#'
#' @export
generateRecombinationDistances = function ( n ) {


    t = findInterval( runif(n), crossoverCDF$vector )
    t = t * 4 / length(crossoverCDF$vector)
    100 * t
}


#hist(generateRecombinationDistances(20000),n=100)


