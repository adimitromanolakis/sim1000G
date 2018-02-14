
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
#' vcf_file = file.path(examples_dir, "region.vcf.gz")
#' vcf = readVCF( vcf_file, maxNumberOfVariants = 100 , min_maf = 0.12 ,max_maf = NA)
#'
#' # For realistic data use the functions downloadGeneticMap / readGeneticMap
#' generateUniformGeneticMap()
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







#' Generates a recombination vector arising from one meiotic event.
#' The origin of segments is coded as (0 - haplotype1 ,  1 - haplotype2 )
#'
#' @param chromosomeLength The length of the region in cm.
#'
#' @examples
#'
#' library("sim1000G")
#'
#' # generate a recombination events for chromosome 4
#' readGeneticMap(4)
#' generateChromosomeRecombinationPositions(500)
#'
#' @export
generateChromosomeRecombinationPositions = function(chromosomeLength = 500) {

    pos = generateRecombinationDistances(10)

    while(sum(pos) < chromosomeLength) pos = c(pos, generateRecombinationDistances(10) )

    pos = cumsum(pos)

    pos = pos[pos <= chromosomeLength]
    pos
}











#' Holds the genetic map information that is used for simulations.
#' @export
geneticMap <- new.env()




#' Downloads a genetic map for a particular chromosome under GRCh37 coordinates for use with sim1000G.
#'
#' @param chromosome Chromosome number to download recombination distances from.
#' @param dir Directory to save the genetic map to (default: extdata)
#'
#'
#' @examples
#'
#'
#'
#' downloadGeneticMap(22, dir=tempdir() )
#'
#'
#' @export
downloadGeneticMap = function(chromosome, dir = NA) {


        fname = sprintf("genetic_map_GRCh37_chr%s.txt.gz" , chromosome )

        url = sprintf(
            "https://github.com/adimitromanolakis/geneticMap-GRCh37/raw/master/genetic_map_GRCh37_chr%s.txt.gz",
            chromosome
        )


        if(is.na(dir)) {
            dest_dir = system.file("datasets", package = "sim1000G")

            if(dest_dir == "") dest_dir = tempdir()
            if(!dir.exists(dest_dir))  dest_dir = tempdir()
        } else {
            dest_dir = dir
        }


        #dest_dir = "./"
        dest_path = file.path(dest_dir, fname)


        if(file.exists(dest_path)) {
            cat(" -> Using genetic map found on ", dest_path, "\n");
            return(dest_path)

        }

        cat(" -> Downloading genetic map from:",url,"\n")
        cat(" -> Saving genetic map to: " , dest_path,"\n")

        download.file(url  ,  destfile = dest_path, quiet=TRUE)
        return(dest_path)

}





#' Reads a genetic map downloaded from the function downloadGeneticMap or reads a genetic map from a specified file. If the argument filename is used
#' then the genetic map is read from the corresponding file. Otherwise, if a chromosome is specified, the genetic map is downloaded for human chromosome
#' using grch37 coordinates.
#'
#' The map must contains a complete chromosome or enough markers to cover the area that
#' will be simulated.
#'
#'
#' @param chromosome Chromosome number to download a genetic map for , or
#' @param filename A filename of an existing genetic map to read from (default NA).
#' @param dir Directory the map file will be saved (only if chromosome is specified).
#'
#'
#'
#' @examples
#'
#'
#'
#'
#' readGeneticMap(chromosome = 22)
#'
#'
#'
#' @export
readGeneticMap = function(chromosome, filename=NA, dir=NA) {


    if( ! is.na(filename) ) { return ( readGeneticMapFromFile(filename) ) }
    #if(is.na(dir)) dir = system.file("extdata", package = "sim1000G")

    fname =   downloadGeneticMap( chromosome ,dir=dir)


    readGeneticMapFromFile(fname)

}




#' Reads a genetic map to be used for simulations. The genetic map should be
#' of a single chromosome and covering the extent of the region to be simulated.
#' Whole chromosome genetic maps can also be used.
#'
#' The file must be contain the following columns in the same order: chromosome, basepaire, rate(not used), centimorgan
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

    makeCDF()
}



#' Generates a uniform genetic map.
#'
#'
#' Generates a uniform genetic map by approximating 1 cm / Mbp. Only used for examples.
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
#' generateUniformGeneticMap()
#' plotRegionalGeneticMap(seq(1e6,100e6,by=1e6/2))
#'
#' @export
generateUniformGeneticMap = function() {

    n = 10000
    bp = seq(0,300e6, by=5000)
    cm = bp / 1e6

    geneticMap$chr = rep(4, length(bp))
    geneticMap$bp = bp
    geneticMap$cm = cm

    makeCDF();
    0;

}



#' Generates a fake genetic map that spans the whole genome.
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
#' generateFakeWholeGenomeGeneticMap(vcf)
#' plotRegionalGeneticMap(seq(1e6,100e6,by=1e6/2))
#'
#' @param vcf A vcf file read by function readVCF.
#' @export
generateFakeWholeGenomeGeneticMap = function(vcf) {

    n = 10000
    bp = vcf$vcf[,2]
    cm = seq(0,4000, l=length(bp))

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
#' generateUniformGeneticMap()
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
#' generateUniformGeneticMap()
#' plotRegionalGeneticMap(seq(1e6,100e6,by=1e6/2))
#'
#'
#' @export
plotRegionalGeneticMap = function(bp) {

    #bp = vcf[,2]
    cm = approx( geneticMap$bp, geneticMap$cm, bp )$y


    len = diff(range(cm))
    print(len)

    plot(bp/1e6,cm, t="l",cex=0.2, xlab="MBp", ylab = "Centimorgan")
    segments( bp/1e6,cm,  bp/1e6,cm - len/10, lwd=0.3,col="blue")


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







#' Generate inter-recombination distances using a chi-square model. Note this are the distances between two succesive recombination events and not
#' the absolute positions of the events. To generate the locations of the recombination events see the example below.
#'
#'
#' @param n Number of distances to generate
#'
#' @return vector of distances between two recombination events.
#'
#' @examples
#'
#' library("sim1000G")
#'
#' distances = generateRecombinationDistances(20)
#'
#'
#' positions_of_recombination = cumsum(distances)
#'
#' if(0) hist(generateRecombinationDistances(20000),n=100)
#'
#' @export
generateRecombinationDistances = function ( n ) {


    if(pkg.opts$recombination == 0) { return(generateRecombinationDistances_noInterference(n))}


    t = findInterval( runif(n), crossoverCDF$vector )
    t = t * 4 / length(crossoverCDF$vector)
    100 * t
}




