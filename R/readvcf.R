

#vcf <- new.env()


#' Creates a regional vcf file using bcftools to extract a region from 1000 genomes vcf files
#'
#'
#' @return none
#' @export
createVCF = function() {

    # node ../fq-thin-vcf-file.js --out vt.vcf --in chr2.vcf --thin 0.7
    # bcftools view -s NA19904,NA19913,NA20317,NA20318,NA20334,NA20355,NA20359,NA20362 vt.vcf|grep -v '^##' > vt2.vcf

    #4 90645249 90759446 strand -1 id SNCA name SNCA bioType protein_coding
    '



CHROM=4
VCF=~/1000genomes/ALL.chr"$CHROM".*.vcf.gz
REGION=4:90645249-90759446
REGION=4:60995249-61569446


# 4 77356252 77704404 strand 1 id SHROOM3 name SHROOM3 bioType protein_coding
#REGION=4:77356252-77704404
#REGION=5:77356252-77704404



bcftools view -S /tmp/CEU.txt --force-samples -r $REGION  $VCF > /tmp/1.vcf
grep -v "^##" /tmp/1.vcf > /tmp/2.vcf
cp /tmp/2.vcf ~/tmp



'

}



#' Read a vcf file, with options to filter out low or high frequency markers.
#'
#'
#' @param filename Input VCF file
#' @param thin How much to thin markers
#' @param maxNumberOfVariants Maximum number of variants to keep from region
#' @param min_maf Minimum allele frequency of markers to keep. If NA skip min_maf filtering.
#' @param max_maf Maximum allele frequency of markers to keep. If NA skip max_maf filtering.
#' @param region_start Extract a region from a vcf files with this starting basepair position
#' @param region_end Extract a region from a vcf files with this ending basepair position
#'
#' @return VCF object to be used by startSimulation function.
#'
#' @examples
#'  \dontrun{
#' library(sim1000G)
#` vcf <- readVCF("region.vcf.gz", maxNumberOfVariants = 100, min_maf = 0.02, max_maf = 0.1)
#` readGeneticMap(chromosome = 4)
#` startSimulation(vcf, totalNumberOfIndividuals = 500)
#'  }
#' @export
readVCF = function(filename = "data.vcf",
                   thin = NA,
                   maxNumberOfVariants = 400,
                   min_maf = 0.02,
                   max_maf = NA,
                   region_start = NA,
                   region_end = NA
                   )
{


    if(  grepl(".rdata", filename) ) {

        cat("[#.......] Loading VCF environment from rdata file\n")

        e1 = new.env()
        load(filename, envir=e1)

        vcf = e1$vcf
        cat("[##......] Number of variants in VCF: ", nrow(vcf$vcf) , "\n");
        cat("[##......] Number of individuals in VCF: ", ncol(vcf$gt1) , "\n");

        return(vcf);


    }




    cat("[#.......] Reading VCF file..\n")



    # vcf = read.table(filename, sep="\t", comment.char="", as.is=T, header=T)

    vcf = read_delim(filename, "\t", comment = "##", progress = FALSE)
    vcf = as.data.frame(vcf)






    # Keep bi-allelic markers only
    vcf = vcf[str_length(vcf[,5]) == 1 ,]


    if( !is.na( region_start ) ) vcf = vcf[ vcf[,2] >= region_start  , ]
    if( !is.na( region_end   ) ) vcf = vcf[ vcf[,2] <= region_end    , ]

    if( nrow(vcf) < 10 ) {
         cat("Too few variants in VCF file\n");
         return(NA);
    }


    # Reduce number of variants
    if( !is.na(thin) )
        vcf = vcf[seq(1,nrow(vcf),by = thin) ,]


    cat("[##......] Chromosome:  ", unique(vcf[,1]),
        " Mbp: " , min(vcf[,2])/1e6,
        " Region Size: ", 1e-3 * ( max(vcf[,2])-min(vcf[,2]) ) ,"kb ",
        "Num of variants:", dim(vcf)[1] ,
        "Num of individuals:", ncol(vcf)-9,
        "\n");


    individual_ids = colnames(vcf)[-(1:9)]
    gt = vcf[,-(1:9) ]
    vcf = vcf[,1:9]


    gt1 = apply( gt, 2,  function(x) as.numeric( str_sub(x,1,1) )  )
    gt2 = apply( gt, 2,  function(x) as.numeric( str_sub(x,3,3) )  )



    #### --- Filter by MAF ---- ####


    cat("[###.....] Filtering and thinning variants\n");

    maf = apply(gt1+gt2,1,function(x) mean(x,na.rm=T)/2)


    #maf[maf>0.5] = 1 - maf[maf>0.5]
    #maf2 = apply(gt2,1,function(x) mean(x,na.rm=T))
    #maf2[maf2>0.5] = 1 - maf2[maf2>0.5]

    ok = apply(gt1,1,function(x) max(x,na.rm=T))
    s = which(ok < 2)


    if( ! is.na(min_maf) ) {  s = intersect( s  ,   which(maf >= min_maf & maf <= 1 - min_maf ) ) }

    if( ! is.na(max_maf) ) {  s = intersect( s  ,   which(maf <= max_maf | maf >= 1 - max_maf) ) }

    if(length(s) > maxNumberOfVariants) s = sort( sample(s,maxNumberOfVariants) )




    gt1 = gt1[s,]
    gt2 = gt2[s,]
    gt = gt[s,]
    vcf = vcf[s,]
    maf = maf[s]




    # cat("Q: flip genotypes if maf > 0.5???????????\n")


    cat("[##......] Chromosome:  ", unique(vcf[,1]),
        " Mbp: " , min(vcf[,2])/1e6,
        " Region Size: ", 1e-3 * ( max(vcf[,2])-min(vcf[,2]) ) ,"kb ",
        "Num of variants:", dim(vcf)[1] , "(after filtering)\n")


    ##


    R = new.env()

    R$vcf = vcf
    R$gt1 = gt1
    R$gt2 = gt2
    R$individual_ids = individual_ids
    R$maf = maf
    # R$gt = gt

    R

}
