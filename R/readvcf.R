

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
REGION=4:60645249-61569446


# 4 77356252 77704404 strand 1 id SHROOM3 name SHROOM3 bioType protein_coding
REGION=4:77356252-77704404



bcftools view -S CEU.txt --force-samples -r $REGION  $VCF > /tmp/1.vcf
grep -v "^##" /tmp/1.vcf > haplosims/2.vcf



'

}



#' Read a vcf file, with options to filter out low or high frequency markers.
#'
#'
#' @param filename Input VCF file
#' @param thin How much to thin markers
#' @param maxNumberOfVariants Maximum number of variants to keep from region
#' @param min_maf Minimum allele frequency of markers to keep
#' @param max_maf Maximum allele frequency of markers to keep
#'
#' @return none
#'
#' @export
readVCF = function(filename = "haplosims/1.vcf", thin = 1, maxNumberOfVariants = 400, min_maf = 0.02, max_maf = NA) {


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



    vcf = read.table(filename, sep="\t", comment.char="", as.is=T, header=T)


    # Keep bi-allelic markers only
    vcf = vcf[str_length(vcf[,5]) == 1 ,]


    # Reduce number of variants
    if(thin > 0)
    vcf = vcf[seq(1,nrow(vcf),by = thin) ,]


    cat("[##......] Chromosome:  ", unique(vcf[,1]),
        " Mbp: " , min(vcf[,2])/1e6,
        " Region Size: ", 1e-3 * ( max(vcf[,2])-min(vcf[,2]) ) ,"kb ",
        "Num of variants:", dim(vcf)[1] ,
        "\n");


    gt = vcf[,-(1:9) ]

    gt1 = apply( gt, 2,  function(x) as.numeric( str_sub(x,1,1) )  )
    gt2 = apply( gt, 2,  function(x) as.numeric( str_sub(x,3,3) )  )


    #### --- Filter by MAF ---- ####


    cat("[###.....] Filtering and thinning variants\n");

    maf = apply(gt1,1,function(x) mean(x,na.rm=T))
    maf[maf>0.5] = 1 - maf[maf>0.5]


    maf2 = apply(gt2,1,function(x) mean(x,na.rm=T))
    maf2[maf2>0.5] = 1 - maf2[maf2>0.5]

    ok = apply(gt1,1,function(x) max(x,na.rm=T))
    s = which(maf > min_maf & maf2 > min_maf & ok < 2)

    if( ! is.na(max_maf) ) {  s = intersect( s  ,   which(maf <= max_maf & maf2 <= max_maf & ok < 2) ) }

    if(length(s)>maxNumberOfVariants) s = sort( sample(s,maxNumberOfVariants) )


    gt1 = gt1[s,]
    gt2 = gt2[s,]
    gt = gt[s,]
    vcf = vcf[s,]

    #


    cat("[##......] Chromosome:  ", unique(vcf[,1]),
        " Mbp: " , min(vcf[,2])/1e6,
        " Region Size: ", 1e-3 * ( max(vcf[,2])-min(vcf[,2]) ) ,"kb ",
        "Num of variants:", dim(vcf)[1] , "(after filtering)\n")


    ##


    R = new.env()

    R$vcf = vcf
    R$gt1 = gt1
    R$gt2 = gt2
    # R$gt = gt

    R

}
