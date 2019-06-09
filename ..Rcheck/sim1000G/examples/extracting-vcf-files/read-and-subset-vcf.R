library(sim1000G)
#library(gplots)


# Read and subset a large region vcf file

vcf_file = "region-chr4-2-r2.0MBp.vcf.gz"

vcf = readVCF( vcf_file, maxNumberOfVariants = 1e6 , min_maf = 1e-6 ,max_maf = 0.01)

plot(vcf$maf)


subsetVCF = function(vcf, var_index = NA) {


    vcf2 = as.list(vcf)

    s = var_index

    vcf2$maf = vcf2$maf [s]
    vcf2$varid = vcf2$varid [s]

    vcf2$gt1 = vcf2$gt1 [s,]
    vcf2$gt2 = vcf2$gt2 [s,]

    vcf2$vcf = vcf2$vcf [s,]

    as.environment(vcf2)


}


str(as.list(vcf))



# find number of variants in the vcf file
num_variants = nrow(vcf$gt1)
num_individuals = ncol(vcf$gt1)


# select first 300 variants

vcf2 = subsetVCF(vcf, 1:400)

str(as.list(vcf2))


dim(vcf$gt1)
dim(vcf2$gt1)

