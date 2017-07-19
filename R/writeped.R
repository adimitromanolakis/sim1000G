
#' @export
writePED = function(fam, filename = "out" ) {
    
    
        filename_ped = sprintf("%s.ped", filename)
        filename_map = sprintf("%s.map", filename)
        
    
        N = SIM$individuals_generated 
        fam[1,]
        
        allele1 = SIM$gt1[fam$gtindex ,]
        allele2 = SIM$gt2[fam$gtindex ,]
        

        getGenotypeFromHaplo = function(i) {
            paste( allele1[i,]+1, allele2[i,]+1, sep=" ")
        }
        
        genotypes = sapply(1:N, getGenotypeFromHaplo)
        
        
        ids = sprintf("IND%s", seq(1,N))
        fam_header = data.frame( fam[,1:4], fam[,5],"1\t",  t(genotypes)) 
        write.table(fam_header, file  =  filename_ped , sep=" ", row.names=F, col.names=F, quote=F)
        
        dim(fam)
        dim(vcf)
        
        map = data.frame( vcf$vcf[,1] , vcf$vcf[,3] , 0 , vcf$vcf[,2])
        write.table(map,   filename_map ,   sep=" ",  row.names=F,  col.names=F,  quote=F)
        
        cat("[] PED file written as ", filename_ped, "\n");
        
        # system("~/tools/prest-plus/prest --geno out.ped --map out.map --wped --ibs1 0.5")

}





