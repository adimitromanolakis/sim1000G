
#' @export
writePED = function(fam) {
    
    
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
        write.table(fam_header,file="out.ped",sep=" ",row=F,col=F,quote=F)
        
        dim(fam)
        dim(vcf)
        
        map = data.frame( vcf$vcf[,1] , vcf$vcf[,3] , 0 , vcf$vcf[,2])
        write.table(map,file="out.map",sep=" ",row=F,col=F,quote=F)
        system("~/tools/prest-plus/prest --geno out.ped --map out.map --wped --ibs1 0.5")

}





