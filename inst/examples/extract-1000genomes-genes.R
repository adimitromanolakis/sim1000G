
## Extract vcf files of genes from selected population in 1000 genomes
##
## Needs:
##   20130606_g1k.ped:  population information file
##   genes.txt: gene list file
##   complete 1000 genomes phased vcf files:
##           files: ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
##   bcftools version 1.3.1

chrom = 4

gene_expand_bp = 250e3


# Read pedigree files from 1000 genomes

ped = read.table("20130606_g1k.ped",h=T,as=T,sep="\t")

pop = ped$Population
names(pop) = ped$Individual.ID

which_populations = c("ASW","LWK","YRI")

    id1 = ped$Individual.ID [ ped$Population %in% which_populations]
    cat(id1,file="sample_subset_afr.txt",sep="\n")


which_populations = c("CEU","TSI","GBR")

    id1 = ped$Individual.ID [ ped$Population %in% which_populations]
    cat(id1,file="sample_subset_eur.txt",sep="\n")






# Read genes table

# Format:
#     V1    V2    V3     V4 V5 V6      V7   V8      V9     V10        V11
#      1 11873 14408 strand  1 id DDX11L1 name DDX11L1 bioType pseudogene


genes_table = read.table("~/fs/genes.txt",h=F,as=T)



s = genes_table$V11 == "protein_coding" & genes_table$V1 == "4"

genes_table = genes_table[s,]



str(genes_table)
table(genes_table$V11)
table(genes_table$V1)





extract_commands = function(num, chrom, start,end,name, output_directory=".") {




    region = sprintf("%s:%d-%d",chrom,start,end)

    cat(region, name, "size=", (end-start)/1000, "kbp\n");

    outfile = sprintf(outfile_template, chrom, num, name);

    outfile = file.path(output_directory, outfile)

    fn = sprintf(fn, chrom)

    cmd = sprintf(cmd, region, fn, outfile)




    cmd = sprintf("echo %s\n%s\n",name,cmd)

    system(cmd)



    cmd


}





i = 10


start = genes_table$V2[i]
end = genes_table$V3[i]
name = genes_table$V7[i]





fn = "ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"







f = function(i) {
    start_bp = genes_table$V2[i]  - gene_expand_bp
    end_bp = genes_table$V3[i]  + gene_expand_bp
    name = genes_table$V7[i]



    extract_commands(i, chrom, start_bp, end_bp, name, output_directory)

}





f_2mbp = function(i) {
    start_bp = 2e6 * i
    end_bp = 2e6*(i+1)
    name = sprintf("r%.1fMBp",i)



    extract_commands(i, chrom, start_bp, end_bp, name, output_directory)

}


# f = f_2mbp




cmd = " bcftools view -r %s -S sample_subset_afr.txt --force-samples %s | bcftools filter -e 'AF==0' | bgzip > %s"
outfile_template = "region-chr%s-%d-%s.vcf.gz"
output_directory ="/tmp/afr"

cmds = sapply(1:500,f)


cmd = " bcftools view -r %s -S sample_subset_eur.txt --force-samples %s | bcftools filter -e 'AF==0' | bgzip > %s"
outfile_template = "region-chr%s-%d-%s.vcf.gz"
output_directory ="/tmp/eur"

cmds2 = sapply(1:500,f)




cat(cmds,cmds2,file="extract-genes.sh",sep="\n")



cat("Extracted the following populations:\n");
table( ped$Population [ ped$Population %in% which_populations] )

