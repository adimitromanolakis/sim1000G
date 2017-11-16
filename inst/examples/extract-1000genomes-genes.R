


# Find population IDs to extract from 1000 genomes
ped = read.table("20130606_g1k.ped",h=T,as=T,sep="\t")
which_populations = c("CEU","TSI","ASW")
#which_populations = c("CEU","TSI","GBR")
which_populations = c("ASW","LWK","YRI")

outfile_template = "ASW-LWK-YRI-region-chr%s-%d-%s.vcf.gz"



id1 = ped$Individual.ID [ ped$Population %in% which_populations]
id1

cat(id1,file="sample_subset.txt",sep="\n")


genes_table = read.table("~/fs/genes.txt",h=F,as=T)


chrom = 4

s = genes_table$V11 == "protein_coding" & genes_table$V1 == "4"
genes_table = genes_table[s,]

str(genes_table)
table(genes_table$V11)
table(genes_table$V1)




fn = "ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"
#fn = "chr%s-ceutsigbr.vcf"
cmd = " bcftools view -r %s -S sample_subset.txt --force-samples %s | bcftools filter -e 'AF==0' | bgzip > %s"




extract_commands = function(num, chrom, start,end,name, output_directory=".") {




    region = sprintf("%s:%d-%d",chrom,start,end)

    cat(region, name, "size=", (end-start)/1000, "kbp\n");

    outfile = sprintf(outfile_template, chrom, num, name);

    outfile = file.path(output_directory, outfile)

    fn = sprintf(fn, chrom)

    cmd = sprintf(cmd, region, fn, outfile)


    cmd = sprintf("echo %s\n%s\n",name,cmd)
    cmd


}





i = 10


start = genes_table$V2[i]
end = genes_table$V3[i]
name = genes_table$V7[i]


extract_commands(i, chrom, start, end, name)

f = function(i) {
    start = genes_table$V2[i]
    end = genes_table$V3[i]
    name = genes_table$V7[i]
    extract_commands(i, chrom, start, end, name, output_directory ="/tmp")



}

dim(genes_table)

cmds = sapply(1:500,f)

cat(cmds,file="extract-genes.sh",sep="\n")

cat("Extracted the following populations:\n");
table( ped$Population [ ped$Population %in% which_populations] )

