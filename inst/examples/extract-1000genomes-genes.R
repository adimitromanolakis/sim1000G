
## Extract vcf files of genes from selected population in 1000 genomes
##
## Needs:
##   20130606_g1k.ped:  population information file
##   genes.txt: gene list file
##   complete 1000 genomes phased vcf files:
##           files: ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz
##   bcftools version 1.3.1

chrom = commandArgs(T)[1]
pops = commandArgs(T)[2]
output_directory = commandArgs(T)[3]

print(commandArgs(T))

if(length(commandArgs(T))<3) stop("Usage: script.R chrom pops output_directory")

pops = strsplit(pops,",")[[1]]
print(pops)



ped_file = "~/1000genomes/20130606_g1k.ped"

which_populations = c("ASW","LWK","YRI")
which_populations = c("CEU","TSI","GBR")
which_populations = pops

gene_expand_bp = 5e3


num_genes_to_extract = 400

run_commands_within_r = 1


sample_subset_file = tempfile()
##


cat("Chromosome ", chrom, "\n")

if(is.na(chrom)) { stop(" no chromosome provided in command line ")}



ped = read.table(ped_file,h=T,as=T,sep="\t")
pop = ped$Population
names(pop) = ped$Individual.ID


id1 = ped$Individual.ID [ ped$Population %in% which_populations]
cat(id1,file=   sample_subset_file ,  sep="\n")

cat("individuals selected: ", length(id1),"\n")

cat("individuals selected: ",id1,"\n")
cat("file=",sample_subset_file,"\n")



# Read genes table

# Format:
#     V1    V2    V3     V4 V5 V6      V7   V8      V9     V10        V11
#      1 11873 14408 strand  1 id DDX11L1 name DDX11L1 bioType pseudogene


genes_table = read.table("~/fs/genes.txt",h=F,as=T)

s = genes_table$V11 == "protein_coding" & genes_table$V1 == chrom

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

    cat("RUN: ",cmd,"\n")

    cmd = sprintf("echo %s\n%s\n",name,cmd)



    if(run_commands_within_r) system(cmd)

    cmd
}









f = function(i) {
    start_bp = genes_table$V2[i]  - gene_expand_bp
    end_bp = genes_table$V3[i]  + gene_expand_bp
    name = genes_table$V7[i]



    print(genes_table[i,])
    extract_commands(i, chrom, start_bp, end_bp, name, output_directory)

}





f_2mbp = function(i) {
    start_bp = 2e6 * i
    end_bp = 2e6*(i+1)
    name = sprintf("r%.1fMBp",i)



    extract_commands(i, chrom, start_bp, end_bp, name, output_directory)

}



#f = f_2mbp




i = 10


start = genes_table$V2[i]
end = genes_table$V3[i]
name = genes_table$V7[i]





fn = "~/1000genomes/ALL.chr%s.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"





cmd = " bcftools view -r %s -S FILE1 --force-samples %s | bcftools filter -e 'AF==0' | ~/truffle/compress --cpu2 10 > %s"
cmd = sub("FILE1",sample_subset_file, cmd)
print(cmd)
#exit()

outfile_template = "region-chr%s-%d-%s.vcf.gz"



if(num_genes_to_extract == -1) num_genes_to_extract = nrow(genes_table)

cmds2 = sapply(1:num_genes_to_extract,f)




cat(cmds2,file="extract-genes.sh",sep="\n")



cat("Extracted the following populations:\n");
table( ped$Population [ ped$Population %in% which_populations] )

