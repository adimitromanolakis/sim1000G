## ----results='asis',echo=F----------------------------------------------------
cat("<style>body { zoom: 1.0; }   .main-content pre> code {  white-space: pre-wrap; } </style>")

## ---- eval=FALSE--------------------------------------------------------------
#  install.packages("sim1000G")

## ----eval=T, message=FALSE, warning=FALSE, collapse=F, results='hide'---------

library(sim1000G)

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")

vcf = readVCF( vcf_file , maxNumberOfVariants = 100 , min_maf = 0.02 , max_maf = NA ) 

# downloadGeneticMap( 4 )
readGeneticMap( chromosome = 4 )

startSimulation( vcf )


## ----eval=T, echo=T, fig.height=12, fig.width=17, message=FALSE, warning=FALSE, results='hide'----



# The following function will erase all the indivudals that have been simulated and start over.
# SIM$reset()


id = c()
for(i in 1:30) id[i] = SIM$addUnrelatedIndividual()

# Show haplotype 1  of first 5 individuals
#print(SIM$gt1[1:5,1:6])

# Show haplotype 2
#print(SIM$gt1[1:5,1:6])



genotypes = SIM$gt1[1:20,] + SIM$gt2[1:20,]

print(dim(genotypes))

str(genotypes)

library(gplots)

heatmap.2(cor(genotypes)^2, col=rev(heat.colors(100)) , trace="none",Rowv=F,Colv=F)



## ---- echo=T, results='hold', eval=T   ,  results='hide'----------------------


# Simulate one family with 2 offspring

fam = newFamilyWithOffspring("fam1",2)
print(fam)



# Simulate 100 families
 
SIM$reset()


## For testing the IBD, we set the cM so that the regions spans to 4000cm
## Remove for normal use
SIM$cm = seq( 0,4000, length = length(SIM$cm) )



time10families = function() {
    
    
    fam = lapply(1:10, function(x) newFamilyWithOffspring(x,2) )
    fam = do.call(rbind, fam)
    fam

}

fam <- time10families() 



# Uncomment the following line to write a plink compatible ped/map file

# writePED(vcf, fam, "out")




## ----  results='hide'---------------------------------------------------------


n = SIM$individuals_generated

IBD1matrix = 
sapply(1:n, function(y) {
        z = sapply(1:n, function(x) computePairIBD12(x,y) [1]) 
        names(z) = 1:n
        z
})

IBD2matrix = 
    sapply(1:n, function(y) {
        z = sapply(1:n, function(x) computePairIBD12(x,y) [2]) 
        names(z) = 1:n
        z
    })








## ---- echo=T, results='asis'--------------------------------------------------
colnames(IBD1matrix) = 1:nrow(IBD1matrix)
rownames(IBD1matrix) = 1:nrow(IBD1matrix)
colnames(IBD2matrix) = 1:nrow(IBD2matrix)
rownames(IBD2matrix) = 1:nrow(IBD2matrix)

knitr::kable(IBD1matrix[1:8,1:8] )


## ---- echo=T, results='asis'--------------------------------------------------
colnames(IBD1matrix) = 1:nrow(IBD1matrix)
rownames(IBD1matrix) = 1:nrow(IBD1matrix)
colnames(IBD2matrix) = 1:nrow(IBD2matrix)
rownames(IBD2matrix) = 1:nrow(IBD2matrix)

knitr::kable(IBD2matrix[1:8,1:8] )


## ---- eval=F------------------------------------------------------------------
#  
#  newFamilyWithOffspring = function(familyid, noffspring = 2) {
#  
#      fam = data.frame(fid = familyid  ,
#                       id = c(1:2) ,
#                       father = c(0,0),
#                       mother = c(0,0),
#                       sex = c(1,2)
#      )
#  
#  
#      j1 = SIM$addUnrelatedIndividual()
#      j2 = SIM$addUnrelatedIndividual()
#  
#      fam$gtindex = c(j1,j2) # Holds the genotype position in the arrays SIM$gt1 and SIM$gt2
#  
#      for(i in 1:noffspring) {
#          j3 = SIM$mate(j1,j2)
#  
#          newFamilyMember = c(familyid, i+10, 1,2, 1 , j3)
#          fam = rbind(fam, newFamilyMember)
#      }
#  
#      return (fam)
#  }
#  
#  

## ---- echo=T, results='hold', eval=T------------------------------------------


# Reset simulation
SIM$reset()



# Set the region size in cM (0-4000cm, for testing the correctness of the function)
SIM$cm = seq(0,4000,l=length(SIM$cm))



A = SIM$addUnrelatedIndividual()
B = SIM$addUnrelatedIndividual()
C = SIM$mate(A,B)
D = SIM$mate(A,B)
G = SIM$addUnrelatedIndividual()
E = SIM$mate(G,C)



computePairIBD12(C,D)
computePairIBD12(E,A)


n = SIM$individuals_generated

IBD1matrix = 
sapply(1:n, function(y) {
        z = sapply(1:n, function(x) computePairIBD12(x,y) [1]) 
        names(z) = 1:n
        z
})

IBD2matrix = 
    sapply(1:n, function(y) {
        z = sapply(1:n, function(x) computePairIBD12(x,y) [2]) 
        names(z) = 1:n
        z
    })



printMatrix(IBD1matrix)


