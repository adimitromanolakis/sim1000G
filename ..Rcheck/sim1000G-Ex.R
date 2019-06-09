pkgname <- "sim1000G"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('sim1000G')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("computePairIBD1")
### * computePairIBD1

flush(stderr()); flush(stdout())

### Name: computePairIBD1
### Title: Computes pairwise IBD1 for a specific pair of individuals. See
###   function computePairIBD12 for description.
### Aliases: computePairIBD1

### ** Examples


library("sim1000G")

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")
vcf = readVCF( vcf_file, maxNumberOfVariants = 100 ,
               min_maf = 0.12 ,max_maf = NA)

# For realistic data use the function downloadGeneticMap
generateUniformGeneticMap()

startSimulation(vcf, totalNumberOfIndividuals = 200)

ped1 = newNuclearFamily(1)

v = computePairIBD1(1, 3)

cat("IBD1 of pair = ", v, "\n");




cleanEx()
nameEx("computePairIBD12")
### * computePairIBD12

flush(stderr()); flush(stdout())

### Name: computePairIBD12
### Title: Computes pairwise IBD1/2 for a specific pair of individuals
### Aliases: computePairIBD12

### ** Examples


library("sim1000G")

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")

vcf = readVCF( vcf_file, maxNumberOfVariants = 100 ,
               min_maf = 0.12 ,max_maf = NA)

generateUniformGeneticMap()

startSimulation(vcf, totalNumberOfIndividuals = 200)

ped1 = newNuclearFamily(1)

v = computePairIBD12(1, 3)

cat("IBD1 of pair = ", v[1], "\n");
cat("IBD2 of pair = ", v[2], "\n");





cleanEx()
nameEx("computePairIBD2")
### * computePairIBD2

flush(stderr()); flush(stdout())

### Name: computePairIBD2
### Title: Computes pairwise IBD2 for a specific pair of individuals
### Aliases: computePairIBD2

### ** Examples


library("sim1000G")

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")
vcf = readVCF( vcf_file, maxNumberOfVariants = 100 ,
               min_maf = 0.12 ,max_maf = NA)

# For realistic data use the function downloadGeneticMap
generateUniformGeneticMap()

startSimulation(vcf, totalNumberOfIndividuals = 200)

ped1 = newNuclearFamily(1)

v = computePairIBD2(1, 3)

cat("IBD2 of pair = ", v, "\n");




cleanEx()
nameEx("downloadGeneticMap")
### * downloadGeneticMap

flush(stderr()); flush(stdout())

### Name: downloadGeneticMap
### Title: Downloads a genetic map for a particular chromosome under GRCh37
###   coordinates for use with sim1000G.
### Aliases: downloadGeneticMap

### ** Examples




downloadGeneticMap(22, dir=tempdir() )





cleanEx()
nameEx("generateChromosomeRecombinationPositions")
### * generateChromosomeRecombinationPositions

flush(stderr()); flush(stdout())

### Name: generateChromosomeRecombinationPositions
### Title: Generates a recombination vector arising from one meiotic event.
###   The origin of segments is coded as (0 - haplotype1 , 1 - haplotype2 )
### Aliases: generateChromosomeRecombinationPositions

### ** Examples


library("sim1000G")

# generate a recombination events for chromosome 4
readGeneticMap(4)
generateChromosomeRecombinationPositions(500)




cleanEx()
nameEx("generateFakeWholeGenomeGeneticMap")
### * generateFakeWholeGenomeGeneticMap

flush(stderr()); flush(stdout())

### Name: generateFakeWholeGenomeGeneticMap
### Title: Generates a fake genetic map that spans the whole genome.
### Aliases: generateFakeWholeGenomeGeneticMap

### ** Examples


library("sim1000G")

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = sprintf("%s/region.vcf.gz", examples_dir)
vcf = readVCF( vcf_file, maxNumberOfVariants = 100 ,
               min_maf = 0.12 ,max_maf = NA)

# For realistic data use the function
# downloadGeneticMap
generateFakeWholeGenomeGeneticMap(vcf)

pdf(file=tempfile())
plotRegionalGeneticMap(seq(1e6,100e6,by=1e6/2))
dev.off()




cleanEx()
nameEx("generateRecombinationDistances")
### * generateRecombinationDistances

flush(stderr()); flush(stdout())

### Name: generateRecombinationDistances
### Title: Generate inter-recombination distances using a chi-square model.
###   Note this are the distances between two succesive recombination
###   events and not the absolute positions of the events. To generate the
###   locations of the recombination events see the example below.
### Aliases: generateRecombinationDistances

### ** Examples


library("sim1000G")

distances = generateRecombinationDistances(20)


positions_of_recombination = cumsum(distances)

if(0) hist(generateRecombinationDistances(20000),n=100)




cleanEx()
nameEx("generateRecombinationDistances_noInterference")
### * generateRecombinationDistances_noInterference

flush(stderr()); flush(stdout())

### Name: generateRecombinationDistances_noInterference
### Title: Generate recombination distances using a no-interference model.
### Aliases: generateRecombinationDistances_noInterference

### ** Examples


library("sim1000G")
mean ( generateRecombinationDistances_noInterference ( 200 ) )




cleanEx()
nameEx("generateSingleRecombinationVector")
### * generateSingleRecombinationVector

flush(stderr()); flush(stdout())

### Name: generateSingleRecombinationVector
### Title: Genetates a recombination vector arising from one meiotic event.
###   The origin of segments is coded as (0 - haplotype1 , 1 - haplotype2 )
### Aliases: generateSingleRecombinationVector

### ** Examples


library("sim1000G")

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")
vcf = readVCF( vcf_file, maxNumberOfVariants = 100 ,
               min_maf = 0.12 ,max_maf = NA)

# For realistic data use the function downloadGeneticMap
generateUniformGeneticMap()
generateSingleRecombinationVector( 1:100 )




cleanEx()
nameEx("generateUniformGeneticMap")
### * generateUniformGeneticMap

flush(stderr()); flush(stdout())

### Name: generateUniformGeneticMap
### Title: Generates a uniform genetic map.
### Aliases: generateUniformGeneticMap

### ** Examples


library("sim1000G")

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = sprintf("%s/region.vcf.gz", examples_dir)
vcf = readVCF( vcf_file, maxNumberOfVariants = 100 ,
               min_maf = 0.12 ,max_maf = NA)

# For realistic data use the function readGeneticMap
generateUniformGeneticMap()

pdf(file=tempfile())
plotRegionalGeneticMap(seq(1e6,100e6,by=1e6/2))
dev.off()




cleanEx()
nameEx("generateUnrelatedIndividuals")
### * generateUnrelatedIndividuals

flush(stderr()); flush(stdout())

### Name: generateUnrelatedIndividuals
### Title: Generates variant data for n unrelated individuals
### Aliases: generateUnrelatedIndividuals

### ** Examples


library("sim1000G")

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")
vcf = readVCF( vcf_file, maxNumberOfVariants = 100 , min_maf = 0.12)

genetic_map_of_region =
   system.file("examples",
     "chr4-geneticmap.txt",
     package = "sim1000G")

readGeneticMapFromFile(genetic_map_of_region)

startSimulation(vcf, totalNumberOfIndividuals = 1200)
ids = generateUnrelatedIndividuals(20)

# See also the documentation on our github page




cleanEx()
nameEx("getCMfromBP")
### * getCMfromBP

flush(stderr()); flush(stdout())

### Name: getCMfromBP
### Title: Converts centimorgan position to base-pair. Return a list of
###   centimorgan positions that correspond to the bp vector (in
###   basepairs).
### Aliases: getCMfromBP

### ** Examples


library("sim1000G")

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = sprintf("%s/region.vcf.gz", examples_dir)
vcf = readVCF( vcf_file, maxNumberOfVariants = 100,
  min_maf = 0.12)

# For realistic data use the function downloadGeneticMap
generateUniformGeneticMap()
getCMfromBP(seq(1e6,100e6,by=1e6))





cleanEx()
nameEx("loadSimulation")
### * loadSimulation

flush(stderr()); flush(stdout())

### Name: loadSimulation
### Title: Load some previously saved simulation data by function
###   saveSimulation
### Aliases: loadSimulation

### ** Examples



examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")

vcf = readVCF( vcf_file, maxNumberOfVariants = 100 ,
           min_maf = 0.12 ,max_maf = NA)

# For a realistic genetic map
# use the function readGeneticMap
generateUniformGeneticMap()

startSimulation(vcf, totalNumberOfIndividuals = 200)

ped1 = newNuclearFamily(1)

saveSimulation("sim1")

vcf = readVCF( vcf_file, maxNumberOfVariants = 100 ,
               min_maf = 0.02 ,max_maf = 0.5)

startSimulation(vcf, totalNumberOfIndividuals = 200)
saveSimulation("sim2")

loadSimulation("sim1")






cleanEx()
nameEx("newFamily3generations")
### * newFamily3generations

flush(stderr()); flush(stdout())

### Name: newFamily3generations
### Title: Generates genotype data for a family of 3 generations
### Aliases: newFamily3generations

### ** Examples


library("sim1000G")

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")
vcf = readVCF( vcf_file, maxNumberOfVariants = 100 ,
               min_maf = 0.12 ,max_maf = NA)

generateUniformGeneticMap()

startSimulation(vcf, totalNumberOfIndividuals = 200)

ped_line = newFamily3generations(12, 3, c(3,3,2) )




cleanEx()
nameEx("newFamilyWithOffspring")
### * newFamilyWithOffspring

flush(stderr()); flush(stdout())

### Name: newFamilyWithOffspring
### Title: Simulates genotypes for 1 family with n offspring
### Aliases: newFamilyWithOffspring

### ** Examples


ped_line = newFamilyWithOffspring(10,3)





cleanEx()
nameEx("newNuclearFamily")
### * newNuclearFamily

flush(stderr()); flush(stdout())

### Name: newNuclearFamily
### Title: Simulates genotypes for 1 family with 1 offspring
### Aliases: newNuclearFamily

### ** Examples


library("sim1000G")

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")
vcf = readVCF( vcf_file, maxNumberOfVariants = 100 ,
   min_maf = 0.12 ,max_maf = NA)

genetic_map_of_region = system.file("examples","chr4-geneticmap.txt",
   package = "sim1000G")
readGeneticMapFromFile(genetic_map_of_region)

startSimulation(vcf, totalNumberOfIndividuals = 1200)
fam1 = newNuclearFamily(1)
fam2 = newNuclearFamily(2)

# See also the documentation on our github page




cleanEx()
nameEx("printMatrix")
### * printMatrix

flush(stderr()); flush(stdout())

### Name: printMatrix
### Title: Utility function that prints a matrix. Useful for IBD12
###   matrices.
### Aliases: printMatrix

### ** Examples


printMatrix (  matrix(runif(16), nrow=4) )



cleanEx()
nameEx("readGeneticMap")
### * readGeneticMap

flush(stderr()); flush(stdout())

### Name: readGeneticMap
### Title: Reads a genetic map downloaded from the function
###   downloadGeneticMap or reads a genetic map from a specified file. If
###   the argument filename is used then the genetic map is read from the
###   corresponding file. Otherwise, if a chromosome is specified, the
###   genetic map is downloaded for human chromosome using grch37
###   coordinates.
### Aliases: readGeneticMap

### ** Examples





readGeneticMap(chromosome = 22)






cleanEx()
nameEx("readGeneticMapFromFile")
### * readGeneticMapFromFile

flush(stderr()); flush(stdout())

### Name: readGeneticMapFromFile
### Title: Reads a genetic map to be used for simulations. The genetic map
###   should be of a single chromosome and covering the extent of the
###   region to be simulated. Whole chromosome genetic maps can also be
###   used.
### Aliases: readGeneticMapFromFile

### ** Examples


## Not run: 
##D 
##D fname = downloadGeneticMap(10)
##D 
##D cat("genetic map downloaded at :", fname, "\n")
##D readGeneticMapFromFile(fname)
##D 
## End(Not run)



cleanEx()
nameEx("readVCF")
### * readVCF

flush(stderr()); flush(stdout())

### Name: readVCF
### Title: Read a vcf file, with options to filter out low or high
###   frequency markers.
### Aliases: readVCF

### ** Examples


examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir,
  "region-chr4-93-TMEM156.vcf.gz")

vcf = readVCF( vcf_file, maxNumberOfVariants = 500 ,
               min_maf = 0.02 ,max_maf = NA)

str(as.list(vcf))



cleanEx()
nameEx("resetSimulation")
### * resetSimulation

flush(stderr()); flush(stdout())

### Name: resetSimulation
### Title: Removes all individuals that have been simulated and resets the
###   simulator.
### Aliases: resetSimulation

### ** Examples


resetSimulation()




cleanEx()
nameEx("retrieveGenotypes")
### * retrieveGenotypes

flush(stderr()); flush(stdout())

### Name: retrieveGenotypes
### Title: Retrieve a matrix of simulated genotypes for a specific set of
###   individual IDs
### Aliases: retrieveGenotypes

### ** Examples


library("sim1000G")

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")
vcf = readVCF( vcf_file, maxNumberOfVariants = 100 ,
               min_maf = 0.12 ,max_maf = NA)

# For realistic data use the function downloadGeneticMap
generateUniformGeneticMap()

startSimulation(vcf, totalNumberOfIndividuals = 200)

ped1 = newNuclearFamily(1)

retrieveGenotypes(ped1$gtindex)




cleanEx()
nameEx("saveSimulation")
### * saveSimulation

flush(stderr()); flush(stdout())

### Name: saveSimulation
### Title: Save the data for a simulation for later use. When simulating
###   multiple populations it allows saving and restoring of simulation
###   data for each population.
### Aliases: saveSimulation

### ** Examples




examples_dir = system.file("examples", package = "sim1000G")

vcf_file = file.path(examples_dir, "region.vcf.gz")
vcf = readVCF( vcf_file, maxNumberOfVariants = 100 ,
               min_maf = 0.12 ,max_maf = NA)


# For realistic data use the functions downloadGeneticMap
generateUniformGeneticMap()

startSimulation(vcf, totalNumberOfIndividuals = 200)

ped1 = newNuclearFamily(1)

saveSimulation("sim1")




cleanEx()
nameEx("setRecombinationModel")
### * setRecombinationModel

flush(stderr()); flush(stdout())

### Name: setRecombinationModel
### Title: Set recombination model to either poisson (no interference) or
###   chi-square.
### Aliases: setRecombinationModel

### ** Examples



generateUniformGeneticMap()

do_plots = 0

setRecombinationModel("chisq")
if(do_plots == 1)
 hist(generateRecombinationDistances(100000),n=200)

setRecombinationModel("poisson")
if(do_plots == 1)
 hist(generateRecombinationDistances(100000),n=200)




cleanEx()
nameEx("startSimulation")
### * startSimulation

flush(stderr()); flush(stdout())

### Name: startSimulation
### Title: Starts and initializes the data structures required for a
###   simulation. A VCF file should be read beforehand with the function
###   readVCF.
### Aliases: startSimulation

### ** Examples

library("sim1000G")
library(gplots)

examples_dir = system.file("examples", package = "sim1000G")
vcf_file = file.path(examples_dir, "region.vcf.gz")

vcf = readVCF( vcf_file, maxNumberOfVariants = 100)


genetic_map_of_region = system.file(
   "examples",
   "chr4-geneticmap.txt",
   package = "sim1000G"
)

readGeneticMapFromFile(genetic_map_of_region)

pdf(file=tempfile())
plotRegionalGeneticMap(vcf$vcf[,2]+1)
dev.off()

startSimulation(vcf, totalNumberOfIndividuals = 200)




cleanEx()
nameEx("subsetVCF")
### * subsetVCF

flush(stderr()); flush(stdout())

### Name: subsetVCF
### Title: Generate a market subset of a vcf file
### Aliases: subsetVCF

### ** Examples


examples_dir = system.file("examples", package = "sim1000G")

vcf_file = file.path(examples_dir, "region-chr4-93-TMEM156.vcf.gz")

vcf = readVCF( vcf_file, maxNumberOfVariants = 500 ,
               min_maf = 0.02 ,max_maf = NA)

vcf2 = subsetVCF(vcf, var_index = 1:50)




### * <FOOTER>
###
options(digits = 7L)
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
