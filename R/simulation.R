






# pkg.env <- new.env()

###set = function(x) { pkg.env$x = x }

###get = function() { pkg.env$x }




#' Holds data necessary for a simulation.
#' @export
SIM = new.env()



# vcf_url = "https://adimitromanolakis.github.io/sim1000G/inst/examples/region.vcf.gz"
# vcf_file = gzcon(url(vcf_url))
# vcf_file = textConnection( readLines(vcf_file) )
# vcf = readVCF(vcf_file, maxNumberOfVariants = 300 , min_maf = 0.12 ,max_maf = NA)
#
# downloadGeneticMap( chromosome  = 4)
# readGeneticMap( chromosome = 4)



#' Start the simulation
#'
#' @param vcf Input vcf file of a region (can be .gz). Must contain phased data.
#' @param totalNumberOfIndividuals Maximum Number of individuals that will ever be generated
#' @param randomdata If 1, disregards the genotypes in the vcf file and generates markers that are not in LD. Generally do not use.
#'
#' @examples
#' library("sim1000G")
#' library(gplots)
#'
#' examples_dir = system.file("examples", package = "sim1000G")
#' vcf_file = sprintf("%s/region.vcf.gz", examples_dir)
#' vcf = readVCF( vcf_file, maxNumberOfVariants = 100 , min_maf = 0.12 ,max_maf = NA)
#'
#' # For a realistic genetic map, use the functions downloadGeneticMap / readGeneticMap
#' generateFakeGeneticMap()
#'
#' plotRegionalGeneticMap(vcf$vcf[,2]+1)
#'
#' startSimulation(vcf, totalNumberOfIndividuals = 200)
#'
#' @export
startSimulation = function(vcf, totalNumberOfIndividuals = 250, randomdata = 0) {

        cat("[#####...] Creating SIM object\n");

        # Original haplotypes from 1000 genomes / other data

        SIM$population_gt1 = vcf$gt1
        SIM$population_gt2 = vcf$gt2



        if(randomdata) {

            str(SIM$population_gt1)
            # indoviduals in columns

            x = SIM$population_gt1*0

            for(i in 1:nrow(x)) x[i,] = rbinom( ncol(x), 1 , 0.25)

            SIM$population_gt1 = x
            SIM$population_gt2 = x

            # SIM$haplodata = haplodata(  t( cbind(SIM$population_gt1, SIM$population_gt2) ) )
            #    ldplot(SIM$haplodata$cor,ld.type="r")

        }




        # Generate haplodata object

        dim( t( cbind(SIM$population_gt1, SIM$population_gt2) )  )


        SIM$haplodata = haplodata(  t( cbind(SIM$population_gt1, SIM$population_gt2) ) )
        cat("[#####...] Haplodata object created\n");



        # Variant information

        SIM$varinfo = vcf$vcf[,1:8]
        SIM$bp = vcf$vcf[,2]
        SIM$cm = approx( geneticMap$bp, geneticMap$cm, SIM$bp )$y

        SIM$N_markers = nrow(vcf$gt1)

        SIM$pool = NA
        SIM$npool = 0

        SIM$total_individuals = totalNumberOfIndividuals
        SIM$individuals_generated = 0

        SIM$gt1 = matrix(nrow = SIM$total_individuals, ncol = SIM$N_markers, -1)
        SIM$gt2 = matrix(nrow = SIM$total_individuals, ncol = SIM$N_markers, -1)


        SIM$origin1 = matrix(nrow = SIM$total_individuals, ncol = SIM$N_markers, -1)
        SIM$origin2 = matrix(nrow = SIM$total_individuals, ncol = SIM$N_markers, -1)


        #SIM$cm = seq(0, 3200, l=dim(SIM$gt1)[2] )
        SIM$npool = 0
        SIM$last_ancestral_index = 10


}



SIM$generateNewHaplotypes = function(n = -1) {

    if(SIM$npool < 2) {
        cat("Generate new individual pool n=200\n");


        SIM$npool = 500
        SIM$pool = haplosim2(SIM$npool, SIM$haplodata, summary = F)$data


    }

    GT = list(gt1 =  SIM$pool[SIM$npool,], gt2 = SIM$pool[SIM$npool-1,] )
    SIM$npool = SIM$npool - 2

    SIM$last_ancestral_index = SIM$last_ancestral_index  + 1

    return(GT)

}


SIM$addUnrelatedIndividual = function() {

    newGenotypes = SIM$generateNewHaplotypes()

    if(SIM$individuals_generated >= SIM$total_individuals) {

        stop("No more space for saving new individual genotypes")

    }


    SIM$individuals_generated = SIM$individuals_generated + 1
    j = SIM$individuals_generated

    cat("Adding individual ",j, " from pool\n");


    SIM$gt1[j,] = newGenotypes$gt1
    SIM$gt2[j,] = newGenotypes$gt2

    SIM$origin1[j,] = SIM$last_ancestral_index
    SIM$origin2[j,] = -SIM$last_ancestral_index


    return(j)
}


SIM$addIndividualFromGenotypes = function(gt1,gt2) {

    if(SIM$individuals_generated >= SIM$total_individuals) {

        stop("No more space for generating new individual genotypes")

    }

    SIM$individuals_generated = SIM$individuals_generated + 1
    j = SIM$individuals_generated

    # cat("Adding individual ",j, " from specified genotypes\n");

    SIM$gt1[j,] = gt1
    SIM$gt2[j,] = gt2


    return(j)
}




debug_flag = 0


SIM$mate = function(i, j) {

    # For now, we hope i,j are of different sex

    recomb1 = generateSingleRecombinationVector(SIM$cm)
    recomb2 = generateSingleRecombinationVector(SIM$cm)


    # FATHER1 = SIM$gt1[i,]
    # FATHER2 = SIM$gt2[i,]
    # MOTHER1 = SIM$gt1[i,]
    # MOTHER2 = SIM$gt2[i,]



    gt1 =  recomb1 * SIM$gt1[i,]  +  (1-recomb1) * SIM$gt2[i,]
    gt2 =  recomb2 * SIM$gt1[j,]  +  (1-recomb2) * SIM$gt2[j,]



    #gt1 = SIM$gt1[i,]
    #gt1[recomb1==1] = SIM$gt2[i,][recomb1==1]

    #gt2 = SIM$gt1[j,]
    #gt2[recomb1==1] = SIM$gt2[j,][recomb1==1]

    #cat("Mate: ",i," " ,j,"\n");

    #p1 = paste(SIM$gt1[i,],collapse="")
    #p2 = paste(SIM$gt2[i,],collapse="")
    #p3 = paste(gt1,collapse="")

    #cat("R1:", paste(recomb1,collapse=""),"\n",sep="")
    #cat("R2:", paste(recomb2,collapse=""),"\n",sep="")

    #cat(p1,"\n",p2,"\n",p3,"\n",sep="")

    SWAP = runif(1) < 0.5
    if(SWAP) { t = gt1; gt1 = gt2; gt2 = t; }


    index = SIM$addIndividualFromGenotypes(gt1, gt2)


    #gt1 =  recomb1 * SIM$gt1[i,]  +  (1-recomb1) * SIM$gt2[i,]
    #gt2 =  recomb2 * SIM$gt1[j,]  +  (1-recomb2) * SIM$gt2[j,]

    if(!SWAP) {
        SIM$origin1[index,] = recomb1 * SIM$origin1[i,] + (1-recomb1) * SIM$origin2[i,];
        SIM$origin2[index,] = recomb2 * SIM$origin1[j,] + (1-recomb2) * SIM$origin2[j,];
    } else {
        SIM$origin2[index,] = recomb1 * SIM$origin1[i,] + (1-recomb1) * SIM$origin2[i,];
        SIM$origin1[index,] = recomb2 * SIM$origin1[j,] + (1-recomb2) * SIM$origin2[j,];
    }

    # last_gt <<- gt1 + gt2

    return (index)

}



SIM$reset = function() {
    SIM$individuals_generated = 0
}






















#' Simulates genotypes for 1 family with 1 offspring
#'
#'
#'
#' @param family_id What will be the family_id (for example: 100)
#'
#' @return family structure object
#'
#' @examples
#'
#' library("sim1000G")
#'
#' examples_dir = system.file("examples", package = "sim1000G")
#' vcf_file = sprintf("%s/region.vcf.gz", examples_dir)
#' vcf = readVCF( vcf_file, maxNumberOfVariants = 100 , min_maf = 0.12 ,max_maf = NA)
#'
#' # For realistic data use the functions downloadGeneticMap / readGeneticMap
#' generateFakeGeneticMap()
#'
#' startSimulation(vcf, totalNumberOfIndividuals = 1200)
#' fam1 = newNuclearFamily(1)
#' fam2 = newNuclearFamily(2)
#'
#' # See also the documentation on our github page
#'
#' @export
newNuclearFamily = function(family_id) {

    fam = data.frame(family_id = family_id  ,
                     id = c(1,2,3) ,
                     father = c(0,0,1) ,
                     mother = c(0,0,2) ,
                     sex = c(1,2,1),
                     generation = c(1,1,2)
    )


    j1 = SIM$addUnrelatedIndividual()
    j2 = SIM$addUnrelatedIndividual()
    j3 = SIM$mate(j1,j2)


    fam$gtindex = c(j1,j2,j3)

    fam
}

#' Simulates genotypes for 1 family with n offspring
#'
#'
#'
#' @param family_id What will be the family_id (for example: 100)
#' @param noffspring Number of offsprings that this family will have
#'
#' @return family structure object
#'
#' @examples
#'
#' ped_line = newFamilyWithOffspring(10,3)
#'
#'
#' @export
newFamilyWithOffspring = function(family_id, noffspring = 2) {

    fam = data.frame(fid = family_id  ,
                     id = c(1:2) ,
                     father = c(0,0),
                     mother = c(0,0),
                     sex = c(1,2),
                     generation = c(1,1)
    )


    j1 = SIM$addUnrelatedIndividual()
    j2 = SIM$addUnrelatedIndividual()

    fam$gtindex = c(j1,j2)

    for(i in 1:noffspring) {
        j3 = SIM$mate(j1,j2)
        newFamLine = c(family_id, i+10, 1,2, sample(c(1,2), 1) ,2 , j3)
        fam = rbind(fam, newFamLine)
    }

    return (fam)
}









#' Generates genotype data for a family of 3 generations
#'
#'
#'
#' @param familyid What will be the family_id (for example: 100)
#' @param noffspring2 Number of offspring in generation 2
#' @param noffspring3 Number of offspring in generation 3 (vector of length noffspring2)
#'
#' @return family structure object
#'
#' @examples
#'
#' library("sim1000G")
#'
#' examples_dir = system.file("examples", package = "sim1000G")
#' vcf_file = sprintf("%s/region.vcf.gz", examples_dir)
#' vcf = readVCF( vcf_file, maxNumberOfVariants = 100 , min_maf = 0.12 ,max_maf = NA)
#'
#' # For realistic data use the functions downloadGeneticMap / readGeneticMap
#' generateFakeGeneticMap()
#'
#' startSimulation(vcf, totalNumberOfIndividuals = 200)
#'
#' ped_line = newFamily3generations(12, 3, c(3,3,2) )
#'
#' @export
newFamily3generations = function(familyid, noffspring2 = 2, noffspring3 = c(1,1)) {

    if(length(noffspring3) == noffspring2){

        fam = data.frame(fid = familyid  ,
                         id = c(1:2) ,
                         father = c(0,0),
                         mother = c(0,0),
                         sex = c(1,2),
                         generation = c(1,1)
        )


        j1 = SIM$addUnrelatedIndividual()
        j2 = SIM$addUnrelatedIndividual()

        fam$gtindex = c(j1,j2)

        id2 <- 2
        for(i in 1:noffspring2) {
            id2 <- id2 + 1
            ids <- id2 + 1
            j3 = SIM$mate(j1,j2)
            gender2 <- sample(c(1,2), 1)
            genders <- ifelse(gender2 == 1, 2, 1)
            newFamLine = c(familyid, id2, 1,2, gender2,2,  j3)
            fam = rbind(fam, newFamLine)

            if(length(noffspring3[i] > 0)){
                j4 <- SIM$addUnrelatedIndividual()
                newFamLine = c(familyid, ids, 0,0, genders ,2,  j4)
                fam = rbind(fam, newFamLine)

                id3 <- ids
                for(j in 1:noffspring3[i]){
                    id3 <- id3 + 1
                    j5 <- SIM$mate(j3, j4)
                    father_id <- ifelse(gender2 == 1, id2, ids)
                    mother_id <- ifelse(gender2 == 2, id2, ids)
                    newFamLine = c(familyid, id3, father_id, mother_id, sample(c(1, 2), 1) ,3,  j5)
                    fam = rbind(fam, newFamLine)


                }
                id2 <- id3
            }
        }

        return (fam)
    }  else  {
        print(noffspring2)
        print(noffspring3)
        stop("The vector for number of offspring in the third generation (possibly zeros) must be equal to the number of offspring in the second generation")
    }
}






#' Computes pairwise IBD1/2 for a specific pair of individuals
#'
#'
#'
#' @param i Index of first individual
#' @param j Index of second individual
#'
#' @return Mean IBD1 and IBD2 as computed from shared haplotypes
#'
#' @examples
#'
#' library("sim1000G")
#'
#' examples_dir = system.file("examples", package = "sim1000G")
#' vcf_file = sprintf("%s/region.vcf.gz", examples_dir)
#' vcf = readVCF( vcf_file, maxNumberOfVariants = 100 , min_maf = 0.12 ,max_maf = NA)
#'
#' # For realistic data use the functions downloadGeneticMap / readGeneticMap
#' generateFakeGeneticMap()
#'
#' startSimulation(vcf, totalNumberOfIndividuals = 200)
#'
#' ped1 = newNuclearFamily(1)
#'
#' v = computePairIBD12(1, 3)
#'
#' cat("IBD1 of pair = ", v[1], "\n");
#' cat("IBD2 of pair = ", v[2], "\n");
#'
#'
#' @export
computePairIBD12 = function(i,j) {

    q1 = SIM$origin1[i,]
    q2 = SIM$origin2[i,]

    r1 = SIM$origin1[j,]
    r2 = SIM$origin2[j,]


    table(q1,q2)

    IBD1H1 = q1 == r1 | q1 == r2
    IBD1H2 = q2 == r1 | q2 == r2


    IBD12 = mean(IBD1H1 | IBD1H2)
    IBD2 = mean(   (q1 == r1 & q2 == r2 )  | (q1 == r2 & q2 == r1) )


    c(IBD1=IBD12-IBD2, IBD2=IBD2)
}






#' Computes pairwise IBD1 for a specific pair of individuals.
#' See function computePairIBD12 for description.
#'
#'
#'
#' @param i Index of first individual
#' @param j Index of second individual
#'
#' @return Mean IBD1 as computed from shared haplotypes
#'
#' @examples
#'
#' library("sim1000G")
#'
#' examples_dir = system.file("examples", package = "sim1000G")
#' vcf_file = sprintf("%s/region.vcf.gz", examples_dir)
#' vcf = readVCF( vcf_file, maxNumberOfVariants = 100 , min_maf = 0.12 ,max_maf = NA)
#'
#' # For realistic data use the functions downloadGeneticMap / readGeneticMap
#' generateFakeGeneticMap()
#'
#' startSimulation(vcf, totalNumberOfIndividuals = 200)
#'
#' ped1 = newNuclearFamily(1)
#'
#' v = computePairIBD1(1, 3)
#'
#' cat("IBD1 of pair = ", v, "\n");
#'
#' @export
computePairIBD1 = function(i,j) {

    q1 = SIM$origin1[i,]
    q2 = SIM$origin2[i,]

    r1 = SIM$origin1[j,]
    r2 = SIM$origin2[j,]



    IBD1H1 = q1 == r1 | q1 == r2
    IBD1H2 = q2 == r1 | q2 == r2

    IBD12 = mean(IBD1H1 | IBD1H2)
    IBD2 = mean(   (q1 == r1 & q2 == r2 )  | (q1 == r2 & q2 == r1) )


    IBD1=IBD12-IBD2

}



#' Computes pairwise IBD2 for a specific pair of individuals
#'
#'
#'
#' @param i Index of first individual
#' @param j Index of second individual
#'
#' @return Mean IBD2 as computed from shared haplotypes
#'
#'
#' @examples
#'
#' library("sim1000G")
#'
#' examples_dir = system.file("examples", package = "sim1000G")
#' vcf_file = sprintf("%s/region.vcf.gz", examples_dir)
#' vcf = readVCF( vcf_file, maxNumberOfVariants = 100 , min_maf = 0.12 ,max_maf = NA)
#'
#' # For realistic data use the functions downloadGeneticMap / readGeneticMap
#' generateFakeGeneticMap()
#'
#' startSimulation(vcf, totalNumberOfIndividuals = 200)
#'
#' ped1 = newNuclearFamily(1)
#'
#' v = computePairIBD2(1, 3)
#'
#' cat("IBD2 of pair = ", v, "\n");
#'
#' @export
computePairIBD2 = function(i,j) {

  q1 = SIM$origin1[i,]
  q2 = SIM$origin2[i,]

  r1 = SIM$origin1[j,]
  r2 = SIM$origin2[j,]

   IBD2 = mean(   (q1 == r1 & q2 == r2 )  | (q1 == r2 & q2 == r1) )


  IBD2

}















#' Utility function that prints a matrix. Useful for IBD12 matrices.
#'
#'
#' @param m Matrix to be printed
#'
#' @examples
#'
#' printMatrix (  matrix(runif(16), nrow=4) )
#' @export
#'
printMatrix = function(m) {
    cat("      " , " " , sprintf("[%4d] ", 1:nrow(m) )  , "\n")

    for(i in 1:nrow(m)) {

        cat(sprintf("[%4d]",i) , " " , sprintf(" %.3f ", m[i,]) , "\n")
    }

}


















