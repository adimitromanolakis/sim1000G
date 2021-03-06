% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simulation.R
\name{saveSimulation}
\alias{saveSimulation}
\title{Save the data for a simulation for later use. When simulating multiple populations it
allows saving and restoring of simulation data for each population.}
\usage{
saveSimulation(id)
}
\arguments{
\item{id}{Name the simulation will be saved as.}
}
\description{
Save the data for a simulation for later use. When simulating multiple populations it
allows saving and restoring of simulation data for each population.
}
\examples{



examples_dir = system.file("examples", package = "sim1000G")

vcf_file = file.path(examples_dir, "region.vcf.gz")
vcf = readVCF( vcf_file, maxNumberOfVariants = 100 ,
               min_maf = 0.12 ,max_maf = NA)


# For realistic data use the functions downloadGeneticMap
generateUniformGeneticMap()

startSimulation(vcf, totalNumberOfIndividuals = 200)

ped1 = newNuclearFamily(1)

saveSimulation("sim1")

}
