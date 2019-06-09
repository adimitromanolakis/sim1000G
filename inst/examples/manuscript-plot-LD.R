

library(sim1000G)
librar(gplots)
simulatePopulation = function(num_individuals = 300)
  
{
  ids = generateUnrelatedIndividuals(num_individuals)
  genotypes = retrieveGenotypes(ids)
}







vcf_file = system.file("example","region.vcf.gz",package="sim1000G")

vcf_file = "~/MEGAsync/work/lunenfeld/sim1000G/inst/examples/region-chr4-93-TMEM156.vcf.gz"
#vcf_file = "~/MEGAsync/work/lunenfeld/sim1000G/inst/examples/region.vcf.gz"


vcf_file = "/tmp/fs/tmp/CEU-TSI-GBR-region-chr4-34-ABLIM2.vcf.gz"
x = list.files(path = "/tmp/fs/tmp/CEU-TSI-GBR/", full.names = T)



vcf_file = "/tmp/fs/tmp/CEU-TSI-GBR/CEU-TSI-GBR-region-chr4-202-NOCT.vcf.gz"
vcf_file = sample(x,1)


vcf_file = "/tmp/fs/tmp/extract-2mbp/eur-2MBp/region-chr1-12-r12.0MBp.vcf.gz"


#vcf_file = "/tmp/fs/tmp/CEU-TSI-GBR//CEU-TSI-GBR-region-chr4-386-GALNTL6.vcf.gz"
vcf = readVCF(vcf_file,
              maxNumberOfVariants = 1000 ,
              min_maf = 0.15 ,
              max_maf = 1-0.15)


str(as.list(vcf))

vcf2 = subsetVCF(vcf,var_id = vcf$varid[100+(1:200) ])
#vcf2= vcf

startSimulation(vcf2, totalNumberOfIndividuals = 3000)






examineMAFandLD = function() {
  
  SIM$reset()
  
  
  ids = generateUnrelatedIndividuals(1900)
  genotypes = retrieveGenotypes(ids)
  
  
  ld_population =   cor   (  t(SIM$population_gt1 + SIM$population_gt2 )   ) ^2
  
  ld_simulated_data = cor(genotypes)^2
  ld_simulated_data[1:10,1:10]
  
  s1 = ld_simulated_data[ lower.tri(ld_simulated_data) ]
  s2 = ld_population[ lower.tri(ld_simulated_data) ]
  
  plot(s1,s2)
  
  maf1 = apply(genotypes,2,mean)/2
  maf2 = apply( t(vcf2$gt1+vcf2$gt2)  , 2 ,mean)/2
  maf1[maf1>0.5] = 1-maf1[maf1>0.5]
  maf2[maf2>0.5] = 1-maf2[maf2>0.5]
  
  plot(maf1,maf2)
  
  s = which(maf2>0.15)
  
  
  n = nrow(ld_simulated_data)
  #ld_simulated_data [ lower.tri(ld_simulated_data)  ] = ld_population[ lower.tri(ld_population) ]
  
  
  library(gplots)
  
  
  pdf(file="/tmp/2.pdf",w=15,h=11)
  heatmap.2(ld_simulated_data[s,s] , col=rev(heat.colors(299)) , trace="none",Rowv=F,Colv=F,
         #   rowsep = 0:10e3, colsep =  0:10e3,           sepwidth=c(0.001,0.001),sepcolor=rgb(0,0,0,0.04),
            offsetRow = 100,offsetCol = 100, breaks=seq(0,1,l=300), 
            
            density.info = "none",keysize = 1,
            
            row.names = "")
  
  heatmap.2(ld_population[s,s] , col=rev(heat.colors(299)) , trace="none",Rowv=F,Colv=F,
           # rowsep = 0:10e3, colsep =  0:10e3,           sepwidth=c(0.001,0.001),sepcolor=rgb(0,0,0,0.04),
            offsetRow = 100,offsetCol = 100, breaks=seq(0,1,l=300),
            density.info = "none",keysize = 1,
            
            row.names = "")
  
  dev.off()
  
  
  
}



examineMAFandLD()



















library(plotly)
packageVersion('plotly')

ld_simulated_data[is.na(ld_simulated_data)] = 0

p <- plot_ly(z = ld_simulated_data, type = "heatmap",  colors = rev(heat.colors(250)) )
p
export(p, file = "/tmp/ldsim.png")

ld_population[is.na(ld_population)] = 0
p <- plot_ly(z = ld_population, type = "heatmap",  colors = rev(heat.colors(250)) )
p
export(p, file = "/tmp/ldpop.png")




plot_ly(z = ~ld_simulated_data^2, type = "surface",  colors = rev(heat.colors(10)))


vals <- unique(scales::rescale(c(ld_simulated_data)))
o <- order(vals, decreasing = FALSE)
cols <- scales::col_numeric("Blues", domain = NULL)(vals)
colz <- setNames(data.frame(vals[o], cols[o]), NULL)

str(colz)

colz[,2] = (redblue(nrow(colz)))
str(colz)
p <- plot_ly(z = ld_simulated_data, colorscale = colz, type = "heatmap")
p

ld_population[is.na(ld_population)] = 0

write.table(ld_population,file="/tmp/t3.csv",quote=F,row=F,col=F,sep=",")
