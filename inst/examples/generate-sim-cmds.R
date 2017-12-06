



nsim=100


for(i in 1:20) {
for(pop1 in c(0,200))
for(eff in c(0,1,2,5,10)) {
 for(pop_strat in c(0,0.5,1,2)) {


         cmd = paste("~/tools/bin/Rscript example2-population-stratification.R ",
                     eff," ",
                     pop1," ",
                     pop_strat," ",
                     nsim,
                     sep="")

         cmd = sprintf("submit sim6-%.2f-%.1f-%d-rep%d %s",
                       eff,pop_strat,pop1,i,
                       cmd)



            cat(cmd,"\n")

     }

 }

}
