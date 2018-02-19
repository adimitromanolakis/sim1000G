

library(plyr)

library(reshape)


if(1) {

    rm(list=ls())

    df = list()
    setwd("~/fs/tmp/sim4-pca//")
    x  = list.files(pattern="*.rdata")
    x = x[ x!="all.rdata" ]

    l = lapply(x,function(fname) { cat(fname,"\n"); load(fname); pz; })

    a = do.call(rbind,l)
    str(a)

}



cast(a, eff + pop1 ~ strat)



length(unique(a$gene))


cast(a, eff+pop1+strat ~ .)



table(a$pop1,a$strat)


d1=cast(a, eff+pop1+strat ~ ., function(x) mean(x < 0.05,na.rm=T) , value="pv1")
d2=cast(a, eff+pop1+strat ~ ., function(x) mean(x < 0.05,na.rm=T) , value="pv2")

merge(d1,d2,by=1:3)

dim(a)
length(unique(a$pv1))
length(unique(a$pv2))




d1=cast(a, eff+pop1+strat ~ ., function(x) mean(x < 0.001,na.rm=T) , value="pv1")
d2=cast(a, eff+pop1+strat ~ ., function(x) mean(x < 0.001,na.rm=T) , value="pv2")

merge(d1,d2,by=1:3)




d1=cast(a, eff+pop1+strat ~ ., function(x) mean(x < 0.001,na.rm=T) , value="pv1")
d2=cast(a, eff+pop1+strat ~ ., function(x) mean(x < 0.001,na.rm=T) , value="pv2")

d1=cast(a[a$eff==0,], eff+pop1+strat ~ ., function(x) quantile(x,p=0.05)
         , value="pv1")
d1




d1=cast(a[a$eff==0,], eff+pop1+strat ~ ., function(x) quantile(x,p=0.05)
        , value="pv1")
d1


efs = sort(unique(a$eff))
efs

l = list()
cnt = 1

for(i in efs)

for(strat in 2)
    for(pop1 in unique(a$pop1))

    {
     s = a$pop1 == pop1 & a$strat == strat & a$eff == i
     snull = a$pop1 == pop1 & a$strat == strat & a$eff == 0

     snull = a[snull,]
     s = a[s,]


     q = quantile(snull$pv2,p=0.01,na.rm=T)
     if(is.na(q)) { power2=NA; } else { power2 = mean(s$pv2<q) }

     q = quantile(snull$pv1,p=0.01,na.rm=T)
     if(is.na(q)) { power1=NA; } else { power1 = mean(s$pv1<q) }


     d = data.frame(exp(i),pop1,strat, q, power1, power2 )


     l[[cnt]] = d
     cnt = cnt+1


    }



l = ldply(l)
l

write.csv(l,"power.csv")

#### Create PDF of qqplots ####




ef = sort(unique(a$eff))

table(a$pop1,exp(a$eff))
exp(ef)

a$pv = a$pv1

pv_null = a [ a$eff == 0 & a$pop1 == 1600 & a$strat == 2,]
pv_null$type = 0


pv_causal = a [ a$eff == ef[3] & a$pop1 == 1600 & a$strat == 2,]
pv_causal$type = 1


mean(pv_null$pv1 < 0.05)
mean(pv_null$pv1 < 0.0001)

mean(pv_causal$pv1 < 0.05)


mean(pv_causal$pv1 < 0.01)
mean(pv_causal$pv1 < 0.001)
mean(pv_causal$pv1 < 0.0001)


set.seed(999)

pdf(w=9,h=5, file="~/fs/tmp/sim4/qq1.pdf");

plot(1)
title(paste("200 var 2 causal pop1=200 strat=1 ef=",exp(ef[4])))

for(i in 1:100) {

par(mfrow=c(1,2))
par(mar=c(5,4,2,1))


s1 =  sample( 1:nrow(pv_null)  ,200-3)
s2 =  sample( 1:nrow(pv_causal)  ,3)

all = rbind( pv_null[s1,]  , pv_causal[s2,] )


mean(all$pv < 0.05/100)
mean(all$pv2 < 0.05/100)



pv = all$pv
pv2 = all$pv2

ylim = max( -log10(c(pv,pv2)) )
ylim



r =  1:length(pv)/length(pv)
#r = ppoints(length(pv))

rsort = sort(-log10(r))

xlim = max(rsort)

all = all[ order(-log10(all$pv)),]

plot(rsort,-log10(all$pv1)   ,xlim=c(0,xlim),ylim=c(0,ylim)    ,
     xlab="Expected",ylab="Observed",cex=0.6     )  ;  abline(0,1);
s = which(all$type != 0)
points(rsort[s], -log10(all$pv1)[s],col= "blue",pch=20,cex=1.5)

title("No covariates")


all = all[ order(-log10(all$pv2)),]

plot(rsort,-log10(all$pv2)   ,xlim=c(0,xlim),ylim=c(0,ylim)  ,
     xlab="Expected",ylab="Observed"  ,cex=0.6               )  ;  abline(0,1);

title("With covariates")

s = which(all$type != 0)

points(rsort[s], -log10(all$pv2)[s],col="blue",pch=20,cex=1.5)



#qqplot( -log10(r), -log10(pv) ,xlim=c(0,3),ylim=c(0,ylim) )
#abline(0,1);



#qqplot( -log10(r), -log10(pv2),xlim=c(0,3) , ylim=c(0,ylim) )
#abline(0,1);
}


dev.off()
