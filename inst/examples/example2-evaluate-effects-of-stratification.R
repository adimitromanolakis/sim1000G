

a = read.table("~/fs/tmp/sim/1",as=T)


load("~/fs/tmp/sim/all.rdata")


a = lb
str(a)
library(plyr)

table(a$V5)
table(a$V6)

if(0) {
colnames(a)[3] = "eff"
colnames(a)[6] = "pop1"
colnames(a)[4] = "strat"
colnames(a)[7] = "pv"
colnames(a)[8] = "pv2"
}

library(reshape)
cast(a, eff + pop1 ~ strat)

table(a$pop1,a$strat)


cast(a, eff+pop1+strat ~ ., function(x) mean(x < 0.05/1,na.rm=T) , value="pv1")



a$pv = a$pv1

pv_null = a [ a$eff == 0 & a$pop1 == 200 & a$strat == 1,]


pv_causal = a [ a$eff == 10 & a$pop1 == 200 & a$strat == 1,]



par(mfrow=c(1,2))
par(mar=c(5,4,2,1))


s1 =  sample( 1:nrow(pv_null)  ,500)
s2 =  sample( 1:nrow(pv_causal)  ,5)

all = rbind( pv_null[s1,]  , pv_causal[s2,] )


mean(all$pv < 0.05/100)
mean(all$pv2 < 0.05/100)



pv = all$pv
r =  1:length(pv)/length(pv)
r = ppoints(length(pv))

qqplot( -log10(r), -log10(pv))
abline(0,1);



pv = all$pv2
qqplot( -log10(r), -log10(pv))
abline(0,1);

