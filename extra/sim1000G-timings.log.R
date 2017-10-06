s = 7

library(ggplot2)
a = read.table("1")

pdf(w=s,h=s/1.4);
 

ggplot(a,aes(V4,V5,pch=factor(V3) )) +geom_point(cex=4)+geom_line() + scale_y_log10() + 
theme_light(10) +labs(x="Number of individuals",y="Time (sec)") + scale_shape_discrete(name="Number of\nmarkers") + scale_x_continuous(limits=c(100,8000)); 


dev.off()

