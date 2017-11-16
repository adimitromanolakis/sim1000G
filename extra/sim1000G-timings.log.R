s = 9

library(ggplot2)
a = read.table("1")
b = read.table("2")

a$V5 = a$V5+b$V5


pdf(w=s,h=s/1.9);


ggplot(a,aes(V4,V5,pch=factor(V3) )) +geom_point(cex=3)+geom_line()  + scale_y_log10() +
theme_light(15) +labs(x="Number of individuals",y="Time (sec)") + scale_shape_discrete(name="Number of\nvariants") + scale_x_continuous(limits=c(100,8000));

# grep Init sim1000G-timings.log.txt  >2

a = read.table("1")


ggplot(a,aes(V4,V5,pch=factor(V3) )) +geom_point(cex=3)+geom_line() + scale_y_log10() +
    theme_light(10) +labs(x="Number of individuals",y="Time (sec)") + scale_shape_discrete(name="Number of\nvariants") + scale_x_continuous(limits=c(100,8000));


dev.off()
