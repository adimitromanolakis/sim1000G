

default:
	Rscript  -e 'library(devtools); build(); document(); '

doc:
	Rscript  -e 'library(devtools); document(); '

tgz:
	cd .. && tar cvzf sim1000G.tar.gz sim1000G/ --exclude sim1000G/.git --exclude sim1000G/inst/doc

check:
	R CMD check .


pdf:
	R CMD Rd2pdf .

install:
	cd .. && R CMD INSTALL BF

package:
	cd .. && R CMD build BF

