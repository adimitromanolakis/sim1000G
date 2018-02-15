

default:
	Rscript  -e 'library(devtools); document(); build(); '

doc:
	Rscript -e 'tools::buildVignettes(dir = ".") ' 
	Rscript  -e 'library(devtools); document(); '

vin:
	Rscript -e 'tools::buildVignettes(package="sim1000G",quiet=F)'


crancheck:
	cd .. && R CMD check --as-cran $(PKGFILE)


tgz:
	cd .. && tar cvzf sim1000G.tar.gz `cat sim1000G/include-filelist.txt | sed -e "s/^/sim1000G\//"`

check:
	R CMD check .


pdf:
	R CMD Rd2pdf .

install:
	cd .. && R CMD INSTALL sim1000G

package:
	cd .. && R CMD build sim1000G






## Development commands:


sync_to_fs:
	cd .. && rsync sim1000G/ -av ~/fs/pedigree/package/sim1000G/


