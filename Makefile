
PKGNAME := $(shell sed -n "s/Package: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGVERS := $(shell sed -n "s/Version: *\([^ ]*\)/\1/p" DESCRIPTION)
PKGSRC  := $(shell basename `pwd`)

default:
	Rscript  -e 'library(devtools); document(); build(); '

build:
	cd ..;\
	R CMD build $(PKGSRC)

document:
	Rscript -e 'tools::buildVignettes(dir = ".") ' 
	#Rscript -e 'devtools::build()' 
	Rscript  -e 'library(devtools); document(); '

checkcitation:
	R CMD INSTALL .
	Rscript -e 'citation("sim1000G")'


tgz:
	cd .. && tar cvzf sim1000G.tar.gz `cat sim1000G/include-filelist.txt | sed -e "s/^/sim1000G\//"`


check: build
	cd ..;\
	R CMD check $(PKGNAME)_$(PKGVERS).tar.gz --as-cran


pdf:
	R CMD Rd2pdf .

install:
	cd .. && R CMD INSTALL sim1000G

package:
	cd .. && R CMD build sim1000G






## Development commands:


sync_to_fs:
	cd .. && rsync sim1000G/ -av ~/fs/pedigree/package/sim1000G/


