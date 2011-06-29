PACKAGE=har_0.1.tar.gz

$(PACKAGE): har/R/* har/man/* har/DESCRIPTION
	R CMD build har
	R CMD check har

.PHONY: install

install: $(PACKAGE)
	R CMD INSTALL har
