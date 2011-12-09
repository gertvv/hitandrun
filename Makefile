PACKAGE=hitandrun_0.1.tar.gz

$(PACKAGE): hitandrun/src/* hitandrun/R/* hitandrun/man/* hitandrun/DESCRIPTION hitandrun/NAMESPACE
	R CMD build hitandrun
	R CMD check hitandrun

.PHONY: install

install: $(PACKAGE)
	R CMD INSTALL hitandrun
