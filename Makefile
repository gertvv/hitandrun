PACKAGE=hitandrun_0.2.tar.gz

$(PACKAGE): hitandrun/src/*.c hitandrun/src/*.h hitandrun/R/* hitandrun/man/* hitandrun/DESCRIPTION hitandrun/NAMESPACE
	rm hitandrun/src/*.o hitandrun/src/*.so
	R CMD build hitandrun
	R CMD check hitandrun

.PHONY: install

install: $(PACKAGE)
	R CMD INSTALL hitandrun
