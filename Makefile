PACKAGE=hitandrun_0.2.tar.gz

$(PACKAGE): hitandrun/src/*.c hitandrun/src/*.h hitandrun/R/* hitandrun/man/* hitandrun/tests/* hitandrun/DESCRIPTION hitandrun/NAMESPACE
	rm -f hitandrun/src/*.o hitandrun/src/*.so
	R CMD build hitandrun
	R CMD check hitandrun

.PHONY: install

clean:
	rm -f hitandrun/src/*.o hitandrun/src/*.so hitandrun/src/symbols.rds

install: $(PACKAGE)
	R CMD INSTALL hitandrun

check:
	R CMD check hitandrun
