read_version = $(shell grep 'Version:' $1/DESCRIPTION | sed 's/Version: //')

PKG_NAME := hitandrun
PKG_VERSION := $(call read_version,$(PKG_NAME))
PACKAGE := $(PKG_NAME)_$(PKG_VERSION).tar.gz

all: $(PACKAGE)

$(PACKAGE): hitandrun/src/*.c hitandrun/src/*.h hitandrun/R/* hitandrun/man/* hitandrun/tests/* hitandrun/DESCRIPTION hitandrun/NAMESPACE
	rm -f hitandrun/src/*.o hitandrun/src/*.so
	R CMD build hitandrun

.PHONY: install clean check test

clean:
	rm -f hitandrun/src/*.o hitandrun/src/*.so hitandrun/src/symbols.rds

install: $(PACKAGE)
	R CMD INSTALL $(PACKAGE)

check: $(PACKAGE)
	R CMD check $(PACKAGE)

test: $(PACKAGE)
	./run-tests.sh $(PACKAGE)
