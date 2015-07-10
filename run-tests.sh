#!/bin/bash

DIR=`mktemp -d -t hitandrun.XXXXXX`
R CMD INSTALL -l $DIR --install-tests $1

R --vanilla --slave -e "library(hitandrun, lib.loc=\"$DIR\"); library(testthat); test_package('hitandrun')"
rm -rf $DIR
