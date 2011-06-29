"Hit and Run" sampler
=====================

This project provides an implementation of the "Hit and Run" algorithm
for sampling from convex shapes. Our motivation for doing this is to
enable the fast generation of arbitrarily constrained weights in
high-dimensional space. This should enable applying e.g. SMAA
(http://smaa.fi/) with more complex weigth constraints.

R implementation
----------------

Currently, an R package is available:

https://github.com/gertvv/hitandrun/downloads

See ??har, ?har, ?har-constraint and ?simplex.sample

Note that the R package is still heavily work-in-progress and we are
preparing to include elaborate examples, validation and computational
tests in the future.

Java implementation
-------------------

A Java implementation for use in JSMAA (http://smaa.fi/jsmaa.php) is
planned.
