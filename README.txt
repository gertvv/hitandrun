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

Dependencies
------------

hitandrun depends on rcdd, which requires GNU MP. On debian derivatives, the
compile-time dependency is provided by libgmp-dev.
