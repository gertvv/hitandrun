#include "har.h"

SEXP hitandrun_hypersphereSample(SEXP _n, SEXP _niter) {
	int const n = asInteger(_n);
	int const niter = asInteger(_niter);

	// allocate output matrix
	SEXP _result = PROTECT(allocMatrix(REALSXP, niter, n));
	Matrix result = { REAL(_result), niter, n };

	GetRNGstate(); // enable use of RNGs

	double x[n];
	for (int i = 0; i < niter; ++i) {
		hitandrun_randDir(x, n);
		writeRow(&result, i, x);
	}

	PutRNGstate();

	UNPROTECT(1);
	return _result;
}
