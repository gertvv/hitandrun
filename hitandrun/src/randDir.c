#include "har.h"

void randDir(double *d, int n) {
	int const inc1 = 1;
	for (int i = 0; i < n; ++i) {
		d[i] = norm_rand();
	}
	double f = 1 / F77_CALL(dnrm2)(&n, d, &inc1);
	F77_CALL(dscal)(&n, &f, d, &inc1);
}

void randDirForR(int *_n, int *_niter, double *_result) {
	int n = *_n;
	int niter = *_niter;
	Matrix result = { _result, niter, n };

	GetRNGstate(); // enable use of RNGs

	double x[n];
	for (int i = 0; i < niter; ++i) {
		randDir(x, n);
		writeRow(&result, i, x);
	}

	PutRNGstate();
}
