#include "har.h"

void hitandrun_har(int *_n, double *_x0, int *_m, double *_constr, double *rhs,
		int *_niter, int *_thin, double *_result) {
	int const n = *_n, m = *_m, niter = *_niter, thin = *_thin;
	const int nh = n + 1; // needed for BLAS
	const int inc1 = 1; // needed for BLAS
	Matrix constr = { _constr, m, n + 1 };
	Matrix result = { _result, niter / thin, n + 1 };

	// Check arguments for sanity
	if (_x0[n] != 1.0) {
		error("The (n + 1)-st component of x0 must be 1");
	}
	if (niter % thin != 0) {
		error("niter % thin != 0");
	}
	if (!hitandrun_hit(&constr, rhs, _x0)) {
		error("The starting point must be inside the region");
	}

	double x[n + 1];
	memcpy(x, _x0, (n + 1) * sizeof(double));
	double d[n + 1];
	d[n] = 0; // homogeneous coordinates -- final direction component always 0.
	double l[2];

	GetRNGstate(); // enable use of RNGs

	for (int i = 0; i < niter; ++i) {
		hitandrun_randDir(d, n); // generate random direction d
		hitandrun_bound(&constr, rhs, x, d, l); // calculate bounds l
		if (!R_FINITE(l[0]) || !R_FINITE(l[1])) {
			error("Bounding function gave NA bounds [%f, %f]", l[0], l[1]);
		}
		if (l[0] == l[1]) { // FIXME: is this an error?
			error("Bounding function gave empty interval");
		}
		double v = l[0] + unif_rand() * (l[1] - l[0]);
		F77_CALL(daxpy)(&nh, &v, d, &inc1, x, &inc1); // x := vd + x

		if ((i + 1) % thin == 0) { // write result
			writeRow(&result, i / thin, x);
		}
	}

	PutRNGstate(); // propagate RNG state back to R (or somewhere)
}
