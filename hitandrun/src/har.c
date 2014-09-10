#include "har.h"

SEXP hitandrun_har(SEXP _x0, SEXP _constr, SEXP _rhs, SEXP _niter, SEXP _thin) {
	// get problem dimensions
	int const niter = asInteger(_niter);
	int const thin = asInteger(_thin);
	int const n = length(_x0) - 1;
	int const m = length(_rhs);
	int const nh = n + 1; // needed for BLAS
	int const inc1 = 1; // needed for BLAS

	// convert input vectors / matrices
	_x0 = PROTECT(coerceVector(_x0, REALSXP));
	_constr = PROTECT(coerceVector(_constr, REALSXP));
	_rhs = PROTECT(coerceVector(_rhs, REALSXP));
	double *x0 = REAL(_x0);
	double *rhs = REAL(_rhs);
	Matrix constr = { REAL(_constr), m, nh };

	// check the starting point
	if (!hitandrun_hit(&constr, rhs, x0, 0.0)) {
		UNPROTECT(3);
		error("The starting point must be inside the region");
	}

	// allocate output matrix
	SEXP _result = PROTECT(allocMatrix(REALSXP, niter / thin, nh));
	Matrix result = { REAL(_result), niter / thin, nh };

	// state variables
	double x[n + 1];
	memcpy(x, x0, (n + 1) * sizeof(double));
	double d[n + 1];
	d[n] = 0; // homogeneous coordinates -- final direction component always 0.
	double l[2];

	GetRNGstate(); // enable use of RNGs

	for (int i = 0; i < niter; ++i) {
		hitandrun_randDir(d, n); // generate random direction d
		hitandrun_bound(&constr, rhs, x, d, l); // calculate bounds l
		if (!R_FINITE(l[0]) || !R_FINITE(l[1])) {
			UNPROTECT(4);
			error("Bounding function gave NA bounds [%f, %f]", l[0], l[1]);
		}
		if (l[0] == l[1]) { // FIXME: is this an error?
			UNPROTECT(4);
			error("Bounding function gave empty interval");
		}
		double v = l[0] + unif_rand() * (l[1] - l[0]);
		F77_CALL(daxpy)(&nh, &v, d, &inc1, x, &inc1); // x := vd + x

		if ((i + 1) % thin == 0) { // write result
			writeRow(&result, i / thin, x);
		}
	}

	PutRNGstate(); // propagate RNG state back to R (or somewhere)

	UNPROTECT(4);
	return _result;
}
