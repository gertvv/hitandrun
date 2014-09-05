#include "har.h"

/**
 * Sample from the boundary of a convex polytope using the "running
 * Shake-and-Bake" algorithm.
 * @param _x0 The starting point.
 * @param _index The index of the boundary (constraint) _x0 lies on.
 * @param _constr The constraint matrix (normalized).
 * @param _rhs The right hand side of the constraints.
 * @param _niter The total number of iterations to run.
 * @param _thin The thinning interval.
 * @return _niter / _thin samples.
 */
SEXP hitandrun_sab(SEXP _x0, SEXP _index, SEXP _constr, SEXP _rhs, SEXP _niter, SEXP _thin) {
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
	int index = asInteger(_index);

	// check the starting point
	if (!hitandrun_hit(&constr, rhs, x0)) {
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
	double l; // length to move

	GetRNGstate(); // enable use of RNGs

	for (int i = 0; i < niter; ++i) {
		hitandrun_rsabDir(d, &constr, index); // generate random direction d

		// calculate the intersection with the boundary
		index = hitandrun_intersect(&constr, rhs, x, d, &l, index);
		F77_CALL(daxpy)(&nh, &l, d, &inc1, x, &inc1); // x := ld + x

		if (!R_FINITE(x[0])) {
			error("Generated NA");
		}

		if ((i + 1) % thin == 0) { // write result
			writeRow(&result, i / thin, x);
		}
	}

	PutRNGstate(); // propagate RNG state back to R (or somewhere)

	UNPROTECT(4);
	return _result;
}
