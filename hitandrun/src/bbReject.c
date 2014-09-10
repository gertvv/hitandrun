#include "har.h"

SEXP hitandrun_bbReject(SEXP _lb, SEXP _ub, SEXP _constr, SEXP _rhs, SEXP _niter) {
	// get problem dimensions
	int const niter = asInteger(_niter);
	int const n = length(_lb);
	int const m = length(_rhs);
	int const nh = n + 1;

	// convert input vectors / matrices
	_lb = PROTECT(coerceVector(_lb, REALSXP));
	_ub = PROTECT(coerceVector(_ub, REALSXP));
	_constr = PROTECT(coerceVector(_constr, REALSXP));
	_rhs = PROTECT(coerceVector(_rhs, REALSXP));
	double *lb = REAL(_lb);
	double *ub = REAL(_ub);
	Matrix constr = { REAL(_constr), m, nh };
	double *rhs = REAL(_rhs);

	// allocate output matrix
	SEXP _result = PROTECT(allocMatrix(REALSXP, niter, nh));
	Matrix result = { REAL(_result), niter, nh };
	double reject = 0.0;

	// internal state variables
	double d[n]; // pre-calculated differences
	for (int j = 0; j < n; ++j) {
		d[j] = ub[j] - lb[j];
	}
	double x[n + 1];
	x[n] = 1.0;

	GetRNGstate(); // enable use of RNGs

	for (int i = 0; i < niter; ++i) {
		int wasHit = 0;
		int miss = 0;
		while (!wasHit) {
			for (int j = 0; j < n; ++j) {
				x[j] = lb[j] + d[j] * unif_rand();
			}
			if (hitandrun_hit(&constr, rhs, x, 0.0)) {
				wasHit = 1;
				writeRow(&result, i, x);
			} else {
				++miss;
			}
		}
		reject += (double)miss / niter;
	}

	PutRNGstate();

	// return samples and rejection rate as a list
	SEXP ans = PROTECT(allocVector(VECSXP, 2));
	SET_VECTOR_ELT(ans, 0, _result);
	SET_VECTOR_ELT(ans, 1, ScalarReal(reject));

	UNPROTECT(6);
	return ans;
}
