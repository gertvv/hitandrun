#include "har.h"

void bbReject(int *_n, double *lb, double *ub,
		int *_m, double *_constr, double *rhs,
		int *_niter, double *_result, double *reject) {
	int const n = *_n, m = *_m, niter = *_niter;
	Matrix constr = { _constr, m, n + 1 };
	Matrix result = { _result, niter, n + 1 };

	double d[n]; // pre-calculate differences
	for (int j = 0; j < n; ++j) {
		d[j] = ub[j] - lb[j];
	}

	double x[n + 1];
	x[n] = 1.0;
	*reject = 0.0;

	GetRNGstate(); // enable use of RNGs

	for (int i = 0; i < niter; ++i) {
		int wasHit = 0;
		int miss = 0;
		while (!wasHit) {
			for (int j = 0; j < n; ++j) {
				x[j] = lb[j] + d[j] * unif_rand();
			}
			if (hit(&constr, rhs, x)) {
				wasHit = 1;
				writeRow(&result, i, x);
			} else {
				++miss;
			}
		}
		*reject += (double)miss / niter;
	}

	PutRNGstate();
}
