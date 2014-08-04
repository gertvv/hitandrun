#include "har.h"

void hitandrun_randDir(double *d, int n) {
	int const inc1 = 1;
	for (int i = 0; i < n; ++i) {
		d[i] = norm_rand();
	}
	double f = 1 / F77_CALL(dnrm2)(&n, d, &inc1);
	F77_CALL(dscal)(&n, &f, d, &inc1);
}
