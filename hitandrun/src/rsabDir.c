#include "har.h"

inline double root(double x, int d) {
	return pow(x, 1.0/d);
}

/**
 * Generate the direction for "running Shake-and-Bake" according to 1.3.3 of
 * Boender et al. (1991)
 */
void hitandrun_rsabDir(double *d, Matrix *constr, int index) {
	const int inc1 = 1; // for BLAS

	int n = constr->nCol - 1;

	double c[n]; // the constraint vector
	for (int i = 0; i < n; ++i) {
		c[i] = *get(constr, index, i);
	}

  if (n == 1) {
    d[0] = -c[0];
    return;
  }

	double r = root(unif_rand(), n - 1);
	hitandrun_randDir(d, n); // \~{u} in the paper

	double cd = F77_CALL(ddot)(&n, c, &inc1, d, &inc1);
	double fd = r / sqrt(1 - cd * cd);
	double fc = -(r * cd / sqrt(1 - cd * cd) + sqrt(1 - r * r));
	F77_CALL(dscal)(&n, &fd, d, &inc1); // d := fd * d
	F77_CALL(daxpy)(&n, &fc, c, &inc1, d, &inc1); // d := fc * c + d
}
