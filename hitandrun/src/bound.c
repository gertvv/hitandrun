#include "har.h"

void hitandrun_bound(Matrix *constr, double *rhs, double *x, double *d, double *l) {
	const int inc1 = 1;
	const double one = 1.0, negone = -1.0, zero = 0.0; // for BLAS
	const char trans = 'N';

	double a[constr->nRow];
	memcpy(a, rhs, constr->nRow * sizeof(double));
	F77_CALL(dgemv)(&trans, &(constr->nRow), &(constr->nCol),
		&negone, constr->data, &(constr->nRow), x, &inc1,
		&one, a, &inc1); // a := -1Ax + 1a (= b - Ax)

	double c[constr->nRow];
	F77_CALL(dgemv)(&trans, &(constr->nRow), &(constr->nCol),
		&one, constr->data, &(constr->nRow), d, &inc1,
		&zero, c, &inc1); // c := 1Ad + 0c

	// we know Ax <= b, now we need to find the values of t such that
	// A(x + td) <= b, i.e. t(Ad) <= b - Ax:
	// T = [
	//   max_{i:(Ad)_i<0} (b - Ax)_i / (Ad)_i,
	//   min_{i:(Ad)_i>0} (b - Ax)_i / (Ad)_i
	// ]
	l[0] = l[1] = NA_REAL;
	for (int i = 0; i < constr->nRow; ++i) {
		if (c[i] < 0.0) {
			double t = a[i] / c[i];
			if (ISNA(l[0]) || t > l[0]) {
				l[0] = t;
			}
		} else if (c[i] > 0.0) {
			double t = a[i] / c[i];
			if (ISNA(l[1]) || t < l[1]) {
				l[1] = t;
			}
		}
	}
}
