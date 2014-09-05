#include "har.h"

/**
 * @param constr: constraint matrix.
 * @param rhs: right hand side of the constraints.
 * @param x: current point.
 * @param d: direction to move in.
 * @param l: (OUTPUT) how far to move.
 * @return the index of the intersected constraint.
 */
int hitandrun_intersect(Matrix *constr, double *rhs, double *x, double *d, double *l, int prev) {
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
	// A(x + td) <= b, i.e. t(Ad) <= b - Ax, and find the maximum:
	//
	//   min_{i:(Ad)_i>0} (b - Ax)_i / (Ad)_i
	*l = NA_REAL;
	int idx;
	for (int i = 0; i < constr->nRow; ++i) {
		if (i != prev && c[i] > 0.0) {
			double t = a[i] / c[i];
			if (ISNA(*l) || t < *l) {
				*l = t;
				idx = i;
			}
		}
	}

	return idx;
}
