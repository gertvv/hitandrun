#include "har.h"

int hitandrun_hit(Matrix *constr, double *rhs, double *x, double epsilon) {
	const int inc1 = 1;
	const double one = 1.0, zero = 0.0; // for BLAS
	const char trans = 'N';

	double a[constr->nRow];
	F77_CALL(dgemv)(&trans, &(constr->nRow), &(constr->nCol),
		&one, constr->data, &(constr->nRow), x, &inc1,
		&zero, a, &inc1); // a := 1Ax + 0a

	for (int i = 0; i < constr->nRow; ++i) {
		if (a[i] - rhs[i] > epsilon) {
			return 0;
		}
	}
	return 1;
}
