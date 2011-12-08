#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>

typedef struct Matrix {
    double * const data;
    int const nRow;
    int const nCol;
} Matrix;

/**
 * @param i Row index.
 * @param j Column index.
 */
inline double *get(Matrix *m, int i, int j) {
    return m->data + j * (m->nRow) + i;
}

/**
 * Generate a random point on the n-dimensional unit sphere.
 * This is done by generating n independent normal variates and then
 * normalizing.
 */
void randDir(double *d, int n) {
	int const inc1 = 1;
	for (int i = 0; i < n; ++i) {
		d[i] = norm_rand();
	}
	double f = 1 / F77_CALL(dnrm2)(&n, d, &inc1);
	F77_CALL(dscal)(&n, &f, d, &inc1);
}

/**
 * Determine if x falls within the space Ax <= b.
 * @param constr: the matrix A.
 * @param rhs: the vector b.
 * @param x: the vector x.
 * @return 1 if x is within the space, 0 otherwise.
 */
int hit(Matrix *constr, double *rhs, double *x) {
	const int inc1 = 1;
	const double one = 1.0, zero = 0.0; // for BLAS
	const char trans = 'N';

	double a[constr->nRow];
	F77_CALL(dgemv)(&trans, &(constr->nRow), &(constr->nCol),
		&one, constr->data, &(constr->nRow), x, &inc1,
		&zero, a, &inc1); // a := 1Ax + 0a

	for (int i = 0; i < constr->nRow; ++i) {
		if (a[i] > rhs[i]) {
			return 0;
		}
	}
	return 1;
}

/**
 * Give bounds for how far we can move from x in the direction of d without
 * exiting the space Ax <= b.
 * @param constr: the matrix A.
 * @param rhs: the vector b.
 * @param x: the vector x.
 * @param d: the direction d.
 * @param l: output -- how far we can move in the negative resp. positive
 * direction.
 */
void bound(Matrix *constr, double *rhs, double *x, double *d, double *l) {
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

#include <stdio.h>

/**
 * Hit-and-Run sampling using exact intersections with the linear constraints.
 * It samples from the convex region Ax <= b.
 * Coordinate vectors are given in homogeneous coordinates (i.e. the last
 * coordinate is set equal to 1) so the constraint matrix can perform
 * translation etc. The generated jump direction will have its last coordinate
 * equal to 0.
 * @param n: the dimension of the sampling space.
 * @param x0: vector representing the seed point (n + 1 components)
 * @param m: the number of constraints
 * @param constr: the m * (n + 1) constraint matrix A
 * @param rhs: the vector b (m components): Ax <= b
 * @param niter: number of iterations to sample for
 * @param thin: thinning interval (will return N=niter/thin samples)
 * @param result: pre-allocated N * (n + 1) result matrix
 */
void har(int *_n, double *_x0, int *_m, double *_constr, double *rhs,
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
	if (!hit(&constr, rhs, _x0)) {
		error("The starting point must be inside the region");
	}

	double x[n + 1];
	memcpy(x, _x0, (n + 1) * sizeof(double));
	double d[n + 1];
	d[n] = 0; // homogeneous coordinates -- final direction component always 0.
	double l[2];

	GetRNGstate(); // enable use of RNGs

	for (int i = 0; i < niter; ++i) {
		randDir(d, n); // generate random direction d
		bound(&constr, rhs, x, d, l); // calculate bounds l
		if (!R_FINITE(l[0]) || !R_FINITE(l[1])) {
			error("Bounding function gave NA bounds [%f, %f]", l[0], l[1]);
		}
		if (l[0] == l[1]) { // FIXME: is this an error?
			error("Bounding function gave empty interval");
		}
		double v = l[0] + unif_rand() * (l[1] - l[0]);
		F77_CALL(daxpy)(&nh, &v, d, &inc1, x, &inc1); // x := vd + x

		if ((i + 1) % thin == 0) { // write result
			for (int j = 0; j < n + 1; ++j) {
				*get(&result, i, j) = x[j];
			}
		}
	}

	PutRNGstate(); // propagate RNG state back to R (or somewhere)
}
