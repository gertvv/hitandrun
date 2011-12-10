#include <R.h>
#include <Rmath.h>
#include <R_ext/RS.h>
#include <R_ext/BLAS.h>

#include <Rinternals.h>
#include <R_ext/Rdynload.h>

typedef struct Matrix {
	double * const data;
	int const nRow;
	int const nCol;
} Matrix;

/**
 * Get an element of a matrix.
 * @param i Row index.
 * @param j Column index.
 */
inline double *get(Matrix *m, int i, int j) {
	return m->data + j * (m->nRow) + i;
}

/**
 * Write a row to a matrix.
 */
inline void writeRow(Matrix *m, int i, double *x) {
	for (int j = 0; j < m->nCol; ++j) {
		*get(m, i, j) = x[j];
	}
}

/**
 * Determine if x falls within the space Ax <= b.
 * @param constr: the matrix A.
 * @param rhs: the vector b.
 * @param x: the vector x.
 * @return 1 if x is within the space, 0 otherwise.
 */
int hit(Matrix *constr, double *rhs, double *x);

/**
 * Generate a random point on the n-dimensional unit sphere.
 * This is done by generating n independent normal variates and then
 * normalizing.
 */
void randDir(double *d, int n);

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
void bound(Matrix *constr, double *rhs, double *x, double *d, double *l);

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
void har(int *_n, double *_x0,
		int *_m, double *_constr, double *rhs,
		int *_niter, int *_thin, double *_result);

/**
 * Rejection sampling from the convex region Ax <= b, by generating uniform
 * variates over a bounding box and rejecting those outside the region.
 * Coordinate vectors are given in homogeneous coordinates (i.e. the last
 * coordinate is set equal to 1) so the constraint matrix can perform
 * translation etc.
 * @param n: the dimension of the sampling space.
 * @param lb: vector of lower bounds for each dimension (n components)
 * @param ub: vector of upper bounds for each dimension (n components)
 * @param m: the number of constraints
 * @param constr: the m * (n + 1) constraint matrix A
 * @param rhs: the vector b (m components): Ax <= b
 * @param niter: number of samples to generate
 * @param result: pre-allocated N * (n + 1) result matrix
 * @param reject: mean rejection rate (output)
 */
void bbReject(int *_n, double *lb, double *ub,
		int *_m, double *_constr, double *rhs,
		int *_niter, double *_result, double *reject);

/**
 * Uniform sampling from the (n-1)-simplex in n-dimensional space.
 * Points on the simplex are vectors x of n components, x_i >=0 and
 * \sum_i x_i = 1.
 * Optionally the x_i can be sorted so that x_i >= x_{i+1}.
 * @param n: the dimension of the sampling space.
 * @param niter: number of samples to generate
 * @param result: pre-allocated N * (n + 1) result matrix
 * @param sort: whether to sort the x_i
 */
void simplexSample(int *_n, int *_sort, int *_niter, double *_result);
