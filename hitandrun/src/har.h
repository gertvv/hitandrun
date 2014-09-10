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
static inline double *get(Matrix *m, int i, int j) {
	return m->data + j * (m->nRow) + i;
}

/**
 * Write a row to a matrix.
 */
static inline void writeRow(Matrix *m, int i, double *x) {
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
int hitandrun_hit(Matrix *constr, double *rhs, double *x, double epsilon);

/**
 * Generate a random point on the n-dimensional unit sphere.
 * This is done by generating n independent normal variates and then
 * normalizing.
 * @param d pointer to an array of size n
 * @param n size of the array d
 */
void hitandrun_randDir(double *d, int n);

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
void hitandrun_bound(Matrix *constr, double *rhs, double *x, double *d, double *l);

/**
 * Hit-and-Run sampling using exact intersections with the linear constraints.
 * It samples from the convex region Ax <= b.
 * Coordinate vectors are given in homogeneous coordinates (i.e. the last
 * coordinate is set equal to 1) so the constraint matrix can perform
 * translation etc. The generated jump direction will have its last coordinate
 * equal to 0.
 * @param x0: vector representing the seed point (n + 1 components)
 * @param constr: the m * (n + 1) constraint matrix A
 * @param rhs: the vector b (m components): Ax <= b
 * @param niter: number of iterations to sample for
 * @param thin: thinning interval (will return N=niter/thin samples)
 * @return N * (n + 1) result matrix
 */
SEXP hitandrun_har(SEXP x0, SEXP constr, SEXP rhs, SEXP niter, SEXP thin);

/**
 * Rejection sampling from the convex region Ax <= b, by generating uniform
 * variates over a bounding box and rejecting those outside the region.
 * Coordinate vectors are given in homogeneous coordinates (i.e. the last
 * coordinate is set equal to 1) so the constraint matrix can perform
 * translation etc.
 * @param lb: vector of lower bounds for each dimension (n components)
 * @param ub: vector of upper bounds for each dimension (n components)
 * @param constr: the m * (n + 1) constraint matrix A
 * @param rhs: the vector b (m components): Ax <= b
 * @param niter: number of samples to generate
 * @return A list of: (1) N * (n + 1) result matrix; (2) rejection rate
 */
SEXP hitandrun_bbReject(SEXP lb, SEXP ub, SEXP constr, SEXP rhs, SEXP niter);

/**
 * Uniform sampling from the (n-1)-simplex in n-dimensional space.
 * Points on the simplex are vectors x of n components, x_i >=0 and
 * \sum_i x_i = 1.
 * Optionally the x_i can be sorted so that x_i >= x_{i+1}.
 * @param n: the dimension of the sampling space.
 * @param sort: whether to sort the x_i
 * @param niter: number of samples to generate
 * @return N * (n + 1) result matrix
 */
SEXP hitandrun_simplexSample(SEXP _n, SEXP _sort, SEXP _niter);

/**
 * Uniform sampling from the boundary of the n-sphere.
 * @param n Dimensionality of the hypersphere
 * @param N The number of samples
 * @return n * N result matrix
 */
SEXP hitandrun_hypersphereSample(SEXP n, SEXP N);

// Shake-and-Bake stuff

int hitandrun_intersect(Matrix *constr, double *rhs, double *x, double *d, double *l, int prev);

void hitandrun_rsabDir(double * d, Matrix *constr, int index);

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
SEXP hitandrun_sab(SEXP _x0, SEXP _index, SEXP _constr, SEXP _rhs, SEXP _niter, SEXP _thin);
