#include "har.h"

void randDir(double *d, int n) {
  int const inc1 = 1;
  for (int i = 0; i < n; ++i) {
    d[i] = norm_rand();
  }
  double f = 1 / F77_CALL(dnrm2)(&n, d, &inc1);
  F77_CALL(dscal)(&n, &f, d, &inc1);
}

void randDirForR(double *result, int *_n, int *_N) {
  int n = *_n;
  int N = *_N;

  Matrix resMat = {result, N, n};

  GetRNGstate(); // enable use of RNGs

  double * curSample = (double *) malloc (sizeof(double) * n);
  for (int i=0;i<N;i++) {
    randDir(curSample, n);
    writeRow(&resMat, i, curSample);
  }
  free(curSample);

  PutRNGstate();
}
