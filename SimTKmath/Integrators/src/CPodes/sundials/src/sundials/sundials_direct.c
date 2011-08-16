/*
 * -----------------------------------------------------------------
 * $Revision: 1.2 $
 * $Date: 2006/11/29 00:05:10 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for operations to be used by a
 * generic direct linear solver.
 * -----------------------------------------------------------------
 */ 

#include <stdio.h>
#include <stdlib.h>

#include <sundials/sundials_direct.h>
#include <sundials/sundials_math.h>

DlsMat NewDenseMat(int M, int N)
{
  DlsMat A;
  int j;

  if ( (M <= 0) || (N <= 0) ) return(NULL);

  A = NULL;
  A = (DlsMat) malloc(sizeof *A);
  if (A==NULL) return (NULL);
  
  A->data = (realtype *) malloc(M * N * sizeof(realtype));
  if (A->data == NULL) {
    free(A); A = NULL;
    return(NULL);
  }
  A->cols = (realtype **) malloc(N * sizeof(realtype *));
  if (A->cols == NULL) {
    free(A->data); A->data = NULL;
    free(A); A = NULL;
    return(NULL);
  }

  for (j=0; j < N; j++) A->cols[j] = A->data + j * M;

  A->M = M;
  A->N = N;
  A->ldim = M;
  A->ldata = M*N;

  A->type = SUNDIALS_DENSE;

  return(A);
}

realtype **newDenseMat(int m, int n)
{
  int j;
  realtype **a;

  if ( (n <= 0) || (m <= 0) ) return(NULL);

  a = NULL;
  a = (realtype **) malloc(n * sizeof(realtype *));
  if (a == NULL) return(NULL);

  a[0] = NULL;
  a[0] = (realtype *) malloc(m * n * sizeof(realtype));
  if (a[0] == NULL) {
    free(a); a = NULL;
    return(NULL);
  }

  for (j=1; j < n; j++) a[j] = a[0] + j * m;

  return(a);
}


DlsMat NewBandMat(int N, int mu, int ml, int smu)
{
  DlsMat A;
  int j, colSize;

  if (N <= 0) return(NULL);
  
  A = NULL;
  A = (DlsMat) malloc(sizeof *A);
  if (A == NULL) return (NULL);

  colSize = smu + ml + 1;
  A->data = NULL;
  A->data = (realtype *) malloc(N * colSize * sizeof(realtype));
  if (A->data == NULL) {
    free(A); A = NULL;
    return(NULL);
  }

  A->cols = NULL;
  A->cols = (realtype **) malloc(N * sizeof(realtype *));
  if (A->cols == NULL) {
    free(A->data);
    free(A); A = NULL;
    return(NULL);
  }

  for (j=0; j < N; j++) A->cols[j] = A->data + j * colSize;

  A->M = N;
  A->N = N;
  A->mu = mu;
  A->ml = ml;
  A->s_mu = smu;
  A->ldim =  colSize;
  A->ldata = N * colSize;

  A->type = SUNDIALS_BAND;

  return(A);
}

realtype **newBandMat(int n, int smu, int ml)
{
  realtype **a;
  int j, colSize;

  if (n <= 0) return(NULL);

  a = NULL;
  a = (realtype **) malloc(n * sizeof(realtype *));
  if (a == NULL) return(NULL);

  colSize = smu + ml + 1;
  a[0] = NULL;
  a[0] = (realtype *) malloc(n * colSize * sizeof(realtype));
  if (a[0] == NULL) {
    free(a); a = NULL;
    return(NULL);
  }

  for (j=1; j < n; j++) a[j] = a[0] + j * colSize;

  return(a);
}

void DestroyMat(DlsMat A)
{
  free(A->data);  A->data = NULL;
  free(A->cols);
  free(A); A = NULL;
}

void destroyMat(realtype **a)
{
  free(a[0]); a[0] = NULL;
  free(a); a = NULL;
}

int *NewIntArray(int N)
{
  int *vec;

  if (N <= 0) return(NULL);

  vec = NULL;
  vec = (int *) malloc(N * sizeof(int));

  return(vec);
}

int *newIntArray(int n)
{
  int *v;

  if (n <= 0) return(NULL);

  v = NULL;
  v = (int *) malloc(n * sizeof(int));

  return(v);
}

realtype *NewRealArray(int N)
{
  realtype *vec;

  if (N <= 0) return(NULL);

  vec = NULL;
  vec = (realtype *) malloc(N * sizeof(realtype));

  return(vec);
}

realtype *newRealArray(int m)
{
  realtype *v;

  if (m <= 0) return(NULL);

  v = NULL;
  v = (realtype *) malloc(m * sizeof(realtype));

  return(v);
}

void DestroyArray(void *V)
{ 
  free(V); 
  V = NULL;
}

void destroyArray(void *v)
{
  free(v); 
  v = NULL;
}

void PrintMat(DlsMat A)
{
  int i, j, start, finish;
  realtype **a;

  switch (A->type) {

  case SUNDIALS_DENSE:

    printf("\n");
    for (i=0; i < A->M; i++) {
      for (j=0; j < A->N; j++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
        printf("%12Lg  ", DENSE_ELEM(A,i,j));
#elif defined(SUNDIALS_DOUBLE_PRECISION)
        printf("%12lg  ", DENSE_ELEM(A,i,j));
#else
        printf("%12g  ", DENSE_ELEM(A,i,j));
#endif
      }
      printf("\n");
    }
    printf("\n");
    
    break;

  case SUNDIALS_BAND:

    a = A->cols;
    printf("\n");
    for (i=0; i < A->N; i++) {
      start = MAX(0,i-A->ml);
      finish = MIN(A->N-1,i+A->mu);
      for (j=0; j < start; j++) printf("%12s  ","");
      for (j=start; j <= finish; j++) {
#if defined(SUNDIALS_EXTENDED_PRECISION)
        printf("%12Lg  ", a[j][i-j+A->s_mu]);
#elif defined(SUNDIALS_DOUBLE_PRECISION)
        printf("%12lg  ", a[j][i-j+A->s_mu]);
#else
        printf("%12g  ", a[j][i-j+A->s_mu]);
#endif
      }
      printf("\n");
    }
    printf("\n");
    
    break;

  }

}


