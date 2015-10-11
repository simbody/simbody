/*
 * -----------------------------------------------------------------
 * $Revision: 1.3 $
 * $Date: 2006/11/22 00:12:51 $
 * -----------------------------------------------------------------
 * Programmer: Radu Serban @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2002, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation file for a generic package of dense
 * matrix operations based on BLAS/LAPACK.
 * -----------------------------------------------------------------
 */

#include <stdio.h>

#include <sundials/sundials_lapack.h>
#include <sundials/sundials_math.h>

#define ZERO RCONST(0.0)
#define ONE  RCONST(1.0)


/*
 * LapackDenseZero sets all elements of the dense matrix A to 0
 */
void LapackDenseZero(DlsMat A)
{
  int i;
  for (i=0; i<A->ldata; i++) A->data[i] = ZERO;
}

/*
 * LapackBandZero sets all elements of the band matrix A to 0
 */
void LapackBandZero(DlsMat A)
{
  int i;
  for (i=0; i<A->ldata; i++) A->data[i] = ZERO;
}

/*
 * LapackDenseAddI overwrites the dense matrix A with I+A
 */
void LapackDenseAddI(DlsMat A)
{
  int j;
  realtype *col_j;
  for (j=0; j<A->N; j++) {
    col_j = A->cols[j];
    col_j[j] += ONE;
  }
}

/*
 * LapackBandAddI overwrites the banded matrix A with I+A
 */
void LapackBandAddI(DlsMat A)
{
  int j;
  realtype *col_j;
  for (j=0; j<A->N; j++) {
    col_j = (A->cols)[j] + A->s_mu;
    col_j[0] += ONE;
  }
}
