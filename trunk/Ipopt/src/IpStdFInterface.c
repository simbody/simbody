/********************************************************************
   Copyright (C) 2004, 2006 International Business Machines and others.
   All Rights Reserved.
   This code is published under the Common Public License.
 
   $Id: IpStdFInterface.c 759 2006-07-07 03:07:08Z andreasw $
 
   Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-03
 ********************************************************************/

#include "IpStdCInterface.h"
#include "IpoptConfig.h"
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

/* ToDo: The following needs to be adapted based on configuration */
typedef FORTRAN_INTEGER_TYPE fint;
typedef double fdouble;
typedef long fptr;

/** Return value for indicating that evaluation could be done without
    problem. */
static const fint OKRetVal = 0;
static const fint NotOKRetVal = 1;

/* Function pointer types for the Fortran callback functions */
typedef void (*FEval_F_CB)(fint* N, fdouble* X, fint* NEW_X,
                           fdouble* OBJVAL, fint* IDAT, fdouble* DDAT,
                           fint* IERR);
typedef void (*FEval_G_CB)(fint* N, fdouble* X, fint* NEW_X,
                           fint* M, fdouble* G, fint* IDAT, fdouble* DDAT,
                           fint* IERR);
typedef void (*FEval_Grad_F_CB)(fint *N, fdouble* X, fint* NEW_X,
                                fdouble* GRAD, fint* IDAT, fdouble* DDAT,
                                fint* IERR);
typedef void (*FEval_Jac_G_CB)(fint* TASK, fint* N, fdouble* X, fint* NEW_X,
                               fint* M, fint* NNZJAC, fint* IROW, fint* JCOL,
                               fdouble* VALUES, fint* IDAT, fdouble* DDAT,
                               fint* IERR);
typedef void (*FEval_Hess_CB)(fint* TASK, fint* N, fdouble* X, fint* NEW_X,
                              fdouble *OBJFACT, fint* M, fdouble* LAMBDA,
                              fint* NEW_LAM, fint* NNZHESS, fint* IROW,
                              fint* JCOL, fdouble* VALUES, fint* IDAT,
                              fdouble* DDAT, fint* IERR);

struct _FUserData
{
  fint* IDAT;
  fdouble* DDAT;
  FEval_F_CB EVAL_F;
  FEval_G_CB EVAL_G;
  FEval_Grad_F_CB EVAL_GRAD_F;
  FEval_Jac_G_CB EVAL_JAC_G;
  FEval_Hess_CB EVAL_HESS;
  IpoptProblem Problem;
};

typedef struct _FUserData FUserData;

static Bool eval_f(Index n, Number* x, Bool new_x,
                   Number* obj_value, UserDataPtr user_data)
{
  fint N = n;
  fint NEW_X = new_x;
  FUserData* fuser_data = (FUserData*)user_data;
  fint* IDAT = fuser_data->IDAT;
  fdouble* DDAT = fuser_data->DDAT;
  fint IERR=0;

  fuser_data->EVAL_F(&N, x, &NEW_X, obj_value, IDAT, DDAT, &IERR);

  return (Bool) (IERR==OKRetVal);
}

static Bool eval_grad_f(Index n, Number* x, Bool new_x,
                        Number* grad_f, UserDataPtr user_data)
{
  fint N = n;
  fint NEW_X = new_x;
  FUserData* fuser_data = (FUserData*)user_data;
  fint* IDAT = fuser_data->IDAT;
  fdouble* DDAT = fuser_data->DDAT;
  fint IERR=0;

  fuser_data->EVAL_GRAD_F(&N, x, &NEW_X, grad_f, IDAT, DDAT, &IERR);

  return (Bool) (IERR==OKRetVal);
}

static Bool eval_g(Index n, Number* x, Bool new_x,
                   Index m, Number* g, UserDataPtr user_data)
{
  fint N = n;
  fint NEW_X = new_x;
  fint M = m;
  FUserData* fuser_data = (FUserData*)user_data;
  fint* IDAT = fuser_data->IDAT;
  fdouble* DDAT = fuser_data->DDAT;
  fint IERR=0;

  fuser_data->EVAL_G(&N, x, &NEW_X, &M, g, IDAT, DDAT, &IERR);

  return (Bool) (IERR==OKRetVal);
}

static Bool eval_jac_g(Index n, Number *x, Bool new_x, Index m, Index nele_jac,
                       Index *iRow, Index *jCol, Number *values,
                       UserDataPtr user_data)
{
  fint N = n;
  fint NEW_X = new_x;
  fint M = m;
  fint NNZJAC = nele_jac;
  fint TASK;
  FUserData* fuser_data = (FUserData*)user_data;
  fint* IDAT = fuser_data->IDAT;
  fdouble* DDAT = fuser_data->DDAT;
  fint IERR=0;

  if(iRow && jCol && !values) {
    /* Only request the structure */
    TASK = 0;
  }
  else if (!iRow && !jCol && values) {
    /* Only request the values */
    TASK = 1;
  }
  else {
    printf("Error in IpStdFInterface eval_jac_g!\n");
    return (Bool) 0;
  }

  fuser_data->EVAL_JAC_G(&TASK, &N, x, &NEW_X, &M, &NNZJAC, iRow, jCol,
                         values, IDAT, DDAT, &IERR);

  return (Bool) (IERR==OKRetVal);
}

static Bool eval_h(Index n, Number *x, Bool new_x, Number obj_factor,
                   Index m, Number *lambda, Bool new_lambda,
                   Index nele_hess, Index *iRow, Index *jCol,
                   Number *values, UserDataPtr user_data)
{
  fint N = n;
  fint NEW_X = new_x;
  fint M = m;
  fint NEW_LAM = new_lambda;
  fint NNZHESS = nele_hess;
  fint TASK;
  FUserData* fuser_data = (FUserData*)user_data;
  fint* IDAT = fuser_data->IDAT;
  fdouble* DDAT = fuser_data->DDAT;
  fint IERR=0;

  if(iRow && jCol && !values) {
    /* Only request the structure */
    TASK = 0;
  }
  else if (!iRow && !jCol && values) {
    /* Only request the values */
    TASK = 1;
  }
  else {
    printf("Error in IpStdFInterface eval_hess!\n");
    return (Bool) 0;
  }

  fuser_data->EVAL_HESS(&TASK, &N, x, &NEW_X, &obj_factor,
                        &M, lambda, &NEW_LAM, &NNZHESS, iRow, jCol,
                        values, IDAT, DDAT, &IERR);

  return (Bool) (IERR==OKRetVal);
}

fptr F77_FUNC(ipcreate,IPCREATE)
(fint* N,
 fdouble* X_L,
 fdouble* X_U,
 fint* M,
 fdouble* G_L,
 fdouble* G_U,
 fint* NELE_JAC,
 fint* NELE_HESS,
 fint* IDX_STY,
 FEval_F_CB EVAL_F,
 FEval_G_CB EVAL_G,
 FEval_Grad_F_CB EVAL_GRAD_F,
 FEval_Jac_G_CB EVAL_JAC_G,
 FEval_Hess_CB EVAL_HESS)
{
  Index n = *N;
  Index m = *M;
  Index nele_jac = *NELE_JAC;
  Index nele_hess = *NELE_HESS;
  Index index_style = *IDX_STY;

  FUserData* fuser_data;

  fuser_data = (FUserData*) malloc(sizeof(FUserData));

  /* First create a new IpoptProblem object; if that fails return 0 */
  fuser_data->Problem =
    CreateIpoptProblem(n, X_L, X_U, m, G_L, G_U, nele_jac, nele_hess, 
		       index_style, eval_f, eval_g, eval_grad_f, 
		       eval_jac_g, eval_h);
  if (fuser_data->Problem == NULL) {
    free(fuser_data);
    return (fptr)NULL;
  }

  /* Store the information for the callback function */
  fuser_data->EVAL_F = EVAL_F;
  fuser_data->EVAL_G = EVAL_G;
  fuser_data->EVAL_GRAD_F = EVAL_GRAD_F;
  fuser_data->EVAL_JAC_G = EVAL_JAC_G;
  fuser_data->EVAL_HESS = EVAL_HESS;

  return (fptr)fuser_data;
}

void F77_FUNC(ipfree,IPFREE)
(fptr* FProblem)
{
  FUserData* fuser_data = (FUserData*) *FProblem;

  FreeIpoptProblem(fuser_data->Problem);
  free(fuser_data);

  *FProblem = (fptr)NULL;
}

fint F77_FUNC(ipsolve,IPSOLVE)
(fptr* FProblem,
 fdouble* X,
 fdouble* G,
 fdouble* OBJ_VAL,
 fdouble* MULT_G,
 fdouble* MULT_X_L,
 fdouble* MULT_X_U,
 fint* IDAT,
 fdouble* DDAT)
{
  FUserData* fuser_data = (FUserData*) *FProblem;
  UserDataPtr user_data;

  fuser_data->IDAT = IDAT;
  fuser_data->DDAT = DDAT;
  user_data = (UserDataPtr) fuser_data;

  return (fint)IpoptSolve(fuser_data->Problem, X, G, OBJ_VAL,
                          MULT_G, MULT_X_L, MULT_X_U, user_data);
}

static char* f2cstr(char* FSTR, int slen)
{
  int len;
  char* cstr;
  for (len=slen;len>0;len--) {
    if (FSTR[len-1]!=' ') {
      break;
    }
  }
  cstr = (char*)malloc(sizeof(char)*(len+1));
  strncpy(cstr, FSTR, len);
  cstr[len]='\0';

  return cstr;
}

/* ToDo make sure position of vlen and klen are at the right place */
fint F77_FUNC(ipaddstroption,IPADDSTROPTION)
(fptr* FProblem,
 char* KEYWORD,
 char* VALUE,
 int klen,
 int vlen)
{
  char* keyword;
  char* val;
  FUserData* fuser_data = (FUserData*) *FProblem;
  fint retval;

  keyword = f2cstr(KEYWORD, klen);
  val = f2cstr(VALUE, vlen);

  retval = AddIpoptStrOption(fuser_data->Problem, keyword, val);

  free(val);
  free(keyword);

  if (retval) {
    return OKRetVal;
  }
  else {
    return NotOKRetVal;
  }
}

fint F77_FUNC(ipaddnumoption,IPADDNUMOPTION)
(fptr* FProblem,
 char* KEYWORD,
 fdouble* VALUE,
 int klen)
{
  char* keyword;
  FUserData* fuser_data = (FUserData*) *FProblem;
  fint retval;

  keyword = f2cstr(KEYWORD, klen);

  retval = AddIpoptNumOption(fuser_data->Problem, keyword, *VALUE);

  free(keyword);

  if (retval) {
    return OKRetVal;
  }
  else {
    return NotOKRetVal;
  }
}

fint F77_FUNC(ipaddintoption,IPADDINTOPTION)
(fptr* FProblem,
 char* KEYWORD,
 fint* VALUE,
 int klen)
{
  char* keyword;
  FUserData* fuser_data = (FUserData*) *FProblem;
  Int value = *VALUE;
  fint retval;

  keyword = f2cstr(KEYWORD, klen);

  retval = AddIpoptIntOption(fuser_data->Problem, keyword, value);

  free(keyword);

  if (retval) {
    return OKRetVal;
  }
  else {
    return NotOKRetVal;
  }
}

fint F77_FUNC(ipopenoutputfile,IPOPENOUTPUTFILE)
(fptr* FProblem,
 char* FILENAME,
 fint* PRINTLEVEL,
 int flen)
{
  char* filename;
  FUserData* fuser_data = (FUserData*) *FProblem;
  Int printlevel = *PRINTLEVEL;
  fint retval;

  filename = f2cstr(FILENAME, flen);

  retval = OpenIpoptOutputFile(fuser_data->Problem, filename, printlevel);

  free(filename);

  if (retval) {
    return OKRetVal;
  }
  else {
    return NotOKRetVal;
  }
}
