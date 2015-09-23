/*
 * -----------------------------------------------------------------
 * $Revision:$
 * $Date:$
 * -----------------------------------------------------------------
 * Programmer: Radu Serban  @ LLNL
 * -----------------------------------------------------------------
 * Copyright (c) 2006, The Regents of the University of California.
 * Produced at the Lawrence Livermore National Laboratory.
 * All rights reserved.
 * For details, see the LICENSE file.
 * -----------------------------------------------------------------
 * This is the implementation for the event detection for CPODES.
 * -----------------------------------------------------------------
 */

/*
 * =================================================================
 * Import Header Files
 * =================================================================
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>

#include "cpodes_private.h"
#include <sundials/sundials_math.h>

/*
 * =================================================================
 * Private Function Prototypes
 * =================================================================
 */

static booleantype cpRootAlloc(CPodeMem cp_mem, int nrt);
static int cpRootfind(CPodeMem cp_mem, realtype ttol);

/*
 * =================================================================
 * EXPORTED FUNCTIONS
 * =================================================================
 */

/*
 * CPodeRootInit
 *
 * CPodeRootInit initializes a rootfinding problem to be solved
 * during the integration of the ODE system.  It loads the root
 * function pointer and the number of root functions, and allocates
 * workspace memory.  The return value is CP_SUCCESS = 0 if no errors
 * occurred, or a negative value otherwise.
 */

int CPodeRootInit(void *cpode_mem, int nrtfn, CPRootFn gfun, void *g_data)
{
  CPodeMem cp_mem;
  booleantype allocOK;
  int i;

  /* Check cpode_mem pointer */
  if (cpode_mem == NULL) {
    cpProcessError(NULL, CP_MEM_NULL, "CPODES", "CPodeRootInit", MSGCP_NO_MEM);
    return(CP_MEM_NULL);
  }
  cp_mem = (CPodeMem) cpode_mem;

  /* If called with nrtfn <= 0, then disable rootfinding and return */
  if (nrtfn <= 0) {
    cp_mem->cp_doRootfinding = FALSE;
    return(CP_SUCCESS);
  }

  /* Check for legal input parameters */
  if (gfun == NULL) {
    cpProcessError(cp_mem, CP_ILL_INPUT, "CPODES", "CPodeRootInit", MSGCP_NULL_G);
    return(CP_ILL_INPUT);
  }

  /* If rerunning CPodeRootInit() with a different number of root
   * functions (changing number of gfun components), then free
   * currently held memory resources */
  if ( (cp_mem->cp_rootMallocDone) && (nrtfn != cp_mem->cp_nrtfn) ) {
    cpRootFree(cp_mem);
    cp_mem->cp_rootMallocDone = FALSE;
  }

  /* Allocate necessary memory and return */
  if (!cp_mem->cp_rootMallocDone) {
    allocOK = cpRootAlloc(cp_mem, nrtfn);
    if (!allocOK) {
      cpProcessError(cp_mem, CP_MEM_FAIL, "CPODES", "CPodeRootInit", MSGCP_MEM_FAIL);
      return(CP_MEM_FAIL);
    }
    cp_mem->cp_rootMallocDone = TRUE;
  }

  /* Set variable values in CPODES memory block */
  cp_mem->cp_nrtfn  = nrtfn;
  cp_mem->cp_gfun   = gfun;
  cp_mem->cp_g_data = g_data;

  /* Set default values for rootdir (both directions)
   * and for gactive (all active) */
  for(i=0; i<nrtfn; i++) {
    cp_mem->cp_rootdir[i] = 0;
    cp_mem->cp_gactive[i] = TRUE;
  }

  /* Rootfinding is now enabled */
  cp_mem->cp_doRootfinding = TRUE;

  return(CP_SUCCESS);
}


/*
 * =================================================================
 * Readibility Constants
 * =================================================================
 */

#define uround         (cp_mem->cp_uround)
#define tn             (cp_mem->cp_tn)
#define h              (cp_mem->cp_h)
#define zn             (cp_mem->cp_zn)
#define y              (cp_mem->cp_y)
#define yp             (cp_mem->cp_yp)

#define lrw            (cp_mem->cp_lrw)
#define liw            (cp_mem->cp_liw)

#define taskc          (cp_mem->cp_taskc)
#define toutc          (cp_mem->cp_toutc)

#define nrtfn          (cp_mem->cp_nrtfn)
#define gfun           (cp_mem->cp_gfun)
#define g_data         (cp_mem->cp_g_data)
#define nge            (cp_mem->cp_nge)

#define gactive        (cp_mem->cp_gactive)
#define iroots         (cp_mem->cp_iroots)
#define rootdir        (cp_mem->cp_rootdir)
#define irfnd          (cp_mem->cp_irfnd)
#define tlo            (cp_mem->cp_tlo)
#define thi            (cp_mem->cp_thi)
#define glo            (cp_mem->cp_glo)
#define ghi            (cp_mem->cp_ghi)
#define trout          (cp_mem->cp_trout)
#define grout          (cp_mem->cp_grout)

/*
 * =================================================================
 * INTERNAL FUNCTIONS
 * =================================================================
 */

/*
 * -----------------------------------------------------------------
 * Root finding functions
 *   cpRcheck1    |
 *   cpRcheck2    |-> interface functions to CPode
 *   cpRcheck3    |
 *   cpRootfind
 * -----------------------------------------------------------------
 * Memory allocation/deallocation for rootfinding
 *   cpRootAlloc
 *   cpRootFree
 * -----------------------------------------------------------------
 */

/*
 * cpRcheck1
 *
 * This routine completes the initialization of rootfinding memory
 * information (it is called only once, at the very first step),
 * and checks whether g has any components that are zero BOTH at AND
 * very near the initial time of the IVP. Those components of g are
 * made inactive and will be later reactivated only when they move
 * away from zero.
 *
 * sherm 111125:
 *  thi (output only)    set to tn
 *  ghi (output only)    set to g(thi)
 * where we may have advanced time by smallh to see whether g's that were zero
 * at tn and deactivated can be reactivated at tn+smallh.
 *
 * The return value will be
 *    CV_RTFUNC_FAIL < 0 if the g function failed
 *    CP_SUCCESS     = 0 otherwise.
 */

int cpRcheck1(CPodeMem cp_mem)
{
  int i, retval;
  realtype ttol, smallh, hratio;
  booleantype zroot;

  for (i = 0; i < nrtfn; i++) iroots[i] = 0;

  thi = tlo = tn;
  ttol = (ABS(tn) + ABS(h))*uround*FUZZ_FACTOR;

  /*
   * Evaluate g at initial t and check for zero values.
   * Note that cpRcheck1 is called at the first step
   * before scaling zn[1] and therefore, y'(t0)=zn[1].
   */

  retval = gfun(tlo, zn[0], zn[1], glo, g_data);
  nge = 1;
  if (retval != 0) return(CP_RTFUNC_FAIL);

  /* Assume we won't find a root at the start. */
  for (i = 0; i < nrtfn; i++) ghi[i] = glo[i];

  zroot = FALSE;
  for (i = 0; i < nrtfn; i++) {
    if (ABS(glo[i]) == ZERO) {
      zroot = TRUE;
      gactive[i] = FALSE;
    }
  }
  if (!zroot) return(CP_SUCCESS);

  /*
   * Some g_i is zero at t0; look at g at t0+(small increment).
   * At the initial time and at order 1, we have:
   * y(t0+smallh) = zn[0] + (smallh/h) * zn[1]
   * y'(t0+smallh) = zn[1]
   */

  hratio = MAX(ttol/ABS(h), PT1);
  smallh = hratio*h;
  thi += smallh;
  N_VLinearSum(ONE, zn[0], hratio, zn[1], y);
  retval = gfun(thi, y, zn[1], ghi, g_data);
  nge++;
  if (retval != 0) return(CP_RTFUNC_FAIL);

  /*
   * We check now only the components of g which were exactly
   * zero at t0 to see if we can 'activate' them.
   */

  for (i = 0; i < nrtfn; i++) {
    if (!gactive[i] && ABS(ghi[i]) != ZERO) {
      gactive[i] = TRUE;
    }
  }

  return(CP_SUCCESS);
}


/*
 * cpRcheck2
 *
 * This routine is called at the beginning of a step to find the beginning tlo
 * of the next root search interval, which is usually the end (thi) of the
 * previous search interval. But it first checks for exact zeros of any active
 * g at thi. It then checks for a close pair of zeros (a condition that
 * would trigger making inactive the corresponding components
 * of g), and for a new root at a nearby point.
 * The endpoint thi of the previous search interval is thus adjusted
 * if necessary to assure that all active g_i are nonzero there,
 * before returning to do a root search in the interval.
 *
 * On entry, thi = tretlast is the last value of tret returned by
 * CPode.  This may be the previous tn, the previous tout value, or
 * the last root location.
 *
 * This routine returns an int equal to:
 *      CP_RTFUNC_FAIL < 0 if gfun failed
 *      RTFOUND        > 0 if a new zero of g was found near tlo, or
 *      CP_SUCCESS     = 0 otherwise.
 */

int cpRcheck2(CPodeMem cp_mem)
{
  int i, retval;
  realtype ttol, smallh, hratio;
  booleantype zroot;

  /* Move tlo up to end of previous search interval in case we find a root
  here. */
  tlo = thi;
  /* Evaluate g(tlo) */
  (void) cpGetSolution(cp_mem, tlo, y, yp);
  retval = gfun(tlo, y, yp, glo, g_data);
  nge++;
  if (retval != 0) return(CP_RTFUNC_FAIL);

  /* Assume we won't find a root at the start. */
  for (i = 0; i < nrtfn; i++) ghi[i] = glo[i];

  /* Check if any active g function is exactly ZERO at tlo.
   * If not, simply return CP_SUCCESS. */

  zroot = FALSE;
  for (i = 0; i < nrtfn; i++) iroots[i] = 0;
  for (i = 0; i < nrtfn; i++) {
    if (!gactive[i]) continue;
    if (ABS(glo[i]) == ZERO) {
      zroot = TRUE;
      iroots[i] = 1;
    }
  }

  /* If no root then tlo==thi, glo==ghi. */
  if (!zroot) return(CP_SUCCESS);

  /* One or more g_i has a zero at tlo.
   * Evaluate g(thi=tlo+smallh). */

  ttol = (ABS(tn) + ABS(h))*uround*FUZZ_FACTOR;
  smallh = (h > ZERO) ? ttol : -ttol;
  thi = tlo+smallh;
  if ( (thi - tn)*h >= ZERO) {
    hratio = smallh/h;
    N_VLinearSum(ONE, y, hratio, zn[1], y);
  } else {
    (void) cpGetSolution(cp_mem, thi, y, yp);
  }
  retval = gfun(thi, y, yp, ghi, g_data);
  nge++;
  if (retval != 0) return(CP_RTFUNC_FAIL);

  /* Check if any active function is ZERO at thi+smallh.
   * Make inactive those that were also ZERO at thi.
   * Report a root for those that only became ZERO at tlo+smallh. */

  zroot = FALSE;
  for (i = 0; i < nrtfn; i++) {
    if (ABS(ghi[i]) == ZERO) {
      if (!gactive[i]) continue;
      if (iroots[i] == 1) { iroots[i] = 0; gactive[i] = FALSE; }
      else                { iroots[i] = 1; zroot = TRUE; }
    }
  }

  if (zroot) return(RTFOUND);

  return(CP_SUCCESS);
}

/*
 * cpRcheck3
 *
 * This routine interfaces to cpRootfind to look for a root of g
 * between tlo and either tn or tout, whichever comes first.
 * Only roots beyond tlo in the direction of integration are sought.
 *
 * On entry, both thi and ghi=g(thi) should have been evaluated. We start by
 * setting tlo=thi and glo=ghi, shiting the search interval to start at the
 * end of the previous one.
 * On return, if there is a root it is in (tlo,thi] which will have been
 * adjusted to a very narrow bracket around the zero crossing. If there is no
 * root then thi and ghi are at the end of the search interval, where they can
 * serve as the start for the next one.
 *
 * This routine returns an int equal to:
 *      CP_RTFUNC_FAIL < 0 if the g function failed,
 *      RTFOUND        > 0 if a root of g was found, or
 *      CP_SUCCESS     = 0 otherwise.
 */

int cpRcheck3(CPodeMem cp_mem)
{
  int i, retval, ier;
  realtype ttol;

  /* Move start of search interval to end of previous one. */
  tlo = thi;
  for (i = 0; i < nrtfn; ++i) glo[i] = ghi[i];

  /* Set thi = tn or tout, whichever comes first. */
  switch (taskc) {
  case CP_ONE_STEP:
    thi = tn;
    break;
  case CP_NORMAL:
    thi = ( (toutc - tn)*h >= ZERO) ? tn : toutc;
    break;
  }

  /* Get y and y' at thi. */
  (void) cpGetSolution(cp_mem, thi, y, yp);

  /* Set ghi = g(thi) */
  retval = gfun(thi, y, yp, ghi, g_data);
  nge++;
  if (retval != 0) return(CP_RTFUNC_FAIL);

  /* Call cpRootfind to search (tlo,thi) for roots, and to modify tlo,thi
  to create a very narrow bracket around the first root. */
  ttol = (ABS(tn) + ABS(h))*uround*FUZZ_FACTOR;
  ier = cpRootfind(cp_mem, ttol);

  /* If the root function g failed, return now */
  if (ier == CP_RTFUNC_FAIL) return(CP_RTFUNC_FAIL);

  /* If any of the inactive components moved away from zero,
   * activate them now. */
  for(i=0; i<nrtfn; i++) {
    if(!gactive[i] && grout[i] != ZERO) gactive[i] = TRUE;
  }

  /* If no root found, return CP_SUCCESS. */
  if (ier == CP_SUCCESS) return(CP_SUCCESS);

  /* If a root was found, interpolate to get y(trout) and return.  */
  (void) cpGetSolution(cp_mem, trout, y, yp);

  return(RTFOUND);
}

/*
 * cpRootFind
 *
 * This routine solves for a root of g(t) between tlo and thi, if
 * one exists.  Only roots of odd multiplicity (i.e. with a change
 * of sign in one of the g_i), or exact zeros, are found.
 * Here the sign of tlo - thi is arbitrary, but if multiple roots
 * are found, the one closest to tlo is returned.
 *
 * The method used is the Illinois algorithm, a modified secant method.
 * Reference: Kathie L. Hiebert and Lawrence F. Shampine, Implicitly
 * Defined Output Points for Solutions of ODEs, Sandia National
 * Laboratory Report SAND80-0180, February 1980.
 *
 * This routine uses the following parameters for communication:
 *
 * nrtfn    = number of functions g_i, or number of components of
 *            the vector-valued function g(t).  Input only.
 *
 * gfun     = user-defined function for g(t).  Its form is
 *            (void) gfun(t, y, gt, g_data)
 *
 * gactive  = array specifying whether a component of g should
 *            or should not be monitored. gactive[i] is initially
 *            set to TRUE for all i=0,...,nrtfn-1, but it may be
 *            reset to FALSE if at the first step g[i] is 0.0
 *            both at the I.C. and at a small perturbation of them.
 *            gactive[i] is then set back on TRUE only after the
 *            corresponding g function moves away from 0.0.
 *
 * rootdir  = array specifying the direction of zero-crossings.
 *            If rootdir[i] > 0, search for roots of g_i only if
 *            g_i is increasing; if rootdir[i] < 0, search for
 *            roots of g_i only if g_i is decreasing; otherwise
 *            always search for roots of g_i.
 *
 * nge      = cumulative counter for gfun calls.
 *
 * ttol     = a convergence tolerance for trout.  Input only.
 *            When a root at trout is found, it is located in time only to
 *            within a tolerance of ttol.  Typically, ttol should
 *            be set to a value on the order of
 *               100 * UROUND * max (ABS(tlo), ABS(thi))
 *            where UROUND is the unit roundoff of the machine.
 *
 * tlo, thi = endpoints of the interval in which roots are sought.
 *            On input, they must be distinct, but tlo - thi may
 *            be of either sign.  The direction of integration is
 *            assumed to be from tlo to thi.  On return, tlo and thi
 *            are the endpoints of the final relevant interval (tlo,thi];
 *            that is, the root has not yet occurred at tlo but has
 *            definitely occurred by thi. The reported root time trout is
 *            always the same as thi.
 *
 * glo, ghi = arrays of length nrtfn containing the vectors g(tlo)
 *            and g(thi) respectively.  Input and output.  On input,
 *            none of the active glo[i] should be zero.
 *
 * trout    = root location (same as thi), if a root was found, or the original
 *            value of thi if not. Output only. trout is the endpoint thi of
 *            the final interval (tlo,thi] bracketing the root, with |thi-tlo|
 *            at most ttol.
 *
 * grout    = array of length nrtfn containing g(trout) (==ghi) on return.
 *
 * iroots   = int array of length nrtfn with root information.
 *            Output only.  If a root was found, iroots indicates
 *            which components g_i have a root at trout.  For
 *            i = 0, ..., nrtfn-1, iroots[i] = 1 if g_i has a root
 *            and g_i is increasing, iroots[i] = -1 if g_i has a
 *            root and g_i is decreasing, and iroots[i] = 0 if g_i
 *            has no roots or g_i varies in the direction opposite
 *            to that indicated by rootdir[i].
 *
 * This routine returns an int equal to:
 *      CP_RTFUNC_FAIL < 0 if gfun faild
 *      RTFOUND        > 0 if a root of g was found, or
 *      CP_SUCCESS     = 0 otherwise.
 */
static int cpRootfind(CPodeMem cp_mem, realtype ttol)
{
  realtype alpha, tmid, my_tmid, tmid_saved, thi_saved, fracint, fracsub;
  int i, retval, side, sideprev;
  booleantype zroot;

  /* alpha is a bias weight in the secant method.
   * On the first two passes, set alpha = 1.  Thereafter, reset alpha
   * according to the side (low vs high) of the subinterval in which
   * the sign change was found in the previous two passes.
   * If the sides were opposite, set alpha = 1.
   * If the sides were the same, then double alpha (if high side),
   * or halve alpha (if low side).
   * The next guess tmid is the secant method value if alpha = 1, but
   * is closer to tlo if alpha < 1, and closer to thi if alpha > 1.
   */
  alpha = ONE;

  /* First, for each active g function, check whether an event occurred in
   * (tlo,thi). Since glo != 0 for an active component, this means we check for
   * a sign change or for ghi = 0 (taking into account rootdir). For each
   * component that triggers an event, we estimate a "proposal" mid point (by
   * bisection if ghi=0 or with secant method otherwise) and select the one
   * closest to tlo. */
  zroot = FALSE;
  tmid = thi;
  for (i = 0;  i < nrtfn; i++) {
    if(!gactive[i]) continue;

    if ( (glo[i]*ghi[i] <= ZERO) && (rootdir[i]*glo[i] <= ZERO) ) {
      zroot = TRUE;
      if (ghi[i] == ZERO) {
        my_tmid = thi - HALF * (thi-tlo);
      } else {
        my_tmid = thi - (thi - tlo)*ghi[i]/(ghi[i] - alpha*glo[i]);
      }
      /* Pick my_tmid if it is closer to tlo than the current tmid. */
      if ( (my_tmid-tmid)*h < ZERO ) tmid = my_tmid;
    }
  }

  /* If no event was detected, set trout to thi and return CP_SUCCESS */
  if (!zroot) {
    trout = thi;
    for (i = 0; i < nrtfn; i++) grout[i] = ghi[i];
    return(CP_SUCCESS);
  }

  /* An event was detected. Loop to locate nearest root. */
  side = 0;

  loop {
    sideprev = side;

    /* If tmid is too close to tlo or thi, adjust it inward,
     * by a fractional distance that is between 0.1 and 0.5. */
    if (ABS(tmid - tlo) < HALF*ttol) {
      fracint = ABS(thi - tlo)/ttol;
      fracsub = (fracint > FIVE) ? PT1 : HALF/fracint;
      tmid = tlo + fracsub*(thi - tlo);
    }
    if (ABS(thi - tmid) < HALF*ttol) {
      fracint = ABS(thi - tlo)/ttol;
      fracsub = (fracint > FIVE) ? PT1 : HALF/fracint;
      tmid = thi - fracsub*(thi - tlo);
    }

    /* Get solution at tmid and evaluate g(tmid). */
    (void) cpGetSolution(cp_mem, tmid, y, yp);
    retval = gfun(tmid, y, yp, grout, g_data);
    nge++;
    if (retval != 0) return(CP_RTFUNC_FAIL);

    /* Check (tlo, tmid) to see if an event occurred in the "low" side
     * First make temporary copies of thi and tmid in case the event
     * turns out to be on the "high" side. */
    tmid_saved = tmid;
    thi_saved  = thi;
    thi = tmid;
    zroot = FALSE;
    for (i = 0;  i < nrtfn; i++) {
      if(!gactive[i]) continue;

      if ( (glo[i]*grout[i] <= ZERO) && (rootdir[i]*glo[i] <= ZERO) ) {
        zroot = TRUE;
        if (grout[i] == ZERO) {
          my_tmid = thi - HALF * (thi-tlo);
        } else {
          my_tmid = thi - (thi - tlo)*grout[i]/(grout[i] - alpha*glo[i]);
        }
        if ( (my_tmid-tmid)*h < ZERO ) tmid = my_tmid;
      }
    }

    /* If we detected an event in the "low" side:
       - accept current value of thi
       - set ghi <- grout
       - test for convergence and break from loop if converged
       - set side=1 (low); if previous side was also 1, scale alpha by 1/2
       - continue looping to refine root location */
    if (zroot) {
      for (i = 0; i < nrtfn; i++) ghi[i] = grout[i];
      if (ABS(thi - tlo) <= ttol) break;

      side = 1;
      if (sideprev == 1) alpha = alpha*HALF;
      else               alpha = ONE;

      continue;
    }

    /* No event detected in "low" side; event must be in "high" side.
       - restore previously saved values for thi and set tlo = tmid_saved
       - set glo <- grout
       - test for convergence and break from loop if converged
       - set side=2 (high); if previous side was also 2, scale alpha by 2
       - continue looping to refine root location */
    thi = thi_saved;
    tlo = tmid_saved;

    for (i = 0; i < nrtfn; i++) glo[i] = grout[i];
    if (ABS(thi - tlo) <= ttol) break;

    tmid = thi;
    for (i = 0;  i < nrtfn; i++) {
      if(!gactive[i]) continue;

      if ( (glo[i]*ghi[i] <= ZERO) && (rootdir[i]*glo[i] <= ZERO) ) {
        if (ghi[i] == ZERO) {
          my_tmid = thi - HALF * (thi-tlo);
        } else {
          my_tmid = thi - (thi - tlo)*ghi[i]/(ghi[i] - alpha*glo[i]);
        }
        if ( (my_tmid-tmid)*h < ZERO ) tmid = my_tmid;
      }
    }

    side = 2;
    if (sideprev == 2) alpha = alpha*TWO;
    else               alpha = ONE;

  } /* End of root-search loop */

  /* Root has been isolated to (tlo,thi] and |thi-tlo| <= ttol. We'll declare
  that the root was found at trout=thi.
     - Reset trout and grout
     - Set iroots
     - Return RTFOUND. */
  trout = thi;
  for (i = 0; i < nrtfn; i++) {
    grout[i] = ghi[i];
    iroots[i] = 0;
    if(!gactive[i]) continue;
    if ( (glo[i]*ghi[i] <= ZERO) && (rootdir[i]*glo[i] <= ZERO) )
      iroots[i] = glo[i] > 0 ? -1:1;
  }

  return(RTFOUND);
}



/*
 * cpRootAlloc allocates memory for the rootfinding algorithm.
 */

static booleantype cpRootAlloc(CPodeMem cp_mem, int nrt)
{
  glo = NULL;
  glo = (realtype *) malloc(nrt*sizeof(realtype));

  ghi = NULL;
  ghi = (realtype *) malloc(nrt*sizeof(realtype));
  if (ghi == NULL) {
    free(glo); glo = NULL;
    return(FALSE);
  }

  grout = NULL;
  grout = (realtype *) malloc(nrt*sizeof(realtype));
  if (grout == NULL) {
    free(glo); glo = NULL;
    free(ghi); ghi = NULL;
    return(FALSE);
  }

  iroots = NULL;
  iroots = (int *) malloc(nrt*sizeof(int));
  if (iroots == NULL) {
    free(glo); glo = NULL;
    free(ghi); ghi = NULL;
    free(grout); grout = NULL;
    return(FALSE);
  }

  rootdir = NULL;
  rootdir = (int *) malloc(nrt*sizeof(int));
  if (rootdir == NULL) {
    free(glo); glo = NULL;
    free(ghi); ghi = NULL;
    free(grout); grout = NULL;
    free(iroots); iroots = NULL;
  }

  gactive = NULL;
  gactive = (booleantype *) malloc(nrt*sizeof(booleantype));
  if (gactive == NULL) {
    free(glo); glo = NULL;
    free(ghi); ghi = NULL;
    free(grout); grout = NULL;
    free(iroots); iroots = NULL;
    free(rootdir); rootdir = NULL;
  }

  lrw += 3*nrt;
  liw += 3*nrt;

  return(TRUE);
}


/*
 * cpRootFree frees the memory allocated in cpRootAlloc.
 */

void cpRootFree(CPodeMem cp_mem)
{
  free(glo); glo = NULL;
  free(ghi); ghi = NULL;
  free(grout); grout = NULL;
  free(gactive); gactive = NULL;
  free(iroots); iroots = NULL;
  free(rootdir); rootdir = NULL;

  lrw -= 3 * nrtfn;
  liw -= 3 * nrtfn;
}

