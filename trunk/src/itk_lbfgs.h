#ifndef lbfgs_h_
#define lbfgs_h_

struct lb3_1_ {
/*
C    GTOL is a DOUBLE PRECISION variable with default value 0.9, which
C        controls the accuracy of the line search routine MCSRCH. If the
C        function and gradient evaluations are inexpensive with respect
C        to the cost of the iteration (which is sometimes the case when
C        solving very large problems) it may be advantageous to set GTOL
C        to a small value. A typical small value is 0.1.  Restriction:
C        GTOL should be greater than 1.D-04.
C
C    STPMIN and STPMAX are non-negative DOUBLE PRECISION variables which
C        specify lower and upper bounds for the step in the line search.
C        Their default values are 1.D-20 and 1.D+20, respectively. These
C        values need not be modified unless the exponents are too large
C        for the machine being used, or unless the problem is extremely
C        badly scaled (in which case the exponents should be increased).
*/
  int mp, lp; /* Fortran i/o stuff.  Unused here. */
  double gtol, stpmin, stpmax;
  double stpawf; /* line search default step length, added by awf */
};

/*#define lb3_1 (*(struct lb3_1_ *) &lb3_)*/
#define lb3_1 lb3_

extern struct lb3_1_ lb3_1;

#endif /* lbfgs_h_ */
