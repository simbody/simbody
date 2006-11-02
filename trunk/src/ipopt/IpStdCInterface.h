/*************************************************************************
   Copyright (C) 2004, 2006 International Business Machines and others.
   All Rights Reserved.
   This code is published under the Common Public License.
 
   $Id: IpStdCInterface.h 759 2006-07-07 03:07:08Z andreasw $
 
   Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-02
 *************************************************************************/

#ifndef __IPSTDCINTERFACE_H__
#define __IPSTDCINTERFACE_H__

#ifdef __cplusplus
extern "C"
{
#endif

  /** Type for all number.  We need to make sure that this is
      identical with what is defined in Common/IpTypes.hpp */
  typedef double Number;

  /** Type for all incides.  We need to make sure that this is
      identical with what is defined in Common/IpTypes.hpp */
  typedef int Index;

  /** Type for all integers.  We need to make sure that this is
      identical with what is defined in Common/IpTypes.hpp */
  typedef int Int;

  /* This includes the SolverReturn enum type */
#include "IpReturnCodes.h"

  /** Structure collecting all information about the problem
   *  definition and solve statistics etc.  This is defined in the
   *  source file. */
  struct IpoptProblemInfo;

  /** Pointer to a Ipopt Problem. */
  typedef struct IpoptProblemInfo* IpoptProblem;

  /** define a boolean type for C */
  typedef int Bool;
#ifndef TRUE
# define TRUE (1)
#endif
#ifndef FALSE
# define FALSE (0)
#endif

  /** A pointer for anything that is to be passed between the called
   *  and individual callback function */
  typedef void * UserDataPtr;

  /** Type defining the callback function for evaluating the value of
   *  the objective function.  Return value should be set to false if
   *  there was a problem doing the evaluation. */
  typedef Bool (*Eval_F_CB)(Index n, Number* x, Bool new_x,
                            Number* obj_value, UserDataPtr user_data);

  /** Type defining the callback function for evaluating the gradient of
   *  the objective function.  Return value should be set to false if
   *  there was a problem doing the evaluation. */
  typedef Bool (*Eval_Grad_F_CB)(Index n, Number* x, Bool new_x,
                                 Number* grad_f, UserDataPtr user_data);

  /** Type defining the callback function for evaluating the value of
   *  the constraint functions.  Return value should be set to false if
   *  there was a problem doing the evaluation. */
  typedef Bool (*Eval_G_CB)(Index n, Number* x, Bool new_x,
                            Index m, Number* g, UserDataPtr user_data);

  /** Type defining the callback function for evaluating the Jacobian of
   *  the constrant functions.  Return value should be set to false if
   *  there was a problem doing the evaluation. */
  typedef Bool (*Eval_Jac_G_CB)(Index n, Number *x, Bool new_x,
                                Index m, Index nele_jac,
                                Index *iRow, Index *jCol, Number *values,
                                UserDataPtr user_data);

  /** Type defining the callback function for evaluating the Hessian of
   *  the Lagrangian function.  Return value should be set to false if
   *  there was a problem doing the evaluation. */
  typedef Bool (*Eval_H_CB)(Index n, Number *x, Bool new_x, Number obj_factor,
                            Index m, Number *lambda, Bool new_lambda,
                            Index nele_hess, Index *iRow, Index *jCol,
                            Number *values, UserDataPtr user_data);

  /** Function for creating a new Ipopt Problem object.  This function
   *  returns an object that can be passed to the IpoptSolve call.  It
   *  contains the basic definition of the optimization problem, such
   *  as number of variables and constraints, bounds on variables and
   *  constraints, information about the derivatives, and the callback
   *  function for the computation of the optimization problem
   *  functions and derivatives.  During this call, the options file
   *  PARAMS.DAT is read as well.
   *
   *  If NULL is returned, there was a problem with one of the inputs
   *  or reading the options file. */
  IpoptProblem CreateIpoptProblem(
      Index n             /** Number of optimization variables */
    , Number* x_L         /** Lower bounds on variables. This array of
                              size n is copied internally, so that the
                              caller can change the incoming data after
                              return without that IpoptProblem is
                              modified.  Any value less or equal than
                              the number specified by option
                              'nlp_lower_bound_inf' is interpreted to
                              be minus infinity. */
    , Number* x_U         /** Upper bounds on variables. This array of
                              size n is copied internally, so that the
                              caller can change the incoming data after
                              return without that IpoptProblem is
                              modified.  Any value greater or equal
                              than the number specified by option
                              'nlp_upper_bound_inf' is interpreted to
                              be plus infinity. */
    , Index m             /** Number of constraints. */
    , Number* g_L         /** Lower bounds on constraints. This array of
                              size m is copied internally, so that the
                              caller can change the incoming data after
                              return without that IpoptProblem is
                              modified.  Any value less or equal than
                              the number specified by option
                              'nlp_lower_bound_inf' is interpreted to
                              be minus infinity. */
    , Number* g_U         /** Upper bounds on constraints. This array of
                              size m is copied internally, so that the
                              caller can change the incoming data after
                              return without that IpoptProblem is
                              modified.  Any value greater or equal
                              than the number specified by option
                              'nlp_upper_bound_inf' is interpreted to
                              be plus infinity. */
    , Index nele_jac      /** Number of non-zero elements in constraint
                              Jacobian. */
    , Index nele_hess     /** Number of non-zero elements in Hessian of
                              Lagrangian. */
    , Index index_style   /** indexing style for iRow & jCol,
				 0 for C style, 1 for Fortran style */
    , Eval_F_CB eval_f    /** Callback function for evaluating
                              objective function */
    , Eval_G_CB eval_g    /** Callback function for evaluating
                              constraint functions */
    , Eval_Grad_F_CB eval_grad_f
                          /** Callback function for evaluating gradient
                              of objective function */
    , Eval_Jac_G_CB eval_jac_g
                          /** Callback function for evaluating Jacobian
                              of constraint functions */
    , Eval_H_CB eval_h    /** Callback function for evaluating Hessian
                              of Lagrangian function */
  );

  /** Method for freeing a previously created IpoptProblem.  After
      freeing an IpoptProblem, it cannot be used anymore. */
  void FreeIpoptProblem(IpoptProblem ipopt_problem);


  /** Function for adding a string option.  Returns FALSE the option
   *  could not be set (e.g., if keyword is unknown) */
  Bool AddIpoptStrOption(IpoptProblem ipopt_problem, char* keyword, char* val);

  /** Function for adding a Number option.  Returns FALSE the option
   *  could not be set (e.g., if keyword is unknown) */
  Bool AddIpoptNumOption(IpoptProblem ipopt_problem, char* keyword, Number val);

  /** Function for adding an Int option.  Returns FALSE the option
   *  could not be set (e.g., if keyword is unknown) */
  Bool AddIpoptIntOption(IpoptProblem ipopt_problem, char* keyword, Int val);

  /** Function for opening an output file for a given name with given
   *  printlevel.  Returns false, if there was a problem opening the
   *  file. */
  Bool OpenIpoptOutputFile(IpoptProblem ipopt_problem, char* file_name,
                           Int print_level);

  /** Function calling the Ipopt optimization algorithm for a problem
      previously defined with CreateIpoptProblem.  The return
      specified outcome of the optimization procedure (e.g., success,
      failure etc).
   */
  enum ApplicationReturnStatus IpoptSolve(
      IpoptProblem ipopt_problem
                         /** Problem that is to be optimized.  Ipopt
                             will use the options previously specified with
                             AddIpoptOption (etc) for this problem. */
    , Number* x          /** Input:  Starting point
                             Output: Optimal solution */
    , Number* g          /** Values of constraint at final point
                             (output only - ignored if set to NULL) */
    , Number* obj_val    /** Final value of objective function
                             (output only - ignored if set to NULL) */
    , Number* mult_g     /** Final multipliers for constraints
                             (output only - ignored if set to NULL) */
    , Number* mult_x_L   /** Final multipliers for lower variable bounds
                             (output only - ignored if set to NULL) */
    , Number* mult_x_U   /** Final multipliers for upper variable bounds
                             (output only - ignored if set to NULL) */
    , UserDataPtr user_data
    /** Pointer to user data.  This will be
    passed unmodified to the callback
    functions. */
  );

  /**
  void IpoptStatisticsCounts;

  void IpoptStatisticsInfeasibilities; */

#ifdef __cplusplus
} /* extern "C" { */
#endif

#endif
