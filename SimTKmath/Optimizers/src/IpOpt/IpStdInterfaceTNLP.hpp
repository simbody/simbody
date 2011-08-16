// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpStdInterfaceTNLP.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPSTDINTERFACETNLP_HPP__
#define __IPSTDINTERFACETNLP_HPP__

#include "IpUtils.hpp"
#include "IpTNLP.hpp"
#include "IpJournalist.hpp"
#include "IpException.hpp"
#include "IpStdCInterface.h"
#include "IpSmartPtr.hpp"

namespace Ipopt
{
  /** Declare excpetion that is thrown when invalid NLP data
   *  is provided */
  DECLARE_STD_EXCEPTION(INVALID_STDINTERFACE_NLP);

  /** Implementation of a TNLP for the Standard C interface.  The
   *  standard C interface is exposed to the user as a single C
   *  function that is given problem dimension, starting points, and
   *  pointers for functions that evaluate objective function etc.
   */
  class StdInterfaceTNLP : public TNLP
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor, given dimensions of problem, function pointers
     *  for evaluation callback functions, and starting points. Note
     *  that the constrctor does not make a copy of any of the Number
     *  arrays, i.e. it is up to the called to keep them around. */
    StdInterfaceTNLP(Index n_var,
                     const Number* x_L, const Number* x_U,
                     Index n_con,
                     const Number* g_L, const Number* g_U,
                     Index nele_jac,
                     Index nele_hess,
                     Index index_style,
                     const Number* start_x,
                     const Number* start_lam,
                     const Number* start_z_L,
                     const Number* start_z_U,
                     Eval_F_CB eval_f,
                     Eval_G_CB eval_g,
                     Eval_Grad_F_CB eval_grad_f,
                     Eval_Jac_G_CB eval_jac_g,
                     Eval_H_CB eval_h,
                     Number* x_sol,
                     Number* z_L_sol,
                     Number* z_U_sol,
                     Number* g_sol,
                     Number* lam_sol,
                     Number* obj_sol,
                     UserDataPtr user_data);

    /** Default destructor */
    virtual ~StdInterfaceTNLP();
    //@}

    /**@name methods to gather information about the NLP. These methods are
     * overloaded from TNLP. See TNLP for their more detailed documentation. */
    //@{
    /** returns dimensions of the nlp. Overloaded from TNLP */
    virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                              Index& nnz_h_lag, IndexStyleEnum& index_style);

    /** returns bounds of the nlp. Overloaded from TNLP */
    virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                                 Index m, Number* g_l, Number* g_u);

    /** provides a starting point for the nlp variables. Overloaded from TNLP */
    virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                    bool init_z, Number* z_L, Number* z_U,
                                    Index m, bool init_lambda, Number* lambda);

    /** evaluates the objective value for the nlp. Overloaded from TNLP */
    virtual bool eval_f(Index n, const Number* x, bool new_x,
                        Number& obj_value);

    /** evaluates the gradient of the objective for the
     *  nlp. Overloaded from TNLP */
    virtual bool eval_grad_f(Index n, const Number* x, bool new_x,
                             Number* grad_f);

    /** evaluates the constraint residuals for the nlp. Overloaded from TNLP */
    virtual bool eval_g(Index n, const Number* x, bool new_x, Index m,
                        Number* g);

    /** specifies the jacobian structure (if values is NULL) and
     *  evaluates the jacobian values (if values is not NULL) for the
     *  nlp. Overloaded from TNLP */
    virtual bool eval_jac_g(Index n, const Number* x, bool new_x, Index m,
                            Index nele_jac, Index* iRow, Index *jCol,
                            Number* values);

    /** specifies the structure of the hessian of the lagrangian (if values is NULL) and
     *  evaluates the values (if values is not NULL). Overloaded from TNLP */
    virtual bool eval_h(Index n, const Number* x, bool new_x,
                        Number obj_factor, Index m, const Number* lambda,
                        bool new_lambda, Index nele_hess, Index* iRow,
                        Index* jCol, Number* values);
    //@}

    /** @name Solution Methods */
    //@{
    virtual void finalize_solution(SolverReturn status,
                                   Index n, const Number* x, const Number* z_L, const Number* z_U,
                                   Index m, const Number* g, const Number* lambda,
                                   Number obj_value);
    //@}

  private:
    /** Journlist */
    SmartPtr<const Journalist> jnlst_;

    /** @name Information about the problem */
    //@{
    /** Number of variables */
    const Index n_var_;
    /** Number of constraints */
    const Index n_con_;
    /** Pointer to Number array containing lower bounds for variables */
    const Number* x_L_;
    /** Pointer to Number array containing upper bounds for variables */
    const Number* x_U_;
    /** Pointer to Number array containing lower bounds for constraints */
    const Number* g_L_;
    /** Pointer to Number array containing upper bounds for constraints */
    const Number* g_U_;
    /** Number of non-zero elements in the constraint Jacobian */
    const Index nele_jac_;
    /** Number of non-zero elements in the Hessian */
    const Index nele_hess_;
    /** Starting value of the iRow and jCol parameters for matrices */
    const Index index_style_;
    /** Pointer to Number array containing starting point for variables */
    const Number* start_x_;
    /** Poitner to Number array containing starting values for
     *  constraint multipliers */
    const Number* start_lam_;
    /** Pointer to Number array containing starting values for lower
     *  bound multipliers */
    const Number* start_z_L_;
    /** Pointer to Number array containing starting values for upper
     *  bound multipliers */
    const Number* start_z_U_;
    /** Pointer to callback function evaluating value of objective function */
    Eval_F_CB eval_f_;
    /**  Pointer to callback function evaluating value of constraints */
    Eval_G_CB eval_g_;
    /** Pointer to callback function evaluating gradient of objective
     *  function */
    Eval_Grad_F_CB eval_grad_f_;
    /** Pointer to callback function evaluating Jacobian of constraints */
    Eval_Jac_G_CB eval_jac_g_;
    /** Pointer to callback function evaluating Hessian of Lagrangian */
    Eval_H_CB eval_h_;
    /** Pointer to user data */
    UserDataPtr user_data_;
    //@}


    /** A non-const copy of x - this is kept up-to-date in apply_new_x */
    Number* non_const_x_;

    /** Pointers to the user provided vectors for solution */
    Number* x_sol_;
    Number* z_L_sol_;
    Number* z_U_sol_;
    Number* g_sol_;
    Number* lambda_sol_;
    Number* obj_sol_;

    /** Internal function to update the internal and ampl state if the
     *  x value changes */
    void apply_new_x(bool new_x, Index n, const Number* x);


    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    StdInterfaceTNLP();

    /** Copy Constructor */
    StdInterfaceTNLP(const StdInterfaceTNLP&);

    /** Overloaded Equals Operator */
    void operator=(const StdInterfaceTNLP&);
    //@}

  };

} // namespace Ipopt

#endif
