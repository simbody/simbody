// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpOptErrorConvCheck.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPOPTERRORCONVCHECK_HPP__
#define __IPOPTERRORCONVCHECK_HPP__

#include "IpConvCheck.hpp"

namespace SimTKIpopt
{

  /** Brief Class Description.
   *  Detailed Class Description.
   */
  class OptimalityErrorConvergenceCheck : public ConvergenceCheck
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    OptimalityErrorConvergenceCheck();

    /** Default destructor */
    virtual ~OptimalityErrorConvergenceCheck();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) override;

    /** Overloaded convergence check */
    virtual ConvergenceStatus
    CheckConvergence(bool call_intermediate_callback = true) override;

    /** Auxilliary function for testing whether current iterate
     *  satisfies the acceptable level of optimality */
    virtual bool CurrentIsAcceptable() override;

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

  protected:
    /** @name Algorithmic parameters */
    //@{
    /** Maximal number of iterations */
    Index max_iterations_;
    /** Tolerance on unscaled dual infeasibility */
    Number dual_inf_tol_;
    /** Tolerance on unscaled constraint violation */
    Number constr_viol_tol_;
    /** Tolerance on unscaled complementarity */
    Number compl_inf_tol_;
    /** Number of iterations with acceptable level of accuracy, after
     *  which the algorithm terminates.  If 0, this heuristic is
     *  disabled. */
    Index acceptable_iter_;
    /** Acceptable tolerance for the problem to terminate earlier if
     *  algorithm seems stuck or cycling */
    Number acceptable_tol_;
    /** Acceptable tolerance on unscaled dual infeasibility */
    Number acceptable_dual_inf_tol_;
    /** Acceptable tolerance on unscaled constraint violation */
    Number acceptable_constr_viol_tol_;
    /** Acceptable tolerance on unscaled complementarity */
    Number acceptable_compl_inf_tol_;
    /** Threshold for primal iterates for divergence test */
    Number diverging_iterates_tol_;
    //@}

  private:
    /**@name Default Compiler Generated Methods (Hidden to avoid
     * implicit creation/calling).  These methods are not implemented
     * and we do not want the compiler to implement them for us, so we
     * declare them private and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    OptimalityErrorConvergenceCheck(const OptimalityErrorConvergenceCheck&);

    /** Overloaded Equals Operator */
    void operator=(const OptimalityErrorConvergenceCheck&);
    //@}

    /** Counter for successive iterations in which acceptability
     *  criteria are met. */
    Index acceptable_counter_;
  };

} // namespace Ipopt

#endif
