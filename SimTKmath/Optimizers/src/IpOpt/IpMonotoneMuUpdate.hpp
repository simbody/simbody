// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpMonotoneMuUpdate.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPMONOTONEMUUPDATE_HPP__
#define __IPMONOTONEMUUPDATE_HPP__

#include "IpMuUpdate.hpp"
#include "IpLineSearch.hpp"
#include "IpRegOptions.hpp"

namespace SimTKIpopt
{

  /** Monotone Mu Update. This class implements the standard monotone mu update
   *  approach.
   */
  class MonotoneMuUpdate : public MuUpdate
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    MonotoneMuUpdate(const SmartPtr<LineSearch>& linesearch);

    /** Default destructor */
    virtual ~MonotoneMuUpdate();
    //@}

    /** Initialize method - overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) override;

    /** Method for determining the barrier parameter for the next
     *  iteration.  When the optimality error for the current barrier
     *  parameter is less than a tolerance, the barrier parameter is
     *  reduced, and the Reset method of the LineSearch object
     *  linesearch is called. */
    virtual bool UpdateBarrierParameter() override;

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(const SmartPtr<RegisteredOptions>& roptions);
    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    MonotoneMuUpdate();

    /** Copy Constructor */
    MonotoneMuUpdate(const MonotoneMuUpdate&);

    /** Overloaded Equals Operator */
    void operator=(const MonotoneMuUpdate&);
    //@}

    /** Internal method for computing the new values for mu and tau */
    void CalcNewMuAndTau(Number &new_mu,
                         Number &new_tau);

    /** @name Algorithmic parameters */
    //@{
    /** Initial value of the barrier parameter */
    Number mu_init_;
    Number barrier_tol_factor_;
    Number mu_linear_decrease_factor_;
    Number mu_superlinear_decrease_power_;
    bool mu_allow_fast_monotone_decrease_;
    /** Tau_min for fraction to boundary rule */
    Number tau_min_;
    Number compl_inf_tol_;
    //@}

    SmartPtr<LineSearch> linesearch_;

    /** Flag indicating whether the method has been called at least once so
     *  far */
    bool initialized_;

    /** If true, no modification of the barrier parameter will be done
     *  at the first call of Update (fix for the restoration phase -
     *  we should clean that up!) */
    bool first_iter_resto_;
  };

} // namespace Ipopt

#endif
