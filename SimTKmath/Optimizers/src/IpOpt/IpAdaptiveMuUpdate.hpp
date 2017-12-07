// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpAdaptiveMuUpdate.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPADAPTIVEMUUPDATE_HPP__
#define __IPADAPTIVEMUUPDATE_HPP__

#include "IpMuUpdate.hpp"
#include "IpLineSearch.hpp"
#include "IpMuOracle.hpp"
#include "IpFilter.hpp"
#include "IpQualityFunctionMuOracle.hpp"

namespace SimTKIpopt
{

  /** Non-monotone mu update.
   */
  class AdaptiveMuUpdate : public MuUpdate
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor */
    AdaptiveMuUpdate(const SmartPtr<LineSearch>& linesearch,
                     const SmartPtr<MuOracle>& free_mu_oracle,
                     const SmartPtr<MuOracle>& fix_mu_oracle=NULL);
    /** Default destructor */
    virtual ~AdaptiveMuUpdate();
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
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
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

    /** Default Constructor */
    AdaptiveMuUpdate();

    /** Copy Constructor */
    AdaptiveMuUpdate(const AdaptiveMuUpdate&);

    /** Overloaded Equals Operator */
    void operator=(const AdaptiveMuUpdate&);
    //@}

    /** @name Algorithmic parameters */
    //@{
    Number mu_max_fact_;
    Number mu_max_;
    Number mu_min_;
    bool mu_min_default_;
    Number tau_min_;
    Number adaptive_mu_safeguard_factor_; //ToDo don't need that?
    Number adaptive_mu_monotone_init_factor_;
    Number barrier_tol_factor_;
    Number mu_linear_decrease_factor_;
    Number mu_superlinear_decrease_power_;
    QualityFunctionMuOracle::NormEnum adaptive_mu_kkt_norm_;
    QualityFunctionMuOracle::CentralityEnum adaptive_mu_kkt_centrality_;
    QualityFunctionMuOracle::BalancingTermEnum adaptive_mu_kkt_balancing_term_;
    /** enumeration for adaptive globalization */
    enum AdaptiveMuGlobalizationEnum
    {
      KKT_ERROR=0,
      FILTER_OBJ_CONSTR,
      FILTER_KKT_ERROR,
      NEVER_MONOTONE_MODE
    };
    /** Flag indicating which globalization strategy should be used. */
    AdaptiveMuGlobalizationEnum adaptive_mu_globalization_;
    /** Maximal margin in filter */
    Number filter_max_margin_;
    /** Factor for filter margin */
    Number filter_margin_fact_;
    /** Unscaled tolerance for complementarity */
    Number compl_inf_tol_;
    //@}

    /** @name Strategy objects */
    //@{
    /** Line search object of the Ipopt algorithm.  */
    SmartPtr<LineSearch> linesearch_;
    /** Pointer to strategy object that is to be used for computing a
     *  suggested value of the barrier parameter in the free mu mode.
     */
    SmartPtr<MuOracle> free_mu_oracle_;
    /** Pointer to strategy object that is to be used for computing a
     *  suggested value for the fixed mu mode.  If NULL, the current
     *  average complementarity is used.
     */
    SmartPtr<MuOracle> fix_mu_oracle_;
    //@}

    /** Dual infeasibility at initial point.  A negative value means
     *  that this quantity has not yet been initialized. */
    Number init_dual_inf_;
    /** Primal infeasibility at initial point.  A negative value means
     *  that this quantity has not yet been initialized. */
    Number init_primal_inf_;

    /** @name Methods and data defining the outer globalization
     *  strategy (might be a strategy object later). */
    //@{
    void InitializeFixedMuGlobalization();
    /** Check whether the point in the "current" fields offers
     *  sufficient reduction in order to remain in or switch to the
     *  free mu mode. */
    bool CheckSufficientProgress();
    /** Include the current point in internal memory to as accepted
     *  point */
    void RememberCurrentPointAsAccepted();
    /** Compute the value of the fixed mu that should be used in a new
     *  fixed mu phase.  This method is called at the beginning of a
     *  new fixed mu phase. */
    Number NewFixedMu();
    /** Compute value for the fraction-to-the-boundary parameter given
     *  mu in the monotone phase */
    Number Compute_tau_monotone(Number mu);

    /** Method for computing the norm of the primal dual system at the
     *  current point.  For consistency, this is computed in the same
     *  way as the quality function is computed.  This is the
     *  quantities used in the nonmonontone KKT reduction
     *  globalization. */
    Number quality_function_pd_system();

    /** Method for computing a lower safeguard bound for the barrier
     *  parameter.  For now, this is related to primal and dual
     *  infeasibility. */
    Number lower_mu_safeguard();

    /** Computer the currently largest reference value. */
    Number max_ref_val();

    /** Computer the currently smallest reference value. */
    Number min_ref_val();

    /** Maximal number of reference values (algorithmic parameter) */
    Index num_refs_max_;
    /** Values of the currently stored reference values (norm of pd
     *  equations) */
    std::list<Number> refs_vals_;
    /** Factor requested to reduce the reference values */
    Number refs_red_fact_;

    /** Alternatively, we might also want to use a filter */
    Filter filter_;
    /** Flag indicating whether the most recent accepted step should
     *  be restored, when switching to the fixed mode. */
    bool restore_accepted_iterate_;
    //@}

    /** Flag indicating whether the problem has any inequality constraints */
    bool no_bounds_;
    /** Flag indicating whether no_bounds_ has been initialized */
    bool check_if_no_bounds_;

    /** @name Most recent accepted point in free mode, from which
     *  fixed mode should be started.
     */
    //@{
    SmartPtr<const IteratesVector> accepted_point_;
    //@}

  };

} // namespace Ipopt

#endif
