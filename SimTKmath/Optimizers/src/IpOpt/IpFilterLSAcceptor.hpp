// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpFilterLSAcceptor.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter                 IBM    2005-10-13
//               derived file from IpFilterLineSearch.hpp

#ifndef __IPFILTERLSACCEPTOR_HPP__
#define __IPFILTERLSACCEPTOR_HPP__

#include "IpFilter.hpp"
#include "IpBacktrackingLSAcceptor.hpp"
#include "IpPDSystemSolver.hpp"

namespace Ipopt
{

  /** Filter line search.  This class implements the filter line
   *  search procedure. 
   */
  class FilterLSAcceptor : public BacktrackingLSAcceptor
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor.  The PDSystemSolver object only needs to be
     *  provided (i.e. not NULL) if second order correction or
     *  corrector steps are to be used. */
    FilterLSAcceptor(const SmartPtr<PDSystemSolver>& pd_solver);

    /** Default destructor */
    virtual ~FilterLSAcceptor();
    //@}

    /** InitializeImpl - overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Reset the acceptor.
     *  This function should be called if all previous information
     *  should be discarded when the line search is performed the
     *  next time.  For example, this method should be called if
     *  the barrier parameter is changed.
     */
    virtual void Reset();

    /** Initialization for the next line search.  The flag in_watchdog
     *  indicates if we are currently in an active watchdog
     *  procedure. */
    virtual void InitThisLineSearch(bool in_watchdog);

    /** Method that is called before the restoration phase is called.
     *  Here, we can set up things that are required in the
     *  termination test for the restoration phase, such as augmenting
     *  a filter. */
    virtual void PrepareRestoPhaseStart();

    /** Method returning the lower bound on the trial step sizes.  If
     *  the backtracking procedure encounters a trial step size below
     *  this value after the first trial set, it swtiches to the
     *  (soft) restoration phase. */
    virtual Number CalculateAlphaMin();

    /** Method for checking if current trial point is acceptable.
     *  It is assumed that the delta information in ip_data is the
     *  search direction used in criteria.  The primal trial point has
     *  to be set before the call.
     */
    virtual bool CheckAcceptabilityOfTrialPoint(Number alpha_primal);

    /** Try a second order correction for the constraints.  If the
     *  first trial step (with incoming alpha_primal) has been reject,
     *  this tries up to max_soc_ second order corrections for the
     *  constraints.  Here, alpha_primal_test is the step size that
     *  has to be used in the filter acceptance tests.  On output
     *  actual_delta_ has been set to the step including the
     *  second order correction if it has been accepted, otherwise it
     *  is unchanged.  If the SOC step has been accepted, alpha_primal
     *  has the fraction-to-the-boundary value for the SOC step on output.
     *  The return value is true, if a SOC step has been accepted.
     */
    virtual bool TrySecondOrderCorrection(Number alpha_primal_test,
                                          Number& alpha_primal,
                                          SmartPtr<IteratesVector>& actual_delta);

    /** Try higher order corrector (for fast local convergence).  In
     *  contrast to a second order correction step, which tries to
     *  make an unacceptable point acceptable by improving constraint
     *  violation, this corrector step is tried even if the regular
     *  primal-dual step is acceptable.
     */
    virtual bool TryCorrector(Number alpha_primal_test,
                              Number& alpha_primal,
                              SmartPtr<IteratesVector>& actual_delta);

    /** Method for ending the current line search.  When it is called,
     *  the internal data should be updates, e.g., the filter might be
     *  augmented.  alpha_primal_test is the value of alpha that has
     *  been used for in the acceptence test ealier. */
    virtual char UpdateForNextIteration(Number alpha_primal_test);

    /** Method for setting internal data if the watchdog procedure is
     *  started. */
    virtual void StartWatchDog();

    /** Method for setting internal data if the watchdog procedure is
     *  stopped. */
    virtual void StopWatchDog();

    /**@name Trial Point Accepting Methods. Used internally to check certain
     * acceptability criteria and used externally (by the restoration phase
     * convergence check object, for instance)
     */
    //@{
    /** Checks if a trial point is acceptable to the current iterate */
    bool IsAcceptableToCurrentIterate(Number trial_barr, Number trial_theta,
                                      bool called_from_restoration=false) const;

    /** Checks if a trial point is acceptable to the current filter */
    bool IsAcceptableToCurrentFilter(Number trial_barr, Number trial_theta) const;
    //@}

    /** Methods for OptionsList */
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
    /** Copy Constructor */
    FilterLSAcceptor(const FilterLSAcceptor&);

    /** Overloaded Equals Operator */
    void operator=(const FilterLSAcceptor&);
    //@}

    /** @name Filter information */
    //@{
    /** Upper bound on infeasibility */
    Number theta_max_;
    Number theta_max_fact_;

    /** Infeasibility switching bound */
    Number theta_min_;
    Number theta_min_fact_;
    //@}

    /** Method for checking if the current step size satisfies the
     *  f-type switching condition.  Here, we use the search direction
     *  stored in ip_data
     */
    bool IsFtype(Number alpha_primal_test);

    /** Method for checking the Armijo condition, given a trial step
     *  size.  The test uses the search direction stored in ip_data,
     *  and the values of the functions at the trial point in ip_data.
     */
    bool ArmijoHolds(Number alpha_primal_test);

    /** Augment the filter used on the current values of the barrier
     *  objective function and the contraint violation */
    void AugmentFilter();

    /** Check comparison "lhs <= rhs", using machine precision based on BasVal */
    //ToDo This should probably not be a static member function if we want to
    //     allow for different relaxation parameters values
    static bool Compare_le(Number lhs, Number rhs, Number BasVal);

    /** @name Parameters for the filter algorithm.  Names as in the paper */
    //@{
    /** \f$ \eta_{\varphi} \f$ */
    Number eta_phi_;
    /** \f$ \delta \f$ */
    Number delta_;
    /** \f$ s_{\varphi} \f$ */
    Number s_phi_;
    /** \f$ s_{\Theta} \f$ */
    Number s_theta_;
    /** \f$ \gamma_{\varphi} \f$ */
    Number gamma_phi_;
    /** \f$ \gamma_{\Theta} \f$ */
    Number gamma_theta_;
    /** \f$ \gamma_{\alpha} \f$ */
    Number alpha_min_frac_;
    /** Maximal number of second order correction steps */
    Index max_soc_;
    /** Required reduction in constraint violation before trying
     *  multiple second order correction steps \f$ \kappa_{soc}\f$.
     */
    Number kappa_soc_;
    /** Maximal increase in objective function in orders of magnitute
     *  (log10).  If the log10(barrier objective function) is
     *  increased more than this compared to the current point, the
     *  trial point is rejected. */
    Number obj_max_inc_;

    /** enumeration for the corrector type */
    enum CorrectorTypeEnum {
      NO_CORRECTOR=0,
      AFFINE_CORRECTOR,
      PRIMAL_DUAL_CORRECTOR
    };
    /** Type of corrector steps that should be tried. */
    CorrectorTypeEnum corrector_type_;
    /** parameter in heurstic that determines whether corrector step
    should be tried. */
    Number corrector_compl_avrg_red_fact_;
    /** Flag indicating whether the corrector should be skipped in an
     *  iteration in which negative curvature is detected */
    bool skip_corr_if_neg_curv_;
    /** Flag indicating whether the corrector should be skipped during
     *  the monotone mu mode. */
    bool skip_corr_in_monotone_mode_;
    /** maximal allowed number of filter resets. */
    Index max_filter_resets_;
    /** interation counter trigger for filter reset.  If the
     *  successive number of iterations in which the last rejected
     *  step was due to the filter, and max_filter_resets is non-zero,
     *  then the filter is reset. */
    Index filter_reset_trigger_;
    //@}

    /** @name Information related to watchdog procedure */
    //@{
    /** Constraint violation at the point with respect to which
     *  progress is to be made */
    Number reference_theta_;
    /** Barrier objective function at the point with respect to which
     *  progress is to be made */
    Number reference_barr_;
    /** Barrier gradient transpose search direction at the point with
     *  respect to which progress is to be made */
    Number reference_gradBarrTDelta_;
    /** Constraint violation at reference point */
    Number watchdog_theta_;
    /** Barrier objective function at reference point */
    Number watchdog_barr_;
    /** Barrier gradient transpose search direction at reference point */
    Number watchdog_gradBarrTDelta_;
    //@}

    /** Filter with entries */
    Filter filter_;

    /** @name Filter reset stuff */
    //@{
    /** True, if last rejected was due to the filter. */
    Number last_rejection_due_to_filter_;
    /** Counter of successive iterations in which filter was reason
     *  for last rejection. */
    Index count_successive_filter_rejections_;
    /** Counter for the filter resets done so far. */
    Index n_filter_resets_;
    //@}

    /** @name Strategy objective that are used */
    //@{
    SmartPtr<PDSystemSolver> pd_solver_;
    //@}
  };

} // namespace Ipopt

#endif
