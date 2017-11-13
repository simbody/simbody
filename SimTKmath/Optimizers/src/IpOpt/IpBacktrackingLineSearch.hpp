// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpBacktrackingLineSearch.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13
//           Andreas Waechter                 IBM    2005-10-13
//               derived file from IpFilterLineSearch.hpp

#ifndef __IPBACKTRACKINGLINESEARCH_HPP__
#define __IPBACKTRACKINGLINESEARCH_HPP__

#include "IpLineSearch.hpp"
#include "IpBacktrackingLSAcceptor.hpp"
#include "IpRestoPhase.hpp"
#include "IpConvCheck.hpp"

namespace SimTKIpopt
{

  /** General implementation of a backtracking line search.  This
   *  class can be used to perform the filter line search procedure or
   *  other procedures.  The BacktrackingLSAcceptor is used to
   *  determine whether trial points are acceptable (e.g., based on a
   *  filter or other methods).
   *
   *  This backtracking line search knows of a restoration phase
   *  (which is called when the trial step size becomes too small or
   *  no search direction could be computed).  It also has the notion
   *  of a "soft restoration phase," which uses the regular steps but
   *  decides on the acceptability based on other measures than the
   *  regular ones (e.g., reduction of the PD error instead of
   *  acceptability to a filter mechanism).
   */
  class BacktrackingLineSearch : public LineSearch
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor.  The acceptor implements the acceptance test for
     *  the line search. The PDSystemSolver object only needs to be
     *  provided (i.e. not NULL) if second order correction is to be
     *  used.  The ConvergenceCheck object is used to determine
     *  whether the current iterate is acceptable (for example, the
     *  restoration phase is not started if the acceptability level
     *  has been reached).  If conv_check is NULL, we assume that the
     *  current iterate is not acceptable (in the sense of the
     *  acceptable_tol option). */
    BacktrackingLineSearch(const SmartPtr<BacktrackingLSAcceptor>& acceptor,
                           const SmartPtr<RestorationPhase>& resto_phase,
                           const SmartPtr<ConvergenceCheck>& conv_check
                          );

    /** Default destructor */
    virtual ~BacktrackingLineSearch();
    //@}

    /** InitializeImpl - overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) override;

    /** Perform the line search.  It is assumed that the search
     *  direction is computed in the data object.
     */
    virtual void FindAcceptableTrialPoint() override;

    /** Reset the line search.
     *  This function should be called if all previous information
     *  should be discarded when the line search is performed the
     *  next time.  For example, this method should be called if
     *  the barrier parameter is changed.
     */
    virtual void Reset() override;

    /** Set flag indicating whether a very rigorous line search should
     *  be performed.  If this flag is set to true, the line search
     *  algorithm might decide to abort the line search and not to
     *  accept a new iterate.  If the line search decided not to
     *  accept a new iterate, the return value of
     *  CheckSkippedLineSearch() is true at the next call.  For
     *  example, in the non-monotone barrier parameter update
     *  procedure, the filter algorithm should not switch to the
     *  restoration phase in the free mode; instead, the algorithm
     *  should swtich to the fixed mode.
     */
    virtual void SetRigorousLineSearch(bool rigorous) override
    {
      rigorous_ = rigorous;
    }

    /** Check if the line search procedure didn't accept a new iterate
     *  during the last call of FindAcceptableTrialPoint().
     *  
     */
    virtual bool CheckSkippedLineSearch() override
    {
      return skipped_line_search_;
    }

    /** Activate fallback mechanism.  Return false, if that is not
     *  possible. */
    virtual bool ActivateFallbackMechanism() override;

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
    BacktrackingLineSearch(const BacktrackingLineSearch&);

    /** Overloaded Equals Operator */
    void operator=(const BacktrackingLineSearch&);
    //@}

    /** Method performing the backtracking line search.  The return
     *  value indicates if the step acceptance criteria are met.  If
     *  the watchdog is active, only one trial step is performed (and
     *  the trial values are set accordingly). */
    bool DoBacktrackingLineSearch(bool skip_first_trial_point,
                                  Number& alpha_primal,
                                  bool& corr_taken,
                                  bool& soc_taken,
                                  Index& n_steps,
                                  bool& evaluation_error,
                                  SmartPtr<IteratesVector>& actual_delta);

    /** Method for starting the watch dog.  Set all appropriate fields
     *  accordingly */
    void StartWatchDog();

    /** Method for stopping the watch dog.  Set all appropriate fields
     *  accordingly. */
    void StopWatchDog(SmartPtr<IteratesVector>& actual_delta);

    /** Method for checking if current trial point is acceptable.
     *  It is assumed that the delta information in ip_data is the
     *  search direction used in criteria.  The primal trial point has
     *  to be set before the call.
     */
    bool CheckAcceptabilityOfTrialPoint(Number alpha_primal);

    /** Check comparison "lhs <= rhs", using machine precision based on BasVal */
    //ToDo This should probably not be a static member function if we want to
    //     allow for different relaxation parameters values
    static bool Compare_le(Number lhs, Number rhs, Number BasVal);

    /** Method for setting the dual variables in the trial fields in
     *  IpData, given the search direction.  The step size for the
     *  bound multipliers is alpha_dual (the fraction-to-the-boundary
     *  step size), and the step size for the equality constraint 
     *  multipliers depends on the choice of alpha_for_y. */
    void PerformDualStep(Number alpha_primal,
                         Number alpha_dual,
                         SmartPtr<IteratesVector>& delta);

    /** Try a step for the soft restoration phase and check if it is
     *  acceptable.  The step size is identical for all variables.  A
     *  point is accepted if it is acceptable for the original
     *  acceptability criterion (in which case
     *  satisfies_original_criterion = true on return), or if the
     *  primal-dual system error was decrease by at least the factor
     *  soft_resto_pderror_reduction_factor_.  The return value is
     *  true, if the trial point was acceptable for the soft
     *  restoration phase. */
    bool TrySoftRestoStep(SmartPtr<IteratesVector>& actual_delta,
                          bool &satisfies_original_criterion);

    /** Try a second order correction for the constraints.  If the
     *  first trial step (with incoming alpha_primal) has been reject,
     *  this tries up to max_soc_ second order corrections for the
     *  constraints.  Here, alpha_primal_test is the step size that
     *  has to be used in the filter acceptance tests.  On output
     *  actual_delta_... has been set to the steps including the
     *  second order correction if it has been accepted, otherwise it
     *  is unchanged.  If the SOC step has been accepted, alpha_primal
     *  has the fraction-to-the-boundary value for the SOC step on output.
     *  The return value is true, if an SOC step has been accepted.
     */
    bool TrySecondOrderCorrection(Number alpha_primal_test,
                                  Number& alpha_primal,
                                  SmartPtr<IteratesVector>& actual_delta);

    /** Try higher order corrector (for fast local convergence).  In
     *  contrast to a second order correction step, which tries to
     *  make an unacceptable point acceptable by improving constraint
     *  violation, this corrector step is tried even if the regular
     *  primal-dual step is acceptable.
     */
    bool TryCorrector(Number alpha_primal_test,
                      Number& alpha_primal,
                      SmartPtr<IteratesVector>& actual_delta);

    /** Perform magic steps.  Take the current values of the slacks in
     *  trial and replace them by better ones that lead to smaller
     *  values of the barrier function and less constraint
     *  violation. */
    void PerformMagicStep();

    /** Detect if the search direction is too small.  This should be
     *  true if the search direction is so small that if makes
     *  numerically no difference. */
    bool DetectTinyStep();

    /** Store current iterate as acceptable point */
    void StoreAcceptablePoint();

    /** Restore acceptable point into the current fields of IpData if
     *  found. Returns true if such as point is available. */
    bool RestoreAcceptablePoint();

    /** Method for determining if the current iterate is acceptable
     *  (in the sense of the acceptable_tol options).  This is a
     *  wrapper for same method from ConvergenceCheck, but returns
     *  false, if no ConvergenceCheck object is provided. */
    bool CurrentIsAcceptable();

    /** @name Parameters for the filter algorithm.  Names as in the paper */
    //@{
    /** factor by which search direction is to be shortened if trial
     *  point is rejected. */
    Number alpha_red_factor_;

    /** enumeration for the different alpha_for_y_ settings */
    enum AlphaForYEnum {
      PRIMAL_ALPHA_FOR_Y=0,
      DUAL_ALPHA_FOR_Y,
      MIN_ALPHA_FOR_Y,
      MAX_ALPHA_FOR_Y,
      FULL_STEP_FOR_Y,
      MIN_DUAL_INFEAS_ALPHA_FOR_Y,
      SAFE_MIN_DUAL_INFEAS_ALPHA_FOR_Y
    };
    /** Flag indicating whether the dual step size is to be used for
     *  the equality constraint multipliers. If 0, the primal step
     *  size is used, if 1 the dual step size, and if 2, the minimum
     *  of both. */
    AlphaForYEnum alpha_for_y_;

    /** Reduction factor for the restoration phase that accepts steps
     *  reducing the optimality error ("soft restoration phase"). If
     *  0., then this restoration phase is not enabled. */
    Number soft_resto_pderror_reduction_factor_;
    /** Maximal number of iterations that can be done in the soft
     *  iteration phase before the algorithm reverts to the regular
     *  restoration phase. */
    Index max_soft_resto_iters_;

    /** Flag indicating whether magic steps should be used. */
    bool magic_steps_;
    /** Flag indicating whether the line search should always accept
     *  the full (fraction-to-the-boundary) step. */
    bool accept_every_trial_step_;
    /** Indicates whether problem can be expected to be infeasible.
     *  This will trigger requesting a tighter reduction in
     *  infeasibility the first time the restoration phase is
     *  called. */
    bool expect_infeasible_problem_;
    /** Tolerance on constraint violation for
     *  expect_infeasible_problem heuristic.  If the constraint
     *  violation becomes that than this value, the heuristic is
     *  disabled for the rest of the optimization run. */
    Number expect_infeasible_problem_ctol_;

    /** Tolerance for detecting tiny steps. */
    Number tiny_step_tol_;

    /** Tolerance for y variables for the tiny step stopping
     *  heuristic.  If repeatedly a tiny step is detected and the step
     *  in the y_c and y_d variables is less than this threshold, we
     *  algorithm will stop. */
    Number tiny_step_y_tol_;

    /** Number of watch dog trial steps. */
    Index watchdog_trial_iter_max_;
    /** Number of shortened iterations that trigger the watchdog. */
    Index watchdog_shortened_iter_trigger_;

    /** Indicates whether the algorithm should start directly with the
     *  restoratin phase */
    bool start_with_resto_;
    //@}

    /** @name Information related to watchdog procedure */
    //@{
    /** Flag indicating if the watchdog is active */
    bool in_watchdog_;
    /** Counter for shortened iterations. */
    Index watchdog_shortened_iter_;
    /** Counter for watch dog iterations */
    Index watchdog_trial_iter_;
    /** Step size for Armijo test in watch dog */
    Number watchdog_alpha_primal_test_;
    /** Watchdog reference iterate */
    SmartPtr<const IteratesVector> watchdog_iterate_;
    /** Watchdog search direction at reference point */
    SmartPtr<const IteratesVector> watchdog_delta_;
    /** Barrier parameter value during last line search */
    Number last_mu_;
    //@}

    /** @name Storage for last iterate that satisfies the acceptable
     *  level of optimality error. */
    //@{
    SmartPtr<const IteratesVector> acceptable_iterate_;
    Index acceptable_iteration_number_;
    //@}

    /** Flag indicating whether the algorithm has asked to immediately
     *  switch to the fallback mechanism (restoration phase) */
    bool fallback_activated_;

    /** Flag indicating whether the line search is to be performed
     * robust (usually this is true, unless SetRigorousLineSearch is
     * called with false).
     */
    bool rigorous_;

    /** Flag indicating whether no acceptable trial point was found
     *  during last line search. */
    bool skipped_line_search_;

    /** Flag indicating whether we are currently in the "soft"
     *  restoration phase mode, in which steps are accepted if they
     *  reduce the optimality error (see
     *  soft_resto_pderror_reduction_factor) */
    bool in_soft_resto_phase_;

    /** Counter for iteration performed in soft restoration phase in a
     *  row */
    Index soft_resto_counter_;

    /** Counter for the number of successive iterations in which the
     *  full step was not accepted. */
    Index count_successive_shortened_steps_;

    /** Flag indicating if a tiny step was detected in previous
     *  iteration */
    bool tiny_step_last_iteration_;

    /** @name Strategy objective that are used */
    //@{
    SmartPtr<BacktrackingLSAcceptor> acceptor_;
    SmartPtr<RestorationPhase> resto_phase_;
    SmartPtr<ConvergenceCheck> conv_check_;
    //@}
  };

} // namespace Ipopt

#endif
