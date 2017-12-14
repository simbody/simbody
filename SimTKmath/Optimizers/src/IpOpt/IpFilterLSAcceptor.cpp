// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpFilterLSAcceptor.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13
//           Andreas Waechter                 IBM    2005-10-13
//               derived file from IpFilterLineSearch.cpp

#include "IpFilterLSAcceptor.hpp"
#include "IpJournalist.hpp"
#include "IpRestoPhase.hpp"
#include "IpAlgTypes.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

namespace SimTKIpopt
{

#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  FilterLSAcceptor::FilterLSAcceptor(const SmartPtr<PDSystemSolver>& pd_solver)
      :
      filter_(2),
      pd_solver_(pd_solver)
  {
    DBG_START_FUN("FilterLSAcceptor::FilterLSAcceptor",
                  dbg_verbosity);
  }

  FilterLSAcceptor::~FilterLSAcceptor()
  {
    DBG_START_FUN("FilterLSAcceptor::~FilterLSAcceptor()",
                  dbg_verbosity);
  }

  void FilterLSAcceptor::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "theta_max_fact",
      "Determines upper bound for constraint violation in the filter.",
      0.0, true, 1e4,
      "The algorithmic parameter theta_max is determined as theta_max_fact "
      "times the maximum of 1 and the constraint violation at initial point.  "
      "Any point with a constraint violation larger than theta_max is "
      "unacceptable to the filter (see Eqn. (21) in implementation paper).");
    roptions->AddLowerBoundedNumberOption(
      "theta_min_fact",
      "Determines constraint violation threshold in the switching rule.",
      0.0, true, 1e-4,
      "The algorithmic parameter theta_min is determined as theta_min_fact "
      "times the maximum of 1 and the constraint violation at initial point.  "
      "The switching rules treats an iteration as an h-type iteration whenever "
      "the current constraint violation is larger than theta_min (see "
      "paragraph before Eqn. (19) in implementation paper).");
    roptions->AddBoundedNumberOption(
      "eta_phi",
      "Relaxation factor in the Armijo condition.",
      0.0, true, 0.5, true, 1e-8,
      "(See Eqn. (20) in implementation paper)");
    roptions->AddLowerBoundedNumberOption(
      "delta", "Multiplier for constraint violation in the switching rule.",
      0.0, true, 1.0,
      "(See Eqn. (19) in implementation paper.)");
    roptions->AddLowerBoundedNumberOption(
      "s_phi",
      "Exponent for linear barrier function model in the switching rule.",
      1.0, true, 2.3,
      "(See Eqn. (19) in implementation paper.)");
    roptions->AddLowerBoundedNumberOption(
      "s_theta",
      "Exponent for current constraint violation in the switching rule.",
      1.0, true, 1.1,
      "(See Eqn. (19) in implementation paper.)");
    roptions->AddBoundedNumberOption(
      "gamma_phi",
      "Relaxation factor in the filter margin for the barrier function.",
      0.0, true, 1.0, true, 1e-8,
      "(See Eqn. (18a) in implementation paper.)");
    roptions->AddBoundedNumberOption(
      "gamma_theta",
      "Relaxation factor in the filter margin for the constraint violation.",
      0.0, true, 1.0, true, 1e-5,
      "(See Eqn. (18b) in implementation paper.)");
    roptions->AddBoundedNumberOption(
      "alpha_min_frac",
      "Safety factor for the minimal step size (before switching to restoration phase).",
      0.0, true, 1.0, true, 0.05,
      "(This is gamma_alpha in Eqn. (20) in implementation paper.)");
    roptions->AddLowerBoundedIntegerOption(
      "max_soc",
      "Maximum number of second order correction trial steps at each iteration.",
      0, 4,
      "Choosing 0 disables the second order "
      "corrections. (This is p^{max} of Step A-5.9 of "
      "Algorithm A in implementation paper.)");
    roptions->AddLowerBoundedNumberOption(
      "kappa_soc",
      "Factor in the sufficient reduction rule for second order correction.",
      0.0, true, 0.99,
      "This option determines how much a second order correction step must reduce the "
      "constraint violation so that further correction steps are attempted.  "
      "(See Step A-5.9 of Algorithm A in implementation paper.)");
    roptions->AddLowerBoundedNumberOption(
      "obj_max_inc",
      "Determines the upper bound on the acceptable increase of barrier objective function.",
      1.0, true, 5.0,
      "Trial points are rejected if they lead to an increase in the "
      "barrier objective function by more than obj_max_inc orders "
      "of magnitude.");

    roptions->AddLowerBoundedIntegerOption(
      "max_filter_resets",
      "Maximal allowed number of filter resets",
      0, 5,
      "A positive number enables a heuristic that resets the filter, whenever "
      "in more than \"filter_reset_trigger\" successive iterations the last "
      "rejected trial steps size was rejected because of the filter.  This "
      "option determine the maximal number of resets that are allowed to take "
      "place.");
    roptions->AddLowerBoundedIntegerOption(
      "filter_reset_trigger",
      "Number of iterations that trigger the filter reset.",
      1, 5,
      "If the filter reset heuristic is active and the number of successive "
      "iterations in which the last rejected trial step size was rejected "
      "because of the filter, the filter is reset.");

    roptions->AddStringOption3(
      "corrector_type",
      "The type of corrector steps that should be taken (unsupported!).",
      "none",
      "none", "no corrector",
      "affine", "corrector step towards mu=0",
      "primal-dual", "corrector step towards current mu",
      "If \"mu_strategy\" is \"adaptive\", this option determines "
      "what kind of corrector steps should be tried.");

    roptions->AddStringOption2(
      "skip_corr_if_neg_curv",
      "Skip the corrector step in negative curvature iteration (unsupported!).",
      "yes",
      "no", "don't skip",
      "yes", "skip",
      "The corrector step is not tried if negative curvature has been "
      "encountered during the computation of the search direction in "
      "the current iteration. This option is only used if \"mu_strategy\" is "
      "\"adaptive\".");

    roptions->AddStringOption2(
      "skip_corr_in_monotone_mode",
      "Skip the corrector step during monotone barrier parameter mode (unsupported!).",
      "yes",
      "no", "don't skip",
      "yes", "skip",
      "The corrector step is not tried if the algorithm is currently in the "
      "monotone mode (see also option \"barrier_strategy\")."
      "This option is only used if \"mu_strategy\" is \"adaptive\".");

    roptions->AddLowerBoundedNumberOption(
      "corrector_compl_avrg_red_fact",
      "Complementarity tolerance factor for accepting corrector step (unsupported!).",
      0.0, true, 1.0,
      "This option determines the factor by which complementarity is allowed to increase "
      "for a corrector step to be accepted.");
  }

  bool FilterLSAcceptor::InitializeImpl(const OptionsList& options,
                                        const std::string& prefix)
  {
    options.GetNumericValue("theta_max_fact", theta_max_fact_, prefix);
    options.GetNumericValue("theta_min_fact", theta_min_fact_, prefix);
    ASSERT_EXCEPTION(theta_min_fact_ < theta_max_fact_, OPTION_INVALID,
                     "Option \"theta_min_fact\": This value must be larger than 0 and less than theta_max_fact.");
    options.GetNumericValue("eta_phi", eta_phi_, prefix);
    options.GetNumericValue("delta", delta_, prefix);
    options.GetNumericValue("s_phi", s_phi_, prefix);
    options.GetNumericValue("s_theta", s_theta_, prefix);
    options.GetNumericValue("gamma_phi", gamma_phi_, prefix);
    options.GetNumericValue("gamma_theta", gamma_theta_, prefix);
    options.GetNumericValue("alpha_min_frac", alpha_min_frac_, prefix);
    options.GetIntegerValue("max_soc", max_soc_, prefix);
    if (max_soc_>0) {
      ASSERT_EXCEPTION(IsValid(pd_solver_), OPTION_INVALID,
                       "Option \"max_soc\": This option is non-negative, but no linear solver for computing the SOC given to FilterLSAcceptor object.");
    }
    options.GetNumericValue("kappa_soc", kappa_soc_, prefix);
    options.GetIntegerValue("max_filter_resets", max_filter_resets_, prefix);
    options.GetIntegerValue("filter_reset_trigger", filter_reset_trigger_,
                            prefix);
    options.GetNumericValue("obj_max_inc", obj_max_inc_, prefix);
    Index enum_int;
    options.GetEnumValue("corrector_type", enum_int, prefix);
    corrector_type_ = CorrectorTypeEnum(enum_int);
    options.GetBoolValue("skip_corr_if_neg_curv", skip_corr_if_neg_curv_, prefix);
    options.GetBoolValue("skip_corr_in_monotone_mode", skip_corr_in_monotone_mode_, prefix);
    options.GetNumericValue("corrector_compl_avrg_red_fact", corrector_compl_avrg_red_fact_, prefix);

    theta_min_ = -1.;
    theta_max_ = -1.;

    n_filter_resets_ = 0;

    Reset();

    return true;
  }

  void FilterLSAcceptor::InitThisLineSearch(bool in_watchdog)
  {
    DBG_START_METH("FilterLSAcceptor::InitThisLineSearch",
                   dbg_verbosity);

    // Set the values for the reference point
    if (!in_watchdog) {
      reference_theta_ = IpCq().curr_constraint_violation();
      reference_barr_ = IpCq().curr_barrier_obj();
      reference_gradBarrTDelta_ = IpCq().curr_gradBarrTDelta();
    }
    else {
      reference_theta_ = watchdog_theta_;
      reference_barr_ = watchdog_barr_;
      reference_gradBarrTDelta_ = watchdog_gradBarrTDelta_;
    }
    filter_.Print(Jnlst());
  }

  bool FilterLSAcceptor::IsFtype(Number alpha_primal_test)
  {
    DBG_START_METH("FilterLSAcceptor::IsFtype",
                   dbg_verbosity);
    Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                   "reference_theta = %e reference_gradBarrTDelta = %e\n",
                   reference_theta_, reference_gradBarrTDelta_);
    Number mach_eps = std::numeric_limits<Number>::epsilon();
    // ToDo find good value
    if (reference_theta_==0. &&  reference_gradBarrTDelta_ > 0. &&
        reference_gradBarrTDelta_ < 100.*mach_eps) {
      reference_gradBarrTDelta_ = -mach_eps;
      Jnlst().Printf(J_WARNING, J_LINE_SEARCH,
                     "reference_theta is slightly positive at feasible point.  Setting it to %e\n",
                     reference_gradBarrTDelta_);
    }
    DBG_ASSERT(reference_theta_>0. || reference_gradBarrTDelta_ < 0.0);
    return (reference_gradBarrTDelta_ < 0.0 &&
            alpha_primal_test*pow(-reference_gradBarrTDelta_,s_phi_) >
            delta_*pow(reference_theta_,s_theta_));
  }

  void FilterLSAcceptor::AugmentFilter()
  {
    DBG_START_METH("FilterLSAcceptor::AugmentFilter",
                   dbg_verbosity);

    Number phi_add = reference_barr_ - gamma_phi_*reference_theta_;
    Number theta_add = (1.-gamma_theta_)*reference_theta_;

    filter_.AddEntry(phi_add, theta_add, IpData().iter_count());
  }

  bool
  FilterLSAcceptor::CheckAcceptabilityOfTrialPoint(Number alpha_primal_test)
  {
    DBG_START_METH("FilterLSAcceptor::CheckAcceptabilityOfTrialPoint",
                   dbg_verbosity);

    bool accept;

    // First compute the barrier function and constraint violation at the
    // current iterate and the trial point

    Number trial_theta = IpCq().trial_constraint_violation();
    // Check if constraint violation is becoming too large
    if (theta_max_ < 0.0) {
      // ToDo should 1.0 be based on dimension? (theta is in 1 norm!!!)
      theta_max_ = theta_max_fact_*Max(1.0, reference_theta_);
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "trial_max is initialized to %e\n",
                     theta_max_);
    }
    if (theta_min_ < 0.0) {
      theta_min_ = theta_min_fact_*Max(1.0, reference_theta_);
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "trial_min is initialized to %e\n",
                     theta_min_);
    }

    if (theta_max_>0 && trial_theta>theta_max_) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "trial_theta = %e is larger than theta_max = %e\n",
                     trial_theta, theta_max_);
      IpData().Append_info_string("Tmax");
      return false;
    }

    Number trial_barr = IpCq().trial_barrier_obj();
    DBG_ASSERT(IsFiniteNumber(trial_barr));

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Checking acceptability for trial step size alpha_primal_test=%13.6e:\n", alpha_primal_test);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  New values of barrier function     = %23.16e  (reference %23.16e):\n", trial_barr, reference_barr_);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  New values of constraint violation = %23.16e  (reference %23.16e):\n", trial_theta, reference_theta_);

    // Check if point is acceptable w.r.t current iterate
    if (alpha_primal_test>0. && IsFtype(alpha_primal_test) &&
        reference_theta_ <= theta_min_) {
      // Armijo condition for the barrier function has to be satisfied
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Checking Armijo Condition...\n");
      accept = ArmijoHolds(alpha_primal_test);
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Checking sufficient reduction...\n");
      accept = IsAcceptableToCurrentIterate(trial_barr, trial_theta);
    }

    if (!accept) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Failed...\n");
      last_rejection_due_to_filter_ = false;
      return accept;
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Succeeded...\n");
    }

    // Now check if that pair is acceptable to the filter
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Checking filter acceptability...\n");
    accept = IsAcceptableToCurrentFilter(trial_barr, trial_theta);
    if (!accept) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Failed...\n");
      last_rejection_due_to_filter_ = true;
      return accept;
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Succeeded...\n");
    }

    // Filter reset heuristic
    if (max_filter_resets_>0) {
      if (n_filter_resets_<max_filter_resets_) {
        if (last_rejection_due_to_filter_) {
          count_successive_filter_rejections_++;
          if (count_successive_filter_rejections_>=filter_reset_trigger_) {
            Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                           "Resetting filter because in %d iterations last rejection was due to filter", count_successive_filter_rejections_);
            IpData().Append_info_string("F+");
            Reset();
          }
        }
        else {
          count_successive_filter_rejections_ = 0;
        }
      }
      else {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "Filter should be reset, but maximal number of resets already exceeded.\n");
        IpData().Append_info_string("F-");
      }
    }
    last_rejection_due_to_filter_= false;

    return accept;
  }

  bool FilterLSAcceptor::ArmijoHolds(Number alpha_primal_test)
  {
    /*
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "ArmijoHolds test with trial_barr = %25.16e reference_barr = %25.16e\n        alpha_primal_test = %25.16e reference_gradBarrTDelta = %25.16e\n", IpCq().trial_barrier_obj(), reference_barr_,alpha_primal_test,reference_gradBarrTDelta_);
    */
    return Compare_le(IpCq().trial_barrier_obj()-reference_barr_,
                      eta_phi_*alpha_primal_test*reference_gradBarrTDelta_,
                      reference_barr_);
  }

  Number FilterLSAcceptor::CalculateAlphaMin()
  {
    Number gBD = IpCq().curr_gradBarrTDelta();
    Number curr_theta = IpCq().curr_constraint_violation();
    Number alpha_min = gamma_theta_;

    if (gBD < 0) {
      alpha_min = Min( gamma_theta_,
                       gamma_phi_*curr_theta/(-gBD));
      if (curr_theta <= theta_min_) {
        alpha_min = Min( alpha_min,
                         delta_*pow(curr_theta,s_theta_)/pow(-gBD,s_phi_)
                       );
      }
    }

    return alpha_min_frac_*alpha_min;
  }

  bool FilterLSAcceptor::IsAcceptableToCurrentIterate(Number trial_barr,
      Number trial_theta,
      bool called_from_restoration /*=false*/) const
  {
    DBG_START_METH("FilterLSAcceptor::IsAcceptableToCurrentIterate",
                   dbg_verbosity);

    // Check if the barrier objective function is increasing to
    // rapidly (according to option obj_max_inc)
    if (!called_from_restoration && trial_barr > reference_barr_) {
      Number basval = 1.;
      if (fabs(reference_barr_)>10.) {
        basval = log10(fabs(reference_barr_));
      }
      if (log10(trial_barr-reference_barr_)>obj_max_inc_+basval) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Rejecting trial point because barrier objective function increasing too rapidly (from %27.15e to %27.15e)\n",reference_barr_,trial_barr);
        return false;
      }
    }

    DBG_PRINT((1,"trial_barr  = %e reference_barr  = %e\n", trial_barr, reference_barr_));
    DBG_PRINT((1,"trial_theta = %e reference_theta = %e\n", trial_theta, reference_theta_));
    return (Compare_le(trial_theta, (1.-gamma_theta_)*reference_theta_, reference_theta_)
            || Compare_le(trial_barr-reference_barr_, -gamma_phi_*reference_theta_, reference_barr_));
  }

  bool FilterLSAcceptor::IsAcceptableToCurrentFilter(Number trial_barr, Number trial_theta) const
  {
    return filter_.Acceptable(trial_barr, trial_theta);
  }

  bool FilterLSAcceptor::Compare_le(Number lhs, Number rhs, Number BasVal)
  {
    DBG_START_FUN("FilterLSAcceptor::Compare_le",
                  dbg_verbosity);
    DBG_PRINT((1,"lhs = %27.16e rhs = %27.16e  BasVal = %27.16e\n",lhs,rhs,BasVal));

    Number mach_eps = std::numeric_limits<Number>::epsilon();
    return (lhs - rhs <= 10.*mach_eps*fabs(BasVal));
  }

  void FilterLSAcceptor::StartWatchDog()
  {
    DBG_START_FUN("FilterLSAcceptor::StartWatchDog", dbg_verbosity);

    watchdog_theta_ = IpCq().curr_constraint_violation();
    watchdog_barr_ = IpCq().curr_barrier_obj();
    watchdog_gradBarrTDelta_ = IpCq().curr_gradBarrTDelta();
  }

  void FilterLSAcceptor::StopWatchDog()
  {
    DBG_START_FUN("FilterLSAcceptor::StopWatchDog", dbg_verbosity);

    reference_theta_ = watchdog_theta_;
    reference_barr_ = watchdog_barr_;
    reference_gradBarrTDelta_ = watchdog_gradBarrTDelta_;
  }

  void FilterLSAcceptor::Reset()
  {
    DBG_START_FUN("FilterLSAcceptor::Reset", dbg_verbosity);

    last_rejection_due_to_filter_ = false;
    count_successive_filter_rejections_ = 0;

    filter_.Clear();
  }

  bool
  FilterLSAcceptor::TrySecondOrderCorrection(
    Number alpha_primal_test,
    Number& alpha_primal,
    SmartPtr<IteratesVector>& actual_delta)
  {
    DBG_START_METH("FilterLSAcceptor::TrySecondOrderCorrection",
                   dbg_verbosity);

    if (max_soc_==0) {
      return false;
    }

    bool accept = false;
    Index count_soc = 0;

    Number theta_soc_old = 0.;
    Number theta_trial = IpCq().trial_constraint_violation();
    Number alpha_primal_soc = alpha_primal;

    SmartPtr<Vector> c_soc = IpCq().curr_c()->MakeNew();
    SmartPtr<Vector> dms_soc = IpCq().curr_d_minus_s()->MakeNew();
    c_soc->Copy(*IpCq().curr_c());
    dms_soc->Copy(*IpCq().curr_d_minus_s());
    while (count_soc<max_soc_ && !accept &&
           (count_soc==0 || theta_trial<=kappa_soc_*theta_soc_old) ) {
      theta_soc_old = theta_trial;

      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Trying second order correction number %d\n",
                     count_soc+1);

      // Compute SOC constraint violation
      c_soc->AddOneVector(1.0, *IpCq().trial_c(), alpha_primal_soc);
      dms_soc->AddOneVector(1.0, *IpCq().trial_d_minus_s(), alpha_primal_soc);

      // Compute the SOC search direction
      SmartPtr<IteratesVector> delta_soc = actual_delta->MakeNewIteratesVector(true);
      SmartPtr<IteratesVector> rhs = actual_delta->MakeNewContainer();
      rhs->Set_x(*IpCq().curr_grad_lag_with_damping_x());
      rhs->Set_s(*IpCq().curr_grad_lag_with_damping_s());
      rhs->Set_y_c(*c_soc);
      rhs->Set_y_d(*dms_soc);
      rhs->Set_z_L(*IpCq().curr_relaxed_compl_x_L());
      rhs->Set_z_U(*IpCq().curr_relaxed_compl_x_U());
      rhs->Set_v_L(*IpCq().curr_relaxed_compl_s_L());
      rhs->Set_v_U(*IpCq().curr_relaxed_compl_s_U());

      bool retval = pd_solver_->Solve(-1.0, 0.0, *rhs, *delta_soc, true);
      if (!retval) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "The linear system could not be solved for the corrector step.\n");
        return false;
      }

      // Compute step size
      alpha_primal_soc =
        IpCq().primal_frac_to_the_bound(IpData().curr_tau(),
                                        *delta_soc->x(),
                                        *delta_soc->s());

      // Check if trial point is acceptable
      try {
        // Compute the primal trial point
        IpData().SetTrialPrimalVariablesFromStep(alpha_primal_soc, *delta_soc->x(), *delta_soc->s());

        // in acceptance tests, use original step size!
        accept = CheckAcceptabilityOfTrialPoint(alpha_primal_test);
      }
      catch(IpoptNLP::Eval_Error& e) {
        e.ReportException(Jnlst(), J_DETAILED);
        Jnlst().Printf(J_WARNING, J_MAIN, "Warning: SOC step rejected due to evaluation error\n");
        IpData().Append_info_string("e");
        accept = false;
        // There is no point in continuing SOC procedure
        break;
      }

      if (accept) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Second order correction step accepted with %d corrections.\n", count_soc+1);
        // Accept all SOC quantities
        alpha_primal = alpha_primal_soc;
        actual_delta = delta_soc;
      }
      else {
        count_soc++;
        theta_trial = IpCq().trial_constraint_violation();
      }
    }
    return accept;
  }

  bool
  FilterLSAcceptor::TryCorrector(
    Number alpha_primal_test,
    Number& alpha_primal,
    SmartPtr<IteratesVector>& actual_delta)
  {
    if (corrector_type_==NO_CORRECTOR ||
        (skip_corr_if_neg_curv_ && IpData().info_regu_x()!=0.) ||
        (skip_corr_in_monotone_mode_ && !IpData().FreeMuMode())) {
      return false;
    }

    DBG_START_METH("FilterLSAcceptor::TryCorrector",
                   dbg_verbosity);

    Index n_bounds = IpData().curr()->z_L()->Dim() + IpData().curr()->z_U()->Dim()
                     + IpData().curr()->v_L()->Dim() + IpData().curr()->v_U()->Dim();
    if (n_bounds==0) {
      // Nothing to be done
      return false;
    }

    IpData().TimingStats().TryCorrector().Start();

    bool accept = false;

    // Compute the corrector step based on corrector_type parameter
    // create a new iterates vector and allocate space for all the entries
    SmartPtr<IteratesVector> delta_corr = actual_delta->MakeNewIteratesVector(true);

    switch (corrector_type_) {
      case AFFINE_CORRECTOR : {
        // 1: Standard MPC corrector

        if (!IpData().HaveAffineDeltas()) {
          Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                         "Solving the Primal Dual System for the affine step\n");
          // First get the right hand side
          SmartPtr<IteratesVector> rhs_aff = delta_corr->MakeNewContainer();

          rhs_aff->Set_x(*IpCq().curr_grad_lag_x());
          rhs_aff->Set_s(*IpCq().curr_grad_lag_s());
          rhs_aff->Set_y_c(*IpCq().curr_c());
          rhs_aff->Set_y_d(*IpCq().curr_d_minus_s());
          rhs_aff->Set_z_L(*IpCq().curr_compl_x_L());
          rhs_aff->Set_z_U(*IpCq().curr_compl_x_U());
          rhs_aff->Set_v_L(*IpCq().curr_compl_s_L());
          rhs_aff->Set_v_U(*IpCq().curr_compl_s_U());

          // create a new iterates vector (with allocated space)
          // for the affine scaling step
          SmartPtr<IteratesVector> step_aff = delta_corr->MakeNewIteratesVector(true);

          // Now solve the primal-dual system to get the step
          pd_solver_->Solve(-1.0, 0.0, *rhs_aff, *step_aff, false);

          DBG_PRINT_VECTOR(2, "step_aff", *step_aff);

          IpData().set_delta_aff(step_aff);
          IpData().SetHaveAffineDeltas(true);
        }

        DBG_ASSERT(IpData().HaveAffineDeltas());

        const SmartPtr<const IteratesVector> delta_aff = IpData().delta_aff();

        delta_corr->Copy(*actual_delta);

        // create a rhs vector and allocate entries
        SmartPtr<IteratesVector> rhs = actual_delta->MakeNewIteratesVector(true);

        rhs->x_NonConst()->Set(0.);
        rhs->s_NonConst()->Set(0.);
        rhs->y_c_NonConst()->Set(0.);
        rhs->y_d_NonConst()->Set(0.);
        IpNLP().Px_L()->TransMultVector(-1., *delta_aff->x(), 0., *rhs->z_L_NonConst());
        rhs->z_L_NonConst()->ElementWiseMultiply(*delta_aff->z_L());
        IpNLP().Px_U()->TransMultVector(1., *delta_aff->x(), 0., *rhs->z_U_NonConst());
        rhs->z_U_NonConst()->ElementWiseMultiply(*delta_aff->z_U());
        IpNLP().Pd_L()->TransMultVector(-1., *delta_aff->s(), 0., *rhs->v_L_NonConst());
        rhs->v_L_NonConst()->ElementWiseMultiply(*delta_aff->v_L());
        IpNLP().Pd_U()->TransMultVector(1., *delta_aff->s(), 0., *rhs->v_U_NonConst());
        rhs->v_U_NonConst()->ElementWiseMultiply(*delta_aff->v_U());

        pd_solver_->Solve(1.0, 1.0, *rhs, *delta_corr, true);

        DBG_PRINT_VECTOR(2, "delta_corr", *delta_corr);
      }
      break;
      case PRIMAL_DUAL_CORRECTOR : {
        // 2: Second order correction for primal-dual step to
        // primal-dual mu

        delta_corr->Copy(*actual_delta);

        // allocate space for the rhs
        SmartPtr<IteratesVector> rhs = actual_delta->MakeNewIteratesVector(true);

        rhs->x_NonConst()->Set(0.);
        rhs->s_NonConst()->Set(0.);
        rhs->y_c_NonConst()->Set(0.);
        rhs->y_d_NonConst()->Set(0.);

        Number mu = IpData().curr_mu();
        SmartPtr<Vector> tmp;

        rhs->z_L_NonConst()->Copy(*IpCq().curr_slack_x_L());
        IpNLP().Px_L()->TransMultVector(-1., *actual_delta->x(),
                                        -1., *rhs->z_L_NonConst());
        tmp = actual_delta->z_L()->MakeNew();
        tmp->AddTwoVectors(1., *IpData().curr()->z_L(), 1., *actual_delta->z_L(), 0.);
        rhs->z_L_NonConst()->ElementWiseMultiply(*tmp);
        rhs->z_L_NonConst()->AddScalar(mu);

        rhs->z_U_NonConst()->Copy(*IpCq().curr_slack_x_U());
        IpNLP().Px_U()->TransMultVector(1., *actual_delta->x(),
                                        -1., *rhs->z_U_NonConst());
        tmp = actual_delta->z_U()->MakeNew();
        tmp->AddTwoVectors(1., *IpData().curr()->z_U(), 1., *actual_delta->z_U(), 0.);
        rhs->z_U_NonConst()->ElementWiseMultiply(*tmp);
        rhs->z_U_NonConst()->AddScalar(mu);

        rhs->v_L_NonConst()->Copy(*IpCq().curr_slack_s_L());
        IpNLP().Pd_L()->TransMultVector(-1., *actual_delta->s(),
                                        -1., *rhs->v_L_NonConst());
        tmp = actual_delta->v_L()->MakeNew();
        tmp->AddTwoVectors(1., *IpData().curr()->v_L(), 1., *actual_delta->v_L(), 0.);
        rhs->v_L_NonConst()->ElementWiseMultiply(*tmp);
        rhs->v_L_NonConst()->AddScalar(mu);

        rhs->v_U_NonConst()->Copy(*IpCq().curr_slack_s_U());
        IpNLP().Pd_U()->TransMultVector(1., *actual_delta->s(),
                                        -1., *rhs->v_U_NonConst());
        tmp = actual_delta->v_U()->MakeNew();
        tmp->AddTwoVectors(1., *IpData().curr()->v_U(), 1., *actual_delta->v_U(), 0.);
        rhs->v_U_NonConst()->ElementWiseMultiply(*tmp);
        rhs->v_U_NonConst()->AddScalar(mu);

        DBG_PRINT_VECTOR(2, "rhs", *rhs);

        pd_solver_->Solve(1.0, 1.0, *rhs, *delta_corr, true);

        DBG_PRINT_VECTOR(2, "delta_corr", *delta_corr);
      }
      break;
      default:
      DBG_ASSERT(false && "Unknown corrector_type value.");
    }

    // Compute step size
    Number alpha_primal_corr =
      IpCq().primal_frac_to_the_bound(IpData().curr_tau(),
                                      *delta_corr->x(),
                                      *delta_corr->s());
    // Set the primal trial point
    IpData().SetTrialPrimalVariablesFromStep(alpha_primal_corr, *delta_corr->x(), *delta_corr->s());

    // Check if we want to not even try the filter criterion
    Number alpha_dual_max =
      IpCq().dual_frac_to_the_bound(IpData().curr_tau(),
                                    *delta_corr->z_L(), *delta_corr->z_U(),
                                    *delta_corr->v_L(), *delta_corr->v_U());

    IpData().SetTrialBoundMultipliersFromStep(alpha_dual_max, *delta_corr->z_L(), *delta_corr->z_U(), *delta_corr->v_L(), *delta_corr->v_U());

    Number trial_avrg_compl = IpCq().trial_avrg_compl();
    Number curr_avrg_compl = IpCq().curr_avrg_compl();
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "avrg_compl(curr) = %e, avrg_compl(trial) = %e\n",
                   curr_avrg_compl, trial_avrg_compl);
    if (corrector_type_==AFFINE_CORRECTOR &&
        trial_avrg_compl>=corrector_compl_avrg_red_fact_*curr_avrg_compl) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Rejecting corrector step, because trial complementarity is too large.\n" );
      IpData().TimingStats().TryCorrector().End();
      return false;
    }

    // Check if trial point is acceptable
    try {
      // in acceptance tests, use original step size!
      accept = CheckAcceptabilityOfTrialPoint(alpha_primal_test);
    }
    catch(IpoptNLP::Eval_Error& e) {
      e.ReportException(Jnlst(), J_DETAILED);
      Jnlst().Printf(J_WARNING, J_MAIN,
                     "Warning: Corrector step rejected due to evaluation error\n");
      IpData().Append_info_string("e");
      accept = false;
    }

    if (accept) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Corrector step accepted with alpha_primal = %e\n",
                     alpha_primal_corr);
      // Accept all SOC quantities
      alpha_primal = alpha_primal_corr;
      actual_delta = delta_corr;

      if (Jnlst().ProduceOutput(J_MOREVECTOR, J_MAIN)) {
        Jnlst().Printf(J_MOREVECTOR, J_MAIN,
                       "*** Accepted corrector for Iteration: %d\n",
                       IpData().iter_count());
        delta_corr->Print(Jnlst(), J_MOREVECTOR, J_MAIN, "delta_corr");
      }
    }

    IpData().TimingStats().TryCorrector().End();
    return accept;
  }

  char FilterLSAcceptor::UpdateForNextIteration(Number alpha_primal_test)
  {
    char info_alpha_primal_char;
    // Augment the filter if required
    if (!IsFtype(alpha_primal_test) ||
        !ArmijoHolds(alpha_primal_test)) {
      AugmentFilter();
      info_alpha_primal_char = 'h';
    }
    else {
      info_alpha_primal_char = 'f';
    }
    return info_alpha_primal_char;
  }

  void FilterLSAcceptor::PrepareRestoPhaseStart()
  {
    AugmentFilter();
  }


} // namespace Ipopt
