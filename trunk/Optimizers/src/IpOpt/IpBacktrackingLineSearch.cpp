// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpBacktrackingLineSearch.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13
//           Andreas Waechter                 IBM    2005-10-13
//               derived file from IpFilterLineSearch.cpp

#include "IpBacktrackingLineSearch.hpp"
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

namespace Ipopt
{

#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  BacktrackingLineSearch::BacktrackingLineSearch(
    const SmartPtr<BacktrackingLSAcceptor>& acceptor,
    const SmartPtr<RestorationPhase>& resto_phase,
    const SmartPtr<ConvergenceCheck>& conv_check)
      :
      LineSearch(),
      acceptor_(acceptor),
      resto_phase_(resto_phase),
      conv_check_(conv_check)
  {
    DBG_START_FUN("BacktrackingLineSearch::BacktrackingLineSearch",
                  dbg_verbosity);
    DBG_ASSERT(IsValid(acceptor_));
  }

  BacktrackingLineSearch::~BacktrackingLineSearch()
  {
    DBG_START_FUN("BacktrackingLineSearch::~BacktrackingLineSearch()",
                  dbg_verbosity);
  }

  void BacktrackingLineSearch::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddBoundedNumberOption(
      "alpha_red_factor",
      "Fractional reduction of the trial step size in the backtracking line search.",
      0.0, true, 1.0, true, 0.5,
      "At every step of the backtracking line search, the trial step size is "
      "reduced by this factor.");

    std::string prev_category = roptions->RegisteringCategory();
    roptions->SetRegisteringCategory("Undocumented");
    roptions->AddStringOption2(
      "magic_steps",
      "Enables magic steps.",
      "no",
      "no", "don't take magic steps",
      "yes", "take magic steps",
      "DOESN'T REALLY WORK YET!");
    roptions->SetRegisteringCategory(prev_category);

    roptions->AddStringOption2(
      "accept_every_trial_step",
      "Always accept the frist trial step.",
      "no",
      "no", "don't arbitrarily accept the full step",
      "yes", "always accept the full step",
      "Setting this option to \"yes\" essentially disables the line search "
      "and makes the algorithm take aggressive steps, without global "
      "convergence guarantees.");

    roptions->AddStringOption7(
      "alpha_for_y",
      "Method to determine the step size for constraint multipliers.",
      "primal",
      "primal", "use primal step size",
      "bound_mult", "use step size for the bound multipliers (good for LPs)",
      "min", "use the min of primal and bound multipliers",
      "max", "use the max of primal and bound multipliers",
      "full", "take a full step of size one",
      "min_dual_infeas", "choose step size minimizing new dual infeasibility",
      "safe_min_dual_infeas", "like \"min_dual_infeas\", but safeguarded by \"min\" and \"max\"",
      "This option determines how the step size (alpha_y) will be calculated when updating the "
      "constraint multipliers.");

    roptions->AddLowerBoundedNumberOption(
      "tiny_step_tol",
      "Tolerance for detecting numerically insignificant steps.",
      0.0, false, 10.0*std::numeric_limits<double>::epsilon(),
      "If the search direction in the primal variables (x and s) is, in "
      "relative terms for each component, less than this value, the "
      "algorithm accepts the full step without line search.  If this happens "
      "repeatedly, the algorithm will terminate with a corresponding exit "
      "message. The default value is 10 times machine precision.");
    roptions->AddLowerBoundedNumberOption(
      "tiny_step_y_tol",
      "Tolerance for quitting because of numerically insignificant steps.",
      0.0, false, 1e-2,
      "If the search direction in the primal variables (x and s) is, in "
      "relative terms for each component, repeatedly less than tiny_step_tol, "
      "and the step in the y variables is smaller than this threshold, the "
      "algorithm will terminate.");
    roptions->AddLowerBoundedIntegerOption(
      "watchdog_shortened_iter_trigger",
      "Number of shortened iterations that trigger the watchdog.",
      0, 10,
      "If the number of successive iterations in which the backtracking line search "
      "did not accept the first trial point exceeds this number, the "
      "watchdog procedure is activated.  Choosing \"0\" here disables the "
      "watchdog procedure.");
    roptions->AddLowerBoundedIntegerOption(
      "watchdog_trial_iter_max",
      "Maximum number of watchdog iterations.",
      1, 3,
      "This option determines the number of trial iterations "
      "allowed before the watchdog "
      "procedure is aborted and the algorithm returns to the stored point.");

    roptions->SetRegisteringCategory("Restoration Phase");
    roptions->AddStringOption2(
      "expect_infeasible_problem",
      "Enable heuristics to quickly detect an infeasible problem.",
      "no",
      "no", "the problem probably be feasible",
      "yes", "the problem has a good chance to be infeasible",
      "This options is meant to activate heuristics that may speed up the "
      "infeasibility determination if you expect that there is a good chance for the problem to be "
      "infeasible.  In the filter line search procedure, the restoration "
      "phase is called more quickly than usually, and more reduction in "
      "the constraint violation is enforced before the restoration phase is "
      "left. If the problem is square, this option is enabled automatically.");
    roptions->AddLowerBoundedNumberOption(
      "expect_infeasible_problem_ctol",
      "Threshold for disabling \"expect_infeasible_problem\" option.",
      0.0, false, 1e-3,
      "If the constraint violation becomes smaller than this threshold, "
      "the \"expect_infeasible_problem\" heuristics in the filter line "
      "search are disabled. If the problem is square, this options is set to "
      "0.");
    roptions->AddStringOption2(
      "start_with_resto",
      "Tells algorithm to switch to restoration phase in first iteration.",
      "no",
      "no", "don't force start in restoration phase",
      "yes", "force start in restoration phase",
      "Setting this option to \"yes\" forces the algorithm to switch to the "
      "feasibility restoration phase in the first iteration. If the initial "
      "point is feasible, the algorithm will abort with a failure.");
    roptions->AddLowerBoundedNumberOption(
      "soft_resto_pderror_reduction_factor",
      "Required reduction in primal-dual error in the soft restoration phase.",
      0.0, false, (1.0 - 1e-4),
      "The soft restoration phase attempts to reduce the "
      "primal-dual error with regular steps. If the damped "
      "primal-dual step (damped only to satisfy the "
      "fraction-to-the-boundary rule) is not decreasing the primal-dual error "
      "by at least this factor, then the regular restoration phase is called. "
      "Choosing \"0\" here disables the soft "
      "restoration phase.");
    roptions->AddLowerBoundedIntegerOption(
      "max_soft_resto_iters",
      "Maximum number of iterations performed successively in soft restoration phase.",
      0, 10,
      "If the soft restoration phase is performed for more than so many "
      "iteratins in a row, the regular restoration phase is called.");
  }

  bool BacktrackingLineSearch::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("alpha_red_factor", alpha_red_factor_, prefix);
    options.GetBoolValue("magic_steps", magic_steps_, prefix);
    options.GetBoolValue("accept_every_trial_step", accept_every_trial_step_, prefix);
    Index enum_int;
    options.GetEnumValue("alpha_for_y", enum_int, prefix);
    alpha_for_y_ = AlphaForYEnum(enum_int);
    options.GetNumericValue("expect_infeasible_problem_ctol", expect_infeasible_problem_ctol_, prefix);
    options.GetBoolValue("expect_infeasible_problem", expect_infeasible_problem_, prefix);

    options.GetBoolValue("start_with_resto", start_with_resto_, prefix);

    options.GetNumericValue("tiny_step_tol", tiny_step_tol_, prefix);
    options.GetNumericValue("tiny_step_y_tol", tiny_step_y_tol_, prefix);
    options.GetIntegerValue("watchdog_trial_iter_max", watchdog_trial_iter_max_, prefix);
    options.GetIntegerValue("watchdog_shortened_iter_trigger", watchdog_shortened_iter_trigger_, prefix);
    options.GetNumericValue("soft_resto_pderror_reduction_factor",
                            soft_resto_pderror_reduction_factor_, prefix);
    options.GetIntegerValue("max_soft_resto_iters", max_soft_resto_iters_,
                            prefix);

    bool retvalue = true;
    if (IsValid(resto_phase_)) {
      if (!resto_phase_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                    options, prefix)) {
        return false;
      }
      if (!acceptor_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                 options, prefix)) {
        return false;
      }
    }

    rigorous_ = true;
    skipped_line_search_ = false;
    tiny_step_last_iteration_ = false;
    fallback_activated_ = false;

    Reset();

    count_successive_shortened_steps_ = 0;

    acceptable_iterate_ = NULL;
    acceptable_iteration_number_ = -1;

    last_mu_ = -1.;

    return retvalue;
  }

  void BacktrackingLineSearch::FindAcceptableTrialPoint()
  {
    DBG_START_METH("BacktrackingLineSearch::FindAcceptableTrialPoint",
                   dbg_verbosity);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "--> Starting filter line search in iteration %d <--\n",
                   IpData().iter_count());

    Number curr_mu = IpData().curr_mu();
    if (last_mu_!=curr_mu) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Mu has changed in line search - resetting watchdog counters.\n");
      // Inactivate the watchdog and release all stored data
      in_watchdog_ = false;
      watchdog_iterate_ = NULL;
      watchdog_delta_ = NULL;
      watchdog_shortened_iter_ = 0;
      last_mu_ = curr_mu;
    }

    // If the problem is square, we want to enable the
    // expect_infeasible_problem option automatically so that the
    // restoration phase is entered soon
    if (IpCq().IsSquareProblem()) {
      expect_infeasible_problem_ = true;
      expect_infeasible_problem_ctol_ = 0.;
    }

    // Store current iterate if the optimality error is on acceptable
    // level to restored if things fail later
    if (CurrentIsAcceptable()) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Storing current iterate as backup acceptable point.\n");
      StoreAcceptablePoint();
    }

    // First assume that line search will find an acceptable trial point
    skipped_line_search_ = false;

    // Get the search directions (this will store the actual search
    // direction, possibly including higher order corrections)
    SmartPtr<IteratesVector> actual_delta;

    bool goto_resto = false;
    if (fallback_activated_) {
      // In this case, the algorithm had trouble to continue and wants
      // to call the restoration phase immediately
      goto_resto = true;
      fallback_activated_ = false; // reset the flag
    }
    else {
      // Initialize the acceptor for this backtracking line search
      acceptor_->InitThisLineSearch(in_watchdog_);
      actual_delta = IpData().delta()->MakeNewContainer();
    }

    if (start_with_resto_) {
      // If the user requested to start with the restoration phase,
      // skip the line search and do exactly that.  Reset the flag so
      // that this happens only once.
      goto_resto = true;
      start_with_resto_= false;
    }

    if (expect_infeasible_problem_ &&
        Max(IpData().curr()->y_c()->Amax(),
            IpData().curr()->y_d()->Amax()) > 1e8) {
      goto_resto = true;
    }

    bool accept = false;
    bool corr_taken = false;
    bool soc_taken = false;
    Index n_steps = 0;
    Number alpha_primal = 0.;

    // Check if search direction becomes too small
    bool tiny_step = (!goto_resto && DetectTinyStep());

    if (in_watchdog_ && (goto_resto || tiny_step)) {
      // If the step could not be computed or is too small and the
      // watchdog is active, stop the watch dog and resume everything
      // from reference point
      StopWatchDog(actual_delta);
      goto_resto = false;
      tiny_step = false;
    }

    // Check if we want to wake up the watchdog
    if (watchdog_shortened_iter_trigger_ > 0 &&
        !in_watchdog_ && !goto_resto && !tiny_step &&
        !in_soft_resto_phase_ && !expect_infeasible_problem_ &&
        watchdog_shortened_iter_ >= watchdog_shortened_iter_trigger_) {
      StartWatchDog();
    }

    // Handle the situation of a tiny step
    if (tiny_step) {
      alpha_primal =
        IpCq().curr_primal_frac_to_the_bound(IpData().curr_tau());
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "Tiny step detected. Use step size alpha = %e unchecked\n",
                     alpha_primal);
      IpData().SetTrialPrimalVariablesFromStep(alpha_primal, *IpData().delta()->x(), *IpData().delta()->s());

      // Evaluate functions at trial point - if that fails, don't use
      // the tiny step and continue with regular line search
      try {
        IpCq().trial_barrier_obj();
        IpCq().trial_constraint_violation();
      }
      catch(IpoptNLP::Eval_Error& e) {
        e.ReportException(Jnlst(), J_DETAILED);
        tiny_step = false;
      }

      if (tiny_step) {
        IpData().Set_info_ls_count(0);

        if (tiny_step_last_iteration_) {
          IpData().Set_info_alpha_primal_char('T');
          IpData().Set_tiny_step_flag(true);
        }
      }
      else {
        IpData().Set_info_alpha_primal_char('t');
      }

      // If the step in the dual variables is also small, we remember
      // that we just did a tiny step so that next time we might
      // decide to quit
      Number delta_y_norm = Max(IpData().curr()->y_c()->Amax(),
                                IpData().curr()->y_d()->Amax());
      if (delta_y_norm < tiny_step_y_tol_) {
        tiny_step_last_iteration_ = true;
      }
      else {
        tiny_step_last_iteration_ = false;
      }
      accept = true;
    }
    else {
      tiny_step_last_iteration_ = false;
    }

    if (!goto_resto && !tiny_step) {

      if (in_soft_resto_phase_) {
        soft_resto_counter_++;
        if (soft_resto_counter_ > max_soft_resto_iters_) {
          accept = false;
        }
        else {
          // If we are currently in the soft restoration phase, continue
          // that way, and switch back if enough progress is made to the
          // original criterion (e.g., the filter)
          bool satisfies_original_criterion = false;
          // ToDo use tiny_step in TrySoftRestoStep?
          accept = TrySoftRestoStep(actual_delta,
                                    satisfies_original_criterion);
          if (accept) {
            IpData().Set_info_alpha_primal_char('s');
            if (satisfies_original_criterion) {
              in_soft_resto_phase_ = false;
              soft_resto_counter_ = 0;
              IpData().Set_info_alpha_primal_char('S');
            }
          }
        }
      }
      else {
        // Start the backtracking line search
        bool done = false;
        bool skip_first_trial_point = false;
        bool evaluation_error;
        while (!done) {
          accept = DoBacktrackingLineSearch(skip_first_trial_point,
                                            alpha_primal,
                                            corr_taken,
                                            soc_taken,
                                            n_steps,
                                            evaluation_error,
                                            actual_delta);
          DBG_PRINT((1, "evaluation_error = %d\n", evaluation_error));
          if (in_watchdog_) {
            if (accept) {
              in_watchdog_ = false;
              IpData().Append_info_string("W");
              Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                             "Watch dog procedure successful!\n");
              done = true;
            }
            else {
              watchdog_trial_iter_++;
              if (evaluation_error ||
                  watchdog_trial_iter_ > watchdog_trial_iter_max_) {
                StopWatchDog(actual_delta);
                skip_first_trial_point = true;
              }
              else {
                done = true;
                accept = true;
              }
            }
          }
          else {
            done = true;
          }
        }
      } /* else: if (in_soft_resto_phase_) { */
    } /* if (!goto_resto && !tiny_step) { */

    // If line search has been aborted because the step size becomes
    // too small, go to the restoration phase or continue with soft
    // restoration phase
    if (!accept) {
      // If we are not asked to do a rigorous line search, do no call
      // the restoration phase.
      if (!rigorous_) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Skipping call of restoration phase...\n");
        skipped_line_search_=true;
      }
      else {
        // Check if we should start the soft restoration phase
        if (!in_soft_resto_phase_
            && !goto_resto && !expect_infeasible_problem_) {
          Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                         "--> Starting soft restoration phase <--\n");
          // Prepare the restoration phase, e.g., augment the filter
          // with the current point.
          acceptor_->PrepareRestoPhaseStart();

          // Try the current search direction for the soft restoration phase
          bool satisfies_original_criterion;
          accept = TrySoftRestoStep(actual_delta,
                                    satisfies_original_criterion);
          // If it has been accepted: If the original criterion is also
          // satisfied, we can just take that step and continue with
          // the regular algorithm, otherwise we stay in the soft
          // restoration phase
          if (accept) {
            if (satisfies_original_criterion) {
              IpData().Set_info_alpha_primal_char('S');
            }
            else {
              in_soft_resto_phase_ = true;
              IpData().Set_info_alpha_primal_char('s');
            }
          }
        }

        if (!accept) {
          // Go to the restoration phase
          if (!in_soft_resto_phase_) {
            // Prepare the restoration phase, e.g., augment the filter
            // with the current point. If we are already in the soft
            // restoration phase, this has been done earlier
            acceptor_->PrepareRestoPhaseStart();
          }
          if (CurrentIsAcceptable()) {
            THROW_EXCEPTION(ACCEPTABLE_POINT_REACHED,
                            "Restoration phase called at acceptable point.");
          }

          if (!IsValid(resto_phase_)) {
            THROW_EXCEPTION(IpoptException, "No Restoration Phase given to this Filter Line Search Object!");
          }
          // ToDo make the 1e-2 below a parameter?
          if (IpCq().curr_constraint_violation()<=
              1e-2*IpData().tol()) {
            bool found_acceptable = RestoreAcceptablePoint();
            if (found_acceptable) {
              Jnlst().Printf(J_WARNING, J_LINE_SEARCH,
                             "Restoration phase is called at almost feasible point,\n  but acceptable point from iteration %d could be restored.\n", acceptable_iteration_number_);
              THROW_EXCEPTION(ACCEPTABLE_POINT_REACHED,
                              "Restoration phase called at almost feasible point, but acceptable point could be restored.\n");
            }
            else {
              // ToDo does that happen too often?
              Jnlst().Printf(J_ERROR, J_LINE_SEARCH,
                             "Restoration phase is called at point that is almost feasible,\n  with constraint violation %e. Abort.\n", IpCq().curr_constraint_violation());
              THROW_EXCEPTION(RESTORATION_FAILED,
                              "Restoration phase called, but point is almost feasible.");
            }
          }

          // Set the info fields for the first output line in the
          // restoration phase which reflects why the restoration phase
          // was called
          IpData().Set_info_alpha_primal(alpha_primal);
          IpData().Set_info_alpha_dual(0.);
          IpData().Set_info_alpha_primal_char('R');
          IpData().Set_info_ls_count(n_steps+1);

          accept = resto_phase_->PerformRestoration();
          if (!accept) {
            bool found_acceptable = RestoreAcceptablePoint();
            if (found_acceptable) {
              THROW_EXCEPTION(ACCEPTABLE_POINT_REACHED,
                              "Restoration phase failed, but acceptable point could be restore.\n");
            }
            else {
              THROW_EXCEPTION(RESTORATION_FAILED,
                              "Failed restoration phase!!!");
            }
          }
          count_successive_shortened_steps_ = 0;
          if (expect_infeasible_problem_) {
            expect_infeasible_problem_ = false;
          }
          in_soft_resto_phase_ = false;
          soft_resto_counter_ = 0;
          watchdog_shortened_iter_ = 0;
        }
      }
    }
    else if (!in_soft_resto_phase_ || tiny_step) {
      // we didn't do the restoration phase and are now updating the
      // dual variables of the trial point
      Number alpha_dual_max =
        IpCq().dual_frac_to_the_bound(IpData().curr_tau(),
                                      *actual_delta->z_L(), *actual_delta->z_U(),
                                      *actual_delta->v_L(), *actual_delta->v_U());

      PerformDualStep(alpha_primal, alpha_dual_max, actual_delta);

      if (n_steps==0) {
        // accepted this if a full step was
        // taken
        count_successive_shortened_steps_ = 0;
        watchdog_shortened_iter_ = 0;
      }
      else {
        count_successive_shortened_steps_++;
        watchdog_shortened_iter_++;
      }

      if (expect_infeasible_problem_ &&
          IpCq().curr_constraint_violation() <= expect_infeasible_problem_ctol_) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "Constraint violation is with %e less than expect_infeasible_problem_ctol.\nDisable expect_infeasible_problem_heuristic.\n", IpCq().curr_constraint_violation());
        expect_infeasible_problem_ = false;
      }
    }
  }

  bool BacktrackingLineSearch::DoBacktrackingLineSearch(bool skip_first_trial_point,
      Number& alpha_primal,
      bool& corr_taken,
      bool& soc_taken,
      Index& n_steps,
      bool& evaluation_error,
      SmartPtr<IteratesVector>& actual_delta)
  {
    evaluation_error = false;
    bool accept = false;

    DBG_START_METH("BacktrackingLineSearch::DoBacktrackingLineSearch",
                   dbg_verbosity);

    // Compute primal fraction-to-the-boundary value
    Number alpha_primal_max =
      IpCq().primal_frac_to_the_bound(IpData().curr_tau(),
                                      *actual_delta->x(),
                                      *actual_delta->s());

    // Compute smallest step size allowed
    Number alpha_min = alpha_primal_max;
    if (!in_watchdog_) {
      alpha_min = acceptor_->CalculateAlphaMin();
    }
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "minimal step size ALPHA_MIN = %E\n", alpha_min);

    // Start line search from maximal step size
    alpha_primal = alpha_primal_max;

    // Step size used in ftype and armijo tests
    Number alpha_primal_test = alpha_primal;
    if (in_watchdog_) {
      alpha_primal_test = watchdog_alpha_primal_test_;
    }

    if (skip_first_trial_point) {
      alpha_primal *= alpha_red_factor_;
    }

    if (!skip_first_trial_point) {
      // Before we do the actual backtracking line search for the
      // regular primal-dual search direction, let's see if a step
      // including a higher-order correctior is already acceptable
      accept = acceptor_->TryCorrector(alpha_primal_test,
                                       alpha_primal,
                                       actual_delta);
    }
    if (accept) {
      corr_taken = true;
    }

    if (!accept) {
      // Loop over decreaseing step sizes until acceptable point is
      // found or until step size becomes too small

      while (alpha_primal>alpha_min ||
             n_steps == 0) { // always allow the "full" step if it is
        // acceptable (even if alpha_primal<=alpha_min)
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "Starting checks for alpha (primal) = %8.2e\n",
                       alpha_primal);

        try {
          // Compute the primal trial point
          IpData().SetTrialPrimalVariablesFromStep(alpha_primal, *actual_delta->x(), *actual_delta->s());

          if (magic_steps_) {
            PerformMagicStep();
          }

          // If it is acceptable, stop the search
          alpha_primal_test = alpha_primal;
          if (accept_every_trial_step_) {
            // We call the evaluation at the trial point here, so that an
            // exception will the thrown if there are problem during the
            // evaluation of the functions (in that case, we want to further
            // reduce the step size
            IpCq().trial_barrier_obj();
            IpCq().trial_constraint_violation();
            accept = true;
          }
          else {
            accept = acceptor_->CheckAcceptabilityOfTrialPoint(alpha_primal_test);
          }
        }
        catch(IpoptNLP::Eval_Error& e) {
          e.ReportException(Jnlst(), J_DETAILED);
          Jnlst().Printf(J_WARNING, J_LINE_SEARCH,
                         "Warning: Cutting back alpha due to evaluation error\n");
          IpData().Append_info_string("e");
          accept = false;
          evaluation_error = true;
        }

        if (accept) {
          break;
        }

        if (in_watchdog_) {
          break;
        }

        // Decide if we want to go to the restoration phase in a
        // short cut to check if the problem is infeasible
        if (expect_infeasible_problem_) {
          if (count_successive_shortened_steps_>=5) {
            break;
          }
        }

        // try second order correction step if the function could
        // be evaluated
        // DoTo: check if we want to do SOC when watchdog is active
        if (!evaluation_error) {
          Number theta_curr = IpCq().curr_constraint_violation();
          Number theta_trial = IpCq().trial_constraint_violation();
          if (alpha_primal==alpha_primal_max &&       // i.e. first trial point
              theta_curr<=theta_trial) {
            // Try second order correction
            accept = acceptor_->TrySecondOrderCorrection(alpha_primal_test,
                     alpha_primal,
                     actual_delta);
          }
          if (accept) {
            soc_taken = true;
            break;
          }
        }

        // Point is not yet acceptable, try a shorter one
        alpha_primal *= alpha_red_factor_;
        n_steps++;
      }
    } /* if (!accept) */

    char info_alpha_primal_char='?';
    if (!accept && in_watchdog_) {
      info_alpha_primal_char = 'w';
    }
    else if (accept) {
      info_alpha_primal_char =
        acceptor_->UpdateForNextIteration(alpha_primal_test);
    }
    if (soc_taken) {
      info_alpha_primal_char = toupper(info_alpha_primal_char);
    }
    IpData().Set_info_alpha_primal_char(info_alpha_primal_char);
    IpData().Set_info_ls_count(n_steps+1);
    if (corr_taken) {
      IpData().Append_info_string("C");
    }

    return accept;
  }

  void BacktrackingLineSearch::StartWatchDog()
  {
    DBG_START_FUN("BacktrackingLineSearch::StartWatchDog", dbg_verbosity);

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Starting Watch Dog\n");

    in_watchdog_ = true;
    watchdog_iterate_ = IpData().curr();
    watchdog_delta_ = IpData().delta();
    watchdog_trial_iter_ = 0;
    watchdog_alpha_primal_test_ =
      IpCq().curr_primal_frac_to_the_bound(IpData().curr_tau());

    acceptor_->StartWatchDog();
  }

  void BacktrackingLineSearch::StopWatchDog(SmartPtr<IteratesVector>& actual_delta)
  {
    DBG_START_FUN("BacktrackingLineSearch::StopWatchDog", dbg_verbosity);

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Stopping Watch Dog\n");

    IpData().Append_info_string("w");

    in_watchdog_ = false;

    // Reset all fields in IpData to reference point
    SmartPtr<IteratesVector> old_trial = watchdog_iterate_->MakeNewContainer();
    IpData().set_trial(old_trial);
    IpData().AcceptTrialPoint();
    actual_delta = watchdog_delta_->MakeNewContainer();
    IpData().SetHaveAffineDeltas(false);

    // reset the stored watchdog iterates
    watchdog_iterate_ = NULL;
    watchdog_delta_ = NULL;

    watchdog_shortened_iter_ = 0;

    acceptor_->StopWatchDog();
  }

  void BacktrackingLineSearch::Reset()
  {
    DBG_START_FUN("BacktrackingLineSearch::Reset", dbg_verbosity);
    in_soft_resto_phase_ = false;
    soft_resto_counter_ = 0;

    acceptor_->Reset();
  }

  void BacktrackingLineSearch::PerformDualStep(Number alpha_primal,
      Number alpha_dual,
      SmartPtr<IteratesVector>& delta)
  {
    DBG_START_FUN("BacktrackingLineSearch::PerformDualStep", dbg_verbosity);

    // set the bound multipliers from the step
    IpData().SetTrialBoundMultipliersFromStep(alpha_dual, *delta->z_L(), *delta->z_U(), *delta->v_L(), *delta->v_U());

    Number alpha_y=-1.;
    switch (alpha_for_y_) {
      case PRIMAL_ALPHA_FOR_Y:
      alpha_y = alpha_primal;
      break;
      case DUAL_ALPHA_FOR_Y:
      alpha_y = alpha_dual;
      break;
      case MIN_ALPHA_FOR_Y:
      alpha_y = Min(alpha_dual, alpha_primal);
      break;
      case MAX_ALPHA_FOR_Y:
      alpha_y = Max(alpha_dual, alpha_primal);
      break;
      case FULL_STEP_FOR_Y:
      alpha_y = 1;
      break;
      case MIN_DUAL_INFEAS_ALPHA_FOR_Y:
      case SAFE_MIN_DUAL_INFEAS_ALPHA_FOR_Y:
      // Here we compute the step size for y so that the dual
      // infeasibility is minimized along delta_y

      // compute the dual infeasibility at new point with old y
      SmartPtr<IteratesVector> temp_trial
      = IpData().trial()->MakeNewContainer();
      temp_trial->Set_y_c(*IpData().curr()->y_c());
      temp_trial->Set_y_d(*IpData().curr()->y_d());
      IpData().set_trial(temp_trial);
      SmartPtr<const Vector> dual_inf_x = IpCq().trial_grad_lag_x();
      SmartPtr<const Vector> dual_inf_s = IpCq().trial_grad_lag_s();

      SmartPtr<Vector> new_jac_times_delta_y =
        IpData().curr()->x()->MakeNew();
      new_jac_times_delta_y->AddTwoVectors(1., *IpCq().trial_jac_cT_times_vec(*delta->y_c()),
                                           1., *IpCq().trial_jac_dT_times_vec(*delta->y_d()),
                                           0.);

      Number a = pow(new_jac_times_delta_y->Nrm2(), 2.) +
                 pow(delta->y_d()->Nrm2(), 2.);
      Number b = dual_inf_x->Dot(*new_jac_times_delta_y) -
                 dual_inf_s->Dot(*delta->y_d());

      Number alpha = - b/a;

      if (alpha_for_y_==SAFE_MIN_DUAL_INFEAS_ALPHA_FOR_Y) {
        alpha_y = Min(Max(alpha_primal, alpha_dual), Max(alpha, Min(alpha_primal, alpha_dual)));
      }
      else {
        alpha_y = Min(1., Max(0., alpha));
      }
      break;
    }

    // Set the eq multipliers from the step now that alpha_y
    // has been calculated.
    DBG_PRINT((1, "alpha_y = %e\n", alpha_y));
    DBG_PRINT_VECTOR(2, "delta_y_c", *delta->y_c());
    DBG_PRINT_VECTOR(2, "delta_y_d", *delta->y_d());
    IpData().SetTrialEqMultipliersFromStep(alpha_y, *delta->y_c(), *delta->y_d());

    // Set some information for iteration summary output
    IpData().Set_info_alpha_primal(alpha_primal);
    IpData().Set_info_alpha_dual(alpha_dual);
  }

  void
  BacktrackingLineSearch::PerformMagicStep()
  {
    DBG_START_METH("BacktrackingLineSearch::PerformMagicStep",
                   2);//dbg_verbosity);

    DBG_PRINT((1,"Incoming barr = %e and constrviol %e\n",
               IpCq().trial_barrier_obj(),
               IpCq().trial_constraint_violation()));
    DBG_PRINT_VECTOR(2, "s in", *IpData().trial()->s());
    DBG_PRINT_VECTOR(2, "d minus s in", *IpCq().trial_d_minus_s());
    DBG_PRINT_VECTOR(2, "slack_s_L in", *IpCq().trial_slack_s_L());
    DBG_PRINT_VECTOR(2, "slack_s_U in", *IpCq().trial_slack_s_U());

    SmartPtr<const Vector> d_L = IpNLP().d_L();
    SmartPtr<const Matrix> Pd_L = IpNLP().Pd_L();
    SmartPtr<Vector> delta_s_magic_L = d_L->MakeNew();
    delta_s_magic_L->Set(0.);
    SmartPtr<Vector> tmp = d_L->MakeNew();
    Pd_L->TransMultVector(1., *IpCq().trial_d_minus_s(), 0., *tmp);
    delta_s_magic_L->ElementWiseMax(*tmp);

    SmartPtr<const Vector> d_U = IpNLP().d_U();
    SmartPtr<const Matrix> Pd_U = IpNLP().Pd_U();
    SmartPtr<Vector> delta_s_magic_U = d_U->MakeNew();
    delta_s_magic_U->Set(0.);
    tmp = d_U->MakeNew();
    Pd_U->TransMultVector(1., *IpCq().trial_d_minus_s(), 0., *tmp);
    delta_s_magic_U->ElementWiseMin(*tmp);

    SmartPtr<Vector> delta_s_magic = IpData().trial()->s()->MakeNew();
    Pd_L->MultVector(1., *delta_s_magic_L, 0., *delta_s_magic);
    Pd_U->MultVector(1., *delta_s_magic_U, 1., *delta_s_magic);
    delta_s_magic_L = NULL; // free memory
    delta_s_magic_U = NULL; // free memory

    // Now find those entries with both lower and upper bounds, there
    // the step is too large
    // ToDo this should only be done if there are inequality
    // constraints with two bounds
    // also this can be done in a smaller space (d_L or d_U whichever
    // is smaller)
    tmp = delta_s_magic->MakeNew();
    tmp->Copy(*IpData().trial()->s());
    Pd_L->MultVector(1., *d_L, -2., *tmp);
    Pd_U->MultVector(1., *d_U, 1., *tmp);
    SmartPtr<Vector> tmp2 = tmp->MakeNew();
    tmp2->Copy(*tmp);
    tmp2->ElementWiseAbs();
    tmp->Axpy(-2., *delta_s_magic);
    tmp->ElementWiseAbs();
    // now, tmp2 = |d_L + d_u - 2*s| and tmp = |d_L + d_u - 2*(s+Delta s)|
    // we want to throw out those for which tmp2 > tmp
    tmp->Axpy(-1., *tmp2);
    tmp->ElementWiseSgn();
    tmp2->Set(0.);
    tmp2->ElementWiseMax(*tmp);
    tmp = d_L->MakeNew();
    Pd_L->TransMultVector(1., *tmp2, 0., *tmp);
    Pd_L->MultVector(1., *tmp, 0., *tmp2);
    tmp = d_U->MakeNew();
    Pd_U->TransMultVector(1., *tmp2, 0., *tmp);
    Pd_U->MultVector(1., *tmp, 0., *tmp2);
    DBG_PRINT_VECTOR(2, "tmp indicator", *tmp2)
    // tmp2 now is one for those entries with both bounds, for which
    // no step should be taken

    tmp = delta_s_magic->MakeNew();
    tmp->Copy(*delta_s_magic);
    tmp->ElementWiseMultiply(*tmp2);
    delta_s_magic->Axpy(-1., *tmp);

    Number delta_s_magic_max = delta_s_magic->Amax();
    Number mach_eps = std::numeric_limits<Number>::epsilon();
    if (delta_s_magic_max>0.) {
      if (delta_s_magic_max > 10*mach_eps*IpData().trial()->s()->Amax()) {
        IpData().Append_info_string("M");
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH, "Magic step with max-norm %.6e taken.\n", delta_s_magic->Amax());
        delta_s_magic->Print(Jnlst(), J_MOREVECTOR, J_LINE_SEARCH,
                             "delta_s_magic");
      }

      // now finally compute the new overall slacks
      delta_s_magic->Axpy(1., *IpData().trial()->s());
      SmartPtr<IteratesVector> trial = IpData().trial()->MakeNewContainer();
      trial->Set_s(*delta_s_magic);

      // also update the set in the dual variables


      IpData().set_trial(trial);
    }

    DBG_PRINT((1,"Outgoing barr = %e and constrviol %e\n", IpCq().trial_barrier_obj(), IpCq().trial_constraint_violation()));
    DBG_PRINT_VECTOR(2, "s out", *IpData().trial()->s());
    DBG_PRINT_VECTOR(2, "d minus s out", *IpCq().trial_d_minus_s());
    DBG_PRINT_VECTOR(2, "slack_s_L out", *IpCq().trial_slack_s_L());
    DBG_PRINT_VECTOR(2, "slack_s_U out", *IpCq().trial_slack_s_U());
  }

  bool BacktrackingLineSearch::TrySoftRestoStep(SmartPtr<IteratesVector>& actual_delta,
      bool &satisfies_original_criterion)
  {
    DBG_START_FUN("FilterLSAcceptor::TrySoftRestoStep", dbg_verbosity);

    if (soft_resto_pderror_reduction_factor_==0.) {
      return false;
    }

    satisfies_original_criterion = false;

    // ToDo: Need to decide if we want to try a corrector step first

    // Compute the maximal step sizes (we use identical step sizes for
    // primal and dual variables
    Number alpha_primal_max =
      IpCq().primal_frac_to_the_bound(IpData().curr_tau(),
                                      *actual_delta->x(),
                                      *actual_delta->s());
    Number alpha_dual_max =
      IpCq().dual_frac_to_the_bound(IpData().curr_tau(),
                                    *actual_delta->z_L(),
                                    *actual_delta->z_U(),
                                    *actual_delta->v_L(),
                                    *actual_delta->v_U());
    Number alpha =  Min(alpha_primal_max, alpha_dual_max);

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Trying soft restoration phase step with step length %13.6e\n",
                   alpha);

    // We allow up to three trials in case there is an evaluation
    // error for the functions
    bool done=false;
    Index count=3;
    while (!done && count>0) {
      // Set the trial point
      IpData().SetTrialPrimalVariablesFromStep(alpha, *actual_delta->x(), *actual_delta->s());
      PerformDualStep(alpha, alpha, actual_delta);

      // Check if that point is acceptable with respect to the current
      // original filter
      try {
        IpCq().trial_barrier_obj();
        IpCq().trial_constraint_violation();
        done=true;
      }
      catch(IpoptNLP::Eval_Error& e) {
        e.ReportException(Jnlst(), J_DETAILED);
        Jnlst().Printf(J_WARNING, J_LINE_SEARCH,
                       "Warning: Evaluation error during soft restoration phase step.\n");
        IpData().Append_info_string("e");
        count--;
      }
    }
    if (!done) {
      return false;
    }

    if (acceptor_->CheckAcceptabilityOfTrialPoint(0.)) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "  Trial step acceptable with respect to original backtracking globalization.\n");
      satisfies_original_criterion = true;
      return true;
    }

    // Evaluate the optimality error at the new point
    Number mu = .0;
    if (!IpData().FreeMuMode()) {
      mu = IpData().curr_mu();
    }
    Number trial_pderror;
    Number curr_pderror;
    try {
      trial_pderror = IpCq().trial_primal_dual_system_error(mu);
      curr_pderror = IpCq().curr_primal_dual_system_error(mu);
    }
    catch(IpoptNLP::Eval_Error& e) {
      e.ReportException(Jnlst(), J_DETAILED);
      Jnlst().Printf(J_WARNING, J_LINE_SEARCH,
                     "Warning: Evaluation error during soft restoration phase step.\n");
      IpData().Append_info_string("e");
      return false;
    }

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  Primal-dual error at current point:  %23.16e\n", curr_pderror);
    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  Primal-dual error at trial point  :  %23.16e\n", trial_pderror);
    // Check if there is sufficient reduction in the optimality error
    if (trial_pderror <= soft_resto_pderror_reduction_factor_*curr_pderror) {
      Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                     "  Trial step accepted.\n");
      return true;
    }

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "  Trial step rejected.\n");
    return false;
  }

  bool
  BacktrackingLineSearch::DetectTinyStep()
  {
    DBG_START_METH("BacktrackingLineSearch::DetectTinyStep",
                   dbg_verbosity);

    Number max_step_x;
    Number max_step_s;

    if (tiny_step_tol_==0.)
      return false;

    // ToDo try to find more efficient implementation
    DBG_PRINT_VECTOR(2, "curr_x", *IpData().curr()->x());
    DBG_PRINT_VECTOR(2, "delta_x", *IpData().delta()->x());

    SmartPtr<Vector> tmp = IpData().curr()->x()->MakeNewCopy();
    tmp->ElementWiseAbs();
    tmp->AddScalar(1.);

    SmartPtr<Vector> tmp2 = IpData().delta()->x()->MakeNewCopy();
    tmp2->ElementWiseDivide(*tmp);
    max_step_x = tmp2->Amax();
    Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                   "Relative step size for delta_x = %e\n",
                   max_step_x);
    if (max_step_x > tiny_step_tol_)
      return false;

    tmp = IpData().curr()->s()->MakeNew();
    tmp->Copy(*IpData().curr()->s());
    tmp->ElementWiseAbs();
    tmp->AddScalar(1.);

    tmp2 = IpData().curr()->s()->MakeNew();
    tmp2->Copy(*IpData().delta()->s());
    tmp2->ElementWiseDivide(*tmp);
    max_step_s = tmp2->Amax();
    Jnlst().Printf(J_MOREDETAILED, J_LINE_SEARCH,
                   "Relative step size for delta_s = %e\n",
                   max_step_s);
    if (max_step_s > tiny_step_tol_)
      return false;

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Tiny step of relative size %e detected.\n",
                   Max(max_step_x, max_step_s));

    return true;
  }

  bool BacktrackingLineSearch::CurrentIsAcceptable()
  {
    return (IsValid(conv_check_) &&
            conv_check_->CurrentIsAcceptable());
  }

  void BacktrackingLineSearch::StoreAcceptablePoint()
  {
    DBG_START_METH("BacktrackingLineSearch::StoreAcceptablePoint",
                   dbg_verbosity);

    acceptable_iterate_ = IpData().curr();
    acceptable_iteration_number_ = IpData().iter_count();
  }

  bool BacktrackingLineSearch::RestoreAcceptablePoint()
  {
    DBG_START_METH("BacktrackingLineSearch::RestoreAcceptablePoint",
                   dbg_verbosity);

    if (!IsValid(acceptable_iterate_)) {
      return false;
    }

    SmartPtr<IteratesVector> prev_iterate = acceptable_iterate_->MakeNewContainer();
    IpData().set_trial(prev_iterate);
    IpData().AcceptTrialPoint();

    return true;
  }

  bool BacktrackingLineSearch::ActivateFallbackMechanism()
  {
    // If we don't have a restoration phase, we don't know what to do
    if (IsNull(resto_phase_)) {
      return false;
    }

    // Reverting to the restoration phase only makes sense if there
    // are constraints
    if (IpData().curr()->y_c()->Dim()+IpData().curr()->y_d()->Dim()==0) {
      return false;
    }

    fallback_activated_ = true;
    rigorous_ = true;

    Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                   "Fallback option activated in BacktrackingLineSearch!\n");

    return true;
  }

} // namespace Ipopt
