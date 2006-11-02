// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpIpoptAlg.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpIpoptAlg.hpp"
#include "IpJournalist.hpp"
#include "IpRestoPhase.hpp"
#include "IpOrigIpoptNLP.hpp"

namespace Ipopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  static bool message_printed = false;

  static void print_message(const Journalist& jnlst)
  {
    jnlst.Printf(J_INSUPPRESSIBLE, J_MAIN,
                 "\n******************************************************************************\n"
                 "This program contains Ipopt, a library for large-scale nonlinear optimization.\n"
                 " Ipopt is released as open source code under the Common Public License (CPL).\n"
                 "         For more information visit http://projects.coin-or.org/Ipopt\n"
                 "******************************************************************************\n\n");
    message_printed = true;
  }

  IpoptAlgorithm::IpoptAlgorithm(const SmartPtr<PDSystemSolver>& pd_solver,
                                 const SmartPtr<LineSearch>& line_search,
                                 const SmartPtr<MuUpdate>& mu_update,
                                 const SmartPtr<ConvergenceCheck>& conv_check,
                                 const SmartPtr<IterateInitializer>& iterate_initializer,
                                 const SmartPtr<IterationOutput>& iter_output,
                                 const SmartPtr<HessianUpdater>& hessian_updater,
                                 const SmartPtr<EqMultiplierCalculator>& eq_multiplier_calculator /* = NULL*/)
      :
      pd_solver_(pd_solver),
      line_search_(line_search),
      mu_update_(mu_update),
      conv_check_(conv_check),
      iterate_initializer_(iterate_initializer),
      iter_output_(iter_output),
      hessian_updater_(hessian_updater),
      eq_multiplier_calculator_(eq_multiplier_calculator)
  {
    DBG_START_METH("IpoptAlgorithm::IpoptAlgorithm",
                   dbg_verbosity);
    DBG_ASSERT(IsValid(pd_solver_));
    DBG_ASSERT(IsValid(line_search_));
    DBG_ASSERT(IsValid(mu_update_));
    DBG_ASSERT(IsValid(conv_check_));
    DBG_ASSERT(IsValid(iterate_initializer_));
    DBG_ASSERT(IsValid(iter_output_));
    DBG_ASSERT(IsValid(hessian_updater_));
  }

  IpoptAlgorithm::~IpoptAlgorithm()
  {
    DBG_START_METH("IpoptAlgorithm::~IpoptAlgorithm()",
                   dbg_verbosity);
  }

  void IpoptAlgorithm::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Line Search");
    roptions->AddLowerBoundedNumberOption(
      "kappa_sigma",
      "Factor limiting the deviation of dual variables from primal estimates.",
      0, true, 1e10,
      "If the dual variables deviate from their primal estimates, a correction "
      "is performed. (See Eqn. (16) in the implementation paper.) "
      "Setting the value to less than 1 disables the correction.");
    roptions->AddStringOption2(
      "recalc_y",
      "Tells the algorithm to recalculate the equality and inequality multipliers as least square estimates.",
      "no",
      "no", "use the Newton step to update the multipliers",
      "yes", "use least-square multiplier estimates",
      "This asks the algorithm to recompute the multipliers, whenever the "
      "current infeasibility is less than recalc_y_feas_tol. "
      "Choosing yes might be helpful in the quasi-Newton option.  However, "
      "each recalculation requires an extra factorization of the linear "
      "system.  If a limited memory quasi-Newton option is chosen, this is "
      "used by default.");
    roptions->AddLowerBoundedNumberOption(
      "recalc_y_feas_tol",
      "Feasibility threshold for recomputation of multipliers.",
      0, true, 1e-6,
      "If recalc_y is chosen and the current infeasibility is less than this "
      "value, then the multipliers are recomputed.");
  }

  bool IpoptAlgorithm::InitializeImpl(const OptionsList& options,
                                      const std::string& prefix)
  {
    DBG_START_METH("IpoptAlgorithm::InitializeImpl",
                   dbg_verbosity);

    // Read the IpoptAlgorithm options
    // Initialize the Data object
    bool retvalue = IpData().Initialize(Jnlst(),
                                        options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the IpIpoptData object failed to initialize.");

    // Initialize the CQ object
    retvalue = IpCq().Initialize(Jnlst(),
                                 options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the IpIpoptCalculatedQuantities object failed to initialize.");

    // Initialize the CQ object
    retvalue = IpNLP().Initialize(Jnlst(),
                                  options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the IpIpoptNLP object failed to initialize.");

    // Initialize all the strategies
    retvalue = iterate_initializer_->Initialize(Jnlst(), IpNLP(), IpData(),
               IpCq(), options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the iterate_initializer strategy failed to initialize.");

    retvalue = mu_update_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                      options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the mu_update strategy failed to initialize.");

    retvalue = pd_solver_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                      options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the pd_solver strategy failed to initialize.");

    retvalue = line_search_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                        options,prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the line_search strategy failed to initialize.");

    retvalue = conv_check_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                       options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the conv_check strategy failed to initialize.");

    retvalue = iter_output_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                        options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the iter_output strategy failed to initialize.");

    retvalue = hessian_updater_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                            options, prefix);
    ASSERT_EXCEPTION(retvalue, FAILED_INITIALIZATION,
                     "the hessian_updater strategy failed to initialize.");

    options.GetNumericValue("kappa_sigma", kappa_sigma_, prefix);
    if (!options.GetBoolValue("recalc_y", recalc_y_, prefix)) {
      Index enum_int;
      if (options.GetEnumValue("hessian_approximation", enum_int, prefix)) {
        HessianApproximationType hessian_approximation =
          HessianApproximationType(enum_int);
        if (hessian_approximation==LIMITED_MEMORY) {
          recalc_y_ = true;
        }
      }
    }
    if (recalc_y_) {
      options.GetNumericValue("recalc_y_feas_tol", recalc_y_feas_tol_, prefix);
    }

    if (prefix=="resto.") {
      skip_print_problem_stats_ = true;
    }
    else {
      skip_print_problem_stats_ = false;
    }

    return true;
  }

  SolverReturn IpoptAlgorithm::Optimize()
  {
    DBG_START_METH("IpoptAlgorithm::Optimize", dbg_verbosity);

    // Start measuring CPU time
    IpData().TimingStats().OverallAlgorithm().Start();

    if (!message_printed) {
      print_message(Jnlst());
    }

    try {
      IpData().TimingStats().InitializeIterates().Start();
      // Initialize the iterates
      InitializeIterates();
      IpData().TimingStats().InitializeIterates().End();

      if (!skip_print_problem_stats_) {
        IpData().TimingStats().PrintProblemStatistics().Start();
        PrintProblemStatistics();
        IpData().TimingStats().PrintProblemStatistics().End();
      }

      IpData().TimingStats().CheckConvergence().Start();
      ConvergenceCheck::ConvergenceStatus conv_status
      = conv_check_->CheckConvergence();
      IpData().TimingStats().CheckConvergence().End();

      // main loop
      while (conv_status == ConvergenceCheck::CONTINUE) {
        // Set the Hessian Matrix
        IpData().TimingStats().UpdateHessian().Start();
        UpdateHessian();
        IpData().TimingStats().UpdateHessian().End();

        // do all the output for this iteration
        IpData().TimingStats().OutputIteration().Start();
        OutputIteration();
        IpData().ResetInfo();
        IpData().TimingStats().OutputIteration().End();

        // initialize the flag that is set to true if the algorithm
        // has to continue with an emergency fallback mode.  For
        // example, when no search direction can be computed, continue
        // with the restoration phase
        bool emergency_mode = false;

        // update the barrier parameter
        IpData().TimingStats().UpdateBarrierParameter().Start();
        emergency_mode = !UpdateBarrierParameter();
        IpData().TimingStats().UpdateBarrierParameter().End();

        if (!emergency_mode) {
          // solve the primal-dual system to get the full step
          IpData().TimingStats().ComputeSearchDirection().Start();
          emergency_mode = !ComputeSearchDirection();
          IpData().TimingStats().ComputeSearchDirection().End();
        }

        // If we are in the emergency mode, ask to line search object
        // to go to the fallback options.  If that isn't possible,
        // issue error message
        if (emergency_mode) {
          bool retval = line_search_->ActivateFallbackMechanism();
          if (retval) {
            Jnlst().Printf(J_WARNING, J_MAIN,
                           "WARNING: Problem in step computation; switching to emergency mode.\n");
          }
          else {
            Jnlst().Printf(J_ERROR, J_MAIN,
                           "ERROR: Problem in step computation, but emergency mode cannot be activated.\n");
            THROW_EXCEPTION(STEP_COMPUTATION_FAILED,
                            "Step computation failed.");
          }
        }

        // Compute the new iterate
        IpData().TimingStats().ComputeAcceptableTrialPoint().Start();
        ComputeAcceptableTrialPoint();
        IpData().TimingStats().ComputeAcceptableTrialPoint().End();

        // Accept the new iterate
        IpData().TimingStats().AcceptTrialPoint().Start();
        AcceptTrialPoint();
        IpData().TimingStats().AcceptTrialPoint().End();

        IpData().Set_iter_count(IpData().iter_count()+1);

        IpData().TimingStats().CheckConvergence().Start();
        conv_status  = conv_check_->CheckConvergence();
        IpData().TimingStats().CheckConvergence().End();
      }

      IpData().TimingStats().OutputIteration().Start();
      OutputIteration();
      IpData().TimingStats().OutputIteration().End();

      IpData().TimingStats().OverallAlgorithm().End();

      if (conv_status == ConvergenceCheck::CONVERGED) {
        return SUCCESS;
      }
      else if (conv_status == ConvergenceCheck::CONVERGED_TO_ACCEPTABLE_POINT) {
        return STOP_AT_ACCEPTABLE_POINT;
      }
      else if (conv_status == ConvergenceCheck::MAXITER_EXCEEDED) {
        return MAXITER_EXCEEDED;
      }
      else if (conv_status == ConvergenceCheck::DIVERGING) {
        return DIVERGING_ITERATES;
      }
      else if (conv_status == ConvergenceCheck::USER_STOP) {
        return USER_REQUESTED_STOP;
      }
    }
    catch(TINY_STEP_DETECTED& exc) {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().UpdateBarrierParameter().EndIfStarted();
      IpData().TimingStats().OverallAlgorithm().End();
      return STOP_AT_TINY_STEP;
    }
    catch(ACCEPTABLE_POINT_REACHED& exc) {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      IpData().TimingStats().OverallAlgorithm().End();
      return STOP_AT_ACCEPTABLE_POINT;
    }
    catch(LOCALLY_INFEASIBLE& exc) {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      IpData().TimingStats().CheckConvergence().EndIfStarted();
      IpData().TimingStats().OverallAlgorithm().End();
      return LOCAL_INFEASIBILITY;
    }
    catch(RESTORATION_CONVERGED_TO_FEASIBLE_POINT& exc) {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      IpData().TimingStats().OverallAlgorithm().End();
      return RESTORATION_FAILURE;
    }
    catch(RESTORATION_FAILED& exc) {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      IpData().TimingStats().OverallAlgorithm().End();
      return RESTORATION_FAILURE;
    }
    catch(RESTORATION_MAXITER_EXCEEDED& exc) {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      IpData().TimingStats().OverallAlgorithm().End();
      return MAXITER_EXCEEDED;
    }
    catch(RESTORATION_USER_STOP& exc) {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      IpData().TimingStats().OverallAlgorithm().End();
      return USER_REQUESTED_STOP;
    }
    catch(STEP_COMPUTATION_FAILED& exc) {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      IpData().TimingStats().OverallAlgorithm().End();
      return ERROR_IN_STEP_COMPUTATION;
    }
    catch(IpoptNLP::Eval_Error& exc) {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      IpData().TimingStats().OverallAlgorithm().End();
      return INVALID_NUMBER_DETECTED;
    }
    catch(FEASIBILITY_PROBLEM_SOLVED& exc) {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      IpData().TimingStats().OverallAlgorithm().End();
      return SUCCESS;
    }
    catch(TOO_FEW_DOF& exc) {
      exc.ReportException(Jnlst(), J_MOREDETAILED);
      IpData().TimingStats().ComputeAcceptableTrialPoint().EndIfStarted();
      IpData().TimingStats().OverallAlgorithm().End();
      return TOO_FEW_DEGREES_OF_FREEDOM;
    }
    catch(INTERNAL_ABORT& exc) {
      exc.ReportException(Jnlst());
      IpData().TimingStats().OverallAlgorithm().End();
      return INTERNAL_ERROR;
    }

    DBG_ASSERT(false && "Unknown return code in the algorithm");

    IpData().TimingStats().OverallAlgorithm().End();
    return INTERNAL_ERROR;
  }

  void IpoptAlgorithm::UpdateHessian()
  {
    Jnlst().Printf(J_DETAILED, J_MAIN, "\n**************************************************\n");
    Jnlst().Printf(J_DETAILED, J_MAIN, "*** Update HessianMatrix for Iteration %d:", IpData().iter_count());
    Jnlst().Printf(J_DETAILED, J_MAIN, "\n**************************************************\n\n");
    hessian_updater_->UpdateHessian();
  }

  bool IpoptAlgorithm::UpdateBarrierParameter()
  {
    Jnlst().Printf(J_DETAILED, J_MAIN, "\n**************************************************\n");
    Jnlst().Printf(J_DETAILED, J_MAIN, "*** Update Barrier Parameter for Iteration %d:", IpData().iter_count());
    Jnlst().Printf(J_DETAILED, J_MAIN, "\n**************************************************\n\n");
    bool retval = mu_update_->UpdateBarrierParameter();

    if (retval) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Barrier Parameter: %e\n", IpData().curr_mu());
    }
    else {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Barrier parameter could not be updated!\n");
    }

    return retval;
  }

  bool IpoptAlgorithm::ComputeSearchDirection()
  {
    DBG_START_METH("IpoptAlgorithm::ComputeSearchDirection", dbg_verbosity);

    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n**************************************************\n");
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "*** Solving the Primal Dual System for Iteration %d:",
                   IpData().iter_count());
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n**************************************************\n\n");

    bool improve_solution = false;
    if (IpData().HaveDeltas()) {
      improve_solution = true;
    }

    SmartPtr<IteratesVector> rhs = IpData().curr()->MakeNewContainer();
    rhs->Set_x(*IpCq().curr_grad_lag_with_damping_x());
    rhs->Set_s(*IpCq().curr_grad_lag_with_damping_s());
    rhs->Set_y_c(*IpCq().curr_c());
    rhs->Set_y_d(*IpCq().curr_d_minus_s());
    rhs->Set_z_L(*IpCq().curr_relaxed_compl_x_L());
    rhs->Set_z_U(*IpCq().curr_relaxed_compl_x_U());
    rhs->Set_v_L(*IpCq().curr_relaxed_compl_s_L());
    rhs->Set_v_U(*IpCq().curr_relaxed_compl_s_U());

    DBG_PRINT_VECTOR(2, "rhs", *rhs);

    // Get space for the search direction
    SmartPtr<IteratesVector> delta =
      IpData().curr()->MakeNewIteratesVector(true);

    if (improve_solution) {
      // We can probably avoid copying and scaling...
      delta->AddOneVector(-1., *IpData().delta(), 0.);
    }

    bool allow_inexact = false;
    bool retval = pd_solver_->Solve(-1.0, 0.0, *rhs, *delta, allow_inexact,
                                    improve_solution);

    if (retval) {
      // Store the search directions in the IpData object
      IpData().set_delta(delta);

      Jnlst().Printf(J_MOREVECTOR, J_MAIN,
                     "*** Step Calculated for Iteration: %d\n",
                     IpData().iter_count());
      IpData().delta()->Print(Jnlst(), J_MOREVECTOR, J_MAIN, "delta");
    }
    else {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "*** Step could not be computed in iteration %d!\n",
                     IpData().iter_count());
    }

    return retval;
  }

  void IpoptAlgorithm::ComputeAcceptableTrialPoint()
  {
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n**************************************************\n");
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "*** Finding Acceptable Trial Point for Iteration %d:",
                   IpData().iter_count());
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n**************************************************\n\n");
    line_search_->FindAcceptableTrialPoint();
  }

  void IpoptAlgorithm::OutputIteration()
  {
    iter_output_->WriteOutput();
  }

  void IpoptAlgorithm::InitializeIterates()
  {
    DBG_START_METH("IpoptAlgorithm::InitializeIterates", dbg_verbosity);

    iterate_initializer_->SetInitialIterates();
  }

  void IpoptAlgorithm::AcceptTrialPoint()
  {
    DBG_START_METH("IpoptAlgorithm::AcceptTrialPoint", dbg_verbosity);
    // If the line search didn't determine a new acceptable trial
    // point, do not accept a new iterate
    if (line_search_->CheckSkippedLineSearch()) {
      Jnlst().Printf(J_SUMMARY, J_MAIN,
                     "Line search didn't find acceptable trial point.\n");
      return;
    }

    // Adjust the bounds if necessary
    Index adjusted_slacks = IpCq().AdjustedTrialSlacks();
    DBG_PRINT((1, "adjusted_slacks = %d\n", adjusted_slacks));
    if (adjusted_slacks>0) {
      IpCq().ResetAdjustedTrialSlacks();
      if (adjusted_slacks==1) {
        Jnlst().Printf(J_WARNING, J_MAIN,
                       "In iteration %d, %d Slack too small, adjusting variable bound\n",
                       IpData().iter_count(), adjusted_slacks);
      }
      else {
        Jnlst().Printf(J_WARNING, J_MAIN,
                       "In iteration %d, %d Slacks too small, adjusting variable bounds\n",
                       IpData().iter_count(), adjusted_slacks);
      }
      if (Jnlst().ProduceOutput(J_VECTOR, J_MAIN)) {
        IpNLP().x_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "old_x_L");
        IpNLP().x_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "old_x_U");
        IpNLP().d_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "old_d_L");
        IpNLP().d_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "old_d_U");
      }

      SmartPtr<Vector> new_x_l = IpNLP().x_L()->MakeNew();
      IpNLP().Px_L()->TransMultVector(1.0, *IpData().trial()->x(),
                                      0.0, *new_x_l);
      new_x_l->Axpy(-1.0, *IpCq().trial_slack_x_L());

      SmartPtr<Vector> new_x_u = IpNLP().x_U()->MakeNew();
      IpNLP().Px_U()->TransMultVector(1.0, *IpData().trial()->x(),
                                      0.0, *new_x_u);
      new_x_u->Axpy(1.0, *IpCq().trial_slack_x_U());

      SmartPtr<Vector> new_d_l = IpNLP().d_L()->MakeNew();
      IpNLP().Pd_L()->TransMultVector(1.0, *IpData().trial()->s(),
                                      0.0, *new_d_l);
      new_d_l->Axpy(-1.0, *IpCq().trial_slack_s_L());

      SmartPtr<Vector> new_d_u = IpNLP().d_U()->MakeNew();
      IpNLP().Pd_U()->TransMultVector(1.0, *IpData().trial()->s(),
                                      0.0, *new_d_u);
      new_d_u->Axpy(1.0, *IpCq().trial_slack_s_U());

      IpNLP().AdjustVariableBounds(*new_x_l, *new_x_u, *new_d_l, *new_d_u);

      if (Jnlst().ProduceOutput(J_VECTOR, J_MAIN)) {
        IpNLP().x_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "new_x_L");
        IpNLP().x_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "new_x_U");
        IpNLP().d_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "new_d_L");
        IpNLP().d_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "new_d_U");
      }

    }

    // Make sure that bound multipliers are not too far from \mu * S^{-1}
    // (see kappa_sigma in paper)
    bool corrected = false;
    Number max_correction;
    SmartPtr<const Vector> new_z_L;
    max_correction = correct_bound_multiplier(
                       *IpData().trial()->z_L(),
                       *IpCq().trial_slack_x_L(),
                       *IpCq().trial_compl_x_L(),
                       new_z_L);
    if (max_correction>0.) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Some value in z_L becomes too large - maximal correction = %8.2e\n",
                     max_correction);
      corrected = true;
    }
    SmartPtr<const Vector> new_z_U;
    max_correction = correct_bound_multiplier(
                       *IpData().trial()->z_U(),
                       *IpCq().trial_slack_x_U(),
                       *IpCq().trial_compl_x_U(),
                       new_z_U);
    if (max_correction>0.) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Some value in z_U becomes too large - maximal correction = %8.2e\n",
                     max_correction);
      corrected = true;
    }
    SmartPtr<const Vector> new_v_L;
    max_correction = correct_bound_multiplier(
                       *IpData().trial()->v_L(),
                       *IpCq().trial_slack_s_L(),
                       *IpCq().trial_compl_s_L(),
                       new_v_L);
    if (max_correction>0.) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Some value in v_L becomes too large - maximal correction = %8.2e\n",
                     max_correction);
      corrected = true;
    }
    SmartPtr<const Vector> new_v_U;
    max_correction = correct_bound_multiplier(
                       *IpData().trial()->v_U(),
                       *IpCq().trial_slack_s_U(),
                       *IpCq().trial_compl_s_U(),
                       new_v_U);
    if (max_correction>0.) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Some value in v_U becomes too large - maximal correction = %8.2e\n",
                     max_correction);
      corrected = true;
    }
    SmartPtr<IteratesVector> trial = IpData().trial()->MakeNewContainer();
    trial->Set_bound_mult(*new_z_L, *new_z_U, *new_v_L, *new_v_U);
    IpData().set_trial(trial);

    if (corrected) {
      IpData().Append_info_string("z");
    }

    // Accept the step
    IpData().AcceptTrialPoint();

    // If we want to recalculate the multipliers (e.g., as least
    // square estimates), call the calculator for that
    if (recalc_y_) {
      // There is no point in doing this if there are no constraints
      if (IpData().curr()->y_c()->Dim()+IpData().curr()->y_d()->Dim()==0) {
        recalc_y_ = false;
      }
    }
    if (recalc_y_ && IpCq().curr_constraint_violation()<recalc_y_feas_tol_) {
      if (Jnlst().ProduceOutput(J_MOREDETAILED, J_MAIN)) {
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "dual infeasisibility before least square multiplier update = %e\n",
                       IpCq().curr_dual_infeasibility(NORM_MAX));
      }
      IpData().Append_info_string("y ");
      DBG_ASSERT(IsValid(eq_multiplier_calculator_));
      if (IpData().curr()->y_c()->Dim()+IpData().curr()->y_d()->Dim()>0) {
        SmartPtr<Vector> y_c = IpData().curr()->y_c()->MakeNew();
        SmartPtr<Vector> y_d = IpData().curr()->y_d()->MakeNew();
        bool retval =
          eq_multiplier_calculator_->CalculateMultipliers(*y_c, *y_d);
        if (retval) {
          SmartPtr<const IteratesVector> curr = IpData().curr();
          SmartPtr<IteratesVector> iterates = curr->MakeNewContainer();
          iterates->Set_x(*curr->x());
          iterates->Set_s(*curr->s());
          iterates->Set_z_L(*curr->z_L());
          iterates->Set_z_U(*curr->z_U());
          iterates->Set_v_L(*curr->v_L());
          iterates->Set_v_U(*curr->v_U());
          iterates->Set_y_c(*y_c);
          iterates->Set_y_d(*y_d);
          IpData().set_trial(iterates);
          IpData().AcceptTrialPoint();
        }
        else {
          Jnlst().Printf(J_DETAILED, J_MAIN,
                         "Recalculation of y multipliers skipped because eq_mult_calc returned false.\n");
        }
      }
    }
  }

  void IpoptAlgorithm::PrintProblemStatistics()
  {
    if (!Jnlst().ProduceOutput(J_SUMMARY, J_STATISTICS)) {
      // nothing to print
      return;
    }

    SmartPtr<const Vector> x = IpData().curr()->x();
    SmartPtr<const Vector> x_L = IpNLP().x_L();
    SmartPtr<const Vector> x_U = IpNLP().x_U();
    SmartPtr<const Matrix> Px_L = IpNLP().Px_L();
    SmartPtr<const Matrix> Px_U = IpNLP().Px_U();

    Index nx_tot, nx_only_lower, nx_both, nx_only_upper;
    calc_number_of_bounds(*IpData().curr()->x(), *IpNLP().x_L(), *IpNLP().x_U(),
                          *IpNLP().Px_L(), *IpNLP().Px_U(),
                          nx_tot, nx_only_lower, nx_both, nx_only_upper);

    Index ns_tot, ns_only_lower, ns_both, ns_only_upper;
    calc_number_of_bounds(*IpData().curr()->s(), *IpNLP().d_L(), *IpNLP().d_U(),
                          *IpNLP().Pd_L(), *IpNLP().Pd_U(),
                          ns_tot, ns_only_lower, ns_both, ns_only_upper);

    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "Total number of variables............................: %8d\n",nx_tot);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "                     variables with only lower bounds: %8d\n",
                   nx_only_lower);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "                variables with lower and upper bounds: %8d\n",nx_both);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "                     variables with only upper bounds: %8d\n",
                   nx_only_upper);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "Total number of equality constraints.................: %8d\n",
                   IpData().curr()->y_c()->Dim());
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "Total number of inequality constraints...............: %8d\n",ns_tot);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "        inequality constraints with only lower bounds: %8d\n",
                   ns_only_lower);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "   inequality constraints with lower and upper bounds: %8d\n",ns_both);
    Jnlst().Printf(J_SUMMARY, J_STATISTICS,
                   "        inequality constraints with only upper bounds: %8d\n\n",
                   ns_only_upper);
  }

  void IpoptAlgorithm::calc_number_of_bounds(
    const Vector& x,
    const Vector& x_L,
    const Vector& x_U,
    const Matrix& Px_L,
    const Matrix& Px_U,
    Index& n_tot,
    Index& n_only_lower,
    Index& n_both,
    Index& n_only_upper)
  {
    DBG_START_METH("IpoptAlgorithm::calc_number_of_bounds",
                   dbg_verbosity);

    n_tot = x.Dim();

    SmartPtr<Vector> tmpx = x.MakeNew();
    SmartPtr<Vector> tmpxL = x_L.MakeNew();
    SmartPtr<Vector> tmpxU = x_U.MakeNew();

    tmpxL->Set(-1.);
    tmpxU->Set(2.);
    Px_L.MultVector(1.0, *tmpxL, 0.0, *tmpx);
    Px_U.MultVector(1.0, *tmpxU, 1.0, *tmpx);
    // Now, x has elements
    //  -1 : if component has only lower bound
    //   0 : if component has no bound
    //   1 : if component has both lower and upper bound
    //   2 : if component has only upper bound
    DBG_PRINT_VECTOR(2, "x-indicator", *tmpx);

    SmartPtr<Vector> tmpx0 = x.MakeNew();
    tmpx0->Set(0.);

    SmartPtr<Vector> tmpx2 = x.MakeNew();
    tmpx2->Set(-1.0);
    tmpx2->Axpy(1.0, *tmpx);
    tmpx2->ElementWiseMax(*tmpx0); // tmpx2 is now 1 in those
    // components with only upper bounds
    n_only_upper = (Index)tmpx2->Asum();

    tmpx->Axpy(-2., *tmpx2);       // now make all those entries for
    // only upper bounds zero in tmpx

    tmpx2->Copy(*tmpx);
    tmpx2->ElementWiseMax(*tmpx0); // tmpx2 is now 1 in those
    // components with both bounds
    n_both = (Index)tmpx2->Asum();

    tmpx->Axpy(-1., *tmpx2);
    tmpx->ElementWiseMin(*tmpx);   // tmpx is now -1 in those with only
    // lower bounds
    n_only_lower = (Index)tmpx->Asum();

  }

  Number IpoptAlgorithm::correct_bound_multiplier(
    const Vector& trial_z,
    const Vector& trial_slack,
    const Vector& trial_compl,
    SmartPtr<const Vector>& new_trial_z)
  {
    DBG_START_METH("IpoptAlgorithm::CorrectBoundMultiplier",
                   dbg_verbosity);

    if (kappa_sigma_<1. || trial_z.Dim()==0) {
      new_trial_z = &trial_z;
      return 0.;
    }

    // We choose as barrier parameter to be used either the current
    // algorithmic barrier parameter (if we are not in the free mode),
    // or the average complementarity (at the trial point)
    Number mu;
    if (IpData().FreeMuMode()) {
      mu = IpCq().trial_avrg_compl();
      mu = Min(mu, 1e3);
    }
    else {
      mu = IpData().curr_mu();
    }
    DBG_PRINT((1,"mu = %8.2e\n", mu));
    DBG_PRINT_VECTOR(2, "trial_z", trial_z);

    // First check quickly if anything need to be corrected, using the
    // trial complementarity directly.  Here, Amax is the same as Max
    // (and we use Amax because that can be used later)
    if (trial_compl.Amax() <= kappa_sigma_*mu &&
        trial_compl.Min() >= 1./kappa_sigma_*mu) {
      new_trial_z = &trial_z;
      return 0.;
    }

    SmartPtr<Vector> one_over_s = trial_z.MakeNew();
    one_over_s->Copy(trial_slack);
    one_over_s->ElementWiseReciprocal();

    SmartPtr<Vector> step_z = trial_z.MakeNew();
    step_z->AddTwoVectors(kappa_sigma_*mu, *one_over_s, -1., trial_z, 0.);

    DBG_PRINT_VECTOR(2, "step_z", *step_z);

    Number max_correction_up = Max(0., -step_z->Min());
    if (max_correction_up>0.) {
      SmartPtr<Vector> tmp = trial_z.MakeNew();
      tmp->Set(0.);
      step_z->ElementWiseMin(*tmp);
      tmp->AddTwoVectors(1., trial_z, 1., *step_z, 0.);
      new_trial_z = GetRawPtr(tmp);
    }
    else {
      new_trial_z = &trial_z;
    }

    step_z->AddTwoVectors(1./kappa_sigma_*mu, *one_over_s, -1., *new_trial_z, 0.);

    Number max_correction_low = Max(0., step_z->Max());
    if (max_correction_low>0.) {
      SmartPtr<Vector> tmp = trial_z.MakeNew();
      tmp->Set(0.);
      step_z->ElementWiseMax(*tmp);
      tmp->AddTwoVectors(1., *new_trial_z, 1., *step_z, 0.);
      new_trial_z = GetRawPtr(tmp);
    }

    DBG_PRINT_VECTOR(2, "new_trial_z", *new_trial_z);

    return Max(max_correction_up, max_correction_low);
  }

} // namespace Ipopt
