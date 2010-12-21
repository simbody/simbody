// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpOptErrorConvCheck.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpOptErrorConvCheck.hpp"

namespace Ipopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  OptimalityErrorConvergenceCheck::OptimalityErrorConvergenceCheck()
  {}

  OptimalityErrorConvergenceCheck::~OptimalityErrorConvergenceCheck()
  {}

  void OptimalityErrorConvergenceCheck::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedIntegerOption(
      "max_iter",
      "Maximum number of iterations.",
      0, 3000,
      "The algorithm terminates with an error message if the number of "
      "iterations exceeded this number.");
    roptions->AddLowerBoundedNumberOption(
      "dual_inf_tol",
      "Desired threshold for the dual infeasibility.",
      0.0, true, 1e-4,
      "Absolute tolerance on the dual infeasibility. Successful termination "
      "requires that the max-norm of the (unscaled) dual infeasibility is less than this "
      "threshold.");
    roptions->AddLowerBoundedNumberOption(
      "constr_viol_tol",
      "Desired threshold for the constraint violation.",
      0.0, true, 1e-4,
      "Absolute tolerance on the constraint violation. Successful termination "
      "requires that the max-norm of the (unscaled) constraint violation is less than this "
      "threshold.");
    roptions->AddLowerBoundedNumberOption(
      "compl_inf_tol",
      "Desired threshold for the complementarity conditions.",
      0.0, true, 1e-4,
      "Absolute tolerance on the complementarity. Successful termination "
      "requires that the max-norm of the (unscaled) complementarity is less than this "
      "threshold.");
    roptions->AddLowerBoundedNumberOption(
      "acceptable_tol",
      "\"Acceptable\" convergence tolerance (relative).",
      0.0, true,  1e-6,
      "Determines which (scaled) overall optimality error is considered to be"
      " \"acceptable.\" There are two levels of termination criteria.  If the "
      "usual \"desired\" tolerances (see tol, dual_inf_tol etc) are satisfied "
      "at an iteration, the algorithm immediately terminates with a success "
      "message.  On the other hand, if the algorithm encounters "
      "\"acceptable_iter\" many iterations in a row that are considered "
      "\"acceptable\", it will terminate before the desired convergence "
      "tolerance is met. This is useful in cases where the algorithm might "
      "not be able to achieve the \"desired\" level of accuracy.");
    roptions->AddLowerBoundedIntegerOption(
      "acceptable_iter",
      "Number of \"acceptable\" iterates before triggering termination.",
      0, 15,
      "If the algorithm encounters this many successive \"acceptable\" iterates "
      "(see \"acceptable_tol\"), it terminates, assuming that the problem "
      "has been solved to best possible accuracy given round-off.  If it is "
      "set to zero, this heuristic is disabled.");
    roptions->AddLowerBoundedNumberOption(
      "acceptable_dual_inf_tol",
      "\"Acceptance\" threshold for the dual infeasibility.",
      0.0, true, 1e-2,
      "Absolute tolerance on the dual infeasibility. \"Acceptable\" termination "
      "requires that the (max-norm of the unscaled) dual infeasibility is less than this "
      "threshold; see also acceptable_tol.");
    roptions->AddLowerBoundedNumberOption(
      "acceptable_constr_viol_tol",
      "\"Acceptance\" threshold for the constraint violation.",
      0.0, true, 1e-2,
      "Absolute tolerance on the constraint violation. \"Acceptable\" termination "
      "requires that the max-norm of the (unscaled) constraint violation is less than this "
      "threshold; see also acceptable_tol.");
    roptions->AddLowerBoundedNumberOption(
      "acceptable_compl_inf_tol",
      "\"Acceptance\" threshold for the complementarity conditions.",
      0.0, true, 1e-2,
      "Absolute tolerance on the complementarity. \"Acceptable\" termination "
      "requires that the max-norm of the (unscaled) complementarity is less than this "
      "threshold; see also acceptable_tol.");
    roptions->AddLowerBoundedNumberOption(
      "diverging_iterates_tol",
      "Threshold for maximal value of primal iterates.",
      0.0, true, 1e20,
      "If any component of the primal iterates exceeded this value (in "
      "absolute terms), the optimization is aborted with the exit message "
      "that the iterates seem to be diverging.");
  }

  bool
  OptimalityErrorConvergenceCheck::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetIntegerValue("max_iter", max_iterations_, prefix);
    options.GetNumericValue("dual_inf_tol", dual_inf_tol_, prefix);
    options.GetNumericValue("constr_viol_tol", constr_viol_tol_, prefix);
    options.GetNumericValue("compl_inf_tol", compl_inf_tol_, prefix);
    options.GetIntegerValue("acceptable_iter", acceptable_iter_, prefix);
    options.GetNumericValue("acceptable_tol", acceptable_tol_, prefix);
    options.GetNumericValue("acceptable_dual_inf_tol", acceptable_dual_inf_tol_, prefix);
    options.GetNumericValue("acceptable_constr_viol_tol", acceptable_constr_viol_tol_, prefix);
    options.GetNumericValue("acceptable_compl_inf_tol", acceptable_compl_inf_tol_, prefix);
    options.GetNumericValue("diverging_iterates_tol", diverging_iterates_tol_, prefix);
    acceptable_counter_ = 0;

    return true;
  }

  ConvergenceCheck::ConvergenceStatus
  OptimalityErrorConvergenceCheck::CheckConvergence(bool call_intermediate_callback /*= true*/)
  {
    DBG_START_METH("OptimalityErrorConvergenceCheck::CheckConvergence", dbg_verbosity);

    if (call_intermediate_callback) {
      // Check if user requested termination by calling the intermediate
      // user callback function
      AlgorithmMode mode = RegularMode;
      // Gather the information also used in the iteration output
      Index iter = IpData().iter_count();
      Number inf_pr = IpCq().curr_primal_infeasibility(NORM_MAX);
      Number inf_du = IpCq().curr_dual_infeasibility(NORM_MAX);
      Number mu = IpData().curr_mu();
      Number dnrm;
      if (IsValid(IpData().delta()) && IsValid(IpData().delta()->x()) && IsValid(IpData().delta()->s())) {
        dnrm = Max(IpData().delta()->x()->Amax(), IpData().delta()->s()->Amax());
      }
      else {
        // This is the first iteration - no search direction has been
        // computed yet.
        dnrm = 0.;
      }
      Number alpha_primal = IpData().info_alpha_primal();
      Number alpha_dual = IpData().info_alpha_dual();
      Number regu_x = IpData().info_regu_x();
      Number unscaled_f = IpCq().unscaled_curr_f();
      Index ls_count = IpData().info_ls_count();
      bool request_stop =
        !IpNLP().IntermediateCallBack(mode, iter, unscaled_f, inf_pr, inf_du,
                                      mu, dnrm, regu_x, alpha_dual,
                                      alpha_primal, ls_count,
                                      &IpData(), &IpCq());

      if (request_stop) {
        return ConvergenceCheck::USER_STOP;
      }
    }

    if (IpData().iter_count() >= max_iterations_) {
      return ConvergenceCheck::MAXITER_EXCEEDED;
    }

    Number overall_error = IpCq().curr_nlp_error();
    Number dual_inf = IpCq().unscaled_curr_dual_infeasibility(NORM_MAX);
    Number constr_viol = IpCq().unscaled_curr_nlp_constraint_violation(NORM_MAX);
    Number compl_inf = IpCq().unscaled_curr_complementarity(0., NORM_MAX);

    if (IpData().curr()->x()->Dim()==IpData().curr()->y_c()->Dim()) {
      // the problem is square, there is no point in looking at dual
      // infeasibility and complementarity as termination criterion
      dual_inf_tol_ = 1e300;
      compl_inf_tol_ = 1e300;
    }

    if (overall_error <= IpData().tol() &&
        dual_inf <= dual_inf_tol_ &&
        constr_viol <= constr_viol_tol_ &&
        compl_inf <= compl_inf_tol_) {
      return ConvergenceCheck::CONVERGED;
    }

    if (acceptable_iter_>0 && CurrentIsAcceptable()) {
      IpData().Append_info_string("A");
      acceptable_counter_++;
      if (acceptable_counter_ >= acceptable_iter_) {
        return ConvergenceCheck::CONVERGED_TO_ACCEPTABLE_POINT;
      }
    }
    else {
      acceptable_counter_ = 0;
    }

    if (IpData().curr()->x()->Amax() > diverging_iterates_tol_) {
      return ConvergenceCheck::DIVERGING;
    }

    return ConvergenceCheck::CONTINUE;
  }

  bool OptimalityErrorConvergenceCheck::CurrentIsAcceptable()
  {
    DBG_START_METH("OptimalityErrorConvergenceCheck::CurrentIsAcceptable",
                   dbg_verbosity);

    Number overall_error = IpCq().curr_nlp_error();
    Number dual_inf = IpCq().unscaled_curr_dual_infeasibility(NORM_MAX);
    Number constr_viol = IpCq().unscaled_curr_nlp_constraint_violation(NORM_MAX);
    Number compl_inf = IpCq().unscaled_curr_complementarity(0., NORM_MAX);

    DBG_PRINT((1, "overall_error = %e\n", overall_error));
    DBG_PRINT((1, "dual_inf = %e\n", dual_inf));
    DBG_PRINT((1, "constr_viol = %e\n", constr_viol));
    DBG_PRINT((1, "compl_inf = %e\n", compl_inf));

    DBG_PRINT((1, "acceptable_tol_ = %e\n", acceptable_tol_));
    DBG_PRINT((1, "acceptable_dual_inf_tol_ = %e\n", acceptable_dual_inf_tol_));
    DBG_PRINT((1, "acceptable_constr_viol_tol_ = %e\n", acceptable_constr_viol_tol_));
    DBG_PRINT((1, "acceptable_compl_inf_tol_ = %e\n", acceptable_compl_inf_tol_));

    if (IpData().curr()->x()->Dim()==IpData().curr()->y_c()->Dim()) {
      // the problem is square, there is no point in looking at dual
      // infeasibility and complementarity as termination criterion
      acceptable_dual_inf_tol_ = 1e300;
      acceptable_compl_inf_tol_ = 1e300;
    }

    return (overall_error <= acceptable_tol_ &&
            dual_inf <= acceptable_dual_inf_tol_ &&
            constr_viol <= acceptable_constr_viol_tol_ &&
            compl_inf <= acceptable_compl_inf_tol_);
  }


} // namespace Ipopt
