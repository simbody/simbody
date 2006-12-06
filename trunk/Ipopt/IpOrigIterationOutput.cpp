// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpOrigIterationOutput.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter              IBM    2004-09-23

#include "IpOrigIterationOutput.hpp"

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
  OrigIterationOutput::OrigIterationOutput()
  {}

  OrigIterationOutput::~OrigIterationOutput()
  {}

  void
  OrigIterationOutput::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    std::string prev_cat = roptions->RegisteringCategory();
    roptions->SetRegisteringCategory("Output");
    roptions->AddStringOption2(
      "print_info_string",
      "Enables printing of additional info string at end of iteration output.",
      "no",
      "no", "don't print string",
      "yes", "print string at end of each iteration output",
      "This string contains some insider information about the current iteration.");
    roptions->SetRegisteringCategory(prev_cat);
  }

  bool OrigIterationOutput::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetBoolValue("print_info_string", print_info_string_, prefix);

    return true;
  }

  void OrigIterationOutput::WriteOutput()
  {
    //////////////////////////////////////////////////////////////////////
    //         First print the summary line for the iteration           //
    //////////////////////////////////////////////////////////////////////

    Index iter = IpData().iter_count();
    std::string header =
      "iter    objective    inf_pr   inf_du lg(mu)  ||d||  lg(rg) alpha_du alpha_pr  ls\n";
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n\n**************************************************\n");
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "*** Summary of Iteration: %d:", IpData().iter_count());
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "\n**************************************************\n\n");
    if (iter%10 == 0 && !IpData().info_skip_output()) {
      // output the header
      Jnlst().Printf(J_ITERSUMMARY, J_MAIN, header.c_str());
    }
    else {
      Jnlst().Printf(J_DETAILED, J_MAIN, header.c_str());
    }
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
    Number unscaled_f = IpCq().unscaled_curr_f();

    // Retrieve some information set in the different parts of the algorithm
    char info_iter=' ';
    Number alpha_primal = IpData().info_alpha_primal();
    char alpha_primal_char = IpData().info_alpha_primal_char();
    Number alpha_dual = IpData().info_alpha_dual();
    Number regu_x = IpData().info_regu_x();
    char regu_x_buf[8];
    char dashes[]="   - ";
    char *regu_x_ptr;
    if (regu_x==.0) {
      regu_x_ptr = dashes;
    }
    else {
      sprintf(regu_x_buf, "%5.1f", log10(regu_x));
      regu_x_ptr = regu_x_buf;
    }
    Index ls_count = IpData().info_ls_count();
    const std::string info_string = IpData().info_string();

    if (!IpData().info_skip_output()) {
      Jnlst().Printf(J_ITERSUMMARY, J_MAIN,
                     "%4d%c%14.7e %7.2e %7.2e %5.1f %7.2e %5s %7.2e %7.2e%c%3d",
                     iter, info_iter, unscaled_f, inf_pr, inf_du, log10(mu), dnrm, regu_x_ptr,
                     alpha_dual, alpha_primal, alpha_primal_char,
                     ls_count);
      if (print_info_string_) {
        Jnlst().Printf(J_ITERSUMMARY, J_MAIN, " %s", info_string.c_str());
      }
      else {
        Jnlst().Printf(J_DETAILED, J_MAIN, " %s", info_string.c_str());
      }
      Jnlst().Printf(J_ITERSUMMARY, J_MAIN, "\n");
    }


    //////////////////////////////////////////////////////////////////////
    //           Now if desired more detail on the iterates             //
    //////////////////////////////////////////////////////////////////////

    if (Jnlst().ProduceOutput(J_DETAILED, J_MAIN)) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "\n**************************************************\n");
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "*** Beginning Iteration %d from the following point:",
                     IpData().iter_count());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "\n**************************************************\n\n");

      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Current barrier parameter mu = %21.16e\n", IpData().curr_mu());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "Current fraction-to-the-boundary parameter tau = %21.16e\n\n",
                     IpData().curr_tau());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_x||_inf   = %.16e\n", IpData().curr()->x()->Amax());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_s||_inf   = %.16e\n", IpData().curr()->s()->Amax());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_y_c||_inf = %.16e\n", IpData().curr()->y_c()->Amax());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_y_d||_inf = %.16e\n", IpData().curr()->y_d()->Amax());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_z_L||_inf = %.16e\n", IpData().curr()->z_L()->Amax());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_z_U||_inf = %.16e\n", IpData().curr()->z_U()->Amax());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_v_L||_inf = %.16e\n", IpData().curr()->v_L()->Amax());
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "||curr_v_U||_inf = %.16e\n", IpData().curr()->v_U()->Amax());
    }
    if (Jnlst().ProduceOutput(J_MOREDETAILED, J_MAIN)) {
      if (IsValid(IpData().delta())) {
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "\n||delta_x||_inf   = %.16e\n", IpData().delta()->x()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "||delta_s||_inf   = %.16e\n", IpData().delta()->s()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "||delta_y_c||_inf = %.16e\n", IpData().delta()->y_c()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "||delta_y_d||_inf = %.16e\n", IpData().delta()->y_d()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "||delta_z_L||_inf = %.16e\n", IpData().delta()->z_L()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "||delta_z_U||_inf = %.16e\n", IpData().delta()->z_U()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "||delta_v_L||_inf = %.16e\n", IpData().delta()->v_L()->Amax());
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "||delta_v_U||_inf = %.16e\n", IpData().delta()->v_U()->Amax());
      }
      else {
        Jnlst().Printf(J_MOREDETAILED, J_MAIN,
                       "\nNo search direction has been computed yet.\n");
      }
    }
    if (Jnlst().ProduceOutput(J_VECTOR, J_MAIN)) {
      IpData().curr()->x()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_x");
      IpData().curr()->s()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_s");

      IpData().curr()->y_c()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_y_c");
      IpData().curr()->y_d()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_y_d");

      IpCq().curr_slack_x_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_slack_x_L");
      IpCq().curr_slack_x_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_slack_x_U");
      IpData().curr()->z_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_z_L");
      IpData().curr()->z_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_z_U");

      IpCq().curr_slack_s_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_slack_s_L");
      IpCq().curr_slack_s_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_slack_s_U");
      IpData().curr()->v_L()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_v_L");
      IpData().curr()->v_U()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_v_U");
    }
    if (Jnlst().ProduceOutput(J_MOREVECTOR, J_MAIN)) {
      IpCq().curr_grad_lag_x()->Print(Jnlst(), J_MOREVECTOR, J_MAIN, "curr_grad_lag_x");
      IpCq().curr_grad_lag_s()->Print(Jnlst(), J_MOREVECTOR, J_MAIN, "curr_grad_lag_s");
      if (IsValid(IpData().delta())) {
        IpData().delta()->Print(Jnlst(), J_MOREVECTOR, J_MAIN, "delta");
      }
    }

    if (Jnlst().ProduceOutput(J_DETAILED, J_MAIN)) {
      Jnlst().Printf(J_DETAILED, J_MAIN,
                     "\n\n***Current NLP Values for Iteration %d:\n",
                     IpData().iter_count());
      Jnlst().Printf(J_DETAILED, J_MAIN, "\n                                   (scaled)                 (unscaled)\n");
      Jnlst().Printf(J_DETAILED, J_MAIN, "Objective...............: %24.16e  %24.16e\n", IpCq().curr_f(), IpCq().unscaled_curr_f());
      Jnlst().Printf(J_DETAILED, J_MAIN, "Dual infeasibility......: %24.16e  %24.16e\n", IpCq().curr_dual_infeasibility(NORM_MAX), IpCq().unscaled_curr_dual_infeasibility(NORM_MAX));
      Jnlst().Printf(J_DETAILED, J_MAIN, "Constraint violation....: %24.16e  %24.16e\n", IpCq().curr_nlp_constraint_violation(NORM_MAX), IpCq().unscaled_curr_nlp_constraint_violation(NORM_MAX));
      Jnlst().Printf(J_DETAILED, J_MAIN, "Complementarity.........: %24.16e  %24.16e\n", IpCq().curr_complementarity(0., NORM_MAX), IpCq().unscaled_curr_complementarity(0., NORM_MAX));
      Jnlst().Printf(J_DETAILED, J_MAIN, "Overall NLP error.......: %24.16e  %24.16e\n\n", IpCq().curr_nlp_error(), IpCq().unscaled_curr_nlp_error());
    }
    if (Jnlst().ProduceOutput(J_VECTOR, J_MAIN)) {
      IpCq().curr_grad_f()->Print(Jnlst(), J_VECTOR, J_MAIN, "grad_f");
      IpCq().curr_c()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_c");
      IpCq().curr_d()->Print(Jnlst(), J_VECTOR, J_MAIN, "curr_d");
      IpCq().curr_d_minus_s()->Print(Jnlst(), J_VECTOR, J_MAIN,
                                     "curr_d - curr_s");
    }

    if (Jnlst().ProduceOutput(J_MATRIX, J_MAIN)) {
      IpCq().curr_jac_c()->Print(Jnlst(), J_MATRIX, J_MAIN, "jac_c");
      IpCq().curr_jac_d()->Print(Jnlst(), J_MATRIX, J_MAIN, "jac_d");
      IpData().W()->Print(Jnlst(), J_MATRIX, J_MAIN, "W");
    }

    Jnlst().Printf(J_DETAILED, J_MAIN, "\n\n");
  }

} // namespace Ipopt
