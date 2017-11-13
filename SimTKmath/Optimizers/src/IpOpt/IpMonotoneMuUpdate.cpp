// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpMonotoneMuUpdate.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpMonotoneMuUpdate.hpp"
#include "IpJournalist.hpp"

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

  MonotoneMuUpdate::MonotoneMuUpdate(const SmartPtr<LineSearch>& linesearch)
      :
      MuUpdate(),
      linesearch_(linesearch),
      initialized_(false)
  {
    DBG_START_METH("MonotoneMuUpdate::MonotoneMuUpdate",
                   dbg_verbosity);
    DBG_ASSERT(IsValid(linesearch_));
  }

  MonotoneMuUpdate::~MonotoneMuUpdate()
  {
    DBG_START_METH("MonotoneMuUpdate::~MonotoneMuUpdate",
                   dbg_verbosity);
  }

  void MonotoneMuUpdate::RegisterOptions(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "mu_init", "Initial value for the barrier parameter.",
      0.0, true,
      0.1,
      "This option determines the initial value for the barrier parameter "
      "(mu).  It is only relevant in the monotone, Fiacco-McCormick "
      "version of the algorithm. (i.e., if \"mu_strategy\" is chosen "
      "as \"monotone\")");
    roptions->AddLowerBoundedNumberOption(
      "barrier_tol_factor",
      "Factor for mu in barrier stop test.",
      0.0, true,
      10.0,
      "The convergence tolerance for each barrier problem in the monotone mode "
      "is the value of the barrier parameter times \"barrier_tol_factor\". "
      "This option is also used in the adaptive mu strategy during the "
      "monotone mode. (This is kappa_epsilon in implementation paper).");
    roptions->AddBoundedNumberOption(
      "mu_linear_decrease_factor",
      "Determines linear decrease rate of barrier parameter.",
      0.0, true, 1.0, true,
      0.2,
      "For the Fiacco-McCormick update procedure the new barrier parameter mu "
      "is obtained by taking the minimum of mu*\"mu_linear_decrease_factor\" "
      "and mu^\"superlinear_decrease_power\".  (This is kappa_mu in "
      "implementation paper.) This option is also used in the adaptive mu "
      "strategy during the monotone mode.");
    roptions->AddBoundedNumberOption(
      "mu_superlinear_decrease_power",
      "Determines superlinear decrease rate of barrier parameter.",
      1.0, true, 2.0, true,
      1.5,
      "For the Fiacco-McCormick update procedure the new barrier parameter mu "
      "is obtained by taking the minimum of mu*\"mu_linear_decrease_factor\" "
      "and mu^\"superlinear_decrease_power\".  (This is theta_mu in "
      "implementation paper.) This option is also used in the adaptive mu "
      "strategy during the monotone mode.");
    roptions->AddStringOption2(
      "mu_allow_fast_monotone_decrease",
      "Allow skipping of barrier problem if barrier test is already met.",
      "yes",
      "no", "Take at least one iteration per barrier problem",
      "yes", "Allow fast decrease of mu if barrier test it met",
      "If set to \"no\", the algorithm enforces at least one iteration per "
      "barrier problem, even if the barrier test is already met for the "
      "updated barrier parameter.");
    roptions->AddBoundedNumberOption(
      "tau_min",
      "Lower bound on fraction-to-the-boundary parameter tau.",
      0.0, true, 1.0, true,
      0.99,
      "(This is tau_min in implementation paper.)  This option is also used "
      "in the adaptive mu strategy during the monotone mode.");
  }

  bool MonotoneMuUpdate::InitializeImpl(const OptionsList& options,
                                        const std::string& prefix)
  {
    options.GetNumericValue("mu_init", mu_init_, prefix);
    options.GetNumericValue("barrier_tol_factor", barrier_tol_factor_, prefix);
    options.GetNumericValue("mu_linear_decrease_factor", mu_linear_decrease_factor_, prefix);
    options.GetNumericValue("mu_superlinear_decrease_power", mu_superlinear_decrease_power_, prefix);
    options.GetBoolValue("mu_allow_fast_monotone_decrease", mu_allow_fast_monotone_decrease_, prefix);
    options.GetNumericValue("tau_min", tau_min_, prefix);
    options.GetNumericValue("compl_inf_tol", compl_inf_tol_, prefix);

    IpData().Set_mu(mu_init_);
    Number tau = Max(tau_min_, 1.0 - mu_init_);
    IpData().Set_tau(tau);

    initialized_ = false;

    //TODO we need to clean up the mu-update for the restoration phase
    if (prefix=="resto.") {
      first_iter_resto_ = true;
    }
    else {
      first_iter_resto_ = false;
    }

    return true;
  }

  bool MonotoneMuUpdate::UpdateBarrierParameter()
  {
    Number mu = IpData().curr_mu();
    Number tau = IpData().curr_tau();

    Number sub_problem_error = IpCq().curr_barrier_error();

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "Optimality Error for Barrier Sub-problem = %e\n",
                   sub_problem_error);
    Number kappa_eps_mu = barrier_tol_factor_ * mu;

    bool done = false;
    bool tiny_step_flag = IpData().tiny_step_flag();
    while ((sub_problem_error <= kappa_eps_mu || tiny_step_flag)
           && !done && !first_iter_resto_) {
      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "  sub_problem_error < kappa_eps * mu (%e)\n", kappa_eps_mu);

      // Compute the new values for mu and tau
      Number new_mu;
      Number new_tau;
      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "Updating mu=%25.16e and tau=%25.16e to ", mu, tau);
      CalcNewMuAndTau(new_mu, new_tau);
      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "new_mu=%25.16e and new_tau=%25.16e\n", new_mu, new_tau);
      bool mu_changed = (mu != new_mu);
      if (!mu_changed && tiny_step_flag) {
        THROW_EXCEPTION(TINY_STEP_DETECTED,
                        "Problem solved to best possible numerical accuracy");
      }

      // Set the new values for mu and tau
      IpData().Set_mu(new_mu);
      IpData().Set_tau(new_tau);
      mu = new_mu;
      tau = new_tau;

      // If this is the first iteration or if
      // mu_allow_fast_monotone_decrease_ is true, we want to check if
      // we can decrease mu even more
      if (initialized_ && !mu_allow_fast_monotone_decrease_) {
        done = true;
      }
      else if (!mu_changed) {
        done = true;
      }
      else {
        sub_problem_error = IpCq().curr_barrier_error();
        kappa_eps_mu = barrier_tol_factor_ * mu;
        done = (sub_problem_error > kappa_eps_mu);
      }

      // Reset the line search
      if (done && mu_changed) {
        linesearch_->Reset();
      }

      tiny_step_flag = false;
    }

    first_iter_resto_ = false;
    initialized_ = true;

    return true;
  }

  void MonotoneMuUpdate::CalcNewMuAndTau(Number &new_mu,
                                         Number &new_tau)
  {
    // update the barrier parameter
    Number mu = IpData().curr_mu();
    Number tol = IpData().tol();

    // Here we need the complementarity tolerance that is posed to the
    // scaled problem
    Number compl_inf_tol =
      IpNLP().NLP_scaling()->apply_obj_scaling(compl_inf_tol_);

    new_mu = Min( mu_linear_decrease_factor_*mu,
                  pow(mu, mu_superlinear_decrease_power_) );
    new_mu = Max(new_mu, Min(tol, compl_inf_tol)/(barrier_tol_factor_+1.));

    // update the fraction to the boundary parameter
    new_tau = Max(tau_min_, 1.-new_mu);
  }

} // namespace Ipopt
