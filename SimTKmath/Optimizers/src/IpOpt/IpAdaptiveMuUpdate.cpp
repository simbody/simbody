// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpAdaptiveMuUpdate.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpAdaptiveMuUpdate.hpp"
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

  AdaptiveMuUpdate::AdaptiveMuUpdate
  (const SmartPtr<LineSearch>& line_search,
   const SmartPtr<MuOracle>& free_mu_oracle,
   const SmartPtr<MuOracle>& fix_mu_oracle)
      :
      MuUpdate(),
      linesearch_(line_search),
      free_mu_oracle_(free_mu_oracle),
      fix_mu_oracle_(fix_mu_oracle),
      filter_(2)
  {
    DBG_ASSERT(IsValid(linesearch_));
    DBG_ASSERT(IsValid(free_mu_oracle_));
    // fix_mu_oracle may be NULL
  }

  AdaptiveMuUpdate::~AdaptiveMuUpdate()
  {}

  void AdaptiveMuUpdate::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "mu_max_fact",
      "Factor for initialization of maximum value for barrier parameter.",
      0.0, true, 1e3,
      "This option determines the upper bound on the barrier parameter.  This "
      "upper bound is computed as the average complementarity at the initial "
      "point times the value of this option. (Only used if option "
      "\"mu_strategy\" is chosen as \"adaptive\".)");
    roptions->AddLowerBoundedNumberOption(
      "mu_max",
      "Maximum value for barrier parameter.",
      0.0, true, 1e5,
      "This option specifies an upper bound on the barrier parameter in the "
      "adaptive mu selection mode.  If this option is set, it overwrites the "
      "effect of mu_max_fact. (Only used if option "
      "\"mu_strategy\" is chosen as \"adaptive\".)");
    roptions->AddLowerBoundedNumberOption(
      "mu_min",
      "Minimum value for barrier parameter.",
      0.0, true, 1e-9,
      "This option specifies the lower bound on the barrier parameter in the "
      "adaptive mu selection mode. By default, it is set to "
      "min(\"tol\",\"compl_inf_tol\")/(\"barrier_tol_factor\"+1), which "
      "should be a reasonable value. (Only used if option "
      "\"mu_strategy\" is chosen as \"adaptive\".)");
    std::string prev_cat = roptions->RegisteringCategory();
    roptions->SetRegisteringCategory("Undocumented");
    roptions->AddLowerBoundedNumberOption(
      "adaptive_mu_safeguard_factor",
      "",
      0.0, false, 0.0);
    roptions->SetRegisteringCategory(prev_cat);

    roptions->AddStringOption3(
      "adaptive_mu_globalization",
      "Globalization strategy for the adaptive mu selection mode.",
      "obj-constr-filter",
      "kkt-error", "nonmonotone decrease of kkt-error",
      "obj-constr-filter", "2-dim filter for objective and constraint violation",
      "never-monotone-mode", "disables globalization",
      "To achieve global convergence of the adaptive version, the algorithm "
      "has to switch to the monotone mode (Fiacco-McCormick approach) when "
      "convergence does not seem to appear.  This option sets the "
      "criterion used to decide when to do this switch. (Only used if option "
      "\"mu_strategy\" is chosen as \"adaptive\".)");

    roptions->AddLowerBoundedIntegerOption(
      "adaptive_mu_kkterror_red_iters",
      "Maximum number of iterations requiring sufficient progress.",
      0, 4,
      "For the \"kkt-error\" based globalization strategy, sufficient "
      "progress must be made for \"adaptive_mu_kkterror_red_iters\" "
      "iterations. If this number of iterations is exceeded, the "
      "globalization strategy switches to the monotone mode.");

    roptions->AddBoundedNumberOption(
      "adaptive_mu_kkterror_red_fact",
      "Sufficient decrease factor for \"kkt-error\" globalization strategy.",
      0.0, true, 1.0, true,
      0.9999,
      "For the \"kkt-error\" based globalization strategy, the error "
      "must decrease by this factor to be deemed sufficient decrease.");

    roptions->AddBoundedNumberOption(
      "filter_margin_fact",
      "Factor determining width of margin for obj-constr-filter adaptive globalization strategy.",
      0.0, true, 1.0, true,
      1e-5,
      "When using the adaptive globalization strategy, \"obj-constr-filter\", "
      "sufficient progress for a filter entry is defined as "
      "follows: (new obj) < (filter obj) - filter_margin_fact*(new "
      "constr-voil) OR (new constr-viol) < (filter constr-viol) - "
      "filter_margin_fact*(new constr-voil).  For the description of "
      "the \"kkt-error-filter\" option see \"filter_max_margin\".");
    roptions->AddLowerBoundedNumberOption(
      "filter_max_margin",
      "Maximum width of margin in obj-constr-filter adaptive globalization strategy.",
      0.0, true,
      1.0,
      // ToDo Detailed description later
      "");
    roptions->AddStringOption2(
      "adaptive_mu_restore_previous_iterate",
      "Indicates if the previous iterate should be restored if the monotone mode is entered.",
      "no",
      "no", "don't restore accepted iterate",
      "yes", "restore accepted iterate",
      "When the globalization strategy for the adaptive barrier algorithm "
      "switches to the monotone mode, it can either start "
      "from the most recent iterate (no), or from the last "
      "iterate that was accepted (yes).");

    roptions->AddLowerBoundedNumberOption(
      "adaptive_mu_monotone_init_factor",
      "Determines the initial value of the barrier parameter when switching to the monotone mode.",
      0.0, true, 0.8,
      "When the globalization strategy for the adaptive barrier algorithm "
      "switches to the monotone mode and fixed_mu_oracle is chosen as "
      "\"average_compl\", the barrier parameter is set to the "
      "current average complementarity times the value of "
      "\"adaptive_mu_monotone_init_factor\".");

    roptions->AddStringOption4(
      "adaptive_mu_kkt_norm_type",
      "Norm used for the KKT error in the adaptive mu globalization strategies.",
      "2-norm-squared",
      "1-norm", "use the 1-norm (abs sum)",
      "2-norm-squared", "use the 2-norm squared (sum of squares)",
      "max-norm", "use the infinity norm (max)",
      "2-norm", "use 2-norm",
      "When computing the KKT error for the globalization strategies, the "
      "norm to be used is specified with this option. Note, this options is also used "
      "in the QualityFunctionMuOracle.");

  }

  bool AdaptiveMuUpdate::InitializeImpl(const OptionsList& options,
                                        const std::string& prefix)
  {
    options.GetNumericValue("mu_max_fact", mu_max_fact_, prefix);
    if (!options.GetNumericValue("mu_max", mu_max_, prefix)) {
      // Set to a negative value as a hint that this value still has
      // to be computed
      mu_max_ = -1.;
    }
    options.GetNumericValue("tau_min", tau_min_, prefix);
    options.GetNumericValue("adaptive_mu_safeguard_factor", adaptive_mu_safeguard_factor_, prefix);
    options.GetNumericValue("adaptive_mu_kkterror_red_fact", refs_red_fact_, prefix);
    options.GetIntegerValue("adaptive_mu_kkterror_red_iters", num_refs_max_, prefix);
    Index enum_int;
    options.GetEnumValue("adaptive_mu_globalization", enum_int, prefix);
    adaptive_mu_globalization_ = AdaptiveMuGlobalizationEnum(enum_int);
    options.GetNumericValue("filter_max_margin", filter_max_margin_, prefix);
    options.GetNumericValue("filter_margin_fact", filter_margin_fact_, prefix);
    options.GetBoolValue("adaptive_mu_restore_previous_iterate", restore_accepted_iterate_, prefix);

    bool retvalue = free_mu_oracle_->Initialize(Jnlst(), IpNLP(), IpData(),
                    IpCq(), options, prefix);
    if (!retvalue) {
      return retvalue;
    }

    if (IsValid(fix_mu_oracle_)) {
      retvalue = fix_mu_oracle_->Initialize(Jnlst(), IpNLP(), IpData(),
                                            IpCq(), options, prefix);
      if (!retvalue) {
        return retvalue;
      }
    }

    options.GetNumericValue("adaptive_mu_monotone_init_factor", adaptive_mu_monotone_init_factor_, prefix);
    options.GetNumericValue("barrier_tol_factor", barrier_tol_factor_, prefix);
    options.GetNumericValue("mu_linear_decrease_factor", mu_linear_decrease_factor_, prefix);
    options.GetNumericValue("mu_superlinear_decrease_power", mu_superlinear_decrease_power_, prefix);

    options.GetEnumValue("quality_function_norm_type", enum_int, prefix);
    adaptive_mu_kkt_norm_ = QualityFunctionMuOracle::NormEnum(enum_int);
    options.GetEnumValue("quality_function_centrality", enum_int, prefix);
    adaptive_mu_kkt_centrality_ = QualityFunctionMuOracle::CentralityEnum(enum_int);
    options.GetEnumValue("quality_function_balancing_term", enum_int, prefix);
    adaptive_mu_kkt_balancing_term_ = QualityFunctionMuOracle::BalancingTermEnum(enum_int);
    options.GetNumericValue("compl_inf_tol", compl_inf_tol_, prefix);
    if (!options.GetNumericValue("mu_min", mu_min_, prefix)) {
      // Compute mu_min based on tolerance (once the NLP scaling is known)
      mu_min_default_ = true;
    }
    else {
      mu_min_default_ = false;
    }

    init_dual_inf_ = -1.;
    init_primal_inf_ = -1.;

    refs_vals_.clear();
    check_if_no_bounds_ = false;
    no_bounds_ = false;
    filter_.Clear();
    IpData().SetFreeMuMode(true);

    accepted_point_ = NULL;

    // The following lines are only here so that
    // IpoptCalculatedQuantities::CalculateSafeSlack and the first
    // output line have something to work with
    IpData().Set_mu(1.);
    IpData().Set_tau(0.);

    return retvalue;
  }

  bool AdaptiveMuUpdate::UpdateBarrierParameter()
  {
    DBG_START_METH("AdaptiveMuUpdate::UpdateBarrierParameter",
                   dbg_verbosity);

    // if min_mu_ has not been given, we now set the default (can't do
    // that earlier, because during call of InitializeImpl, the
    // scaling in the NLP is not yet determined).  We compute this
    // here in every iteration, since the tolerance might be changed
    // (e.g. in the restoration phase)
    if (mu_min_default_) {
      mu_min_ = Min(IpData().tol(),
                    IpNLP().NLP_scaling()->apply_obj_scaling(compl_inf_tol_))/
                (barrier_tol_factor_+1);
    }

    // if mu_max has not yet been computed, do so now, based on the
    // current average complementarity
    if (mu_max_<0.) {
      mu_max_ = mu_max_fact_*IpCq().curr_avrg_compl();
      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "Setting mu_max to %e.\n", mu_max_);
    }

    // if there are not bounds, we always return the minimum MU value
    if (!check_if_no_bounds_) {
      Index n_bounds = IpData().curr()->z_L()->Dim() + IpData().curr()->z_U()->Dim()
                       + IpData().curr()->v_L()->Dim() + IpData().curr()->v_U()->Dim();

      if (n_bounds==0) {
        no_bounds_ = true;
        IpData().Set_mu(mu_min_);
        IpData().Set_tau(tau_min_);
      }

      check_if_no_bounds_ = true;
    }

    if (no_bounds_)
      return true;

    bool tiny_step_flag = IpData().tiny_step_flag();
    if (!IpData().FreeMuMode()) {
      // if we are in the fixed mu mode, we need to check if the
      // current iterate is good enough to continue with the free mode
      bool sufficient_progress = CheckSufficientProgress();
      if (sufficient_progress && !tiny_step_flag) {
        Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                       "Switching back to free mu mode.\n");
        IpData().SetFreeMuMode(true);
        // Skipping Restoration phase?
        RememberCurrentPointAsAccepted();
      }
      else {
        Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                       "Remaining in fixed mu mode.\n");

        // ToDo decide whether we want this for all options
        Number sub_problem_error = IpCq().curr_barrier_error();
        Number mu = IpData().curr_mu();
        if (sub_problem_error <= barrier_tol_factor_ * mu ||
            tiny_step_flag) {
          // If the current barrier problem has been solved sufficiently
          // well, decrease mu
          // ToDo combine this code with MonotoneMuUpdate
          Number tol = IpData().tol();
          Number compl_inf_tol =
            IpNLP().NLP_scaling()->apply_obj_scaling(compl_inf_tol_);

          Number new_mu = Min( mu_linear_decrease_factor_*mu,
                               pow(mu, mu_superlinear_decrease_power_) );
          DBG_PRINT((1,"new_mu = %e, compl_inf_tol = %e tol = %e\n", new_mu, compl_inf_tol, tol));
          new_mu = Max(new_mu,
                       Min(compl_inf_tol, tol)/(barrier_tol_factor_+1));
          if (tiny_step_flag && new_mu == mu) {
            THROW_EXCEPTION(TINY_STEP_DETECTED,
                            "Problem solved to best possible numerical accuracy");
          }
          Number new_tau = Compute_tau_monotone(mu);
          IpData().Set_mu(new_mu);
          IpData().Set_tau(new_tau);
          Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                         "Reducing mu to %24.16e in fixed mu mode. Tau becomes %24.16e\n", new_mu, new_tau);
          linesearch_->Reset();
        }
      }
    }
    else {
      // Here we are in the free mu mode.
      bool sufficient_progress = CheckSufficientProgress();
      if (linesearch_->CheckSkippedLineSearch() || tiny_step_flag ) {
        sufficient_progress = false;
      }
      if (sufficient_progress) {
        Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                       "Staying in free mu mode.\n");
        RememberCurrentPointAsAccepted();
      }
      else {
        IpData().SetFreeMuMode(false);

        if (restore_accepted_iterate_) {
          // Restore most recent accepted iterate to start fixed mode from
          Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                         "Restoring most recent accepted point.\n");
          SmartPtr<IteratesVector> prev_iter = accepted_point_->MakeNewContainer();
          IpData().set_trial(prev_iter);
          IpData().AcceptTrialPoint();
        }

        // Set the new values for mu and tau and tell the linesearch
        // to reset its memory
        Number mu = NewFixedMu();
        Number tau = Compute_tau_monotone(mu);

        if (tiny_step_flag && mu==IpData().curr_mu()) {
          THROW_EXCEPTION(TINY_STEP_DETECTED,
                          "Problem solved to best possible numerical accuracy");
        }

        IpData().Set_mu(mu);
        IpData().Set_tau(tau);
        Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                       "Switching to fixed mu mode with mu = %24.16e and tau = %24.16e.\n", mu, tau);
        linesearch_->Reset();
        // Skipping Restoration phase?
      }
    }

    if (IpData().FreeMuMode()) {

      // Choose the fraction-to-the-boundary parameter for the current
      // iteration
      // ToDo: Is curr_nlp_error really what we should use here?
      Number tau = Max(tau_min_, 1-IpCq().curr_nlp_error());
      IpData().Set_tau(tau);

      // Compute the new barrier parameter via the oracle
      Number mu;
      bool retval = free_mu_oracle_->CalculateMu(mu_min_, mu_max_, mu);
      if (!retval) {
        Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                       "The mu oracle could not compute a new value of the barrier parameter.\n");
        return false;
      }

      mu = Max(mu, mu_min_);
      Number mu_lower_safe = lower_mu_safeguard();
      if (mu < mu_lower_safe) {
        Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                       "mu = %e smaller than safeguard = %e. Increasing mu.\n",
                       mu, mu_lower_safe);
        mu = mu_lower_safe;
        IpData().Append_info_string("m");
      }

      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "Barrier parameter mu computed by oracle is %e\n",
                     mu);

      // Apply safeguards if appropriate
      mu = Min(mu, mu_max_);
      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "Barrier parameter mu after safeguards is %e\n",
                     mu);

      // Set the new values
      IpData().Set_mu(mu);

      linesearch_->Reset();
      // Uncomment the next line if the line search should not switch to
      // the restoration phase in the free mode

      // linesearch_->SetRigorousLineSearch(false);
    }
    else {
      IpData().Append_info_string("F");
      linesearch_->SetRigorousLineSearch(true);
    }

    return true;
  }

  bool
  AdaptiveMuUpdate::CheckSufficientProgress()
  {
    bool retval = true;

    switch (adaptive_mu_globalization_) {
      case KKT_ERROR : {
        Index num_refs = (Index)refs_vals_.size();
        if (num_refs >= num_refs_max_) {
          retval = false;
          Number curr_error = quality_function_pd_system();
          std::list<Number>::iterator iter;
          for (iter = refs_vals_.begin(); iter != refs_vals_.end();
               ++iter) {
            if ( curr_error <= refs_red_fact_*(*iter) ) {
              retval = true;
            }
          }
        }
      }
      break;
      case FILTER_OBJ_CONSTR : {
        /*
               retval = filter_.Acceptable(IpCq().curr_f(),
                                           IpCq().curr_constraint_violation());
        */
        // ToDo: Is curr_nlp_error really what we should use here?
        Number curr_error = IpCq().curr_nlp_error();
        Number margin = filter_margin_fact_*Min(filter_max_margin_, curr_error);
        retval = filter_.Acceptable(IpCq().curr_f() + margin,
                                    IpCq().curr_constraint_violation() + margin);
      }
      break;
      case FILTER_KKT_ERROR : {
        DBG_ASSERT(false && "Unknown adaptive_mu_globalization value.");
      }
      break;
      case NEVER_MONOTONE_MODE :
      retval = true;
      break;
      default:
      DBG_ASSERT(false && "Unknown adaptive_mu_globalization value.");
    }

    return retval;
  }

  void
  AdaptiveMuUpdate::RememberCurrentPointAsAccepted()
  {
    switch (adaptive_mu_globalization_) {
      case KKT_ERROR : {
        Number curr_error = quality_function_pd_system();
        Index num_refs = (Index)refs_vals_.size();
        if (num_refs >= num_refs_max_) {
          refs_vals_.pop_front();
        }
        refs_vals_.push_back(curr_error);

        if (Jnlst().ProduceOutput(J_MOREDETAILED, J_BARRIER_UPDATE)) {
          Index num_refs = 0;
          std::list<Number>::iterator iter;
          for (iter = refs_vals_.begin(); iter != refs_vals_.end();
               ++iter) {
            num_refs++;
            Jnlst().Printf(J_MOREDETAILED, J_BARRIER_UPDATE,
                           "pd system reference[%2d] = %.6e\n", num_refs, *iter);
          }
        }
      }
      break;
      case FILTER_OBJ_CONSTR : {
        /*
               Number theta = IpCq().curr_constraint_violation();
               filter_.AddEntry(IpCq().curr_f() - filter_margin_fact_*theta,
                                IpCq().curr_constraint_violation() - filter_margin_fact_*theta,
                                IpData().iter_count());
               filter_.Print(Jnlst());
        */
        filter_.AddEntry(IpCq().curr_f(),
                         IpCq().curr_constraint_violation(),
                         IpData().iter_count());
        filter_.Print(Jnlst());
      }
      break;
      case FILTER_KKT_ERROR : {
        DBG_ASSERT(false && "Unknown corrector_type value.");
      }
      break;
      default:
      DBG_ASSERT(false && "Unknown corrector_type value.");
    }

    if (restore_accepted_iterate_) {
      // Keep pointers to this iterate so that it could be restored
      accepted_point_ = IpData().curr();
    }
  }

  Number
  AdaptiveMuUpdate::Compute_tau_monotone(Number mu)
  {
    return Max(tau_min_, 1-mu);
  }

  Number
  AdaptiveMuUpdate::min_ref_val()
  {
    DBG_ASSERT(adaptive_mu_globalization_==KKT_ERROR);
    Number min_ref;
    DBG_ASSERT(refs_vals_.size()>0);
    std::list<Number>::iterator iter = refs_vals_.begin();
    min_ref = *iter;
    ++iter;
    while (iter != refs_vals_.end()) {
      min_ref = Min(min_ref, *iter);
      ++iter;
    }
    return min_ref;
  }

  Number
  AdaptiveMuUpdate::max_ref_val()
  {
    DBG_ASSERT(adaptive_mu_globalization_==KKT_ERROR);
    Number max_ref;
    DBG_ASSERT(refs_vals_.size()>0);
    std::list<Number>::iterator iter = refs_vals_.begin();
    max_ref = *iter;
    ++iter;
    while (iter != refs_vals_.end()) {
      max_ref = Max(max_ref, *iter);
      ++iter;
    }
    return max_ref;
  }

  Number
  AdaptiveMuUpdate::NewFixedMu()
  {
    Number max_ref;
    // ToDo: Decide whether we should impose an upper bound on
    // mu based on the smallest reference value.  For now, don't
    // impose one.
    max_ref = 1e20;
    /*
    switch (adaptive_mu_globalization_) {
      case 1 :
      max_ref = max_ref_val();
      break;
      case 2 : {
        max_ref = 1e20;
      }
      break;
      default:
      DBG_ASSERT("Unknown corrector_type value.");
    }
    */

    Number new_mu;
    bool have_mu = false;
    ;
    if (IsValid(fix_mu_oracle_)) {
      have_mu = fix_mu_oracle_->CalculateMu(mu_min_, mu_max_, new_mu);
      if (!have_mu) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "New fixed value for mu could not be computed from the mu_oracle.\n");
      }
    }
    if (!have_mu) {
      new_mu = adaptive_mu_monotone_init_factor_*IpCq().curr_avrg_compl();
    }
    new_mu = Max(new_mu, lower_mu_safeguard());
    new_mu = Min(new_mu, Number(0.1) * max_ref);

    new_mu = Max(new_mu, mu_min_);
    new_mu = Min(new_mu, mu_max_);

    return new_mu;
  }

  Number
  AdaptiveMuUpdate::quality_function_pd_system()
  {
    Index n_dual = IpData().curr()->x()->Dim() + IpData().curr()->s()->Dim();
    Index n_pri = IpData().curr()->y_c()->Dim() + IpData().curr()->y_d()->Dim();
    Index n_comp = IpData().curr()->z_L()->Dim() + IpData().curr()->z_U()->Dim() +
                   IpData().curr()->v_L()->Dim() + IpData().curr()->v_U()->Dim();

    Number dual_inf=0.;
    Number primal_inf=0.;
    Number complty=0.;
    switch (adaptive_mu_kkt_norm_) {
      case QualityFunctionMuOracle::NM_NORM_1:
      dual_inf =
        IpCq().curr_dual_infeasibility(NORM_1);
      primal_inf =
        IpCq().curr_primal_infeasibility(NORM_1);
      complty =
        IpCq().curr_complementarity(0., NORM_1);
      dual_inf /= (Number)n_dual;
      DBG_ASSERT(n_pri>0 || primal_inf==0.);
      if (n_pri>0) {
        primal_inf /= (Number)n_pri;
      }
      DBG_ASSERT(n_comp>0 || complty==0.);
      if (n_comp>0) {
        complty /= (Number)n_comp;
      }
      break;
      case QualityFunctionMuOracle::NM_NORM_2_SQUARED:
      dual_inf =
        IpCq().curr_dual_infeasibility(NORM_2);
      dual_inf *= dual_inf;
      primal_inf =
        IpCq().curr_primal_infeasibility(NORM_2);
      primal_inf *= primal_inf;
      complty =
        IpCq().curr_complementarity(0., NORM_2);
      complty *= complty;
      dual_inf /= (Number)n_dual;
      DBG_ASSERT(n_pri>0 || primal_inf==0.);
      if (n_pri>0) {
        primal_inf /= (Number)n_pri;
      }
      DBG_ASSERT(n_comp>0 || complty==0.);
      if (n_comp>0) {
        complty /= (Number)n_comp;
      }
      break;
      case QualityFunctionMuOracle::NM_NORM_MAX:
      dual_inf =
        IpCq().curr_dual_infeasibility(NORM_MAX);
      primal_inf =
        IpCq().curr_primal_infeasibility(NORM_MAX);
      complty =
        IpCq().curr_complementarity(0., NORM_MAX);
      break;
      case QualityFunctionMuOracle::NM_NORM_2:
      dual_inf =
        IpCq().curr_dual_infeasibility(NORM_2);
      primal_inf =
        IpCq().curr_primal_infeasibility(NORM_2);
      complty =
        IpCq().curr_complementarity(0., NORM_2);
      dual_inf /= sqrt((Number)n_dual);
      DBG_ASSERT(n_pri>0 || primal_inf==0.);
      if (n_pri>0) {
        primal_inf /= sqrt((Number)n_pri);
      }
      DBG_ASSERT(n_comp>0 || complty==0.);
      if (n_comp>0) {
        complty /= sqrt((Number)n_comp);
      }
      break;
    }

    Number centrality = 0.;
    if (adaptive_mu_kkt_centrality_!=0) {
      Number xi = IpCq().curr_centrality_measure();
      switch (adaptive_mu_kkt_centrality_) {
        case 1:
        centrality = -complty*log(xi);
        break;
        case 2:
        centrality = complty/xi;
        break;
        case 3:
        centrality = complty/pow(xi,3);
        break;
        default:
        DBG_ASSERT(false && "Unknown value for adaptive_mu_kkt_centrality_");
      }
    }

    Number balancing_term=0.;
    switch (adaptive_mu_kkt_balancing_term_) {
      case 0:
      //Nothing
      break;
      case 1:
      balancing_term = pow(Max(0., Max(dual_inf,primal_inf)-complty),3);
      break;
      default:
      DBG_ASSERT(false && "Unknown value for adaptive_mu_kkt_balancing_term");
    }

    DBG_ASSERT(centrality>=0.);
    DBG_ASSERT(balancing_term>=0);
    Number kkt_error = primal_inf + dual_inf + complty +
                       centrality + balancing_term;

    Jnlst().Printf(J_MOREDETAILED, J_BARRIER_UPDATE,
                   "KKT error in barrier update check:\n"
                   "  primal infeasibility: %15.6e\n"
                   "    dual infeasibility: %15.6e\n"
                   "       complementarity: %15.6e\n"
                   "            centrality: %15.6e\n"
                   "             kkt error: %15.6e\n",
                   primal_inf, dual_inf, complty, centrality, kkt_error);

    return kkt_error;
  }

  Number
  AdaptiveMuUpdate::lower_mu_safeguard()
  {
    DBG_START_METH("AdaptiveMuUpdate::lower_mu_safeguard",
                   dbg_verbosity);
    if (adaptive_mu_safeguard_factor_ == 0.)
      return 0.;

    Number dual_inf =
      IpCq().curr_dual_infeasibility(NORM_1);
    Number primal_inf =
      IpCq().curr_primal_infeasibility(NORM_1);
    Index n_dual = IpData().curr()->x()->Dim() + IpData().curr()->s()->Dim();
    dual_inf /= (Number)n_dual;
    Index n_pri = IpData().curr()->y_c()->Dim() + IpData().curr()->y_d()->Dim();
    DBG_ASSERT(n_pri>0 || primal_inf==0.);
    if (n_pri>0) {
      primal_inf /= (Number)n_pri;
    }

    if (init_dual_inf_ < 0.) {
      init_dual_inf_ = Max(1., dual_inf);
    }
    if (init_primal_inf_ < 0.) {
      init_primal_inf_ = Max(1., primal_inf);
    }

    Number lower_mu_safeguard =
      Max(adaptive_mu_safeguard_factor_ * (dual_inf/init_dual_inf_),
          adaptive_mu_safeguard_factor_ * (primal_inf/init_primal_inf_));
    DBG_PRINT((1,"dual_inf=%e init_dual_inf_=%e primal_inf=%e init_primal_inf_=%e\n", dual_inf, init_dual_inf_, primal_inf, init_primal_inf_));

    if (adaptive_mu_globalization_==KKT_ERROR) {
      lower_mu_safeguard = Min(lower_mu_safeguard, min_ref_val());
    }

    return lower_mu_safeguard;
  }

} // namespace Ipopt
