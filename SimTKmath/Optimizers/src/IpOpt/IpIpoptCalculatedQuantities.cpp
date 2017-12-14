// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpIpoptCalculatedQuantities.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpIpoptCalculatedQuantities.hpp"
#include "IpSumSymMatrix.hpp"
#include "IpLowRankUpdateSymMatrix.hpp"
#include "IpRestoIpoptNLP.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

#include <limits>

namespace SimTKIpopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  IpoptCalculatedQuantities::IpoptCalculatedQuantities
  (const SmartPtr<IpoptNLP>& ip_nlp,
   const SmartPtr<IpoptData>& ip_data)
      :
      ip_nlp_(ip_nlp),
      ip_data_(ip_data),

      curr_slack_x_L_cache_(1),
      curr_slack_x_U_cache_(1),
      curr_slack_s_L_cache_(1),
      curr_slack_s_U_cache_(1),
      trial_slack_x_L_cache_(1),
      trial_slack_x_U_cache_(1),
      trial_slack_s_L_cache_(1),
      trial_slack_s_U_cache_(1),
      num_adjusted_slack_x_L_(0),
      num_adjusted_slack_x_U_(0),
      num_adjusted_slack_s_L_(0),
      num_adjusted_slack_s_U_(0),

      curr_f_cache_(2),
      trial_f_cache_(5),
      curr_grad_f_cache_(2),
      trial_grad_f_cache_(1),

      curr_barrier_obj_cache_(2),
      trial_barrier_obj_cache_(5),
      curr_grad_barrier_obj_x_cache_(1),
      curr_grad_barrier_obj_s_cache_(1),
      grad_kappa_times_damping_x_cache_(1),
      grad_kappa_times_damping_s_cache_(1),

      curr_c_cache_(1),
      trial_c_cache_(2),
      curr_d_cache_(1),
      trial_d_cache_(2),
      curr_d_minus_s_cache_(1),
      trial_d_minus_s_cache_(1),
      curr_jac_c_cache_(1),
      trial_jac_c_cache_(1),
      curr_jac_d_cache_(1),
      trial_jac_d_cache_(1),
      curr_jac_cT_times_vec_cache_(2),
      trial_jac_cT_times_vec_cache_(1),
      curr_jac_dT_times_vec_cache_(2),
      trial_jac_dT_times_vec_cache_(1),
      curr_jac_c_times_vec_cache_(1),
      curr_jac_d_times_vec_cache_(1),
      curr_constraint_violation_cache_(2),
      trial_constraint_violation_cache_(5),
      curr_nlp_constraint_violation_cache_(3),
      unscaled_curr_nlp_constraint_violation_cache_(3),

      curr_exact_hessian_cache_(1),

      curr_grad_lag_x_cache_(1),
      trial_grad_lag_x_cache_(1),
      curr_grad_lag_s_cache_(1),
      trial_grad_lag_s_cache_(1),
      curr_grad_lag_with_damping_x_cache_(0),
      curr_grad_lag_with_damping_s_cache_(0),
      curr_compl_x_L_cache_(1),
      curr_compl_x_U_cache_(1),
      curr_compl_s_L_cache_(1),
      curr_compl_s_U_cache_(1),
      trial_compl_x_L_cache_(1),
      trial_compl_x_U_cache_(1),
      trial_compl_s_L_cache_(1),
      trial_compl_s_U_cache_(1),
      curr_relaxed_compl_x_L_cache_(1),
      curr_relaxed_compl_x_U_cache_(1),
      curr_relaxed_compl_s_L_cache_(1),
      curr_relaxed_compl_s_U_cache_(1),
      curr_primal_infeasibility_cache_(3),
      trial_primal_infeasibility_cache_(3),
      curr_dual_infeasibility_cache_(3),
      trial_dual_infeasibility_cache_(3),
      unscaled_curr_dual_infeasibility_cache_(3),
      curr_complementarity_cache_(6),
      trial_complementarity_cache_(6),
      curr_centrality_measure_cache_(1),
      curr_nlp_error_cache_(1),
      unscaled_curr_nlp_error_cache_(1),
      curr_barrier_error_cache_(1),
      curr_primal_dual_system_error_cache_(1),
      trial_primal_dual_system_error_cache_(3),

      primal_frac_to_the_bound_cache_(5),
      dual_frac_to_the_bound_cache_(5),

      curr_sigma_x_cache_(1),
      curr_sigma_s_cache_(1),

      curr_avrg_compl_cache_(1),
      trial_avrg_compl_cache_(1),
      curr_gradBarrTDelta_cache_(1),

      dampind_x_L_(NULL),
      dampind_x_U_(NULL),
      dampind_s_L_(NULL),
      dampind_s_U_(NULL),

      initialize_called_(false)
  {
    DBG_START_METH("IpoptCalculatedQuantities::IpoptCalculatedQuantities",
                   dbg_verbosity);
    DBG_ASSERT(IsValid(ip_nlp_) && IsValid(ip_data_));
  }

  IpoptCalculatedQuantities::~IpoptCalculatedQuantities()
  {}

  void IpoptCalculatedQuantities::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Convergence");
    roptions->AddLowerBoundedNumberOption(
      "s_max",
      "Scaling threshold for the NLP error.",
      0.0, true, 100.0,
      "(See paragraph after Eqn. (6) in the implementation paper.)");

    roptions->SetRegisteringCategory("NLP");
    roptions->AddLowerBoundedNumberOption(
      "kappa_d",
      "Weight for linear damping term (to handle one-sided bounds).",
      0.0, false, 1e-5,
      "(see Section 3.7 in implementation paper.)");

    roptions->SetRegisteringCategory("Line Search");
    roptions->AddLowerBoundedNumberOption(
      "slack_move",
      "Correction size for very small slacks.",
      0.0, false,
      pow(std::numeric_limits<double>::epsilon(), 0.75),
      "Due to numerical issues or the lack of an interior, the slack variables might "
      "become very small.  If a slack becomes very small compared to machine "
      "precision, the corresponding bound is moved slightly.  This parameter "
      "determines how large the move should be.  Its default value is "
      "mach_eps^{3/4}.  (See also end of Section 3.5 in implementation paper "
      "- but actual implementation might be somewhat different.)");
    roptions->SetRegisteringCategory("Line search");
    roptions->AddStringOption3(
      "constraint_violation_norm_type",
      "Norm to be used for the constraint violation in the line search.",
      "1-norm",
      "1-norm", "use the 1-norm",
      "2-norm", "use the 2-norm",
      "max-norm", "use the infinity norm",
      "Determines which norm should be used when the algorithm computes the "
      "constraint violation in the line search.");
  }

  bool IpoptCalculatedQuantities::Initialize(const Journalist& jnlst,
      const OptionsList& options,
      const std::string& prefix)
  {
    std::string svalue;
    Index enum_int;

    options.GetNumericValue("s_max", s_max_, prefix);
    options.GetNumericValue("kappa_d", kappa_d_, prefix);
    options.GetNumericValue("slack_move", slack_move_, prefix);
    options.GetEnumValue("constraint_violation_norm_type", enum_int, prefix);
    constr_viol_normtype_ = ENormType(enum_int);
    // The following option is registered by OrigIpoptNLP
    options.GetBoolValue("warm_start_same_structure",
                         warm_start_same_structure_, prefix);

    if (!warm_start_same_structure_) {
      dampind_x_L_ = NULL;
      dampind_x_U_ = NULL;
      dampind_s_L_ = NULL;
      dampind_s_U_ = NULL;

      tmp_x_ = NULL;
      tmp_s_ = NULL;
      tmp_c_ = NULL;
      tmp_d_ = NULL;
      tmp_x_L_ = NULL;
      tmp_x_U_ = NULL;
      tmp_s_L_ = NULL;
      tmp_s_U_ = NULL;
    }

    num_adjusted_slack_x_L_ = 0;
    num_adjusted_slack_x_U_ = 0;
    num_adjusted_slack_s_L_ = 0;
    num_adjusted_slack_s_U_ = 0;

    initialize_called_ = true;
    return true;
  }

  ///////////////////////////////////////////////////////////////////////////
  //                         Slack Calculations                            //
  ///////////////////////////////////////////////////////////////////////////

  SmartPtr<Vector>
  IpoptCalculatedQuantities::CalcSlack_L(const Matrix& P,
                                         const Vector& x,
                                         const Vector& x_bound)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalcSlack_L",
                   dbg_verbosity);
    SmartPtr<Vector> result;
    result = x_bound.MakeNew();
    result->Copy(x_bound);
    P.TransMultVector(1.0, x, -1.0, *result);
    return result;
  }

  SmartPtr<Vector>
  IpoptCalculatedQuantities::CalcSlack_U(const Matrix& P,
                                         const Vector& x,
                                         const Vector& x_bound)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalcSlack_U",
                   dbg_verbosity);
    SmartPtr<Vector> result;
    result = x_bound.MakeNew();
    result->Copy(x_bound);
    P.TransMultVector(-1.0, x, 1.0, *result);
    return result;
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::curr_slack_x_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_slack_x_L()",
                   dbg_verbosity);
    SmartPtr<Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> x_bound = ip_nlp_->x_L();
    if (!curr_slack_x_L_cache_.GetCachedResult1Dep(result, *x)) {
      if (!trial_slack_x_L_cache_.GetCachedResult1Dep(result, *x)) {
        SmartPtr<const Matrix> P = ip_nlp_->Px_L();
        DBG_PRINT_VECTOR(2,"x_L", *x_bound);
        result = CalcSlack_L(*P, *x, *x_bound);
        DBG_ASSERT(num_adjusted_slack_x_L_==0);
        num_adjusted_slack_x_L_ =
          CalculateSafeSlack(result, x_bound, x, ip_data_->curr()->z_L());
      }
      curr_slack_x_L_cache_.AddCachedResult1Dep(result, *x);
    }
    return ConstPtr(result);
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::curr_slack_x_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_slack_x_U()",
                   dbg_verbosity);

    SmartPtr<Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> x_bound = ip_nlp_->x_U();
    if (!curr_slack_x_U_cache_.GetCachedResult1Dep(result, *x)) {
      if (!trial_slack_x_U_cache_.GetCachedResult1Dep(result, *x)) {
        SmartPtr<const Matrix> P = ip_nlp_->Px_U();
        result = CalcSlack_U(*P, *x, *x_bound);
        DBG_ASSERT(num_adjusted_slack_x_U_==0);
        num_adjusted_slack_x_U_ =
          CalculateSafeSlack(result, x_bound, x, ip_data_->curr()->z_U());
      }
      curr_slack_x_U_cache_.AddCachedResult1Dep(result, *x);
    }
    return ConstPtr(result);
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::curr_slack_s_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_slack_s_L()",
                   dbg_verbosity);

    SmartPtr<Vector> result;
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    SmartPtr<const Vector> s_bound = ip_nlp_->d_L();
    if (!curr_slack_s_L_cache_.GetCachedResult1Dep(result, *s)) {
      if (!trial_slack_s_L_cache_.GetCachedResult1Dep(result, *s)) {
        SmartPtr<const Matrix> P = ip_nlp_->Pd_L();
        result = CalcSlack_L(*P, *s, *s_bound);
        DBG_ASSERT(num_adjusted_slack_s_L_==0);
        num_adjusted_slack_s_L_ =
          CalculateSafeSlack(result, s_bound, s, ip_data_->curr()->v_L());
      }
      curr_slack_s_L_cache_.AddCachedResult1Dep(result, *s);
    }
    return ConstPtr(result);
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::curr_slack_s_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_slack_s_U()",
                   dbg_verbosity);

    SmartPtr<Vector> result;
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    SmartPtr<const Vector> s_bound = ip_nlp_->d_U();
    if (!curr_slack_s_U_cache_.GetCachedResult1Dep(result, *s)) {
      if (!trial_slack_s_U_cache_.GetCachedResult1Dep(result, *s)) {
        SmartPtr<const Matrix> P = ip_nlp_->Pd_U();
        result = CalcSlack_U(*P, *s, *s_bound);
        DBG_ASSERT(num_adjusted_slack_s_U_==0);
        num_adjusted_slack_s_U_ =
          CalculateSafeSlack(result, s_bound, s, ip_data_->curr()->v_U());
        DBG_PRINT_VECTOR(2, "result", *result);
        DBG_PRINT((1, "num_adjusted_slack_s_U = %d\n", num_adjusted_slack_s_U_));
      }
      curr_slack_s_U_cache_.AddCachedResult1Dep(result, *s);
    }
    return ConstPtr(result);
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::trial_slack_x_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_slack_x_L()",
                   dbg_verbosity);

    num_adjusted_slack_x_L_ = 0;
    SmartPtr<Vector> result;
    SmartPtr<const Vector> x = ip_data_->trial()->x();
    SmartPtr<const Vector> x_bound = ip_nlp_->x_L();
    if (!trial_slack_x_L_cache_.GetCachedResult1Dep(result, *x)) {
      if (!curr_slack_x_L_cache_.GetCachedResult1Dep(result, *x)) {
        SmartPtr<const Matrix> P = ip_nlp_->Px_L();
        result = CalcSlack_L(*P, *x, *x_bound);
        DBG_ASSERT(num_adjusted_slack_x_L_==0);
        num_adjusted_slack_x_L_ =
          CalculateSafeSlack(result, x_bound, x, ip_data_->curr()->z_L());
      }
      trial_slack_x_L_cache_.AddCachedResult1Dep(result, *x);
    }
    return ConstPtr(result);
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::trial_slack_x_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_slack_x_U()",
                   dbg_verbosity);

    num_adjusted_slack_x_U_ = 0;
    SmartPtr<Vector> result;
    SmartPtr<const Vector> x = ip_data_->trial()->x();
    SmartPtr<const Vector> x_bound = ip_nlp_->x_U();
    if (!trial_slack_x_U_cache_.GetCachedResult1Dep(result, *x)) {
      if (!curr_slack_x_U_cache_.GetCachedResult1Dep(result, *x)) {
        SmartPtr<const Matrix> P = ip_nlp_->Px_U();
        result = CalcSlack_U(*P, *x, *x_bound);
        DBG_ASSERT(num_adjusted_slack_x_U_==0);
        num_adjusted_slack_x_U_ =
          CalculateSafeSlack(result, x_bound, x, ip_data_->curr()->z_U());
      }
      trial_slack_x_U_cache_.AddCachedResult1Dep(result, *x);
    }
    return ConstPtr(result);
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::trial_slack_s_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_slack_s_L()",
                   dbg_verbosity);

    num_adjusted_slack_s_L_ = 0;
    SmartPtr<Vector> result;
    SmartPtr<const Vector> s = ip_data_->trial()->s();
    SmartPtr<const Vector> s_bound = ip_nlp_->d_L();
    if (!trial_slack_s_L_cache_.GetCachedResult1Dep(result, *s)) {
      if (!curr_slack_s_L_cache_.GetCachedResult1Dep(result, *s)) {
        SmartPtr<const Matrix> P = ip_nlp_->Pd_L();
        result = CalcSlack_L(*P, *s, *s_bound);
        DBG_ASSERT(num_adjusted_slack_s_L_==0);
        num_adjusted_slack_s_L_ =
          CalculateSafeSlack(result, s_bound, s, ip_data_->curr()->v_L());
      }
      trial_slack_s_L_cache_.AddCachedResult1Dep(result, *s);
    }
    return ConstPtr(result);
  }

  SmartPtr<const Vector> IpoptCalculatedQuantities::trial_slack_s_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_slack_s_U()",
                   dbg_verbosity);

    num_adjusted_slack_s_U_ = 0;
    SmartPtr<Vector> result;
    SmartPtr<const Vector> s = ip_data_->trial()->s();
    SmartPtr<const Vector> s_bound = ip_nlp_->d_U();
    if (!trial_slack_s_U_cache_.GetCachedResult1Dep(result, *s)) {
      if (!curr_slack_s_U_cache_.GetCachedResult1Dep(result, *s)) {
        SmartPtr<const Matrix> P = ip_nlp_->Pd_U();
        result = CalcSlack_U(*P, *s, *s_bound);
        DBG_ASSERT(num_adjusted_slack_s_U_==0);
        num_adjusted_slack_s_U_ =
          CalculateSafeSlack(result, s_bound, s, ip_data_->curr()->v_U());
        DBG_PRINT((1, "num_adjusted_slack_s_U = %d\n", num_adjusted_slack_s_U_));
        DBG_PRINT_VECTOR(2, "trial_slack_s_U", *result);
      }
      trial_slack_s_U_cache_.AddCachedResult1Dep(result, *s);
    }
    return ConstPtr(result);
  }

  Index IpoptCalculatedQuantities::
  CalculateSafeSlack(SmartPtr<Vector>& slack,
                     const SmartPtr<const Vector>& bound,
                     const SmartPtr<const Vector>& curr_point,
                     const SmartPtr<const Vector>& multiplier)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalculateSafeSlack", dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    Index retval = 0;
    if (slack->Dim() > 0) {
      Number min_slack = slack->Min();
      // TODO we need to make sure that this also works for non-monotone MUs
      Number s_min = std::numeric_limits<Number>::epsilon()
                     * Min(1., ip_data_->curr_mu());
      DBG_PRINT((1,"s_min = %g, min_slack=%g\n", s_min, min_slack));
      if (min_slack < s_min) {
        // Need to correct the slacks and calculate new bounds...
        SmartPtr<Vector> t = slack->MakeNew();
        t->Copy(*slack);
        t->AddScalar(-s_min);
        t->ElementWiseSgn();

        SmartPtr<Vector> zero_vec = t->MakeNew();
        zero_vec->Set(0.0);
        t->ElementWiseMin(*zero_vec);
        t->Scal(-1.0);
        retval = (Index)t->Asum();
        DBG_PRINT((1,"Number of slack corrections = %d\n", retval));
        DBG_PRINT_VECTOR(2, "t(sgn)", *t);

        // ToDo AW: I added the follwing line b/c I found a case where
        // slack was negative and this correction produced 0
        slack->ElementWiseMax(*zero_vec);

        SmartPtr<Vector> t2 = t->MakeNew();
        t2->Set(ip_data_->curr_mu());
        t2->ElementWiseDivide(*multiplier);

        SmartPtr<Vector> s_min_vec = t2->MakeNew();
        s_min_vec->Set(s_min);

        t2->ElementWiseMax(*s_min_vec);
        t2->Axpy(-1.0, *slack);
        DBG_PRINT_VECTOR(2, "tw(smin,mu/mult)", *t2);

        t->ElementWiseMultiply(*t2);
        t->Axpy(1.0, *slack);

        SmartPtr<Vector> t_max = t2;
        t_max->Set(1.0);
        SmartPtr<Vector> abs_bound = bound->MakeNew();
        abs_bound->Copy(*bound);
        abs_bound->ElementWiseAbs();
        t_max->ElementWiseMax(*abs_bound);
        DBG_PRINT_VECTOR(2, "t_max1", *t_max);
        DBG_PRINT_VECTOR(2, "slack", *slack);
        t_max->AddOneVector(1.0, *slack, slack_move_);
        DBG_PRINT_VECTOR(2, "t_max2", *t_max);

        t->ElementWiseMin(*t_max);
        DBG_PRINT_VECTOR(2, "new_slack", *t);

        slack = t;
        return retval;
      }
    }

    return retval;
  }

  Index
  IpoptCalculatedQuantities::AdjustedTrialSlacks()
  {
    DBG_START_METH("IpoptCalculatedQuantities::AdjustedTrialSlacks()",
                   dbg_verbosity);
    Index result =  (num_adjusted_slack_x_L_ +
                     num_adjusted_slack_x_U_ +
                     num_adjusted_slack_s_L_ +
                     num_adjusted_slack_s_U_);
    DBG_PRINT((1,"result = %d\n", result));
    return result;
  }

  void
  IpoptCalculatedQuantities::ResetAdjustedTrialSlacks()
  {
    DBG_START_METH("IpoptCalculatedQuantities::ResetAdjustedTrialSlacks()",
                   dbg_verbosity);
    num_adjusted_slack_x_L_
    = num_adjusted_slack_x_U_
      = num_adjusted_slack_s_L_
        = num_adjusted_slack_s_U_ = 0;
  }

  ///////////////////////////////////////////////////////////////////////////
  //                          Objective Function                           //
  ///////////////////////////////////////////////////////////////////////////

  Number
  IpoptCalculatedQuantities::curr_f()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_f()",
                   dbg_verbosity);
    Number result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();
    DBG_PRINT_VECTOR(2,"curr_x",*x);
    DBG_PRINT((1, "curr_x tag = %d\n", x->GetTag()));

    bool objective_depends_on_mu = ip_nlp_->objective_depends_on_mu();
    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(x);
    std::vector<Number> sdeps(1);
    if (objective_depends_on_mu) {
      sdeps[0] = ip_data_->curr_mu();
    }
    else {
      sdeps[0] = -1.;
    }

    if (!curr_f_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (!trial_f_cache_.GetCachedResult(result, tdeps, sdeps)) {
        DBG_PRINT((2,"evaluate curr f\n"));
        if (objective_depends_on_mu) {
          result = ip_nlp_->f(*x, ip_data_->curr_mu());
        }
        else {
          result = ip_nlp_->f(*x);
        }
      }
      curr_f_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    DBG_PRINT((1,"result (curr_f) = %e\n", result));
    return result;
  }

  Number
  IpoptCalculatedQuantities::unscaled_curr_f()
  {
    return ip_nlp_->NLP_scaling()->unapply_obj_scaling(curr_f());
  }

  Number
  IpoptCalculatedQuantities::trial_f()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_f()",
                   dbg_verbosity);
    Number result;
    SmartPtr<const Vector> x = ip_data_->trial()->x();
    DBG_PRINT_VECTOR(2,"trial_x",*x);
    DBG_PRINT((1, "trial_x tag = %d\n", x->GetTag()));

    bool objective_depends_on_mu = ip_nlp_->objective_depends_on_mu();
    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(x);
    std::vector<Number> sdeps(1);
    if (objective_depends_on_mu) {
      sdeps[0] = ip_data_->curr_mu();
    }
    else {
      sdeps[0] = -1.;
    }

    if (!trial_f_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (!curr_f_cache_.GetCachedResult(result, tdeps, sdeps)) {
        DBG_PRINT((2,"evaluate trial f\n"));
        if (objective_depends_on_mu) {
          result = ip_nlp_->f(*x, ip_data_->curr_mu());
        }
        else {
          result = ip_nlp_->f(*x);
        }
      }
      trial_f_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    DBG_PRINT((1,"result (trial_f) = %e\n", result));
    return result;
  }

  Number
  IpoptCalculatedQuantities::unscaled_trial_f()
  {
    return ip_nlp_->NLP_scaling()->unapply_obj_scaling(trial_f());
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_grad_f()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_grad_f()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();

    bool objective_depends_on_mu = ip_nlp_->objective_depends_on_mu();
    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(x);
    std::vector<Number> sdeps(1);
    if (objective_depends_on_mu) {
      sdeps[0] = ip_data_->curr_mu();
    }
    else {
      sdeps[0] = -1.;
    }

    if (!curr_grad_f_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (!trial_grad_f_cache_.GetCachedResult(result, tdeps, sdeps)) {
        if (objective_depends_on_mu) {
          result = ip_nlp_->grad_f(*x, ip_data_->curr_mu());
        }
        else {
          result = ip_nlp_->grad_f(*x);
        }
      }
      curr_grad_f_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_grad_f()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_grad_f()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->trial()->x();

    bool objective_depends_on_mu = ip_nlp_->objective_depends_on_mu();
    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(x);
    std::vector<Number> sdeps(1);
    if (objective_depends_on_mu) {
      sdeps[0] = ip_data_->curr_mu();
    }
    else {
      sdeps[0] = -1.;
    }

    if (!trial_grad_f_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (!curr_grad_f_cache_.GetCachedResult(result, tdeps, sdeps)) {
        if (objective_depends_on_mu) {
          result = ip_nlp_->grad_f(*x, ip_data_->curr_mu());
        }
        else {
          result = ip_nlp_->grad_f(*x);
        }
      }
      trial_grad_f_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  ///////////////////////////////////////////////////////////////////////////
  //                    Barrier Objective Function                         //
  ///////////////////////////////////////////////////////////////////////////
  Number
  IpoptCalculatedQuantities::CalcBarrierTerm(Number mu,
      const Vector& slack_x_L,
      const Vector& slack_x_U,
      const Vector& slack_s_L,
      const Vector& slack_s_U)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalcBarrierTerm",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);

    DBG_PRINT_VECTOR(2, "slack_x_L", slack_x_L);
    DBG_PRINT_VECTOR(2, "slack_x_U", slack_x_U);
    DBG_PRINT_VECTOR(2, "slack_s_L", slack_s_L);
    DBG_PRINT_VECTOR(2, "slack_s_U", slack_s_U);

    Number retval=0.;
    retval += slack_x_L.SumLogs();
    DBG_PRINT((1, "BarrierTerm after x_L = %25.16e\n", retval));
    retval += slack_x_U.SumLogs();
    DBG_PRINT((1, "BarrierTerm after x_U = %25.16e\n", retval));
    retval += slack_s_L.SumLogs();
    DBG_PRINT((1, "BarrierTerm after s_L = %25.16e\n", retval));
    retval += slack_s_U.SumLogs();
    DBG_PRINT((1, "BarrierTerm after s_U = %25.16e\n", retval));
    retval *= -mu;

    DBG_PRINT((1, "BarrierTerm without damping = %25.16e\n", retval));

    // Include the linear damping term if kappa_d is nonzero.
    if (kappa_d_>0) {
      SmartPtr<const Vector> dampind_x_L;
      SmartPtr<const Vector> dampind_x_U;
      SmartPtr<const Vector> dampind_s_L;
      SmartPtr<const Vector> dampind_s_U;
      ComputeDampingIndicators(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

      Tmp_x_L().Copy(slack_x_L);
      Tmp_x_L().ElementWiseMultiply(*dampind_x_L);
      retval += kappa_d_ * mu * Tmp_x_L().Asum();
      Tmp_x_U().Copy(slack_x_U);
      Tmp_x_U().ElementWiseMultiply(*dampind_x_U);
      retval += kappa_d_ * mu * Tmp_x_U().Asum();
      Tmp_s_L().Copy(slack_s_L);
      Tmp_s_L().ElementWiseMultiply(*dampind_s_L);
      retval += kappa_d_ * mu * Tmp_s_L().Asum();
      Tmp_s_U().Copy(slack_s_U);
      Tmp_s_U().ElementWiseMultiply(*dampind_s_U);
      retval += kappa_d_ * mu * Tmp_s_U().Asum();
    }

    DBG_PRINT((1, "BarrierTerm with damping = %25.16e\n", retval));

    DBG_ASSERT(IsFiniteNumber(retval));
    return retval;
  }

  Number
  IpoptCalculatedQuantities::curr_barrier_obj()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_barrier_obj()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    DBG_PRINT_VECTOR(2,"curr_x",*x);
    DBG_PRINT_VECTOR(2,"curr_s",*s);
    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);

    Number mu = ip_data_->curr_mu();
    DBG_PRINT((1,"curr_mu=%e\n",mu));
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;

    if (!curr_barrier_obj_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (!trial_barrier_obj_cache_.GetCachedResult(result, tdeps, sdeps)) {
        result = curr_f();
        DBG_PRINT((1,"curr_F=%e\n",result));
        result += CalcBarrierTerm(mu,
                                  *curr_slack_x_L(),
                                  *curr_slack_x_U(),
                                  *curr_slack_s_L(),
                                  *curr_slack_s_U());
      }
      curr_barrier_obj_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    DBG_ASSERT(IsFiniteNumber(result));
    return result;
  }

  Number
  IpoptCalculatedQuantities::trial_barrier_obj()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_barrier_obj()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->trial()->x();
    SmartPtr<const Vector> s = ip_data_->trial()->s();
    DBG_PRINT_VECTOR(2,"trial_x",*x);
    DBG_PRINT_VECTOR(2,"trial_s",*s);
    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);

    Number mu = ip_data_->curr_mu();
    DBG_PRINT((1,"trial_mu=%e\n",mu));
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;

    if (!trial_barrier_obj_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (!curr_barrier_obj_cache_.GetCachedResult(result, tdeps, sdeps)) {
        result = trial_f();
        DBG_PRINT((1,"trial_F=%e\n",result));
        DBG_PRINT_VECTOR(2, "trial_slack_s_U", *trial_slack_s_U());
        result += CalcBarrierTerm(ip_data_->curr_mu(),
                                  *trial_slack_x_L(),
                                  *trial_slack_x_U(),
                                  *trial_slack_s_L(),
                                  *trial_slack_s_U());
      }
      trial_barrier_obj_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    DBG_ASSERT(IsFiniteNumber(result));
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_grad_barrier_obj_x()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_grad_barrier_obj_x()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(x);
    Number mu = ip_data_->curr_mu();
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;
    DBG_PRINT((1,"curr_mu=%e\n",mu));

    if (!curr_grad_barrier_obj_x_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp1 = x->MakeNew();
      tmp1->Copy(*curr_grad_f());

      Tmp_x_L().Set(1.);
      ip_nlp_->Px_L()->AddMSinvZ(-mu, *curr_slack_x_L(), Tmp_x_L(), *tmp1);

      Tmp_x_U().Set(1.);
      ip_nlp_->Px_U()->AddMSinvZ(mu, *curr_slack_x_U(), Tmp_x_U(), *tmp1);

      DBG_PRINT_VECTOR(2, "Barrier_Grad_x without damping", *tmp1);

      // Take care of linear damping terms
      if (kappa_d_>0.) {
        SmartPtr<const Vector> dampind_x_L;
        SmartPtr<const Vector> dampind_x_U;
        SmartPtr<const Vector> dampind_s_L;
        SmartPtr<const Vector> dampind_s_U;
        ComputeDampingIndicators(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

        DBG_PRINT((1, "kappa_d*mu = %e\n", kappa_d_*mu));
        DBG_PRINT_VECTOR(2, "dampind_x_L", *dampind_x_L);
        ip_nlp_->Px_L()->MultVector(kappa_d_*mu, *dampind_x_L, 1., *tmp1);
        ip_nlp_->Px_U()->MultVector(-kappa_d_*mu, *dampind_x_U, 1., *tmp1);
      }

      DBG_PRINT_VECTOR(2, "Barrier_Grad_x with damping", *tmp1);

      result = ConstPtr(tmp1);

      curr_grad_barrier_obj_x_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::grad_kappa_times_damping_x()
  {
    DBG_START_METH("IpoptCalculatedQuantities::grad_kappa_times_damping_x()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();

    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(ip_nlp_->Px_L());
    tdeps[1] = GetRawPtr(ip_nlp_->Px_U());
    std::vector<Number> sdeps(1);
    sdeps[0] = kappa_d_;
    if (!grad_kappa_times_damping_x_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp1 = x->MakeNew();
      if (kappa_d_>0.) {
        SmartPtr<const Vector> dampind_x_L;
        SmartPtr<const Vector> dampind_x_U;
        SmartPtr<const Vector> dampind_s_L;
        SmartPtr<const Vector> dampind_s_U;
        ComputeDampingIndicators(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

        ip_nlp_->Px_L()->MultVector(kappa_d_, *dampind_x_L, 0., *tmp1);
        ip_nlp_->Px_U()->MultVector(-kappa_d_, *dampind_x_U, 1., *tmp1);
      }
      else {
        tmp1->Set(0.);
      }
      result = ConstPtr(tmp1);

      grad_kappa_times_damping_x_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_grad_barrier_obj_s()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_grad_barrier_obj_s()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> s = ip_data_->curr()->s();
    std::vector<const TaggedObject*> tdeps(1);
    tdeps[0] = GetRawPtr(s);
    Number mu = ip_data_->curr_mu();
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;
    DBG_PRINT((1,"curr_mu=%e\n",mu));

    if (!curr_grad_barrier_obj_s_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp1 = s->MakeNew();

      Tmp_s_L().Set(-mu);
      Tmp_s_L().ElementWiseDivide(*curr_slack_s_L());
      ip_nlp_->Pd_L()->MultVector(1., Tmp_s_L(), 0., *tmp1);

      Tmp_s_U().Set(1.);
      ip_nlp_->Pd_U()->AddMSinvZ(mu, *curr_slack_s_U(), Tmp_s_U(), *tmp1);

      DBG_PRINT_VECTOR(2, "Barrier_Grad_s without damping", *tmp1);

      // Take care of linear damping terms
      if (kappa_d_>0.) {
        SmartPtr<const Vector> dampind_x_L;
        SmartPtr<const Vector> dampind_x_U;
        SmartPtr<const Vector> dampind_s_L;
        SmartPtr<const Vector> dampind_s_U;
        ComputeDampingIndicators(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

        DBG_PRINT((1, "kappa_d*mu = %e\n", kappa_d_*mu));
        DBG_PRINT_VECTOR(2, "dampind_s_L", *dampind_s_L);
        DBG_PRINT_VECTOR(2, "dampind_s_U", *dampind_s_U);
        ip_nlp_->Pd_L()->MultVector(kappa_d_*mu, *dampind_s_L, 1., *tmp1);
        ip_nlp_->Pd_U()->MultVector(-kappa_d_*mu, *dampind_s_U, 1., *tmp1);
      }

      DBG_PRINT_VECTOR(2, "Barrier_Grad_s with damping", *tmp1);

      result = ConstPtr(tmp1);

      curr_grad_barrier_obj_s_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::grad_kappa_times_damping_s()
  {
    DBG_START_METH("IpoptCalculatedQuantities::grad_kappa_times_damping_s()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> s = ip_data_->curr()->s();

    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(ip_nlp_->Pd_L());
    tdeps[1] = GetRawPtr(ip_nlp_->Pd_U());
    std::vector<Number> sdeps(1);
    sdeps[0] = kappa_d_;
    if (!grad_kappa_times_damping_s_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp1 = s->MakeNew();
      if (kappa_d_>0.) {
        SmartPtr<const Vector> dampind_x_L;
        SmartPtr<const Vector> dampind_x_U;
        SmartPtr<const Vector> dampind_s_L;
        SmartPtr<const Vector> dampind_s_U;
        ComputeDampingIndicators(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

        ip_nlp_->Pd_L()->MultVector(kappa_d_, *dampind_s_L, 0., *tmp1);
        ip_nlp_->Pd_U()->MultVector(-kappa_d_, *dampind_s_U, 1., *tmp1);
      }
      else {
        tmp1->Set(0.);
      }
      result = ConstPtr(tmp1);

      grad_kappa_times_damping_s_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  void
  IpoptCalculatedQuantities::ComputeDampingIndicators(
    SmartPtr<const Vector>& dampind_x_L,
    SmartPtr<const Vector>& dampind_x_U,
    SmartPtr<const Vector>& dampind_s_L,
    SmartPtr<const Vector>& dampind_s_U)
  {
    DBG_START_METH("IpoptCalculatedQuantities::ComputeDampingFilters()",
                   dbg_verbosity);

    // Assume that all indicators have to be computed if one of the
    // SmartPtrs is still zero.
    if (IsNull(dampind_x_L_)) {
      // First for x
      Tmp_x_L().Set(1.0);
      ip_nlp_->Px_L()->MultVector(1.0, Tmp_x_L(), 0.0, Tmp_x());
      Tmp_x_U().Set(1.0);
      ip_nlp_->Px_U()->MultVector(-1.0, Tmp_x_U(), 1.0, Tmp_x());

      dampind_x_L_ = ip_nlp_->x_L()->MakeNew();
      ip_nlp_->Px_L()->TransMultVector(1.0, Tmp_x(), 0.0, *dampind_x_L_);

      dampind_x_U_ = ip_nlp_->x_U()->MakeNew();
      ip_nlp_->Px_U()->TransMultVector(-1.0, Tmp_x(), 0.0, *dampind_x_U_);

      // Now for s
      Tmp_s_L().Set(1.0);
      ip_nlp_->Pd_L()->MultVector(1.0, Tmp_s_L(), 0.0, Tmp_s());
      Tmp_s_U().Set(1.0);
      ip_nlp_->Pd_U()->MultVector(-1.0, Tmp_s_U(), 1.0, Tmp_s());

      dampind_s_L_ = ip_nlp_->d_L()->MakeNew();
      ip_nlp_->Pd_L()->TransMultVector(1.0, Tmp_s(), 0.0, *dampind_s_L_);

      dampind_s_U_ = ip_nlp_->d_U()->MakeNew();
      ip_nlp_->Pd_U()->TransMultVector(-1.0, Tmp_s(), 0.0, *dampind_s_U_);

      DBG_PRINT_VECTOR(2, "dampind_x_L_", *dampind_x_L_);
      DBG_PRINT_VECTOR(2, "dampind_x_U_", *dampind_x_U_);
      DBG_PRINT_VECTOR(2, "dampind_s_L_", *dampind_s_L_);
      DBG_PRINT_VECTOR(2, "dampind_s_U_", *dampind_s_U_);
    }

    dampind_x_L = ConstPtr(dampind_x_L_);
    dampind_x_U = ConstPtr(dampind_x_U_);
    dampind_s_L = ConstPtr(dampind_s_L_);
    dampind_s_U = ConstPtr(dampind_s_U_);
  }

  ///////////////////////////////////////////////////////////////////////////
  //                                Constraints                            //
  ///////////////////////////////////////////////////////////////////////////

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_c()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_c()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();

    if (!curr_c_cache_.GetCachedResult1Dep(result, *x)) {
      if (!trial_c_cache_.GetCachedResult1Dep(result, *x)) {
        result = ip_nlp_->c(*x);
      }
      curr_c_cache_.AddCachedResult1Dep(result, *x);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::unscaled_curr_c()
  {
    return ip_nlp_->NLP_scaling()->unapply_vector_scaling_c(curr_c());
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_c()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_c()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->trial()->x();

    if (!trial_c_cache_.GetCachedResult1Dep(result, *x)) {
      if (!curr_c_cache_.GetCachedResult1Dep(result, *x)) {
        result = ip_nlp_->c(*x);
      }
      trial_c_cache_.AddCachedResult1Dep(result, *x);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_d()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_d()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();

    if (!curr_d_cache_.GetCachedResult1Dep(result, *x)) {
      if (!trial_d_cache_.GetCachedResult1Dep(result, *x)) {
        result = ip_nlp_->d(*x);
      }
      curr_d_cache_.AddCachedResult1Dep(result, *x);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::unscaled_curr_d()
  {
    return ip_nlp_->NLP_scaling()->unapply_vector_scaling_d(curr_d());
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_d()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_d()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->trial()->x();

    if (!trial_d_cache_.GetCachedResult1Dep(result, *x)) {
      if (!curr_d_cache_.GetCachedResult1Dep(result, *x)) {
        result = ip_nlp_->d(*x);
      }
      trial_d_cache_.AddCachedResult1Dep(result, *x);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_d_minus_s()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_d_minus_s()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();

    if (!curr_d_minus_s_cache_.GetCachedResult2Dep(result, *x, *s)) {
      if (!trial_d_minus_s_cache_.GetCachedResult2Dep(result, *x, *s)) {
        SmartPtr<Vector> tmp = s->MakeNew();
        tmp->AddTwoVectors(1., *curr_d(), -1., *s, 0.);
        result = ConstPtr(tmp);
      }
      curr_d_minus_s_cache_.AddCachedResult2Dep(result, *x, *s);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_d_minus_s()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_d_minus_s()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_->trial()->x();
    SmartPtr<const Vector> s = ip_data_->trial()->s();

    if (!trial_d_minus_s_cache_.GetCachedResult2Dep(result, *x, *s)) {
      if (!curr_d_minus_s_cache_.GetCachedResult2Dep(result, *x, *s)) {
        SmartPtr<Vector> tmp = s->MakeNew();
        tmp->AddTwoVectors(1., *trial_d(), -1., *s, 0.);
        result = ConstPtr(tmp);
      }
      trial_d_minus_s_cache_.AddCachedResult2Dep(result, *x, *s);
    }

    return result;
  }

  SmartPtr<const Matrix>
  IpoptCalculatedQuantities::curr_jac_c()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_c()",
                   dbg_verbosity);
    SmartPtr<const Matrix> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();

    if (!curr_jac_c_cache_.GetCachedResult1Dep(result, *x)) {
      if (!trial_jac_c_cache_.GetCachedResult1Dep(result, *x)) {
        result = ip_nlp_->jac_c(*x);
      }
      curr_jac_c_cache_.AddCachedResult1Dep(result, *x);
    }
    return result;
  }

  SmartPtr<const Matrix>
  IpoptCalculatedQuantities::trial_jac_c()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_jac_c()",
                   dbg_verbosity);
    SmartPtr<const Matrix> result;
    SmartPtr<const Vector> x = ip_data_->trial()->x();

    if (!trial_jac_c_cache_.GetCachedResult1Dep(result, *x)) {
      if (!curr_jac_c_cache_.GetCachedResult1Dep(result, *x)) {
        result = ip_nlp_->jac_c(*x);
      }
      trial_jac_c_cache_.AddCachedResult1Dep(result, *x);
    }
    return result;
  }

  SmartPtr<const Matrix>
  IpoptCalculatedQuantities::curr_jac_d()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_d()",
                   dbg_verbosity);
    SmartPtr<const Matrix> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();

    if (!curr_jac_d_cache_.GetCachedResult1Dep(result, *x)) {
      if (!trial_jac_d_cache_.GetCachedResult1Dep(result, *x)) {
        result = ip_nlp_->jac_d(*x);
      }
      curr_jac_d_cache_.AddCachedResult1Dep(result, *x);
    }
    return result;
  }

  SmartPtr<const Matrix>
  IpoptCalculatedQuantities::trial_jac_d()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_jac_d()",
                   dbg_verbosity);
    SmartPtr<const Matrix> result;
    SmartPtr<const Vector> x = ip_data_->trial()->x();

    if (!trial_jac_d_cache_.GetCachedResult1Dep(result, *x)) {
      if (!curr_jac_d_cache_.GetCachedResult1Dep(result, *x)) {
        result = ip_nlp_->jac_d(*x);
      }
      trial_jac_d_cache_.AddCachedResult1Dep(result, *x);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_jac_c_times_vec(const Vector& vec)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_c_times_vec",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();

    if (!curr_jac_c_times_vec_cache_.GetCachedResult2Dep(result, *x, vec)) {
      SmartPtr<Vector> tmp = ip_data_->curr()->y_c()->MakeNew();
      curr_jac_c()->MultVector(1.0, vec, 0., *tmp);
      result = ConstPtr(tmp);
      curr_jac_c_times_vec_cache_.AddCachedResult2Dep(result, *x, vec);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_jac_d_times_vec(const Vector& vec)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_d_times_vec()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();

    if (!curr_jac_d_times_vec_cache_.GetCachedResult2Dep(result, *x, vec)) {
      SmartPtr<Vector> tmp = ip_data_->curr()->s()->MakeNew();
      curr_jac_d()->MultVector(1.0, vec, 0., *tmp);
      result = ConstPtr(tmp);
      curr_jac_d_times_vec_cache_.AddCachedResult2Dep(result, *x, vec);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_jac_cT_times_curr_y_c()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_cT_times_curr_y_c()",
                   dbg_verbosity);
    return curr_jac_cT_times_vec(*ip_data_->curr()->y_c());
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_jac_cT_times_trial_y_c()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_jac_cT_times_trial_y_c()",
                   dbg_verbosity);
    return trial_jac_cT_times_vec(*ip_data_->trial()->y_c());
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_jac_dT_times_curr_y_d()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_dT_times_curr_y_d()",
                   dbg_verbosity);
    return curr_jac_dT_times_vec(*ip_data_->curr()->y_d());
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_jac_dT_times_trial_y_d()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_jac_dT_times_trial_y_d()",
                   dbg_verbosity);
    return trial_jac_dT_times_vec(*ip_data_->trial()->y_d());
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_jac_cT_times_vec(const Vector& vec)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_cT_times_vec",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();

    if (!curr_jac_cT_times_vec_cache_.GetCachedResult2Dep(result, *x, vec)) {
      if (!trial_jac_cT_times_vec_cache_.GetCachedResult2Dep(result, *x, vec)) {
        SmartPtr<Vector> tmp = x->MakeNew();
        curr_jac_c()->TransMultVector(1.0, vec, 0., *tmp);
        result = ConstPtr(tmp);
      }
      curr_jac_cT_times_vec_cache_.AddCachedResult2Dep(result, *x, vec);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_jac_cT_times_vec(const Vector& vec)
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_jac_cT_times_vec",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->trial()->x();

    if (!trial_jac_cT_times_vec_cache_.GetCachedResult2Dep(result, *x, vec)) {
      if (!curr_jac_cT_times_vec_cache_.GetCachedResult2Dep(result, *x, vec)) {
        SmartPtr<Vector> tmp = x->MakeNew();
        trial_jac_c()->TransMultVector(1.0, vec, 0., *tmp);
        result = ConstPtr(tmp);
      }
      trial_jac_cT_times_vec_cache_.AddCachedResult2Dep(result, *x, vec);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_jac_dT_times_vec(const Vector& vec)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_jac_dT_times_vec()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();

    if (!curr_jac_dT_times_vec_cache_.GetCachedResult2Dep(result, *x, vec)) {
      if (!trial_jac_dT_times_vec_cache_.GetCachedResult2Dep(result, *x, vec)) {
        SmartPtr<Vector> tmp = x->MakeNew();
        curr_jac_d()->TransMultVector(1.0, vec, 0., *tmp);
        result = ConstPtr(tmp);
      }
      curr_jac_dT_times_vec_cache_.AddCachedResult2Dep(result, *x, vec);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_jac_dT_times_vec(const Vector& vec)
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_jac_dT_times_vec()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->trial()->x();

    if (!trial_jac_dT_times_vec_cache_.GetCachedResult2Dep(result, *x, vec)) {
      if (!curr_jac_dT_times_vec_cache_.GetCachedResult2Dep(result, *x, vec)) {
        SmartPtr<Vector> tmp = x->MakeNew();
        trial_jac_d()->TransMultVector(1.0, vec, 0., *tmp);
        result = ConstPtr(tmp);
      }
      trial_jac_dT_times_vec_cache_.AddCachedResult2Dep(result, *x, vec);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_constraint_violation()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_constraint_violation()",
                   dbg_verbosity);
    return curr_primal_infeasibility(constr_viol_normtype_);
  }

  Number
  IpoptCalculatedQuantities::trial_constraint_violation()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_constraint_violation()",
                   dbg_verbosity);
    return trial_primal_infeasibility(constr_viol_normtype_);
  }

  Number
  IpoptCalculatedQuantities::curr_nlp_constraint_violation
  (ENormType NormType)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_nlp_constraint_violation()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();

    std::vector<const TaggedObject*> deps(1);
    deps[0] = GetRawPtr(x);
    std::vector<Number> sdeps(1);
    sdeps[0] = (Number)NormType;

    if (!curr_nlp_constraint_violation_cache_.GetCachedResult(result, deps, sdeps)) {
      SmartPtr<const Vector> c = curr_c();
      SmartPtr<const Vector> d = curr_d();

      SmartPtr<Vector> d_viol_L = ip_nlp_->d_L()->MakeNewCopy();
      ip_nlp_->Pd_L()->TransMultVector(-1., *d, 1., *d_viol_L);
      SmartPtr<Vector> tmp = d_viol_L->MakeNew();
      tmp->Set(0.);
      d_viol_L->ElementWiseMax(*tmp);
      DBG_PRINT_VECTOR(2, "d_viol_L", *d_viol_L);

      SmartPtr<Vector> d_viol_U = ip_nlp_->d_U()->MakeNewCopy();
      ip_nlp_->Pd_U()->TransMultVector(-1., *d, 1., *d_viol_U);
      tmp = d_viol_U->MakeNew();
      tmp->Set(0.);
      d_viol_U->ElementWiseMin(*tmp);
      DBG_PRINT_VECTOR(2, "d_viol_U", *d_viol_U);

      std::vector<SmartPtr<const Vector> > vecs(3);
      vecs[0] = c;
      vecs[1] = GetRawPtr(d_viol_L);
      vecs[2] = GetRawPtr(d_viol_U);
      result = CalcNormOfType(NormType, vecs);
      curr_nlp_constraint_violation_cache_.AddCachedResult(result, deps, sdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::unscaled_curr_nlp_constraint_violation
  (ENormType NormType)
  {
    DBG_START_METH("IpoptCalculatedQuantities::unscaled_curr_nlp_constraint_violation()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();

    std::vector<const TaggedObject*> deps(1);
    deps[0] = GetRawPtr(x);
    std::vector<Number> sdeps(1);
    sdeps[0] = (Number)NormType;

    if (!unscaled_curr_nlp_constraint_violation_cache_.GetCachedResult(result, deps, sdeps)) {
      SmartPtr<const Vector> c = unscaled_curr_c();

      SmartPtr<const Vector> d = curr_d();

      SmartPtr<Vector> d_viol_L = ip_nlp_->d_L()->MakeNewCopy();
      ip_nlp_->Pd_L()->TransMultVector(-1., *d, 1., *d_viol_L);
      SmartPtr<Vector> tmp = d_viol_L->MakeNew();
      tmp->Set(0.);
      d_viol_L->ElementWiseMax(*tmp);
      DBG_PRINT_VECTOR(2, "d_viol_L", *d_viol_L);

      SmartPtr<Vector> d_viol_U = ip_nlp_->d_U()->MakeNewCopy();
      ip_nlp_->Pd_U()->TransMultVector(-1., *d, 1., *d_viol_U);
      tmp = d_viol_U->MakeNew();
      tmp->Set(0.);
      d_viol_U->ElementWiseMin(*tmp);
      DBG_PRINT_VECTOR(2, "d_viol_U", *d_viol_U);

      std::vector<SmartPtr<const Vector> > vecs(3);
      vecs[0] = c;
      vecs[1] = GetRawPtr(d_viol_L);
      vecs[2] = GetRawPtr(d_viol_U);
      result = CalcNormOfType(NormType, vecs);
      unscaled_curr_nlp_constraint_violation_cache_.AddCachedResult(result, deps, sdeps);
    }

    return result;
  }

  ///////////////////////////////////////////////////////////////////////////
  //                Exact Hessian using second derivatives                 //
  ///////////////////////////////////////////////////////////////////////////

  SmartPtr<const SymMatrix>
  IpoptCalculatedQuantities::curr_exact_hessian()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_exact_hessian()",
                   dbg_verbosity);

    SmartPtr<const SymMatrix> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> y_c = ip_data_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();

    bool objective_depends_on_mu = ip_nlp_->objective_depends_on_mu();
    std::vector<const TaggedObject*> tdeps(3);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(y_c);
    tdeps[2] = GetRawPtr(y_d);
    std::vector<Number> sdeps(1);
    if (objective_depends_on_mu) {
      sdeps[0] = ip_data_->curr_mu();
    }
    else {
      sdeps[0] = -1.;
    }

    if (!curr_exact_hessian_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (objective_depends_on_mu) {
        result = ip_nlp_->h(*x, 1.0, *y_c, *y_d, ip_data_->curr_mu());
      }
      else {
        result = ip_nlp_->h(*x, 1.0, *y_c, *y_d);
      }
      curr_exact_hessian_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  ///////////////////////////////////////////////////////////////////////////
  //                  Optimality Error and its components                  //
  ///////////////////////////////////////////////////////////////////////////

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_grad_lag_x()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_grad_lag_x()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> y_c = ip_data_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    SmartPtr<const Vector> z_L = ip_data_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr()->z_U();

    std::vector<const TaggedObject*> deps(5);
    deps[0] = GetRawPtr(x);
    deps[1] = GetRawPtr(y_c);
    deps[2] = GetRawPtr(y_d);
    deps[3] = GetRawPtr(z_L);
    deps[4] = GetRawPtr(z_U);

    if (!curr_grad_lag_x_cache_.GetCachedResult(result, deps)) {
      if (!trial_grad_lag_x_cache_.GetCachedResult(result, deps)) {
        SmartPtr<Vector> tmp = x->MakeNew();
        DBG_PRINT_VECTOR(2,"curr_grad_f",*curr_grad_f());
        tmp->Copy(*curr_grad_f());
        tmp->AddTwoVectors(1., *curr_jac_cT_times_curr_y_c(),
                           1., *curr_jac_dT_times_curr_y_d(), 1.);
        DBG_PRINT_VECTOR(2,"jac_cT*y_c",*curr_jac_cT_times_curr_y_c());
        DBG_PRINT_VECTOR(2,"jac_dT*y_d",*curr_jac_dT_times_curr_y_d());
        ip_nlp_->Px_L()->MultVector(-1., *z_L, 1., *tmp);
        ip_nlp_->Px_U()->MultVector(1., *z_U, 1., *tmp);
        result = ConstPtr(tmp);
      }
      curr_grad_lag_x_cache_.AddCachedResult(result, deps);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_grad_lag_x()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_grad_lag_x()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_->trial()->x();
    SmartPtr<const Vector> y_c = ip_data_->trial()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->trial()->y_d();
    SmartPtr<const Vector> z_L = ip_data_->trial()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->trial()->z_U();

    std::vector<const TaggedObject*> deps(5);
    deps[0] = GetRawPtr(x);
    deps[1] = GetRawPtr(y_c);
    deps[2] = GetRawPtr(y_d);
    deps[3] = GetRawPtr(z_L);
    deps[4] = GetRawPtr(z_U);

    if (!trial_grad_lag_x_cache_.GetCachedResult(result, deps)) {
      if (!curr_grad_lag_x_cache_.GetCachedResult(result, deps)) {
        SmartPtr<Vector> tmp = x->MakeNew();
        DBG_PRINT_VECTOR(2,"trial_grad_f",*trial_grad_f());
        tmp->Copy(*trial_grad_f());
        tmp->AddTwoVectors(1., *trial_jac_cT_times_trial_y_c(),
                           1., *trial_jac_dT_times_trial_y_d(), 1.);
        ip_nlp_->Px_L()->MultVector(-1., *z_L, 1., *tmp);
        ip_nlp_->Px_U()->MultVector(1., *z_U, 1., *tmp);
        result = ConstPtr(tmp);
      }
      trial_grad_lag_x_cache_.AddCachedResult(result, deps);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_grad_lag_s()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_grad_lag_s()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    SmartPtr<const Vector> v_L = ip_data_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr()->v_U();

    std::vector<const TaggedObject*> deps(3);
    deps[0] = GetRawPtr(y_d);
    deps[1] = GetRawPtr(v_L);
    deps[2] = GetRawPtr(v_U);

    if (!curr_grad_lag_s_cache_.GetCachedResult(result, deps)) {
      if (!trial_grad_lag_s_cache_.GetCachedResult(result, deps)) {
        SmartPtr<Vector> tmp = y_d->MakeNew();
        ip_nlp_->Pd_U()->MultVector(1., *v_U, 0., *tmp);
        ip_nlp_->Pd_L()->MultVector(-1., *v_L, 1., *tmp);
        tmp->Axpy(-1., *y_d);
        result = ConstPtr(tmp);
      }
      curr_grad_lag_s_cache_.AddCachedResult(result, deps);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_grad_lag_s()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_grad_lag_s()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> y_d = ip_data_->trial()->y_d();
    SmartPtr<const Vector> v_L = ip_data_->trial()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->trial()->v_U();

    std::vector<const TaggedObject*> deps(3);
    deps[0] = GetRawPtr(y_d);
    deps[1] = GetRawPtr(v_L);
    deps[2] = GetRawPtr(v_U);

    if (!trial_grad_lag_s_cache_.GetCachedResult(result, deps)) {
      if (!curr_grad_lag_s_cache_.GetCachedResult(result, deps)) {
        SmartPtr<Vector> tmp = y_d->MakeNew();
        ip_nlp_->Pd_U()->MultVector(1., *v_U, 0., *tmp);
        ip_nlp_->Pd_L()->MultVector(-1., *v_L, 1., *tmp);
        tmp->Axpy(-1., *y_d);
        result = ConstPtr(tmp);
      }
      trial_grad_lag_s_cache_.AddCachedResult(result, deps);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_grad_lag_with_damping_x()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_grad_lag_with_damping_x()",
                   dbg_verbosity);

    /* If no damping is used, just return the gradient of the regular
       Lagrangian function */
    if (kappa_d_==0.) {
      return curr_grad_lag_x();
    }

    SmartPtr<const Vector> result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> y_c = ip_data_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    SmartPtr<const Vector> z_L = ip_data_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr()->z_U();
    Number mu = ip_data_->curr_mu();

    std::vector<const TaggedObject*> deps(5);
    deps[0] = GetRawPtr(x);
    deps[1] = GetRawPtr(y_c);
    deps[2] = GetRawPtr(y_d);
    deps[3] = GetRawPtr(z_L);
    deps[4] = GetRawPtr(z_U);
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;

    if (!curr_grad_lag_with_damping_x_cache_.GetCachedResult(result, deps, sdeps)) {
      SmartPtr<Vector> tmp = x->MakeNew();
      tmp->Copy(*curr_grad_lag_x());

      SmartPtr<const Vector> dampind_x_L;
      SmartPtr<const Vector> dampind_x_U;
      SmartPtr<const Vector> dampind_s_L;
      SmartPtr<const Vector> dampind_s_U;
      ComputeDampingIndicators(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

      ip_nlp_->Px_L()->MultVector(kappa_d_*mu, *dampind_x_L, 1., *tmp);
      ip_nlp_->Px_U()->MultVector(-kappa_d_*mu, *dampind_x_U, 1., *tmp);

      result = ConstPtr(tmp);
      curr_grad_lag_with_damping_x_cache_.AddCachedResult(result, deps, sdeps);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_grad_lag_with_damping_s()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_grad_lag_with_damping_s()",
                   dbg_verbosity);

    /* If no damping is used, just return the gradient of the regular
       Lagrangian function */
    if (kappa_d_==0.) {
      return curr_grad_lag_s();
    }

    SmartPtr<const Vector> result;

    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    SmartPtr<const Vector> v_L = ip_data_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr()->v_U();
    Number mu = ip_data_->curr_mu();

    std::vector<const TaggedObject*> deps(3);
    deps[0] = GetRawPtr(y_d);
    deps[1] = GetRawPtr(v_L);
    deps[2] = GetRawPtr(v_U);
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;

    if (!curr_grad_lag_with_damping_s_cache_.GetCachedResult(result, deps, sdeps)) {
      SmartPtr<Vector> tmp = y_d->MakeNew();
      tmp->Copy(*curr_grad_lag_s());

      SmartPtr<const Vector> dampind_x_L;
      SmartPtr<const Vector> dampind_x_U;
      SmartPtr<const Vector> dampind_s_L;
      SmartPtr<const Vector> dampind_s_U;
      ComputeDampingIndicators(dampind_x_L, dampind_x_U, dampind_s_L, dampind_s_U);

      ip_nlp_->Pd_L()->MultVector(kappa_d_*mu, *dampind_s_L, 1., *tmp);
      ip_nlp_->Pd_U()->MultVector(-kappa_d_*mu, *dampind_s_U, 1., *tmp);

      result = ConstPtr(tmp);
      curr_grad_lag_with_damping_s_cache_.AddCachedResult(result, deps, sdeps);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::CalcCompl(const Vector& slack,
                                       const Vector& mult)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalcCompl()",
                   dbg_verbosity);
    SmartPtr<Vector> result = slack.MakeNew();
    result->Copy(slack);
    result->ElementWiseMultiply(mult);
    return ConstPtr(result);
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_compl_x_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_compl_x_L()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_x_L();
    SmartPtr<const Vector> mult = ip_data_->curr()->z_L();
    DBG_PRINT_VECTOR(2, "slack_x_L", *slack);
    DBG_PRINT_VECTOR(2, "z_L", *mult);

    if (!curr_compl_x_L_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
      if (!trial_compl_x_L_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
        result = CalcCompl(*slack, *mult);
      }
      curr_compl_x_L_cache_.AddCachedResult2Dep(result, *slack, *mult);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_compl_x_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_compl_x_L()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = trial_slack_x_L();
    SmartPtr<const Vector> mult = ip_data_->trial()->z_L();
    DBG_PRINT_VECTOR(2, "slack_x_L", *slack);
    DBG_PRINT_VECTOR(2, "z_L", *mult);

    if (!trial_compl_x_L_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
      if (!curr_compl_x_L_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
        result = CalcCompl(*slack, *mult);
      }
      trial_compl_x_L_cache_.AddCachedResult2Dep(result, *slack, *mult);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_compl_x_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_compl_x_U()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_x_U();
    SmartPtr<const Vector> mult = ip_data_->curr()->z_U();

    if (!curr_compl_x_U_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
      if (!trial_compl_x_U_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
        result = CalcCompl(*slack, *mult);
      }
      curr_compl_x_U_cache_.AddCachedResult2Dep(result, *slack, *mult);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_compl_x_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_compl_x_U()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = trial_slack_x_U();
    SmartPtr<const Vector> mult = ip_data_->trial()->z_U();

    if (!trial_compl_x_U_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
      if (!curr_compl_x_U_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
        result = CalcCompl(*slack, *mult);
      }
      trial_compl_x_U_cache_.AddCachedResult2Dep(result, *slack, *mult);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_compl_s_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_compl_s_L()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_s_L();
    SmartPtr<const Vector> mult = ip_data_->curr()->v_L();

    if (!curr_compl_s_L_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
      if (!trial_compl_s_L_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
        result = CalcCompl(*slack, *mult);
      }
      curr_compl_s_L_cache_.GetCachedResult2Dep(result, *slack, *mult);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_compl_s_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_compl_s_L()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = trial_slack_s_L();
    SmartPtr<const Vector> mult = ip_data_->trial()->v_L();

    if (!trial_compl_s_L_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
      if (!curr_compl_s_L_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
        result = CalcCompl(*slack, *mult);
      }
      trial_compl_s_L_cache_.GetCachedResult2Dep(result, *slack, *mult);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_compl_s_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_compl_s_U()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_s_U();
    SmartPtr<const Vector> mult = ip_data_->curr()->v_U();

    if (!curr_compl_s_U_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
      if (!trial_compl_s_U_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
        result = CalcCompl(*slack, *mult);
      }
      curr_compl_s_U_cache_.GetCachedResult2Dep(result, *slack, *mult);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::trial_compl_s_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_compl_s_U()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = trial_slack_s_U();
    SmartPtr<const Vector> mult = ip_data_->trial()->v_U();

    if (!trial_compl_s_U_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
      if (!curr_compl_s_U_cache_.GetCachedResult2Dep(result, *slack, *mult)) {
        result = CalcCompl(*slack, *mult);
      }
      trial_compl_s_U_cache_.GetCachedResult2Dep(result, *slack, *mult);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_relaxed_compl_x_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_relaxed_compl_x_L()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_x_L();
    SmartPtr<const Vector> mult = ip_data_->curr()->z_L();
    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(slack);
    tdeps[1] = GetRawPtr(mult);

    Number mu = ip_data_->curr_mu();
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;

    if (!curr_relaxed_compl_x_L_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp = slack->MakeNew();
      tmp->Copy(*curr_compl_x_L());
      tmp->AddScalar(-mu);
      result = ConstPtr(tmp);
      curr_relaxed_compl_x_L_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_relaxed_compl_x_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_relaxed_compl_x_U()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_x_U();
    SmartPtr<const Vector> mult = ip_data_->curr()->z_U();
    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(slack);
    tdeps[1] = GetRawPtr(mult);

    Number mu = ip_data_->curr_mu();
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;

    if (!curr_relaxed_compl_x_U_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp = slack->MakeNew();
      tmp->Copy(*curr_compl_x_U());
      tmp->AddScalar(-mu);
      result = ConstPtr(tmp);
      curr_relaxed_compl_x_U_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_relaxed_compl_s_L()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_relaxed_compl_s_L()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_s_L();
    SmartPtr<const Vector> mult = ip_data_->curr()->v_L();
    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(slack);
    tdeps[1] = GetRawPtr(mult);

    Number mu = ip_data_->curr_mu();
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;

    if (!curr_relaxed_compl_s_L_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp = slack->MakeNew();
      tmp->Copy(*curr_compl_s_L());
      tmp->AddScalar(-mu);
      result = ConstPtr(tmp);
      curr_relaxed_compl_s_L_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_relaxed_compl_s_U()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_relaxed_compl_s_U()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;

    SmartPtr<const Vector> slack = curr_slack_s_U();
    SmartPtr<const Vector> mult = ip_data_->curr()->v_U();
    std::vector<const TaggedObject*> tdeps(2);
    tdeps[0] = GetRawPtr(slack);
    tdeps[1] = GetRawPtr(mult);

    Number mu = ip_data_->curr_mu();
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;

    if (!curr_relaxed_compl_s_U_cache_.GetCachedResult(result, tdeps, sdeps)) {
      SmartPtr<Vector> tmp = slack->MakeNew();
      tmp->Copy(*curr_compl_s_U());
      tmp->AddScalar(-mu);
      result = ConstPtr(tmp);
      curr_relaxed_compl_s_U_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  Number
  IpoptCalculatedQuantities::CalcNormOfType
  (ENormType NormType,
   const Vector& vec1, const Vector& vec2)
  {
    switch (NormType) {
      case NORM_1 :
      return vec1.Asum() + vec2.Asum();
      case NORM_2 :
      return sqrt(pow(vec1.Nrm2(),2) + pow(vec2.Nrm2(),2));
      case NORM_MAX :
      return Max(vec1.Amax(), vec2.Amax());
      default:
      DBG_ASSERT(false && "Unknown NormType.");
      return 0.0;
    }
  }

  Number
  IpoptCalculatedQuantities::CalcNormOfType
  (ENormType NormType,
   std::vector<SmartPtr<const Vector> > vecs)
  {
    Number result=0.;

    switch (NormType) {
      case NORM_1 :
      for (Index i=0; i<(Index)vecs.size(); i++) {
        result += vecs[i]->Asum();
      }
      break;
      case NORM_2 :
      for (Index i=0; i<(Index)vecs.size(); i++) {
        Number nrm = vecs[i]->Nrm2();
        result += nrm*nrm;
      }
      result = sqrt(result);
      break;
      case NORM_MAX :
      for (Index i=0; i<(Index)vecs.size(); i++) {
        result = Max(result, vecs[i]->Amax());
      }
      break;
      default:
      DBG_ASSERT(false && "Unknown NormType.");
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_primal_infeasibility
  (ENormType NormType)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_primal_infeasibility()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();

    DBG_PRINT_VECTOR(2, "x to eval", *x);
    DBG_PRINT_VECTOR(2, "s to eval", *s);
    DBG_PRINT((1,"NormType = %d\n", NormType))

    std::vector<const TaggedObject*> deps(2);
    deps[0] = GetRawPtr(x);
    deps[1] = GetRawPtr(s);
    std::vector<Number> sdeps(1);
    sdeps[0] = (Number)NormType;

    if (!curr_primal_infeasibility_cache_.GetCachedResult(result, deps, sdeps)) {
      if (!trial_primal_infeasibility_cache_.GetCachedResult(result, deps, sdeps)) {
        DBG_PRINT((1,"Recomputing recomputing infeasibility.\n"));
        SmartPtr<const Vector> c = curr_c();
        SmartPtr<const Vector> d_minus_s = curr_d_minus_s();

        DBG_PRINT_VECTOR(2,"c", *c);
        DBG_PRINT_VECTOR(2,"d_minus_s", *d_minus_s);

        result = CalcNormOfType(NormType, *c, *d_minus_s);

      }
      curr_primal_infeasibility_cache_.AddCachedResult(result, deps, sdeps);
    }

    DBG_PRINT((1,"result = %e\n",result));
    return result;
  }

  Number
  IpoptCalculatedQuantities::trial_primal_infeasibility
  (ENormType NormType)
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_primal_infeasibility()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->trial()->x();
    SmartPtr<const Vector> s = ip_data_->trial()->s();

    DBG_PRINT_VECTOR(2, "x to eval", *x);
    DBG_PRINT_VECTOR(2, "s to eval", *s);
    DBG_PRINT((1,"NormType = %d\n", NormType))

    std::vector<const TaggedObject*> deps(2);
    deps[0] = GetRawPtr(x);
    deps[1] = GetRawPtr(s);
    std::vector<Number> sdeps(1);
    sdeps[0] = (Number)NormType;

    if (!trial_primal_infeasibility_cache_.GetCachedResult(result, deps, sdeps)) {
      if (!curr_primal_infeasibility_cache_.GetCachedResult(result, deps, sdeps)) {
        DBG_PRINT((1,"Recomputing recomputing infeasibility.\n"));
        SmartPtr<const Vector> c = trial_c();
        SmartPtr<const Vector> d_minus_s = trial_d_minus_s();

        DBG_PRINT_VECTOR(2,"c", *c);
        DBG_PRINT_VECTOR(2,"d_minus_s", *d_minus_s);

        result = CalcNormOfType(NormType, *c, *d_minus_s);
      }
      trial_primal_infeasibility_cache_.AddCachedResult(result, deps, sdeps);
    }

    DBG_PRINT((1,"result = %e\n",result));
    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_dual_infeasibility
  (ENormType NormType)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_dual_infeasibility()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    SmartPtr<const Vector> y_c = ip_data_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    SmartPtr<const Vector> z_L = ip_data_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr()->z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr()->v_U();

    std::vector<const TaggedObject*> deps(8);
    deps[0] = GetRawPtr(x);
    deps[1] = GetRawPtr(s);
    deps[2] = GetRawPtr(y_c);
    deps[3] = GetRawPtr(y_d);
    deps[4] = GetRawPtr(z_L);
    deps[5] = GetRawPtr(z_U);
    deps[6] = GetRawPtr(v_L);
    deps[7] = GetRawPtr(v_U);
    std::vector<Number> sdeps(1);
    sdeps[0] = (Number)NormType;

    if (!curr_dual_infeasibility_cache_.GetCachedResult(result, deps, sdeps)) {
      if (!trial_dual_infeasibility_cache_.GetCachedResult(result, deps, sdeps)) {
        SmartPtr<const Vector> grad_lag_x = curr_grad_lag_x();
        SmartPtr<const Vector> grad_lag_s = curr_grad_lag_s();
        DBG_PRINT_VECTOR(2,"grad_lag_x", *grad_lag_x);
        DBG_PRINT_VECTOR(2,"grad_lag_s", *grad_lag_s);

        result = CalcNormOfType(NormType, *grad_lag_x, *grad_lag_s);
      }
      curr_dual_infeasibility_cache_.AddCachedResult(result, deps, sdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::trial_dual_infeasibility
  (ENormType NormType)
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_dual_infeasibility()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->trial()->x();
    SmartPtr<const Vector> s = ip_data_->trial()->s();
    SmartPtr<const Vector> y_c = ip_data_->trial()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->trial()->y_d();
    SmartPtr<const Vector> z_L = ip_data_->trial()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->trial()->z_U();
    SmartPtr<const Vector> v_L = ip_data_->trial()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->trial()->v_U();

    std::vector<const TaggedObject*> deps(8);
    deps[0] = GetRawPtr(x);
    deps[1] = GetRawPtr(s);
    deps[2] = GetRawPtr(y_c);
    deps[3] = GetRawPtr(y_d);
    deps[4] = GetRawPtr(z_L);
    deps[5] = GetRawPtr(z_U);
    deps[6] = GetRawPtr(v_L);
    deps[7] = GetRawPtr(v_U);
    std::vector<Number> sdeps(1);
    sdeps[0] = (Number)NormType;

    if (!trial_dual_infeasibility_cache_.GetCachedResult(result, deps, sdeps)) {
      if (!curr_dual_infeasibility_cache_.GetCachedResult(result, deps, sdeps)) {
        SmartPtr<const Vector> grad_lag_x = trial_grad_lag_x();
        SmartPtr<const Vector> grad_lag_s = trial_grad_lag_s();
        DBG_PRINT_VECTOR(2,"grad_lag_x", *grad_lag_x);
        DBG_PRINT_VECTOR(2,"grad_lag_s", *grad_lag_s);

        result = CalcNormOfType(NormType, *grad_lag_x, *grad_lag_s);

      }
      trial_dual_infeasibility_cache_.AddCachedResult(result, deps, sdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::unscaled_curr_dual_infeasibility
  (ENormType NormType)
  {
    DBG_START_METH("IpoptCalculatedQuantities::unscaled_curr_dual_infeasibility()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    SmartPtr<const Vector> y_c = ip_data_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    SmartPtr<const Vector> z_L = ip_data_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr()->z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr()->v_U();

    std::vector<const TaggedObject*> deps(8);
    deps[0] = GetRawPtr(x);
    deps[1] = GetRawPtr(s);
    deps[2] = GetRawPtr(y_c);
    deps[3] = GetRawPtr(y_d);
    deps[4] = GetRawPtr(z_L);
    deps[5] = GetRawPtr(z_U);
    deps[6] = GetRawPtr(v_L);
    deps[7] = GetRawPtr(v_U);
    std::vector<Number> sdeps(1);
    sdeps[0] = (Number)NormType;

    if (!unscaled_curr_dual_infeasibility_cache_.GetCachedResult(result, deps, sdeps)) {
      SmartPtr<const Vector> grad_lag_x =
        ip_nlp_->NLP_scaling()->unapply_grad_obj_scaling(curr_grad_lag_x());

      Number obj_unscal = ip_nlp_->NLP_scaling()->unapply_obj_scaling(1.);
      SmartPtr<const Vector> grad_lag_s;
      if (obj_unscal != 1.) {
        SmartPtr<Vector> tmp =
          ip_nlp_->NLP_scaling()->apply_vector_scaling_d_NonConst(ConstPtr(curr_grad_lag_s()));
        tmp->Scal(obj_unscal);
        grad_lag_s = ConstPtr(tmp);
      }
      else {
        grad_lag_s = ip_nlp_->NLP_scaling()->apply_vector_scaling_d(curr_grad_lag_s());
      }

      result = CalcNormOfType(NormType, *grad_lag_x, *grad_lag_s);
      unscaled_curr_dual_infeasibility_cache_.AddCachedResult(result, deps, sdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_complementarity
  (Number mu, ENormType NormType)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_complementarity()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    SmartPtr<const Vector> z_L = ip_data_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr()->z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr()->v_U();

    std::vector<const TaggedObject*> deps(6);
    deps[0] = GetRawPtr(x);
    deps[1] = GetRawPtr(s);
    deps[2] = GetRawPtr(z_L);
    deps[3] = GetRawPtr(z_U);
    deps[4] = GetRawPtr(v_L);
    deps[5] = GetRawPtr(v_U);
    std::vector<Number> sdeps(2);
    sdeps[0] = (Number)NormType;
    sdeps[1] = mu;

    if (!curr_complementarity_cache_.GetCachedResult(result, deps, sdeps)) {
      if (!trial_complementarity_cache_.GetCachedResult(result, deps, sdeps)) {

        std::vector<SmartPtr<const Vector> > vecs(4);
        SmartPtr<const Vector> compl_x_L = curr_compl_x_L();
        SmartPtr<const Vector> compl_x_U = curr_compl_x_U();
        SmartPtr<const Vector> compl_s_L = curr_compl_s_L();
        SmartPtr<const Vector> compl_s_U = curr_compl_s_U();

        if (mu==.0) {
          vecs[0] = GetRawPtr(compl_x_L);
          vecs[1] = GetRawPtr(compl_x_U);
          vecs[2] = GetRawPtr(compl_s_L);
          vecs[3] = GetRawPtr(compl_s_U);
        }
        else {
          SmartPtr<Vector> tmp = compl_x_L->MakeNew();
          tmp->Copy(*compl_x_L);
          tmp->AddScalar(-mu);
          vecs[0] = GetRawPtr(tmp);
          tmp = compl_x_U->MakeNew();
          tmp->Copy(*compl_x_U);
          tmp->AddScalar(-mu);
          vecs[1] = GetRawPtr(tmp);
          tmp = compl_s_L->MakeNew();
          tmp->Copy(*compl_s_L);
          tmp->AddScalar(-mu);
          vecs[2] = GetRawPtr(tmp);
          tmp = compl_s_U->MakeNew();
          tmp->Copy(*compl_s_U);
          tmp->AddScalar(-mu);
          vecs[3] = GetRawPtr(tmp);
        }

        result = CalcNormOfType(NormType, vecs);
      }

      curr_complementarity_cache_.AddCachedResult(result, deps, sdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::trial_complementarity
  (Number mu, ENormType NormType)
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_complementarity()",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x = ip_data_->trial()->x();
    SmartPtr<const Vector> s = ip_data_->trial()->s();
    SmartPtr<const Vector> z_L = ip_data_->trial()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->trial()->z_U();
    SmartPtr<const Vector> v_L = ip_data_->trial()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->trial()->v_U();

    std::vector<const TaggedObject*> deps(6);
    deps[0] = GetRawPtr(x);
    deps[1] = GetRawPtr(s);
    deps[2] = GetRawPtr(z_L);
    deps[3] = GetRawPtr(z_U);
    deps[4] = GetRawPtr(v_L);
    deps[5] = GetRawPtr(v_U);
    std::vector<Number> sdeps(2);
    sdeps[0] = (Number)NormType;
    sdeps[1] = mu;

    if (!trial_complementarity_cache_.GetCachedResult(result, deps, sdeps)) {
      if (!curr_complementarity_cache_.GetCachedResult(result, deps, sdeps)) {

        std::vector<SmartPtr<const Vector> > vecs(4);
        SmartPtr<const Vector> compl_x_L = trial_compl_x_L();
        SmartPtr<const Vector> compl_x_U = trial_compl_x_U();
        SmartPtr<const Vector> compl_s_L = trial_compl_s_L();
        SmartPtr<const Vector> compl_s_U = trial_compl_s_U();

        if (mu==.0) {
          vecs[0] = GetRawPtr(compl_x_L);
          vecs[1] = GetRawPtr(compl_x_U);
          vecs[2] = GetRawPtr(compl_s_L);
          vecs[3] = GetRawPtr(compl_s_U);
        }
        else {
          SmartPtr<Vector> tmp = compl_x_L->MakeNew();
          tmp->Copy(*compl_x_L);
          tmp->AddScalar(-mu);
          vecs[0] = GetRawPtr(tmp);
          tmp = compl_x_U->MakeNew();
          tmp->Copy(*compl_x_U);
          tmp->AddScalar(-mu);
          vecs[1] = GetRawPtr(tmp);
          tmp = compl_s_L->MakeNew();
          tmp->Copy(*compl_s_L);
          tmp->AddScalar(-mu);
          vecs[2] = GetRawPtr(tmp);
          tmp = compl_s_U->MakeNew();
          tmp->Copy(*compl_s_U);
          tmp->AddScalar(-mu);
          vecs[3] = GetRawPtr(tmp);
        }

        result = CalcNormOfType(NormType, vecs);
      }

      trial_complementarity_cache_.AddCachedResult(result, deps, sdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::unscaled_curr_complementarity
  (Number mu, ENormType NormType)
  {
    return ip_nlp_->NLP_scaling()->unapply_obj_scaling(curr_complementarity(mu, NormType));
  }

  Number
  IpoptCalculatedQuantities::CalcCentralityMeasure(const Vector& compl_x_L,
      const Vector& compl_x_U,
      const Vector& compl_s_L,
      const Vector& compl_s_U)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalcCentralityMeasure()",
                   dbg_verbosity);

    Number MinCompl = std::numeric_limits<Number>::max();
    bool have_bounds = false;

    Index n_compl_x_L = compl_x_L.Dim();
    Index n_compl_x_U = compl_x_U.Dim();
    Index n_compl_s_L = compl_s_L.Dim();
    Index n_compl_s_U = compl_s_U.Dim();

    // Compute the Minimum of all complementarities
    if( n_compl_x_L>0 ) {
      if( have_bounds ) {
        MinCompl = Min(MinCompl, compl_x_L.Min());
      }
      else {
        MinCompl = compl_x_L.Min();
      }
      have_bounds = true;
    }
    if( n_compl_x_U>0 ) {
      if( have_bounds ) {
        MinCompl = Min(MinCompl, compl_x_U.Min());
      }
      else {
        MinCompl = compl_x_U.Min();
      }
      have_bounds = true;
    }
    if( n_compl_s_L>0 ) {
      if( have_bounds ) {
        MinCompl = Min(MinCompl, compl_s_L.Min());
      }
      else {
        MinCompl = compl_s_L.Min();
      }
      have_bounds = true;
    }
    if( n_compl_s_U>0 ) {
      if( have_bounds ) {
        MinCompl = Min(MinCompl, compl_s_U.Min());
      }
      else {
        MinCompl = compl_s_U.Min();
      }
      have_bounds = true;
    }

    // If there are no bounds, just return 0.;
    if (!have_bounds) {
      return 0.;
    }

    DBG_PRINT_VECTOR(2, "compl_x_L", compl_x_L);
    DBG_PRINT_VECTOR(2, "compl_x_U", compl_x_U);
    DBG_PRINT_VECTOR(2, "compl_s_L", compl_s_L);
    DBG_PRINT_VECTOR(2, "compl_s_U", compl_s_U);

    DBG_ASSERT(MinCompl>0. && "There is a zero complementarity entry");

    Number avrg_compl = (compl_x_L.Asum() + compl_x_U.Asum() +
                         compl_s_L.Asum() + compl_s_U.Asum());
    DBG_PRINT((1,"sum_compl = %25.16e\n", avrg_compl));
    avrg_compl /= (n_compl_x_L + n_compl_x_U + n_compl_s_L + n_compl_s_U);
    DBG_PRINT((1,"avrg_compl = %25.16e\n", avrg_compl));
    DBG_PRINT((1,"MinCompl = %25.16e\n", MinCompl));

    Number xi = MinCompl/avrg_compl;
    // The folloking line added for the case that avrg_compl is
    // slightly smaller than MinCompl, due to numerical roundoff
    xi = Min(1., xi);

    return xi;
  }

  Number
  IpoptCalculatedQuantities::curr_centrality_measure()
  {
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    SmartPtr<const Vector> z_L = ip_data_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr()->z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr()->v_U();

    std::vector<const TaggedObject*> tdeps(6);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    tdeps[2] = GetRawPtr(z_L);
    tdeps[3] = GetRawPtr(z_U);
    tdeps[4] = GetRawPtr(v_L);
    tdeps[5] = GetRawPtr(v_U);

    if (!curr_centrality_measure_cache_.GetCachedResult(result, tdeps)) {
      SmartPtr<const Vector> compl_x_L = curr_compl_x_L();
      SmartPtr<const Vector> compl_x_U = curr_compl_x_U();
      SmartPtr<const Vector> compl_s_L = curr_compl_s_L();
      SmartPtr<const Vector> compl_s_U = curr_compl_s_U();

      result = CalcCentralityMeasure(*compl_x_L, *compl_x_U,
                                     *compl_s_L, *compl_s_U);

      curr_centrality_measure_cache_.AddCachedResult(result, tdeps);
    }
    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_nlp_error()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_nlp_error()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    SmartPtr<const Vector> y_c = ip_data_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    SmartPtr<const Vector> z_L = ip_data_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr()->z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr()->v_U();

    std::vector<const TaggedObject*> tdeps(8);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    tdeps[2] = GetRawPtr(y_c);
    tdeps[3] = GetRawPtr(y_d);
    tdeps[4] = GetRawPtr(z_L);
    tdeps[5] = GetRawPtr(z_U);
    tdeps[6] = GetRawPtr(v_L);
    tdeps[7] = GetRawPtr(v_U);

    if (!curr_nlp_error_cache_.GetCachedResult(result, tdeps)) {
      if (ip_data_->curr()->x()->Dim()==ip_data_->curr()->y_c()->Dim()) {
        // This is a square problem, we only need to consider the
        // infeasibility
        result = curr_nlp_constraint_violation(NORM_MAX);
      }
      else {
        Number s_d = 0;
        Number s_c = 0;
        ComputeOptimalityErrorScaling(*ip_data_->curr()->y_c(), *ip_data_->curr()->y_d(),
                                      *ip_data_->curr()->z_L(), *ip_data_->curr()->z_U(),
                                      *ip_data_->curr()->v_L(), *ip_data_->curr()->v_U(),
                                      s_max_,
                                      s_d, s_c);
        DBG_PRINT((1, "s_d = %lf, s_c = %lf\n", s_d, s_c));

        // Dual infeasibility
        DBG_PRINT((1, "curr_dual_infeasibility(NORM_MAX) = %8.2e\n",
                   curr_dual_infeasibility(NORM_MAX)));
        result = curr_dual_infeasibility(NORM_MAX)/s_d;
        /*
        // Primal infeasibility
        DBG_PRINT((1, "curr_primal_infeasibility(NORM_MAX) = %8.2e\n",
        curr_primal_infeasibility(NORM_MAX)));
        result = Max(result, curr_primal_infeasibility(NORM_MAX));
        */
        result = Max(result, curr_nlp_constraint_violation(NORM_MAX));
        // Complementarity
        DBG_PRINT((1, "curr_complementarity(0., NORM_MAX) = %8.2e\n",
                   curr_complementarity(0., NORM_MAX)));
        result = Max(result, curr_complementarity(0., NORM_MAX)/s_c);
      }

      curr_nlp_error_cache_.AddCachedResult(result, tdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::unscaled_curr_nlp_error()
  {
    DBG_START_METH("IpoptCalculatedQuantities::unscaled_curr_nlp_error()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    SmartPtr<const Vector> y_c = ip_data_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    SmartPtr<const Vector> z_L = ip_data_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr()->z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr()->v_U();

    std::vector<const TaggedObject*> tdeps(8);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    tdeps[2] = GetRawPtr(y_c);
    tdeps[3] = GetRawPtr(y_d);
    tdeps[4] = GetRawPtr(z_L);
    tdeps[5] = GetRawPtr(z_U);
    tdeps[6] = GetRawPtr(v_L);
    tdeps[7] = GetRawPtr(v_U);

    if (!unscaled_curr_nlp_error_cache_.GetCachedResult(result, tdeps)) {

      // Dual infeasibility
      result = unscaled_curr_dual_infeasibility(NORM_MAX);
      // Constraint violation
      result = Max(result, unscaled_curr_nlp_constraint_violation(NORM_MAX));
      // Complementarity (ToDo use unscaled?)
      DBG_PRINT((1, "curr_complementarity(0., NORM_MAX) = %8.2e\n",
                 curr_complementarity(0., NORM_MAX)));
      result = Max(result, unscaled_curr_complementarity(0., NORM_MAX));

      unscaled_curr_nlp_error_cache_.AddCachedResult(result, tdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_barrier_error()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_barrier_error()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    SmartPtr<const Vector> y_c = ip_data_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    SmartPtr<const Vector> z_L = ip_data_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr()->z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr()->v_U();
    Number mu = ip_data_->curr_mu();

    std::vector<const TaggedObject*> tdeps(8);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    tdeps[2] = GetRawPtr(y_c);
    tdeps[3] = GetRawPtr(y_d);
    tdeps[4] = GetRawPtr(z_L);
    tdeps[5] = GetRawPtr(z_U);
    tdeps[6] = GetRawPtr(v_L);
    tdeps[7] = GetRawPtr(v_U);
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;

    if (!curr_barrier_error_cache_.GetCachedResult(result, tdeps, sdeps)) {
      Number s_d = 0;
      Number s_c = 0;
      ComputeOptimalityErrorScaling(*ip_data_->curr()->y_c(), *ip_data_->curr()->y_d(),
                                    *ip_data_->curr()->z_L(), *ip_data_->curr()->z_U(),
                                    *ip_data_->curr()->v_L(), *ip_data_->curr()->v_U(),
                                    s_max_,
                                    s_d, s_c);
      DBG_PRINT((1, "s_d = %lf, s_c = %lf\n", s_d, s_c));

      // Primal infeasibility
      result = curr_dual_infeasibility(NORM_MAX)/s_d;
      // Dual infeasibility
      result = Max(result, curr_primal_infeasibility(NORM_MAX));
      // Complementarity
      result = Max(result, curr_complementarity(mu, NORM_MAX)/s_c);

      curr_barrier_error_cache_.AddCachedResult(result, tdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_primal_dual_system_error(Number mu)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_primal_dual_system_error()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    SmartPtr<const Vector> y_c = ip_data_->curr()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->curr()->y_d();
    SmartPtr<const Vector> z_L = ip_data_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr()->z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr()->v_U();

    std::vector<const TaggedObject*> tdeps(8);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    tdeps[2] = GetRawPtr(y_c);
    tdeps[3] = GetRawPtr(y_d);
    tdeps[4] = GetRawPtr(z_L);
    tdeps[5] = GetRawPtr(z_U);
    tdeps[6] = GetRawPtr(v_L);
    tdeps[7] = GetRawPtr(v_U);
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;

    if (!curr_primal_dual_system_error_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (!trial_primal_dual_system_error_cache_.GetCachedResult(result, tdeps, sdeps)) {
        // For now we use the 1 norm, and scale each component by the number of entries...
        Index n_dual = x->Dim() + s->Dim();
        Number dual_inf = curr_dual_infeasibility(NORM_1)/((Number)n_dual);

        Index n_primal = y_c->Dim() + y_d->Dim();
        Number primal_inf = 0.;
        if (n_primal>0) {
          primal_inf = curr_primal_infeasibility(NORM_1)/((Number)n_primal);
        }

        Index n_cmpl = z_L->Dim() + z_U->Dim() + v_L->Dim() + v_U->Dim();
        Number cmpl = 0.;
        if (n_cmpl>0) {
          cmpl = curr_complementarity(mu, NORM_1)/((Number)n_cmpl);
        }

        result = dual_inf + primal_inf + cmpl;
      }
      curr_primal_dual_system_error_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::trial_primal_dual_system_error(Number mu)
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_primal_dual_system_error()",
                   dbg_verbosity);
    DBG_ASSERT(initialize_called_);
    Number result;

    SmartPtr<const Vector> x = ip_data_->trial()->x();
    SmartPtr<const Vector> s = ip_data_->trial()->s();
    SmartPtr<const Vector> y_c = ip_data_->trial()->y_c();
    SmartPtr<const Vector> y_d = ip_data_->trial()->y_d();
    SmartPtr<const Vector> z_L = ip_data_->trial()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->trial()->z_U();
    SmartPtr<const Vector> v_L = ip_data_->trial()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->trial()->v_U();

    std::vector<const TaggedObject*> tdeps(8);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    tdeps[2] = GetRawPtr(y_c);
    tdeps[3] = GetRawPtr(y_d);
    tdeps[4] = GetRawPtr(z_L);
    tdeps[5] = GetRawPtr(z_U);
    tdeps[6] = GetRawPtr(v_L);
    tdeps[7] = GetRawPtr(v_U);
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;

    if (!trial_primal_dual_system_error_cache_.GetCachedResult(result, tdeps, sdeps)) {
      if (!curr_primal_dual_system_error_cache_.GetCachedResult(result, tdeps, sdeps)) {
        // For now we use the 1 norm, and scale each component by the number of entries...
        Index n_dual = x->Dim() + s->Dim();
        Number dual_inf = trial_dual_infeasibility(NORM_1)/((Number)n_dual);

        Index n_primal = y_c->Dim() + y_d->Dim();
        Number primal_inf = 0.;
        if (n_primal>0) {
          primal_inf = trial_primal_infeasibility(NORM_1)/((Number)n_primal);
        }

        Index n_cmpl = z_L->Dim() + z_U->Dim() + v_L->Dim() + v_U->Dim();
        Number cmpl = 0.;
        if (n_cmpl>0) {
          cmpl = trial_complementarity(mu, NORM_1)/((Number)n_cmpl);
        }

        result = dual_inf + primal_inf + cmpl;
      }
      trial_primal_dual_system_error_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  ///////////////////////////////////////////////////////////////////////////
  //                Fraction-to-the-boundary step sizes                    //
  ///////////////////////////////////////////////////////////////////////////

  Number
  IpoptCalculatedQuantities::CalcFracToBound(const Vector& slack_L,
      Vector& tmp_L,
      const Matrix& P_L,
      const Vector& slack_U,
      Vector& tmp_U,
      const Matrix& P_U,
      const Vector& delta,
      Number tau)
  {
    DBG_START_METH("IpoptCalculatedQuantities::CalcFracToBound",
                   dbg_verbosity);

    Number alpha_L = 1.0;
    Number alpha_U = 1.0;
    if (slack_L.Dim() > 0) {
      P_L.TransMultVector(1.0, delta, 0.0, tmp_L);
      alpha_L = slack_L.FracToBound(tmp_L, tau);
    }

    if (slack_U.Dim() > 0) {
      P_U.TransMultVector(-1.0, delta, 0.0, tmp_U);
      alpha_U = slack_U.FracToBound(tmp_U, tau);
    }

    DBG_PRINT((1,"alpha_L = %lf, alpha_U = %lf\n", alpha_L, alpha_U));
    DBG_ASSERT(alpha_L >= 0.0 && alpha_L <= 1.0
               && alpha_U >=0.0 && alpha_U <= 1.0);

    return Min(alpha_L, alpha_U);
  }

  Number
  IpoptCalculatedQuantities::primal_frac_to_the_bound(Number tau,
      const Vector& delta_x,
      const Vector& delta_s)
  {
    DBG_START_METH("IpoptCalculatedQuantities::primal_frac_to_the_bound",
                   dbg_verbosity);
    Number result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    std::vector<const TaggedObject*> tdeps(4);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    tdeps[2] = &delta_x;
    tdeps[3] = &delta_s;

    std::vector<Number> sdeps(1);
    sdeps[0] = tau;

    if (!primal_frac_to_the_bound_cache_.GetCachedResult(result, tdeps, sdeps)) {
      result = Min(CalcFracToBound(*curr_slack_x_L(), Tmp_x_L(), *ip_nlp_->Px_L(),
                                   *curr_slack_x_U(), Tmp_x_U(), *ip_nlp_->Px_U(),
                                   delta_x, tau),
                   CalcFracToBound(*curr_slack_s_L(), Tmp_s_L(), *ip_nlp_->Pd_L(),
                                   *curr_slack_s_U(), Tmp_s_U(), *ip_nlp_->Pd_U(),
                                   delta_s, tau));

      primal_frac_to_the_bound_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_primal_frac_to_the_bound(Number tau)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_primal_frac_to_the_bound()",
                   dbg_verbosity);
    return primal_frac_to_the_bound(tau, *ip_data_->delta()->x(),
                                    *ip_data_->delta()->s());
  }

  Number
  IpoptCalculatedQuantities::uncached_dual_frac_to_the_bound(
    Number tau,
    const Vector& delta_z_L,
    const Vector& delta_z_U,
    const Vector& delta_v_L,
    const Vector& delta_v_U)
  {
    DBG_START_METH("IpoptCalculatedQuantities::uncached_dual_frac_to_the_bound",
                   dbg_verbosity);
    Number result;

    result = ip_data_->curr()->z_L()->FracToBound(delta_z_L, tau);
    result = Min(result, ip_data_->curr()->z_U()->FracToBound(delta_z_U, tau));
    result = Min(result, ip_data_->curr()->v_L()->FracToBound(delta_v_L, tau));
    result = Min(result, ip_data_->curr()->v_U()->FracToBound(delta_v_U, tau));

    return result;
  }

  Number
  IpoptCalculatedQuantities::dual_frac_to_the_bound(
    Number tau,
    const Vector& delta_z_L,
    const Vector& delta_z_U,
    const Vector& delta_v_L,
    const Vector& delta_v_U)
  {
    DBG_START_METH("IpoptCalculatedQuantities::dual_frac_to_the_bound",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> z_L = ip_data_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr()->z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr()->v_U();
    std::vector<const TaggedObject*> tdeps(8);
    tdeps[0] = GetRawPtr(z_L);
    tdeps[1] = GetRawPtr(z_U);
    tdeps[2] = GetRawPtr(v_L);
    tdeps[3] = GetRawPtr(v_U);
    tdeps[4] = &delta_z_L;
    tdeps[5] = &delta_z_U;
    tdeps[6] = &delta_v_L;
    tdeps[7] = &delta_v_U;

    std::vector<Number> sdeps(1);
    sdeps[0] = tau;

    if (!dual_frac_to_the_bound_cache_.GetCachedResult(result, tdeps, sdeps)) {
      result = z_L->FracToBound(delta_z_L, tau);
      result = Min(result, z_U->FracToBound(delta_z_U, tau));
      result = Min(result, v_L->FracToBound(delta_v_L, tau));
      result = Min(result, v_U->FracToBound(delta_v_U, tau));

      dual_frac_to_the_bound_cache_.AddCachedResult(result, tdeps, sdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_dual_frac_to_the_bound(Number tau)
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_dual_frac_to_the_bound()",
                   dbg_verbosity);
    return dual_frac_to_the_bound(tau, *ip_data_->delta()->z_L(),
                                  *ip_data_->delta()->z_U(),
                                  *ip_data_->delta()->v_L(),
                                  *ip_data_->delta()->v_U());
  }

  Number
  IpoptCalculatedQuantities::uncached_slack_frac_to_the_bound(
    Number tau,
    const Vector& delta_x_L,
    const Vector& delta_x_U,
    const Vector& delta_s_L,
    const Vector& delta_s_U)
  {
    DBG_START_METH("IpoptCalculatedQuantities::slack_frac_to_the_bound",
                   dbg_verbosity);
    Number result;

    SmartPtr<const Vector> x_L = curr_slack_x_L();
    SmartPtr<const Vector> x_U = curr_slack_x_U();
    SmartPtr<const Vector> s_L = curr_slack_s_L();
    SmartPtr<const Vector> s_U = curr_slack_s_U();

    result = x_L->FracToBound(delta_x_L, tau);
    result = Min(result, x_U->FracToBound(delta_x_U, tau));
    result = Min(result, s_L->FracToBound(delta_s_L, tau));
    result = Min(result, s_U->FracToBound(delta_s_U, tau));

    return result;
  }

  ///////////////////////////////////////////////////////////////////////////
  //                             Sigma Matrices                            //
  ///////////////////////////////////////////////////////////////////////////

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_sigma_x()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_sigma_x()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> z_L = ip_data_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr()->z_U();

    if (!curr_sigma_x_cache_.GetCachedResult3Dep(result, *x, *z_L, *z_U)) {
      SmartPtr<Vector> sigma = x->MakeNew();

      sigma->Set(0.);
      ip_nlp_->Px_L()->AddMSinvZ(1., *curr_slack_x_L(), *z_L, *sigma);
      ip_nlp_->Px_U()->AddMSinvZ(1., *curr_slack_x_U(), *z_U, *sigma);

      DBG_PRINT_VECTOR(2,"sigma_x", *sigma);

      result = ConstPtr(sigma);
      curr_sigma_x_cache_.AddCachedResult3Dep(result, *x, *z_L, *z_U);
    }

    return result;
  }

  SmartPtr<const Vector>
  IpoptCalculatedQuantities::curr_sigma_s()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_sigma_s()",
                   dbg_verbosity);
    SmartPtr<const Vector> result;
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    SmartPtr<const Vector> v_L = ip_data_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr()->v_U();

    if (!curr_sigma_s_cache_.GetCachedResult3Dep(result, *s, *v_L, *v_U)) {
      SmartPtr<Vector> sigma = s->MakeNew();

      sigma->Set(0.);
      ip_nlp_->Pd_L()->AddMSinvZ(1., *curr_slack_s_L(), *v_L, *sigma);
      ip_nlp_->Pd_U()->AddMSinvZ(1., *curr_slack_s_U(), *v_U, *sigma);

      result = ConstPtr(sigma);
      curr_sigma_s_cache_.AddCachedResult3Dep(result, *s, *v_L, *v_U);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::curr_avrg_compl()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_avrg_compl()",
                   dbg_verbosity);

    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    SmartPtr<const Vector> z_L = ip_data_->curr()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->curr()->z_U();
    SmartPtr<const Vector> v_L = ip_data_->curr()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->curr()->v_U();

    std::vector<const TaggedObject*> tdeps(6);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    tdeps[2] = GetRawPtr(z_L);
    tdeps[3] = GetRawPtr(z_U);
    tdeps[4] = GetRawPtr(v_L);
    tdeps[5] = GetRawPtr(v_U);

    if (!curr_avrg_compl_cache_.GetCachedResult(result, tdeps)) {
      if (!trial_avrg_compl_cache_.GetCachedResult(result, tdeps)) {

        SmartPtr<const Vector> slack_x_L = curr_slack_x_L();
        SmartPtr<const Vector> slack_x_U = curr_slack_x_U();
        SmartPtr<const Vector> slack_s_L = curr_slack_s_L();
        SmartPtr<const Vector> slack_s_U = curr_slack_s_U();

        Index ncomps = z_L->Dim() + z_U->Dim() + v_L->Dim() + v_U->Dim();

        if (ncomps>0) {
          result = z_L->Dot(*slack_x_L);
          result += z_U->Dot(*slack_x_U);
          result += v_L->Dot(*slack_s_L);
          result += v_U->Dot(*slack_s_U);

          result /= (Number)ncomps;
        }
        else {
          result = 0.;
        }
      }

      curr_avrg_compl_cache_.AddCachedResult(result, tdeps);
    }

    return result;
  }

  Number
  IpoptCalculatedQuantities::trial_avrg_compl()
  {
    DBG_START_METH("IpoptCalculatedQuantities::trial_avrg_compl()",
                   dbg_verbosity);

    Number result;

    SmartPtr<const Vector> x = ip_data_->trial()->x();
    SmartPtr<const Vector> s = ip_data_->trial()->s();
    SmartPtr<const Vector> z_L = ip_data_->trial()->z_L();
    SmartPtr<const Vector> z_U = ip_data_->trial()->z_U();
    SmartPtr<const Vector> v_L = ip_data_->trial()->v_L();
    SmartPtr<const Vector> v_U = ip_data_->trial()->v_U();

    std::vector<const TaggedObject*> tdeps(6);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    tdeps[2] = GetRawPtr(z_L);
    tdeps[3] = GetRawPtr(z_U);
    tdeps[4] = GetRawPtr(v_L);
    tdeps[5] = GetRawPtr(v_U);

    if (!trial_avrg_compl_cache_.GetCachedResult(result, tdeps)) {
      if (!curr_avrg_compl_cache_.GetCachedResult(result, tdeps)) {

        SmartPtr<const Vector> slack_x_L = trial_slack_x_L();
        SmartPtr<const Vector> slack_x_U = trial_slack_x_U();
        SmartPtr<const Vector> slack_s_L = trial_slack_s_L();
        SmartPtr<const Vector> slack_s_U = trial_slack_s_U();

        Index ncomps = z_L->Dim() + z_U->Dim() + v_L->Dim() + v_U->Dim();

        if (ncomps>0) {
          result = z_L->Dot(*slack_x_L);
          result += z_U->Dot(*slack_x_U);
          result += v_L->Dot(*slack_s_L);
          result += v_U->Dot(*slack_s_U);

          result /= (Number)ncomps;
        }
        else {
          result = 0.;
        }
      }

      trial_avrg_compl_cache_.AddCachedResult(result, tdeps);
    }

    return result;
  }

  void IpoptCalculatedQuantities::ComputeOptimalityErrorScaling(const Vector& y_c, const Vector& y_d,
      const Vector& z_L, const Vector& z_U,
      const Vector& v_L, const Vector& v_U,
      Number s_max,
      Number& s_d, Number& s_c)
  {
    DBG_ASSERT(initialize_called_);

    s_c = z_L.Asum() + z_U.Asum() + v_L.Asum() + v_U.Asum();
    Number n = (z_L.Dim() + z_U.Dim() + v_L.Dim() + v_U.Dim());
    if (n == 0) {
      s_c = 1.0;
    }
    else {
      s_c = s_c / n;
      s_c = Max(s_max, s_c)/s_max;
    }

    s_d = y_c.Asum() + y_d.Asum() + z_L.Asum() + z_U.Asum() + v_L.Asum() + v_U.Asum();
    n = (y_c.Dim() + y_d.Dim() + z_L.Dim() + z_U.Dim() + v_L.Dim() + v_U.Dim());
    if ( n == 0 ) {
      s_d = 1.0;
    }
    else {
      s_d = s_d / n;
      s_d = Max(s_max, s_d)/s_max;
    }
  }

  Number IpoptCalculatedQuantities::curr_gradBarrTDelta()
  {
    DBG_START_METH("IpoptCalculatedQuantities::curr_gradBarrTDelta()",
                   dbg_verbosity);

    Number result;

    SmartPtr<const Vector> x = ip_data_->curr()->x();
    SmartPtr<const Vector> s = ip_data_->curr()->s();
    SmartPtr<const Vector> delta_x = ip_data_->delta()->x();
    SmartPtr<const Vector> delta_s = ip_data_->delta()->s();
    std::vector<const TaggedObject*> tdeps(4);
    tdeps[0] = GetRawPtr(x);
    tdeps[1] = GetRawPtr(s);
    tdeps[2] = GetRawPtr(delta_x);
    tdeps[3] = GetRawPtr(delta_s);
    Number mu = ip_data_->curr_mu();
    std::vector<Number> sdeps(1);
    sdeps[0] = mu;
    DBG_PRINT((1,"curr_mu=%e\n",mu));

    if (!curr_gradBarrTDelta_cache_.GetCachedResult(result, tdeps, sdeps)) {
      result = curr_grad_barrier_obj_x()->Dot(*delta_x) +
               curr_grad_barrier_obj_s()->Dot(*delta_s);

      curr_gradBarrTDelta_cache_.AddCachedResult(result, tdeps, sdeps);
    }
    return result;
  }

  bool IpoptCalculatedQuantities::IsSquareProblem() const
  {
    return (ip_data_->curr()->x()->Dim() == ip_data_->curr()->y_c()->Dim());
  }

  Vector& IpoptCalculatedQuantities::Tmp_x()
  {
    if (!IsValid(tmp_x_)) {
      tmp_x_ = ip_data_->curr()->x()->MakeNew();
    }
    return *tmp_x_;
  }

  Vector& IpoptCalculatedQuantities::Tmp_s()
  {
    if (!IsValid(tmp_s_)) {
      tmp_s_ = ip_data_->curr()->s()->MakeNew();
    }
    return *tmp_s_;
  }

  Vector& IpoptCalculatedQuantities::Tmp_c()
  {
    if (!IsValid(tmp_c_)) {
      tmp_c_ = ip_data_->curr()->y_c()->MakeNew();
    }
    return *tmp_c_;
  }

  Vector& IpoptCalculatedQuantities::Tmp_d()
  {
    if (!IsValid(tmp_d_)) {
      tmp_d_ = ip_data_->curr()->y_d()->MakeNew();
    }
    return *tmp_d_;
  }

  Vector& IpoptCalculatedQuantities::Tmp_x_L()
  {
    if (!IsValid(tmp_x_L_)) {
      tmp_x_L_ = ip_nlp_->x_L()->MakeNew();
    }
    return *tmp_x_L_;
  }

  Vector& IpoptCalculatedQuantities::Tmp_x_U()
  {
    if (!IsValid(tmp_x_U_)) {
      tmp_x_U_ = ip_nlp_->x_U()->MakeNew();
    }
    return *tmp_x_U_;
  }

  Vector& IpoptCalculatedQuantities::Tmp_s_L()
  {
    if (!IsValid(tmp_s_L_)) {
      tmp_s_L_ = ip_nlp_->d_L()->MakeNew();
    }
    return *tmp_s_L_;
  }

  Vector& IpoptCalculatedQuantities::Tmp_s_U()
  {
    if (!IsValid(tmp_s_U_)) {
      tmp_s_U_ = ip_nlp_->d_U()->MakeNew();
    }
    return *tmp_s_U_;
  }


} // namespace Ipopt
