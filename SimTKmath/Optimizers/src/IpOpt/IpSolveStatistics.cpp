// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpSolveStatistics.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter          IBM    2005-08-15

#include "IpSolveStatistics.hpp"
#include "IpIpoptCalculatedQuantities.hpp"

namespace Ipopt
{

  SolveStatistics::SolveStatistics(
    const SmartPtr<IpoptNLP>& ip_nlp,
    const SmartPtr<IpoptData>& ip_data,
    const SmartPtr<IpoptCalculatedQuantities>& ip_cq)
      :
      num_iters_(ip_data->iter_count()),
      total_cpu_time_(ip_data->TimingStats().OverallAlgorithm().TotalTime()),
      num_obj_evals_(ip_nlp->f_evals()),
      num_constr_evals_(Max(ip_nlp->c_evals(), ip_nlp->d_evals())),
      num_obj_grad_evals_(ip_nlp->grad_f_evals()),
      num_constr_jac_evals_(Max(ip_nlp->jac_c_evals(),ip_nlp->jac_c_evals())),
      num_hess_evals_(ip_nlp->h_evals()),

      scaled_obj_val_(ip_cq->curr_f()),
      obj_val_(ip_cq->unscaled_curr_f()),
      scaled_dual_inf_(ip_cq->curr_dual_infeasibility(NORM_MAX)),
      dual_inf_(ip_cq->unscaled_curr_dual_infeasibility(NORM_MAX)),
      scaled_constr_viol_(ip_cq->curr_nlp_constraint_violation(NORM_MAX)),
      constr_viol_(ip_cq->unscaled_curr_nlp_constraint_violation(NORM_MAX)),
      scaled_compl_(ip_cq->curr_complementarity(0., NORM_MAX)),
      compl_(ip_cq->unscaled_curr_complementarity(0., NORM_MAX)),
      scaled_kkt_error_(ip_cq->curr_nlp_error()),
      kkt_error_(ip_cq->unscaled_curr_nlp_error())
  {}

  Index SolveStatistics::IterationCount() const
  {
    return num_iters_;
  }

  Number SolveStatistics::TotalCPUTime() const
  {
    return total_cpu_time_;
  }

  void SolveStatistics::NumberOfEvaluations(
    Index& num_obj_evals,
    Index& num_constr_evals,
    Index& num_obj_grad_evals,
    Index& num_constr_jac_evals,
    Index& num_hess_evals) const
  {
    num_obj_evals = num_obj_evals_;
    num_constr_evals = num_constr_evals_;
    num_obj_grad_evals = num_obj_grad_evals_;
    num_constr_jac_evals = num_constr_jac_evals_;
    num_hess_evals = num_hess_evals_;
  }

  void SolveStatistics::Infeasibilities(Number& dual_inf,
                                        Number& constr_viol,
                                        Number& complementarity,
                                        Number& kkt_error) const
  {
    dual_inf = dual_inf_;
    constr_viol = constr_viol_;
    complementarity = compl_;
    kkt_error = kkt_error_;
  }

  void SolveStatistics::ScaledInfeasibilities(Number& scaled_dual_inf,
      Number& scaled_constr_viol,
      Number& scaled_complementarity,
      Number& scaled_kkt_error) const
  {
    scaled_dual_inf = scaled_dual_inf_;
    scaled_constr_viol = scaled_constr_viol_;
    scaled_complementarity = scaled_compl_;
    scaled_kkt_error = scaled_kkt_error_;
  }

  Number SolveStatistics::FinalObjective() const
  {
    return obj_val_;
  }

  Number SolveStatistics::FinalScaledObjective() const
  {
    return scaled_obj_val_;
  }

} // namespace Ipopt
