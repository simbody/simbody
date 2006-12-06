// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpProbingMuOracle.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpProbingMuOracle.hpp"

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

  ProbingMuOracle::ProbingMuOracle(const SmartPtr<PDSystemSolver>& pd_solver)
      :
      MuOracle(),
      pd_solver_(pd_solver)
  {
    DBG_ASSERT(IsValid(pd_solver_));
  }

  ProbingMuOracle::~ProbingMuOracle()
  {}

  void ProbingMuOracle::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    // None to register...
  }

  bool ProbingMuOracle::InitializeImpl(const OptionsList& options,
                                       const std::string& prefix)
  {
    options.GetNumericValue("sigma_max", sigma_max_, prefix);

    return true;
  }

  bool ProbingMuOracle::CalculateMu(Number mu_min, Number mu_max,
                                    Number& new_mu)
  {
    DBG_START_METH("ProbingMuOracle::CalculateMu",
                   dbg_verbosity);

    /////////////////////////////////////
    // Compute the affine scaling step //
    /////////////////////////////////////

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "Solving the Primal Dual System for the affine step\n");
    // First get the right hand side
    SmartPtr<IteratesVector> rhs = IpData().curr()->MakeNewContainer();

    rhs->Set_x(*IpCq().curr_grad_lag_x());
    rhs->Set_s(*IpCq().curr_grad_lag_s());
    rhs->Set_y_c(*IpCq().curr_c());
    rhs->Set_y_d(*IpCq().curr_d_minus_s());
    rhs->Set_z_L(*IpCq().curr_compl_x_L());
    rhs->Set_z_U(*IpCq().curr_compl_x_U());
    rhs->Set_v_L(*IpCq().curr_compl_s_L());
    rhs->Set_v_U(*IpCq().curr_compl_s_U());

    // Get space for the affine scaling step
    SmartPtr<IteratesVector> step = rhs->MakeNewIteratesVector(true);

    // Now solve the primal-dual system to get the affine step.  We
    // allow a somewhat inexact solution here
    bool allow_inexact = true;
    bool retval = pd_solver_->Solve(-1.0, 0.0,
                                    *rhs,
                                    *step,
                                    allow_inexact
                                   );
    if (!retval) {
      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "The linear system could not be solved for the affine step!\n");
      return false;
    }

    DBG_PRINT_VECTOR(2, "step", *step);

    /////////////////////////////////////////////////////////////
    // Use Mehrotra's formula to compute the barrier parameter //
    /////////////////////////////////////////////////////////////

    // First compute the fraction-to-the-boundary step sizes
    Number alpha_primal_aff = IpCq().primal_frac_to_the_bound(1.0,
                              *step->x(),
                              *step->s());

    Number alpha_dual_aff = IpCq().dual_frac_to_the_bound(1.0,
                            *step->z_L(),
                            *step->z_U(),
                            *step->v_L(),
                            *step->v_U());

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "  The affine maximal step sizes are\n"
                   "   alpha_primal_aff = %23.16e\n"
                   "   alpha_dual_aff = %23.16e\n",
                   alpha_primal_aff,
                   alpha_dual_aff);

    // now compute the average complementarity at the affine step
    // ToDo shoot for mu_min instead of 0?
    Number mu_aff = CalculateAffineMu(alpha_primal_aff, alpha_dual_aff, *step);
    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "  The average complementariy at the affine step is %23.16e\n",
                   mu_aff);

    // get the current average complementarity
    Number mu_curr = IpCq().curr_avrg_compl();
    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "  The average complementariy at the current point is %23.16e\n",
                   mu_curr);
    DBG_ASSERT(mu_curr>0.);

    // Apply Mehrotra's rule
    Number sigma = pow((mu_aff/mu_curr),3);
    // Make sure, sigma is not too large
    sigma = Min(sigma, sigma_max_);

    Number mu = sigma*mu_curr;

    // Store the affine search direction (in case it is needed in the
    // line search for a corrector step)
    IpData().set_delta_aff(step);
    IpData().SetHaveAffineDeltas(true);

    char ssigma[40];
    sprintf(ssigma, " sigma=%8.2e", sigma);
    IpData().Append_info_string(ssigma);
    //sprintf(ssigma, " xi=%8.2e ", IpCq().curr_centrality_measure());
    //IpData().Append_info_string(ssigma);

    new_mu = Max(Min(mu, mu_max), mu_min);
    return true;
  }

  Number ProbingMuOracle::CalculateAffineMu
  (
    Number alpha_primal,
    Number alpha_dual,
    const IteratesVector& step)
  {
    // Get the current values of the slack variables and bound multipliers
    SmartPtr<const Vector> slack_x_L = IpCq().curr_slack_x_L();
    SmartPtr<const Vector> slack_x_U = IpCq().curr_slack_x_U();
    SmartPtr<const Vector> slack_s_L = IpCq().curr_slack_s_L();
    SmartPtr<const Vector> slack_s_U = IpCq().curr_slack_s_U();

    SmartPtr<const Vector> z_L = IpData().curr()->z_L();
    SmartPtr<const Vector> z_U = IpData().curr()->z_U();
    SmartPtr<const Vector> v_L = IpData().curr()->v_L();
    SmartPtr<const Vector> v_U = IpData().curr()->v_U();

    SmartPtr<Vector> tmp_slack;
    SmartPtr<Vector> tmp_mult;
    SmartPtr<const Matrix> P;
    Index ncomp = 0;
    Number sum =0.;

    // For each combination of slack and multiplier, compute the new
    // values and their dot products.

    // slack_x_L
    if (slack_x_L->Dim()>0) {
      ncomp += slack_x_L->Dim();

      P = IpNLP().Px_L();
      tmp_slack = slack_x_L->MakeNew();
      tmp_slack->Copy(*slack_x_L);
      P->TransMultVector(alpha_primal, *step.x(), 1.0, *tmp_slack);

      tmp_mult = z_L->MakeNew();
      tmp_mult->Copy(*z_L);
      tmp_mult->Axpy(alpha_dual, *step.z_L());

      sum += tmp_slack->Dot(*tmp_mult);
    }

    // slack_x_U
    if (slack_x_U->Dim()>0) {
      ncomp += slack_x_U->Dim();

      P = IpNLP().Px_U();
      tmp_slack = slack_x_U->MakeNew();
      tmp_slack->Copy(*slack_x_U);
      P->TransMultVector(-alpha_primal, *step.x(), 1.0, *tmp_slack);

      tmp_mult = z_U->MakeNew();
      tmp_mult->Copy(*z_U);
      tmp_mult->Axpy(alpha_dual, *step.z_U());

      sum += tmp_slack->Dot(*tmp_mult);
    }

    // slack_s_L
    if (slack_s_L->Dim()>0) {
      ncomp += slack_s_L->Dim();

      P = IpNLP().Pd_L();
      tmp_slack = slack_s_L->MakeNew();
      tmp_slack->Copy(*slack_s_L);
      P->TransMultVector(alpha_primal, *step.s(), 1.0, *tmp_slack);

      tmp_mult = v_L->MakeNew();
      tmp_mult->Copy(*v_L);
      tmp_mult->Axpy(alpha_dual, *step.v_L());

      sum += tmp_slack->Dot(*tmp_mult);
    }

    // slack_s_U
    if (slack_s_U->Dim()>0) {
      ncomp += slack_s_U->Dim();

      P = IpNLP().Pd_U();
      tmp_slack = slack_s_U->MakeNew();
      tmp_slack->Copy(*slack_s_U);
      P->TransMultVector(-alpha_primal, *step.s(), 1.0, *tmp_slack);

      tmp_mult = v_U->MakeNew();
      tmp_mult->Copy(*v_U);
      tmp_mult->Axpy(alpha_dual, *step.v_U());

      sum += tmp_slack->Dot(*tmp_mult);
    }

    DBG_ASSERT(ncomp>0);

    return sum/((Number)ncomp);
  }

} // namespace Ipopt
