// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpQualityFunctionMuOracle.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter            IBM    2004-11-12

#include "IpQualityFunctionMuOracle.hpp"

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

  QualityFunctionMuOracle::QualityFunctionMuOracle(const SmartPtr<PDSystemSolver>& pd_solver)
      :
      MuOracle(),
      pd_solver_(pd_solver),

      tmp_step_x_L_(NULL),
      tmp_step_x_U_(NULL),
      tmp_step_s_L_(NULL),
      tmp_step_s_U_(NULL),
      tmp_step_z_L_(NULL),
      tmp_step_z_U_(NULL),
      tmp_step_v_L_(NULL),
      tmp_step_v_U_(NULL),

      tmp_slack_x_L_(NULL),
      tmp_slack_x_U_(NULL),
      tmp_slack_s_L_(NULL),
      tmp_slack_s_U_(NULL),
      tmp_z_L_(NULL),
      tmp_z_U_(NULL),
      tmp_v_L_(NULL),
      tmp_v_U_(NULL),

      count_qf_evals_(0)
  {
    DBG_ASSERT(IsValid(pd_solver_));
  }

  QualityFunctionMuOracle::~QualityFunctionMuOracle()
  {}

  void QualityFunctionMuOracle::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "sigma_max",
      "Maximum value of the centering parameter.",
      0.0, true, 1e2,
      "This is the upper bound for the centering parameter chosen by the "
      "quality function based barrier parameter update. (Only used if option "
      "\"mu_oracle\" is set to \"quality-function\".)");
    roptions->AddLowerBoundedNumberOption(
      "sigma_min",
      "Minimum value of the centering parameter.",
      0.0, false, 1e-6,
      "This is the lower bound for the centering parameter chosen by the "
      "quality function based barrier parameter update. (Only used if option "
      "\"mu_oracle\" is set to \"quality-function\".)");
    roptions->AddStringOption4(
      "quality_function_norm_type",
      "Norm used for components of the quality function.",
      "2-norm-squared",
      "1-norm", "use the 1-norm (abs sum)",
      "2-norm-squared", "use the 2-norm squared (sum of squares)",
      "max-norm", "use the infinity norm (max)",
      "2-norm", "use 2-norm",
      "(Only used if option \"mu_oracle\" is set to \"quality-function\".)");
    roptions->AddStringOption4(
      "quality_function_centrality",
      "The penalty term for centrality that is included in quality function.",
      "none",
      "none", "no penalty term is added",
      "log", "complementarity * the log of the centrality measure",
      "reciprocal", "complementarity * the reciprocal of the centrality measure",
      "cubed-reciprocal", "complementarity * the reciprocal of the centrality measure cubed",
      "This determines whether a term is added to the quality function to "
      "penalize deviation from centrality with respect to complementarity.  The "
      "complementarity measure here is the xi in the Loqo update rule. (Only used if option "
      "\"mu_oracle\" is set to \"quality-function\".)");
    roptions->AddStringOption2(
      "quality_function_balancing_term",
      "The balancing term included in the quality function for centrality.",
      "none",
      "none", "no balancing term is added",
      "cubic", "Max(0,Max(dual_inf,primal_inf)-compl)^3",
      "This determines whether a term is added to the quality function that "
      "penalizes situations where the complementarity is much smaller "
      "than dual and primal infeasibilities. (Only used if option "
      "\"mu_oracle\" is set to \"quality-function\".)");
    roptions->AddLowerBoundedIntegerOption(
      "quality_function_max_section_steps",
      "Maximum number of search steps during direct search procedure "
      "determining the optimal centering parameter.",
      0, 8,
      "The golden section search is performed for the quality function based "
      "mu oracle. (Only used if option "
      "\"mu_oracle\" is set to \"quality-function\".)");
    roptions->AddBoundedNumberOption(
      "quality_function_section_sigma_tol",
      "Tolerance for the section search procedure determining "
      "the optimal centering parameter (in sigma space).",
      0.0, false, 1.0, true,
      1e-2,
      "The golden section search is performed for the quality function based "
      "mu oractle. (Only used if option "
      "\"mu_oracle\" is set to \"quality-function\".)");
    roptions->AddBoundedNumberOption(
      "quality_function_section_qf_tol",
      "Tolerance for the golden section search procedure determining "
      "the optimal centering parameter (in the function value space).",
      0.0, false, 1.0, true,
      0e-2,
      "The golden section search is performed for the quality function based mu "
      "oractle. (Only used if option "
      "\"mu_oracle\" is set to \"quality-function\".)");
  }


  bool QualityFunctionMuOracle::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    Index enum_int;

    options.GetNumericValue("sigma_max", sigma_max_, prefix);
    options.GetNumericValue("sigma_min", sigma_min_, prefix);

    options.GetEnumValue("quality_function_norm_type", enum_int, prefix);
    quality_function_norm_ = NormEnum(enum_int);
    options.GetEnumValue("quality_function_centrality", enum_int, prefix);
    quality_function_centrality_ = CentralityEnum(enum_int);
    options.GetEnumValue("quality_function_balancing_term", enum_int, prefix);
    quality_function_balancing_term_ = BalancingTermEnum(enum_int);
    options.GetIntegerValue("quality_function_max_section_steps",
                            quality_function_max_section_steps_, prefix);
    options.GetNumericValue("quality_function_section_sigma_tol",
                            quality_function_section_sigma_tol_, prefix);
    options.GetNumericValue("quality_function_section_qf_tol",
                            quality_function_section_qf_tol_, prefix);

    initialized_ = false;

    return true;
  }

  bool QualityFunctionMuOracle::CalculateMu(Number mu_min, Number mu_max,
      Number& new_mu)
  {
    DBG_START_METH("QualityFunctionMuOracle::CalculateMu",
                   dbg_verbosity);

    ///////////////////////////////////////////////////////////////////////////
    // Reserve memory for temporary vectors used in CalculateQualityFunction //
    ///////////////////////////////////////////////////////////////////////////

    tmp_step_x_L_ = IpNLP().x_L()->MakeNew();
    tmp_step_x_U_ = IpNLP().x_U()->MakeNew();
    tmp_step_s_L_ = IpNLP().d_L()->MakeNew();
    tmp_step_s_U_ = IpNLP().d_U()->MakeNew();
    tmp_step_z_L_ = IpNLP().x_L()->MakeNew();
    tmp_step_z_U_ = IpNLP().x_U()->MakeNew();
    tmp_step_v_L_ = IpNLP().d_L()->MakeNew();
    tmp_step_v_U_ = IpNLP().d_U()->MakeNew();

    tmp_slack_x_L_ = IpNLP().x_L()->MakeNew();
    tmp_slack_x_U_ = IpNLP().x_U()->MakeNew();
    tmp_slack_s_L_ = IpNLP().d_L()->MakeNew();
    tmp_slack_s_U_ = IpNLP().d_U()->MakeNew();
    tmp_z_L_ = IpNLP().x_L()->MakeNew();
    tmp_z_U_ = IpNLP().x_U()->MakeNew();
    tmp_v_L_ = IpNLP().d_L()->MakeNew();
    tmp_v_U_ = IpNLP().d_U()->MakeNew();

    /////////////////////////////////////
    // Compute the affine scaling step //
    /////////////////////////////////////

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "Solving the Primal Dual System for the affine step\n");
    // First get the right hand side
    SmartPtr<IteratesVector> rhs_aff = IpData().curr()->MakeNewIteratesVector(false);
    rhs_aff->Set_x(*IpCq().curr_grad_lag_x());
    rhs_aff->Set_s(*IpCq().curr_grad_lag_s());
    rhs_aff->Set_y_c(*IpCq().curr_c());
    rhs_aff->Set_y_d(*IpCq().curr_d_minus_s());
    rhs_aff->Set_z_L(*IpCq().curr_compl_x_L());
    rhs_aff->Set_z_U(*IpCq().curr_compl_x_U());
    rhs_aff->Set_v_L(*IpCq().curr_compl_s_L());
    rhs_aff->Set_v_U(*IpCq().curr_compl_s_U());

    // Get space for the affine scaling step
    SmartPtr<IteratesVector> step_aff = IpData().curr()->MakeNewIteratesVector(true);

    // Now solve the primal-dual system to get the step.  We allow a
    // somewhat inexact solution, iterative refinement will be done
    // after mu is known
    bool allow_inexact = true;
    bool retval = pd_solver_->Solve(-1.0, 0.0,
                                    *rhs_aff,
                                    *step_aff,
                                    allow_inexact
                                   );
    if (!retval) {
      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "The linear system could not be solved for the affine step!\n");
      return false;
    }

    DBG_PRINT_VECTOR(2, "step_aff", *step_aff);

    /////////////////////////////////////
    // Compute the pure centering step //
    /////////////////////////////////////

    Number avrg_compl = IpCq().curr_avrg_compl();

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "Solving the Primal Dual System for the centering step\n");
    // First get the right hand side
    SmartPtr<IteratesVector> rhs_cen = IpData().curr()->MakeNewIteratesVector(true);
    rhs_cen->x_NonConst()->AddOneVector(-avrg_compl,
                                        *IpCq().grad_kappa_times_damping_x(),
                                        0.);
    rhs_cen->s_NonConst()->AddOneVector(-avrg_compl,
                                        *IpCq().grad_kappa_times_damping_s(),
                                        0.);

    rhs_cen->y_c_NonConst()->Set(0.);
    rhs_cen->y_d_NonConst()->Set(0.);
    rhs_cen->z_L_NonConst()->Set(avrg_compl);
    rhs_cen->z_U_NonConst()->Set(avrg_compl);
    rhs_cen->v_L_NonConst()->Set(avrg_compl);
    rhs_cen->v_U_NonConst()->Set(avrg_compl);

    // Get space for the centering step
    SmartPtr<IteratesVector> step_cen = IpData().curr()->MakeNewIteratesVector(true);

    // Now solve the primal-dual system to get the step
    allow_inexact = true;
    retval = pd_solver_->Solve(1.0, 0.0,
                               *rhs_cen,
                               *step_cen,
                               allow_inexact
                              );
    if (!retval) {
      Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                     "The linear system could not be solved for the centering step!\n");
      return false;
    }

    DBG_PRINT_VECTOR(2, "step_cen", *step_cen);

    // Start the timing for the quality function search here
    IpData().TimingStats().QualityFunctionSearch().Start();

    // Some initializations
    if (!initialized_) {
      n_dual_ = IpData().curr()->x()->Dim() + IpData().curr()->s()->Dim();
      n_pri_ = IpData().curr()->y_c()->Dim() + IpData().curr()->y_d()->Dim();
      n_comp_ = IpData().curr()->z_L()->Dim() + IpData().curr()->z_U()->Dim() +
                IpData().curr()->v_L()->Dim() + IpData().curr()->v_U()->Dim();
      initialized_ = true;
    }

    count_qf_evals_ = 0;

    // Compute some quantities used for the quality function evaluations
    // (This way we try to avoid retrieving numbers from cache...

    curr_slack_x_L_ = IpCq().curr_slack_x_L();
    curr_slack_x_U_ = IpCq().curr_slack_x_U();
    curr_slack_s_L_ = IpCq().curr_slack_s_L();
    curr_slack_s_U_ = IpCq().curr_slack_s_U();

    curr_z_L_ = IpData().curr()->z_L();
    curr_z_U_ = IpData().curr()->z_U();
    curr_v_L_ = IpData().curr()->v_L();
    curr_v_U_ = IpData().curr()->v_U();

    IpData().TimingStats().Task5().Start();
    switch (quality_function_norm_) {
      case NM_NORM_1:
      curr_grad_lag_x_asum_ = IpCq().curr_grad_lag_x()->Asum();
      curr_grad_lag_s_asum_ = IpCq().curr_grad_lag_s()->Asum();
      curr_c_asum_ = IpCq().curr_c()->Asum();
      curr_d_minus_s_asum_ = IpCq().curr_d_minus_s()->Asum();
      break;
      case NM_NORM_2_SQUARED:
      case NM_NORM_2:
      curr_grad_lag_x_nrm2_ = IpCq().curr_grad_lag_x()->Nrm2();
      curr_grad_lag_s_nrm2_ = IpCq().curr_grad_lag_s()->Nrm2();
      curr_c_nrm2_ = IpCq().curr_c()->Nrm2();
      curr_d_minus_s_nrm2_ = IpCq().curr_d_minus_s()->Nrm2();
      break;
      case NM_NORM_MAX:
      curr_grad_lag_x_amax_ = IpCq().curr_grad_lag_x()->Amax();
      curr_grad_lag_s_amax_ = IpCq().curr_grad_lag_s()->Amax();
      curr_c_amax_ = IpCq().curr_c()->Amax();
      curr_d_minus_s_amax_ = IpCq().curr_d_minus_s()->Amax();
      break;
      default:
      DBG_ASSERT(false && "Unknown value for quality_function_norm_");
    }
    IpData().TimingStats().Task5().End();

    // Some initializations
    if (!initialized_) {
      n_dual_ = IpData().curr()->x()->Dim() + IpData().curr()->s()->Dim();
      n_pri_ = IpData().curr()->y_c()->Dim() + IpData().curr()->y_d()->Dim();
      n_comp_ = IpData().curr()->z_L()->Dim() + IpData().curr()->z_U()->Dim() +
                IpData().curr()->v_L()->Dim() + IpData().curr()->v_U()->Dim();
      initialized_ = true;
    }

    count_qf_evals_ = 0;

    // Compute some quantities used for the quality function evaluations
    // (This way we try to avoid retrieving numbers from cache...

    curr_slack_x_L_ = IpCq().curr_slack_x_L();
    curr_slack_x_U_ = IpCq().curr_slack_x_U();
    curr_slack_s_L_ = IpCq().curr_slack_s_L();
    curr_slack_s_U_ = IpCq().curr_slack_s_U();

    curr_z_L_ = IpData().curr()->z_L();
    curr_z_U_ = IpData().curr()->z_U();
    curr_v_L_ = IpData().curr()->v_L();
    curr_v_U_ = IpData().curr()->v_U();

    IpData().TimingStats().Task5().Start();
    switch (quality_function_norm_) {
      case NM_NORM_1:
      curr_grad_lag_x_asum_ = IpCq().curr_grad_lag_x()->Asum();
      curr_grad_lag_s_asum_ = IpCq().curr_grad_lag_s()->Asum();
      curr_c_asum_ = IpCq().curr_c()->Asum();
      curr_d_minus_s_asum_ = IpCq().curr_d_minus_s()->Asum();
      break;
      case NM_NORM_2_SQUARED:
      case NM_NORM_2:
      curr_grad_lag_x_nrm2_ = IpCq().curr_grad_lag_x()->Nrm2();
      curr_grad_lag_s_nrm2_ = IpCq().curr_grad_lag_s()->Nrm2();
      curr_c_nrm2_ = IpCq().curr_c()->Nrm2();
      curr_d_minus_s_nrm2_ = IpCq().curr_d_minus_s()->Nrm2();
      break;
      case NM_NORM_MAX:
      curr_grad_lag_x_amax_ = IpCq().curr_grad_lag_x()->Amax();
      curr_grad_lag_s_amax_ = IpCq().curr_grad_lag_s()->Amax();
      curr_c_amax_ = IpCq().curr_c()->Amax();
      curr_d_minus_s_amax_ = IpCq().curr_d_minus_s()->Amax();
      break;
      default:
      DBG_ASSERT(false && "Unknown value for quality_function_norm_");
    }
    IpData().TimingStats().Task5().End();

    // We now compute the step for the slack variables.  This safes
    // time, because we then don't have to do this any more for each
    // evaluation of the quality function
    SmartPtr<Vector> step_aff_x_L = step_aff->z_L()->MakeNew();
    SmartPtr<Vector> step_aff_x_U = step_aff->z_U()->MakeNew();
    SmartPtr<Vector> step_aff_s_L = step_aff->v_L()->MakeNew();
    SmartPtr<Vector> step_aff_s_U = step_aff->v_U()->MakeNew();
    IpNLP().Px_L()->TransMultVector(1., *step_aff->x(), 0., *step_aff_x_L);
    IpNLP().Px_U()->TransMultVector(-1., *step_aff->x(), 0., *step_aff_x_U);
    IpNLP().Pd_L()->TransMultVector(1., *step_aff->s(), 0., *step_aff_s_L);
    IpNLP().Pd_U()->TransMultVector(-1., *step_aff->s(), 0., *step_aff_s_U);
    SmartPtr<Vector> step_cen_x_L = step_cen->z_L()->MakeNew();
    SmartPtr<Vector> step_cen_x_U = step_cen->z_U()->MakeNew();
    SmartPtr<Vector> step_cen_s_L = step_cen->v_L()->MakeNew();
    SmartPtr<Vector> step_cen_s_U = step_cen->v_U()->MakeNew();
    IpNLP().Px_L()->TransMultVector(1., *step_cen->x(), 0., *step_cen_x_L);
    IpNLP().Px_U()->TransMultVector(-1., *step_cen->x(), 0., *step_cen_x_U);
    IpNLP().Pd_L()->TransMultVector(1., *step_cen->s(), 0., *step_cen_s_L);
    IpNLP().Pd_U()->TransMultVector(-1., *step_cen->s(), 0., *step_cen_s_U);

    Number sigma;

    // First we determine whether we want to search for a value of
    // sigma larger or smaller than 1.  For this, we estimate the
    // slope of the quality function at sigma=1.
    Number qf_1 = CalculateQualityFunction(1.,
                                           *step_aff_x_L,
                                           *step_aff_x_U,
                                           *step_aff_s_L,
                                           *step_aff_s_U,
                                           *step_aff->y_c(),
                                           *step_aff->y_d(),
                                           *step_aff->z_L(),
                                           *step_aff->z_U(),
                                           *step_aff->v_L(),
                                           *step_aff->v_U(),
                                           *step_cen_x_L,
                                           *step_cen_x_U,
                                           *step_cen_s_L,
                                           *step_cen_s_U,
                                           *step_cen->y_c(),
                                           *step_cen->y_d(),
                                           *step_cen->z_L(),
                                           *step_cen->z_U(),
                                           *step_cen->v_L(),
                                           *step_cen->v_U());

    Number sigma_1minus = 1.-Max(1e-4, quality_function_section_sigma_tol_);
    Number qf_1minus = CalculateQualityFunction(sigma_1minus,
                       *step_aff_x_L,
                       *step_aff_x_U,
                       *step_aff_s_L,
                       *step_aff_s_U,
                       *step_aff->y_c(),
                       *step_aff->y_d(),
                       *step_aff->z_L(),
                       *step_aff->z_U(),
                       *step_aff->v_L(),
                       *step_aff->v_U(),
                       *step_cen_x_L,
                       *step_cen_x_U,
                       *step_cen_s_L,
                       *step_cen_s_U,
                       *step_cen->y_c(),
                       *step_cen->y_d(),
                       *step_cen->z_L(),
                       *step_cen->z_U(),
                       *step_cen->v_L(),
                       *step_cen->v_U());

    if (qf_1minus > qf_1) {
      // It seems that the quality function decreases for values
      // larger than sigma, so perform golden section search for sigma
      // > 1.
      Number sigma_up = Min(sigma_max_, mu_max/avrg_compl);
      Number sigma_lo = 1.;
      if (sigma_lo >= sigma_up) {
        sigma = sigma_up;
      }
      else {
        // ToDo maybe we should use different tolerances for sigma>1
        sigma = PerformGoldenSection(sigma_up, -100., sigma_lo, qf_1,
                                     quality_function_section_sigma_tol_,
                                     quality_function_section_qf_tol_,
                                     *step_aff_x_L,
                                     *step_aff_x_U,
                                     *step_aff_s_L,
                                     *step_aff_s_U,
                                     *step_aff->y_c(),
                                     *step_aff->y_d(),
                                     *step_aff->z_L(),
                                     *step_aff->z_U(),
                                     *step_aff->v_L(),
                                     *step_aff->v_U(),
                                     *step_cen_x_L,
                                     *step_cen_x_U,
                                     *step_cen_s_L,
                                     *step_cen_s_U,
                                     *step_cen->y_c(),
                                     *step_cen->y_d(),
                                     *step_cen->z_L(),
                                     *step_cen->z_U(),
                                     *step_cen->v_L(),
                                     *step_cen->v_U());
      }
    }
    else {
      // Search for sigma less than 1

      Number sigma_lo = Max(sigma_min_, mu_min/avrg_compl);
      Number sigma_up = Min(Max(sigma_lo, sigma_1minus), mu_max/avrg_compl);
      if (sigma_lo >= sigma_up) {
        // Skip the search, we are already at the minimum
        sigma = sigma_lo;
      }
      else {
        sigma = PerformGoldenSection(sigma_up, qf_1minus, sigma_lo, -100.,
                                     quality_function_section_sigma_tol_,
                                     quality_function_section_qf_tol_,
                                     *step_aff_x_L,
                                     *step_aff_x_U,
                                     *step_aff_s_L,
                                     *step_aff_s_U,
                                     *step_aff->y_c(),
                                     *step_aff->y_d(),
                                     *step_aff->z_L(),
                                     *step_aff->z_U(),
                                     *step_aff->v_L(),
                                     *step_aff->v_U(),
                                     *step_cen_x_L,
                                     *step_cen_x_U,
                                     *step_cen_s_L,
                                     *step_cen_s_U,
                                     *step_cen->y_c(),
                                     *step_cen->y_d(),
                                     *step_cen->z_L(),
                                     *step_cen->z_U(),
                                     *step_cen->v_L(),
                                     *step_cen->v_U());
      }
    }

    //#define tracequalityfunction
#ifdef tracequalityfunction
    char fname[100];
    sprintf(fname, "qf_values_%d.dat", IpData().iter_count());
    FILE* fid = fopen(fname, "w");

    Number sigma_1 = sigma_max_;
    Number sigma_2 = 1e-9/avrg_compl;
    Number sigma_trace = sigma_1;
    while(sigma_trace > sigma_2) {
      Number qf = CalculateQualityFunction(sigma_trace,
                                           *step_aff_x_L,
                                           *step_aff_x_U,
                                           *step_aff_s_L,
                                           *step_aff_s_U,
                                           *step_aff->y_c(),
                                           *step_aff->y_d(),
                                           *step_aff->z_L(),
                                           *step_aff->z_U(),
                                           *step_aff->v_L(),
                                           *step_aff->v_U(),
                                           *step_cen_x_L,
                                           *step_cen_x_U,
                                           *step_cen_s_L,
                                           *step_cen_s_U,
                                           *step_cen->y_c(),
                                           *step_cen->y_d(),
                                           *step_cen->z_L(),
                                           *step_cen->z_U(),
                                           *step_cen->v_L(),
                                           *step_cen->v_U());
      fprintf(fid, "%9.2e %25.16e\n", sigma_trace, qf);
      sigma_trace /= 1.1;
    }
    fclose(fid);
#endif

    // End timing of quality function search
    IpData().TimingStats().QualityFunctionSearch().End();

    Jnlst().Printf(J_DETAILED, J_BARRIER_UPDATE,
                   "Sigma = %e\n", sigma);
    Number mu = sigma*avrg_compl;

    // Store the affine search direction (in case it is needed in the
    // line search for a corrector step)
    IpData().set_delta_aff(step_aff);
    IpData().SetHaveAffineDeltas(true);

    // Now construct the overall search direction here
    SmartPtr<IteratesVector> step = IpData().curr()->MakeNewIteratesVector(true);
    step->AddTwoVectors(sigma, *step_cen, 1.0, *IpData().delta_aff(), 0.0);

    DBG_PRINT_VECTOR(2, "step", *step);
    IpData().set_delta(step);
    IpData().SetHaveDeltas(true);

    ///////////////////////////////////////////////////////////////////////////
    // Release memory for temporary vectors used in CalculateQualityFunction //
    ///////////////////////////////////////////////////////////////////////////

    tmp_step_x_L_ = NULL;
    tmp_step_x_U_ = NULL;
    tmp_step_s_L_ = NULL;
    tmp_step_s_U_ = NULL;
    tmp_step_z_L_ = NULL;
    tmp_step_z_U_ = NULL;
    tmp_step_v_L_ = NULL;
    tmp_step_v_U_ = NULL;

    tmp_slack_x_L_ = NULL;
    tmp_slack_x_U_ = NULL;
    tmp_slack_s_L_ = NULL;
    tmp_slack_s_U_ = NULL;
    tmp_z_L_ = NULL;
    tmp_z_U_ = NULL;
    tmp_v_L_ = NULL;
    tmp_v_U_ = NULL;

    curr_slack_x_L_ = NULL;
    curr_slack_x_U_ = NULL;
    curr_slack_s_L_ = NULL;
    curr_slack_s_U_ = NULL;

    // DELETEME
    char ssigma[40];
    sprintf(ssigma, " sigma=%8.2e", sigma);
    IpData().Append_info_string(ssigma);
    sprintf(ssigma, " qf=%d", count_qf_evals_);
    IpData().Append_info_string(ssigma);
    /*
    sprintf(ssigma, " xi=%8.2e ", IpCq().curr_centrality_measure());
    IpData().Append_info_string(ssigma);
    if (sigma>1.) {
      IpData().Append_info_string("LARGESIGMA");
    }
    */

    new_mu = mu;
    return true;
  }

  Number QualityFunctionMuOracle::CalculateQualityFunction
  (Number sigma,
   const Vector& step_aff_x_L,
   const Vector& step_aff_x_U,
   const Vector& step_aff_s_L,
   const Vector& step_aff_s_U,
   const Vector& step_aff_y_c,
   const Vector& step_aff_y_d,
   const Vector& step_aff_z_L,
   const Vector& step_aff_z_U,
   const Vector& step_aff_v_L,
   const Vector& step_aff_v_U,
   const Vector& step_cen_x_L,
   const Vector& step_cen_x_U,
   const Vector& step_cen_s_L,
   const Vector& step_cen_s_U,
   const Vector& step_cen_y_c,
   const Vector& step_cen_y_d,
   const Vector& step_cen_z_L,
   const Vector& step_cen_z_U,
   const Vector& step_cen_v_L,
   const Vector& step_cen_v_U
  )
  {
    DBG_START_METH("QualityFunctionMuOracle::CalculateQualityFunction",
                   dbg_verbosity);
    count_qf_evals_++;

    IpData().TimingStats().Task1().Start();
    tmp_step_x_L_->AddTwoVectors(1., step_aff_x_L, sigma, step_cen_x_L, 0.);
    tmp_step_x_U_->AddTwoVectors(1., step_aff_x_U, sigma, step_cen_x_U, 0.);
    tmp_step_s_L_->AddTwoVectors(1., step_aff_s_L, sigma, step_cen_s_L, 0.);
    tmp_step_s_U_->AddTwoVectors(1., step_aff_s_U, sigma, step_cen_s_U, 0.);
    tmp_step_z_L_->AddTwoVectors(1., step_aff_z_L, sigma, step_cen_z_L, 0.);
    tmp_step_z_U_->AddTwoVectors(1., step_aff_z_U, sigma, step_cen_z_U, 0.);
    tmp_step_v_L_->AddTwoVectors(1., step_aff_v_L, sigma, step_cen_v_L, 0.);
    tmp_step_v_U_->AddTwoVectors(1., step_aff_v_U, sigma, step_cen_v_U, 0.);
    IpData().TimingStats().Task1().End();

    // Compute the fraction-to-the-boundary step sizes
    IpData().TimingStats().Task2().Start();
    Number tau = IpData().curr_tau();
    Number alpha_primal = IpCq().uncached_slack_frac_to_the_bound(tau,
                          *tmp_step_x_L_,
                          *tmp_step_x_U_,
                          *tmp_step_s_L_,
                          *tmp_step_s_U_);

    Number alpha_dual = IpCq().uncached_dual_frac_to_the_bound(tau,
                        *tmp_step_z_L_,
                        *tmp_step_z_U_,
                        *tmp_step_v_L_,
                        *tmp_step_v_U_);
    IpData().TimingStats().Task2().End();

    Number xi = 0.; // centrality measure

    IpData().TimingStats().Task1().Start();
    tmp_slack_x_L_->AddTwoVectors(1., *curr_slack_x_L_,
                                  alpha_primal, *tmp_step_x_L_, 0.);
    tmp_slack_x_U_->AddTwoVectors(1., *curr_slack_x_U_,
                                  alpha_primal, *tmp_step_x_U_, 0.);
    tmp_slack_s_L_->AddTwoVectors(1., *curr_slack_s_L_,
                                  alpha_primal, *tmp_step_s_L_, 0.);
    tmp_slack_s_U_->AddTwoVectors(1., *curr_slack_s_U_,
                                  alpha_primal, *tmp_step_s_U_, 0.);

    tmp_z_L_->AddTwoVectors(1., *curr_z_L_,
                            alpha_dual, *tmp_step_z_L_, 0.);
    tmp_z_U_->AddTwoVectors(1., *curr_z_U_,
                            alpha_dual, *tmp_step_z_U_, 0.);
    tmp_v_L_->AddTwoVectors(1., *curr_v_L_,
                            alpha_dual, *tmp_step_v_L_, 0.);
    tmp_v_U_->AddTwoVectors(1., *curr_v_U_,
                            alpha_dual, *tmp_step_v_U_, 0.);
    IpData().TimingStats().Task1().End();

    IpData().TimingStats().Task3().Start();
    tmp_slack_x_L_->ElementWiseMultiply(*tmp_z_L_);
    tmp_slack_x_U_->ElementWiseMultiply(*tmp_z_U_);
    tmp_slack_s_L_->ElementWiseMultiply(*tmp_v_L_);
    tmp_slack_s_U_->ElementWiseMultiply(*tmp_v_U_);
    IpData().TimingStats().Task3().End();

    DBG_PRINT_VECTOR(2, "compl_x_L", *tmp_slack_x_L_);
    DBG_PRINT_VECTOR(2, "compl_x_U", *tmp_slack_x_U_);
    DBG_PRINT_VECTOR(2, "compl_s_L", *tmp_slack_s_L_);
    DBG_PRINT_VECTOR(2, "compl_s_U", *tmp_slack_s_U_);

    Number dual_inf=-1.;
    Number primal_inf=-1.;
    Number compl_inf=-1.;

    IpData().TimingStats().Task5().Start();
    switch (quality_function_norm_) {
      case NM_NORM_1:
      dual_inf = (1.-alpha_dual)*(curr_grad_lag_x_asum_ +
                                  curr_grad_lag_s_asum_);

      primal_inf = (1.-alpha_primal)*(curr_c_asum_ +
                                      curr_d_minus_s_asum_);

      compl_inf = tmp_slack_x_L_->Asum() + tmp_slack_x_U_->Asum() +
                  tmp_slack_s_L_->Asum() + tmp_slack_s_U_->Asum();

      dual_inf /= n_dual_;
      if (n_pri_>0) {
        primal_inf /= n_pri_;
      }
      DBG_ASSERT(n_comp_>0);
      compl_inf /= n_comp_;
      break;
      case NM_NORM_2_SQUARED:
      dual_inf =
        pow(1.-alpha_dual, 2)*(pow(curr_grad_lag_x_nrm2_, 2) +
                               pow(curr_grad_lag_s_nrm2_, 2));
      primal_inf =
        pow(1.-alpha_primal, 2)*(pow(curr_c_nrm2_, 2) +
                                 pow(curr_d_minus_s_nrm2_, 2));
      compl_inf =
        pow(tmp_slack_x_L_->Nrm2(), 2) + pow(tmp_slack_x_U_->Nrm2(), 2) +
        pow(tmp_slack_s_L_->Nrm2(), 2) + pow(tmp_slack_s_U_->Nrm2(), 2);

      dual_inf /= n_dual_;
      if (n_pri_>0) {
        primal_inf /= n_pri_;
      }
      DBG_ASSERT(n_comp_>0);
      compl_inf /= n_comp_;
      break;
      case NM_NORM_MAX:
      dual_inf =
        (1.-alpha_dual)*Max(curr_grad_lag_x_amax_,
                            curr_grad_lag_s_amax_);
      primal_inf =
        (1.-alpha_primal)*Max(curr_c_amax_,
                              curr_d_minus_s_amax_);
      compl_inf =
        Max(tmp_slack_x_L_->Amax(), tmp_slack_x_U_->Amax(),
            tmp_slack_s_L_->Amax(), tmp_slack_s_U_->Amax());
      break;
      case NM_NORM_2:
      dual_inf =
        (1.-alpha_dual)*sqrt(pow(curr_grad_lag_x_nrm2_, 2) +
                             pow(curr_grad_lag_s_nrm2_, 2));
      primal_inf =
        (1.-alpha_primal)*sqrt(pow(curr_c_nrm2_, 2) +
                               pow(curr_d_minus_s_nrm2_, 2));
      compl_inf =
        sqrt(pow(tmp_slack_x_L_->Nrm2(), 2) + pow(tmp_slack_x_U_->Nrm2(), 2) +
             pow(tmp_slack_s_L_->Nrm2(), 2) + pow(tmp_slack_s_U_->Nrm2(), 2));

      dual_inf /= sqrt((Number)n_dual_);
      if (n_pri_>0) {
        primal_inf /= sqrt((Number)n_pri_);
      }
      DBG_ASSERT(n_comp_>0);
      compl_inf /= sqrt((Number)n_comp_);
      break;
      default:
      DBG_ASSERT(false && "Unknown value for quality_function_norm_");
    }
    IpData().TimingStats().Task5().End();

    Number quality_function = dual_inf + primal_inf + compl_inf;

    if (quality_function_centrality_!=CEN_NONE) {
      IpData().TimingStats().Task4().Start();
      xi = IpCq().CalcCentralityMeasure(*tmp_slack_x_L_, *tmp_slack_x_U_,
                                        *tmp_slack_s_L_, *tmp_slack_s_U_);
      IpData().TimingStats().Task4().End();
    }
    switch (quality_function_centrality_) {
      case CEN_NONE:
      //Nothing
      break;
      case CEN_LOG:
      quality_function -= compl_inf*log(xi);
      break;
      case CEN_RECIPROCAL:
      quality_function += compl_inf/xi;
      case CEN_CUBED_RECIPROCAL:
      quality_function += compl_inf/pow(xi,3);
      break;
      default:
      DBG_ASSERT(false && "Unknown value for quality_function_centrality_");
    }

    switch (quality_function_balancing_term_) {
      case BT_NONE:
      //Nothing
      break;
      case BT_CUBIC:
      quality_function += pow(Max(0., Max(dual_inf,primal_inf)-compl_inf),3);
      break;
      default:
      DBG_ASSERT(false && "Unknown value for quality_function_balancing term_");
    }

    Jnlst().Printf(J_MOREDETAILED, J_BARRIER_UPDATE,
                   "sigma = %8.2e d_inf = %18.12e p_inf = %18.12e cmpl = %18.12e q = %18.12e a_pri = %8.2e a_dual = %8.2e xi = %8.2e\n", sigma, dual_inf, primal_inf, compl_inf, quality_function, alpha_primal, alpha_dual, xi);


    return quality_function;
    //return compl_inf;
  }

  Number
  QualityFunctionMuOracle::PerformGoldenSection
  (Number sigma_up_in,
   Number q_up,
   Number sigma_lo_in,
   Number q_lo,
   Number sigma_tol,
   Number qf_tol,
   const Vector& step_aff_x_L,
   const Vector& step_aff_x_U,
   const Vector& step_aff_s_L,
   const Vector& step_aff_s_U,
   const Vector& step_aff_y_c,
   const Vector& step_aff_y_d,
   const Vector& step_aff_z_L,
   const Vector& step_aff_z_U,
   const Vector& step_aff_v_L,
   const Vector& step_aff_v_U,
   const Vector& step_cen_x_L,
   const Vector& step_cen_x_U,
   const Vector& step_cen_s_L,
   const Vector& step_cen_s_U,
   const Vector& step_cen_y_c,
   const Vector& step_cen_y_d,
   const Vector& step_cen_z_L,
   const Vector& step_cen_z_U,
   const Vector& step_cen_v_L,
   const Vector& step_cen_v_U
  )
  {
    Number sigma_up = ScaleSigma(sigma_up_in);
    Number sigma_lo = ScaleSigma(sigma_lo_in);

    Number sigma;
    Number gfac = (3.-sqrt(5.))/2.;
    Number sigma_mid1 = sigma_lo + gfac*(sigma_up-sigma_lo);
    Number sigma_mid2 = sigma_lo + (1.-gfac)*(sigma_up-sigma_lo);

    Number qmid1 = CalculateQualityFunction(UnscaleSigma(sigma_mid1),
                                            step_aff_x_L,
                                            step_aff_x_U,
                                            step_aff_s_L,
                                            step_aff_s_U,
                                            step_aff_y_c,
                                            step_aff_y_d,
                                            step_aff_z_L,
                                            step_aff_z_U,
                                            step_aff_v_L,
                                            step_aff_v_U,
                                            step_cen_x_L,
                                            step_cen_x_U,
                                            step_cen_s_L,
                                            step_cen_s_U,
                                            step_cen_y_c,
                                            step_cen_y_d,
                                            step_cen_z_L,
                                            step_cen_z_U,
                                            step_cen_v_L,
                                            step_cen_v_U);
    Number qmid2 = CalculateQualityFunction(UnscaleSigma(sigma_mid2),
                                            step_aff_x_L,
                                            step_aff_x_U,
                                            step_aff_s_L,
                                            step_aff_s_U,
                                            step_aff_y_c,
                                            step_aff_y_d,
                                            step_aff_z_L,
                                            step_aff_z_U,
                                            step_aff_v_L,
                                            step_aff_v_U,
                                            step_cen_x_L,
                                            step_cen_x_U,
                                            step_cen_s_L,
                                            step_cen_s_U,
                                            step_cen_y_c,
                                            step_cen_y_d,
                                            step_cen_z_L,
                                            step_cen_z_U,
                                            step_cen_v_L,
                                            step_cen_v_U);

    Index nsections = 0;
    while ((sigma_up-sigma_lo)>=sigma_tol*sigma_up &&
           //while ((sigma_up-sigma_lo)>=sigma_tol &&  // Note we are using the non-relative criterion here for sigma
           (1.-Min(q_lo, q_up, qmid1, qmid2)/Max(q_lo, q_up, qmid1, qmid2))>=qf_tol &&
           nsections<quality_function_max_section_steps_) {
      nsections++;
      //printf("sigma_lo=%e sigma_mid1=%e sigma_mid2=%e sigma_up=%e\n",sigma_lo, sigma_mid1, sigma_mid2, sigma_up);
      if (qmid1 > qmid2) {
        sigma_lo = sigma_mid1;
        q_lo = qmid1;
        sigma_mid1 = sigma_mid2;
        qmid1 = qmid2;
        sigma_mid2 = sigma_lo + (1.-gfac)*(sigma_up-sigma_lo);
        qmid2 = CalculateQualityFunction(UnscaleSigma(sigma_mid2),
                                         step_aff_x_L,
                                         step_aff_x_U,
                                         step_aff_s_L,
                                         step_aff_s_U,
                                         step_aff_y_c,
                                         step_aff_y_d,
                                         step_aff_z_L,
                                         step_aff_z_U,
                                         step_aff_v_L,
                                         step_aff_v_U,
                                         step_cen_x_L,
                                         step_cen_x_U,
                                         step_cen_s_L,
                                         step_cen_s_U,
                                         step_cen_y_c,
                                         step_cen_y_d,
                                         step_cen_z_L,
                                         step_cen_z_U,
                                         step_cen_v_L,
                                         step_cen_v_U);
      }
      else {
        sigma_up = sigma_mid2;
        q_up = qmid2;
        sigma_mid2 = sigma_mid1;
        qmid2 = qmid1;
        sigma_mid1 = sigma_lo + gfac*(sigma_up-sigma_lo);
        qmid1 = CalculateQualityFunction(UnscaleSigma(sigma_mid1),
                                         step_aff_x_L,
                                         step_aff_x_U,
                                         step_aff_s_L,
                                         step_aff_s_U,
                                         step_aff_y_c,
                                         step_aff_y_d,
                                         step_aff_z_L,
                                         step_aff_z_U,
                                         step_aff_v_L,
                                         step_aff_v_U,
                                         step_cen_x_L,
                                         step_cen_x_U,
                                         step_cen_s_L,
                                         step_cen_s_U,
                                         step_cen_y_c,
                                         step_cen_y_d,
                                         step_cen_z_L,
                                         step_cen_z_U,
                                         step_cen_v_L,
                                         step_cen_v_U);
      }
    }

    if ((sigma_up-sigma_lo)>=sigma_tol*sigma_up &&
        (1.-Min(q_lo, q_up, qmid1, qmid2)/Max(q_lo, q_up, qmid1, qmid2))<qf_tol) {
      // The qf tolerance make it stop
      IpData().Append_info_string("qf_tol ");
      Number qf_min = Min(q_lo, q_up, qmid1, qmid2);
      DBG_ASSERT(qf_min>-100.);
      if (qf_min == q_lo) {
        sigma = sigma_lo;
      }
      else if (qf_min == qmid1) {
        sigma = sigma_mid1;
      }
      else if (qf_min == qmid2) {
        sigma = sigma_mid2;
      }
      else {
        sigma = sigma_up;
      }
    }
    else {
      Number q;
      if (qmid1 < qmid2) {
        sigma = sigma_mid1;
        q = qmid1;
      }
      else {
        sigma = sigma_mid2;
        q = qmid2;
      }
      if (sigma_up == ScaleSigma(sigma_up_in)) {
        Number qtmp;
        if (q_up<0.) {
          qtmp = CalculateQualityFunction(UnscaleSigma(sigma_up),
                                          step_aff_x_L,
                                          step_aff_x_U,
                                          step_aff_s_L,
                                          step_aff_s_U,
                                          step_aff_y_c,
                                          step_aff_y_d,
                                          step_aff_z_L,
                                          step_aff_z_U,
                                          step_aff_v_L,
                                          step_aff_v_U,
                                          step_cen_x_L,
                                          step_cen_x_U,
                                          step_cen_s_L,
                                          step_cen_s_U,
                                          step_cen_y_c,
                                          step_cen_y_d,
                                          step_cen_z_L,
                                          step_cen_z_U,
                                          step_cen_v_L,
                                          step_cen_v_U);
        }
        else {
          qtmp = q_up;
        }
        if (qtmp < q) {
          sigma = sigma_up;
          q = qtmp;
        }
      }
      else if (sigma_lo == ScaleSigma(sigma_lo_in)) {
        Number qtmp;
        if (q_lo<0.) {
          qtmp = CalculateQualityFunction(UnscaleSigma(sigma_lo),
                                          step_aff_x_L,
                                          step_aff_x_U,
                                          step_aff_s_L,
                                          step_aff_s_U,
                                          step_aff_y_c,
                                          step_aff_y_d,
                                          step_aff_z_L,
                                          step_aff_z_U,
                                          step_aff_v_L,
                                          step_aff_v_U,
                                          step_cen_x_L,
                                          step_cen_x_U,
                                          step_cen_s_L,
                                          step_cen_s_U,
                                          step_cen_y_c,
                                          step_cen_y_d,
                                          step_cen_z_L,
                                          step_cen_z_U,
                                          step_cen_v_L,
                                          step_cen_v_U);
        }
        else {
          qtmp = q_lo;
        }
        if (qtmp < q) {
          sigma = sigma_lo;
          q = qtmp;
        }
      }
    }

    return UnscaleSigma(sigma);
  }

  /*
  Number QualityFunctionMuOracle::ScaleSigma(Number sigma) {return log(sigma);}
  Number QualityFunctionMuOracle::UnscaleSigma(Number scaled_sigma) {return exp(scaled_sigma);}
  */
  Number QualityFunctionMuOracle::ScaleSigma(Number sigma)
  {
    return sigma;
  }
  Number QualityFunctionMuOracle::UnscaleSigma(Number scaled_sigma)
  {
    return scaled_sigma;
  }

  /* AW: Tried search in the log space, but that was even worse than
     search in unscaled space */
  /*
  Number
  QualityFunctionMuOracle::PerformGoldenSectionLog
  (Number sigma_up,
   Number sigma_lo,
   Number tol,
   const Vector& step_aff_x_L,
   const Vector& step_aff_x_U,
   const Vector& step_aff_s_L,
   const Vector& step_aff_s_U,
   const Vector& step_aff_y_c,
   const Vector& step_aff_y_d,
   const Vector& step_aff_z_L,
   const Vector& step_aff_z_U,
   const Vector& step_aff_v_L,
   const Vector& step_aff_v_U,
   const Vector& step_cen_x_L,
   const Vector& step_cen_x_U,
   const Vector& step_cen_s_L,
   const Vector& step_cen_s_U,
   const Vector& step_cen_y_c,
   const Vector& step_cen_y_d,
   const Vector& step_cen_z_L,
   const Vector& step_cen_z_U,
   const Vector& step_cen_v_L,
   const Vector& step_cen_v_U
  )
  {
    Number log_sigma;
    Number log_sigma_up = log(sigma_up);
    Number log_sigma_lo = log(sigma_lo);

    Number log_sigma_up_in = log_sigma_up;
    Number log_sigma_lo_in = log_sigma_lo;
    Number gfac = (3.-sqrt(5.))/2.;
    Number log_sigma_mid1 = log_sigma_lo + gfac*(log_sigma_up-log_sigma_lo);
    Number log_sigma_mid2 = log_sigma_lo + (1.-gfac)*(log_sigma_up-log_sigma_lo);

    Number qmid1 = CalculateQualityFunction(exp(log_sigma_mid1),
                                            step_aff_x_L,
                                            step_aff_x_U,
                                            step_aff_s_L,
                                            step_aff_s_U,
                                            step_aff_y_c,
                                            step_aff_y_d,
                                            step_aff_z_L,
                                            step_aff_z_U,
                                            step_aff_v_L,
                                            step_aff_v_U,
                                            step_cen_x_L,
                                            step_cen_x_U,
                                            step_cen_s_L,
                                            step_cen_s_U,
                                            step_cen_y_c,
                                            step_cen_y_d,
                                            step_cen_z_L,
                                            step_cen_z_U,
                                            step_cen_v_L,
                                            step_cen_v_U);
    Number qmid2 = CalculateQualityFunction(exp(log_sigma_mid2),
                                            step_aff_x_L,
                                            step_aff_x_U,
                                            step_aff_s_L,
                                            step_aff_s_U,
                                            step_aff_y_c,
                                            step_aff_y_d,
                                            step_aff_z_L,
                                            step_aff_z_U,
                                            step_aff_v_L,
                                            step_aff_v_U,
                                            step_cen_x_L,
                                            step_cen_x_U,
                                            step_cen_s_L,
                                            step_cen_s_U,
                                            step_cen_y_c,
                                            step_cen_y_d,
                                            step_cen_z_L,
                                            step_cen_z_U,
                                            step_cen_v_L,
                                            step_cen_v_U);

    Index nsections = 0;
    while ((log_sigma_up-log_sigma_lo)>=tol*log_sigma_up && nsections<quality_function_max_section_steps_) {
      nsections++;
      if (qmid1 > qmid2) {
        log_sigma_lo = log_sigma_mid1;
        log_sigma_mid1 = log_sigma_mid2;
        qmid1 = qmid2;
        log_sigma_mid2 = log_sigma_lo + (1.-gfac)*(log_sigma_up-log_sigma_lo);
        qmid2 = CalculateQualityFunction(exp(log_sigma_mid2),
                                         step_aff_x_L,
                                         step_aff_x_U,
                                         step_aff_s_L,
                                         step_aff_s_U,
                                         step_aff_y_c,
                                         step_aff_y_d,
                                         step_aff_z_L,
                                         step_aff_z_U,
                                         step_aff_v_L,
                                         step_aff_v_U,
                                         step_cen_x_L,
                                         step_cen_x_U,
                                         step_cen_s_L,
                                         step_cen_s_U,
                                         step_cen_y_c,
                                         step_cen_y_d,
                                         step_cen_z_L,
                                         step_cen_z_U,
                                         step_cen_v_L,
                                         step_cen_v_U);
      }
      else {
        log_sigma_up = log_sigma_mid2;
        log_sigma_mid2 = log_sigma_mid1;
        qmid2 = qmid1;
        log_sigma_mid1 = log_sigma_lo + gfac*(log_sigma_up-log_sigma_lo);
        qmid1 = CalculateQualityFunction(exp(log_sigma_mid1),
                                         step_aff_x_L,
                                         step_aff_x_U,
                                         step_aff_s_L,
                                         step_aff_s_U,
                                         step_aff_y_c,
                                         step_aff_y_d,
                                         step_aff_z_L,
                                         step_aff_z_U,
                                         step_aff_v_L,
                                         step_aff_v_U,
                                         step_cen_x_L,
                                         step_cen_x_U,
                                         step_cen_s_L,
                                         step_cen_s_U,
                                         step_cen_y_c,
                                         step_cen_y_d,
                                         step_cen_z_L,
                                         step_cen_z_U,
                                         step_cen_v_L,
                                         step_cen_v_U);
      }
    }

    Number q;
    if (qmid1 < qmid2) {
      log_sigma = log_sigma_mid1;
      q = qmid1;
    }
    else {
      log_sigma = log_sigma_mid2;
      q = qmid2;
    }
    if (log_sigma_up == log_sigma_up_in) {
      Number qtmp = CalculateQualityFunction(exp(log_sigma_up),
                                             step_aff_x_L,
                                             step_aff_x_U,
                                             step_aff_s_L,
                                             step_aff_s_U,
                                             step_aff_y_c,
                                             step_aff_y_d,
                                             step_aff_z_L,
                                             step_aff_z_U,
                                             step_aff_v_L,
                                             step_aff_v_U,
                                             step_cen_x_L,
                                             step_cen_x_U,
                                             step_cen_s_L,
                                             step_cen_s_U,
                                             step_cen_y_c,
                                             step_cen_y_d,
                                             step_cen_z_L,
                                             step_cen_z_U,
                                             step_cen_v_L,
                                             step_cen_v_U);
      if (qtmp < q) {
        log_sigma = log_sigma_up;
        q = qtmp;
      }
    }
    else if (log_sigma_lo == log_sigma_lo_in) {
      Number qtmp = CalculateQualityFunction(exp(log_sigma_lo),
                                             step_aff_x_L,
                                             step_aff_x_U,
                                             step_aff_s_L,
                                             step_aff_s_U,
                                             step_aff_y_c,
                                             step_aff_y_d,
                                             step_aff_z_L,
                                             step_aff_z_U,
                                             step_aff_v_L,
                                             step_aff_v_U,
                                             step_cen_x_L,
                                             step_cen_x_U,
                                             step_cen_s_L,
                                             step_cen_s_U,
                                             step_cen_y_c,
                                             step_cen_y_d,
                                             step_cen_z_L,
                                             step_cen_z_U,
                                             step_cen_v_L,
                                             step_cen_v_U);
      if (qtmp < q) {
        log_sigma = log_sigma_lo;
        q = qtmp;
      }
    }

    return exp(log_sigma);
  }
  */


} // namespace Ipopt
