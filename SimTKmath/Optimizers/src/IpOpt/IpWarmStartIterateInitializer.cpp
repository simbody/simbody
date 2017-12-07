// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpWarmStartIterateInitializer.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter              IBM    2005-04-01

#include "IpWarmStartIterateInitializer.hpp"
#include "IpDefaultIterateInitializer.hpp"

// ToDo make independent of DenseVector
#include "IpDenseVector.hpp"

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

  WarmStartIterateInitializer::WarmStartIterateInitializer()
      :
      IterateInitializer()
  {}

  void WarmStartIterateInitializer::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "warm_start_bound_push",
      "same as bound_push for the regular initializer.",
      0.0, true, 1e-3);
    roptions->AddBoundedNumberOption(
      "warm_start_bound_frac",
      "same as bound_frac for the regular initializer.",
      0.0, true, 0.5, false, 1e-3);
    roptions->AddLowerBoundedNumberOption(
      "warm_start_mult_bound_push",
      "same as mult_bound_push for the regular initializer.",
      0.0, true, 1e-3);
    roptions->AddNumberOption(
      "warm_start_mult_init_max",
      "Maximum initial value for the equality multipliers.",
      1e6);
    roptions->AddNumberOption(
      "warm_start_target_mu",
      "Unsupported!",
      0e-3);
    roptions->AddStringOption2(
      "warm_start_entire_iterate",
      "Tells algorithm whether to use the GetWarmStartIterate method in the NLP.",
      "no",
      "no", "call GetStartingPoint in the NLP",
      "yes", "call GetWarmStartIterate in the NLP",
      "");
  }

  bool WarmStartIterateInitializer::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("warm_start_bound_push",
                            warm_start_bound_push_, prefix);
    options.GetNumericValue("warm_start_bound_frac",
                            warm_start_bound_frac_, prefix);
    options.GetNumericValue("warm_start_mult_bound_push",
                            warm_start_mult_bound_push_, prefix);
    options.GetNumericValue("warm_start_mult_init_max",
                            warm_start_mult_init_max_, prefix);
    options.GetNumericValue("warm_start_target_mu",
                            warm_start_target_mu_, prefix);
    options.GetBoolValue("warm_start_entire_iterate",
                         warm_start_entire_iterate_, prefix);

    return true;
  }

  bool WarmStartIterateInitializer::SetInitialIterates()
  {
    DBG_START_METH("WarmStartIterateInitializer::SetInitialIterates",
                   dbg_verbosity);

    // Get the starting values provided by the NLP and store them
    // in the ip_data current fields.

    SmartPtr<IteratesVector> init_vec;
    bool have_iterate = false;

    if (warm_start_entire_iterate_) {
      IpData().InitializeDataStructures(IpNLP(), false, false, false,
                                        false, false);

      init_vec = IpData().curr()->MakeNewIteratesVector(true);

      have_iterate = IpNLP().GetWarmStartIterate(*init_vec);

      if (!have_iterate) {
        Jnlst().Printf(J_DETAILED, J_LINE_SEARCH,
                       "Tried to obtain entire warm start iterate from NLP, but it returned false.\n");
        IpData().Append_info_string("NW");
      }

      // Make sure given bounds are respected
      if (have_iterate && warm_start_mult_init_max_>0.) {
        SmartPtr<Vector> y_c = init_vec->create_new_y_c_copy();
        SmartPtr<Vector> tmp = y_c->MakeNew();
        tmp->Set(warm_start_mult_init_max_);
        y_c->ElementWiseMin(*tmp);
        tmp->Set(-warm_start_mult_init_max_);
        y_c->ElementWiseMax(*tmp);

        SmartPtr<Vector> y_d = init_vec->create_new_y_d_copy();
        tmp = y_d->MakeNew();
        tmp->Set(warm_start_mult_init_max_);
        y_d->ElementWiseMin(*tmp);
        tmp->Set(-warm_start_mult_init_max_);
        y_d->ElementWiseMax(*tmp);

        SmartPtr<Vector> z_L = init_vec->create_new_z_L_copy();
        tmp = z_L->MakeNew();
        tmp->Set(warm_start_mult_init_max_);
        z_L->ElementWiseMin(*tmp);

        SmartPtr<Vector> z_U = init_vec->create_new_z_U_copy();
        tmp = z_U->MakeNew();
        tmp->Set(warm_start_mult_init_max_);
        z_U->ElementWiseMin(*tmp);

        SmartPtr<Vector> v_L = init_vec->create_new_v_L_copy();
        tmp = v_L->MakeNew();
        tmp->Set(warm_start_mult_init_max_);
        v_L->ElementWiseMin(*tmp);

        SmartPtr<Vector> v_U = init_vec->create_new_v_U_copy();
        tmp = v_U->MakeNew();
        tmp->Set(warm_start_mult_init_max_);
        v_U->ElementWiseMin(*tmp);
      }
    }

    if (!have_iterate) {

      /////////////////////////////////////////////////////////////////////
      //                   Initialize primal variables                   //
      /////////////////////////////////////////////////////////////////////

      // Get the initial values for x, y_c, y_d, z_L, z_U,
      IpData().InitializeDataStructures(IpNLP(), true, true, true, true, true);

      IpData().curr()->x()->Print(Jnlst(), J_VECTOR, J_INITIALIZATION,
                                  "user-provided x");
      IpData().curr()->y_c()->Print(Jnlst(), J_VECTOR, J_INITIALIZATION,
                                    "user-provided y_c");
      IpData().curr()->y_d()->Print(Jnlst(), J_VECTOR, J_INITIALIZATION,
                                    "user-provided y_d");
      IpData().curr()->z_L()->Print(Jnlst(), J_VECTOR, J_INITIALIZATION,
                                    "user-provided z_L");
      IpData().curr()->z_U()->Print(Jnlst(), J_VECTOR, J_INITIALIZATION,
                                    "user-provided z_U");
      if (Jnlst().ProduceOutput(J_MOREVECTOR, J_INITIALIZATION)) {
        IpCq().curr_d()->Print(Jnlst(), J_MOREVECTOR, J_INITIALIZATION,
                               "d at user-provided x");
      }

      SmartPtr<Vector> tmp;

      init_vec = IpData().curr()->MakeNewContainer();

      // If requested, make sure that the multipliers are not too large
      if (warm_start_mult_init_max_>0.) {
        SmartPtr<Vector> y_c = init_vec->create_new_y_c_copy();
        tmp = y_c->MakeNew();
        tmp->Set(warm_start_mult_init_max_);
        y_c->ElementWiseMin(*tmp);
        tmp->Set(-warm_start_mult_init_max_);
        y_c->ElementWiseMax(*tmp);

        SmartPtr<Vector> y_d = init_vec->create_new_y_d_copy();
        tmp = y_d->MakeNew();
        tmp->Set(warm_start_mult_init_max_);
        y_d->ElementWiseMin(*tmp);
        tmp->Set(-warm_start_mult_init_max_);
        y_d->ElementWiseMax(*tmp);

        SmartPtr<Vector> z_L = init_vec->create_new_z_L_copy();
        tmp = z_L->MakeNew();
        tmp->Set(warm_start_mult_init_max_);
        z_L->ElementWiseMin(*tmp);

        SmartPtr<Vector> z_U = init_vec->create_new_z_U_copy();
        tmp = z_U->MakeNew();
        tmp->Set(warm_start_mult_init_max_);
        z_U->ElementWiseMin(*tmp);
      }

      // Get the initial values for v_L and v_U out of y_d
      SmartPtr<Vector> v_L = init_vec->create_new_v_L();
      IpNLP().Pd_L()->TransMultVector(-1., *init_vec->y_d(), 0., *v_L);
      tmp = v_L->MakeNew();
      tmp->Set(warm_start_mult_bound_push_);
      v_L->ElementWiseMax(*tmp);

      SmartPtr<Vector> v_U = init_vec->create_new_v_U();
      IpNLP().Pd_U()->TransMultVector(1., *init_vec->y_d(), 0., *v_U);
      tmp = v_U->MakeNew();
      tmp->Set(warm_start_mult_bound_push_);
      v_U->ElementWiseMax(*tmp);

      // Initialize slack variables
      init_vec->Set_s(*IpCq().curr_d());
    }

    // Make the corrected values current (and initialize s)
    IpData().set_trial(init_vec);
    IpData().AcceptTrialPoint();

    // Now apply the target mu heuristic if required
    if (warm_start_target_mu_>0.) {
      SmartPtr<const Vector> new_x;
      SmartPtr<const Vector> new_z_L;

      SmartPtr<const IteratesVector> curr = IpData().curr();
      process_target_mu(1., *curr->x(), *IpCq().curr_slack_x_L(),
                        *curr->z_L(), *IpNLP().Px_L(),
                        new_x, new_z_L);
      SmartPtr<const Vector> new_s;
      SmartPtr<const Vector> new_v_L;
      process_target_mu(1., *curr->s(), *IpCq().curr_slack_s_L(),
                        *curr->v_L(), *IpNLP().Pd_L(),
                        new_s, new_v_L);

      // Set the trial pointers to new_x and new_s. The process_target_mu
      // methods below create new vectors in new_x and new_s and do not alter
      // the existing ones.
      init_vec->Set_x(*new_x);
      init_vec->Set_s(*new_s);
      IpData().set_trial(init_vec);

      SmartPtr<const Vector> new_z_U;
      process_target_mu(-1., *IpData().trial()->x(), *IpCq().trial_slack_x_U(),
                        *IpData().curr()->z_U(), *IpNLP().Px_U(),
                        new_x, new_z_U);
      SmartPtr<const Vector> new_v_U;
      process_target_mu(-1., *IpData().trial()->s(), *IpCq().trial_slack_s_U(),
                        *IpData().curr()->v_U(), *IpNLP().Pd_U(),
                        new_s, new_v_U);

      // Now submit the full modified point
      init_vec->Set_x(*new_x);
      init_vec->Set_s(*new_s);
      // y_c and y_d currently contain a copy of curr()->y_c...
      // we set them back to the actual pointer to reuse the tags
      init_vec->Set_y_c(*IpData().curr()->y_c());
      init_vec->Set_y_d(*IpData().curr()->y_d());
      init_vec->Set_z_L(*new_z_L);
      init_vec->Set_z_U(*new_z_U);
      init_vec->Set_v_L(*new_v_L);
      init_vec->Set_v_U(*new_v_U);
      IpData().set_trial(init_vec);
      IpData().AcceptTrialPoint();

      // We need to call this to make sure that we don't get an error
      // message because at some point a slack became too small
      IpCq().ResetAdjustedTrialSlacks();
    }

    SmartPtr<const Vector> new_x;
    SmartPtr<const Vector> new_s;
    // Push the primal x variables
    DefaultIterateInitializer::push_variables(Jnlst(),
        warm_start_bound_push_,
        warm_start_bound_frac_,
        "x",
        *IpData().curr()->x(),
        new_x,
        *IpNLP().x_L(),
        *IpNLP().x_U(),
        *IpNLP().Px_L(),
        *IpNLP().Px_U());

    // Push the primal s variables
    DefaultIterateInitializer::push_variables(Jnlst(),
        warm_start_bound_push_,
        warm_start_bound_frac_,
        "s",
        *IpData().curr()->s(),
        new_s,
        *IpNLP().d_L(),
        *IpNLP().d_U(),
        *IpNLP().Pd_L(),
        *IpNLP().Pd_U());

    // Push the multipliers
    SmartPtr<Vector> new_z_L = IpData().curr()->z_L()->MakeNewCopy();
    SmartPtr<Vector> tmp = IpData().curr()->z_L()->MakeNew();
    tmp->Set(warm_start_mult_bound_push_);
    new_z_L->ElementWiseMax(*tmp);

    SmartPtr<Vector> new_z_U = IpData().curr()->z_U()->MakeNewCopy();
    tmp = IpData().curr()->z_U()->MakeNew();
    tmp->Set(warm_start_mult_bound_push_);
    new_z_U->ElementWiseMax(*tmp);

    SmartPtr<Vector> new_v_L = IpData().curr()->v_L()->MakeNewCopy();
    tmp = IpData().curr()->v_L()->MakeNew();
    tmp->Set(warm_start_mult_bound_push_);
    new_v_L->ElementWiseMax(*tmp);

    SmartPtr<Vector> new_v_U = IpData().curr()->v_U()->MakeNewCopy();
    tmp = IpData().curr()->v_U()->MakeNew();
    tmp->Set(warm_start_mult_bound_push_);
    new_v_U->ElementWiseMax(*tmp);

    // Make sure the new variables are current
    init_vec = IpData().curr()->MakeNewContainer();
    init_vec->Set_x(*new_x);
    init_vec->Set_s(*new_s);
    init_vec->Set_z_L(*new_z_L);
    init_vec->Set_z_U(*new_z_U);
    init_vec->Set_v_L(*new_v_L);
    init_vec->Set_v_U(*new_v_U);
    IpData().set_trial(init_vec);
    IpData().AcceptTrialPoint();

    IpData().curr()->x()->Print(Jnlst(), J_VECTOR, J_INITIALIZATION,
                                "initial x");
    IpData().curr()->s()->Print(Jnlst(), J_VECTOR, J_INITIALIZATION,
                                "initial s");
    IpData().curr()->y_c()->Print(Jnlst(), J_VECTOR, J_INITIALIZATION,
                                  "initial y_c");
    IpData().curr()->y_d()->Print(Jnlst(), J_VECTOR, J_INITIALIZATION,
                                  "initial y_d");
    IpData().curr()->z_L()->Print(Jnlst(), J_VECTOR, J_INITIALIZATION,
                                  "initial z_L");
    IpData().curr()->z_U()->Print(Jnlst(), J_VECTOR, J_INITIALIZATION,
                                  "initial z_U");
    IpData().curr()->v_L()->Print(Jnlst(), J_VECTOR, J_INITIALIZATION,
                                  "initial v_L");
    IpData().curr()->v_U()->Print(Jnlst(), J_VECTOR, J_INITIALIZATION,
                                  "initial v_U");
    if (Jnlst().ProduceOutput(J_MOREVECTOR, J_INITIALIZATION)) {
      IpCq().curr_slack_x_L()->Print(Jnlst(), J_MOREVECTOR, J_INITIALIZATION,
                                     "initial slack_x_L");
      IpCq().curr_slack_x_U()->Print(Jnlst(), J_MOREVECTOR, J_INITIALIZATION,
                                     "initial slack_x_U");
      IpCq().curr_slack_s_L()->Print(Jnlst(), J_MOREVECTOR, J_INITIALIZATION,
                                     "initial slack_s_L");
      IpCq().curr_slack_s_U()->Print(Jnlst(), J_MOREVECTOR, J_INITIALIZATION,
                                     "initial slack_s_U");
    }

    return true;
  }

  void WarmStartIterateInitializer::process_target_mu(Number factor,
      const Vector& curr_vars,
      const Vector& curr_slacks,
      const Vector& curr_mults,
      const Matrix& P,
      SmartPtr<const Vector>& ret_vars,
      SmartPtr<const Vector>& ret_mults)
  {
    SmartPtr<Vector> new_slacks = curr_slacks.MakeNewCopy();
    SmartPtr<Vector> new_mults = curr_mults.MakeNewCopy();
    adapt_to_target_mu(*new_slacks, *new_mults, warm_start_target_mu_);
    new_slacks->Axpy(-1, curr_slacks); // this is now correction step
    SmartPtr<Vector> new_vars = curr_vars.MakeNew();
    new_vars->Copy(curr_vars);
    P.MultVector(factor, *new_slacks, 1., *new_vars);

    ret_vars = ConstPtr(new_vars);
    ret_mults = ConstPtr(new_mults);
  }

  void WarmStartIterateInitializer::adapt_to_target_mu(Vector& new_s,
      Vector& new_z,
      Number target_mu)
  {
    DBG_ASSERT(new_s.Dim() == new_z.Dim());

    DenseVector* dnew_s = dynamic_cast<DenseVector*>(&new_s);
    assert(dnew_s);
    DenseVector* dnew_z = dynamic_cast<DenseVector*>(&new_z);
    assert(dnew_z);
    Number* values_s = dnew_s->Values();
    Number* values_z = dnew_z->Values();

    for (Index i=0; i<new_s.Dim(); i++) {
      if (values_s[i] > 1e4*values_z[i]) {
        values_z[i] = target_mu/values_s[i];
        if (values_z[i]>values_s[i]) {
          values_s[i] = values_z[i] = sqrt(target_mu);
        }
      }
      else if (values_z[i] > 1e4*values_s[i]) {
        values_s[i] = target_mu/values_z[i];
        if (values_s[i]>values_z[i]) {
          values_s[i] = values_z[i] = sqrt(target_mu);
        }
      }
      else {
        values_s[i] = values_z[i] = sqrt(target_mu);
      }
    }
  }

} // namespace Ipopt
