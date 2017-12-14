// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpDefaultIterateInitializer.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter              IBM    2004-09-23

#include "IpDefaultIterateInitializer.hpp"

namespace SimTKIpopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  DefaultIterateInitializer::DefaultIterateInitializer
  (const SmartPtr<EqMultiplierCalculator>& eq_mult_calculator,
   const SmartPtr<IterateInitializer>& warm_start_initializer)
      :
      IterateInitializer(),
      eq_mult_calculator_(eq_mult_calculator),
      warm_start_initializer_(warm_start_initializer)
  {}

  void DefaultIterateInitializer::RegisterOptions(SmartPtr<RegisteredOptions> reg_options)
  {
    reg_options->AddLowerBoundedNumberOption(
      "bound_push",
      "Desired minimum absolute distance from the initial point to bound.",
      0.0, true, 0.01,
      "Determines how much the initial point might have to "
      "be modified in order to be sufficiently inside "
      "the bounds (together with \"bound_frac\").  (This is kappa_1 in "
      "Section 3.6 of implementation paper.)");
    reg_options->AddBoundedNumberOption(
      "bound_frac",
      "Desired minimum relative distance from the initial point to bound.",
      0, true, 0.5, false, 0.01,
      "Determines how much the initial point might have to "
      "be modified in order to be sufficiently inside "
      "the bounds (together with \"bound_push\").  (This is kappa_2 in "
      "Section 3.6 of implementation paper.)");
    reg_options->AddLowerBoundedNumberOption(
      "constr_mult_init_max",
      "Maximum allowed least-square guess of constraint multipliers.",
      0, false, 1e3,
      "Determines how large the initial least-square guesses of the constraint "
      "multipliers are allowed to be (in max-norm). If the guess is larger "
      "than this value, it is discarded and all constraint multipliers are "
      "set to zero.  This options is also used when initializing the "
      "restoration phase. By default, \"resto.constr_mult_init_max\" (the one "
      "used in RestoIterateInitializer) is set to zero.");
    reg_options->AddLowerBoundedNumberOption(
      "bound_mult_init_val",
      "Initial value for the bound multipliers.",
      0, true, 1.0,
      "All dual variables corresponding to bound constraints are "
      "initialized to this value.");
    reg_options->SetRegisteringCategory("Warm Start");
    reg_options->AddStringOption2(
      "warm_start_init_point",
      "Warm-start for initial point", "no",
      "no", "do not use the warm start initialization",
      "yes", "use the warm start initialization",
      "Indicates whether this optimization should use a warm start "
      "initialization, where values of primal and dual variables are "
      "given (e.g., from a previous optimization of a related problem.)");
  }

  bool DefaultIterateInitializer::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    // Check for the algorithm options
    options.GetNumericValue("bound_push", bound_push_, prefix);
    options.GetNumericValue("bound_frac", bound_frac_, prefix);
    options.GetNumericValue("constr_mult_init_max",
                            constr_mult_init_max_, prefix);
    options.GetNumericValue("bound_mult_init_val",
                            bound_mult_init_val_, prefix);
    options.GetBoolValue("warm_start_init_point",
                         warm_start_init_point_, prefix);

    bool retvalue = true;
    if (IsValid(eq_mult_calculator_)) {
      retvalue = eq_mult_calculator_->Initialize(Jnlst(), IpNLP(), IpData(),
                 IpCq(), options, prefix);
      if (!retvalue) {
        return retvalue;
      }
    }
    if (IsValid(warm_start_initializer_)) {
      retvalue =
        warm_start_initializer_->Initialize(Jnlst(), IpNLP(), IpData(),
                                            IpCq(), options, prefix);
    }
    return retvalue;
  }

  bool DefaultIterateInitializer::SetInitialIterates()
  {
    DBG_START_METH("DefaultIterateInitializer::SetInitialIterates",
                   dbg_verbosity);

    if (warm_start_init_point_) {
      DBG_ASSERT(IsValid(warm_start_initializer_));
      return warm_start_initializer_->SetInitialIterates();
    }

    // Get the starting values provided by the NLP and store them
    // in the ip_data current fields.  The following line only requests
    // initial values for the primal variables x, but later we might
    // make this more flexible based on user options.

    /////////////////////////////////////////////////////////////////////
    //                   Initialize primal variables                   //
    /////////////////////////////////////////////////////////////////////

    IpData().InitializeDataStructures(IpNLP(), true, false, false,
                                      false, false);

    // get a container of the current point. We will modify parts of this
    // IteratesVector to set the trial point.
    SmartPtr<IteratesVector> iterates = IpData().curr()->MakeNewContainer();

    DBG_PRINT_VECTOR(2, "curr_x", *iterates->x());

    // Now we compute the initial values that the algorithm is going to
    // actually use.  We first store them in the trial fields in ip_data.

    // Push the x iterates sufficiently inside the bounds
    // Calculate any required shift in x0 and s0

    SmartPtr<const Vector> new_x;
    push_variables(Jnlst(), bound_push_, bound_frac_,
                   "x", *iterates->x(), new_x, *IpNLP().x_L(),
                   *IpNLP().x_U(), *IpNLP().Px_L(), *IpNLP().Px_U());

    iterates->Set_x(*new_x);
    IpData().set_trial(iterates);

    // Calculate the shift in s...
    SmartPtr<const Vector> s = IpCq().trial_d();
    DBG_PRINT_VECTOR(2, "s", *s);

    SmartPtr<const Vector> new_s;
    push_variables(Jnlst(), bound_push_, bound_frac_,
                   "s", *s, new_s, *IpNLP().d_L(),
                   *IpNLP().d_U(), *IpNLP().Pd_L(), *IpNLP().Pd_U());

    iterates = IpData().trial()->MakeNewContainer();
    iterates->Set_s(*new_s);

    /////////////////////////////////////////////////////////////////////
    //                   Initialize bound multipliers                  //
    /////////////////////////////////////////////////////////////////////

    // Initialize the bound multipliers to bound_mult_init_val.
    iterates->create_new_z_L();
    iterates->create_new_z_U();
    iterates->create_new_v_L();
    iterates->create_new_v_U();
    iterates->z_L_NonConst()->Set(bound_mult_init_val_);
    iterates->z_U_NonConst()->Set(bound_mult_init_val_);
    iterates->v_L_NonConst()->Set(bound_mult_init_val_);
    iterates->v_U_NonConst()->Set(bound_mult_init_val_);

    IpData().set_trial(iterates);

    /////////////////////////////////////////////////////////////////////
    //           Initialize equality constraint multipliers            //
    /////////////////////////////////////////////////////////////////////

    least_square_mults(Jnlst(), IpNLP(), IpData(), IpCq(),
                       eq_mult_calculator_, constr_mult_init_max_);

    // upgrade the trial to the current point
    IpData().AcceptTrialPoint();

    return true;
  }

  void DefaultIterateInitializer::push_variables(
    const Journalist& jnlst,
    Number bound_push,
    Number bound_frac,
    std::string name,
    const Vector& orig_x,
    SmartPtr<const Vector>& new_x,
    const Vector& x_L,
    const Vector& x_U,
    const Matrix& Px_L,
    const Matrix& Px_U)
  {
    DBG_START_FUN("DefaultIterateInitializer::push_variables",
                  dbg_verbosity);
    // Calculate any required shift in x0 and s0
    const double dbl_min = std::numeric_limits<double>::min();
    const double tiny_double = 100.0*dbl_min;

    SmartPtr<Vector> tmp = orig_x.MakeNew();
    SmartPtr<Vector> tmp_l = x_L.MakeNew();
    SmartPtr<Vector> tmp_u = x_U.MakeNew();
    SmartPtr<Vector> tiny_l = x_L.MakeNew();
    tiny_l->Set(tiny_double);

    // Calculate p_l
    SmartPtr<Vector> q_l = x_L.MakeNew();
    SmartPtr<Vector> p_l = x_L.MakeNew();

    DBG_PRINT_VECTOR(2,"orig_x", orig_x);
    DBG_PRINT_MATRIX(2,"Px_L", Px_L);
    DBG_PRINT_VECTOR(2, "x_L", x_L);
    DBG_PRINT_MATRIX(2,"Px_U", Px_U);
    DBG_PRINT_VECTOR(2, "x_U", x_U);

    Px_L.MultVector(1.0, x_L, 0.0, *tmp);
    Px_U.TransMultVector(1.0, *tmp, 0.0, *tmp_u);
    tmp_u->AddOneVector(1., x_U, -1.);
    Px_U.MultVector(1.0, *tmp_u, 0.0, *tmp);
    Px_L.TransMultVector(1.0, *tmp, 0.0, *q_l);
    q_l->AddOneVector(-1.0, *tiny_l, bound_frac);
    DBG_PRINT_VECTOR(2, "q_l", *q_l);
    // At this point, q_l is
    // bound_frac * Px_L^T Px_U(x_U - Px_U^T Px_L x_L)  -  tiny_double
    // i.e., it is bound_frac*(x_U - x_L) for those components that have
    //          two bounds
    //       and - tiny_double for those that have only one or no bounds

    tmp_l->Set(bound_push);
    p_l->AddOneVector(bound_push, x_L, 0.);
    p_l->ElementWiseAbs();
    p_l->ElementWiseMax(*tmp_l);
    // now p_l is bound_push * max(|x_L|,1)

    q_l->ElementWiseReciprocal();
    p_l->ElementWiseReciprocal();

    p_l->ElementWiseMax(*q_l);
    p_l->ElementWiseReciprocal();
    //    p_l->Axpy(1.0, *tiny_l);  we shouldn't need this

    // At this point, p_l is
    //  min(bound_push * max(|x_L|,1), bound_frac*(x_U-x_L)  for components
    //                                                       with two bounds
    //  bound_push * max(|x_L|,1)                            otherwise
    // This is the margin we want to the lower bound
    DBG_PRINT_VECTOR(1, "p_l", *p_l);

    // Calculate p_u
    SmartPtr<Vector> q_u = x_U.MakeNew();
    SmartPtr<Vector> p_u = x_U.MakeNew();
    SmartPtr<Vector> tiny_u = x_U.MakeNew();
    tiny_u->Set(tiny_double);

    Px_U.MultVector(1.0, x_U, 0.0, *tmp);
    Px_L.TransMultVector(1.0, *tmp, 0.0, *tmp_l);
    tmp_l->Axpy(-1.0, x_L);
    Px_L.MultVector(1.0, *tmp_l, 0.0, *tmp);
    Px_U.TransMultVector(1.0, *tmp, 0.0, *q_u);
    q_u->AddOneVector(-1.0, *tiny_u, bound_frac);
    DBG_PRINT_VECTOR(2,"q_u",*q_u);
    // q_u is now the same as q_l above, but of the same dimension as x_L

    tmp_u->Set(bound_push);
    p_u->Copy(x_U);
    p_u->AddOneVector(bound_push, x_U, 0.);
    p_u->ElementWiseAbs();
    p_u->ElementWiseMax(*tmp_u);
    DBG_PRINT_VECTOR(2,"p_u",*p_u);

    q_u->ElementWiseReciprocal();
    p_u->ElementWiseReciprocal();

    p_u->ElementWiseMax(*q_u);
    p_u->ElementWiseReciprocal();
    p_u->Axpy(1.0, *tiny_u);
    // At this point, p_l is
    //  min(bound_push * max(|x_U|,1), bound_frac*(x_U-x_L)  for components
    //                                                       with two bounds
    //  bound_push * max(|x_U|,1)                            otherwise
    // This is the margin we want to the upper bound
    DBG_PRINT_VECTOR(2,"actual_p_u",*p_u);

    // Calculate the new x
    SmartPtr<Vector> delta_x = orig_x.MakeNew();

    SmartPtr<Vector> zero_l = x_L.MakeNew();
    zero_l->Set(0.0);
    SmartPtr<Vector> zero_u = x_U.MakeNew();
    zero_u->Set(0.0);

    Px_L.TransMultVector(-1.0, orig_x, 0.0, *tmp_l);
    tmp_l->AddTwoVectors(1.0, x_L, 1.0, *p_l, 1.);
    tmp_l->ElementWiseMax(*zero_l);
    // tmp_l is now max(x_L + p_l - x, 0), i.e., the amount by how
    // much need to correct the variable
    Number nrm_l = tmp_l->Amax();
    if (nrm_l>0.) {
      Px_L.MultVector(1.0, *tmp_l, 0.0, *delta_x);
    }
    else {
      delta_x->Set(0.);
    }

    Px_U.TransMultVector(1.0, orig_x, 0.0, *tmp_u);
    tmp_u->AddTwoVectors(-1.0, x_U, 1.0, *p_u, 1.);
    tmp_u->ElementWiseMax(*zero_u);
    // tmp_u is now max(x - (x_U-p_u), 0), i.e., the amount by how
    // much need to correct the variable
    Number nrm_u = tmp_u->Amax();
    if (nrm_u>0.) {
      Px_U.MultVector(-1.0, *tmp_u, 1.0, *delta_x);
    }

    if (nrm_l > 0 || nrm_u > 0) {
      delta_x->Axpy(1.0, orig_x);
      new_x = ConstPtr(delta_x);
      jnlst.Printf(J_DETAILED, J_INITIALIZATION, "Moved initial values of %s sufficiently inside the bounds.\n", name.c_str());
      orig_x.Print(jnlst, J_VECTOR, J_INITIALIZATION, "original vars");
      new_x->Print(jnlst, J_VECTOR, J_INITIALIZATION, "new vars");
    }
    else {
      new_x = &orig_x;
      jnlst.Printf(J_DETAILED, J_INITIALIZATION, "Initial values of %s sufficiently inside the bounds.\n", name.c_str());
    }
  }

  void DefaultIterateInitializer::least_square_mults(
    const Journalist& jnlst,
    IpoptNLP& ip_nlp,
    IpoptData& ip_data,
    IpoptCalculatedQuantities& ip_cq,
    const SmartPtr<EqMultiplierCalculator>& eq_mult_calculator,
    Number constr_mult_init_max)
  {
    DBG_START_FUN("DefaultIterateInitializer::least_square_mults",
                  dbg_verbosity);

    SmartPtr<IteratesVector> iterates = ip_data.trial()->MakeNewContainer();
    iterates->create_new_y_c();
    iterates->create_new_y_d();

    if (iterates->y_c_NonConst()->Dim()==iterates->x()->Dim()) {
      // This problem is square, we just set the multipliers to zero
      iterates->y_c_NonConst()->Set(0.0);
      iterates->y_d_NonConst()->Set(0.0);
      ip_data.Append_info_string("s");
    }
    else if (IsValid(eq_mult_calculator) && constr_mult_init_max>0. &&
             iterates->y_c_NonConst()->Dim()+iterates->y_d_NonConst()->Dim()>0) {
      // First move all the trial data into the current fields, since
      // those values are needed to compute the initial values for
      // the multipliers
      ip_data.CopyTrialToCurrent();
      SmartPtr<Vector> y_c = iterates->y_c_NonConst();
      SmartPtr<Vector> y_d = iterates->y_d_NonConst();
      bool retval = eq_mult_calculator->CalculateMultipliers(*y_c, *y_d);
      if (!retval) {
        y_c->Set(0.0);
        y_d->Set(0.0);
      }
      else {
        /*
        {
          ip_data.set_trial(iterates);
          printf("grad_x = %e grad_s = %e y_c = %e y_d = %e\n",
             ip_cq.trial_grad_lag_x()->Amax(),
             ip_cq.trial_grad_lag_s()->Amax(),
             y_c->Amax(),
             y_d->Amax());
          iterates = ip_data.trial()->MakeNewContainer();
        }
        */
        jnlst.Printf(J_DETAILED, J_INITIALIZATION,
                     "Least square estimates max(y_c) = %e, max(y_d) = %e\n",
                     y_c->Amax(), y_d->Amax());
        Number yinitnrm = Max(y_c->Amax(), y_d->Amax());
        if (yinitnrm > constr_mult_init_max) {
          y_c->Set(0.0);
          y_d->Set(0.0);
        }
        else {
          ip_data.Append_info_string("y");
        }
      }
    }
    else {
      iterates->y_c_NonConst()->Set(0.0);
      iterates->y_d_NonConst()->Set(0.0);
    }
    ip_data.set_trial(iterates);

    DBG_PRINT_VECTOR(2, "y_c", *ip_data.trial()->y_c());
    DBG_PRINT_VECTOR(2, "y_d", *ip_data.trial()->y_d());
  }

} // namespace Ipopt
