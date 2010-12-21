// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpOrigIpoptNLP.cpp 765 2006-07-14 18:03:23Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpOrigIpoptNLP.hpp"
#include "IpLowRankUpdateSymMatrix.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"

namespace Ipopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  OrigIpoptNLP::OrigIpoptNLP(const SmartPtr<const Journalist>& jnlst,
                             const SmartPtr<NLP>& nlp,
                             const SmartPtr<NLPScalingObject>& nlp_scaling)
      :
      IpoptNLP(nlp_scaling),
      jnlst_(jnlst),
      nlp_(nlp),
      x_space_(NULL),
      f_cache_(1),
      grad_f_cache_(1),
      c_cache_(1),
      jac_c_cache_(1),
      d_cache_(1),
      jac_d_cache_(1),
      h_cache_(1),
      f_evals_(0),
      grad_f_evals_(0),
      c_evals_(0),
      jac_c_evals_(0),
      d_evals_(0),
      jac_d_evals_(0),
      h_evals_(0),
      initialized_(false)
  {}

  OrigIpoptNLP::~OrigIpoptNLP()
  {}

  void OrigIpoptNLP::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "bound_relax_factor",
      "Factor for initial relaxation of the bounds.",
      0, false,
      1e-8,
      "Before start of the optimization, the bounds given by the user are "
      "relaxed.  This option sets the factor for this relaxation.  If it "
      "is set to zero, then then bounds relaxation is disabled. "
      "(See Eqn.(35) in implementation paper.)");
    roptions->AddStringOption2(
      "honor_original_bounds",
      "Indicates whether final points should be projected into original bounds.",
      "yes",
      "no", "Leave final point unchanged",
      "yes", "Project final point back into original bounds",
      "Ipopt might relax the bounds during the optimization (see, e.g., option "
      "\"bound_relax_factor\").  This option determines whether the final "
      "point should be projected back into the user-provide original bounds "
      "after the optimization.");
    roptions->SetRegisteringCategory("Warm Start");
    roptions->AddStringOption2(
      "warm_start_same_structure",
      "Indicates whether a problem with a structure identical to the previous one is to be solved.",
      "no",
      "no", "Assume this is a new problem.",
      "yes", "Assume this is problem has known structure",
      "If \"yes\" is chosen, then the algorithm assumes that an NLP is now to "
      "be solved, whose strcture is identical to one that already was "
      "considered (with the same NLP object).");
    roptions->SetRegisteringCategory("NLP");
    roptions->AddStringOption2(
      "check_derivatives_for_naninf",
      "Indicates whether it is desired to check for Nan/Inf in derivative matrices",
      "no",
      "no", "Don't check (faster).",
      "yes", "Check Jacobians and Hessian for Nan and Inf.",
      "Activating this option will cause an error if an invalid number is "
      "detected in the constraint Jacobians or the Lagrangian Hessian.  If "
      "this is not activated, the test is skipped, and the algorithm might "
      "proceed with invalid numbers and fail.");
    roptions->SetRegisteringCategory("Hessian Approximation");
    roptions->AddStringOption2(
      "hessian_approximation",
      "Indicates what Hessian information is to be used.",
      "exact",
      "exact", "Use second derivatives provided by the NLP.",
      "limited-memory", "Perform a limited-memory quasi-Newton  approximation",
      "This determines which kind of information for the Hessian of the "
      "Lagrangian function is used by the algorithm.");
  }

  bool OrigIpoptNLP::Initialize(const Journalist& jnlst,
                                const OptionsList& options,
                                const std::string& prefix)
  {
    options.GetNumericValue("bound_relax_factor", bound_relax_factor_, prefix);
    options.GetBoolValue("honor_original_bounds",
                         honor_original_bounds_, prefix);
    options.GetBoolValue("warm_start_same_structure",
                         warm_start_same_structure_, prefix);
    options.GetBoolValue("check_derivatives_for_naninf",
                         check_derivatives_for_naninf_, prefix);
    Index enum_int;
    options.GetEnumValue("hessian_approximation", enum_int, prefix);
    hessian_approximation_ = HessianApproximationType(enum_int);

    // Reset the function evaluation counters (for warm start)
    f_evals_=0;
    grad_f_evals_=0;
    c_evals_=0;
    jac_c_evals_=0;
    d_evals_=0;
    jac_d_evals_=0;
    h_evals_=0;

    // Reset the cache entries belonging to a dummy dependency.  This
    // is required for repeated solve, since the cache is not updated
    // if a dimension is zero
    std::vector<const TaggedObject*> deps(1);
    deps[0] = NULL;
    std::vector<Number> sdeps(0);
    c_cache_.InvalidateResult(deps, sdeps);
    d_cache_.InvalidateResult(deps, sdeps);
    jac_c_cache_.InvalidateResult(deps, sdeps);
    jac_d_cache_.InvalidateResult(deps, sdeps);

    if (!nlp_->ProcessOptions(options, prefix)) {
      return false;
    }

    initialized_ = true;
    return IpoptNLP::Initialize(jnlst, options, prefix);
  }

  bool OrigIpoptNLP::InitializeStructures(SmartPtr<Vector>& x,
                                          bool init_x,
                                          SmartPtr<Vector>& y_c,
                                          bool init_y_c,
                                          SmartPtr<Vector>& y_d,
                                          bool init_y_d,
                                          SmartPtr<Vector>& z_L,
                                          bool init_z_L,
                                          SmartPtr<Vector>& z_U,
                                          bool init_z_U,
                                          SmartPtr<Vector>& v_L,
                                          SmartPtr<Vector>& v_U
                                         )
  {
    DBG_START_METH("OrigIpoptNLP::InitializeStructures", dbg_verbosity);
    DBG_ASSERT(initialized_);
    bool retValue;

    if (!warm_start_same_structure_) {

      retValue = nlp_->GetSpaces(x_space_, c_space_, d_space_,
                                 x_l_space_, px_l_space_,
                                 x_u_space_, px_u_space_,
                                 d_l_space_, pd_l_space_,
                                 d_u_space_, pd_u_space_,
                                 jac_c_space_, jac_d_space_,
                                 h_space_);

      if (!retValue) {
        return false;
      }

      // Check if the Hessian space is actually a limited-memory
      // approximation.  If so, get the required information from the
      // NLP and create an appropreate h_space
      if (hessian_approximation_==LIMITED_MEMORY) {
        SmartPtr<VectorSpace> approx_vecspace;
        SmartPtr<Matrix> P_approx;
        nlp_->GetQuasiNewtonApproximationSpaces(approx_vecspace,
                                                P_approx);
        if (IsValid(approx_vecspace)) {
          DBG_ASSERT(IsValid(P_approx));
          h_space_ = new LowRankUpdateSymMatrixSpace(x_space_->Dim(),
                     ConstPtr(P_approx),
                     ConstPtr(approx_vecspace),
                     true);
          jnlst_->Printf(J_DETAILED, J_INITIALIZATION,
                         "Hessian approximation will be done in smaller space of dimension %d (instead of %d)\n\n",
                         P_approx->NCols(), P_approx->NRows());
        }
        else {
          DBG_ASSERT(IsNull(P_approx));
          h_space_ = new LowRankUpdateSymMatrixSpace(x_space_->Dim(),
                     ConstPtr(P_approx),
                     ConstPtr(x_space_),
                     true);
          jnlst_->Printf(J_DETAILED, J_INITIALIZATION,
                         "Hessian approximation will be done in the space of all %d x variables.\n\n",
                         x_space_->Dim());
        }
      }

      NLP_scaling()->DetermineScaling(x_space_,
                                      c_space_, d_space_,
                                      jac_c_space_, jac_d_space_,
                                      h_space_,
                                      scaled_jac_c_space_, scaled_jac_d_space_,
                                      scaled_h_space_);

      ASSERT_EXCEPTION(x_space_->Dim() >= c_space_->Dim(), TOO_FEW_DOF,
                       "Too few degrees of freedom!");
      ASSERT_EXCEPTION(x_space_->Dim() > 0, TOO_FEW_DOF,
                       "Too few degrees of freedom (no free variables)!");

      // cannot have any null pointers, want zero length vectors
      // instead of null - this will later need to be changed for _h;
      retValue = (IsValid(x_space_) && IsValid(c_space_) && IsValid(d_space_)
                  && IsValid(x_l_space_) && IsValid(px_l_space_)
                  && IsValid(x_u_space_) && IsValid(px_u_space_)
                  && IsValid(d_u_space_) && IsValid(pd_u_space_)
                  && IsValid(d_l_space_) && IsValid(pd_l_space_)
                  && IsValid(jac_c_space_) && IsValid(jac_d_space_)
                  && IsValid(h_space_)
                  && IsValid(scaled_jac_c_space_)
                  && IsValid(scaled_jac_d_space_)
                  && IsValid(scaled_h_space_));

      DBG_ASSERT(retValue && "Model cannot return null vector or matrix prototypes or spaces,"
                 " please return zero length vectors instead");
    }
    else {
      ASSERT_EXCEPTION(IsValid(x_space_), INVALID_WARMSTART,
                       "OrigIpoptNLP called with warm_start_same_structure, but the problem is solved for the first time.");
    }

    // Create the bounds structures
    SmartPtr<Vector> x_L = x_l_space_->MakeNew();
    SmartPtr<Matrix> Px_L = px_l_space_->MakeNew();
    SmartPtr<Vector> x_U = x_u_space_->MakeNew();
    SmartPtr<Matrix> Px_U = px_u_space_->MakeNew();
    SmartPtr<Vector> d_L = d_l_space_->MakeNew();
    SmartPtr<Matrix> Pd_L = pd_l_space_->MakeNew();
    SmartPtr<Vector> d_U = d_u_space_->MakeNew();
    SmartPtr<Matrix> Pd_U = pd_u_space_->MakeNew();

    retValue = nlp_->GetBoundsInformation(*Px_L, *x_L, *Px_U, *x_U,
                                          *Pd_L, *d_L, *Pd_U, *d_U);

    if (!retValue) {
      return false;
    }

    x_L->Print(*jnlst_, J_MOREVECTOR, J_INITIALIZATION,
               "original x_L unscaled");
    x_U->Print(*jnlst_, J_MOREVECTOR, J_INITIALIZATION,
               "original x_U unscaled");
    d_L->Print(*jnlst_, J_MOREVECTOR, J_INITIALIZATION,
               "original d_L unscaled");
    d_U->Print(*jnlst_, J_MOREVECTOR, J_INITIALIZATION,
               "original d_U unscaled");

    if (honor_original_bounds_) {
      SmartPtr<Vector> tmp;
      tmp = x_L->MakeNewCopy();
      orig_x_L_ = ConstPtr(tmp);
      tmp = x_U->MakeNewCopy();
      orig_x_U_ = ConstPtr(tmp);
    }

    relax_bounds(-bound_relax_factor_, *x_L);
    relax_bounds( bound_relax_factor_, *x_U);
    relax_bounds(-bound_relax_factor_, *d_L);
    relax_bounds( bound_relax_factor_, *d_U);

    x_L_ = ConstPtr(x_L);
    Px_L_ = ConstPtr(Px_L);
    x_U_ = ConstPtr(x_U);
    Px_U_ = ConstPtr(Px_U);
    d_L_ = ConstPtr(d_L);
    Pd_L_ = ConstPtr(Pd_L);
    d_U_ = ConstPtr(d_U);
    Pd_U_ = ConstPtr(Pd_U);

    // now create and store the scaled bounds
    x_L_ = NLP_scaling()->apply_vector_scaling_x_LU(*Px_L_, x_L_, *x_space_);
    x_U_ = NLP_scaling()->apply_vector_scaling_x_LU(*Px_U_, x_U_, *x_space_);
    d_L_ = NLP_scaling()->apply_vector_scaling_d_LU(*Pd_L_, d_L_, *d_space_);
    d_U_ = NLP_scaling()->apply_vector_scaling_d_LU(*Pd_U_, d_U_, *d_space_);

    x_L->Print(*jnlst_, J_VECTOR, J_INITIALIZATION,
               "modified x_L scaled");
    x_U->Print(*jnlst_, J_VECTOR, J_INITIALIZATION,
               "modified x_U scaled");
    d_L->Print(*jnlst_, J_VECTOR, J_INITIALIZATION,
               "modified d_L scaled");
    d_U->Print(*jnlst_, J_VECTOR, J_INITIALIZATION,
               "modified d_U scaled");

    // Create the iterates structures
    x = x_space_->MakeNew();
    y_c = c_space_->MakeNew();
    y_d = d_space_->MakeNew();
    z_L = x_l_space_->MakeNew();
    z_U = x_u_space_->MakeNew();
    v_L = d_l_space_->MakeNew();
    v_U = d_u_space_->MakeNew();

    retValue = nlp_->GetStartingPoint(GetRawPtr(x), init_x,
                                      GetRawPtr(y_c), init_y_c,
                                      GetRawPtr(y_d), init_y_d,
                                      GetRawPtr(z_L), init_z_L,
                                      GetRawPtr(z_U), init_z_U);

    if (!retValue) {
      return false;
    }


    Number obj_scal = NLP_scaling()->apply_obj_scaling(1.);
    if (init_x) {
      x->Print(*jnlst_, J_VECTOR, J_INITIALIZATION, "initial x unscaled");
      if (NLP_scaling()->have_x_scaling()) {
        x = NLP_scaling()->apply_vector_scaling_x_NonConst(ConstPtr(x));
      }
    }
    if (init_y_c) {
      y_c->Print(*jnlst_, J_VECTOR, J_INITIALIZATION, "initial y_c unscaled");
      if (NLP_scaling()->have_c_scaling()) {
        y_c = NLP_scaling()->unapply_vector_scaling_c_NonConst(ConstPtr(y_c));
      }
      if (obj_scal!=1.) {
        y_c->Scal(obj_scal);
      }
    }
    if (init_y_d) {
      y_d->Print(*jnlst_, J_VECTOR, J_INITIALIZATION, "initial y_d unscaled");
      if (NLP_scaling()->have_d_scaling()) {
        y_d = NLP_scaling()->unapply_vector_scaling_d_NonConst(ConstPtr(y_d));
      }
      if (obj_scal!=1.) {
        y_d->Scal(obj_scal);
      }
    }
    if (init_z_L) {
      z_L->Print(*jnlst_, J_VECTOR, J_INITIALIZATION, "initial z_L unscaled");
      if (NLP_scaling()->have_x_scaling()) {
        z_L = NLP_scaling()->apply_vector_scaling_x_LU_NonConst(*Px_L_, ConstPtr(z_L), *x_space_);
      }
      if (obj_scal!=1.) {
        z_L->Scal(obj_scal);
      }
    }
    if (init_z_U) {
      z_U->Print(*jnlst_, J_VECTOR, J_INITIALIZATION, "initial z_U unscaled");
      if (NLP_scaling()->have_x_scaling()) {
        z_U = NLP_scaling()->apply_vector_scaling_x_LU_NonConst(*Px_U_, ConstPtr(z_U), *x_space_);
      }
      if (obj_scal!=1.) {
        z_U->Scal(obj_scal);
      }
    }

    return true;
  }

  void
  OrigIpoptNLP::relax_bounds(Number bound_relax_factor, Vector& bounds)
  {
    DBG_START_METH("OrigIpoptNLP::relax_bounds", dbg_verbosity);
    if (bound_relax_factor!=0.) {
      SmartPtr<Vector> tmp = bounds.MakeNew();
      tmp->Copy(bounds);
      tmp->ElementWiseAbs();
      SmartPtr<Vector> ones = bounds.MakeNew();
      ones->Set(1.);
      tmp->ElementWiseMax(*ones);
      DBG_PRINT((1, "bound_relax_factor = %e", bound_relax_factor));
      DBG_PRINT_VECTOR(2, "tmp", *tmp);
      DBG_PRINT_VECTOR(2, "bounds before", bounds);
      bounds.Axpy(bound_relax_factor, *tmp);
      DBG_PRINT_VECTOR(2, "bounds after", bounds);
    }
  }

  Number OrigIpoptNLP::f(const Vector& x)
  {
    x.Dot(x);
    x.Dot(x);
    DBG_START_METH("OrigIpoptNLP::f", dbg_verbosity);
    Number ret = 0.0;
    DBG_PRINT((2, "x.Tag = %d\n", x.GetTag()));
    if (!f_cache_.GetCachedResult1Dep(ret, &x)) {
      f_evals_++;
      SmartPtr<const Vector> unscaled_x = NLP_scaling()->unapply_vector_scaling_x(&x);
      f_eval_time_.Start();
      bool success = nlp_->Eval_f(*unscaled_x, ret);
      f_eval_time_.End();
      DBG_PRINT((1, "success = %d ret = %e\n", success, ret));
      ASSERT_EXCEPTION(success && IsFiniteNumber(ret), Eval_Error,
                       "Error evaluating the objective function");
      ret = NLP_scaling()->apply_obj_scaling(ret);
      f_cache_.AddCachedResult1Dep(ret, &x);
    }

    return ret;
  }

  Number OrigIpoptNLP::f(const Vector& x, Number mu)
  {
    assert(false && "ERROR: This method is only a placeholder for f(mu) and should not be called");
    return 0.;
  }

  SmartPtr<const Vector> OrigIpoptNLP::grad_f(const Vector& x)
  {
    SmartPtr<Vector> unscaled_grad_f;
    SmartPtr<const Vector> retValue;
    if (!grad_f_cache_.GetCachedResult1Dep(retValue, &x)) {
      grad_f_evals_++;
      unscaled_grad_f = x_space_->MakeNew();

      SmartPtr<const Vector> unscaled_x = NLP_scaling()->unapply_vector_scaling_x(&x);
      grad_f_eval_time_.Start();
      bool success = nlp_->Eval_grad_f(*unscaled_x, *unscaled_grad_f);
      grad_f_eval_time_.End();
      ASSERT_EXCEPTION(success && IsFiniteNumber(unscaled_grad_f->Nrm2()),
                       Eval_Error, "Error evaluating the gradient of the objective function");
      retValue = NLP_scaling()->apply_grad_obj_scaling(ConstPtr(unscaled_grad_f));
      grad_f_cache_.AddCachedResult1Dep(retValue, &x);
    }

    return retValue;
  }

  SmartPtr<const Vector> OrigIpoptNLP::grad_f(const Vector& x, Number mu)
  {
    THROW_EXCEPTION(INTERNAL_ABORT,
                    "ERROR: This method is only a placeholder for grad_f(mu) and should not be called");
    return NULL;
  }

  /** Equality constraint residual */
  SmartPtr<const Vector> OrigIpoptNLP::c(const Vector& x)
  {
    SmartPtr<const Vector> retValue;
    if (c_space_->Dim()==0) {
      // We do this caching of an empty vector so that the returned
      // Vector has always the same tag (this might make a difference
      // in cases where only the constraints are supposed to change...
      SmartPtr<const Vector> dep = NULL;
      if (!c_cache_.GetCachedResult1Dep(retValue, GetRawPtr(dep))) {
        retValue = c_space_->MakeNew();
        c_cache_.AddCachedResult1Dep(retValue, GetRawPtr(dep));
      }
    }
    else {
      if (!c_cache_.GetCachedResult1Dep(retValue, x)) {
        SmartPtr<Vector> unscaled_c = c_space_->MakeNew();
        c_evals_++;
        SmartPtr<const Vector> unscaled_x = NLP_scaling()->unapply_vector_scaling_x(&x);
        c_eval_time_.Start();
        bool success = nlp_->Eval_c(*unscaled_x, *unscaled_c);
        c_eval_time_.End();
        ASSERT_EXCEPTION(success && IsFiniteNumber(unscaled_c->Nrm2()),
                         Eval_Error, "Error evaluating the equality constraints");
        retValue = NLP_scaling()->apply_vector_scaling_c(ConstPtr(unscaled_c));
        c_cache_.AddCachedResult1Dep(retValue, x);
      }
    }

    return retValue;
  }

  SmartPtr<const Vector> OrigIpoptNLP::d(const Vector& x)
  {
    DBG_START_METH("OrigIpoptNLP::d", dbg_verbosity);
    SmartPtr<const Vector> retValue;
    if (d_space_->Dim()==0) {
      // We do this caching of an empty vector so that the returned
      // Vector has always the same tag (this might make a difference
      // in cases where only the constraints are supposed to change...
      SmartPtr<const Vector> dep = NULL;
      if (!d_cache_.GetCachedResult1Dep(retValue, GetRawPtr(dep))) {
        retValue = d_space_->MakeNew();
        d_cache_.AddCachedResult1Dep(retValue, GetRawPtr(dep));
      }
    }
    else {
      if (!d_cache_.GetCachedResult1Dep(retValue, x)) {
        d_evals_++;
        SmartPtr<Vector> unscaled_d = d_space_->MakeNew();

        DBG_PRINT_VECTOR(2, "scaled_x", x);
        SmartPtr<const Vector> unscaled_x = NLP_scaling()->unapply_vector_scaling_x(&x);
        d_eval_time_.Start();
        bool success = nlp_->Eval_d(*unscaled_x, *unscaled_d);
        d_eval_time_.End();
        DBG_PRINT_VECTOR(2, "unscaled_d", *unscaled_d);
        ASSERT_EXCEPTION(success && IsFiniteNumber(unscaled_d->Nrm2()),
                         Eval_Error, "Error evaluating the inequality constraints");
        retValue = NLP_scaling()->apply_vector_scaling_d(ConstPtr(unscaled_d));
        d_cache_.AddCachedResult1Dep(retValue, x);
      }
    }

    return retValue;
  }

  SmartPtr<const Matrix> OrigIpoptNLP::jac_c(const Vector& x)
  {
    SmartPtr<const Matrix> retValue;
    if (c_space_->Dim()==0) {
      // We do this caching of an empty vector so that the returned
      // Matrix has always the same tag
      SmartPtr<const Vector> dep = NULL;
      if (!jac_c_cache_.GetCachedResult1Dep(retValue, GetRawPtr(dep))) {
        SmartPtr<Matrix> unscaled_jac_c = jac_c_space_->MakeNew();
        retValue = NLP_scaling()->apply_jac_c_scaling(ConstPtr(unscaled_jac_c));
        jac_c_cache_.AddCachedResult1Dep(retValue, GetRawPtr(dep));
      }
    }
    else {
      if (!jac_c_cache_.GetCachedResult1Dep(retValue, x)) {
        jac_c_evals_++;
        SmartPtr<Matrix> unscaled_jac_c = jac_c_space_->MakeNew();

        SmartPtr<const Vector> unscaled_x = NLP_scaling()->unapply_vector_scaling_x(&x);
        jac_c_eval_time_.Start();
        bool success = nlp_->Eval_jac_c(*unscaled_x, *unscaled_jac_c);
        jac_c_eval_time_.End();
        ASSERT_EXCEPTION(success, Eval_Error, "Error evaluating the jacobian of the equality constraints");
        if (check_derivatives_for_naninf_) {
          if (!unscaled_jac_c->HasValidNumbers()) {
            jnlst_->Printf(J_WARNING, J_NLP,
                           "The Jacobian for the equality constraints contains an invalid number\n");
            THROW_EXCEPTION(Eval_Error, "The Jacobian for the equality constraints contains an invalid number");
          }
        }
        retValue = NLP_scaling()->apply_jac_c_scaling(ConstPtr(unscaled_jac_c));
        jac_c_cache_.AddCachedResult1Dep(retValue, x);
      }
    }

    return retValue;
  }

  SmartPtr<const Matrix> OrigIpoptNLP::jac_d(const Vector& x)
  {
    SmartPtr<const Matrix> retValue;
    if (d_space_->Dim()==0) {
      // We do this caching of an empty vector so that the returned
      // Matrix has always the same tag
      SmartPtr<const Vector> dep = NULL;
      if (!jac_d_cache_.GetCachedResult1Dep(retValue, GetRawPtr(dep))) {
        SmartPtr<Matrix> unscaled_jac_d = jac_d_space_->MakeNew();
        retValue = NLP_scaling()->apply_jac_d_scaling(ConstPtr(unscaled_jac_d));
        jac_d_cache_.AddCachedResult1Dep(retValue, GetRawPtr(dep));
      }
    }
    else {
      if (!jac_d_cache_.GetCachedResult1Dep(retValue, x)) {
        jac_d_evals_++;
        SmartPtr<Matrix> unscaled_jac_d = jac_d_space_->MakeNew();

        SmartPtr<const Vector> unscaled_x = NLP_scaling()->unapply_vector_scaling_x(&x);
        jac_d_eval_time_.Start();
        bool success = nlp_->Eval_jac_d(*unscaled_x, *unscaled_jac_d);
        jac_d_eval_time_.End();
        ASSERT_EXCEPTION(success, Eval_Error, "Error evaluating the jacobian of the inequality constraints");
        retValue = NLP_scaling()->apply_jac_d_scaling(ConstPtr(unscaled_jac_d));
        jac_d_cache_.AddCachedResult1Dep(retValue, x);
      }
    }

    return retValue;
  }

  SmartPtr<const SymMatrix> OrigIpoptNLP::uninitialized_h()
  {
    return h_space_->MakeNewSymMatrix();
  }

  SmartPtr<const SymMatrix> OrigIpoptNLP::h(const Vector& x,
      Number obj_factor,
      const Vector& yc,
      const Vector& yd)
  {
    std::vector<const TaggedObject*> deps(3);
    deps[0] = &x;
    deps[1] = &yc;
    deps[2] = &yd;
    std::vector<Number> scalar_deps(1);
    scalar_deps[0] = obj_factor;

    SmartPtr<SymMatrix> unscaled_h;
    SmartPtr<const SymMatrix> retValue;
    if (!h_cache_.GetCachedResult(retValue, deps, scalar_deps)) {
      h_evals_++;
      unscaled_h = h_space_->MakeNewSymMatrix();

      SmartPtr<const Vector> unscaled_x = NLP_scaling()->unapply_vector_scaling_x(&x);
      SmartPtr<const Vector> unscaled_yc = NLP_scaling()->apply_vector_scaling_c(&yc);
      SmartPtr<const Vector> unscaled_yd = NLP_scaling()->apply_vector_scaling_d(&yd);
      Number scaled_obj_factor = NLP_scaling()->apply_obj_scaling(obj_factor);
      h_eval_time_.Start();
      bool success = nlp_->Eval_h(*unscaled_x, scaled_obj_factor, *unscaled_yc, *unscaled_yd, *unscaled_h);
      h_eval_time_.End();
      ASSERT_EXCEPTION(success, Eval_Error, "Error evaluating the hessian of the lagrangian");
      retValue = NLP_scaling()->apply_hessian_scaling(ConstPtr(unscaled_h));
      h_cache_.AddCachedResult(retValue, deps, scalar_deps);
    }

    return retValue;
  }

  SmartPtr<const SymMatrix> OrigIpoptNLP::h(const Vector& x,
      Number obj_factor,
      const Vector& yc,
      const Vector& yd,
      Number mu)
  {
    THROW_EXCEPTION(INTERNAL_ABORT,
                    "ERROR: This method is only a for h(mu) and should not be called");
    return NULL;
  }


  void OrigIpoptNLP::GetSpaces(SmartPtr<const VectorSpace>& x_space,
                               SmartPtr<const VectorSpace>& c_space,
                               SmartPtr<const VectorSpace>& d_space,
                               SmartPtr<const VectorSpace>& x_l_space,
                               SmartPtr<const MatrixSpace>& px_l_space,
                               SmartPtr<const VectorSpace>& x_u_space,
                               SmartPtr<const MatrixSpace>& px_u_space,
                               SmartPtr<const VectorSpace>& d_l_space,
                               SmartPtr<const MatrixSpace>& pd_l_space,
                               SmartPtr<const VectorSpace>& d_u_space,
                               SmartPtr<const MatrixSpace>& pd_u_space,
                               SmartPtr<const MatrixSpace>& Jac_c_space,
                               SmartPtr<const MatrixSpace>& Jac_d_space,
                               SmartPtr<const SymMatrixSpace>& Hess_lagrangian_space)
  {
    // Make sure that we already have all the pointers
    DBG_ASSERT(IsValid(x_space_) &&
               IsValid(c_space_) &&
               IsValid(d_space_) &&
               IsValid(x_l_space_) &&
               IsValid(px_l_space_) &&
               IsValid(x_u_space_) &&
               IsValid(px_u_space_) &&
               IsValid(d_l_space_) &&
               IsValid(pd_l_space_) &&
               IsValid(d_u_space_) &&
               IsValid(pd_u_space_) &&
               IsValid(scaled_jac_c_space_) &&
               IsValid(scaled_jac_d_space_) &&
               IsValid(scaled_h_space_));

    DBG_ASSERT(IsValid(NLP_scaling()));

    x_space = x_space_;
    c_space = c_space_;
    d_space = d_space_;
    x_l_space = x_l_space_;
    px_l_space = px_l_space_;
    x_u_space = x_u_space_;
    px_u_space = px_u_space_;
    d_l_space = d_l_space_;
    pd_l_space = pd_l_space_;
    d_u_space = d_u_space_;
    pd_u_space = pd_u_space_;
    Jac_c_space = scaled_jac_c_space_;
    Jac_d_space = scaled_jac_d_space_;
    Hess_lagrangian_space = scaled_h_space_;
  }

  void OrigIpoptNLP::FinalizeSolution(SolverReturn status,
                                      const Vector& x, const Vector& z_L, const Vector& z_U,
                                      const Vector& c, const Vector& d,
                                      const Vector& y_c, const Vector& y_d,
                                      Number obj_value)
  {
    DBG_START_METH("OrigIpoptNLP::FinalizeSolution", dbg_verbosity);
    // need to submit the unscaled solution back to the nlp
    SmartPtr<const Vector> unscaled_x =
      NLP_scaling()->unapply_vector_scaling_x(&x);
    SmartPtr<const Vector> unscaled_c =
      NLP_scaling()->unapply_vector_scaling_c(&c);
    SmartPtr<const Vector> unscaled_d =
      NLP_scaling()->unapply_vector_scaling_d(&d);
    const Number unscaled_obj = NLP_scaling()->unapply_obj_scaling(obj_value);

    SmartPtr<const Vector> unscaled_z_L;
    SmartPtr<const Vector> unscaled_z_U;
    SmartPtr<const Vector> unscaled_y_c;
    SmartPtr<const Vector> unscaled_y_d;

    // The objective function scaling factor also appears in the constraints
    Number obj_unscale_factor = NLP_scaling()->unapply_obj_scaling(1.);
    if (obj_unscale_factor!=1.) {
      SmartPtr<Vector> tmp = NLP_scaling()->apply_vector_scaling_x_LU_NonConst(*Px_L_, &z_L, *x_space_);
      tmp->Scal(obj_unscale_factor);
      unscaled_z_L = ConstPtr(tmp);

      tmp = NLP_scaling()->apply_vector_scaling_x_LU_NonConst(*Px_U_, &z_U, *x_space_);
      tmp->Scal(obj_unscale_factor);
      unscaled_z_U = ConstPtr(tmp);

      tmp = NLP_scaling()->apply_vector_scaling_c_NonConst(&y_c);
      tmp->Scal(obj_unscale_factor);
      unscaled_y_c = ConstPtr(tmp);

      tmp = NLP_scaling()->apply_vector_scaling_d_NonConst(&y_d);
      tmp->Scal(obj_unscale_factor);
      unscaled_y_d = ConstPtr(tmp);
    }
    else {
      unscaled_z_L = NLP_scaling()->apply_vector_scaling_x_LU(*Px_L_, &z_L, *x_space_);
      unscaled_z_U = NLP_scaling()->apply_vector_scaling_x_LU(*Px_U_, &z_U, *x_space_);
      unscaled_y_c = NLP_scaling()->apply_vector_scaling_c(&y_c);
      unscaled_y_d = NLP_scaling()->apply_vector_scaling_d(&y_d);
    }

    if (honor_original_bounds_ && (Px_L_->NCols()>0 || Px_U_->NCols()>0)) {
      // Make sure the user specified bounds are satisfied
      SmartPtr<Vector> tmp;
      SmartPtr<Vector> un_x = unscaled_x->MakeNewCopy();
      if (Px_L_->NCols()>0) {
        tmp = orig_x_L_->MakeNewCopy();
        Px_L_->TransMultVector(1., *un_x, 0., *tmp);
        Px_L_->MultVector(-1., *tmp, 1., *un_x);
        tmp->ElementWiseMax(*orig_x_L_);
        Px_L_->MultVector(1., *tmp, 1., *un_x);
      }
      if (Px_U_->NCols()>0) {
        tmp = orig_x_U_->MakeNewCopy();
        Px_U_->TransMultVector(1., *un_x, 0., *tmp);
        Px_U_->MultVector(-1., *tmp, 1., *un_x);
        tmp->ElementWiseMin(*orig_x_U_);
        Px_U_->MultVector(1., *tmp, 1., *un_x);
      }
      unscaled_x = ConstPtr(un_x);
    }

    unscaled_x->Print(*jnlst_, J_VECTOR, J_SOLUTION, "final x unscaled");
    unscaled_y_c->Print(*jnlst_, J_VECTOR, J_SOLUTION, "final y_c unscaled");
    unscaled_y_d->Print(*jnlst_, J_VECTOR, J_SOLUTION, "final y_d unscaled");
    unscaled_z_L->Print(*jnlst_, J_VECTOR, J_SOLUTION, "final z_L unscaled");
    unscaled_z_U->Print(*jnlst_, J_VECTOR, J_SOLUTION, "final z_U unscaled");

    nlp_->FinalizeSolution(status, *unscaled_x,
                           *unscaled_z_L, *unscaled_z_U,
                           *unscaled_c, *unscaled_d,
                           *unscaled_y_c, *unscaled_y_d,
                           unscaled_obj);
  }

  bool OrigIpoptNLP::IntermediateCallBack(AlgorithmMode mode,
                                          Index iter, Number obj_value,
                                          Number inf_pr, Number inf_du,
                                          Number mu, Number d_norm,
                                          Number regularization_size,
                                          Number alpha_du, Number alpha_pr,
                                          Index ls_trials,
                                          SmartPtr<const IpoptData> ip_data,
                                          SmartPtr<IpoptCalculatedQuantities> ip_cq)
  {
    return nlp_->IntermediateCallBack(mode, iter, obj_value, inf_pr, inf_du,
                                      mu, d_norm, regularization_size,
                                      alpha_du, alpha_pr, ls_trials,
                                      GetRawPtr(ip_data), GetRawPtr(ip_cq));
  }

  void OrigIpoptNLP::AdjustVariableBounds(const Vector& new_x_L, const Vector& new_x_U,
                                          const Vector& new_d_L, const Vector& new_d_U)
  {
    x_L_ = new_x_L.MakeNewCopy();
    x_U_ = new_x_U.MakeNewCopy();
    d_L_ = new_d_L.MakeNewCopy();
    d_U_ = new_d_U.MakeNewCopy();
  }

  void
  OrigIpoptNLP::PrintTimingStatistics(
    Journalist& jnlst,
    EJournalLevel level,
    EJournalCategory category) const
  {
    if (!jnlst.ProduceOutput(level, category))
      return;

    jnlst.Printf(level, category,
                 "Function Evaluations................: %10.3f\n",
                 TotalFunctionEvaluationCPUTime());
    jnlst.Printf(level, category,
                 " Objective function.................: %10.3f\n",
                 f_eval_time_.TotalTime());
    jnlst.Printf(level, category,
                 " Equality constraints...............: %10.3f\n",
                 c_eval_time_.TotalTime());
    jnlst.Printf(level, category,
                 " Inequality constraints.............: %10.3f\n",
                 d_eval_time_.TotalTime());
    jnlst.Printf(level, category,
                 " Equality constraint Jacobian.......: %10.3f\n",
                 jac_c_eval_time_.TotalTime());
    jnlst.Printf(level, category,
                 " Inequality constraint Jacobian.....: %10.3f\n",
                 jac_d_eval_time_.TotalTime());
    jnlst.Printf(level, category,
                 " Lagrangian Hessian.................: %10.3f\n",
                 h_eval_time_.TotalTime());
  }

  Number
  OrigIpoptNLP::TotalFunctionEvaluationCPUTime() const
  {
    return f_eval_time_.TotalTime()+
           c_eval_time_.TotalTime()+
           d_eval_time_.TotalTime()+
           jac_c_eval_time_.TotalTime()+
           jac_d_eval_time_.TotalTime()+
           h_eval_time_.TotalTime();
  }

} // namespace Ipopt
