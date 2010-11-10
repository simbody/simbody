// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpNLPScaling.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-06-25

#include "IpNLPScaling.hpp"

namespace Ipopt
{

#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  SmartPtr<Vector> NLPScalingObject::apply_vector_scaling_x_LU_NonConst(
    const Matrix& Px_LU,
    const SmartPtr<const Vector>& lu,
    const VectorSpace& x_space)
  {
    SmartPtr<Vector> scaled_x_LU = lu->MakeNew();
    if (have_x_scaling()) {
      SmartPtr<Vector> tmp_x = x_space.MakeNew();

      // move to full x space
      Px_LU.MultVector(1.0, *lu, 0.0, *tmp_x);

      // scale in full x space
      tmp_x = apply_vector_scaling_x_NonConst(ConstPtr(tmp_x));

      // move back to x_L space
      Px_LU.TransMultVector(1.0, *tmp_x, 0.0, *scaled_x_LU);
    }
    else {
      scaled_x_LU->Copy(*lu);
    }

    return scaled_x_LU;
  }

  SmartPtr<const Vector> NLPScalingObject::apply_vector_scaling_x_LU(
    const Matrix& Px_LU,
    const SmartPtr<const Vector>& lu,
    const VectorSpace& x_space)
  {
    if (have_x_scaling()) {
      return ConstPtr(apply_vector_scaling_x_LU_NonConst(Px_LU, lu, x_space));
    }
    else {
      return lu;
    }
  }

  SmartPtr<Vector> NLPScalingObject::apply_vector_scaling_d_LU_NonConst(
    const Matrix& Pd_LU,
    const SmartPtr<const Vector>& lu,
    const VectorSpace& d_space)
  {
    SmartPtr<Vector> scaled_d_LU = lu->MakeNew();
    if (have_d_scaling()) {
      SmartPtr<Vector> tmp_d = d_space.MakeNew();

      // move to full d space
      Pd_LU.MultVector(1.0, *lu, 0.0, *tmp_d);

      // scale in full x space
      tmp_d = apply_vector_scaling_d_NonConst(ConstPtr(tmp_d));

      // move back to x_L space
      Pd_LU.TransMultVector(1.0, *tmp_d, 0.0, *scaled_d_LU);
    }
    else {
      scaled_d_LU->Copy(*lu);
    }

    return scaled_d_LU;
  }

  SmartPtr<const Vector> NLPScalingObject::apply_vector_scaling_d_LU(
    const Matrix& Pd_LU,
    const SmartPtr<const Vector>& lu,
    const VectorSpace& d_space)
  {
    if (have_d_scaling()) {
      return ConstPtr(apply_vector_scaling_d_LU_NonConst(Pd_LU, lu, d_space));
    }
    else {
      return lu;
    }
  }

  SmartPtr<Vector> NLPScalingObject::unapply_vector_scaling_d_LU_NonConst(
    const Matrix& Pd_LU,
    const SmartPtr<const Vector>& lu,
    const VectorSpace& d_space)
  {
    SmartPtr<Vector> unscaled_d_LU = lu->MakeNew();
    if (have_d_scaling()) {
      SmartPtr<Vector> tmp_d = d_space.MakeNew();

      // move to full d space
      Pd_LU.MultVector(1.0, *lu, 0.0, *tmp_d);

      // scale in full x space
      tmp_d = unapply_vector_scaling_d_NonConst(ConstPtr(tmp_d));

      // move back to x_L space
      Pd_LU.TransMultVector(1.0, *tmp_d, 0.0, *unscaled_d_LU);
    }
    else {
      unscaled_d_LU->Copy(*lu);
    }

    return unscaled_d_LU;
  }

  SmartPtr<const Vector> NLPScalingObject::unapply_vector_scaling_d_LU(
    const Matrix& Pd_LU,
    const SmartPtr<const Vector>& lu,
    const VectorSpace& d_space)
  {
    if (have_d_scaling()) {
      return ConstPtr(unapply_vector_scaling_d_LU_NonConst(Pd_LU, lu, d_space));
    }
    else {
      return lu;
    }
  }

  SmartPtr<Vector> NLPScalingObject::apply_grad_obj_scaling_NonConst(
    const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> scaled_v = unapply_vector_scaling_x_NonConst(v);
    Number df = apply_obj_scaling(1.0);
    if (df != 1.) {
      scaled_v->Scal(df);
    }
    return scaled_v;
  }

  SmartPtr<const Vector> NLPScalingObject::apply_grad_obj_scaling(
    const SmartPtr<const Vector>& v)
  {
    Number df = apply_obj_scaling(1.);
    if (df != 1.) {
      SmartPtr<Vector> scaled_v = apply_grad_obj_scaling_NonConst(v);
      return ConstPtr(scaled_v);
    }
    else {
      SmartPtr<const Vector> scaled_v = unapply_vector_scaling_x(v);
      return scaled_v;
    }
  }

  SmartPtr<Vector> NLPScalingObject::unapply_grad_obj_scaling_NonConst(
    const SmartPtr<const Vector>& v)
  {
    SmartPtr<Vector> unscaled_v = apply_vector_scaling_x_NonConst(v);
    Number df = unapply_obj_scaling(1.);
    if (df != 1.) {
      unscaled_v->Scal(df);
    }
    return unscaled_v;
  }

  SmartPtr<const Vector> NLPScalingObject::unapply_grad_obj_scaling(
    const SmartPtr<const Vector>& v)
  {
    Number df = unapply_obj_scaling(1.);
    if (df != 1.) {
      SmartPtr<Vector> unscaled_v = unapply_grad_obj_scaling_NonConst(v);
      return ConstPtr(unscaled_v);
    }
    else {
      SmartPtr<const Vector> scaled_v = apply_vector_scaling_x(v);
      return scaled_v;
    }
  }

  void
  StandardScalingBase::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddNumberOption(
      "obj_scaling_factor",
      "Scaling factor for the objective function.",
      1.,
      "This option sets a scaling factor for the objective function. "
      "The scaling is seen internally by Ipopt but the unscaled objective is "
      "reported in the console output. "
      "If additional scaling parameters are computed "
      "(e.g. user-scaling or gradient-based), both factors are multiplied. "
      "If this value is chosen to be negative, Ipopt will "
      "maximize the objective function instead of minimizing it.");
  }

  bool StandardScalingBase::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    options.GetNumericValue("obj_scaling_factor", obj_scaling_factor_, prefix);
    return true;
  }

  void StandardScalingBase::DetermineScaling(
    const SmartPtr<const VectorSpace> x_space,
    const SmartPtr<const VectorSpace> c_space,
    const SmartPtr<const VectorSpace> d_space,
    const SmartPtr<const MatrixSpace> jac_c_space,
    const SmartPtr<const MatrixSpace> jac_d_space,
    const SmartPtr<const SymMatrixSpace> h_space,
    SmartPtr<const MatrixSpace>& new_jac_c_space,
    SmartPtr<const MatrixSpace>& new_jac_d_space,
    SmartPtr<const SymMatrixSpace>& new_h_space)
  {
    SmartPtr<Vector> dc;
    SmartPtr<Vector> dd;
    DetermineScalingParametersImpl(x_space, c_space, d_space,
                                   jac_c_space, jac_d_space,
                                   h_space, df_, dx_, dc, dd);

    df_ *= obj_scaling_factor_;

    if (Jnlst().ProduceOutput(J_VECTOR, J_MAIN)) {
      Jnlst().Printf(J_VECTOR, J_MAIN, "objective scaling factor = %g\n", df_);
      if (IsValid(dx_)) {
        dx_->Print(Jnlst(), J_VECTOR, J_MAIN, "x scaling vector");
      }
      else {
        Jnlst().Printf(J_VECTOR, J_MAIN, "No x scaling provided\n");
      }
      if (IsValid(dc)) {
        dc->Print(Jnlst(), J_VECTOR, J_MAIN, "c scaling vector");
      }
      else {
        Jnlst().Printf(J_VECTOR, J_MAIN, "No c scaling provided\n");
      }
      if (IsValid(dd)) {
        dd->Print(Jnlst(), J_VECTOR, J_MAIN, "d scaling vector");
      }
      else {
        Jnlst().Printf(J_VECTOR, J_MAIN, "No d scaling provided\n");
      }
    }

    // create the scaling matrix spaces
    if (IsValid(dx_) || IsValid(dc)) {
      scaled_jac_c_space_ =
        new ScaledMatrixSpace(ConstPtr(dc), false, jac_c_space,
                              ConstPtr(dx_), true);
      new_jac_c_space = GetRawPtr(scaled_jac_c_space_);
    }
    else {
      scaled_jac_c_space_ = NULL;
      new_jac_c_space = jac_c_space;
    }

    if (IsValid(dx_) || IsValid(dc)) {
      scaled_jac_d_space_ =
        new ScaledMatrixSpace(ConstPtr(dd), false, jac_d_space,
                              ConstPtr(dx_), true);
      new_jac_d_space = GetRawPtr(scaled_jac_d_space_);
    }
    else {
      scaled_jac_d_space_ = NULL;
      new_jac_d_space =jac_d_space ;
    }

    if (IsValid(h_space)) {
      if (IsValid(dx_)) {
        scaled_h_space_ = new SymScaledMatrixSpace(ConstPtr(dx_), true, h_space);
        new_h_space = GetRawPtr(scaled_h_space_);
      }
      else {
        scaled_h_space_ = NULL;
        new_h_space = h_space;
      }
    }
    else {
      new_h_space = NULL;
    }
  }

  Number StandardScalingBase::apply_obj_scaling(const Number& f)
  {
    return df_*f;
  }

  Number StandardScalingBase::unapply_obj_scaling(const Number& f)
  {
    return f/df_;
  }

  SmartPtr<Vector> StandardScalingBase::apply_vector_scaling_x_NonConst(
    const SmartPtr<const Vector>& v)
  {
    DBG_START_METH("StandardScalingBase::apply_vector_scaling_x_NonConst",
                   dbg_verbosity);
    SmartPtr<Vector> scaled_x = v->MakeNewCopy();
    if (IsValid(dx_)) {
      scaled_x->ElementWiseMultiply(*dx_);
    }
    else {
      DBG_PRINT((1, "Creating copy in apply_vector_scaling_x_NonConst!"));
    }
    return scaled_x;
  }

  SmartPtr<const Vector> StandardScalingBase::apply_vector_scaling_x(
    const SmartPtr<const Vector>& v)
  {
    if (IsValid(dx_)) {
      return ConstPtr(apply_vector_scaling_x_NonConst(v));
    }
    else {
      return v;
    }
  }

  SmartPtr<Vector> StandardScalingBase::unapply_vector_scaling_x_NonConst(
    const SmartPtr<const Vector>& v)
  {
    DBG_START_METH("StandardScalingBase::unapply_vector_scaling_x_NonConst",
                   dbg_verbosity);
    SmartPtr<Vector> unscaled_x = v->MakeNewCopy();
    if (IsValid(dx_)) {
      unscaled_x->ElementWiseDivide(*dx_);
    }
    else {
      DBG_PRINT((1, "Creating copy in unapply_vector_scaling_x_NonConst!"));
    }
    return unscaled_x;
  }

  SmartPtr<const Vector> StandardScalingBase::unapply_vector_scaling_x(
    const SmartPtr<const Vector>& v)
  {
    if (IsValid(dx_)) {
      return ConstPtr(unapply_vector_scaling_x_NonConst(v));
    }
    else {
      return v;
    }
  }

  SmartPtr<Vector> StandardScalingBase::apply_vector_scaling_c_NonConst(
    const SmartPtr<const Vector>& v)
  {
    DBG_START_METH("StandardScalingBase::apply_vector_scaling_c_NonConst",
                   dbg_verbosity);
    SmartPtr<Vector> scaled_c = v->MakeNewCopy();
    if (IsValid(scaled_jac_c_space_) &&
        IsValid(scaled_jac_c_space_->RowScaling())) {
      scaled_c->ElementWiseMultiply(*scaled_jac_c_space_->RowScaling());
    }
    else {
      DBG_PRINT((1,"Creating copy in apply_vector_scaling_c_NonConst!"));
    }
    return scaled_c;
  }

  SmartPtr<const Vector> StandardScalingBase::apply_vector_scaling_c(
    const SmartPtr<const Vector>& v)
  {
    if (IsValid(scaled_jac_c_space_) &&
        IsValid(scaled_jac_c_space_->RowScaling())) {
      return ConstPtr(apply_vector_scaling_c_NonConst(v));
    }
    else {
      return v;
    }
  }

  SmartPtr<Vector> StandardScalingBase::unapply_vector_scaling_c_NonConst(
    const SmartPtr<const Vector>& v)
  {
    DBG_START_METH("StandardScalingBase::unapply_vector_scaling_c_NonConst",
                   dbg_verbosity);
    SmartPtr<Vector> scaled_c = v->MakeNewCopy();
    if (IsValid(scaled_jac_c_space_) &&
        IsValid(scaled_jac_c_space_->RowScaling())) {
      scaled_c->ElementWiseDivide(*scaled_jac_c_space_->RowScaling());
    }
    else {
      DBG_PRINT((1,"Creating copy in unapply_vector_scaling_c_NonConst!"));
    }
    return scaled_c;
  }

  SmartPtr<const Vector> StandardScalingBase::unapply_vector_scaling_c(
    const SmartPtr<const Vector>& v)
  {
    if (IsValid(scaled_jac_c_space_) &&
        IsValid(scaled_jac_c_space_->RowScaling())) {
      return ConstPtr(unapply_vector_scaling_c_NonConst(v));
    }
    else {
      return v;
    }
  }

  SmartPtr<Vector> StandardScalingBase::apply_vector_scaling_d_NonConst(
    const SmartPtr<const Vector>& v)
  {
    DBG_START_METH("StandardScalingBase::apply_vector_scaling_d_NonConst",
                   dbg_verbosity);
    SmartPtr<Vector> scaled_d = v->MakeNewCopy();
    if (IsValid(scaled_jac_d_space_) &&
        IsValid(scaled_jac_d_space_->RowScaling())) {
      scaled_d->ElementWiseMultiply(*scaled_jac_d_space_->RowScaling());
    }
    else {
      DBG_PRINT((1,"Creating copy in apply_vector_scaling_d_NonConst!"));
    }
    return scaled_d;
  }

  SmartPtr<const Vector> StandardScalingBase::apply_vector_scaling_d(
    const SmartPtr<const Vector>& v)
  {
    if (IsValid(scaled_jac_d_space_) &&
        IsValid(scaled_jac_d_space_->RowScaling())) {
      return ConstPtr(apply_vector_scaling_d_NonConst(v));
    }
    else {
      return v;
    }
  }

  SmartPtr<Vector> StandardScalingBase::unapply_vector_scaling_d_NonConst(
    const SmartPtr<const Vector>& v)
  {
    DBG_START_METH("StandardScalingBase::unapply_vector_scaling_d_NonConst",
                   dbg_verbosity);
    SmartPtr<Vector> scaled_d = v->MakeNewCopy();
    if (IsValid(scaled_jac_d_space_) &&
        IsValid(scaled_jac_d_space_->RowScaling())) {
      scaled_d->ElementWiseDivide(*scaled_jac_d_space_->RowScaling());
    }
    else {
      DBG_PRINT((1,"Creating copy in unapply_vector_scaling_d_NonConst!"));
    }
    return scaled_d;
  }

  SmartPtr<const Vector> StandardScalingBase::unapply_vector_scaling_d(
    const SmartPtr<const Vector>& v)
  {
    if (IsValid(scaled_jac_d_space_) &&
        IsValid(scaled_jac_d_space_->RowScaling())) {
      return ConstPtr(unapply_vector_scaling_d_NonConst(v));
    }
    else {
      return v;
    }
  }

  // ToDo: matrix not passed by reference, so setting to NULL doesn't make difference
  SmartPtr<const Matrix> StandardScalingBase::apply_jac_c_scaling(
    SmartPtr<const Matrix> matrix)
  {
    if (IsValid(scaled_jac_c_space_)) {
      SmartPtr<ScaledMatrix> ret = scaled_jac_c_space_->MakeNewScaledMatrix(false);
      ret->SetUnscaledMatrix(matrix);
      return GetRawPtr(ret);
    }
    else {
      SmartPtr<const Matrix> ret = matrix;
      matrix = NULL;
      return ret;
    }
  }

  SmartPtr<const Matrix> StandardScalingBase::apply_jac_d_scaling(
    SmartPtr<const Matrix> matrix)
  {
    if (IsValid(scaled_jac_d_space_)) {
      SmartPtr<ScaledMatrix> ret = scaled_jac_d_space_->MakeNewScaledMatrix(false);
      ret->SetUnscaledMatrix(matrix);
      return GetRawPtr(ret);
    }
    else {
      SmartPtr<const Matrix> ret = matrix;
      matrix = NULL;
      return ret;
    }
  }

  SmartPtr<const SymMatrix> StandardScalingBase::apply_hessian_scaling(
    SmartPtr<const SymMatrix> matrix)
  {
    if (IsValid(scaled_h_space_)) {
      SmartPtr<SymScaledMatrix> ret = scaled_h_space_->MakeNewSymScaledMatrix(false);
      ret->SetUnscaledMatrix(matrix);
      return GetRawPtr(ret);
    }
    else {
      SmartPtr<const SymMatrix> ret = matrix;
      matrix = NULL;
      return ret;
    }
  }

  bool StandardScalingBase::have_x_scaling()
  {
    return IsValid(dx_);
  }

  bool StandardScalingBase::have_c_scaling()
  {
    return (IsValid(scaled_jac_c_space_) &&
            IsValid(scaled_jac_c_space_->RowScaling()));
  }

  bool StandardScalingBase::have_d_scaling()
  {
    return (IsValid(scaled_jac_d_space_) &&
            IsValid(scaled_jac_d_space_->RowScaling()));
  }

  void NoNLPScalingObject::DetermineScalingParametersImpl(
    const SmartPtr<const VectorSpace> x_space,
    const SmartPtr<const VectorSpace> c_space,
    const SmartPtr<const VectorSpace> d_space,
    const SmartPtr<const MatrixSpace> jac_c_space,
    const SmartPtr<const MatrixSpace> jac_d_space,
    const SmartPtr<const SymMatrixSpace> h_space,
    Number& df,
    SmartPtr<Vector>& dx,
    SmartPtr<Vector>& dc,
    SmartPtr<Vector>& dd)
  {
    df = 1.;
    dx = NULL;
    dc = NULL;
    dd = NULL;
  }

} // namespace Ipopt
