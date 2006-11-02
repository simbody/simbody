// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpGradientScaling.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-07-13

#include "IpGradientScaling.hpp"
#include "IpTripletHelper.hpp"

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

  void GradientScaling::RegisterOptions(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->AddLowerBoundedNumberOption(
      "nlp_scaling_max_gradient", "Maximum gradient after NLP scaling.",
      0, true, 100.0,
      "This is the gradient scaling cut-off. If the maximum"
      " gradient is above this value, then gradient based scaling"
      " will be performed. Scaling parameters are calculated to"
      " scale the maximum gradient back to this value. (This is g_max in "
      "Section 3.8 of the implementation paper.) Note: This"
      " option is only used if \"nlp_scaling_method\" is chosen as"
      " \"gradient-based\".");
  }

  bool GradientScaling::InitializeImpl(const OptionsList& options,
                                       const std::string& prefix)
  {
    options.GetNumericValue("nlp_scaling_max_gradient", scaling_max_gradient_, prefix);
    return StandardScalingBase::InitializeImpl(options, prefix);
  }


  void GradientScaling::DetermineScalingParametersImpl(
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
    DBG_ASSERT(IsValid(nlp_));

    SmartPtr<Vector> x = x_space->MakeNew();
    if (!nlp_->GetStartingPoint(GetRawPtr(x), true,
                                NULL, false,
                                NULL, false,
                                NULL, false,
                                NULL, false)) {
      THROW_EXCEPTION(FAILED_INITIALIZATION,
                      "Error getting initial point from NLP in GradientScaling.\n");
    }

    //
    // Calculate grad_f scaling
    //
    SmartPtr<Vector> grad_f = x_space->MakeNew();
    if (nlp_->Eval_grad_f(*x, *grad_f)) {
      double max_grad_f = grad_f->Amax();
      df = 1.;
      if (max_grad_f > scaling_max_gradient_) {
        df = scaling_max_gradient_ / max_grad_f;
      }
      Jnlst().Printf(J_DETAILED, J_INITIALIZATION,
                     "Scaling parameter for objective function = %e\n", df);
    }
    else {
      Jnlst().Printf(J_WARNING, J_INITIALIZATION,
                     "Error evaluating objective gradient at user provided starting point.\n  No scaling factor for objective function computed!\n");
      df = 1.;
    }
    //
    // No x scaling
    //
    dx = NULL;

    //
    // Calculate c scaling
    //
    SmartPtr<Matrix> jac_c = jac_c_space->MakeNew();
    if (nlp_->Eval_jac_c(*x, *jac_c)) {
      // ToDo: Don't use TripletHelper, have special methods on matrices instead
      Index nnz = TripletHelper::GetNumberEntries(*jac_c);
      Index* irow = new Index[nnz];
      Index* jcol = new Index[nnz];
      Number* values = new Number[nnz];
      TripletHelper::FillRowCol(nnz, *jac_c, irow, jcol);
      TripletHelper::FillValues(nnz, *jac_c, values);
      Number* c_scaling = new Number[jac_c->NRows()];

      for (Index r=0; r<jac_c->NRows(); r++) {
        c_scaling[r] = 0;
      }

      // put the max of each row into c_scaling...
      bool need_c_scale = false;
      for (Index i=0; i<nnz; i++) {
        if (fabs(values[i]) > scaling_max_gradient_) {
          c_scaling[irow[i]-1] = Max(c_scaling[irow[i]-1], fabs(values[i]));
          need_c_scale = true;
        }
      }

      if (need_c_scale) {
        // now compute the scaling factors for each row
        for (Index r=0; r<jac_c->NRows(); r++) {
          Number scaling = 1.0;
          if (c_scaling[r] > scaling_max_gradient_) {
            scaling = scaling_max_gradient_/c_scaling[r];
          }
          c_scaling[r] = scaling;
        }

        dc = c_space->MakeNew();
        TripletHelper::PutValuesInVector(jac_c->NRows(), c_scaling, *dc);
        if (Jnlst().ProduceOutput(J_DETAILED, J_INITIALIZATION)) {
          Jnlst().Printf(J_DETAILED, J_INITIALIZATION,
                         "Equality constraints are scaled with smallest scaling parameter is %e\n", dc->Min());
        }
      }
      else {
        Jnlst().Printf(J_DETAILED, J_INITIALIZATION,
                       "Equality constraints are not scaled.\n");
        dc = NULL;
      }

      delete [] irow;
      irow = NULL;
      delete [] jcol;
      jcol = NULL;
      delete [] values;
      values = NULL;
      delete [] c_scaling;
      c_scaling = NULL;
    }
    else {
      Jnlst().Printf(J_WARNING, J_INITIALIZATION,
                     "Error evaluating Jacobian of equality constraints at user provided starting point.\n  No scaling factors for equality constraints computed!\n");
      dc = NULL;
    }

    //
    // Calculate d scaling
    //
    SmartPtr<Matrix> jac_d = jac_d_space->MakeNew();
    if (nlp_->Eval_jac_d(*x, *jac_d)) {
      Index nnz = TripletHelper::GetNumberEntries(*jac_d);
      Index* irow = new Index[nnz];
      Index* jcol = new Index[nnz];
      Number* values = new Number[nnz];
      TripletHelper::FillRowCol(nnz, *jac_d, irow, jcol);
      TripletHelper::FillValues(nnz, *jac_d, values);
      Number* d_scaling = new Number[jac_d->NRows()];

      for (Index r=0; r<jac_d->NRows(); r++) {
        d_scaling[r] = 0;
      }

      // put the max of each row into c_scaling...
      bool need_d_scale = false;
      for (Index i=0; i<nnz; i++) {
        if (fabs(values[i]) > scaling_max_gradient_) {
          d_scaling[irow[i]-1] = Max(d_scaling[irow[i]-1], fabs(values[i]));
          need_d_scale = true;
        }
      }

      if (need_d_scale) {
        // now compute the scaling factors for each row
        for (Index r=0; r<jac_d->NRows(); r++) {
          Number scaling = 1.0;
          if (d_scaling[r] > scaling_max_gradient_) {
            scaling = scaling_max_gradient_/d_scaling[r];
          }
          d_scaling[r] = scaling;
        }

        dd = d_space->MakeNew();
        TripletHelper::PutValuesInVector(jac_d->NRows(), d_scaling, *dd);
        if (Jnlst().ProduceOutput(J_DETAILED, J_INITIALIZATION)) {
          Jnlst().Printf(J_DETAILED, J_INITIALIZATION,
                         "Inequality constraints are scaled with smallest scaling parameter is %e\n", dd->Min());
        }
      }
      else {
        dd = NULL;
        Jnlst().Printf(J_DETAILED, J_INITIALIZATION,
                       "Inequality constraints are not scaled.\n");
      }

      delete [] irow;
      irow = NULL;
      delete [] jcol;
      jcol = NULL;
      delete [] values;
      values = NULL;
      delete [] d_scaling;
      d_scaling = NULL;
    }
    else {
      Jnlst().Printf(J_WARNING, J_INITIALIZATION,
                     "Error evaluating Jacobian of inequality constraints at user provided starting point.\n  No scaling factors for inequality constraints computed!\n");
      dd = NULL;
    }
  }

} // namespace Ipopt
