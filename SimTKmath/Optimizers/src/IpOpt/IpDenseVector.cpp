// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpDenseVector.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpDenseVector.hpp"
#include "IpBlas.hpp"
#include "IpUtils.hpp"
#include "IpDebug.hpp"

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

  DenseVector::DenseVector(const DenseVectorSpace* owner_space)
      :
      Vector(owner_space),
      owner_space_(owner_space),
      values_(NULL),
      initialized_(false)
  {
    DBG_START_METH("DenseVector::DenseVector(Index dim)", dbg_verbosity);
    if (Dim() == 0) {
      initialized_ = true;
      homogeneous_ = false;
    }
  }

  DenseVector::~DenseVector()
  {
    DBG_START_METH("DenseVector::~DenseVector()", dbg_verbosity);
    if (values_) {
      owner_space_->FreeInternalStorage(values_);
    }
  }

  void DenseVector::SetValues(const Number* x)
  {
    initialized_ = true;
    IpBlasDcopy(Dim(), x, 1, values_allocated(), 1);
    homogeneous_ = false;
    // This is not an overloaded method from
    // Vector. Here, we must call ObjectChanged()
    // manually.
    ObjectChanged();
  }

  void DenseVector::set_values_from_scalar()
  {
    DBG_ASSERT(homogeneous_);
    initialized_ = true;
    homogeneous_ = false;
    Number* vals = values_allocated();
    IpBlasDcopy(Dim(), &scalar_, 0, vals, 1);
  }

  void DenseVector::CopyImpl(const Vector& x)
  {
    DBG_START_METH("DenseVector::CopyImpl(const Vector& x)", dbg_verbosity);
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x);
    if (dense_x) {
      DBG_ASSERT(dense_x->initialized_);
      DBG_ASSERT(Dim() == dense_x->Dim());
      homogeneous_ = dense_x->homogeneous_;
      if (homogeneous_) {
        scalar_ = dense_x->scalar_;
      }
      else {
        IpBlasDcopy(Dim(), dense_x->values_, 1, values_allocated(), 1);
      }
    }
    initialized_=true;
  }

  void DenseVector::ScalImpl(Number alpha)
  {
    DBG_ASSERT(initialized_);
    if (homogeneous_) {
      scalar_ *= alpha;
    }
    else {
      IpBlasDscal(Dim(), alpha, values_, 1);
    }
  }

  void DenseVector::AxpyImpl(Number alpha, const Vector &x)
  {
    DBG_ASSERT(initialized_);
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x);
    if (dense_x) {
      DBG_ASSERT(dense_x->initialized_);
      DBG_ASSERT(Dim() == dense_x->Dim());
      if (homogeneous_) {
        if (dense_x->homogeneous_) {
          scalar_ += alpha * dense_x->scalar_;
        }
        else {
          homogeneous_ = false;
          Number* vals = values_allocated();
          for (Index i=0; i<Dim(); i++) {
            vals[i] = scalar_ + alpha*dense_x->values_[i];
          }
        }
      }
      else {
        if (dense_x->homogeneous_) {
          if (dense_x->scalar_!=0.) {
            IpBlasDaxpy(Dim(), alpha, &dense_x->scalar_, 0, values_, 1);
          }
        }
        else {
          IpBlasDaxpy(Dim(), alpha, dense_x->values_, 1, values_, 1);
        }
      }
    }
  }

  Number DenseVector::DotImpl(const Vector &x) const
  {
    DBG_ASSERT(initialized_);
    Number retValue;
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x);
    DBG_ASSERT(dense_x->initialized_);
    DBG_ASSERT(Dim() == dense_x->Dim());
    if (homogeneous_) {
      if (dense_x->homogeneous_) {
        retValue = Dim() * scalar_ * dense_x->scalar_;
      }
      else {
        retValue = IpBlasDdot(Dim(), dense_x->values_, 1, &scalar_, 0);
      }
    }
    else {
      if (dense_x->homogeneous_) {
        retValue = IpBlasDdot(Dim(), &dense_x->scalar_, 0, values_, 1);
      }
      else {
        retValue = IpBlasDdot(Dim(), dense_x->values_, 1, values_, 1);
      }
    }
    return retValue;
  }

  Number DenseVector::Nrm2Impl() const
  {
    DBG_ASSERT(initialized_);
    if (homogeneous_) {
      return sqrt((double)Dim()) * fabs(scalar_);
    }
    else {
      return IpBlasDnrm2(Dim(), values_, 1);
    }
  }

  Number DenseVector::AsumImpl() const
  {
    DBG_ASSERT(initialized_);
    if (homogeneous_) {
      return Dim() * fabs(scalar_);
    }
    else {
      return IpBlasDasum(Dim(), values_, 1);
    }
  }

  Number DenseVector::AmaxImpl() const
  {
    DBG_ASSERT(initialized_);
    if (Dim()==0) {
      return 0.;
    }
    else {
      if (homogeneous_) {
        return fabs(scalar_);
      }
      else {
        return fabs(values_[IpBlasIdamax(Dim(), values_, 1)-1]);
      }
    }
  }

  void DenseVector::SetImpl(Number value)
  {
    initialized_ = true;
    homogeneous_ = true;
    scalar_ = value;
    // ToDo decide if we want this here:
    if (values_) {
      owner_space_->FreeInternalStorage(values_);
      values_ = NULL;
    }
  }

  void DenseVector::ElementWiseDivideImpl(const Vector& x)
  {
    DBG_ASSERT(initialized_);
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x);
    if (dense_x) {
      DBG_ASSERT(dense_x->initialized_);
      const Number* values_x = dense_x->values_;
      DBG_ASSERT(Dim() == dense_x->Dim());
      if (homogeneous_) {
        if (dense_x->homogeneous_) {
          scalar_ /= dense_x->scalar_;
        }
        else {
          homogeneous_ = false;
          Number* vals = values_allocated();
          for (Index i=0; i<Dim(); i++) {
            vals[i] = scalar_/values_x[i];
          }
        }
      }
      else {
        if (dense_x->homogeneous_) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] /= dense_x->scalar_;
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] /= values_x[i];
          }
        }
      }
    }
  }

  void DenseVector::ElementWiseMultiplyImpl(const Vector& x)
  {
    DBG_ASSERT(initialized_);
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x);
    if (dense_x) {
      DBG_ASSERT(dense_x->initialized_);
      const Number* values_x = dense_x->values_;
      DBG_ASSERT(Dim() == dense_x->Dim());
      if (homogeneous_) {
        if (dense_x->homogeneous_) {
          scalar_ *= dense_x->scalar_;
        }
        else {
          homogeneous_ = false;
          Number* vals = values_allocated();
          for (Index i=0; i<Dim(); i++) {
            vals[i] = scalar_*values_x[i];
          }
        }
      }
      else {
        if (dense_x->homogeneous_) {
          if (dense_x->scalar_ != 1.0) {
            for (Index i=0; i<Dim(); i++) {
              values_[i] *= dense_x->scalar_;
            }
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] *= values_x[i];
          }
        }
      }
    }
  }

  void DenseVector::ElementWiseMaxImpl(const Vector& x)
  {
    DBG_ASSERT(initialized_);
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    assert(dense_x); // ToDo: Implement Others
    if (dense_x) {
      DBG_ASSERT(dense_x->initialized_);
      const Number* values_x = dense_x->values_;
      DBG_ASSERT(Dim() == dense_x->Dim());
      if (homogeneous_) {
        if (dense_x->homogeneous_) {
          scalar_ = Ipopt::Max(scalar_, dense_x->scalar_);
        }
        else {
          homogeneous_ = false;
          Number* vals = values_allocated();
          for (Index i=0; i<Dim(); i++) {
            vals[i] = Ipopt::Max(scalar_, values_x[i]);
          }
        }
      }
      else {
        if (dense_x->homogeneous_) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = Ipopt::Max(values_[i], dense_x->scalar_);
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = Ipopt::Max(values_[i], values_x[i]);
          }
        }
      }
    }
  }

  void DenseVector::ElementWiseMinImpl(const Vector& x)
  {
    DBG_ASSERT(initialized_);
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x);
    if (dense_x) {
      DBG_ASSERT(dense_x->initialized_);
      const Number* values_x = dense_x->values_;
      DBG_ASSERT(Dim() == dense_x->Dim());
      if (homogeneous_) {
        if (dense_x->homogeneous_) {
          scalar_ = Ipopt::Min(scalar_, dense_x->scalar_);
        }
        else {
          homogeneous_ = false;
          Number* vals = values_allocated();
          for (Index i=0; i<Dim(); i++) {
            vals[i] = Ipopt::Min(scalar_, values_x[i]);
          }
        }
      }
      else {
        if (dense_x->homogeneous_) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = Ipopt::Min(values_[i], dense_x->scalar_);
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = Ipopt::Min(values_[i], values_x[i]);
          }
        }
      }
    }
  }

  void DenseVector::ElementWiseReciprocalImpl()
  {
    DBG_ASSERT(initialized_);
    if (homogeneous_) {
      scalar_ = 1.0/scalar_;
    }
    else {
      for (Index i=0; i<Dim(); i++) {
        values_[i] = 1.0/values_[i];
      }
    }
  }

  void DenseVector::ElementWiseAbsImpl()
  {
    DBG_ASSERT(initialized_);
    if (homogeneous_) {
      scalar_ = fabs(scalar_);
    }
    else {
      for (Index i=0; i<Dim(); i++) {
        values_[i] = fabs(values_[i]);
      }
    }
  }

  void DenseVector::ElementWiseSqrtImpl()
  {
    DBG_ASSERT(initialized_);
    if (homogeneous_) {
      scalar_ = sqrt(scalar_);
    }
    else {
      for (Index i=0; i<Dim(); i++) {
        values_[i] = sqrt(values_[i]);
      }
    }
  }

  void DenseVector::AddScalarImpl(Number scalar)
  {
    DBG_ASSERT(initialized_);
    if (homogeneous_) {
      scalar_ += scalar;
    }
    else {
      IpBlasDaxpy(Dim(), 1., &scalar, 0, values_, 1);
    }
  }

  Number DenseVector::MaxImpl() const
  {
    DBG_ASSERT(initialized_);
    DBG_ASSERT(Dim() > 0 && "There is no Max of a zero length vector (no reasonable default can be returned)");
    Number max;
    if (homogeneous_) {
      max = scalar_;
    }
    else {
      max = values_[0];
      for (Index i=1; i<Dim(); i++) {
        max = Ipopt::Max(values_[i], max);
      }
    }
    return max;
  }

  Number DenseVector::MinImpl() const
  {
    DBG_ASSERT(initialized_);
    DBG_ASSERT(Dim() > 0 && "There is no Min of a zero length vector"
               "(no reasonable default can be returned) - "
               "Check for zero length vector before calling");
    Number min;
    if (homogeneous_) {
      min = scalar_;
    }
    else {
      min = values_[0];
      for (Index i=1; i<Dim(); i++) {
        min = Ipopt::Min(values_[i], min);
      }
    }
    return min;
  }

  Number DenseVector::SumImpl() const
  {
    DBG_ASSERT(initialized_);
    Number sum;
    if (homogeneous_) {
      sum = Dim()*scalar_;
    }
    else {
      sum = 0.;
      for (Index i=0; i<Dim(); i++) {
        sum += values_[i];
      }
    }
    return sum;
  }

  Number DenseVector::SumLogsImpl() const
  {
    DBG_ASSERT(initialized_);
    Number sum;
    if (homogeneous_) {
      sum = Dim() * log(scalar_);
    }
    else {
      sum = 0.0;
      for (Index i=0; i<Dim(); i++) {
        sum += log(values_[i]);
      }
    }
    return sum;
  }

  void DenseVector::ElementWiseSgnImpl()
  {
    DBG_ASSERT(initialized_);
    if (homogeneous_) {
      if (scalar_ > 0.) {
        scalar_ = 1.;
      }
      else if (scalar_ < 0.) {
        scalar_ = -1.;
      }
      else {
        scalar_ = 0.;
      }
    }
    else {
      for (Index i=0; i<Dim(); i++) {
        if (values_[i] > 0.) {
          values_[i] = 1.;
        }
        else if (values_[i] < 0.) {
          values_[i] = -1.;
        }
        else {
          values_[i] = 0;
        }
      }
    }
  }

  // Specialized Functions
  void DenseVector::AddTwoVectorsImpl(Number a, const Vector& v1,
                                      Number b, const Vector& v2, Number c)
  {
    Number* values_v1=NULL;
    bool homogeneous_v1=false;
    Number scalar_v1 = 0;
    if (a!=0.) {
      const DenseVector* dense_v1 = dynamic_cast<const DenseVector*>(&v1);
      DBG_ASSERT(dense_v1);
      DBG_ASSERT(dense_v1->initialized_);
      DBG_ASSERT(Dim() == dense_v1->Dim());
      values_v1=dense_v1->values_;
      homogeneous_v1=dense_v1->homogeneous_;
      if (homogeneous_v1)
        scalar_v1 = dense_v1->scalar_;
    }
    Number* values_v2=NULL;
    bool homogeneous_v2=false;
    Number scalar_v2 = 0;
    if (b!=0.) {
      const DenseVector* dense_v2 = dynamic_cast<const DenseVector*>(&v2);
      DBG_ASSERT(dense_v2);
      DBG_ASSERT(dense_v2->initialized_);
      DBG_ASSERT(Dim() == dense_v2->Dim());
      values_v2=dense_v2->values_;
      homogeneous_v2=dense_v2->homogeneous_;
      if (homogeneous_v2)
        scalar_v2 = dense_v2->scalar_;
    }
    DBG_ASSERT(c==0. || initialized_);
    if ((c==0. || homogeneous_) && homogeneous_v1 && homogeneous_v2 ) {
      homogeneous_ = true;
      Number val = 0;
      if (c!=0.) {
        val = c*scalar_;
      }
      scalar_ = val + a*scalar_v1 + b*scalar_v2;
      initialized_ = true;
      return;
    }
    if (c==0.) {
      // make sure we have memory allocated for this vector
      values_allocated();
      homogeneous_ = false;
    }

    // If any of the vectors is homogeneous, call the default implementation
    if ( homogeneous_ || homogeneous_v1 || homogeneous_v2) {
      // ToDo:Should we implement specialized methods here too?
      Vector::AddTwoVectorsImpl(a, v1, b, v2, c);
      return;
    }

    // I guess I'm going over board here, but it might be best to
    // capture all cases for a, b, and c separately...
    if (c==0 ) {
      if (a==1.) {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v1[i];
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v1[i] + values_v2[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v1[i] - values_v2[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v1[i] + b*values_v2[i];
          }
        }
      }
      else if (a==-1.) {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v1[i];
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v1[i] + values_v2[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v1[i] - values_v2[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v1[i] + b*values_v2[i];
          }
        }
      }
      else if (a==0.) {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = 0.;
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v2[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v2[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = b*values_v2[i];
          }
        }
      }
      else {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = a*values_v1[i];
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = a*values_v1[i] + values_v2[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = a*values_v1[i] - values_v2[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = a*values_v1[i] + b*values_v2[i];
          }
        }
      }
    }
    else if (c==1.) {
      if (a==1.) {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] += values_v1[i];
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] += values_v1[i] + values_v2[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] += values_v1[i] - values_v2[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] += values_v1[i] + b*values_v2[i];
          }
        }
      }
      else if (a==-1.) {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] -= values_v1[i];
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] += -values_v1[i] + values_v2[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] += -values_v1[i] - values_v2[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] += -values_v1[i] + b*values_v2[i];
          }
        }
      }
      else if (a==0.) {
        if (b==0.) {
          /* Nothing */
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] += values_v2[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] -= values_v2[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] += b*values_v2[i];
          }
        }
      }
      else {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] += a*values_v1[i];
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] += a*values_v1[i] + values_v2[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] += a*values_v1[i] - values_v2[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] += a*values_v1[i] + b*values_v2[i];
          }
        }
      }
    }
    else if (c==-1.) {
      if (a==1.) {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v1[i] - values_[i];
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v1[i] + values_v2[i] - values_[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v1[i] - values_v2[i] - values_[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v1[i] + b*values_v2[i] - values_[i];
          }
        }
      }
      else if (a==-1.) {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v1[i] - values_[i];
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v1[i] + values_v2[i] - values_[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v1[i] - values_v2[i] - values_[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v1[i] + b*values_v2[i] - values_[i];
          }
        }
      }
      else if (a==0.) {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] *= -1.;
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v2[i] - values_[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v2[i] - values_[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = b*values_v2[i] - values_[i];
          }
        }
      }
      else {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = a*values_v1[i] - values_[i];
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = a*values_v1[i] + values_v2[i] - values_[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = a*values_v1[i] - values_v2[i] - values_[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = a*values_v1[i] + b*values_v2[i] - values_[i];
          }
        }
      }
    }
    else {
      if (a==1.) {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v1[i] + c*values_[i];
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v1[i] + values_v2[i] + c*values_[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v1[i] - values_v2[i] + c*values_[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v1[i] + b*values_v2[i] + c*values_[i];
          }
        }
      }
      else if (a==-1.) {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v1[i] + c*values_[i];
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v1[i] + values_v2[i] + c*values_[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v1[i] - values_v2[i] + c*values_[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v1[i] + b*values_v2[i] + c*values_[i];
          }
        }
      }
      else if (a==0.) {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] *= c;
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = values_v2[i] + c*values_[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = -values_v2[i] + c*values_[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = b*values_v2[i] + c*values_[i];
          }
        }
      }
      else {
        if (b==0.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = a*values_v1[i] + c*values_[i];
          }
        }
        else if (b==1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = a*values_v1[i] + values_v2[i] + c*values_[i];
          }
        }
        else if (b==-1.) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = a*values_v1[i] - values_v2[i] + c*values_[i];
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = a*values_v1[i] + b*values_v2[i] + c*values_[i];
          }
        }
      }
    }
    initialized_=true;
  }

  Number
  DenseVector::FracToBoundImpl(const Vector& delta, Number tau) const
  {
    DBG_ASSERT(Dim()==delta.Dim());
    DBG_ASSERT(tau>=0.);
    const DenseVector* dense_delta = dynamic_cast<const DenseVector*>(&delta);
    DBG_ASSERT(dense_delta);

    Number alpha = 1.;
    Number* values_x = values_;
    Number* values_delta = dense_delta->values_;
    if (homogeneous_) {
      if (dense_delta->homogeneous_) {
        if (dense_delta->scalar_<0.) {
          alpha = Ipopt::Min(alpha, -tau/dense_delta->scalar_ * scalar_);
        }
      }
      else {
        for (Index i=0; i<Dim(); i++) {
          if (values_delta[i]<0.) {
            alpha = Ipopt::Min(alpha, -tau/values_delta[i] * scalar_);
          }
        }
      }
    }
    else {
      if (dense_delta->homogeneous_) {
        if (dense_delta->scalar_<0.) {
          for (Index i=0; i<Dim(); i++) {
            alpha = Ipopt::Min(alpha, -tau/dense_delta->scalar_ * values_x[i]);
          }
        }
      }
      else {
        for (Index i=0; i<Dim(); i++) {
          if (values_delta[i]<0.) {
            alpha = Ipopt::Min(alpha, -tau/values_delta[i] * values_x[i]);
          }
        }
      }
    }

    DBG_ASSERT(alpha>=0.);
    return alpha;
  }

  void DenseVector::AddVectorQuotientImpl(Number a, const Vector& z,
                                          const Vector& s, Number c)
  {
    DBG_ASSERT(Dim()==z.Dim());
    DBG_ASSERT(Dim()==s.Dim());
    const DenseVector* dense_z = dynamic_cast<const DenseVector*>(&z);
    DBG_ASSERT(dense_z);
    DBG_ASSERT(dense_z->initialized_);
    const DenseVector* dense_s = dynamic_cast<const DenseVector*>(&s);
    DBG_ASSERT(dense_s);
    DBG_ASSERT(dense_s->initialized_);

    DBG_ASSERT(c==0. || initialized_);
    bool homogeneous_z = dense_z->homogeneous_;
    bool homogeneous_s = dense_s->homogeneous_;

    if ((c==0. || homogeneous_) && homogeneous_z && homogeneous_s) {
      if (c==0.) {
        scalar_ = a * dense_z->scalar_ / dense_s->scalar_;
      }
      else {
        scalar_ = c * scalar_ + a * dense_z->scalar_ / dense_s->scalar_;
      }
      initialized_ = true;
      homogeneous_ = true;
      if (values_) {
        owner_space_->FreeInternalStorage(values_);
        values_ = NULL;
      }
      return;
    }

    // At least one is not homogeneous
    // Make sure we have memory to store a non-homogeneous vector
    values_allocated();

    Number* values_z = dense_z->values_;
    Number* values_s = dense_s->values_;

    if (c==0.) {
      if (homogeneous_z) {
        // then s is not homogeneous
        for (Index i=0; i<Dim(); i++) {
          values_[i] = a * dense_z->scalar_ / values_s[i];
        }
      }
      else if (homogeneous_s) {
        // then z is not homogeneous
        for (Index i=0; i<Dim(); i++) {
          values_[i] = values_z[i] * a / dense_s->scalar_;
        }
      }
      else {
        for (Index i=0; i<Dim(); i++) {
          values_[i] = a * values_z[i] / values_s[i];
        }
      }
    }
    else if (homogeneous_) {
      Number val = c*scalar_;
      if (homogeneous_z) {
        // then s is not homogeneous
        for (Index i=0; i<Dim(); i++) {
          values_[i] = val + a * dense_z->scalar_ / values_s[i];
        }
      }
      else if (homogeneous_s) {
        // then z is not homogeneous
        for (Index i=0; i<Dim(); i++) {
          values_[i] = val + values_z[i] * a / dense_s->scalar_;
        }
      }
      else {
        for (Index i=0; i<Dim(); i++) {
          values_[i] = val + a * values_z[i] / values_s[i];
        }
      }
    }
    else {
      // ToDo could distinguish c = 1
      if (homogeneous_z) {
        if (homogeneous_s) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = c*values_[i] + a * dense_z->scalar_/dense_s->scalar_;
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = c*values_[i] + a * dense_z->scalar_/values_s[i];
          }
        }
      }
      else {
        if (homogeneous_s) {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = c*values_[i] + values_z[i] * a /dense_s->scalar_;
          }
        }
        else {
          for (Index i=0; i<Dim(); i++) {
            values_[i] = c*values_[i] + a * values_z[i]/values_s[i];
          }
        }
      }
    }

    initialized_ = true;
    homogeneous_ = false;
  }

  void DenseVector::CopyToPos(Index Pos, const Vector& x)
  {
    Index dim_x = x.Dim();
    DBG_ASSERT(dim_x+Pos<=Dim());
    const DenseVector* dense_x = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dense_x);

    Number* vals = values_allocated();
    homogeneous_ = false;

    if (dense_x->homogeneous_) {
      IpBlasDcopy(dim_x, &scalar_, 1, vals+Pos, 1);
    }
    else {
      IpBlasDcopy(dim_x, dense_x->values_, 1, vals+Pos, 1);
    }
    initialized_ = true;
    ObjectChanged();
  }

  void DenseVector::CopyFromPos(Index Pos, Vector& x) const
  {
    Index dim_x = x.Dim();
    DBG_ASSERT(dim_x+Pos<=Dim());
    DenseVector* dense_x = dynamic_cast<DenseVector*>(&x);
    DBG_ASSERT(dense_x);
    DBG_ASSERT(dense_x->homogeneous_); // This might have to made more general

    if (homogeneous_) {
      IpBlasDcopy(dim_x, &scalar_, 1, dense_x->values_, 1);
    }
    else {
      IpBlasDcopy(dim_x, values_+Pos, 1, dense_x->values_, 1);
    }
    // We need to tell X that it has changed!
    dense_x->ObjectChanged();
    dense_x->initialized_=true;
  }

  void DenseVector::PrintImpl(const Journalist& jnlst,
                              EJournalLevel level,
                              EJournalCategory category,
                              const std::string& name,
                              Index indent,
                              const std::string& prefix) const
  {
    jnlst.PrintfIndented(level, category, indent,
                         "%sDenseVector \"%s\" with %d elements:\n",
                         prefix.c_str(), name.c_str(), Dim());
    if (initialized_) {
      if (homogeneous_) {
        jnlst.PrintfIndented(level, category, indent,
                             "%sHomogeneous vector, all elements have value %23.16e\n",
                             prefix.c_str(), scalar_);
      }
      else {
        for (Index i=0; i<Dim(); i++) {
          jnlst.PrintfIndented(level, category, indent,
                               "%s%s[%5d]=%23.16e\n",
                               prefix.c_str(), name.c_str(), i+1, values_[i]);
        }
      }
    }
    else {
      jnlst.PrintfIndented(level, category, indent,
                           "%sUninitialized!\n",
                           prefix.c_str());
    }
  }
} // namespace Ipopt
