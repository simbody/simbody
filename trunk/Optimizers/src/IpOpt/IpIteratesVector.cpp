// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpIteratesVector.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-06-06

#include "IpIteratesVector.hpp"

namespace Ipopt
{

  IteratesVector::IteratesVector(const IteratesVectorSpace* owner_space, bool create_new)
      :
      CompoundVector(owner_space, create_new),
      owner_space_(owner_space)
  {
    DBG_ASSERT(owner_space_);
  }

  IteratesVector::~IteratesVector()
  {}

  SmartPtr<IteratesVector> IteratesVector::MakeNewIteratesVector(bool create_new) const
  {
    return owner_space_->MakeNewIteratesVector(create_new);
  }

  SmartPtr<IteratesVector> IteratesVector::MakeNewContainer() const
  {
    SmartPtr<IteratesVector> ret = MakeNewIteratesVector(false);

    if (IsValid(x())) {
      ret->Set_x(*x());
    }
    if (IsValid(s())) {
      ret->Set_s(*s());
    }
    if (IsValid(y_c())) {
      ret->Set_y_c(*y_c());
    }
    if (IsValid(y_d())) {
      ret->Set_y_d(*y_d());
    }
    if (IsValid(z_L())) {
      ret->Set_z_L(*z_L());
    }
    if (IsValid(z_U())) {
      ret->Set_z_U(*z_U());
    }
    if (IsValid(v_L())) {
      ret->Set_v_L(*v_L());
    }
    if (IsValid(v_U())) {
      ret->Set_v_U(*v_U());
    }

    return ret;

    // We may need a non const version
    //     if (IsCompConst(0)) {
    //       ret->Set_x(*x());
    //     }
    //     else {
    //       ret->Set_x_NonConst(*x_NonConst());
    //     }

    //     if (IsCompConst(1)) {
    //       ret->Set_s(*s());
    //     }
    //     else {
    //       ret->Set_s_NonConst(*s_NonConst());
    //     }

    //     if (IsCompConst(2)) {
    //       ret->Set_y_c(*y_c());
    //     }
    //     else {
    //       ret->Set_y_c_NonConst(*y_c_NonConst());
    //     }

    //     if (IsCompConst(3)) {
    //       ret->Set_y_d(*y_d());
    //     }
    //     else {
    //       ret->Set_y_d_NonConst(*y_d_NonConst());
    //     }

    //     if (IsCompConst(4)) {
    //       ret->Set_z_L(*z_L());
    //     }
    //     else {
    //       ret->Set_z_L_NonConst(*z_L_NonConst());
    //     }

    //     if (IsCompConst(5)) {
    //       ret->Set_z_U(*z_U());
    //     }
    //     else {
    //       ret->Set_z_U_NonConst(*z_U_NonConst());
    //     }

    //     if (IsCompConst(6)) {
    //       ret->Set_v_L(*v_L());
    //     }
    //     else {
    //       ret->Set_v_L_NonConst(*v_L_NonConst());
    //     }

    //     if (IsCompConst(7)) {
    //       ret->Set_v_U(*v_U());
    //     }
    //     else {
    //       ret->Set_v_U_NonConst(*v_U_NonConst());
    //     }

    //    return ret;
  }

  IteratesVectorSpace::IteratesVectorSpace(const VectorSpace& x_space, const VectorSpace& s_space,
      const VectorSpace& y_c_space, const VectorSpace& y_d_space,
      const VectorSpace& z_L_space, const VectorSpace& z_U_space,
      const VectorSpace& v_L_space, const VectorSpace& v_U_space
                                          )
      :
      CompoundVectorSpace(8, x_space.Dim() + s_space.Dim()
                          + y_c_space.Dim() + y_d_space.Dim()
                          + z_L_space.Dim() + z_U_space.Dim()
                          + v_L_space.Dim() + v_U_space.Dim()
                         )
  {
    x_space_ = &x_space;
    s_space_ = &s_space;
    y_c_space_ = &y_c_space;
    y_d_space_ = &y_d_space;
    z_L_space_ = &z_L_space;
    z_U_space_ = &z_U_space;
    v_L_space_ = &v_L_space;
    v_U_space_ = &v_U_space;

    this->CompoundVectorSpace::SetCompSpace(0, *x_space_);
    this->CompoundVectorSpace::SetCompSpace(1, *s_space_);
    this->CompoundVectorSpace::SetCompSpace(2, *y_c_space_);
    this->CompoundVectorSpace::SetCompSpace(3, *y_d_space_);
    this->CompoundVectorSpace::SetCompSpace(4, *z_L_space_);
    this->CompoundVectorSpace::SetCompSpace(5, *z_U_space_);
    this->CompoundVectorSpace::SetCompSpace(6, *v_L_space_);
    this->CompoundVectorSpace::SetCompSpace(7, *v_U_space_);
  }

  IteratesVectorSpace::~IteratesVectorSpace()
  {}

} // namespae Ipopt
