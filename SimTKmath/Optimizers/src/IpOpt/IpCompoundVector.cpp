// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpCompoundVector.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpoptConfig.h"
#include "IpCompoundVector.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

// Keeps MS VC++ 8 quiet about sprintf, strcpy, etc.
#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif





#ifdef HAVE_CSTDIO
# include <cstdio>
#else
# ifdef HAVE_STDIO_H
#  include <stdio.h>
# else
#  error "don't have header file for stdio"
# endif
#endif

#include <limits>

namespace SimTKIpopt
{

#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  CompoundVector::CompoundVector(const CompoundVectorSpace* owner_space, bool create_new)
      :
      Vector(owner_space),
      comps_(owner_space->NCompSpaces()),
      const_comps_(owner_space->NCompSpaces()),
      owner_space_(owner_space),
      vectors_valid_(false)
  {
    Index dim_check = 0;
    for (Index i=0; i<NComps(); i++) {
      SmartPtr<const VectorSpace> space = owner_space_->GetCompSpace(i);
      DBG_ASSERT(IsValid(space));
      dim_check += space->Dim();

      if (create_new) {
        comps_[i] = space->MakeNew();
      }
    }

    DBG_ASSERT(dim_check == Dim());

    if (create_new) {
      vectors_valid_ = VectorsValid();
    }
  }

  CompoundVector::~CompoundVector()
  {
    // ToDo: Do we need an empty here?
  }

  void CompoundVector::SetComp(Index icomp, const Vector& vec)
  {
    DBG_ASSERT(icomp<NComps());
    comps_[icomp] = NULL;
    const_comps_[icomp] = &vec;

    vectors_valid_ = VectorsValid();
    ObjectChanged();
  }

  void CompoundVector::SetCompNonConst(Index icomp, Vector& vec)
  {
    DBG_ASSERT(icomp < NComps());
    comps_[icomp] = &vec;
    const_comps_[icomp] = NULL;

    vectors_valid_ = VectorsValid();
    ObjectChanged();
  }

  void CompoundVector::CopyImpl(const Vector& x)
  {
    DBG_START_METH("CompoundVector::CopyImpl(const Vector& x)", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    const CompoundVector* comp_x = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(comp_x);
    DBG_ASSERT(NComps() == comp_x->NComps());
    for(Index i=0; i<NComps(); i++) {
      Comp(i)->Copy(*comp_x->GetComp(i));
    }
  }

  void CompoundVector::ScalImpl(Number alpha)
  {
    DBG_START_METH("CompoundVector::ScalImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    for(Index i=0; i<NComps(); i++) {
      DBG_ASSERT(Comp(i));
      Comp(i)->Scal(alpha);
    }
  }

  void CompoundVector::AxpyImpl(Number alpha, const Vector &x)
  {
    DBG_START_METH("CompoundVector::AxpyImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    const CompoundVector* comp_x = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(comp_x);
    DBG_ASSERT(NComps() == comp_x->NComps());
    for(Index i=0; i<NComps(); i++) {
      DBG_ASSERT(Comp(i));
      Comp(i)->Axpy(alpha, *comp_x->GetComp(i));
    }
  }

  Number CompoundVector::DotImpl(const Vector &x) const
  {
    DBG_START_METH("CompoundVector::DotImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    const CompoundVector* comp_x = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(comp_x);
    DBG_ASSERT(NComps() == comp_x->NComps());
    Number dot = 0.;
    for(Index i=0; i<NComps(); i++) {
      DBG_ASSERT(ConstComp(i));
      dot += ConstComp(i)->Dot(*comp_x->GetComp(i));
    }
    return dot;
  }

  Number CompoundVector::Nrm2Impl() const
  {
    DBG_START_METH("CompoundVector::Nrm2Impl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    Number sum=0.;
    for(Index i=0; i<NComps(); i++) {
      Number nrm2 = ConstComp(i)->Nrm2();
      sum += nrm2*nrm2;
    }
    return sqrt(sum);
  }

  Number CompoundVector::AsumImpl() const
  {
    DBG_START_METH("CompoundVector::AsumImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    Number sum=0.;
    for(Index i=0; i<NComps(); i++) {
      sum += ConstComp(i)->Asum();
    }
    return sum;
  }

  Number CompoundVector::AmaxImpl() const
  {
    DBG_START_METH("CompoundVector::AmaxImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    Number max=0.;
    for(Index i=0; i<NComps(); i++) {
      max = SimTKIpopt::Max(max, ConstComp(i)->Amax());
    }
    return max;
  }

  void CompoundVector::SetImpl(Number value)
  {
    DBG_START_METH("CompoundVector::SetImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    for(Index i=0; i<NComps(); i++) {
      Comp(i)->Set(value);
    }
  }

  void CompoundVector::ElementWiseDivideImpl(const Vector& x)
  {
    DBG_START_METH("CompoundVector::ElementWiseDivideImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    const CompoundVector* comp_x = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(comp_x);
    DBG_ASSERT(NComps() == comp_x->NComps());
    for(Index i=0; i<NComps(); i++) {
      Comp(i)->ElementWiseDivide(*comp_x->GetComp(i));
    }
  }

  void CompoundVector::ElementWiseMultiplyImpl(const Vector& x)
  {
    DBG_START_METH("CompoundVector::ElementWiseMultiplyImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    const CompoundVector* comp_x = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(comp_x);
    DBG_ASSERT(NComps() == comp_x->NComps());
    for(Index i=0; i<NComps(); i++) {
      Comp(i)->ElementWiseMultiply(*comp_x->GetComp(i));
    }
  }

  void CompoundVector::ElementWiseMaxImpl(const Vector& x)
  {
    DBG_START_METH("CompoundVector::ElementWiseMaxImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    const CompoundVector* comp_x = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(comp_x);
    DBG_ASSERT(NComps() == comp_x->NComps());
    for(Index i=0; i<NComps(); i++) {
      Comp(i)->ElementWiseMax(*comp_x->GetComp(i));
    }
  }

  void CompoundVector::ElementWiseMinImpl(const Vector& x)
  {
    DBG_START_METH("CompoundVector::ElementWiseMinImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    const CompoundVector* comp_x = dynamic_cast<const CompoundVector*>(&x);
    DBG_ASSERT(comp_x);
    DBG_ASSERT(NComps() == comp_x->NComps());
    for(Index i=0; i<NComps(); i++) {
      Comp(i)->ElementWiseMin(*comp_x->GetComp(i));
    }
  }

  void CompoundVector::ElementWiseReciprocalImpl()
  {
    DBG_START_METH("CompoundVector::ElementWiseReciprocalImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    for(Index i=0; i<NComps(); i++) {
      Comp(i)->ElementWiseReciprocal();
    }
  }

  void CompoundVector::ElementWiseAbsImpl()
  {
    DBG_START_METH("CompoundVector::ElementWiseAbsImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    for(Index i=0; i<NComps(); i++) {
      Comp(i)->ElementWiseAbs();
    }
  }

  void CompoundVector::ElementWiseSqrtImpl()
  {
    DBG_START_METH("CompoundVector::ElementWiseSqrtImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    for(Index i=0; i<NComps(); i++) {
      Comp(i)->ElementWiseSqrt();
    }
  }

  void CompoundVector::AddScalarImpl(Number scalar)
  {
    DBG_START_METH("CompoundVector::AddScalarImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    for(Index i=0; i<NComps(); i++) {
      Comp(i)->AddScalar(scalar);
    }
  }

  Number CompoundVector::MaxImpl() const
  {
    DBG_START_METH("CompoundVector::MaxImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    DBG_ASSERT(NComps() > 0 && Dim() > 0 && "There is no Max of a zero length vector (no reasonable default can be returned)");
    Number max = -std::numeric_limits<Number>::max();
    for(Index i=0; i<NComps(); i++) {
      if (ConstComp(i)->Dim() != 0) {
        max = SimTKIpopt::Max(max, ConstComp(i)->Max());
      }
    }
    return max;
  }

  Number CompoundVector::MinImpl() const
  {
    DBG_START_METH("CompoundVector::MinImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    DBG_ASSERT(NComps() > 0 && Dim() > 0 && "There is no Min of a zero length vector (no reasonable default can be returned)");
    Number min = std::numeric_limits<Number>::max();
    for (Index i=0; i<NComps(); i++) {
      if (ConstComp(i)->Dim() != 0) {
        min = SimTKIpopt::Min(min, ConstComp(i)->Min());
      }
    }
    return min;
  }

  Number CompoundVector::SumImpl() const
  {
    DBG_START_METH("CompoundVector::SumImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    Number sum=0.;
    for(Index i=0; i<NComps(); i++) {
      sum += ConstComp(i)->Sum();
    }
    return sum;
  }

  Number CompoundVector::SumLogsImpl() const
  {
    DBG_START_METH("CompoundVector::SumLogsImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    Number sum=0.;
    for(Index i=0; i<NComps(); i++) {
      sum += ConstComp(i)->SumLogs();
    }
    return sum;
  }

  void CompoundVector::ElementWiseSgnImpl()
  {
    DBG_START_METH("CompoundVector::ElementWiseSgnImpl", dbg_verbosity);
    DBG_ASSERT(vectors_valid_);
    for(Index i=0; i<NComps(); i++) {
      Comp(i)->ElementWiseSgn();
    }
  }

  // Specialized Functions
  void CompoundVector::AddTwoVectorsImpl(Number a, const Vector& v1,
                                         Number b, const Vector& v2, Number c)
  {
    DBG_ASSERT(vectors_valid_);
    const CompoundVector* comp_v1 = dynamic_cast<const CompoundVector*>(&v1);
    DBG_ASSERT(comp_v1);
    DBG_ASSERT(NComps() == comp_v1->NComps());
    const CompoundVector* comp_v2 = dynamic_cast<const CompoundVector*>(&v2);
    DBG_ASSERT(comp_v2);
    DBG_ASSERT(NComps() == comp_v2->NComps());

    for(Index i=0; i<NComps(); i++) {
      Comp(i)->AddTwoVectors(a, *comp_v1->GetComp(i), b, *comp_v2->GetComp(i), c);
    }
  }

  Number
  CompoundVector::FracToBoundImpl(const Vector& delta, Number tau) const
  {
    DBG_ASSERT(vectors_valid_);
    const CompoundVector* comp_delta =
      dynamic_cast<const CompoundVector*>(&delta);
    DBG_ASSERT(comp_delta);
    DBG_ASSERT(NComps() == comp_delta->NComps());

    Number alpha = 1.;
    for(Index i=0; i<NComps(); i++) {
      alpha = SimTKIpopt::Min(alpha,
                         ConstComp(i)->FracToBound(*comp_delta->GetComp(i), tau));
    }
    return alpha;
  }

  void CompoundVector::AddVectorQuotientImpl(Number a, const Vector& z,
      const Vector& s, Number c)
  {
    DBG_ASSERT(vectors_valid_);
    const CompoundVector* comp_z =
      dynamic_cast<const CompoundVector*>(&z);
    DBG_ASSERT(comp_z);
    DBG_ASSERT(NComps() == comp_z->NComps());
    const CompoundVector* comp_s =
      dynamic_cast<const CompoundVector*>(&s);
    DBG_ASSERT(comp_s);
    DBG_ASSERT(NComps() == comp_s->NComps());

    for(Index i=0; i<NComps(); i++) {
      Comp(i)->AddVectorQuotient(a, *comp_z->GetComp(i),
                                 *comp_s->GetComp(i), c);
    }
  }

  bool CompoundVector::HasValidNumbersImpl() const
  {
    DBG_ASSERT(vectors_valid_);
    for(Index i=0; i<NComps(); i++) {
      if (!ConstComp(i)->HasValidNumbers()) {
        return false;
      }
    }
    return true;
  }

  void CompoundVector::PrintImpl(const Journalist& jnlst,
                                 EJournalLevel level,
                                 EJournalCategory category,
                                 const std::string& name,
                                 Index indent,
                                 const std::string& prefix) const
  {
    DBG_START_METH("CompoundVector::PrintImpl", dbg_verbosity);
    jnlst.Printf(level, category, "\n");
    jnlst.PrintfIndented(level, category, indent,
                         "%sCompoundVector \"%s\" with %d components:\n",
                         prefix.c_str(), name.c_str(), NComps());
    for(Index i=0; i<NComps(); i++) {
      jnlst.Printf(level, category, "\n");
      jnlst.PrintfIndented(level, category, indent,
                           "%sComponent %d:\n", prefix.c_str(), i+1);
      if (ConstComp(i)) {
        DBG_ASSERT(name.size()<200);
        char buffer[256];
        sprintf(buffer, "%s[%2d]", name.c_str(), i);
        std::string term_name = buffer;
        ConstComp(i)->Print(&jnlst, level, category, term_name,
                            indent+1, prefix);
      }
      else {
        jnlst.PrintfIndented(level, category, indent,
                             "%sComponent %d is not yet set!\n",
                             prefix.c_str(), i+1);
      }
    }
  }

  bool CompoundVector::VectorsValid()
  {
    bool retVal = true;
    for (Index i=0; i<NComps(); i++) {
      // Better not have an entry in both (sanity check)
      DBG_ASSERT(IsNull(comps_[i]) || IsNull(const_comps_[i]));
      if (IsNull(comps_[i])  && IsNull(const_comps_[i])) {
        retVal = false;
        break;
      }
    }
    return retVal;
  }

  CompoundVectorSpace::CompoundVectorSpace(Index ncomp_spaces, Index total_dim)
      :
      VectorSpace(total_dim),
      ncomp_spaces_(ncomp_spaces),
      comp_spaces_(ncomp_spaces)
  {}

  void CompoundVectorSpace::SetCompSpace(Index icomp, const VectorSpace& vec_space)
  {
    DBG_ASSERT(icomp<ncomp_spaces_);
    DBG_ASSERT(IsNull(comp_spaces_[icomp]));
    comp_spaces_[icomp] = &vec_space;
  }

  SmartPtr<const VectorSpace> CompoundVectorSpace::GetCompSpace(Index icomp) const
  {
    DBG_ASSERT(icomp<ncomp_spaces_);
    return comp_spaces_[icomp];
  }

}
