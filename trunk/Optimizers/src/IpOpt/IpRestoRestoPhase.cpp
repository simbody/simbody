// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpRestoRestoPhase.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-02-11

#include "IpRestoRestoPhase.hpp"
#include "IpRestoIpoptNLP.hpp"

namespace Ipopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  RestoRestorationPhase::RestoRestorationPhase()
  {}

  RestoRestorationPhase::~RestoRestorationPhase()
  {}

  bool RestoRestorationPhase::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    return true;
  }

  bool
  RestoRestorationPhase::PerformRestoration()
  {
    DBG_START_METH("RestoRestorationPhase::PerformRestoration",
                   dbg_verbosity);
    Jnlst().Printf(J_DETAILED, J_MAIN,
                   "Performing second level restoration phase for current constriant violation %8.2e\n", IpCq().curr_constraint_violation());

    DBG_ASSERT(IpCq().curr_constraint_violation()>0.);

    // Get a grip on the restoration phase NLP and obtain the pointers
    // to the original NLP data
    SmartPtr<RestoIpoptNLP> resto_ip_nlp =
      dynamic_cast<RestoIpoptNLP*> (&IpNLP());
    DBG_ASSERT(IsValid(resto_ip_nlp));
    SmartPtr<IpoptNLP> orig_ip_nlp =
      dynamic_cast<IpoptNLP*> (&resto_ip_nlp->OrigIpNLP());
    DBG_ASSERT(IsValid(orig_ip_nlp));

    // Get the current point and create a new vector for the result
    SmartPtr<const CompoundVector> Ccurr_x =
      dynamic_cast<const CompoundVector*> (GetRawPtr(IpData().curr()->x()));
    SmartPtr<Vector> new_x = IpData().curr()->x()->MakeNew();
    SmartPtr<CompoundVector> Cnew_x =
      dynamic_cast<CompoundVector*> (GetRawPtr(new_x));

    // The x values remain unchanged
    SmartPtr<Vector> x = Cnew_x->GetCompNonConst(0);
    x->Copy(*Ccurr_x->GetComp(0));

    // ToDo in free mu mode - what to do here?
    Number mu = IpData().curr_mu();

    // Compute the initial values for the n and p variables for the
    // equality constraints
    Number rho = resto_ip_nlp->Rho();
    SmartPtr<Vector> nc = Cnew_x->GetCompNonConst(1);
    SmartPtr<Vector> pc = Cnew_x->GetCompNonConst(2);
    SmartPtr<const Vector> cvec = orig_ip_nlp->c(*Ccurr_x->GetComp(0));
    SmartPtr<Vector> a = nc->MakeNew();
    SmartPtr<Vector> b = nc->MakeNew();
    a->Set(mu/(2.*rho));
    a->Axpy(-0.5, *cvec);
    b->Copy(*cvec);
    b->Scal(mu/(2.*rho));
    solve_quadratic(*a, *b, *nc);
    pc->Copy(*cvec);
    pc->Axpy(1., *nc);
    DBG_PRINT_VECTOR(2, "nc", *nc);
    DBG_PRINT_VECTOR(2, "pc", *pc);

    // initial values for the n and p variables for the inequality
    // constraints
    SmartPtr<Vector> nd = Cnew_x->GetCompNonConst(3);
    SmartPtr<Vector> pd = Cnew_x->GetCompNonConst(4);
    SmartPtr<Vector> dvec = pd->MakeNew();
    dvec->Copy(*orig_ip_nlp->d(*Ccurr_x->GetComp(0)));
    dvec->Axpy(-1., *IpData().curr()->s());
    a = nd->MakeNew();
    b = nd->MakeNew();
    a->Set(mu/(2.*rho));
    a->Axpy(-0.5, *dvec);
    b->Copy(*dvec);
    b->Scal(mu/(2.*rho));
    solve_quadratic(*a, *b, *nd);
    pd->Copy(*dvec);
    pd->Axpy(1., *nd);
    DBG_PRINT_VECTOR(2, "nd", *nd);
    DBG_PRINT_VECTOR(2, "pd", *pd);

    // Now set the trial point to the solution of the restoration phase
    // s and all multipliers remain unchanged
    SmartPtr<IteratesVector> new_trial = IpData().curr()->MakeNewContainer();
    new_trial->Set_x(*new_x);
    IpData().set_trial(new_trial);

    IpData().Append_info_string("R");

    return true;
  }

  void RestoRestorationPhase::solve_quadratic(const Vector& a,
      const Vector& b,
      Vector& v)
  {
    v.Copy(a);
    v.ElementWiseMultiply(a);

    v.Axpy(1., b);
    v.ElementWiseSqrt();

    v.Axpy(1., a);
  }

} // namespace Ipopt
