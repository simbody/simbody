// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpLeastSquareMults.cpp 765 2006-07-14 18:03:23Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter              IBM    2004-09-23

#include "IpLeastSquareMults.hpp"

namespace SimTKIpopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  LeastSquareMultipliers::LeastSquareMultipliers(AugSystemSolver& augSysSolver)
      :
      EqMultiplierCalculator(),
      augsyssolver_(&augSysSolver)
  {}

  bool LeastSquareMultipliers::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    return augsyssolver_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                     options, prefix);
  }

  bool LeastSquareMultipliers::CalculateMultipliers
  (Vector& y_c,
   Vector& y_d)
  {
    DBG_START_METH("LeastSquareMultipliers::CalculateMultipliers",
                   dbg_verbosity);

    SmartPtr<const SymMatrix> zeroW = IpNLP().uninitialized_h();
    DBG_PRINT_MATRIX(2, "zeroW", *zeroW);
    SmartPtr<const Matrix> J_c = IpCq().curr_jac_c();
    SmartPtr<const Matrix> J_d = IpCq().curr_jac_d();
    SmartPtr<const Vector> grad_f = IpCq().curr_grad_f();
    SmartPtr<const Matrix> Px_L = IpNLP().Px_L();
    SmartPtr<const Matrix> Px_U = IpNLP().Px_U();
    SmartPtr<const Matrix> Pd_L = IpNLP().Pd_L();
    SmartPtr<const Matrix> Pd_U = IpNLP().Pd_U();
    SmartPtr<const Vector> z_L = IpData().curr()->z_L();
    SmartPtr<const Vector> z_U = IpData().curr()->z_U();
    SmartPtr<const Vector> v_L = IpData().curr()->v_L();
    SmartPtr<const Vector> v_U = IpData().curr()->v_U();

    // Compute the right hand side
    SmartPtr<Vector> rhs_x = grad_f->MakeNew();
    rhs_x->Copy(*grad_f);
    DBG_PRINT_VECTOR(2, "rhs_x grad_f", *rhs_x);
    Px_L->MultVector(1., *z_L, -1., *rhs_x);
    Px_U->MultVector(-1., *z_U, 1., *rhs_x);

    SmartPtr<Vector> rhs_s = IpData().curr()->s()->MakeNew();
    Pd_L->MultVector(1., *v_L, 0., *rhs_s);
    Pd_U->MultVector(-1., *v_U, 1., *rhs_s);

    SmartPtr<Vector> rhs_c = y_c.MakeNew();
    rhs_c->Set(0.);
    SmartPtr<Vector> rhs_d = y_d.MakeNew();
    rhs_d->Set(0.);

    SmartPtr<Vector> sol_x = rhs_x->MakeNew();
    SmartPtr<Vector> sol_s = rhs_s->MakeNew();

    DBG_PRINT_VECTOR(2, "rhs_x", *rhs_x);
    DBG_PRINT_VECTOR(2, "rhs_s", *rhs_s);
    DBG_PRINT_VECTOR(2, "rhs_c", *rhs_c);
    DBG_PRINT_VECTOR(2, "rhs_d", *rhs_d);

    enum ESymSolverStatus retval;
    Index numberOfEVals=rhs_c->Dim()+rhs_d->Dim();
    retval = augsyssolver_->Solve(GetRawPtr(zeroW), 0.0, NULL, 1.0, NULL,
                                  1.0, GetRawPtr(J_c), NULL, 0.,
                                  GetRawPtr(J_d), NULL, 0., *rhs_x, *rhs_s,
                                  *rhs_c, *rhs_d, *sol_x, *sol_s,
                                  y_c, y_d, true, numberOfEVals);
    if (retval!=SYMSOLVER_SUCCESS) {
      return false;
    }

    DBG_PRINT_VECTOR(2, "sol_x", *sol_x);
    DBG_PRINT_VECTOR(2, "sol_s", *sol_s);
    DBG_PRINT_VECTOR(2, "sol_c", y_c);
    DBG_PRINT_VECTOR(2, "sol_d", y_d);

    return true;
  }

} // namespace Ipopt
