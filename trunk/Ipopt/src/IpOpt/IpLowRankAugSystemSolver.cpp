// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpLowRankAugSystemSolver.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter                IBM    2005-12-27

#include "IpLowRankAugSystemSolver.hpp"
#include "IpLowRankUpdateSymMatrix.hpp"

namespace Ipopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  LowRankAugSystemSolver::LowRankAugSystemSolver(
    AugSystemSolver& aug_system_solver)
      :
      AugSystemSolver(),
      aug_system_solver_(&aug_system_solver),
      w_tag_(0),
      w_factor_(0.),
      d_x_tag_(0),
      delta_x_(0.),
      d_s_tag_(0),
      delta_s_(0.),
      j_c_tag_(0),
      d_c_tag_(0),
      delta_c_(0.),
      j_d_tag_(0),
      d_d_tag_(0),
      delta_d_(0.)
  {
    DBG_START_METH("LowRankAugSystemSolver::LowRankAugSystemSolver()",dbg_verbosity);
    DBG_ASSERT(IsValid(aug_system_solver_));
  }

  LowRankAugSystemSolver::~LowRankAugSystemSolver()
  {
    DBG_START_METH("LowRankAugSystemSolver::~LowRankAugSystemSolver()",dbg_verbosity);
  }

  bool LowRankAugSystemSolver::InitializeImpl(const OptionsList& options,
      const std::string& prefix)
  {
    first_call_ = true;
    J1_ = NULL;
    J2_ = NULL;
    Vtilde1_ = NULL;
    Utilde2_ = NULL;
    Wdiag_ = NULL;
    compound_sol_vecspace_ = NULL;

    return aug_system_solver_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                          options, prefix);
  }

  ESymSolverStatus LowRankAugSystemSolver::Solve(
    const SymMatrix* W,
    double W_factor,
    const Vector* D_x,
    double delta_x,
    const Vector* D_s,
    double delta_s,
    const Matrix* J_c,
    const Vector* D_c,
    double delta_c,
    const Matrix* J_d,
    const Vector* D_d,
    double delta_d,
    const Vector& rhs_x,
    const Vector& rhs_s,
    const Vector& rhs_c,
    const Vector& rhs_d,
    Vector& sol_x,
    Vector& sol_s,
    Vector& sol_c,
    Vector& sol_d,
    bool check_NegEVals,
    Index numberOfNegEVals)
  {
    DBG_START_METH("LowRankAugSystemSolver::Solve",dbg_verbosity);

    ESymSolverStatus retval;

    if (first_call_) {
      DBG_ASSERT(IsNull(Wdiag_));
      // Set up the diagonal matrix Wdiag_
      Index dimx = rhs_x.Dim();
      SmartPtr<DiagMatrixSpace> Wdiag_space = new DiagMatrixSpace(dimx);
      Wdiag_ = Wdiag_space->MakeNewDiagMatrix();
    }

    if (first_call_ ||
        AugmentedSystemRequiresChange(W, W_factor, D_x, delta_x, D_s, delta_s,
                                      *J_c, D_c, delta_c, *J_d, D_d,
                                      delta_d) ) {
      retval = UpdateFactorization(W, W_factor, D_x, delta_x, D_s, delta_s,
                                   *J_c, D_c, delta_c, *J_d, D_d, delta_d,
                                   rhs_x, rhs_s, rhs_c, rhs_d,
                                   check_NegEVals, numberOfNegEVals);
      if (retval != SYMSOLVER_SUCCESS) {
        return retval;
      }

      // Store the tags
      w_tag_ = W->GetTag();
      w_factor_ = W_factor;
      if (D_x) {
        d_x_tag_ = D_x->GetTag();
      }
      else {
        d_x_tag_ = 0;
      }
      delta_x_ = delta_x;
      if (D_s) {
        d_s_tag_ = D_s->GetTag();
      }
      else {
        d_s_tag_ = 0;
      }
      delta_s_ = delta_s;
      if (J_c) {
        j_c_tag_ = J_c->GetTag();
      }
      else {
        j_c_tag_ = 0;
      }
      if (D_c) {
        d_c_tag_ = D_c->GetTag();
      }
      else {
        d_c_tag_ = 0;
      }
      delta_c_ = delta_c;
      if (J_d) {
        j_d_tag_ = J_d->GetTag();
      }
      else {
        j_d_tag_ = 0;
      }
      if (D_d) {
        d_d_tag_ = D_d->GetTag();
      }
      else {
        d_d_tag_ = 0;
      }
      delta_d_ = delta_d;

      first_call_ = false;
    }

    // Now solve the system for the given right hand side, using the
    // Sherman-Morrison formula with factorization information already
    // computed.
    retval = aug_system_solver_->Solve(GetRawPtr(Wdiag_), W_factor,
                                       D_x, delta_x, D_s, delta_s,
                                       J_c, D_c, delta_c, J_d, D_d, delta_d,
                                       rhs_x, rhs_s, rhs_c, rhs_d,
                                       sol_x, sol_s, sol_c, sol_d,
                                       check_NegEVals, numberOfNegEVals);
    if (aug_system_solver_->ProvidesInertia()) {
      num_neg_evals_ = aug_system_solver_->NumberOfNegEVals();
    }
    if (retval != SYMSOLVER_SUCCESS) {
      Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                     "LowRankAugSystemSolver: AugSystemSolver returned retval = %d for right hand side.\n", retval);
      return retval;
    }

    if (IsValid(Vtilde1_) || IsValid(Utilde2_)) {
      // Create a CompoundVectors to store the right hand side and
      // solutions
      SmartPtr<CompoundVector> crhs =
        compound_sol_vecspace_->MakeNewCompoundVector(false);
      crhs->SetComp(0, rhs_x);
      crhs->SetComp(1, rhs_s);
      crhs->SetComp(2, rhs_c);
      crhs->SetComp(3, rhs_d);
      SmartPtr<CompoundVector> csol =
        compound_sol_vecspace_->MakeNewCompoundVector(false);
      csol->SetCompNonConst(0, sol_x);
      csol->SetCompNonConst(1, sol_s);
      csol->SetCompNonConst(2, sol_c);
      csol->SetCompNonConst(3, sol_d);

      if (IsValid(Utilde2_)) {
        Index nU = Utilde2_->NCols();
        SmartPtr<DenseVectorSpace> bUspace =
          new DenseVectorSpace(nU);
        SmartPtr<DenseVector> bU = bUspace->MakeNewDenseVector();
        Utilde2_->TransMultVector(1., *crhs, 0., *bU);
        J2_->CholeskySolveVector(*bU);
        Utilde2_->MultVector(1., *bU, 1., *csol);
      }
      if (IsValid(Vtilde1_)) {
        Index nV = Vtilde1_->NCols();
        SmartPtr<DenseVectorSpace> bVspace =
          new DenseVectorSpace(nV);
        SmartPtr<DenseVector> bV = bVspace->MakeNewDenseVector();
        Vtilde1_->TransMultVector(1., *crhs, 0., *bV);
        J1_->CholeskySolveVector(*bV);
        Vtilde1_->MultVector(-1., *bV, 1., *csol);
      }
    }

    return retval;
  }

  ESymSolverStatus LowRankAugSystemSolver::UpdateFactorization(
    const SymMatrix* W,
    double W_factor,
    const Vector* D_x,
    double delta_x,
    const Vector* D_s,
    double delta_s,
    const Matrix& J_c,
    const Vector* D_c,
    double delta_c,
    const Matrix& J_d,
    const Vector* D_d,
    double delta_d,
    const Vector& proto_rhs_x,
    const Vector& proto_rhs_s,
    const Vector& proto_rhs_c,
    const Vector& proto_rhs_d,
    bool check_NegEVals,
    Index numberOfNegEVals)
  {
    DBG_START_METH("LowRankAugSystemSolver::UpdateFactorization",
                   dbg_verbosity);

    DBG_ASSERT(W_factor == 0.0 || W_factor == 1.0);
    ESymSolverStatus retval = SYMSOLVER_SUCCESS;

    // Get the low update information out of W
    const LowRankUpdateSymMatrix* LR_W =
      dynamic_cast<const LowRankUpdateSymMatrix*> (W);
    DBG_ASSERT(LR_W);
    DBG_PRINT_MATRIX(2, "LR_W", *LR_W);

    SmartPtr<const Vector> B0;
    SmartPtr<const MultiVectorMatrix> V;
    SmartPtr<const MultiVectorMatrix> U;
    if (W_factor == 1.0) {
      V = LR_W->GetV();
      U = LR_W->GetU();
      B0 = LR_W->GetDiag();
    }
    SmartPtr<const Matrix> P_LM = LR_W->P_LowRank();
    SmartPtr<const VectorSpace> LR_VecSpace = LR_W->LowRankVectorSpace();

    if (IsNull(B0)) {
      SmartPtr<Vector> zero_B0 = (IsValid(P_LM)) ? LR_VecSpace->MakeNew() : proto_rhs_x.MakeNew();
      zero_B0->Set(0.0);
      B0 = GetRawPtr(zero_B0);
    }

    // set up the Hessian for the underlying augmented system solver
    // without the low-rank update
    if (IsValid(P_LM) && LR_W->ReducedDiag()) {
      DBG_ASSERT(IsValid(B0));
      SmartPtr<Vector> fullx = proto_rhs_x.MakeNew();
      P_LM->MultVector(1., *B0, 0., *fullx);
      Wdiag_->SetDiag(*fullx);
    }
    else {
      Wdiag_->SetDiag(*B0);
    }

    SmartPtr<MultiVectorMatrix> Vtilde1_x;
    if (IsValid(V)) {
      SmartPtr<MultiVectorMatrix> V_x;
      Index nV = V->NCols();
      //DBG_PRINT((1, "delta_x  = %e\n", delta_x));
      //DBG_PRINT_MATRIX(2, "V", *V);
      retval = SolveMultiVector(D_x, delta_x, D_s, delta_s, J_c,
                                D_c, delta_c, J_d, D_d, delta_d,
                                proto_rhs_x, proto_rhs_s, proto_rhs_c,
                                proto_rhs_d, *V, P_LM, V_x, Vtilde1_,
                                Vtilde1_x, check_NegEVals, numberOfNegEVals);
      if (retval != SYMSOLVER_SUCCESS) {
        Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                       "LowRankAugSystemSolver: SolveMultiVector returned retval = %d for V.\n", retval);
        return retval;
      }
      //DBG_PRINT_MATRIX(2, "Vtilde1_x", *Vtilde1_x);

      SmartPtr<DenseSymMatrixSpace> M1space =
        new DenseSymMatrixSpace(nV);
      SmartPtr<DenseSymMatrix> M1 = M1space->MakeNewDenseSymMatrix();
      M1->FillIdentity();
      M1->HighRankUpdateTranspose(1., *Vtilde1_x, *V_x, 1.);
      //DBG_PRINT_MATRIX(2, "M1", *M1);
      SmartPtr<DenseGenMatrixSpace> J1space =
        new DenseGenMatrixSpace(nV, nV);
      J1_ = J1space->MakeNewDenseGenMatrix();
      bool retchol = J1_->ComputeCholeskyFactor(*M1);
      // M1 must be positive definite!
      //DBG_ASSERT(retchol);
      if (!retchol) {
        Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                       "LowRankAugSystemSolver: Cholesky for M1 returned error!\n");
        retval = SYMSOLVER_WRONG_INERTIA;
        num_neg_evals_++;
        return retval;
      }
    }
    else {
      Vtilde1_ = NULL;
      J1_ = NULL;
    }

    if (IsValid(U)) {
      Index nU = U->NCols();
      SmartPtr<MultiVectorMatrix> U_x;
      SmartPtr<MultiVectorMatrix> Utilde1;
      SmartPtr<MultiVectorMatrix> Utilde1_x;
      SmartPtr<MultiVectorMatrix> Utilde2_x;
      retval = SolveMultiVector(D_x, delta_x, D_s, delta_s, J_c,
                                D_c, delta_c, J_d, D_d, delta_d,
                                proto_rhs_x, proto_rhs_s, proto_rhs_c,
                                proto_rhs_d, *U, P_LM, U_x, Utilde1,
                                Utilde1_x, check_NegEVals, numberOfNegEVals);
      if (retval != SYMSOLVER_SUCCESS) {
        Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                       "LowRankAugSystemSolver: SolveMultiVector returned retval = %d for U.\n", retval);
        return retval;
      }

      if (IsNull(Vtilde1_)) {
        Utilde2_ = Utilde1;
        Utilde2_x = Utilde1_x;
      }
      else {
        Index nV = Vtilde1_->NCols();
        SmartPtr<DenseGenMatrixSpace> Cspace =
          new DenseGenMatrixSpace(nV, nU);
        SmartPtr<DenseGenMatrix> C = Cspace->MakeNewDenseGenMatrix();
        C->HighRankUpdateTranspose(1., *Vtilde1_x, *U_x, 0.);
        J1_->CholeskySolveMatrix(*C);
        Utilde2_ = Utilde1;
        Utilde2_->AddRightMultMatrix(-1, *Vtilde1_, *C, 1.);
        Utilde2_x = Utilde1_x->MakeNewMultiVectorMatrix();
        for (Index i=0; i<Utilde1_x->NCols(); i++) {
          const CompoundVector* cvec =
            dynamic_cast<const CompoundVector*> (GetRawPtr(Utilde2_->GetVector(i)));
          DBG_ASSERT(cvec);
          Utilde2_x->SetVector(i, *cvec->GetComp(0));
        }
      }

      SmartPtr<DenseSymMatrixSpace> M2space =
        new DenseSymMatrixSpace(nU);
      SmartPtr<DenseSymMatrix> M2 = M2space->MakeNewDenseSymMatrix();
      M2->FillIdentity();
      M2->HighRankUpdateTranspose(-1., *Utilde2_x, *U_x, 1.);
      SmartPtr<DenseGenMatrixSpace> J2space =
        new DenseGenMatrixSpace(nU, nU);
      J2_ = J2space->MakeNewDenseGenMatrix();
      //DBG_PRINT_MATRIX(2, "M2", *M2);
      bool retchol = J2_->ComputeCholeskyFactor(*M2);
      if (!retchol) {
        Jnlst().Printf(J_DETAILED, J_SOLVE_PD_SYSTEM,
                       "LowRankAugSystemSolver: Cholesky for M2 returned error.\n");
        retval = SYMSOLVER_WRONG_INERTIA;
        num_neg_evals_++;
        return retval;
      }
    }
    else {
      J2_ = NULL;
      Utilde2_ = NULL;
    }

    return retval;
  }

  ESymSolverStatus LowRankAugSystemSolver::SolveMultiVector(
    const Vector* D_x,
    double delta_x,
    const Vector* D_s,
    double delta_s,
    const Matrix& J_c,
    const Vector* D_c,
    double delta_c,
    const Matrix& J_d,
    const Vector* D_d,
    double delta_d,
    const Vector& proto_rhs_x,
    const Vector& proto_rhs_s,
    const Vector& proto_rhs_c,
    const Vector& proto_rhs_d,
    const MultiVectorMatrix& V,
    const SmartPtr<const Matrix>& P_LM,
    SmartPtr<MultiVectorMatrix>& V_x,
    SmartPtr<MultiVectorMatrix>& Vtilde,
    SmartPtr<MultiVectorMatrix>& Vtilde_x,
    bool check_NegEVals,
    Index numberOfNegEVals)
  {
    DBG_START_METH("LowRankAugSystemSolver::SolveMultiVector",
                   dbg_verbosity);

    ESymSolverStatus retval;

    Index nrhs = V.NCols();
    DBG_ASSERT(nrhs>0);

    SmartPtr<MultiVectorMatrixSpace> V_xspace =
      new MultiVectorMatrixSpace(nrhs, *proto_rhs_x.OwnerSpace());
    V_x = V_xspace->MakeNewMultiVectorMatrix();

    // Create the right hand sides
    std::vector<SmartPtr<const Vector> > rhs_xV(nrhs);
    std::vector<SmartPtr<const Vector> > rhs_sV(nrhs);
    std::vector<SmartPtr<const Vector> > rhs_cV(nrhs);
    std::vector<SmartPtr<const Vector> > rhs_dV(nrhs);

    for (Index i=0; i<nrhs; i++) {
      if (IsNull(P_LM)) {
        rhs_xV[i] = V.GetVector(i);
        DBG_ASSERT(rhs_xV[i]->Dim() == proto_rhs_x.Dim());
      }
      else {
        SmartPtr<Vector> fullx = proto_rhs_x.MakeNew();
        P_LM->MultVector(1., *V.GetVector(i), 0., *fullx);
        rhs_xV[i] = ConstPtr(fullx);
      }
      V_x->SetVector(i, *rhs_xV[i]);
      SmartPtr<Vector> tmp;
      tmp =  proto_rhs_s.MakeNew();
      tmp->Set(0.);
      rhs_sV[i] = ConstPtr(tmp);
      tmp =  proto_rhs_c.MakeNew();
      tmp->Set(0.);
      rhs_cV[i] = ConstPtr(tmp);
      tmp =  proto_rhs_d.MakeNew();
      tmp->Set(0.);
      rhs_dV[i] = ConstPtr(tmp);
    }

    // now get space for the solution
    std::vector<SmartPtr<Vector> > sol_xV(nrhs);
    std::vector<SmartPtr<Vector> > sol_sV(nrhs);
    std::vector<SmartPtr<Vector> > sol_cV(nrhs);
    std::vector<SmartPtr<Vector> > sol_dV(nrhs);
    for (Index i=0; i<nrhs; i++) {
      sol_xV[i] = proto_rhs_x.MakeNew();
      sol_sV[i] = proto_rhs_s.MakeNew();
      sol_cV[i] = proto_rhs_c.MakeNew();
      sol_dV[i] = proto_rhs_d.MakeNew();
    }

    // Call the actual augmented system solver to obtain Vtilde
    retval = aug_system_solver_->MultiSolve(GetRawPtr(Wdiag_), 1.0, D_x, delta_x, D_s, delta_s,
                                            &J_c, D_c, delta_c, &J_d, D_d, delta_d,
                                            rhs_xV, rhs_sV, rhs_cV, rhs_dV,
                                            sol_xV, sol_sV, sol_cV, sol_dV,
                                            check_NegEVals, numberOfNegEVals);

    if (aug_system_solver_->ProvidesInertia()) {
      num_neg_evals_ = aug_system_solver_->NumberOfNegEVals();
    }
    if (retval != SYMSOLVER_SUCCESS) {
      return retval;
    }

    // Pack the results into Vtilde
    if (IsNull(compound_sol_vecspace_)) {
      Index dimx = proto_rhs_x.Dim();
      Index dims = proto_rhs_s.Dim();
      Index dimc = proto_rhs_c.Dim();
      Index dimd = proto_rhs_d.Dim();
      Index dimtot = dimx+dims+dimc+dimd;
      SmartPtr<CompoundVectorSpace> vecspace =
        new CompoundVectorSpace(4, dimtot);
      vecspace->SetCompSpace(0, *proto_rhs_x.OwnerSpace());
      vecspace->SetCompSpace(1, *proto_rhs_s.OwnerSpace());
      vecspace->SetCompSpace(2, *proto_rhs_c.OwnerSpace());
      vecspace->SetCompSpace(3, *proto_rhs_d.OwnerSpace());
      compound_sol_vecspace_ = ConstPtr(vecspace);
    }
    SmartPtr<MultiVectorMatrixSpace> V1space =
      new MultiVectorMatrixSpace(nrhs, *compound_sol_vecspace_);
    Vtilde = V1space->MakeNewMultiVectorMatrix();
    Vtilde_x = V_xspace->MakeNewMultiVectorMatrix();
    for (Index i=0; i<nrhs; i++) {
      Vtilde_x->SetVector(i, *sol_xV[i]);
      SmartPtr<CompoundVector> cvec =
        compound_sol_vecspace_->MakeNewCompoundVector(false);
      cvec->SetCompNonConst(0, *sol_xV[i]);
      cvec->SetCompNonConst(1, *sol_sV[i]);
      cvec->SetCompNonConst(2, *sol_cV[i]);
      cvec->SetCompNonConst(3, *sol_dV[i]);
      Vtilde->SetVectorNonConst(i, *cvec);
    }

    return retval;
  }

  bool LowRankAugSystemSolver::AugmentedSystemRequiresChange(
    const SymMatrix* W,
    double W_factor,
    const Vector* D_x,
    double delta_x,
    const Vector* D_s,
    double delta_s,
    const Matrix& J_c,
    const Vector* D_c,
    double delta_c,
    const Matrix& J_d,
    const Vector* D_d,
    double delta_d)
  {
    DBG_START_METH("LowRankAugSystemSolver::AugmentedSystemRequiresChange",
                   dbg_verbosity);

#ifdef IP_DEBUG

    bool Wtest = (W && W->GetTag() != w_tag_);
    bool iWtest = (!W && w_tag_ != 0);
    bool wfactor_test = (W_factor != w_factor_);
    bool D_xtest = (D_x && D_x->GetTag() != d_x_tag_);
    bool iD_xtest = (!D_x && d_x_tag_ != 0);
    bool delta_xtest = (delta_x != delta_x_);
    bool D_stest = (D_s && D_s->GetTag() != d_s_tag_);
    bool iD_stest = (!D_s && d_s_tag_ != 0);
    bool delta_stest = (delta_s != delta_s_);
    bool J_ctest = (J_c.GetTag() != j_c_tag_);
    bool D_ctest = (D_c && D_c->GetTag() != d_c_tag_);
    bool iD_ctest = (!D_c && d_c_tag_ != 0);
    bool delta_ctest = (delta_c != delta_c_);
    bool J_dtest = (J_d.GetTag() != j_d_tag_);
    bool D_dtest = (D_d && D_d->GetTag() != d_d_tag_);
    bool iD_dtest = (!D_d && d_d_tag_ != 0);
    bool delta_dtest = (delta_d != delta_d_);
#endif

    DBG_PRINT((2,"Wtest = %d\n", Wtest));
    DBG_PRINT((2,"iWtest = %d\n", iWtest));
    DBG_PRINT((2,"wfactor_test = %d\n", wfactor_test));
    DBG_PRINT((2,"D_xtest = %d\n", D_xtest));
    DBG_PRINT((2,"iD_xtest = %d\n", iD_xtest));
    DBG_PRINT((2,"delta_xtest = %d\n", delta_xtest));
    DBG_PRINT((2,"D_stest = %d\n", D_stest));
    DBG_PRINT((2,"iD_stest = %d\n", iD_stest));
    DBG_PRINT((2,"delta_stest = %d\n", delta_stest));
    DBG_PRINT((2,"J_ctest = %d\n", J_ctest));
    DBG_PRINT((2,"D_ctest = %d\n", D_ctest));
    DBG_PRINT((2,"iD_ctest = %d\n", iD_ctest));
    DBG_PRINT((2,"delta_ctest = %d\n", delta_ctest));
    DBG_PRINT((2,"J_dtest = %d\n", J_dtest));
    DBG_PRINT((2,"D_dtest = %d\n", D_dtest));
    DBG_PRINT((2,"iD_dtest = %d\n", iD_dtest));
    DBG_PRINT((2,"delta_dtest = %d\n", delta_dtest));

    if ( (W && W->GetTag() != w_tag_)
         || (!W && w_tag_ != 0)
         || (W_factor != w_factor_)
         || (D_x && D_x->GetTag() != d_x_tag_)
         || (!D_x && d_x_tag_ != 0)
         || (delta_x != delta_x_)
         || (D_s && D_s->GetTag() != d_s_tag_)
         || (!D_s && d_s_tag_ != 0)
         || (delta_s != delta_s_)
         || (J_c.GetTag() != j_c_tag_)
         || (D_c && D_c->GetTag() != d_c_tag_)
         || (!D_c && d_c_tag_ != 0)
         || (delta_c != delta_c_)
         || (J_d.GetTag() != j_d_tag_)
         || (D_d && D_d->GetTag() != d_d_tag_)
         || (!D_d && d_d_tag_ != 0)
         || (delta_d != delta_d_) ) {
      return true;
    }

    return false;
  }

  Index LowRankAugSystemSolver::NumberOfNegEVals() const
  {
    DBG_ASSERT(!first_call_);
    return num_neg_evals_;
  }

  bool LowRankAugSystemSolver::ProvidesInertia() const
  {
    return aug_system_solver_->ProvidesInertia();
  }

  bool LowRankAugSystemSolver::IncreaseQuality()
  {
    return aug_system_solver_->IncreaseQuality();
  }

} // namespace Ipopt
