// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpStdAugSystemSolver.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpStdAugSystemSolver.hpp"
#include "IpDebug.hpp"

#include "IpCompoundSymMatrix.hpp"
#include "IpCompoundVector.hpp"
#include "IpSumSymMatrix.hpp"
#include "IpDiagMatrix.hpp"
#include "IpIdentityMatrix.hpp"
// ToDo: Remove below here - for debug only
#include "IpTripletHelper.hpp"
// ToDo: Remove above here

namespace Ipopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  StdAugSystemSolver::StdAugSystemSolver(SymLinearSolver& linSolver)
      :
      AugSystemSolver(),
      linsolver_(&linSolver),
      augmented_system_space_(NULL),
      sumsym_space_x_(NULL),
      diag_space_x_(NULL),
      diag_space_s_(NULL),
      diag_space_c_(NULL),
      ident_space_ds_(NULL),
      diag_space_d_(NULL),
      w_tag_(0),
      d_x_tag_(0),
      delta_x_(0.),
      d_s_tag_(0),
      delta_s_(0.),
      j_c_tag_(0),
      d_c_tag_(0),
      delta_c_(0.),
      j_d_tag_(0),
      d_d_tag_(0),
      delta_d_(0.),
      old_w_(NULL)
  {
    DBG_START_METH("StdAugSystemSolver::StdAugSystemSolver()",dbg_verbosity);
    DBG_ASSERT(IsValid(linsolver_));
  }

  StdAugSystemSolver::~StdAugSystemSolver()
  {
    DBG_START_METH("StdAugSystemSolver::~StdAugSystemSolver()",dbg_verbosity);
  }


  bool StdAugSystemSolver::InitializeImpl(const OptionsList& options,
                                          const std::string& prefix)
  {
    // This option is registered by OrigIpoptNLP
    options.GetBoolValue("warm_start_same_structure",
                         warm_start_same_structure_, prefix);

    if (!warm_start_same_structure_) {
      augsys_tag_ = 0;
      augmented_system_ = NULL;
    }
    else {
      ASSERT_EXCEPTION(IsValid(augmented_system_), INVALID_WARMSTART,
                       "StdAugSystemSolver called with warm_start_same_structure, but augmented system is not initialized.");
    }

    return linsolver_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                  options, prefix);
  }

  ESymSolverStatus StdAugSystemSolver::MultiSolve(
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
    std::vector<SmartPtr<const Vector> >& rhs_xV,
    std::vector<SmartPtr<const Vector> >& rhs_sV,
    std::vector<SmartPtr<const Vector> >& rhs_cV,
    std::vector<SmartPtr<const Vector> >& rhs_dV,
    std::vector<SmartPtr<Vector> >& sol_xV,
    std::vector<SmartPtr<Vector> >& sol_sV,
    std::vector<SmartPtr<Vector> >& sol_cV,
    std::vector<SmartPtr<Vector> >& sol_dV,
    bool check_NegEVals,
    Index numberOfNegEVals)
  {
    DBG_START_METH("StdAugSystemSolver::MultiSolve",dbg_verbosity);
    DBG_ASSERT(J_c && J_d && "Currently, you MUST specify J_c and J_d in the augmented system");

    DBG_ASSERT(W_factor == 0.0 || W_factor == 1.0);

    Index nrhs = (Index)rhs_xV.size();
    DBG_ASSERT(nrhs>0);
    DBG_ASSERT(nrhs==(Index)rhs_sV.size());
    DBG_ASSERT(nrhs==(Index)rhs_cV.size());
    DBG_ASSERT(nrhs==(Index)rhs_dV.size());
    DBG_ASSERT(nrhs==(Index)sol_xV.size());
    DBG_ASSERT(nrhs==(Index)sol_sV.size());
    DBG_ASSERT(nrhs==(Index)sol_cV.size());
    DBG_ASSERT(nrhs==(Index)sol_dV.size());

    // Create the compound matrix of the augmented system if it has not
    // yet been created - It is assumed that the structure will not change
    // after this call
    DBG_DO(bool debug_first_time_through = false;)
    if (!IsValid(augmented_system_)) {
      // pass in the information to form the structure of the augmented system
      // rhs_? are passed in to provide a prototype vector
      // for D_? (since these may be NULL)
      DBG_ASSERT(W && J_c && J_d); // W must exist during the first call to setup the structure!
      CreateAugmentedSpace(*W, *J_c, *J_d, *rhs_xV[0], *rhs_sV[0],
                           *rhs_cV[0], *rhs_dV[0]);
      CreateAugmentedSystem(W, W_factor, D_x, delta_x, D_s, delta_s,
                            *J_c, D_c, delta_c, *J_d, D_d, delta_d,
                            *rhs_xV[0], *rhs_sV[0], *rhs_cV[0], *rhs_dV[0]);
      DBG_DO(debug_first_time_through = true;)
    }


    // Check if anything that was just passed in is different from
    // what is currently in the compound matrix of the augmented
    // system. If anything is different, then update the augmented
    // system
    if ( AugmentedSystemRequiresChange(W, W_factor, D_x, delta_x, D_s, delta_s,
                                       *J_c, D_c, delta_c, *J_d, D_d, delta_d) ) {
      DBG_ASSERT(!debug_first_time_through);
      CreateAugmentedSystem(W, W_factor, D_x, delta_x, D_s, delta_s,
                            *J_c, D_c, delta_c, *J_d, D_d, delta_d,
                            *rhs_xV[0], *rhs_sV[0], *rhs_cV[0], *rhs_dV[0]);
    }

    // Sanity checks
    DBG_ASSERT(rhs_xV[0]->Dim() == sol_xV[0]->Dim());
    DBG_ASSERT(rhs_sV[0]->Dim() == sol_sV[0]->Dim());
    DBG_ASSERT(rhs_cV[0]->Dim() == sol_cV[0]->Dim());
    DBG_ASSERT(rhs_dV[0]->Dim() == sol_dV[0]->Dim());

    // Now construct the overall right hand side vector that will be passed
    // to the linear solver
    std::vector<SmartPtr<const Vector> > augmented_rhsV(nrhs);
    for (Index i=0; i<nrhs; i++) {
      SmartPtr<CompoundVector> augrhs =
        augmented_vector_space_->MakeNewCompoundVector();
      augrhs->SetComp(0, *rhs_xV[i]);
      augrhs->SetComp(1, *rhs_sV[i]);
      augrhs->SetComp(2, *rhs_cV[i]);
      augrhs->SetComp(3, *rhs_dV[i]);
      char buffer[16];
      sprintf(buffer, "RHS[%2d]", i);
      augrhs->Print(Jnlst(), J_MOREVECTOR, J_LINEAR_ALGEBRA, buffer);
      augmented_rhsV[i] = GetRawPtr(augrhs);
    }

    augmented_system_->Print(Jnlst(), J_MATRIX, J_LINEAR_ALGEBRA, "KKT");
    if (Jnlst().ProduceOutput(J_MOREMATRIX, J_LINEAR_ALGEBRA)) {
      // ToDo: remove below here - for debug only
      Index dbg_nz = TripletHelper::GetNumberEntries(*augmented_system_);
      Index* dbg_iRows = new Index[dbg_nz];
      Index* dbg_jCols = new Index[dbg_nz];
      Number* dbg_values = new Number[dbg_nz];
      TripletHelper::FillRowCol(dbg_nz, *augmented_system_, dbg_iRows, dbg_jCols);
      TripletHelper::FillValues(dbg_nz, *augmented_system_, dbg_values);
      Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA, "******* KKT SYSTEM *******\n");
      for (Index dbg_i=0; dbg_i<dbg_nz; dbg_i++) {
        Jnlst().Printf(J_MOREMATRIX, J_LINEAR_ALGEBRA, "(%d) KKT[%d][%d] = %23.15e\n", dbg_i, dbg_iRows[dbg_i], dbg_jCols[dbg_i], dbg_values[dbg_i]);
      }
      delete [] dbg_iRows;
      dbg_iRows = NULL;
      delete [] dbg_jCols;
      dbg_jCols = NULL;
      delete [] dbg_values;
      dbg_values = NULL;
      // ToDo: remove above here
    }

    // Call the linear solver
    std::vector<SmartPtr<Vector> > augmented_solV(nrhs);
    for (Index i=0; i<nrhs; i++) {
      SmartPtr<CompoundVector> augsol =
        augmented_vector_space_->MakeNewCompoundVector();
      augsol->SetCompNonConst(0, *sol_xV[i]);
      augsol->SetCompNonConst(1, *sol_sV[i]);
      augsol->SetCompNonConst(2, *sol_cV[i]);
      augsol->SetCompNonConst(3, *sol_dV[i]);
      augmented_solV[i] = GetRawPtr(augsol);
    }
    ESymSolverStatus retval;
    retval = linsolver_->MultiSolve(*augmented_system_, augmented_rhsV,
                                    augmented_solV, check_NegEVals,
                                    numberOfNegEVals);
    if (retval==SYMSOLVER_SUCCESS) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Factorization successful.\n");
      for (Index i=0; i<nrhs; i++) {
        char buffer[16];
        sprintf(buffer, "SOL[%2d]", i);
        augmented_solV[i]->Print(Jnlst(), J_MOREVECTOR, J_LINEAR_ALGEBRA,
                                 buffer);
      }
    }
    else if (retval==SYMSOLVER_FATAL_ERROR) {
      THROW_EXCEPTION(FATAL_ERROR_IN_LINEAR_SOLVER,"A fatal error occured in the linear solver.");
    }
    else {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA, "Factorization failed with retval = %d\n", retval);
    }

    return retval;
  }

  void StdAugSystemSolver::CreateAugmentedSpace(const SymMatrix& W,
      const Matrix& J_c,
      const Matrix& J_d,
      const Vector& proto_x,
      const Vector& proto_s,
      const Vector& proto_c,
      const Vector& proto_d)
  {
    DBG_ASSERT(!IsValid(augmented_system_));

    old_w_ = &W;

    //===
    // Setup the augmented system matrix (described in IpAugSystemSolver.hpp")
    //===

    // created the compound symmetric matrix space
    Index n_x = J_c.NCols();
    Index n_s = J_d.NRows();
    Index n_c = J_c.NRows();
    Index n_d = n_s;

    Index total_nRows = n_x + n_s + n_c + n_d;
    augmented_system_space_ = new CompoundSymMatrixSpace(4, total_nRows);
    augmented_system_space_->SetBlockDim(0, n_x);
    augmented_system_space_->SetBlockDim(1, n_s);
    augmented_system_space_->SetBlockDim(2, n_c);
    augmented_system_space_->SetBlockDim(3, n_d);

    // (1,1) block
    // create the spaces and sum matrix for the upper left corner (W + D_x delta_x*I)
    // of the hessian part for the 1,1 block
    diag_space_x_ = new DiagMatrixSpace(n_x);

    sumsym_space_x_ = new SumSymMatrixSpace(n_x, 2);
    sumsym_space_x_->SetTermSpace(0, *W.OwnerSymMatrixSpace());
    sumsym_space_x_->SetTermSpace(1, *diag_space_x_);
    augmented_system_space_->SetCompSpace(0,0, *sumsym_space_x_);


    // (2,2) block
    // create the spaces and diag matrix for the lower right corner (D_s + delta_s*I)
    // of the hessian part, the 2,2 block
    diag_space_s_ = new DiagMatrixSpace(n_s);
    augmented_system_space_->SetCompSpace(1,1, *diag_space_s_);

    // (3,1) block
    augmented_system_space_->SetCompSpace(2,0, *J_c.OwnerSpace());

    // (3,3) block
    // create the matrix space and matrix for the 3,3 block
    diag_space_c_ = new DiagMatrixSpace(n_c);
    augmented_system_space_->SetCompSpace(2,2, *diag_space_c_);

    // (4,1) block
    augmented_system_space_->SetCompSpace(3,0, *J_d.OwnerSpace());

    // (4,2) block
    // create the identity matrix space and matrix for the 4,2 block
    ident_space_ds_ = new IdentityMatrixSpace(n_s);
    augmented_system_space_->SetCompSpace(3,1, *ident_space_ds_);

    // (4,4) block
    // create the sum matrix space and matrix for the 4,4 block
    diag_space_d_ = new DiagMatrixSpace(n_d);
    augmented_system_space_->SetCompSpace(3,3, *diag_space_d_);

    // Create the space for the vectors
    augmented_vector_space_ = new CompoundVectorSpace(4, n_x + n_s + n_c + n_d);
    augmented_vector_space_->SetCompSpace(0, *proto_x.OwnerSpace());
    augmented_vector_space_->SetCompSpace(1, *proto_s.OwnerSpace());
    augmented_vector_space_->SetCompSpace(2, *proto_c.OwnerSpace());
    augmented_vector_space_->SetCompSpace(3, *proto_d.OwnerSpace());

  }

  void StdAugSystemSolver::CreateAugmentedSystem(
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
    const Vector& proto_x,
    const Vector& proto_s,
    const Vector& proto_c,
    const Vector& proto_d)
  {
    augmented_system_ = augmented_system_space_->MakeNewCompoundSymMatrix();

    // (1,1) block
    SmartPtr<SumSymMatrix> sumsym_x = sumsym_space_x_->MakeNewSumSymMatrix();

    if (W) {
      sumsym_x->SetTerm(0, W_factor, *W);
      old_w_ = W;
      w_tag_ = W->GetTag();
    }
    else {
      sumsym_x->SetTerm(0, 0.0, *old_w_);
      w_tag_ = 0;
    }
    w_factor_ = W_factor;

    SmartPtr<DiagMatrix> diag_x = diag_space_x_->MakeNewDiagMatrix();
    if (D_x) {
      if (delta_x==0.) {
        diag_x->SetDiag(*D_x);
      }
      else {
        SmartPtr<Vector> tmp = D_x->MakeNewCopy();
        tmp->AddScalar(delta_x);
        diag_x->SetDiag(*tmp);
      }
      d_x_tag_ = D_x->GetTag();
    }
    else {
      SmartPtr<Vector> tmp = proto_x.MakeNew();
      tmp->Set(delta_x);
      diag_x->SetDiag(*tmp);
      d_x_tag_ = 0;
    }
    sumsym_x->SetTerm(1, 1.0, *diag_x);
    delta_x_ = delta_x;

    augmented_system_->SetComp(0,0, *sumsym_x);

    // (2,2) block
    SmartPtr<DiagMatrix> diag_s = diag_space_s_->MakeNewDiagMatrix();
    if (D_s) {
      if (delta_s==0.) {
        diag_s->SetDiag(*D_s);
      }
      else {
        SmartPtr<Vector> tmp = D_s->MakeNewCopy();
        tmp->AddScalar(delta_s);
        diag_s->SetDiag(*tmp);
      }
      d_s_tag_ = D_s->GetTag();
    }
    else {
      SmartPtr<Vector> tmp = proto_s.MakeNew();
      tmp->Set(delta_s);
      diag_s->SetDiag(*tmp);
      d_s_tag_ = 0;
    }
    delta_s_ = delta_s;

    augmented_system_->SetComp(1, 1, *diag_s);

    // (3,1) block
    augmented_system_->SetComp(2,0, J_c);
    j_c_tag_ = J_c.GetTag();

    // (3,3) block
    SmartPtr<DiagMatrix> diag_c = diag_space_c_->MakeNewDiagMatrix();
    if (D_c) {
      if (delta_c==0.) {
        diag_c->SetDiag(*D_c);
      }
      else {
        SmartPtr<Vector> tmp = D_c->MakeNewCopy();
        tmp->AddScalar(-delta_c);
        diag_c->SetDiag(*tmp);
      }
      d_c_tag_ = D_c->GetTag();
    }
    else {
      SmartPtr<Vector> tmp = proto_c.MakeNew();
      tmp->Set(-delta_c);
      diag_c->SetDiag(*tmp);
      d_c_tag_ = 0;
    }
    delta_c_ = delta_c;

    augmented_system_->SetComp(2,2, *diag_c);

    // (4,1) block
    augmented_system_->SetComp(3,0, J_d);
    j_d_tag_ = J_d.GetTag();

    // (4,2) block
    SmartPtr<IdentityMatrix> ident_ds = ident_space_ds_->MakeNewIdentityMatrix();
    ident_ds->SetFactor(-1.0);
    augmented_system_->SetComp(3,1, *ident_ds);

    // (4,4) block
    SmartPtr<DiagMatrix> diag_d = diag_space_d_->MakeNewDiagMatrix();
    if (D_d) {
      if (delta_d==0.) {
        diag_d->SetDiag(*D_d);
      }
      else {
        SmartPtr<Vector> tmp = D_d->MakeNewCopy();
        tmp->AddScalar(-delta_d);
        diag_d->SetDiag(*tmp);
      }
      d_d_tag_ = D_d->GetTag();
    }
    else {
      SmartPtr<Vector> tmp = proto_d.MakeNew();
      tmp->Set(-delta_d);
      diag_d->SetDiag(*tmp);
      d_d_tag_ = 0;
    }
    delta_d_ = delta_d;

    augmented_system_->SetComp(3,3, *diag_d);

    augsys_tag_ = augmented_system_->GetTag();
  }


  bool StdAugSystemSolver::AugmentedSystemRequiresChange(
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
    DBG_START_METH("StdAugSystemSolver::AugmentedSystemRequiresChange",dbg_verbosity);
    DBG_ASSERT(augsys_tag_ == augmented_system_->GetTag() && "Someone has changed the augmented system outside of the AugSystemSolver. This should NOT happen.");

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

  Index StdAugSystemSolver::NumberOfNegEVals() const
  {
    DBG_ASSERT(IsValid(augmented_system_));
    return linsolver_->NumberOfNegEVals();
  }

  bool StdAugSystemSolver::ProvidesInertia() const
  {
    return linsolver_->ProvidesInertia();
  }

  bool StdAugSystemSolver::IncreaseQuality()
  {
    return linsolver_->IncreaseQuality();
  }

} // namespace Ipopt
