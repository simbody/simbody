// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpTSymLinearSolver.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-03-17

#include "IpTSymLinearSolver.hpp"
#include "IpTripletHelper.hpp"

namespace SimTKIpopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  TSymLinearSolver::TSymLinearSolver
  (SmartPtr<SparseSymLinearSolverInterface> solver_interface,
   SmartPtr<TSymScalingMethod> scaling_method)
      :
      SymLinearSolver(),
      atag_(0),
      dim_(0),
      nonzeros_triplet_(0),
      nonzeros_compressed_(0),
      have_structure_(false),
      initialized_(false),

      solver_interface_(solver_interface),
      scaling_method_(scaling_method),
      scaling_factors_(NULL),
      airn_(NULL),
      ajcn_(NULL)
  {
    DBG_START_METH("TSymLinearSolver::TSymLinearSolver()",dbg_verbosity);
    DBG_ASSERT(IsValid(solver_interface));
  }

  TSymLinearSolver::~TSymLinearSolver()
  {
    DBG_START_METH("TSymLinearSolver::~TSymLinearSolver()",
                   dbg_verbosity);
    delete [] airn_;
    delete [] ajcn_;
    delete [] scaling_factors_;
  }

  void TSymLinearSolver::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddStringOption2(
      "linear_scaling_on_demand",
      "Flag indicating that linear scaling is only done if it seems required.",
      "yes",
      "no", "Always scale the linear system.",
      "yes", "Start using linear system scaling if solutions seem not good.",
      "This option is only important if a linear scaling method (e.g., mc19) "
      "is used.  If you choose \"no\", then the scaling factors are computed "
      "for every linear system from the start.  This can be quite expensive. "
      "Choosing \"yes\" means that the algorithm will start the scaling "
      "method only when the solutions to the linear system seem not good, and "
      "then use it until the end.");
  }

  bool TSymLinearSolver::InitializeImpl(const OptionsList& options,
                                        const std::string& prefix)
  {
    if (IsValid(scaling_method_)) {
      options.GetBoolValue("linear_scaling_on_demand",
                           linear_scaling_on_demand_, prefix);
    }
    else {
      linear_scaling_on_demand_ = false;
    }
    // This option is registered by OrigIpoptNLP
    options.GetBoolValue("warm_start_same_structure",
                         warm_start_same_structure_, prefix);

    if (!solver_interface_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                       options, prefix)) {
      return false;
    }

    if (!warm_start_same_structure_) {
      // Reset all private data
      atag_=0;
      dim_=0;
      nonzeros_triplet_=0;
      nonzeros_compressed_=0;
      have_structure_=false;

      matrix_format_ = solver_interface_->MatrixFormat();
      switch (matrix_format_) {
        case SparseSymLinearSolverInterface::Dense_Format:
        triplet_to_dense_converter_ = new TripletToDenseConverter(0);
        break;
        case SparseSymLinearSolverInterface::CSR_Format_0_Offset:
        triplet_to_csr_converter_ = new TripletToCSRConverter(0);
        break;
        case SparseSymLinearSolverInterface::CSR_Format_1_Offset:
        triplet_to_csr_converter_ = new TripletToCSRConverter(1);
        break;
        case SparseSymLinearSolverInterface::Triplet_Format:
        triplet_to_csr_converter_ = NULL;
        break;
        default:
        DBG_ASSERT(false && "Invalid MatrixFormat returned from solver interface.");
        return false;
      }
    }
    else {
      ASSERT_EXCEPTION(have_structure_, INVALID_WARMSTART,
                       "TSymLinearSolver called with warm_start_same_structure, but the internal structures are not initialized.");
    }

    // reset the initialize flag to make sure that InitializeStructure
    // is called for the linear solver
    initialized_=false;

    if (IsValid(scaling_method_) && !linear_scaling_on_demand_) {
      use_scaling_ = true;
    }
    else {
      use_scaling_ = false;
    }
    just_switched_on_scaling_ = false;

    bool retval = true;
    if (IsValid(scaling_method_)) {
      IpData().TimingStats().LinearSystemScaling().Start();
      retval = scaling_method_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                           options, prefix);
      IpData().TimingStats().LinearSystemScaling().End();
    }
    return retval;
  }

  ESymSolverStatus
  TSymLinearSolver::MultiSolve(const SymMatrix& sym_A,
                               std::vector<SmartPtr<const Vector> >& rhsV,
                               std::vector<SmartPtr<Vector> >& solV,
                               bool check_NegEVals,
                               Index numberOfNegEVals)
  {
    DBG_START_METH("TSymLinearSolver::MultiSolve",dbg_verbosity);
    DBG_ASSERT(!check_NegEVals || ProvidesInertia());

    // Check if this object has ever seen a matrix If not,
    // allocate memory of the matrix structure and copy the nonzeros
    // structure (it is assumed that this will never change).
    if (!initialized_) {
      ESymSolverStatus retval = InitializeStructure(sym_A);
      if (retval != SYMSOLVER_SUCCESS) {
        return retval;
      }
    }

    DBG_ASSERT(nonzeros_triplet_== TripletHelper::GetNumberEntries(sym_A));

    // Check if the matrix has been changed
    DBG_PRINT((1,"atag_=%d sym_A->GetTag()=%d\n",atag_,sym_A.GetTag()));
    bool new_matrix = sym_A.HasChanged(atag_);
    atag_ = sym_A.GetTag();

    // If a new matrix is encountered, get the array for storing the
    // entries from the linear solver interface, fill in the new
    // values, compute the new scaling factors (if required), and
    // scale the matrix
    if (new_matrix || just_switched_on_scaling_) {
      GiveMatrixToSolver(true, sym_A);
      new_matrix = true;
    }

    // Retrieve the right hand sides and scale if required
    Index nrhs = (Index)rhsV.size();
    Number* rhs_vals = new Number[dim_*nrhs];
    for (Index irhs=0; irhs<nrhs; irhs++) {
      TripletHelper::FillValuesFromVector(dim_, *rhsV[irhs],
                                          &rhs_vals[irhs*(dim_)]);
      if (use_scaling_) {
        IpData().TimingStats().LinearSystemScaling().Start();
        for (Index i=0; i<dim_; i++) {
          rhs_vals[irhs*(dim_)+i] *= scaling_factors_[i];
        }
        IpData().TimingStats().LinearSystemScaling().End();
      }
    }

    bool done = false;
    // Call the linear solver through the interface to solve the
    // system.  This is repeated, if the return values is S_CALL_AGAIN
    // after the values have been restored (this might be necessary
    // for MA27 if the size of the work space arrays was not large
    // enough).
    ESymSolverStatus retval;
    while (!done) {
      const Index* ia;
      const Index* ja;
      if (matrix_format_==SparseSymLinearSolverInterface::Triplet_Format) {
        ia = airn_;
        ja = ajcn_;
      } else if (matrix_format_==SparseSymLinearSolverInterface::Dense_Format) {
        ia = triplet_to_dense_converter_->IA();
        ja = triplet_to_dense_converter_->JA();
      }else {
        IpData().TimingStats().LinearSystemStructureConverter().Start();
        ia = triplet_to_csr_converter_->IA();
        ja = triplet_to_csr_converter_->JA();
        IpData().TimingStats().LinearSystemStructureConverter().End();
      }

      retval = solver_interface_->MultiSolve(new_matrix, ia, ja,
                                             nrhs, rhs_vals, check_NegEVals,
                                             numberOfNegEVals);
      if (retval==SYMSOLVER_CALL_AGAIN) {
        DBG_PRINT((1, "Solver interface asks to be called again.\n"));
        GiveMatrixToSolver(false, sym_A);
      }
      else {
        done = true;
      }
    }

    // If the solve was successful, unscale the solution (if required)
    // and transfer the result into the Vectors
    if (retval==SYMSOLVER_SUCCESS) {
      for (Index irhs=0; irhs<nrhs; irhs++) {
        if (use_scaling_) {
          IpData().TimingStats().LinearSystemScaling().Start();
          for (Index i=0; i<dim_; i++) {
            rhs_vals[irhs*(dim_)+i] *= scaling_factors_[i];
          }
          IpData().TimingStats().LinearSystemScaling().End();
        }
        TripletHelper::PutValuesInVector(dim_, &rhs_vals[irhs*(dim_)],
                                         *solV[irhs]);
      }
    }

    delete[] rhs_vals;

    return retval;
  }

  // Initialize the local copy of the positions of the nonzero
  // elements
  ESymSolverStatus
  TSymLinearSolver::InitializeStructure(const SymMatrix& sym_A)
  {
    DBG_START_METH("TSymLinearSolver::InitializeStructure",
                   dbg_verbosity);
    DBG_ASSERT(!initialized_);

    ESymSolverStatus retval;

    // have_structure_ is already true if this is a warm start for a
    // problem with identical structure
    if (!have_structure_) {

      dim_ = sym_A.Dim();
      nonzeros_triplet_ = TripletHelper::GetNumberEntries(sym_A);

      delete [] airn_;
      delete [] ajcn_;
      airn_ = new Index[nonzeros_triplet_];
      ajcn_ = new Index[nonzeros_triplet_];

      TripletHelper::FillRowCol(nonzeros_triplet_, sym_A, airn_, ajcn_);

      // If the solver wants the compressed format, the converter has to
      // be initialized
      const Index *ia;
      const Index *ja;
      Index nonzeros;
      if (matrix_format_ == SparseSymLinearSolverInterface::Triplet_Format) {
        ia = airn_;
        ja = ajcn_;
        nonzeros = nonzeros_triplet_;
      } else if(matrix_format_ == SparseSymLinearSolverInterface::Dense_Format) {
        nonzeros_compressed_ =
          triplet_to_dense_converter_->InitializeConverter(dim_, nonzeros_triplet_,
              airn_, ajcn_);
        ia = triplet_to_dense_converter_->IA();
        ja = triplet_to_dense_converter_->JA();
        nonzeros = nonzeros_compressed_;
      } else {
        IpData().TimingStats().LinearSystemStructureConverter().Start();
        IpData().TimingStats().LinearSystemStructureConverterInit().Start();
        nonzeros_compressed_ =
          triplet_to_csr_converter_->InitializeConverter(dim_, nonzeros_triplet_,
              airn_, ajcn_);
        IpData().TimingStats().LinearSystemStructureConverterInit().End();
        ia = triplet_to_csr_converter_->IA();
        ja = triplet_to_csr_converter_->JA();
        IpData().TimingStats().LinearSystemStructureConverter().End();
        nonzeros = nonzeros_compressed_;
      }

      retval = solver_interface_->InitializeStructure(dim_, nonzeros, ia, ja);
      if (retval != SYMSOLVER_SUCCESS) {
        return retval;
      }

      // Get space for the scaling factors
      delete [] scaling_factors_;
      if (IsValid(scaling_method_)) {
        IpData().TimingStats().LinearSystemScaling().Start();
        scaling_factors_ = new Number[dim_];
        IpData().TimingStats().LinearSystemScaling().End();
      }

      have_structure_ = true;
    }
    else {
      ASSERT_EXCEPTION(dim_==sym_A.Dim(), INVALID_WARMSTART,
                       "TSymLinearSolver called with warm_start_same_structure, but the problem is solved for the first time.");
      // This is a warm start for identical structure, so we don't need to
      // recompute the nonzeros location arrays
      const Index *ia;
      const Index *ja;
      Index nonzeros;
      if (matrix_format_ == SparseSymLinearSolverInterface::Triplet_Format) {
        ia = airn_;
        ja = ajcn_;
        nonzeros = nonzeros_triplet_;
      } else if (matrix_format_ == SparseSymLinearSolverInterface::Dense_Format) {
        ia = triplet_to_dense_converter_->IA();
        ja = triplet_to_dense_converter_->JA();
        nonzeros = nonzeros_compressed_;
      } else {
        IpData().TimingStats().LinearSystemStructureConverter().Start();
        ia = triplet_to_csr_converter_->IA();
        ja = triplet_to_csr_converter_->JA();
        IpData().TimingStats().LinearSystemStructureConverter().End();
        nonzeros = nonzeros_compressed_;
      }
      retval = solver_interface_->InitializeStructure(dim_, nonzeros, ia, ja);
    }
    initialized_=true;
    return retval;
  }

  Index TSymLinearSolver::NumberOfNegEVals() const
  {
    DBG_START_METH("TSymLinearSolver::NumberOfNegEVals",dbg_verbosity);
    return solver_interface_->NumberOfNegEVals();
  }

  bool TSymLinearSolver::IncreaseQuality()
  {
    DBG_START_METH("TSymLinearSolver::IncreaseQuality",dbg_verbosity);

    if (IsValid(scaling_method_) && !use_scaling_ &&
        linear_scaling_on_demand_) {
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Switching on scaling of the linear system (on demand).\n");
      IpData().Append_info_string("Mc");
      use_scaling_ = true;
      just_switched_on_scaling_ = true;
      return true;
    }

    return solver_interface_->IncreaseQuality();
  }

  bool TSymLinearSolver::ProvidesInertia() const
  {
    DBG_START_METH("TSymLinearSolver::ProvidesInertia",dbg_verbosity);

    return solver_interface_->ProvidesInertia();
  }

  void TSymLinearSolver::GiveMatrixToSolver(bool new_matrix,
      const SymMatrix& sym_A)
  {
    DBG_START_METH("TSymLinearSolver::GiveMatrixToSolver",dbg_verbosity);
    DBG_PRINT((1,"new_matrix = %d\n",new_matrix));

    Number* pa = solver_interface_->GetValuesArrayPtr();
    Number* atriplet;

    if (matrix_format_!=SparseSymLinearSolverInterface::Triplet_Format) {
      atriplet = new Number[nonzeros_triplet_];
    }
    else {
      atriplet = pa;
    }

    //DBG_PRINT_MATRIX(3, "Aunscaled", sym_A);
    TripletHelper::FillValues(nonzeros_triplet_, sym_A, atriplet);
    if (DBG_VERBOSITY()>=3) {
      for (Index i=0; i<nonzeros_triplet_; i++) {
        DBG_PRINT((3, "KKTunscaled(%6d,%6d) = %24.16e\n", airn_[i], ajcn_[i], atriplet[i]));
      }
    }

    if (use_scaling_) {
      IpData().TimingStats().LinearSystemScaling().Start();
      DBG_ASSERT(scaling_factors_);
      if (new_matrix || just_switched_on_scaling_) {
        // only compute scaling factors if the matrix has not been
        // changed since the last call to this method
        bool retval =
          scaling_method_->ComputeSymTScalingFactors(dim_, nonzeros_triplet_,
              airn_, ajcn_,
              atriplet, scaling_factors_);
        if (!retval) {
          Jnlst().Printf(J_ERROR, J_LINEAR_ALGEBRA,
                         "Error during computation of scaling factors.\n");
          THROW_EXCEPTION(ERROR_IN_LINEAR_SCALING_METHOD, "scaling_method_->ComputeSymTScalingFactors returned false.")
        }
        // complain if not in debug mode
        if (Jnlst().ProduceOutput(J_MOREVECTOR, J_LINEAR_ALGEBRA)) {
          for (Index i=0; i<dim_; i++) {
            Jnlst().Printf(J_MOREVECTOR, J_LINEAR_ALGEBRA,
                           "scaling factor[%6d] = %22.17e\n",
                           i, scaling_factors_[i]);
          }
        }
        just_switched_on_scaling_ = false;
      }
      for (Index i=0; i<nonzeros_triplet_; i++) {
        atriplet[i] *=
          scaling_factors_[airn_[i]-1] * scaling_factors_[ajcn_[i]-1];
      }
      if (DBG_VERBOSITY()>=3) {
        for (Index i=0; i<nonzeros_triplet_; i++) {
          DBG_PRINT((3, "KKTscaled(%6d,%6d) = %24.16e\n", airn_[i], ajcn_[i], atriplet[i]));
        }
      }
      IpData().TimingStats().LinearSystemScaling().End();
    }

    if (matrix_format_!=SparseSymLinearSolverInterface::Triplet_Format) {
      IpData().TimingStats().LinearSystemStructureConverter().Start();
      if( matrix_format_ == SparseSymLinearSolverInterface::Dense_Format ){
         triplet_to_dense_converter_->ConvertValues(nonzeros_triplet_, atriplet,
          dim_, pa);
      } else {
         triplet_to_csr_converter_->ConvertValues(nonzeros_triplet_, atriplet,
          nonzeros_compressed_, pa);
      }
      IpData().TimingStats().LinearSystemStructureConverter().End();
      delete[] atriplet;
    }

  }

} // namespace Ipopt
