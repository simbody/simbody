// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpSparseSymLinearSolverInterface.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-03-17

#ifndef __IPSPARSESYMLINEARSOLVERINTERFACE_HPP__
#define __IPSPARSESYMLINEARSOLVERINTERFACE_HPP__

#include "IpUtils.hpp"
#include "IpAlgStrategy.hpp"
#include "IpSymLinearSolver.hpp"

namespace SimTKIpopt
{

  /** Base class for interfaces to symmetric indefinite linear solvers
   *  for sparse matrices.
   *
   *  This defines the general interface to linear solvers for sparse
   *  symmetric indefinite matrices.  The matrices can be provided
   *  either in "triplet format" (like for Harwell's MA27 solver), or
   *  in compressed sparse row (CSR) format for the lower triangular
   *  part of the symmetric matrix.
   *
   *  The solver should be able to compute the interia of the matrix,
   *  or more specifically, the number of negative eigenvalues in the
   *  factorized matrix.
   *
   *  This interface is used by the calling objective in the following
   *  way: 
   *
   *  1. The InitializeImpl method is called at the very beginning
   *  (for every optimization run), which allows the linear solver
   *  object to retrieve options given in the OptionsList (such as
   *  pivot tolerances etc).  At this point, some internal data can
   *  also be initialized.
   *
   *  2. The calling class calls MatrixFormat to find out which matrix
   *  representation the linear solver requires.  The possible options
   *  are Triplet_Format, as well as CSR_Format_0_Offset and
   *  CSR_Format_1_Offset.  The difference between the last two is
   *  that for CSR_Format_0_Offset the couning of the element position
   *  in the ia and ja arrays starts are 0 (C-style numbering),
   *  whereas for the other one it starts at 1 (Fortran-style
   *  numbering).
   *
   *  3. After this, the InitializeStructure method is called (once).
   *  Here, the structure of the matrix is provided.  If the linear
   *  solver requires a symbolic preprocessing phase that can be done
   *  without knowledge of the matrix element values, it can be done
   *  here.
   *
   *  4. The calling class will request an array for storing the
   *  actual values for a matrix using the GetValuesArrayPtr method.
   *  This array must be at least as large as the number of nonzeros
   *  in the matrix (as given to this class by the InitializeStructure
   *  method call).  After a call of this method, the calling class
   *  will fill this array with the actual values of the matrix.
   *
   *  5. Every time lateron, when actual solves of a linear system is
   *  requested, the calling class will call the MultiSolve to request
   *  the solve, possibly for mulitple right-hand sides.  The flag
   *  new_matrix then indicates if the values of the matrix have
   *  changed and if a factorization is required, or if an old
   *  factorization can be used to do the solve.
   *
   *  Note that the GetValuesArrayPtr method will be called before
   *  every call of MultiSolve with new_matrix=true, or before a
   *  renewed call of MultiSolve if the most previous return value was
   *  SYMSOLV_CALL_AGAIN.
   *
   *  6. The calling class might request with NumberOfNegEVals the
   *  number of the negative eigenvalues for the original matrix that
   *  were detected during the most recently performed factorization.
   *
   *  7. The calling class might ask the linear solver to increase the
   *  quality of the solution.  For example, if the linear solver uses
   *  a pivot tolerance, a larger value should be used for the next
   *  solve (which might require a refactorization).
   *
   *  8. Finally, when the destructor is called, the internal storage,
   *  also in the linear solver, should be released.
   *
   *  Note, if the matrix is given in triplet format, entries might be
   *  listed multiple times, in which case the corresponsing elements
   *  have to be added.
   *
   *  A note for warm starts: If the option
   *  "warm_start_same_structure" is specified with "yes", the
   *  algorithm assumes that a problem with the same sparsity
   *  structure is solved for a repeated time.  In that case, the
   *  linear solver might reuse information from the previous
   *  optimization.  See Ma27TSolverInterface for an example.
  */
  class SparseSymLinearSolverInterface: public AlgorithmStrategyObject
  {
  public:
    /** Enum to specify sparse matrix format. */
    enum EMatrixFormat {
      /** Triplet (MA27) format. */
      Triplet_Format,
      /** Compressed sparse row format for lower triangular part, with
       *  0 offset. */
      CSR_Format_0_Offset,
      /** Compressed sparse row format for lower triangular part, with
       *  1 offset. */
      CSR_Format_1_Offset,
      /** Dense format 
       *  1 offset. */
      Dense_Format 
    };
    /** @name Constructor/Destructor */
    //@{
    SparseSymLinearSolverInterface()
    {}

    virtual ~SparseSymLinearSolverInterface()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix) = 0;

    /** @name Methods for requesting solution of the linear system. */
    //@{
    /** Method for initializing internal stuctures.  Here, ndim gives
     *  the number of rows and columns of the matrix, nonzeros give
     *  the number of nonzero elements, and ia and ja give the
     *  positions of the nonzero elements, given in the matrix format
     *  determined by MatrixFormat.
     */
    virtual ESymSolverStatus InitializeStructure(Index dim, Index nonzeros,
        const Index* ia,
        const Index* ja) = 0;

    /** Method returning an internal array into which the nonzero
     *  elements (in the same order as ja) will be stored by the
     *  calling routine before a call to MultiSolve with a
     *  new_matrix=true (or after a return of MultiSolve with
     *  SYMSOLV_CALL_AGAIN). The returned array must have space for at
     *  least nonzero elements. */
    virtual Number* GetValuesArrayPtr() = 0;

    /** Solve operation for multiple right hand sides.  Solves the
     *  linear system A * x = b with multiple right hand sides, where
     *  A is the symmtric indefinite matrix.  Here, ia and ja give the
     *  positions of the values (in the required matrix data format).
     *  The actual values of the matrix will have been given to this
     *  object by copying them into the array provided by
     *  GetValuesArrayPtr. ia and ja are identical to the ones given
     *  to InitializeStructure.  The flag new_matrix is set to true,
     *  if the values of the matrix has changed, and a refactorzation
     *  is required.
     *
     *  The return code is SYMSOLV_SUCCESS if the factorization and
     *  solves were successful, SYMSOLV_SINGULAR if the linear system
     *  is singular, and SYMSOLV_WRONG_INERTIA if check_NegEVals is
     *  true and the number of negative eigenvalues in the matrix does
     *  not match numberOfNegEVals.  If SYMSOLV_CALL_AGAIN is
     *  returned, then the calling function will request the pointer
     *  for the array for storing a again (with GetValuesPtr), write
     *  the values of the nonzero elements into it, and call this
     *  MultiSolve method again with the same right-hand sides.  (This
     *  can be done, for example, if the linear solver realized it
     *  does not have sufficient memory and needs to redo the
     *  factorization; e.g., for MA27.)
     *
     *  The number of right-hand sides is given by nrhs, the values of
     *  the right-hand sides are given in rhs_vals (one full right-hand
     *  side stored immediately after the other), and solutions are
     *  to be returned in the same array.
     *
     *  check_NegEVals will not be chosen true, if ProvidesInertia()
     *  returns false.
     */
    virtual ESymSolverStatus MultiSolve(bool new_matrix,
                                        const Index* ia,
                                        const Index* ja,
                                        Index nrhs,
                                        Number* rhs_vals,
                                        bool check_NegEVals,
                                        Index numberOfNegEVals)=0;

    /** Number of negative eigenvalues detected during last
     *  factorization.  Returns the number of negative eigenvalues of
     *  the most recent factorized matrix.  This must not be called if
     *  the linear solver does not compute this quantities (see
     *  ProvidesInertia).
     */
    virtual Index NumberOfNegEVals() const =0;
    //@}

    //* @name Options of Linear solver */
    //@{
    /** Request to increase quality of solution for next solve.  The
     *  calling class asks linear solver to increase quality of
     *  solution for the next solve (e.g. increase pivot tolerance).
     *  Returns false, if this is not possible (e.g. maximal pivot
     *  tolerance already used.)
     */
    virtual bool IncreaseQuality() =0;

    /** Query whether inertia is computed by linear solver.  Returns
     *  true, if linear solver provides inertia.
     */
    virtual bool ProvidesInertia() const =0;

    /** Query of requested matrix type that the linear solver
     *  understands.
     */
    virtual EMatrixFormat MatrixFormat() const =0;
    //@}
  };


} // namespace Ipopt

#endif
