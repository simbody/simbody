// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpStdAugSystemSolver.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IP_STDAUGSYSTEMSOLVER_HPP__
#define __IP_STDAUGSYSTEMSOLVER_HPP__

#include "IpAugSystemSolver.hpp"
#include "IpCompoundMatrix.hpp"
#include "IpCompoundSymMatrix.hpp"
#include "IpCompoundVector.hpp"
#include "IpSumSymMatrix.hpp"
#include "IpDiagMatrix.hpp"
#include "IpIdentityMatrix.hpp"

namespace SimTKIpopt
{

  DECLARE_STD_EXCEPTION(FATAL_ERROR_IN_LINEAR_SOLVER);

  /** Solver for the augmented system for triple type matrices.
   *
   *  The current implemetation assumes that all matrices are of the
   *  type SymTMatrix, and all vectors are of the type DenseVector.
   */
  class StdAugSystemSolver : public AugSystemSolver
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor using only a linear solver object */
    StdAugSystemSolver(SymLinearSolver& LinSolver);

    /** Default destructor */
    virtual ~StdAugSystemSolver();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix) override;

    /** Set up the augmented system and solve it for a set of given
     *  right hand side - implementation for GenTMatrices and
     *  SymTMatrices.
     */
    virtual ESymSolverStatus MultiSolve(
      const SymMatrix* W,
      Number W_factor,
      const Vector* D_x,
      Number delta_x,
      const Vector* D_s,
      Number delta_s,
      const Matrix* J_c,
      const Vector* D_c,
      Number delta_c,
      const Matrix* J_d,
      const Vector* D_d,
      Number delta_d,
      std::vector<SmartPtr<const Vector> >& rhs_xV,
      std::vector<SmartPtr<const Vector> >& rhs_sV,
      std::vector<SmartPtr<const Vector> >& rhs_cV,
      std::vector<SmartPtr<const Vector> >& rhs_dV,
      std::vector<SmartPtr<Vector> >& sol_xV,
      std::vector<SmartPtr<Vector> >& sol_sV,
      std::vector<SmartPtr<Vector> >& sol_cV,
      std::vector<SmartPtr<Vector> >& sol_dV,
      bool check_NegEVals,
      Index numberOfNegEVals) override;

    /** Number of negative eigenvalues detected during last
     * solve.  Returns the number of negative eigenvalues of
     * the most recent factorized matrix.  This must not be called if
     * the linear solver does not compute this quantities (see
     * ProvidesInertia).
     */
    virtual Index NumberOfNegEVals() const override;

    /** Query whether inertia is computed by linear solver.
     * Returns true, if linear solver provides inertia.
     */
    virtual bool ProvidesInertia() const override;

    /** Request to increase quality of solution for next solve.  Ask
     *  underlying linear solver to increase quality of solution for
     *  the next solve (e.g. increase pivot tolerance).  Returns
     *  false, if this is not possible (e.g. maximal pivot tolerance
     *  already used.)
     */
    virtual bool IncreaseQuality() override;

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default constructor. */
    StdAugSystemSolver();
    /** Copy Constructor */
    StdAugSystemSolver(const StdAugSystemSolver&);

    /** Overloaded Equals Operator */
    void operator=(const StdAugSystemSolver&);
    //@}


    /** Create the matrix space for the Compound Sym Matrix that
     *  represents the augmented system. This signifies the "first"
     *  time through and requires all structural knowledge */
    void CreateAugmentedSpace(const SymMatrix& W,
                              const Matrix& J_c,
                              const Matrix& J_d,
                              const Vector& proto_x,
                              const Vector& proto_s,
                              const Vector& proto_c,
                              const Vector& proto_d);


    /** Create the new compound sym matrix that represents the
     *  augmented system. This is done EVERY time Solve is called
     *  with ANY different information */
    void CreateAugmentedSystem(const SymMatrix* W,
                               Number W_factor,
                               const Vector* D_x,
                               Number delta_x,
                               const Vector* D_s,
                               Number delta_s,
                               const Matrix& J_c,
                               const Vector* D_c,
                               Number delta_c,
                               const Matrix& J_d,
                               const Vector* D_d,
                               Number delta_d,
                               const Vector& proto_x,
                               const Vector& proto_s,
                               const Vector& proto_c,
                               const Vector& proto_d);

    /** Check the internal tags and decide if the passed variables are
     *  different from what is in the augmented_system_ */
    bool AugmentedSystemRequiresChange(const SymMatrix* W,
                                       Number W_factor,
                                       const Vector* D_x,
                                       Number delta_x,
                                       const Vector* D_s,
                                       Number delta_s,
                                       const Matrix& J_c,
                                       const Vector* D_c,
                                       Number delta_c,
                                       const Matrix& J_d,
                                       const Vector* D_d,
                                       Number delta_d);

    /** The linear solver object that is to be used to solve the
     *  linear systems.
     */
    SmartPtr<SymLinearSolver> linsolver_;

    /** Spaces for piecing together the augmented system */
    SmartPtr<CompoundSymMatrixSpace> augmented_system_space_;
    SmartPtr<SumSymMatrixSpace> sumsym_space_x_;
    SmartPtr<DiagMatrixSpace> diag_space_x_;
    SmartPtr<DiagMatrixSpace> diag_space_s_;
    SmartPtr<DiagMatrixSpace> diag_space_c_;
    SmartPtr<IdentityMatrixSpace> ident_space_ds_;
    SmartPtr<DiagMatrixSpace> diag_space_d_;

    SmartPtr<CompoundVectorSpace> augmented_vector_space_;

    /**@name Tags and values to track in order to decide whether the
       matrix has to be updated compared to the most recent call of
       the Set method.
     */
    //@{
    /** Tag for W matrix.  If W has been given to Set as NULL, then
     *  this tag is set to 0
     */
    TaggedObject::Tag w_tag_;
    /** Most recent value of W_factor */
    double w_factor_;
    /** Tag for D_x vector, representing the diagonal matrix D_x.  If
     *  D_x has been given to Set as NULL, then this tag is set to 0
     */
    TaggedObject::Tag d_x_tag_;
    /** Most recent value of delta_x from Set method */
    double delta_x_;
    /** Tag for D_s vector, representing the diagonal matrix D_s.  If
     *  D_s has been given to Set as NULL, then this tag is set to 0
     */
    TaggedObject::Tag d_s_tag_;
    /** Most recent value of delta_s from Set method */
    double delta_s_;
    /** Tag for J_c matrix.  If J_c has been given to Set as NULL, then
     *  this tag is set to 0
     */
    TaggedObject::Tag j_c_tag_;
    /** Tag for D_c vector, representing the diagonal matrix D_c.  If
     *  D_c has been given to Set as NULL, then this tag is set to 0
     */
    TaggedObject::Tag d_c_tag_;
    /** Most recent value of delta_c from Set method */
    double delta_c_;
    /** Tag for J_d matrix.  If J_d has been given to Set as NULL, then
     *  this tag is set to 0
     */
    TaggedObject::Tag j_d_tag_;
    /** Tag for D_d vector, representing the diagonal matrix D_d.  If
     *  D_d has been given to Set as NULL, then this tag is set to 0
     */
    TaggedObject::Tag d_d_tag_;
    /** Most recent value of delta_d from Set method */
    double delta_d_;

    /** This is the tag of the matrix storing the augmented system.  Since
     *  this object owns this matrix, no changes should happen outside.
     *  However, since it is given away as a smart pointer, someone outside
     *  might change it.  For debugging purposes, we now track its tag as
     *  well.
     */
    TaggedObject::Tag augsys_tag_;
    //@}

    /** The resulting augmented matrix.
     *  This matrix is stored as follows:  First we have the diagonal elements
     *  for the upper left block (for D_W and delta_W), then the elements for
     *  the Hessian W, then the Jacobian A, and finally the diagonal elements
     *  for the lower right block (for D_C and delta_C).
     */
    SmartPtr<CompoundSymMatrix> augmented_system_;

    /** A copy of a previous W used in the augmented_system_. Since Solve can
     *  be called with a NULL W, we keep a copy of the last W passed to keep
     *  the nonzero structure of the augmented_system_ consistent */
    SmartPtr<const SymMatrix> old_w_;

    /** @name Algorithmic parameters */
    //@{
    /** Flag indicating whether the TNLP with identical structure has
     *  already been solved before. */
    bool warm_start_same_structure_;
    //@}
  };

} // namespace Ipopt

#endif
