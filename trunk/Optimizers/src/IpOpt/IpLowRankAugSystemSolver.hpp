// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpLowRankAugSystemSolver.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter                IBM    2005-12-27

#ifndef __IP_LOWRANKAUGSYSTEMSOLVER_HPP__
#define __IP_LOWRANKAUGSYSTEMSOLVER_HPP__

#include "IpAugSystemSolver.hpp"
#include "IpDenseGenMatrix.hpp"
#include "IpMultiVectorMatrix.hpp"
#include "IpDiagMatrix.hpp"

namespace Ipopt
{

  /** Solver for the augmented system with LowRankUpdateSymMatrix
   *  Hessian matrices.
   */
  class LowRankAugSystemSolver : public AugSystemSolver
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor using only a linear solver object */
    LowRankAugSystemSolver(AugSystemSolver& aug_system_solver);

    /** Default destructor */
    virtual ~LowRankAugSystemSolver();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix);

    /** Set up the augmented system and solve it for a given right hand
     *  side.
     */
    virtual ESymSolverStatus Solve(
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
      Index numberOfNegEVals);

    /** Number of negative eigenvalues detected during last
     * solve.  Returns the number of negative eigenvalues of
     * the most recent factorized matrix.  This must not be called if
     * the linear solver does not compute this quantities (see
     * ProvidesInertia).
     */
    virtual Index NumberOfNegEVals() const;

    /** Query whether inertia is computed by linear solver.
     * Returns true, if linear solver provides inertia.
     */
    virtual bool ProvidesInertia() const;

    /** Request to increase quality of solution for next solve.  Ask
     *  underlying linear solver to increase quality of solution for
     *  the next solve (e.g. increase pivot tolerance).  Returns
     *  false, if this is not possible (e.g. maximal pivot tolerance
     *  already used.)
     */
    virtual bool IncreaseQuality();

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
    LowRankAugSystemSolver();
    /** Copy Constructor */
    LowRankAugSystemSolver(const LowRankAugSystemSolver&);

    /** Overloaded Equals Operator */
    void operator=(const LowRankAugSystemSolver&);
    //@}

    /** The augmented system solver object that should be used for the
     *  factorization of the augmented system without the low-rank
     *  update.
     */
    SmartPtr<AugSystemSolver> aug_system_solver_;

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
    //@}

    /** @name Information to be stored in order to resolve for the
     *  same matrix with a different right hand side. */
    //@{
    bool first_call_;
    SmartPtr<DenseGenMatrix> J1_;
    SmartPtr<DenseGenMatrix> J2_;
    SmartPtr<MultiVectorMatrix> Vtilde1_;
    SmartPtr<MultiVectorMatrix> Utilde2_;
    /** Hessian Matrix passed to the augmented system solver solving
     *  the matrix without the low-rank update. */
    SmartPtr<DiagMatrix> Wdiag_;
    /** Vector space for Compound vectors that capture the entire
     *  right hand side and solution vectors .*/
    SmartPtr<const CompoundVectorSpace> compound_sol_vecspace_;
    //@}

    /** Stores the number of negative eigenvalues detected during most
     *  recent factorization.  This is what is returned by
     *  NumberOfNegEVals() of this class.  It usually is the number of
     *  negative eigenvalues retured from the aug_system_solver solve,
     *  but if a Cholesky factorization could not be performed, the
     *  returned value is one more than this what the
     *  aug_system_solver returned. */
    Index num_neg_evals_;

    /** @name Internal functions */
    //@{
    /** Method for updating the factorization, including J1_, J2_,
     *  Vtilde1_, Utilde2, Wdiag_, compound_sol_vecspace_ */
    ESymSolverStatus UpdateFactorization(
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
      Index numberOfNegEVals);

    /** Method for solving the augmented system without low-rank
     *  update for multiple right hand sides that are provided as
     *  MultiVectorMatrix.  The result is returned as a
     *  MultiVectorMatrix in Vtilde1.  V_x and Vtilde1_x are V and
     *  Vtilde1 in the x-space. */
    ESymSolverStatus SolveMultiVector(
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
      SmartPtr<MultiVectorMatrix>& Vtilde1,
      SmartPtr<MultiVectorMatrix>& Vtilde1_x,
      bool check_NegEVals,
      Index numberOfNegEVals);

    /** Method that compares the tags of the data for the matrix with
     *  those from the previous call.  Returns true, if there was a
     *  change and the factorization has to be updated. */
    bool AugmentedSystemRequiresChange(
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
      double delta_d);
    //@}

  };

} // namespace Ipopt

#endif
