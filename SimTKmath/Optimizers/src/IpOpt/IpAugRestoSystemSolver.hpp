// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpAugRestoSystemSolver.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IP_AUGRESTOSYSTEMSOLVER_HPP__
#define __IP_AUGRESTOSYSTEMSOLVER_HPP__

#include "IpAugSystemSolver.hpp"

namespace SimTKIpopt
{

  /** Class that converts the an augmented system with compound restoration
   *  pieces into a smaller "pivoted" system to be solved with an
   *  existing AugSystemSolver. This is really a decorator that changes
   *  the behavior of the AugSystemSolver to account for the known structure
   *  of the restoration phase.
   */
  class AugRestoSystemSolver: public AugSystemSolver
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor. Here, orig_aug_solver is the object for solving
     *  the original augmented system.  The flag
     *  skip_orig_aug_solver_init indicates, if the initialization
     *  call (to Initialize) should be skipped; this flag will usually
     *  be true, so that the symbolic factorization of the main
     *  algorithm will be used. */
    AugRestoSystemSolver(AugSystemSolver& orig_aug_solver,
                         bool skip_orig_aug_solver_init=true);

    /** Default destructor */
    virtual ~AugRestoSystemSolver();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix) override;

    /** Translate the augmented system (in the full space of the
     *  restoration variables) into the smaller space of the original
     *  variables
     */
    virtual ESymSolverStatus Solve(
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
      const Vector& rhs_x,
      const Vector& rhs_s,
      const Vector& rhs_c,
      const Vector& rhs_d,
      Vector& sol_x,
      Vector& sol_s,
      Vector& sol_c,
      Vector& sol_d,
      bool check_NegEVals,
      Index numberOfNegEVals) override;

    /** Returns the number of negative eigenvalues from the original
     *  augmented system call
     */
    virtual Index NumberOfNegEVals() const override
    {
      return orig_aug_solver_->NumberOfNegEVals();
    }

    /** Query whether inertia is computed by linear solver.
     * Returns true, if linear solver provides inertia.
     */
    virtual bool ProvidesInertia() const override
    {
      return orig_aug_solver_->ProvidesInertia();
    }

    /** Request to increase quality of solution for next solve.  Ask
     *  underlying linear solver to increase quality of solution for
     *  the next solve (e.g. increase pivot tolerance).  Returns
     *  false, if this is not possible (e.g. maximal pivot tolerance
     *  already used.)
     */
    virtual bool IncreaseQuality() override
    {
      return orig_aug_solver_->IncreaseQuality();
    }

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    AugRestoSystemSolver();

    /** Copy Constructor */
    AugRestoSystemSolver(const AugRestoSystemSolver&);

    /** Overloaded Equals Operator */
    void operator=(const AugRestoSystemSolver&);
    //@}

    /**@name Caches for some of the necessary calculated quantities */
    //@{
    CachedResults< SmartPtr<Vector> > neg_omega_c_plus_D_c_cache_;
    CachedResults< SmartPtr<Vector> > neg_omega_d_plus_D_d_cache_;
    CachedResults< SmartPtr<Vector> > sigma_tilde_n_c_inv_cache_;
    CachedResults< SmartPtr<Vector> > sigma_tilde_p_c_inv_cache_;
    CachedResults< SmartPtr<Vector> > sigma_tilde_n_d_inv_cache_;
    CachedResults< SmartPtr<Vector> > sigma_tilde_p_d_inv_cache_;
    CachedResults< SmartPtr<Vector> > d_x_plus_wr_d_cache_;
    CachedResults< SmartPtr<Vector> > rhs_cR_cache_;
    CachedResults< SmartPtr<Vector> > rhs_dR_cache_;
    //@}

    /**@name Methods to calculate the cached quantities */
    //@{
    SmartPtr<const Vector> Neg_Omega_c_plus_D_c(
      const SmartPtr<const Vector>& sigma_tilde_n_c_inv,
      const SmartPtr<const Vector>& sigma_tilde_p_c_inv,
      const Vector* D_c,
      const Vector& any_vec_in_c);

    SmartPtr<const Vector> Neg_Omega_d_plus_D_d(
      const Matrix& Pd_L,
      const SmartPtr<const Vector>& sigma_tilde_n_d_inv,
      const Matrix& neg_Pd_U,
      const SmartPtr<const Vector>& sigma_tilde_p_d_inv,
      const Vector* D_d,
      const Vector& any_vec_in_d);

    /** Sigma tilde is the sum of Sigma and delta_x times the identity */
    SmartPtr<const Vector> Sigma_tilde_n_c_inv(
      const SmartPtr<const Vector>& sigma_tilde_n_c,
      Number delta_x,
      const Vector& any_vec_in_n_c);

    SmartPtr<const Vector> Sigma_tilde_p_c_inv(
      const SmartPtr<const Vector>& sigma_tilde_p_c,
      Number delta_x,
      const Vector& any_vec_in_p_c);

    SmartPtr<const Vector> Sigma_tilde_n_d_inv(
      const SmartPtr<const Vector>& sigma_tilde_n_d,
      Number delta_x,
      const Vector& any_vec_in_n_d);

    SmartPtr<const Vector> Sigma_tilde_p_d_inv(
      const SmartPtr<const Vector>& sigma_tilde_p_d,
      Number delta_x,
      const Vector& any_vec_in_p_d);

    SmartPtr<const Vector> D_x_plus_wr_d(
      const SmartPtr<const Vector>& CD_x0,
      Number factor,
      const Vector& wr_d);

    SmartPtr<const Vector> Rhs_cR(
      const Vector& rhs_c,
      const SmartPtr<const Vector>& sigma_tilde_n_c_inv,
      const Vector& rhs_n_c,
      const SmartPtr<const Vector>& sigma_tilde_p_c_inv,
      const Vector& rhs_p_c);

    SmartPtr<const Vector> Rhs_dR(
      const Vector& rhs_d,
      const SmartPtr<const Vector>& sigma_tilde_n_d_inv,
      const Vector& rhs_n_d,
      const Matrix& pd_L,
      const SmartPtr<const Vector>& sigma_tilde_p_d_inv,
      const Vector& rhs_p_d,
      const Matrix& pd_U);
    //@}

    SmartPtr<AugSystemSolver> orig_aug_solver_;
    bool skip_orig_aug_solver_init_;
  };

} // namespace Ipopt

#endif
