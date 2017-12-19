// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpPDFullSpaceSolver.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPPDFULLSPACESOLVER_HPP__
#define __IPPDFULLSPACESOLVER_HPP__

#include "IpPDSystemSolver.hpp"
#include "IpAugSystemSolver.hpp"
#include "IpPDPerturbationHandler.hpp"

namespace SimTKIpopt
{

  /** This is the implemetation of the Primal-Dual System, using the
   *  full space approach with a direct linear solver.
   *
   *  A note on the iterative refinement: We perform at least
   *  min_refinement_steps number of iterative refinement steps.  If after
   *  one iterative refinement the quality of the solution (defined in
   *  ResidualRatio) does not improve or the maximal number of
   *  iterative refinement steps is exceeded before the tolerance
   *  residual_ratio_max_ is satisfied, we first ask the linear solver
   *  to solve the system more accurately (e.g. by increasing the
   *  pivot tolerance).  If that doesn't help or is not possible, we
   *  treat the system, as if it is singular (i.e. increase delta's).
   */
  class PDFullSpaceSolver: public PDSystemSolver
  {
  public:
    /** @name /Destructor */
    //@{
    /** Constructor that takes in the Augmented System solver that
     *  is to be used inside
     */
    PDFullSpaceSolver(AugSystemSolver& augSysSolver,
                      PDPerturbationHandler& perturbHandler);

    /** Default destructor */
    virtual ~PDFullSpaceSolver();
    //@}

    /* overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix) override;

    /** Solve the primal dual system, given one right hand side.
     */
    virtual bool Solve(Number alpha,
                       Number beta,
                       const IteratesVector& rhs,
                       IteratesVector& res,
                       bool allow_inexact=false,
                       bool improve_solution=false) override;

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

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
    PDFullSpaceSolver();
    /** Overloaded Equals Operator */
    PDFullSpaceSolver& operator=(const PDFullSpaceSolver&);
    //@}

    /** @name Strategy objects to hold on to. */
    //@{
    /** Pointer to the Solver for the augmented system */
    SmartPtr<AugSystemSolver> augSysSolver_;
    /** Pointer to the Perturbation Handler. */
    SmartPtr<PDPerturbationHandler> perturbHandler_;
    //@}

    /**@name Data about the correction made to the system */
    //@{
    /** A dummy cache to figure out if the deltas are still up to date*/
    CachedResults<void*> dummy_cache_;
    /** Flag indicating if for the current matrix the solution quality
     *  of the augmented system solver has already been increased. */
    bool augsys_improved_;
    //@}

    /** @name Parameters */
    //@{
    /** Minimal number of iterative refinement performed per backsolve */
    Index min_refinement_steps_;
    /** Maximal number of iterative refinement performed per backsolve */
    Index max_refinement_steps_;
    /** Maximal allowed ratio of the norm of the residual over the
     *  norm of the right hand side and solution. */
    Number residual_ratio_max_;
    /** If the residual_ratio is larger than this value after trying
     *  to improve the solution, the linear system is assumed to be
     *  singular and modified. */
    Number residual_ratio_singular_;
    /** Factor defining require improvement to consider iterative
     *  refinement successful. */
    Number residual_improvement_factor_;
    //@}

    /** Internal function for a single backsolve (which will be used
     *  for iterative refinement on the outside).  This method returns
     *  false, if for some reason the linear system could not be
     *  solved (e.g. when the regularization parameter becomes too
     *  large.)
     */
    bool SolveOnce(bool resolve_unmodified,
                   bool pretend_singular,
                   const SymMatrix& W,
                   const Matrix& J_c,
                   const Matrix& J_d,
                   const Matrix& Px_L,
                   const Matrix& Px_U,
                   const Matrix& Pd_L,
                   const Matrix& Pd_U,
                   const Vector& z_L,
                   const Vector& z_U,
                   const Vector& v_L,
                   const Vector& v_U,
                   const Vector& slack_x_L,
                   const Vector& slack_x_U,
                   const Vector& slack_s_L,
                   const Vector& slack_s_U,
                   const Vector& sigma_x,
                   const Vector& sigma_s,
                   Number alpha,
                   Number beta,
                   const IteratesVector& rhs,
                   IteratesVector& res);

    /** Internal function for computing the residual (resid) given the
     * right hand side (rhs) and the solution of the system (res).
     */
    void ComputeResiduals(const SymMatrix& W,
                          const Matrix& J_c,
                          const Matrix& J_d,
                          const Matrix& Px_L,
                          const Matrix& Px_U,
                          const Matrix& Pd_L,
                          const Matrix& Pd_U,
                          const Vector& z_L,
                          const Vector& z_U,
                          const Vector& v_L,
                          const Vector& v_U,
                          const Vector& slack_x_L,
                          const Vector& slack_x_U,
                          const Vector& slack_s_L,
                          const Vector& slack_s_U,
                          const Vector& sigma_x,
                          const Vector& sigma_s,
                          Number alpha,
                          Number beta,
                          const IteratesVector& rhs,
                          const IteratesVector& res,
                          IteratesVector& resid);

    /** Internal function for computing the ratio of the residual
     *  compared to the right hand side and solution.  The smaller
     *  this value, the better the solution. */
    Number ComputeResidualRatio(const IteratesVector& rhs,
                                const IteratesVector& res,
                                const IteratesVector& resid);

    /** @name Auxilliary functions */
    //@{
    /** Compute \f$ x = S^{-1}(r + \alpha Z P^T d)\f$ */
    void SinvBlrmZPTdBr(Number alpha, const Vector& S,
                        const Vector& R, const Vector& Z,
                        const Matrix& P, const Vector&g, Vector& X);
    //@}
  };

} // namespace Ipopt

#endif
