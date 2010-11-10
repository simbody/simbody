// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpPDFullSpaceSolver.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpPDFullSpaceSolver.hpp"
#include "IpDebug.hpp"

namespace Ipopt
{

#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  PDFullSpaceSolver::PDFullSpaceSolver(AugSystemSolver& augSysSolver,
                                       PDPerturbationHandler& perturbHandler)
      :
      PDSystemSolver(),
      augSysSolver_(&augSysSolver),
      perturbHandler_(&perturbHandler),
      dummy_cache_(1)
  {
    DBG_START_METH("PDFullSpaceSolver::PDFullSpaceSolver",dbg_verbosity);
  }

  PDFullSpaceSolver::~PDFullSpaceSolver()
  {
    DBG_START_METH("PDFullSpaceSolver::~PDFullSpaceSolver()",dbg_verbosity);
  }

  void PDFullSpaceSolver::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedIntegerOption(
      "min_refinement_steps",
      "Minimum number of iterative refinement steps per linear system solve.",
      0, 1,
      "Iterative refinement (on the full unsymmetric system) is performed for "
      "each right hand side.  This option determines the minimum number "
      "of iterative refinements (i.e. at least \"min_refinement_steps\" "
      "iterative refinement steps are enforced per right hand side.)");
    roptions->AddLowerBoundedIntegerOption(
      "max_refinement_steps",
      "Maximum number of iterative refinement steps per linear system solve.",
      0, 10,
      "Iterative refinement (on the full unsymmetric system) is performed for "
      "each right hand side.  This option determines the maximum number "
      "of iterative refinement steps.");
    roptions->AddLowerBoundedNumberOption(
      "residual_ratio_max",
      "Iterative refinement tolerance",
      0.0, true, 1e-10,
      "Iterative refinement is performed until the residual test ratio is "
      "less than this tolerance (or until \"max_refinement_steps\" refinement "
      "steps are performed).");
    roptions->AddLowerBoundedNumberOption(
      "residual_ratio_singular",
      "Threshold for declaring linear system singular after failed iterative refinement.",
      0.0, true, 1e-5,
      "If the residual test ratio is larger than this value after failed "
      "iterative refinement, the algorithm pretends that the linear system is "
      "singular.");
    // ToDo Think about following option - are the correct norms used?
    roptions->AddLowerBoundedNumberOption(
      "residual_improvement_factor",
      "Minimal required reduction of residual test ratio in iterative refinement.",
      0.0, true, 0.999999999,
      "If the improvement of the residual test ratio made by one iterative "
      "refinement step is not better than this factor, iterative refinement "
      "is aborted.");
  }


  bool PDFullSpaceSolver::InitializeImpl(const OptionsList& options,
                                         const std::string& prefix)
  {
    // Check for the algorithm options
    options.GetIntegerValue("min_refinement_steps", min_refinement_steps_, prefix);
    options.GetIntegerValue("max_refinement_steps", max_refinement_steps_, prefix);
    ASSERT_EXCEPTION(max_refinement_steps_ >= min_refinement_steps_, OPTION_INVALID,
                     "Option \"max_refinement_steps\": This value must be larger than or equal to min_refinement_steps (default 1)");

    options.GetNumericValue("residual_ratio_max", residual_ratio_max_, prefix);
    options.GetNumericValue("residual_ratio_singular", residual_ratio_singular_, prefix);
    ASSERT_EXCEPTION(residual_ratio_singular_ >= residual_ratio_max_, OPTION_INVALID,
                     "Option \"residual_ratio_singular\": This value must be not smaller than residual_ratio_max.");
    options.GetNumericValue("residual_improvement_factor", residual_improvement_factor_, prefix);

    // Reset internal flags and data
    augsys_improved_ = false;

    if (!augSysSolver_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                   options, prefix)) {
      return false;
    }

    return perturbHandler_->Initialize(Jnlst(), IpNLP(), IpData(), IpCq(),
                                       options, prefix);
  }

  bool PDFullSpaceSolver::Solve(Number alpha,
                                Number beta,
                                const IteratesVector& rhs,
                                IteratesVector& res,
                                bool allow_inexact,
                                bool improve_solution /* = false */)
  {
    DBG_START_METH("PDFullSpaceSolver::Solve",dbg_verbosity);
    DBG_ASSERT(!allow_inexact || !improve_solution);
    DBG_ASSERT(!improve_solution || beta==0.);

    // Timing of PDSystem solver starts here
    IpData().TimingStats().PDSystemSolverTotal().Start();

    DBG_PRINT_VECTOR(2, "rhs_x", *rhs.x());
    DBG_PRINT_VECTOR(2, "rhs_s", *rhs.s());
    DBG_PRINT_VECTOR(2, "rhs_c", *rhs.y_c());
    DBG_PRINT_VECTOR(2, "rhs_d", *rhs.y_d());
    DBG_PRINT_VECTOR(2, "rhs_zL", *rhs.z_L());
    DBG_PRINT_VECTOR(2, "rhs_zU", *rhs.z_U());
    DBG_PRINT_VECTOR(2, "rhs_vL", *rhs.v_L());
    DBG_PRINT_VECTOR(2, "rhs_vU", *rhs.v_U());
    DBG_PRINT_VECTOR(2, "res_x in", *res.x());
    DBG_PRINT_VECTOR(2, "res_s in", *res.s());
    DBG_PRINT_VECTOR(2, "res_c in", *res.y_c());
    DBG_PRINT_VECTOR(2, "res_d in", *res.y_d());
    DBG_PRINT_VECTOR(2, "res_zL in", *res.z_L());
    DBG_PRINT_VECTOR(2, "res_zU in", *res.z_U());
    DBG_PRINT_VECTOR(2, "res_vL in", *res.v_L());
    DBG_PRINT_VECTOR(2, "res_vU in", *res.v_U());

    // if beta is nonzero, keep a copy of the incoming values in res_ */
    SmartPtr<IteratesVector> copy_res;
    if (beta != 0.) {
      copy_res = res.MakeNewIteratesVectorCopy();
    }

    // Receive data about matrix
    SmartPtr<const Vector> x = IpData().curr()->x();
    SmartPtr<const Vector> s = IpData().curr()->s();
    SmartPtr<const SymMatrix> W = IpData().W();
    SmartPtr<const Matrix> J_c = IpCq().curr_jac_c();
    SmartPtr<const Matrix> J_d = IpCq().curr_jac_d();
    SmartPtr<const Matrix> Px_L = IpNLP().Px_L();
    SmartPtr<const Matrix> Px_U = IpNLP().Px_U();
    SmartPtr<const Matrix> Pd_L = IpNLP().Pd_L();
    SmartPtr<const Matrix> Pd_U = IpNLP().Pd_U();
    SmartPtr<const Vector> z_L = IpData().curr()->z_L();
    SmartPtr<const Vector> z_U = IpData().curr()->z_U();
    SmartPtr<const Vector> v_L = IpData().curr()->v_L();
    SmartPtr<const Vector> v_U = IpData().curr()->v_U();
    SmartPtr<const Vector> slack_x_L = IpCq().curr_slack_x_L();
    SmartPtr<const Vector> slack_x_U = IpCq().curr_slack_x_U();
    SmartPtr<const Vector> slack_s_L = IpCq().curr_slack_s_L();
    SmartPtr<const Vector> slack_s_U = IpCq().curr_slack_s_U();
    SmartPtr<const Vector> sigma_x = IpCq().curr_sigma_x();
    SmartPtr<const Vector> sigma_s = IpCq().curr_sigma_s();
    DBG_PRINT_VECTOR(2, "Sigma_x", *sigma_x);
    DBG_PRINT_VECTOR(2, "Sigma_s", *sigma_s);

    bool done = false;
    // The following flag is set to true, if we asked the linear
    // solver to improve the quality of the solution in
    // the next solve
    bool resolve_with_better_quality = false;
    // the following flag is set to true, if iterative refinement
    // failed and we want to try if a modified system is able to
    // remedy that problem by pretending the matrix is singular
    bool pretend_singular = false;
    bool pretend_singular_last_time = false;

    // Beginning of loop for solving the system (including all
    // modifications for the linear system to ensure good solution
    // quality)
    while (!done) {

      // if improve_solution is true, we are given already a solution
      // from the calling function, so we can skip the first solve
      bool solve_retval = true;
      if (!improve_solution) {
        solve_retval =
          SolveOnce(resolve_with_better_quality, pretend_singular,
                    *W, *J_c, *J_d, *Px_L, *Px_U, *Pd_L, *Pd_U, *z_L, *z_U,
                    *v_L, *v_U, *slack_x_L, *slack_x_U, *slack_s_L, *slack_s_U,
                    *sigma_x, *sigma_s, 1., 0., rhs, res);
        resolve_with_better_quality = false;
        pretend_singular = false;
      }
      improve_solution = false;

      if (!solve_retval) {
        // If system seems not to be solvable, we return with false
        // and let the calling routine deal with it.
        IpData().TimingStats().PDSystemSolverTotal().End();
        return false;
      }

      if (allow_inexact) {
        // no safety checks required
        break;
      }

      // Get space for the residual
      SmartPtr<IteratesVector> resid = res.MakeNewIteratesVector(true);

      // ToDo don't to that after max refinement?
      ComputeResiduals(*W, *J_c, *J_d, *Px_L, *Px_U, *Pd_L, *Pd_U,
                       *z_L, *z_U, *v_L, *v_U, *slack_x_L, *slack_x_U,
                       *slack_s_L, *slack_s_U, *sigma_x, *sigma_s,
                       alpha, beta, rhs, res, *resid);

      Number residual_ratio =
        ComputeResidualRatio(rhs, res, *resid);
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "residual_ratio = %e\n", residual_ratio);
      Number residual_ratio_old = residual_ratio;

      // Beginning of loop for iterative refinement
      Index num_iter_ref = 0;
      bool quit_refinement = false;
      while (!allow_inexact && !quit_refinement &&
             (num_iter_ref < min_refinement_steps_ ||
              residual_ratio > residual_ratio_max_) ) {

        // To the next back solve
        solve_retval =
          SolveOnce(resolve_with_better_quality, false,
                    *W, *J_c, *J_d, *Px_L, *Px_U, *Pd_L, *Pd_U, *z_L, *z_U,
                    *v_L, *v_U, *slack_x_L, *slack_x_U, *slack_s_L, *slack_s_U,
                    *sigma_x, *sigma_s, -1., 1., *resid, res);
        ASSERT_EXCEPTION(solve_retval, INTERNAL_ABORT,
                         "SolveOnce returns false during iterative refinement.");

        ComputeResiduals(*W, *J_c, *J_d, *Px_L, *Px_U, *Pd_L, *Pd_U,
                         *z_L, *z_U, *v_L, *v_U, *slack_x_L, *slack_x_U,
                         *slack_s_L, *slack_s_U, *sigma_x, *sigma_s,
                         alpha, beta, rhs, res, *resid);

        residual_ratio =
          ComputeResidualRatio(rhs, res, *resid);
        Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                       "residual_ratio = %e\n", residual_ratio);

        num_iter_ref++;
        // Check if we have to give up on iterative refinement
        if (num_iter_ref>min_refinement_steps_ &&
            (num_iter_ref>max_refinement_steps_ ||
             residual_ratio>residual_improvement_factor_*residual_ratio_old)) {

          Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                         "Iterative refinement failed with residual_ratio = %e\n", residual_ratio);
          quit_refinement = true;

          // Pretend singularity only once - if it didn't help, we
          // have to live with what we got so far
          resolve_with_better_quality = false;
          DBG_PRINT((1, "pretend_singular = %d\n", pretend_singular));
          if (!pretend_singular_last_time) {
            // First try if we can ask the augmented system solver to
            // improve the quality of the solution (only if that hasn't
            // been done before for this linear system)
            if (!augsys_improved_) {
              Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                             "Asking augmented system solver to improve quality of its solutions.\n");
              augsys_improved_ = augSysSolver_->IncreaseQuality();
              if (augsys_improved_) {
                IpData().Append_info_string("q");
                resolve_with_better_quality = true;
              }
              else {
                // solver said it cannot improve quality, so let
                // possibly conclude that the current modification is
                // singular
                pretend_singular = true;
              }
            }
            else {
              // we had already asked the solver before to improve the
              // quality of the solution, so let's now pretend that the
              // modification is possibly singular
              pretend_singular = true;
            }
            pretend_singular_last_time = pretend_singular;
            if (pretend_singular) {
              // let's only conclude that the current linear system
              // including modifications is singular, if the residual is
              // quite bad
              if (residual_ratio < residual_ratio_singular_) {
                pretend_singular = false;
                IpData().Append_info_string("S");
                Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                               "Just accept current solution.\n");
              }
              else {
                IpData().Append_info_string("s");
                Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                               "Pretend that the current system (including modifications) is singular.\n");
              }
            }
          }
          else {
            pretend_singular = false;
            DBG_PRINT((1,"Resetting pretend_singular to false.\n"));
          }
        }

        residual_ratio_old = residual_ratio;
      } // End of loop for iterative refinement

      done = !(resolve_with_better_quality) && !(pretend_singular);

    } // End of loop for solving the linear system (incl. modifications)

    // Finally let's assemble the res result vectors
    if (alpha != 0.) {
      res.Scal(alpha);
    }

    if (beta != 0.) {
      res.Axpy(beta, *copy_res);
    }

    DBG_PRINT_VECTOR(2, "res_x", *res.x());
    DBG_PRINT_VECTOR(2, "res_s", *res.s());
    DBG_PRINT_VECTOR(2, "res_c", *res.y_c());
    DBG_PRINT_VECTOR(2, "res_d", *res.y_d());
    DBG_PRINT_VECTOR(2, "res_zL", *res.z_L());
    DBG_PRINT_VECTOR(2, "res_zU", *res.z_U());
    DBG_PRINT_VECTOR(2, "res_vL", *res.v_L());
    DBG_PRINT_VECTOR(2, "res_vU", *res.v_U());

    IpData().TimingStats().PDSystemSolverTotal().End();

    return true;
  }

  bool PDFullSpaceSolver::SolveOnce(bool resolve_with_better_quality,
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
                                    IteratesVector& res)
  {
    // TO DO LIST:
    //
    // 1. decide for reasonable return codes (e.g. fatal error, too
    //    ill-conditioned...)
    // 2. Make constants parameters that can be set from the outside
    // 3. Get Information out of Ipopt structures
    // 4. add heuristic for structurally singular problems
    // 5. see if it makes sense to distinguish delta_x and delta_s,
    //    or delta_c and delta_d
    // 6. increase pivot tolerance if number of get evals so too small
    DBG_START_METH("PDFullSpaceSolver::SolveOnce",dbg_verbosity);

    IpData().TimingStats().PDSystemSolverSolveOnce().Start();

    // Compute the right hand side for the augmented system formulation
    SmartPtr<Vector> augRhs_x = rhs.x()->MakeNewCopy();
    Px_L.AddMSinvZ(1.0, slack_x_L, *rhs.z_L(), *augRhs_x);
    Px_U.AddMSinvZ(-1.0, slack_x_U, *rhs.z_U(), *augRhs_x);

    SmartPtr<Vector> augRhs_s = rhs.s()->MakeNewCopy();
    Pd_L.AddMSinvZ(1.0, slack_s_L, *rhs.v_L(), *augRhs_s);
    Pd_U.AddMSinvZ(-1.0, slack_s_U, *rhs.v_U(), *augRhs_s);

    // Get space into which we can put the solution of the augmented system
    SmartPtr<IteratesVector> sol = res.MakeNewIteratesVector(true);

    // Now check whether any data has changed
    std::vector<const TaggedObject*> deps(13);
    deps[0] = &W;
    deps[1] = &J_c;
    deps[2] = &J_d;
    deps[3] = &z_L;
    deps[4] = &z_U;
    deps[5] = &v_L;
    deps[6] = &v_U;
    deps[7] = &slack_x_L;
    deps[8] = &slack_x_U;
    deps[9] = &slack_s_L;
    deps[10] = &slack_s_U;
    deps[11] = &sigma_x;
    deps[12] = &sigma_s;
    void* dummy;
    bool uptodate = dummy_cache_.GetCachedResult(dummy, deps);
    if (!uptodate) {
      dummy_cache_.AddCachedResult(dummy, deps);
      augsys_improved_ = false;
    }
    // improve_current_solution can only be true, if that system has
    // been solved before
    DBG_ASSERT((!resolve_with_better_quality && !pretend_singular) || uptodate);

    ESymSolverStatus retval;
    if (uptodate && !pretend_singular) {

      // Get the perturbation values
      Number delta_x;
      Number delta_s;
      Number delta_c;
      Number delta_d;
      perturbHandler_->CurrentPerturbation(delta_x, delta_s, delta_c, delta_d);

      // No need to go throught the pain of finding the appropriate
      // values for the deltas, because the matrix hasn't changed since
      // the last call.  So, just call the Solve Method
      //
      // Note: resolve_with_better_quality is true, then the Solve
      // method has already asked the augSysSolver to increase the
      // quality at the end solve, and we are now getting the solution
      // with that better quality
      retval = augSysSolver_->Solve(&W, 1.0, &sigma_x, delta_x,
                                    &sigma_s, delta_s, &J_c, NULL,
                                    delta_c, &J_d, NULL, delta_d,
                                    *augRhs_x, *augRhs_s, *rhs.y_c(), *rhs.y_d(),
                                    *sol->x_NonConst(), *sol->s_NonConst(),
                                    *sol->y_c_NonConst(), *sol->y_d_NonConst(),
                                    false, 0);
      if (retval!=SYMSOLVER_SUCCESS) {
        IpData().TimingStats().PDSystemSolverSolveOnce().End();
        return false;
      }
    }
    else {
      Index numberOfEVals=rhs.y_c()->Dim()+rhs.y_d()->Dim();
      // counter for the number of trial evaluations
      // (ToDo is not at the correct place)
      Index count = 0;

      // Get the very first perturbation values from the perturbation
      // Handler
      Number delta_x;
      Number delta_s;
      Number delta_c;
      Number delta_d;
      perturbHandler_->ConsiderNewSystem(delta_x, delta_s, delta_c, delta_d);

      retval = SYMSOLVER_SINGULAR;
      bool fail = false;

      while (retval!= SYMSOLVER_SUCCESS && !fail) {

        if (pretend_singular) {
          retval = SYMSOLVER_SINGULAR;
          pretend_singular = false;
        }
        else {
          count++;
          Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                         "Solving system with delta_x=%e delta_s=%e\n                    delta_c=%e delta_d=%e\n",
                         delta_x, delta_s, delta_c, delta_d);
          retval = augSysSolver_->Solve(&W, 1.0, &sigma_x, delta_x,
                                        &sigma_s, delta_s, &J_c, NULL,
                                        delta_c, &J_d, NULL, delta_d,
                                        *augRhs_x, *augRhs_s, *rhs.y_c(), *rhs.y_d(),
                                        *sol->x_NonConst(), *sol->s_NonConst(),
                                        *sol->y_c_NonConst(), *sol->y_d_NonConst(),
                                        true, numberOfEVals);
        }
        assert(retval!=SYMSOLVER_FATAL_ERROR); //TODO make return code
        if (retval==SYMSOLVER_SINGULAR &&
            (rhs.y_c()->Dim()+rhs.y_d()->Dim() > 0) ) {

          // Get new perturbation factors from the perturbation
          // handlers for the singular case
          bool pert_return =
                             perturbHandler_->PerturbForSingularity(delta_x, delta_s,
                                                                    delta_c, delta_d);
          if (!pert_return) {
            Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                           "PerturbForSingularity can't be done\n");
            IpData().TimingStats().PDSystemSolverSolveOnce().End();
            return false;
          }
        }
        else if (retval==SYMSOLVER_WRONG_INERTIA &&
                 augSysSolver_->NumberOfNegEVals() < numberOfEVals) {
          Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                         "Number of negative eigenvalues too small!\n");
          // If the number of negative eigenvalues is too small, then
          // we first try to remedy this by asking for better quality
          // solution (e.g. increasing pivot tolerance), and if that
          // doesn't help, we assume that the system is singular
          bool assume_singular = true;
          if (!augsys_improved_) {
            Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                           "Asking augmented system solver to improve quality of its solutions.\n");
            augsys_improved_ = augSysSolver_->IncreaseQuality();
            if (augsys_improved_) {
              IpData().Append_info_string("q");
              resolve_with_better_quality = true;
              assume_singular = false;
            }
            else {
              Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                             "Quality could not be improved\n");
            }
          }
          if (assume_singular) {
            bool pert_return =
                               perturbHandler_->PerturbForSingularity(delta_x, delta_s,
                                                                      delta_c, delta_d);
            if (!pert_return) {
              Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                             "PerturbForSingularity can't be done for assume singular.\n");
              IpData().TimingStats().PDSystemSolverSolveOnce().End();
              return false;
            }
            IpData().Append_info_string("a");
          }
        }
        else if (retval==SYMSOLVER_WRONG_INERTIA ||
                 retval==SYMSOLVER_SINGULAR) {
          // Get new perturbation factors from the perturbation
          // handlers for the case of wrong inertia
          bool pert_return =
                             perturbHandler_->PerturbForWrongInertia(delta_x, delta_s,
                                                                     delta_c, delta_d);
          if (!pert_return) {
            Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                           "PerturbForWrongInertia can't be done for wrong interia or singular.\n");
            IpData().TimingStats().PDSystemSolverSolveOnce().End();
            return false;
          }
        }
      } // while (retval!=SYMSOLVER_SUCCESS && !fail) {

      // Some output
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Number of trial factorizations performed: %d\n",
                     count);
      Jnlst().Printf(J_DETAILED, J_LINEAR_ALGEBRA,
                     "Perturbation parameters: delta_x=%e delta_s=%e\n                         delta_c=%e delta_d=%e\n",
                     delta_x, delta_s, delta_c, delta_d);
    }

    // Compute the remaining sol Vectors
    Px_L.SinvBlrmZMTdBr(-1., slack_x_L, *rhs.z_L(), z_L, *sol->x_NonConst(), *sol->z_L_NonConst());
    Px_U.SinvBlrmZMTdBr(1., slack_x_U, *rhs.z_U(), z_U, *sol->x_NonConst(), *sol->z_U_NonConst());
    Pd_L.SinvBlrmZMTdBr(-1., slack_s_L, *rhs.v_L(), v_L, *sol->s_NonConst(), *sol->v_L_NonConst());
    Pd_U.SinvBlrmZMTdBr(1., slack_s_U, *rhs.v_U(), v_U, *sol->s_NonConst(), *sol->v_U_NonConst());

    // Finally let's assemble the res result vectors
    res.AddOneVector(alpha, *sol, beta);

    IpData().TimingStats().PDSystemSolverSolveOnce().End();

    return true;
  }

  void PDFullSpaceSolver::ComputeResiduals(
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
    const IteratesVector& res,
    IteratesVector& resid)
  {
    DBG_START_METH("PDFullSpaceSolver::ComputeResiduals", dbg_verbosity);

    DBG_PRINT_VECTOR(2, "res", res);
    IpData().TimingStats().ComputeResiduals().Start();

    // Get the current sizes of the perturbation factors
    Number delta_x;
    Number delta_s;
    Number delta_c;
    Number delta_d;
    perturbHandler_->CurrentPerturbation(delta_x, delta_s, delta_c, delta_d);

    SmartPtr<Vector> tmp;

    // x
    W.MultVector(1., *res.x(), 0., *resid.x_NonConst());
    J_c.TransMultVector(1., *res.y_c(), 1., *resid.x_NonConst());
    J_d.TransMultVector(1., *res.y_d(), 1., *resid.x_NonConst());
    Px_L.MultVector(-1., *res.z_L(), 1., *resid.x_NonConst());
    Px_U.MultVector(1., *res.z_U(), 1., *resid.x_NonConst());
    resid.x_NonConst()->AddTwoVectors(delta_x, *res.x(), -1., *rhs.x(), 1.);

    // s
    Pd_U.MultVector(1., *res.v_U(), 0., *resid.s_NonConst());
    Pd_L.MultVector(-1., *res.v_L(), 1., *resid.s_NonConst());
    resid.s_NonConst()->AddTwoVectors(-1., *res.y_d(), -1., *rhs.s(), 1.);
    if (delta_s!=0.) {
      resid.s_NonConst()->Axpy(delta_s, *res.s());
    }

    // c
    J_c.MultVector(1., *res.x(), 0., *resid.y_c_NonConst());
    resid.y_c_NonConst()->AddTwoVectors(-delta_c, *res.y_c(), -1., *rhs.y_c(), 1.);

    // d
    J_d.MultVector(1., *res.x(), 0., *resid.y_d_NonConst());
    resid.y_d_NonConst()->AddTwoVectors(-1., *res.s(), -1., *rhs.y_d(), 1.);
    if (delta_d!=0.) {
      resid.y_d_NonConst()->Axpy(-delta_d, *res.y_d());
    }

    // zL
    resid.z_L_NonConst()->Copy(*res.z_L());
    resid.z_L_NonConst()->ElementWiseMultiply(slack_x_L);
    tmp = z_L.MakeNew();
    Px_L.TransMultVector(1., *res.x(), 0., *tmp);
    tmp->ElementWiseMultiply(z_L);
    resid.z_L_NonConst()->AddTwoVectors(1., *tmp, -1., *rhs.z_L(), 1.);

    // zU
    resid.z_U_NonConst()->Copy(*res.z_U());
    resid.z_U_NonConst()->ElementWiseMultiply(slack_x_U);
    tmp = z_U.MakeNew();
    Px_U.TransMultVector(1., *res.x(), 0., *tmp);
    tmp->ElementWiseMultiply(z_U);
    resid.z_U_NonConst()->AddTwoVectors(-1., *tmp, -1., *rhs.z_U(), 1.);

    // vL
    resid.v_L_NonConst()->Copy(*res.v_L());
    resid.v_L_NonConst()->ElementWiseMultiply(slack_s_L);
    tmp = v_L.MakeNew();
    Pd_L.TransMultVector(1., *res.s(), 0., *tmp);
    tmp->ElementWiseMultiply(v_L);
    resid.v_L_NonConst()->AddTwoVectors(1., *tmp, -1., *rhs.v_L(), 1.);

    // vU
    resid.v_U_NonConst()->Copy(*res.v_U());
    resid.v_U_NonConst()->ElementWiseMultiply(slack_s_U);
    tmp = v_U.MakeNew();
    Pd_U.TransMultVector(1., *res.s(), 0., *tmp);
    tmp->ElementWiseMultiply(v_U);
    resid.v_U_NonConst()->AddTwoVectors(-1., *tmp, -1., *rhs.v_U(), 1.);

    DBG_PRINT_VECTOR(2, "resid", resid);

    if (Jnlst().ProduceOutput(J_MOREVECTOR, J_LINEAR_ALGEBRA)) {
      resid.Print(Jnlst(), J_MOREVECTOR, J_LINEAR_ALGEBRA, "resid");
    }

    if (Jnlst().ProduceOutput(J_MOREDETAILED, J_LINEAR_ALGEBRA)) {
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_x  %e\n", resid.x()->Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_s  %e\n", resid.s()->Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_c  %e\n", resid.y_c()->Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_d  %e\n", resid.y_d()->Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_zL %e\n", resid.z_L()->Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_zU %e\n", resid.z_U()->Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_vL %e\n", resid.v_L()->Amax());
      Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                     "max-norm resid_vU %e\n", resid.v_U()->Amax());
    }
    IpData().TimingStats().ComputeResiduals().End();
  }

  Number PDFullSpaceSolver::ComputeResidualRatio(const IteratesVector& rhs,
      const IteratesVector& res,
      const IteratesVector& resid)
  {
    DBG_START_METH("PDFullSpaceSolver::ComputeResidualRatio", dbg_verbosity);

    Number nrm_rhs = rhs.Amax();
    Number nrm_res = res.Amax();
    Number nrm_resid = resid.Amax();
    Jnlst().Printf(J_MOREDETAILED, J_LINEAR_ALGEBRA,
                   "nrm_rhs = %8.2e nrm_sol = %8.2e nrm_resid = %8.2e\n",
                   nrm_rhs, nrm_res, nrm_resid);

    if (nrm_rhs+nrm_res == 0.) {
      return nrm_resid;  // this should be zero
    }
    else {
      // ToDo: determine how to include norm of matrix, and what
      // safeguard to use against incredibly large solution vectors
      Number max_cond = 1e6;
      return nrm_resid/(Min(nrm_res, max_cond*nrm_rhs)+nrm_rhs);
    }
  }

} // namespace Ipopt
