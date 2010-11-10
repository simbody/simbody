// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpLimMemQuasiNewtonUpdater.cpp 795 2006-10-11 19:01:37Z andreasw $
//
// Authors:  Andreas Waechter                 IBM    2005-12-26

#include "IpLimMemQuasiNewtonUpdater.hpp"
#include "IpRestoIpoptNLP.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

#include <limits>

namespace Ipopt
{

#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  LimMemQuasiNewtonUpdater::LimMemQuasiNewtonUpdater(
    bool update_for_resto)
      :
      sigma_safe_min_(1e-8),
      sigma_safe_max_(1e+8),
      update_for_resto_(update_for_resto)
  {}

  void LimMemQuasiNewtonUpdater::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->AddLowerBoundedIntegerOption(
      "limited_memory_max_history",
      "Maximum size of the history for the limited quasi-Newton Hessian approximation.",
      0, 6,
      "This option determines the number of most recent iterations that are "
      "taken into account for the limited-memory quasi-Newton approximation.");

    roptions->AddStringOption2(
      "limited_memory_update_type",
      "Quasi-Newton update formula for the limited memory approximation.",
      "bfgs",
      "bfgs", "BFGS update (with skipping)",
      "sr1", "SR1 (not working well)",
      "Determines which update formula is to be used for the limited-memory "
      "quasi-Newton approximation.");

    roptions->AddStringOption3(
      "limited_memory_initialization",
      "Initialization strategy for the limited memory quasi-Newton approximation.",
      "scalar1",
      "scalar1", "sigma = s^Ty/s^Ts",
      "scalar2", "sigma = y^Ty/s^Ty",
      "constant", "sigma = limited_memory_init_val",
      "Determines how the diagonal Matrix B_0 as the first term in the "
      "limited memory approximation should be computed.");

    roptions->AddLowerBoundedNumberOption(
      "limited_memory_init_val",
      "Value for B0 in low-rank update.",
      0, true, 1.,
      "The starting matrix in the low rank update, B0, is chosen to be this "
      "multiple of the identity in the first iteration (when no updates have "
      "been performed yet), and is constantly chosen as this value, if "
      "\"limited_memory_initialization\" is \"constant\".");

    roptions->AddLowerBoundedIntegerOption(
      "limited_memory_max_skipping",
      "Threshold for successive iterations where update is skipped.",
      1, 2,
      "If the update is skipped more than this number of successive "
      "iterations, we quasi-Newton approximation is reset.");
  }

  bool LimMemQuasiNewtonUpdater::InitializeImpl(
    const OptionsList& options,
    const std::string& prefix)
  {
    options.GetIntegerValue("limited_memory_max_history",
                            limited_memory_max_history_, prefix);
    Index enum_int;
    options.GetEnumValue("limited_memory_update_type", enum_int, prefix);
    limited_memory_update_type_ = LMUpdateType(enum_int);
    options.GetEnumValue("limited_memory_initialization", enum_int, prefix);
    limited_memory_initialization_ = LMInitialization(enum_int);
    options.GetNumericValue("limited_memory_init_val",
                            limited_memory_init_val_, prefix);
    options.GetIntegerValue("limited_memory_max_skipping",
                            limited_memory_max_skipping_, prefix);

    h_space_ = NULL;
    curr_lm_memory_ = 0;
    S_ = NULL;
    Y_ = NULL;
    Ypart_ = NULL;
    D_ = NULL;
    L_ = NULL;
    sigma_ = -1;
    V_ = NULL;
    U_ = NULL;
    SdotS_ = NULL;
    SdotS_uptodate_ = false;
    STDRS_ = NULL;
    DRS_ = NULL;
    curr_DR_x_tag_ = 0;

    last_x_ = NULL;
    last_grad_f_ = NULL;
    last_jac_c_ = NULL;
    last_jac_d_ = NULL;
    lm_skipped_iter_ = 0;

    last_eta_ = -1.;

    return true;
  }

  void LimMemQuasiNewtonUpdater::UpdateHessian()
  {
    DBG_START_METH("LimMemQuasiNewtonUpdater::UpdateHessian",
                   dbg_verbosity);

    // First check if this is the first call - it is if the h_space_
    // has not been set yet.
    if (IsNull(h_space_)) {
      if (update_for_resto_) {
        SmartPtr<const SymMatrixSpace> sp = IpNLP().HessianMatrixSpace();
        const CompoundSymMatrixSpace* csp =
          dynamic_cast<const CompoundSymMatrixSpace*> (GetRawPtr(sp));
        DBG_ASSERT(csp);
        h_space_ = dynamic_cast<const LowRankUpdateSymMatrixSpace*>
                   (GetRawPtr(csp->GetCompSpace(0,0)));
      }
      else {
        // ToDo don't need that?!? Can always get the space from the NLP?
        SmartPtr<const SymMatrixSpace> sp = IpNLP().HessianMatrixSpace();
        h_space_ = dynamic_cast<const LowRankUpdateSymMatrixSpace*>(GetRawPtr(sp));
        ASSERT_EXCEPTION(IsValid(h_space_), OPTION_INVALID,
                         "Limited-memory quasi-Newton option chosen, but NLP doesn't provide LowRankUpdateSymMatrixSpace.");
      }
      DBG_ASSERT((h_space_->ReducedDiag() && !update_for_resto_) ||
                 (!h_space_->ReducedDiag() && update_for_resto_));
    }

    SmartPtr<const Matrix> P_LM = h_space_->P_LowRank();
    SmartPtr<const VectorSpace> LM_vecspace =
      h_space_->LowRankVectorSpace();
    DBG_ASSERT(IsValid(LM_vecspace));

    // If we are in the restoration phase, get some additional data
    // for the structured update
    if (update_for_resto_) {
      DBG_ASSERT(IpNLP().objective_depends_on_mu());
      RestoIpoptNLP* resto_nlp = dynamic_cast<RestoIpoptNLP*>
                                 (&IpNLP());
      DBG_ASSERT(resto_nlp);
      curr_DR_x_ = resto_nlp->DR_x();
      DBG_ASSERT(IsValid(curr_DR_x_));
      DBG_ASSERT(curr_DR_x_tag_==0 || curr_DR_x_tag_==curr_DR_x_->GetTag());
      curr_DR_x_tag_ = curr_DR_x_->GetTag();
      if (IsNull(P_LM)) {
        curr_red_DR_x_ = curr_DR_x_;
      }
      else {
        SmartPtr<Vector> tmp = LM_vecspace->MakeNew();
        P_LM->TransMultVector(1, *curr_DR_x_, 0., *tmp);
        curr_red_DR_x_ = ConstPtr(tmp);
      }
      curr_eta_ = resto_nlp->Eta(IpData().curr_mu());
      eta_changed_ = (curr_eta_!=last_eta_);
      Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                     "curr_eta (for B0) is %e\n", curr_eta_);
    }

    // Get the current values
    SmartPtr<const Vector> curr_x;
    SmartPtr<const Vector> curr_grad_f;
    SmartPtr<const Matrix> curr_jac_c;
    SmartPtr<const Matrix> curr_jac_d;
    if (update_for_resto_) {
      const CompoundVector* cv =
        dynamic_cast<const CompoundVector*> (GetRawPtr(IpData().curr()->x()));
      DBG_ASSERT(cv);
      //DBG_PRINT_VECTOR(2, "cv", *cv);
      curr_x = cv->GetComp(0);
      const CompoundMatrix* cm =
        dynamic_cast<const CompoundMatrix*> (GetRawPtr(IpCq().curr_jac_c()));
      curr_jac_c = cm->GetComp(0,0);
      cm = dynamic_cast<const CompoundMatrix*> (GetRawPtr(IpCq().curr_jac_d()));
      curr_jac_d = cm->GetComp(0,0);
    }
    else {
      curr_x = IpData().curr()->x();
      curr_grad_f = IpCq().curr_grad_f();
      curr_jac_c = IpCq().curr_jac_c();
      curr_jac_d = IpCq().curr_jac_d();
    }
    //DBG_PRINT_VECTOR(2, "curr_x", *curr_x);

    // If this is the first iteration, we just gather information and
    // set W to be the identity matrix.
    if (IsNull(last_x_) ||
        lm_skipped_iter_ >= limited_memory_max_skipping_) {
      if (IsNull(last_x_)) {
        DBG_ASSERT(IsNull(last_grad_f_));
        DBG_ASSERT(IsNull(last_jac_c_));
        DBG_ASSERT(IsNull(last_jac_d_));
        Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                       "Limited-Memory approximation started; store data at current iterate.\n");
      }
      else {
        Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                       "Resetting Limited-Memory Update.\n");
        IpData().Append_info_string("Wr");
      }

      last_x_ = curr_x;
      last_grad_f_ = curr_grad_f;
      last_jac_c_ = curr_jac_c;
      last_jac_d_ = curr_jac_d;

      curr_lm_memory_ = 0;
      lm_skipped_iter_ = 0;
      S_ = NULL;
      Y_ = NULL;
      Ypart_ = NULL;
      D_ = NULL;
      L_ = NULL;
      V_ = NULL;
      U_ = NULL;
      SdotS_ = NULL;
      SdotS_uptodate_ = false;
      STDRS_ = NULL;
      DRS_ = NULL;

      if (update_for_resto_) {
        sigma_ = -1;
        last_eta_ = -1.;
      }
      else {
        // Set up W to be multiple of I
        sigma_ = limited_memory_init_val_; // not for resto
      }
      SetW();
      return;
    }

    // Compute s and y for the current iteration

    // s = x_k - x_{k-1}
    SmartPtr<Vector> s_full_new = curr_x->MakeNew();
    DBG_PRINT_VECTOR(2, "last_x", *last_x_);
    DBG_PRINT_VECTOR(2, "curr_x", *curr_x);
    s_full_new->AddTwoVectors(1, *curr_x, -1, *last_x_, 0.);

    // y = grad_lag(x_k,y_k) - grad_lag(x_{k-1},y_k)
    SmartPtr<Vector> y_full_new = curr_x->MakeNew();
    if (update_for_resto_) {
      curr_jac_c->TransMultVector(1., *IpData().curr()->y_c(),
                                  0., *y_full_new);
      curr_jac_d->TransMultVector(1., *IpData().curr()->y_d(),
                                  1., *y_full_new);
    }
    else {
      y_full_new->AddTwoVectors(1., *IpCq().curr_grad_f(),
                                -1, *last_grad_f_, 0.);
      y_full_new->AddTwoVectors(1., *IpCq().curr_jac_cT_times_curr_y_c(),
                                1., *IpCq().curr_jac_dT_times_curr_y_d(), 1.);
    }
    last_jac_c_->TransMultVector(-1., *IpData().curr()->y_c(),
                                 1., *y_full_new);
    last_jac_d_->TransMultVector(-1., *IpData().curr()->y_d(),
                                 1., *y_full_new);

    SmartPtr<Vector> s_new;
    SmartPtr<Vector> y_new;

    if (IsNull(P_LM)) {
      // Then the approximation is in the same space as x
      s_new = s_full_new;
      y_new = y_full_new;
    }
    else {
      s_new = LM_vecspace->MakeNew();
      y_new = LM_vecspace->MakeNew();
      P_LM->TransMultVector(1., *s_full_new, 0., *s_new);
      P_LM->TransMultVector(1., *y_full_new, 0., *y_new);
    }
    s_full_new = NULL;
    y_full_new = NULL;
    DBG_PRINT_VECTOR(2, "s_new", *s_new);
    DBG_PRINT_VECTOR(2, "ypart_new", *y_new);

    // In the restoration phase case, y_new is only the y without the
    // objective part, so now we add the explicitly known objective
    // part
    SmartPtr<Vector> ypart_new;
    if (update_for_resto_) {
      ypart_new = y_new;
      y_new = ypart_new->MakeNew();
      y_new->AddOneVector(curr_eta_, *s_new, 0.);
      y_new->ElementWiseMultiply(*curr_red_DR_x_);
      y_new->AddOneVector(1., *ypart_new, 1.);
    }

    bool skipping = false;
    bool retroactive_skip = false;

    // Sometimes the change in the primal variables is very small, and
    // we should then skip the update...
    // ToDo: Find good number or do relative test?
    Number s_new_max = s_new->Amax();
    Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                   "In limited-memory update, s_new_max is %e\n", s_new_max);
    if (s_new_max<100.*std::numeric_limits<Number>::epsilon()) {
      skipping = true;
      IpData().Append_info_string("WS");
    }

    if (!skipping) {
      switch (limited_memory_update_type_) {
        case BFGS: {
          skipping = CheckSkippingBFGS(*s_new, *y_new);
          DBG_PRINT_VECTOR(2, "y_new", *y_new);
          if (skipping) {
            break;
          }

          if (update_for_resto_) {
            // In the restoration phase we don't know if the update has
            // to be skipped yet, since a change in curr_eta_ might
            // cause some of the s^Ty pairs to become negative.  For
            // now, we then just skip the update and use what we had
            // before.
            //
            // ToDo: It would probably be better just to reset the
            // entire update.
            StoreInternalDataBackup();
          }

          bool augment_memory = UpdateInternalData(*s_new, *y_new, ypart_new);

          if (update_for_resto_) {
            Number dmin = D_->Min();
            if (dmin <= 0.) {
              Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                             "Skipping BFGS update in restoration phase because Dmin is %e\n", dmin);
              // Check if any of the s^Ty pairs are non-positive.  If so, skip
              IpData().Append_info_string("We");
              skipping = true;
              RestoreInternalDataBackup();
              break;
            }
          }

          Number sTy_new = s_new->Dot(*y_new);
          if (!update_for_resto_) {
            // Compute the initial matrix B_0
            switch (limited_memory_initialization_) {
              case SCALAR1:
              sigma_ = sTy_new/pow(s_new->Nrm2(),2);
              break;
              case SCALAR2:
              sigma_ = pow(y_new->Nrm2(),2)/sTy_new;
              break;
              case CONSTANT:
              sigma_ = limited_memory_init_val_;
              break;
            }
            Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                           "sigma (for B0) is %e\n", sigma_);
            if (sigma_ < sigma_safe_min_ ||
                sigma_ > sigma_safe_max_) {
              sigma_ = Max(Min(sigma_safe_max_, sigma_), sigma_safe_min_);
              IpData().Append_info_string("Wp");
              Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                             "Projecting sigma into safeguards to be %e!\n", sigma_);
            }
          }

          if (limited_memory_max_history_ == 0 ) {
            break;
          }

          // First update V - here only the last column is updated
          DBG_ASSERT(sTy_new > 0.);
          SmartPtr<Vector> v_new = y_new->MakeNewCopy();
          v_new->Scal(1./sqrt(sTy_new));
          if (augment_memory) {
            AugmentMultiVector(V_, *v_new);
          }
          else {
            ShiftMultiVector(V_, *v_new);
          }

          // Compute Ltilde = L * diag(D^{-1/2});
          SmartPtr<DenseVector> Dtilde = D_->MakeNewDenseVector();
          Dtilde->Copy(*D_);
          Dtilde->ElementWiseSqrt();
          Dtilde->ElementWiseReciprocal();
          SmartPtr<DenseGenMatrix> Ltilde = L_->MakeNewDenseGenMatrix();
          DBG_PRINT_MATRIX(3, "D", *D_);
          DBG_PRINT_MATRIX(3, "L", *L_);
          Ltilde->Copy(*L_);
          Ltilde->ScaleColumns(*Dtilde);
          DBG_PRINT_MATRIX(3, "Ltilde", *Ltilde);

          // M = Ltilde * Ltilde^T
          SmartPtr<DenseSymMatrixSpace> Mspace =
            new DenseSymMatrixSpace(curr_lm_memory_);
          SmartPtr<DenseSymMatrix> M = Mspace->MakeNewDenseSymMatrix();
          M->HighRankUpdate(false, 1., *Ltilde, 0.);

          // M += S^T B_0 S
          if (!update_for_resto_) {
            // For now, we assume that B_0 is sigma*I
            DBG_ASSERT(SdotS_uptodate_);
            DBG_PRINT_MATRIX(3, "SdotS", *SdotS_);
            M->AddMatrix(sigma_, *SdotS_, 1.);
          }
          else {
            DBG_PRINT_MATRIX(3, "STDRS", *STDRS_);
            M->AddMatrix(curr_eta_, *STDRS_, 1.);
          }

          // Compute Cholesky factor J with M = J J^T
          DBG_PRINT_MATRIX(3, "M", *M);
          SmartPtr<DenseGenMatrix> J = L_->MakeNewDenseGenMatrix();
          bool cholesky_retval = J->ComputeCholeskyFactor(*M);
          DBG_PRINT_MATRIX(3, "J", *J);
          if (!cholesky_retval) {
            Jnlst().Printf(J_WARNING, J_HESSIAN_APPROXIMATION,
                           "Cholesky factorization failed for LBFGS update! Skipping update.\n");
            skipping = true;
            break;
          }

          // Compute C = J^{-T}
          SmartPtr<DenseGenMatrix> C = J->MakeNewDenseGenMatrix();
          C->FillIdentity();
          J->CholeskyBackSolveMatrix(true, 1., *C);

          // Compute U = B_0 * S * C
          U_ = S_->MakeNewMultiVectorMatrix();
          if (!update_for_resto_) {
            DBG_ASSERT(sigma_>0.);
            U_->AddRightMultMatrix(sigma_, *S_, *C, 0.);
          }
          else {
            DBG_ASSERT(sigma_<0.);
            U_->AddRightMultMatrix(curr_eta_, *DRS_, *C, 0.);
          }

          // Compute Lbar = Ltilde^T * C
          SmartPtr<DenseGenMatrix> Lbar = Ltilde->MakeNewDenseGenMatrix();
          Lbar->AddMatrixProduct(1., *Ltilde, true, *C, false, 0.);

          // Compute U += V * Lbar;
          U_->AddRightMultMatrix(1., *V_, *Lbar, 1.);
          break;
        }
        case SR1:
        // TODO IMPLEMENT WELL!
        if (IpData().info_regu_x()>0.) {
          RestoreInternalDataBackup();
          Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                         "Undoing most recent SR1 update.\n");
          IpData().Append_info_string("Wb");
          retroactive_skip = true;
        }

        // We don't know if the update has to be skipped, so we store a
        // backup of all internal data that can be restored in case we
        // have to skip the update.
        StoreInternalDataBackup();

        // Update the internal stuff
        UpdateInternalData(*s_new, *y_new, ypart_new);

        if (!update_for_resto_) {
          // Set B0 for now as we do for BFGS - except that we take the
          // abs value?
          //
          // It seems that it is not a good idea to use that update if
          // this is the first contribution to the limited-memory history,
          // since otherwise the update will be skipped (and then all
          // updates will be skipped)
          if (curr_lm_memory_==1) {
            sigma_ = limited_memory_init_val_;
          }
          else {
            // ToDo: What lower bound to use?
            Number sTy_new = Max(1e-8, fabs(s_new->Dot(*y_new)));
            DBG_ASSERT(sTy_new!=0.);
            switch (limited_memory_initialization_) {
              case SCALAR1:
              sigma_ = sTy_new/pow(s_new->Nrm2(),2);
              break;
              case SCALAR2:
              sigma_ = pow(y_new->Nrm2(),2)/sTy_new;
              break;
              case CONSTANT:
              sigma_ = limited_memory_init_val_;
              break;
            }
          }
          Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                         "sigma (for B0) is %e\n", sigma_);
          // ToDo Decide what to use for SR1 - or at least use different
          // skipping rule
        }

        if (limited_memory_max_history_ == 0 ) {
          break;
        }

        DBG_PRINT_VECTOR(2,"S",*S_);
        DBG_PRINT_VECTOR(2,"Y",*Y_);

        // Compute Z as D + L + L^T - S^TB_0S
        SmartPtr<DenseSymMatrix> Z;
        if (!update_for_resto_) {
          Z = SdotS_->MakeNewDenseSymMatrix();
          DBG_PRINT_MATRIX(3, "SdotS", *SdotS_);
          DBG_PRINT((1, "sigma_ = %e\n", sigma_));
          Z->AddMatrix(-sigma_, *SdotS_, 0.);
        }
        else {
          Z = STDRS_->MakeNewDenseSymMatrix();
          Z->AddMatrix(-curr_eta_, *STDRS_, 0.);
        }

        DBG_PRINT_MATRIX(3, "L", *L_);
        DBG_PRINT_VECTOR(3, "D", *D_);
        Z->SpecialAddForLMSR1(*D_, *L_);
        // Compute the eigenvectors Q and eignevalues E for Z
        SmartPtr<DenseGenMatrix> Q = L_->MakeNewDenseGenMatrix();
        SmartPtr<DenseVector> E = D_->MakeNewDenseVector();
        DBG_PRINT_MATRIX(3, "Z", *Z);
        bool retval = Q->ComputeEigenVectors(*Z, *E);
        ASSERT_EXCEPTION(retval, INTERNAL_ABORT,
                         "Eigenvalue decomposition failed for limited-memory SR1 update.");
        DBG_PRINT_VECTOR(2, "E", *E);
        DBG_PRINT_MATRIX(3, "Q", *Q);
        SmartPtr<DenseGenMatrix> Qminus;
        SmartPtr<DenseGenMatrix> Qplus;
        // Split the eigenvectors and scale them
        skipping = SplitEigenvalues(*Q, *E, Qminus, Qplus);
        if (skipping) {
          RestoreInternalDataBackup();
          break;
        }

        // Compute Vtilde = Y - B_0*S
        SmartPtr<MultiVectorMatrix> Vtilde = Y_->MakeNewMultiVectorMatrix();
        Vtilde->AddOneMultiVectorMatrix(1., *Y_, 0.);
        if (!update_for_resto_) {
          Vtilde->AddOneMultiVectorMatrix(-sigma_, *S_, 1.);
        }
        else {
          DBG_ASSERT(sigma_<0.);
          Vtilde->AddOneMultiVectorMatrix(-curr_eta_, *DRS_, 1.);
        }

        // Now get U as the Vtilde * Qminus
        if (IsValid(Qminus)) {
          SmartPtr<MultiVectorMatrixSpace> U_space =
            new MultiVectorMatrixSpace(Qminus->NCols(), *s_new->OwnerSpace());
          U_ = U_space->MakeNewMultiVectorMatrix();
          U_->AddRightMultMatrix(1., *Vtilde, *Qminus, 0.);
          DBG_PRINT_MATRIX(3, "U", *U_);
        }
        else {
          U_ = NULL;
        }

        // Now get V as the Vtilde * Qplus
        if (IsValid(Qplus)) {
          SmartPtr<MultiVectorMatrixSpace> V_space =
            new MultiVectorMatrixSpace(Qplus->NCols(), *s_new->OwnerSpace());
          V_ = V_space->MakeNewMultiVectorMatrix();
          V_->AddRightMultMatrix(1., *Vtilde, *Qplus, 0.);
          DBG_PRINT_MATRIX(3, "V", *V_);
        }
        else {
          V_ = NULL;
        }
        break;
      }
    }

    if (!skipping) {
      // Put together W
      SetW();
      if (retroactive_skip) {
        lm_skipped_iter_++;
      }
      else {
        lm_skipped_iter_ = 0;
      }
    }
    else { // if (!skipping) {
      IpData().Append_info_string("Ws");
      lm_skipped_iter_++;
    }
    Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                   "Number of successive iterations with skipping: %d\n",
                   lm_skipped_iter_);

    // Keep stuff around in case we want to skip SR1 retroactively
    // because of negative curvature!
    //
    //ReleaseInternalDataBackup();

    last_x_ = curr_x;
    last_grad_f_ = curr_grad_f;
    last_jac_c_ = curr_jac_c;
    last_jac_d_ = curr_jac_d;

    if (update_for_resto_) {
      last_eta_ = curr_eta_;
      curr_DR_x_ = NULL;
    }
  }

  void LimMemQuasiNewtonUpdater::
  StoreInternalDataBackup()
  {
    DBG_START_METH("LimMemQuasiNewtonUpdater::StoreInternalDataBackup",
                   dbg_verbosity);
    curr_lm_memory_old_ = curr_lm_memory_;
    S_old_ = S_;
    Y_old_ = Y_;
    Ypart_old_ = Ypart_;
    D_old_ = D_;
    L_old_ = L_;
    SdotS_old_ = SdotS_;
    SdotS_uptodate_old_ = SdotS_uptodate_;
    STDRS_old_ = STDRS_;
    DRS_old_ = DRS_;
    sigma_old_ = sigma_;
    V_old_ = V_;
    U_old_ = U_;
  }

  void LimMemQuasiNewtonUpdater::
  RestoreInternalDataBackup()
  {
    DBG_START_METH("LimMemQuasiNewtonUpdater::RestoreInternalDataBackup",
                   dbg_verbosity);
    curr_lm_memory_ = curr_lm_memory_old_;
    S_ = S_old_;
    Y_ = Y_old_;
    Ypart_ = Ypart_old_;
    D_ = D_old_;
    L_ = L_old_;
    SdotS_ = SdotS_old_;
    SdotS_uptodate_ = SdotS_uptodate_old_;
    STDRS_ = STDRS_old_;
    DRS_ = DRS_old_;
    sigma_ = sigma_old_;
    V_ = V_old_;
    U_ = U_old_;
  }

  void LimMemQuasiNewtonUpdater::
  ReleaseInternalDataBackup()
  {
    DBG_START_METH("LimMemQuasiNewtonUpdater::ReleaseInternalDataBackup",
                   dbg_verbosity);
    S_old_ = NULL;
    Y_old_ = NULL;
    Ypart_old_ = NULL;
    D_old_ = NULL;
    L_old_ = NULL;
    SdotS_old_ = NULL;
    SdotS_uptodate_old_ = false;
    STDRS_old_ = NULL;
    DRS_old_ = NULL;
    V_old_ = NULL;
    U_old_ = NULL;
  }

  bool LimMemQuasiNewtonUpdater::
  UpdateInternalData(const Vector& s_new, const Vector& y_new,
                     SmartPtr<Vector> ypart_new)
  {
    DBG_START_METH("LimMemQuasiNewtonUpdater::UpdateInternalData",
                   dbg_verbosity);

    if (limited_memory_max_history_==0) {
      return false;
    }

    bool augment_memory;
    if (curr_lm_memory_ < limited_memory_max_history_) {
      augment_memory = true;
      curr_lm_memory_++;
    }
    else {
      augment_memory = false;
    }

    if (!update_for_resto_) {
      // Update the internal information
      if (augment_memory) {
        // If the memory is still
        // growing, increase the vector spaces etc
        AugmentMultiVector(S_, s_new);
        AugmentMultiVector(Y_, y_new);
        AugmentDenseVector(D_, s_new.Dot(y_new));
        AugmentLMatrix(L_, *S_, *Y_);
        DBG_ASSERT(SdotS_uptodate_ || S_->NCols()==1);
        AugmentSdotSMatrix(SdotS_, *S_);
        SdotS_uptodate_ = true;
      }
      else {
        // Otherwise, we shift the internal data
        ShiftMultiVector(S_, s_new);
        ShiftMultiVector(Y_, y_new);
        ShiftDenseVector(D_, s_new.Dot(y_new));
        ShiftLMatrix(L_, *S_, *Y_);
        DBG_ASSERT(SdotS_uptodate_);
        ShiftSdotSMatrix(SdotS_, *S_);
      }
    }
    else {
      // Compute DR*s_new;
      SmartPtr<Vector> DRs_new = s_new.MakeNewCopy();
      DRs_new->ElementWiseMultiply(*curr_red_DR_x_);
      if (augment_memory) {
        AugmentMultiVector(S_, s_new);
        AugmentMultiVector(DRS_, *DRs_new);
        AugmentMultiVector(Ypart_, *ypart_new);
        AugmentSTDRSMatrix(STDRS_, *S_, *DRS_);
      }
      else {
        ShiftMultiVector(S_, s_new);
        ShiftMultiVector(DRS_, *DRs_new);
        ShiftMultiVector(Ypart_, *ypart_new);
        ShiftSTDRSMatrix(STDRS_, *S_, *DRS_);
      }
      DBG_PRINT((1,"curr_eta = %e\n", curr_eta_));
      DBG_PRINT_VECTOR(2,"curr_red_DR_x", *curr_red_DR_x_);
      DBG_PRINT_VECTOR(2,"S", *S_);
      DBG_PRINT_VECTOR(2,"Ypart", *Ypart_);
      RecalcY(curr_eta_, *curr_red_DR_x_, *S_, *Ypart_, Y_);
      DBG_PRINT_VECTOR(2,"Y", *Y_);
      RecalcD(*S_, *Y_, D_);
      RecalcL(*S_, *Y_, L_);
    }

    return augment_memory;
  }

  bool LimMemQuasiNewtonUpdater::
  SplitEigenvalues(DenseGenMatrix& Q, const DenseVector& E,
                   SmartPtr<DenseGenMatrix>& Qminus,
                   SmartPtr<DenseGenMatrix>& Qplus)
  {
    DBG_START_METH("LimMemQuasiNewtonUpdater::SplitEigenvalues",
                   dbg_verbosity);

    Index dim = E.Dim();
    DBG_ASSERT(dim==Q.NCols());
    DBG_ASSERT(dim==Q.NRows());

    const Number* Evals = E.Values();
    const Number* Qvals = Q.Values();

    // Determine number of negative eigenvalues
    Index nneg=0;
    for (Index i=0; i<dim; i++) {
      if (Evals[i]<0.) {
        nneg++;
      }
    }

    // Determine the ratio of smallest over the largest eigenvalue
    Number emax = Max(fabs(Evals[0]), fabs(Evals[dim-1]));
    if (emax==0.) {
      return true;
    }
    Number emin;
    if (nneg==0) {
      emin = Evals[0];
    }
    else if (nneg==dim) {
      emin = -Evals[dim-1];
    }
    else {
      emin = Min(-Evals[nneg-1],Evals[nneg]);
    }
    Number ratio = emin/emax;
    Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                   "Eigenvalues in SR1 update: emin=%e emax=%e ratio=%e\n",
                   emin, emax, ratio);
    DBG_ASSERT(ratio>=0.);
    // ToDo make the following an option?
    const Number tol = 1e-12;
    if (ratio<tol) {
      return true;
    }

    // Consider the special cases where there are only positive or
    // only negative eigenvalues
    if (nneg==0) {
      SmartPtr<DenseVector> tmp = E.MakeNewDenseVector();
      tmp->Copy(E);
      tmp->ElementWiseSqrt();
      tmp->ElementWiseReciprocal();
      Q.ScaleColumns(*tmp);
      Qplus = &Q;
      Qminus = NULL;
      return false;
    }
    else if (nneg==E.Dim()) {
      SmartPtr<DenseVector> tmp = E.MakeNewDenseVector();
      tmp->AddOneVector(-1., E, 0.);
      tmp->ElementWiseSqrt();
      tmp->ElementWiseReciprocal();
      Q.ScaleColumns(*tmp);
      Qminus = &Q;
      Qplus =  NULL;
      return false;
    }

    // Create Qminus
    SmartPtr<DenseGenMatrixSpace> Qminus_space =
      new DenseGenMatrixSpace(dim, nneg);
    Qminus = Qminus_space->MakeNewDenseGenMatrix();
    Number* Qminus_vals = Qminus->Values();
    for (Index j=0; j<nneg; j++) {
      Number esqrt = sqrt(-Evals[j]);
      for (Index i=0; i<dim; i++) {
        Qminus_vals[i+j*dim] = Qvals[i+j*dim]/esqrt;
      }
    }

    // Create Qplus
    SmartPtr<DenseGenMatrixSpace> Qplus_space =
      new DenseGenMatrixSpace(dim, dim-nneg);
    Qplus = Qplus_space->MakeNewDenseGenMatrix();
    Number* Qplus_vals = Qplus->Values();
    for (Index j=0; j<dim-nneg; j++) {
      DBG_ASSERT(Evals[j+nneg]>0.);
      Number esqrt = sqrt(Evals[j+nneg]);
      for (Index i=0; i<dim; i++) {
        Qplus_vals[i+j*dim] = Qvals[i+(j+nneg)*dim]/esqrt;
      }
    }

    return false;
  }

  bool LimMemQuasiNewtonUpdater::
  CheckSkippingBFGS(Vector& s_new, Vector& y_new)
  {
    Number sTy = s_new.Dot(y_new);
    Number snrm = s_new.Nrm2();
    Number ynrm = y_new.Nrm2();

    // ToDo make a parameter?
    Number tol = sqrt(std::numeric_limits<Number>::epsilon());

    Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                   "Limited-Memory test for skipping:\n");
    Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                   "     s^Ty = %e snrm = %e ynrm = %e\n",
                   sTy, snrm, ynrm);

    bool skipping;

    DBG_ASSERT(limited_memory_update_type_==BFGS);
    skipping = (sTy <= tol*snrm*ynrm);

    if (skipping) {
      Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                     "     Skip the update.\n");
    }
    else {
      Jnlst().Printf(J_DETAILED, J_HESSIAN_APPROXIMATION,
                     "     Perform the update.\n");
    }

    return skipping;
  }

  void LimMemQuasiNewtonUpdater::
  AugmentMultiVector(SmartPtr<MultiVectorMatrix>& V,
                     const Vector& v_new)
  {
    Index ncols;
    if (IsValid(V)) {
      ncols = V->NCols();
    }
    else {
      ncols = 0;
    }

    SmartPtr<const VectorSpace> vec_space = v_new.OwnerSpace();
    SmartPtr<MultiVectorMatrixSpace> new_Vspace =
      new MultiVectorMatrixSpace(ncols+1, *vec_space);
    SmartPtr<MultiVectorMatrix> new_V =
      new_Vspace->MakeNewMultiVectorMatrix();
    for (Index i=0; i<ncols; i++) {
      new_V->SetVector(i, *V->GetVector(i));
    }
    new_V->SetVector(ncols, v_new);

    V = new_V;
  }

  void LimMemQuasiNewtonUpdater::
  AugmentDenseVector(SmartPtr<DenseVector>& V,
                     Number v_new)
  {
    Index ndim;
    if (IsValid(V)) {
      ndim = V->Dim();
    }
    else {
      ndim = 0;
    }

    SmartPtr<DenseVectorSpace> new_Vspace =
      new DenseVectorSpace(ndim+1);
    SmartPtr<DenseVector> new_V =
      new_Vspace->MakeNewDenseVector();
    Number* newVvalues = new_V->Values();
    if (IsValid(V)) {
      DBG_ASSERT(!V->IsHomogeneous());
      const Number* Vvalues = V->Values();
      for (Index i=0; i<ndim; i++) {
        newVvalues[i] = Vvalues[i];
      }
    }
    newVvalues[ndim] = v_new;

    V = new_V;
  }

  void LimMemQuasiNewtonUpdater::
  AugmentLMatrix(SmartPtr<DenseGenMatrix>& V,
                 const MultiVectorMatrix& S,
                 const MultiVectorMatrix& Y)
  {
    Index ndim;
    if (IsValid(V)) {
      ndim = V->NCols();
      DBG_ASSERT(ndim==V->NRows());
    }
    else {
      ndim = 0;
    }
    DBG_ASSERT(S.NCols()==ndim+1);
    DBG_ASSERT(Y.NCols()==ndim+1);

    SmartPtr<DenseGenMatrixSpace> new_Vspace =
      new DenseGenMatrixSpace(ndim+1, ndim+1);
    SmartPtr<DenseGenMatrix> new_V =
      new_Vspace->MakeNewDenseGenMatrix();
    Number* newVvalues = new_V->Values();
    if (IsValid(V)) {
      const Number* Vvalues = V->Values();
      for (Index j=0; j<ndim; j++) {
        for (Index i=0; i<ndim; i++) {
          newVvalues[i+j*(ndim+1)] = Vvalues[i+j*ndim];
        }
      }
    }

    for (Index j=0; j<ndim; j++) {
      newVvalues[ndim + j*(ndim+1)] =
        S.GetVector(ndim)->Dot(*Y.GetVector(j));
    }

    for (Index i=0; i<ndim+1; i++) {
      newVvalues[i + ndim*(ndim+1)] = 0.;
    }

    V = new_V;
  }

  void LimMemQuasiNewtonUpdater::
  AugmentSdotSMatrix(SmartPtr<DenseSymMatrix>& V,
                     const MultiVectorMatrix& S)
  {
    Index ndim;
    if (IsValid(V)) {
      ndim = V->Dim();
    }
    else {
      ndim = 0;
    }
    DBG_ASSERT(S.NCols()==ndim+1);

    SmartPtr<DenseSymMatrixSpace> new_Vspace =
      new DenseSymMatrixSpace(ndim+1);
    SmartPtr<DenseSymMatrix> new_V =
      new_Vspace->MakeNewDenseSymMatrix();
    Number* newVvalues = new_V->Values();
    if (IsValid(V)) {
      const Number* Vvalues = V->Values();
      for (Index j=0; j<ndim; j++) {
        for (Index i=j; i<ndim; i++) {
          newVvalues[i+j*(ndim+1)] = Vvalues[i+j*ndim];
        }
      }
    }

    for (Index j=0; j<ndim+1; j++) {
      newVvalues[ndim + j*(ndim+1)] =
        S.GetVector(ndim)->Dot(*S.GetVector(j));
    }

    V = new_V;
  }

  void LimMemQuasiNewtonUpdater::
  AugmentSTDRSMatrix(SmartPtr<DenseSymMatrix>& V,
                     const MultiVectorMatrix& S,
                     const MultiVectorMatrix& DRS)
  {
    Index ndim;
    if (IsValid(V)) {
      ndim = V->Dim();
    }
    else {
      ndim = 0;
    }
    DBG_ASSERT(S.NCols()==ndim+1);

    SmartPtr<DenseSymMatrixSpace> new_Vspace =
      new DenseSymMatrixSpace(ndim+1);
    SmartPtr<DenseSymMatrix> new_V =
      new_Vspace->MakeNewDenseSymMatrix();
    Number* newVvalues = new_V->Values();
    if (IsValid(V)) {
      const Number* Vvalues = V->Values();
      for (Index j=0; j<ndim; j++) {
        for (Index i=j; i<ndim; i++) {
          newVvalues[i+j*(ndim+1)] = Vvalues[i+j*ndim];
        }
      }
    }

    for (Index j=0; j<ndim+1; j++) {
      newVvalues[ndim + j*(ndim+1)] =
        S.GetVector(ndim)->Dot(*DRS.GetVector(j));
    }

    V = new_V;
  }

  void LimMemQuasiNewtonUpdater::
  ShiftMultiVector(SmartPtr<MultiVectorMatrix>& V, const Vector& v_new)
  {
    Index ncols = V->NCols();

    SmartPtr<MultiVectorMatrix> new_V = V->MakeNewMultiVectorMatrix();

    for (Index i=0; i<ncols-1; i++) {
      new_V->SetVector(i, *V->GetVector(i+1));
    }
    new_V->SetVector(ncols-1, v_new);

    V = new_V;
  }

  void LimMemQuasiNewtonUpdater::
  ShiftDenseVector(SmartPtr<DenseVector>& V, Number v_new)
  {
    Index ndim = V->Dim();

    SmartPtr<DenseVector> new_V = V->MakeNewDenseVector();

    DBG_ASSERT(!V->IsHomogeneous());
    Number* Vvalues = V->Values();
    Number* new_Vvalues = new_V->Values();
    for (Index i=0; i<ndim-1; i++) {
      new_Vvalues[i] = Vvalues[i+1];
    }
    new_Vvalues[ndim-1] = v_new;

    V = new_V;
  }

  void LimMemQuasiNewtonUpdater::
  ShiftLMatrix(SmartPtr<DenseGenMatrix>& V,
               const MultiVectorMatrix& S,
               const MultiVectorMatrix& Y)
  {
    Index ndim = V->NCols();
    DBG_ASSERT(ndim==V->NRows());
    DBG_ASSERT(S.NCols()==ndim);
    DBG_ASSERT(Y.NCols()==ndim);

    SmartPtr<DenseGenMatrix> new_V = V->MakeNewDenseGenMatrix();

    Number* Vvalues = V->Values();
    Number* new_Vvalues = new_V->Values();
    for (Index j=0; j<ndim-1; j++) {
      for (Index i=0; i<ndim-1; i++) {
        new_Vvalues[i+j*ndim] = Vvalues[i+1+(j+1)*ndim];
      }
    }

    for (Index j=0; j<ndim-1; j++) {
      new_Vvalues[ndim-1 + j*ndim] =
        S.GetVector(ndim-1)->Dot(*Y.GetVector(j));
    }

    for (Index i=0; i<ndim; i++) {
      new_Vvalues[i + ndim*(ndim-1)] = 0.;
    }

    V = new_V;
  }

  void LimMemQuasiNewtonUpdater::
  ShiftSdotSMatrix(SmartPtr<DenseSymMatrix>& V,
                   const MultiVectorMatrix& S)
  {
    Index ndim = V->Dim();
    DBG_ASSERT(S.NCols()==ndim);

    SmartPtr<DenseSymMatrix> new_V = V->MakeNewDenseSymMatrix();

    Number* Vvalues = V->Values();
    Number* new_Vvalues = new_V->Values();
    for (Index j=0; j<ndim-1; j++) {
      for (Index i=j; i<ndim-1; i++) {
        new_Vvalues[i+j*ndim] = Vvalues[i+1+(j+1)*ndim];
      }
    }

    for (Index j=0; j<ndim; j++) {
      new_Vvalues[ndim-1 + j*ndim] =
        S.GetVector(ndim-1)->Dot(*S.GetVector(j));
    }

    V = new_V;
  }

  void LimMemQuasiNewtonUpdater::
  ShiftSTDRSMatrix(SmartPtr<DenseSymMatrix>& V,
                   const MultiVectorMatrix& S,
                   const MultiVectorMatrix& DRS)
  {
    Index ndim = V->Dim();
    DBG_ASSERT(S.NCols()==ndim);

    SmartPtr<DenseSymMatrix> new_V = V->MakeNewDenseSymMatrix();

    Number* Vvalues = V->Values();
    Number* new_Vvalues = new_V->Values();
    for (Index j=0; j<ndim-1; j++) {
      for (Index i=j; i<ndim-1; i++) {
        new_Vvalues[i+j*ndim] = Vvalues[i+1+(j+1)*ndim];
      }
    }

    for (Index j=0; j<ndim; j++) {
      new_Vvalues[ndim-1 + j*ndim] =
        S.GetVector(ndim-1)->Dot(*DRS.GetVector(j));
    }

    V = new_V;
  }

  void LimMemQuasiNewtonUpdater::SetW()
  {
    DBG_START_METH("LimMemQuasiNewtonUpdater::SetW",
                   dbg_verbosity);

    SmartPtr<Vector> B0;
    if (update_for_resto_) {
      B0 = curr_DR_x_->MakeNew();
      B0->AddOneVector(curr_eta_, *curr_DR_x_, 0.);
    }
    else {
      SmartPtr<const VectorSpace> LM_vecspace =
        h_space_->LowRankVectorSpace();
      DBG_ASSERT(IsValid(LM_vecspace));
      B0 = LM_vecspace->MakeNew();
      B0->Set(sigma_);
    }
    DBG_PRINT_VECTOR(2, "B0", *B0);

    SmartPtr<LowRankUpdateSymMatrix> W =
      h_space_->MakeNewLowRankUpdateSymMatrix();
    W->SetDiag(*B0);
    if (IsValid(V_)) {
      W->SetV(*V_);
    }
    if (IsValid(U_)) {
      W->SetU(*U_);
    }
    if (update_for_resto_) {
      SmartPtr<const SymMatrixSpace> sp = IpNLP().HessianMatrixSpace();
      const CompoundSymMatrixSpace* csp =
        dynamic_cast<const CompoundSymMatrixSpace*> (GetRawPtr(sp));
      SmartPtr<CompoundSymMatrix> CW =
        csp->MakeNewCompoundSymMatrix();
      CW->SetComp(0,0,*W);
      IpData().Set_W(GetRawPtr(CW));
    }
    else {
      IpData().Set_W(GetRawPtr(W));
    }

#ifdef PRINT_W
    // DELETEME
    const DenseVector* dx = dynamic_cast<const DenseVector*>
                            (GetRawPtr(IpData().curr()->x()));
    DBG_ASSERT(dx);
    SmartPtr<DenseVector> tmpx = dx->MakeNewDenseVector();
    SmartPtr<DenseVector> tmpy = dx->MakeNewDenseVector();
    for (Index i=0; i<dx->Dim(); i++) {
      Number* tmpx_vals = tmpx->Values();
      for (Index j=0; j<dx->Dim(); j++) {
        tmpx_vals[j] = 0.;
      }
      tmpx_vals[i] = 1.;
      W->MultVector(1., *tmpx, 0., *tmpy);
      tmpx->Print(Jnlst(), J_DETAILED, J_MAIN, "tmpx");
      tmpy->Print(Jnlst(), J_DETAILED, J_MAIN, "tmpy");
    }
    // ENDDELETEME
#endif

  }

  void LimMemQuasiNewtonUpdater::RecalcY(Number eta, const Vector& DR_x,
                                         MultiVectorMatrix& DRS,
                                         MultiVectorMatrix& Ypart,
                                         SmartPtr<MultiVectorMatrix>& Y)
  {
    SmartPtr<const MultiVectorMatrixSpace> mvspace =
      Ypart.MultiVectorMatrixOwnerSpace();
    Y = mvspace->MakeNewMultiVectorMatrix();
    Y->AddOneMultiVectorMatrix(eta, DRS, 0.);
    Y->AddOneMultiVectorMatrix(1., Ypart, 1.);
  }

  void LimMemQuasiNewtonUpdater::RecalcD(MultiVectorMatrix& S,
                                         MultiVectorMatrix& Y,
                                         SmartPtr<DenseVector>& D)
  {
    SmartPtr<DenseVectorSpace> space =
      new DenseVectorSpace(S.NCols());
    D = space->MakeNewDenseVector();
    Number* Dvalues = D->Values();
    for (Index i=0; i<S.NCols(); i++) {
      Dvalues[i] = S.GetVector(i)->Dot(*Y.GetVector(i));
    }
  }

  void LimMemQuasiNewtonUpdater::RecalcL(MultiVectorMatrix& S,
                                         MultiVectorMatrix& Y,
                                         SmartPtr<DenseGenMatrix>& L)
  {
    Index dim = S.NCols();
    SmartPtr<DenseGenMatrixSpace> space =
      new DenseGenMatrixSpace(dim, dim);
    L = space->MakeNewDenseGenMatrix();
    Number* Lvalues = L->Values();
    for (Index j=0; j<dim; j++) {
      for (Index i=0; i<=j; i++) {
        Lvalues[i+j*dim] = 0.;
      }
      for (Index i=j+1; i<dim; i++) {
        Lvalues[i+j*dim] = S.GetVector(i)->Dot(*Y.GetVector(j));
        ;
      }
    }
  }

} // namespace Ipopt
