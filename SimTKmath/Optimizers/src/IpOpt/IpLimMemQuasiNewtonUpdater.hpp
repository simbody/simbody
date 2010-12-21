// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpLimMemQuasiNewtonUpdater.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Andreas Waechter            IBM    2005-12-26

#ifndef __IPLIMMEMQUASINEWTONUPDATER_HPP__
#define __IPLIMMEMQUASINEWTONUPDATER_HPP__

#include "IpHessianUpdater.hpp"
#include "IpLowRankUpdateSymMatrix.hpp"
#include "IpMultiVectorMatrix.hpp"
#include "IpDenseVector.hpp"
#include "IpDenseGenMatrix.hpp"
#include "IpDenseSymMatrix.hpp"

namespace Ipopt
{

  /** Implementation of the HessianUpdater for limit-memory
   *  quasi-Newton approximation of the Lagrangian Hessian.
   */
  class LimMemQuasiNewtonUpdater : public HessianUpdater
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    LimMemQuasiNewtonUpdater(bool update_for_resto);

    /** Default destructor */
    virtual ~LimMemQuasiNewtonUpdater()
    {}
    //@}

    /** overloaded from AlgorithmStrategyObject */
    virtual bool InitializeImpl(const OptionsList& options,
                                const std::string& prefix);

    /** Update the Hessian based on the current information in IpData.
     */
    virtual void UpdateHessian();

    /** Methods for OptionsList */
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
    /** Copy Constructor */
    LimMemQuasiNewtonUpdater(const LimMemQuasiNewtonUpdater&);

    /** Overloaded Equals Operator */
    void operator=(const LimMemQuasiNewtonUpdater&);
    //@}

    /** Matrix space for the low-rank Hessian approximation. */
    SmartPtr<const LowRankUpdateSymMatrixSpace> h_space_;

    /** @name Algorithmic parameters */
    //@{
    /** Size of memory for limited memory update. */
    Index limited_memory_max_history_;
    /** enumeration for the Hessian update type. */
    enum LMUpdateType {
      BFGS=0,
      SR1
    };
    /** Type of Hessian update. */
    LMUpdateType limited_memory_update_type_;
    /** enumeration for the Hessian initialization. */
    enum LMInitialization {
      SCALAR1=0,
      SCALAR2,
      CONSTANT
    };
    /** How to choose B0 in the low-rank update. */
    LMInitialization limited_memory_initialization_;
    /** Value of B0 (as this multiple of the identity in certain
     *  situations.)  */
    Number limited_memory_init_val_;
    /** Number of successive iterations of skipped updates after which
     *  the approximation is reset. */
    Index limited_memory_max_skipping_;
    /** Minimal safeguard value for sigma */
    Number sigma_safe_min_;
    /** Maximal safeguard value for sigma */
    Number sigma_safe_max_;
    //@}

    /** Flag indicating if the update is to be done for the original
     *  NLP or for the restoration phase NLP.  In the latter case, we
     *  are performing a "structured" update, taking into account the
     *  first explicit term in the objective function of the form
     *  eta*D_r*x_k */
    const bool update_for_resto_;
    /** Most recent value for eta in the restoration phase objective
     *  function (only for update_for_resto_ = true) */
    Number last_eta_;
    /** Current DR_x scaling factors in the restoration phase
     *  objective function (only for update_for_resto_ = true).  This
     *  should not change throughout one restoration phase. */
    SmartPtr<const Vector> curr_DR_x_;
    /** Tag  for curr_DR_x_ */
    TaggedObject::Tag curr_DR_x_tag_;
    /** Current DR_x scaling factors in the restoration phase
     *  objective function in the smaller space for the approximation -
     *  this is only computed if the space is indeed smaller than the
     *  x space (only for update_for_resto_ = true) */
    SmartPtr<const Vector> curr_red_DR_x_;
    /** Current value of weighing factor eta in the restoration phase
     *  objective function (only for update_for_resto_ = true) */
    Number curr_eta_;
    /** Flag inidicating whether DR_x or eta have changed since the
     *  last update */
    bool eta_changed_;

    /** Counter for successive iterations in which the update was
     *  skipped */
    Index lm_skipped_iter_;

    /** @name Information for the limited memory update */
    //@{
    /** current size of limited memory */
    Index curr_lm_memory_;
    /** s pairs for the recent iterations */
    SmartPtr<MultiVectorMatrix> S_;
    /** y pairs for the recent iterations.  If update_for_resto is
     *  true, then this includes only the information for the
     *  constraints. */
    SmartPtr<MultiVectorMatrix> Y_;
    /** For restoration phase update: Y without the quadratic
     *  objective function part */
    SmartPtr<MultiVectorMatrix> Ypart_;
    /** Diagonal elements D_k for compact formulation from last
     *  update. */
    SmartPtr<DenseVector> D_;
    /** Matrix L_k for compact formulation from last update.*/
    SmartPtr<DenseGenMatrix> L_;
    /** First term (starting matrix) for the approximation. */
    SmartPtr<Vector> B0_;
    /** First term (starting matrix) for the approximation.  If that
     *  first terms is a multiple of the identy, sigma give that
     *  factor.  Otherwise sigma = -1.  */
    Number sigma_;
    /** V in LowRankUpdateMatrix from last update */
    SmartPtr<MultiVectorMatrix> V_;
    /** U in LowRankUpdateMatrix from last update */
    SmartPtr<MultiVectorMatrix> U_;
    /** For efficient implementation, we store the pairwise products
     *  for s's. */
    SmartPtr<DenseSymMatrix> SdotS_;
    /** Flag indicating whether SdotS_ is update to date from most
     * recent update. */
    bool SdotS_uptodate_;
    /** DR * S (only for restoration phase) */
    SmartPtr<MultiVectorMatrix> DRS_;
    /** For efficient implementation, we store the S^T S DR * S. Only
     *  for restoration phase. */
    SmartPtr<DenseSymMatrix> STDRS_;
    /** Primal variables x from most recent update */
    SmartPtr<const Vector> last_x_;
    /** Gradient of objective function w.r.t. x at x_last_ */
    SmartPtr<const Vector> last_grad_f_;
    /** Jacobian for equality constraints w.r.t x at x_last */
    SmartPtr<const Matrix> last_jac_c_;
    /** Jacobian for inequality constraints w.r.t x at x_last */
    SmartPtr<const Matrix> last_jac_d_;
    /** current size of limited memory */
    Index curr_lm_memory_old_;
    /** s pairs for the recent iterations (backup) */
    SmartPtr<MultiVectorMatrix> S_old_;
    /** y pairs for the recent iterations.  If update_for_resto is
     *  true, then this includes only the information for the
     *  constraints. (backup) */
    SmartPtr<MultiVectorMatrix> Y_old_;
    /** For restoration phase update: Y without the quadratic
     *  objective function part (backup) */
    SmartPtr<MultiVectorMatrix> Ypart_old_;
    /** Diagonal elements D_k for compact formulation from last
     *  update (backup). */
    SmartPtr<DenseVector> D_old_;
    /** Matrix L_k for compact formulation from last update (backup).*/
    SmartPtr<DenseGenMatrix> L_old_;
    /** First term (starting matrix) for the approximation (backup). */
    SmartPtr<Vector> B0_old_;
    /** First term (starting matrix) for the approximation.  If that
     *  first terms is a multiple of the identy, sigma give that
     *  factor.  Otherwise sigma = -1.  (backup) */
    Number sigma_old_;
    /** V in LowRankUpdateMatrix from last update (backup) */
    SmartPtr<MultiVectorMatrix> V_old_;
    /** U in LowRankUpdateMatrix from last update (backup) */
    SmartPtr<MultiVectorMatrix> U_old_;
    /** For efficient implementation, we store the pairwise products
     *  for s's (backup). */
    SmartPtr<DenseSymMatrix> SdotS_old_;
    /** Flag indicating whether SdotS_ is update to date from most
     *  recent update (backup). */
    bool SdotS_uptodate_old_;
    /** DR * S (only for restoration phase) (backup) */
    SmartPtr<MultiVectorMatrix> DRS_old_;
    /** For efficient implementation, we store the S^T S DR * S. Only
     *  for restoration phase. (backup) */
    SmartPtr<DenseSymMatrix> STDRS_old_;
    //@}

    /** @name Auxilliary function */
    //@{
    /** Method deciding whether the BFGS update should be skipped.  It
     *  returns true, if no update is to be performed this time. If
     *  Powell-damping is performed, the Vectors s_new and y_new,
     *  might be adapted. */
    bool CheckSkippingBFGS(Vector& s_new, Vector& y_new);
    /** Update the internal data, such as the S, Y, L, D etc matrices
     *  and vectors that are required for computing the compact
     *  representation.  The method returns true if the limited memory
     *  history grew (i.e., curr_lm_memory_ was increased). */
    bool UpdateInternalData(const Vector& s_new, const Vector& y_new,
                            SmartPtr<Vector> ypart_new);
    /** Given a MutliVector V, create a new MultiVectorSpace with one
     *  more column, and return V as a member of that space,
     *  consisting of all previous vectors, and in addition v_new in
     *  the last column.  If V is NULL, then a new MatrixSpace with
     *  one column is created. */
    void AugmentMultiVector(SmartPtr<MultiVectorMatrix>& V,
                            const Vector& v_new);
    /** Given a DenseVector V, create a new DenseVectorSpace with one
     *  more row, and return V as a member of that space,
     *  consisting of all previous elements, and in addition v_new in
     *  the last row.  If V is NULL, then a new DenseVectorSpace with
     *  dimension one is created. */
    void AugmentDenseVector(SmartPtr<DenseVector>& V,
                            Number v_new);
    /** Given a strictly-lower triangular square DenseGenMatrix V,
     *  create a new DenseGenMatrixSpace with one more dimension, and
     *  return V as a member of that space, consisting of all previous
     *  elements, and in addition elements s_i^Ty_j for (i<j), where s
     *  and y are the vectors in the MultiVectors S and Y.  If V is
     *  NULL, then a new DenseGenMatrixSpace with dimension one is
     *  created. */
    void AugmentLMatrix(SmartPtr<DenseGenMatrix>& V,
                        const MultiVectorMatrix& S,
                        const MultiVectorMatrix& Y);
    /** Given a DenseSymMatrix V, create a new DenseGenMatrixSpace
     *  with one more dimension, and return V as a member of that
     *  space, consisting of all previous elements, and in addition
     *  elements s_i^Ts_j for the new entries, where s are the vectors
     *  in the MultiVector S.  If V is NULL, then a new
     *  DenseGenMatrixSpace with dimension one is created. */
    void AugmentSdotSMatrix(SmartPtr<DenseSymMatrix>& V,
                            const MultiVectorMatrix& S);
    /** Given a DenseSymMatrix V, create a new DenseGenMatrixSpace
     *  with one more dimension, and return V as a member of that
     *  space, consisting of all previous elements, and in addition
     *  elements s_i^TDRs_j for the new entries, where s are the
     *  vectors in the MultiVector S, and DRs are the vectors in DRS.
     *  If V is NULL, then a new DenseGenMatrixSpace with dimension
     *  one is created. */
    void AugmentSTDRSMatrix(SmartPtr<DenseSymMatrix>& V,
                            const MultiVectorMatrix& S,
                            const MultiVectorMatrix& DRS);

    /** Given a MutliVector V, get rid of the first column, shift all
     *  other columns to the left, and make v_new the last column.
     *  The entity that V points to at the call, is not changed - a
     *  new entity is created in the method and returned as V. */
    void ShiftMultiVector(SmartPtr<MultiVectorMatrix>& V, const Vector& v_new);
    /** Given a DenseVector V, get rid of the first element, shift all
     *  other elements one position to the top, and make v_new the
     *  last entry. The entity that V points to at the call, is not
     *  changed - a new entity is created in the method and returned
     *  as V. */
    void ShiftDenseVector(SmartPtr<DenseVector>& V, Number v_new);
    /** Given a strictly-lower triangular square DenseGenMatrix V,
     *  shift everything one row and column up, and fill the new
     *  strictly lower triangular entries as s_i^Ty_j for (i<j), where
     *  s and y are the vectors in the MultiVectors S and Y. The
     *  entity that V points to at the call, is not changed - a new
     *  entity is created in the method and returned as V. */
    void ShiftLMatrix(SmartPtr<DenseGenMatrix>& V,
                      const MultiVectorMatrix& S,
                      const MultiVectorMatrix& Y);
    /** Given a DenseSymMatrix V, shift everything up one row and
     *  column, and fill the new entries as s_i^Ts_j, where s are the
     *  vectors in the MultiVector S. The entity that V points to at
     *  the call, is not changed - a new entity is created in the
     *  method and returned as V. */
    void ShiftSdotSMatrix(SmartPtr<DenseSymMatrix>& V,
                          const MultiVectorMatrix& S);
    /** Given a DenseSymMatrix V, shift everything up one row and
     *  column, and fill the new entries as s_i^TDRs_j, where s are
     *  the vectors in the MultiVector S, and DRs are the vectors in
     *  DRS. The entity that V points to at the call, is not changed -
     *  a new entity is created in the method and returned as V. */
    void ShiftSTDRSMatrix(SmartPtr<DenseSymMatrix>& V,
                          const MultiVectorMatrix& S,
                          const MultiVectorMatrix& DRS);
    /** Method for recomputing Y from scratch, using Ypart (only for
     *  restoration phase) */
    void RecalcY(Number eta, const Vector& DR_x,
                 MultiVectorMatrix& S,
                 MultiVectorMatrix& Ypart,
                 SmartPtr<MultiVectorMatrix>& Y);
    /** Method for recomputing D from S and Y */
    void RecalcD(MultiVectorMatrix& S,
                 MultiVectorMatrix& Y,
                 SmartPtr<DenseVector>& D);
    /** Method for recomputing L from S and Y */
    void RecalcL(MultiVectorMatrix& S,
                 MultiVectorMatrix& Y,
                 SmartPtr<DenseGenMatrix>& L);
    /** Split the eigenvectors into negative and positive ones.  Given
     *  the eigenvectors in Q and the eigenvalues (in ascending order)
     *  in, this returns Qminus as the negative eigenvectors times
     *  sqrt(-eval), and Qplus as the positive eigenvectors times
     *  sqrt(eval). If Qminus or Qplus is NULL, it means that there
     *  are not negetive or positive eigenvalues.  Q might be changed
     *  during this call.  The return value is true, if the ratio of
     *  the smallest over the largest eigenvalue (in absolute values)
     *  is too small; in that case, the update should be skipped.  */
    bool SplitEigenvalues(DenseGenMatrix& Q, const DenseVector& E,
                          SmartPtr<DenseGenMatrix>& Qminus,
                          SmartPtr<DenseGenMatrix>& Qplus);
    /** Store a copy of the pointers to the internal data (S, Y, D, L,
     *  SdotS, curr_lm_memory) This is called in case the update is
     *  started but skipped during the process. */
    void StoreInternalDataBackup();
    /** Restore the copy of the pointers to the internal data most
     *  recently stored with StoreInternalDataBackup(). */
    void RestoreInternalDataBackup();
    /** Release anything that we allocated for
     *  StoreInternalDataBackup and is no longer needed. */
    void ReleaseInternalDataBackup();
    /** Set the W field in IpData based on the current values of
     *  B0_, V_, and U_ */
    void SetW();
    //@}

  };

} // namespace Ipopt

#endif
