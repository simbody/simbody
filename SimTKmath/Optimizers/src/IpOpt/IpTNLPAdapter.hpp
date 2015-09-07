// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpTNLPAdapter.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPTNLPADAPTER_HPP__
#define __IPTNLPADAPTER_HPP__

#include "IpNLP.hpp"
#include "IpTNLP.hpp"
#include "IpOrigIpoptNLP.hpp"

namespace Ipopt
{

  // forward declarations
  class ExpansionMatrix;
  class ExpansionMatrixSpace;
  class IteratesVector;

  /** This class Adapts the TNLP interface so it looks like an NLP interface.
   *  This is an Adapter class (Design Patterns) that converts  a TNLP to an
   *  NLP. This allows users to write to the "more convenient" TNLP interface.
   */
  class TNLPAdapter : public NLP
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default constructor */
    TNLPAdapter(const SmartPtr<TNLP> tnlp,
                const SmartPtr<const Journalist> jnlst = NULL);

    /** Default destructor */
    virtual ~TNLPAdapter();
    //@}

    /**@name Exceptions */
    //@{
    DECLARE_STD_EXCEPTION(INVALID_TNLP);
    DECLARE_STD_EXCEPTION(ERROR_IN_TNLP_DERIVATIVE_TEST);
    //@}

    /** @name TNLPAdapter Initialization. */
    //@{
    virtual bool ProcessOptions(const OptionsList& options,
                                const std::string& prefix) override;

    /** Method for creating the derived vector / matrix types
     *  (Do not delete these, the ). */
    virtual bool GetSpaces(SmartPtr<const VectorSpace>& x_space,
                           SmartPtr<const VectorSpace>& c_space,
                           SmartPtr<const VectorSpace>& d_space,
                           SmartPtr<const VectorSpace>& x_l_space,
                           SmartPtr<const MatrixSpace>& px_l_space,
                           SmartPtr<const VectorSpace>& x_u_space,
                           SmartPtr<const MatrixSpace>& px_u_space,
                           SmartPtr<const VectorSpace>& d_l_space,
                           SmartPtr<const MatrixSpace>& pd_l_space,
                           SmartPtr<const VectorSpace>& d_u_space,
                           SmartPtr<const MatrixSpace>& pd_u_space,
                           SmartPtr<const MatrixSpace>& Jac_c_space,
                           SmartPtr<const MatrixSpace>& Jac_d_space,
                           SmartPtr<const SymMatrixSpace>& Hess_lagrangian_space) override;

    /** Method for obtaining the bounds information */
    virtual bool GetBoundsInformation(const Matrix& Px_L,
                                      Vector& x_L,
                                      const Matrix& Px_U,
                                      Vector& x_U,
                                      const Matrix& Pd_L,
                                      Vector& d_L,
                                      const Matrix& Pd_U,
                                      Vector& d_U) override;

    /** Method for obtaining the starting point
     *  for all the iterates. */
    virtual bool GetStartingPoint(
      SmartPtr<Vector> x,
      bool need_x,
      SmartPtr<Vector> y_c,
      bool need_y_c,
      SmartPtr<Vector> y_d,
      bool need_y_d,
      SmartPtr<Vector> z_L,
      bool need_z_L,
      SmartPtr<Vector> z_U,
      bool need_z_U
    ) override;

    /** Method for obtaining an entire iterate as a warmstart point.
     *  The incoming IteratesVector has to be filled. */
    virtual bool GetWarmStartIterate(IteratesVector& warm_start_iterate) override;
    //@}

    /** @name TNLPAdapter evaluation routines. */
    //@{
    virtual bool Eval_f(const Vector& x, Number& f) override;

    virtual bool Eval_grad_f(const Vector& x, Vector& g_f) override;

    virtual bool Eval_c(const Vector& x, Vector& c) override;

    virtual bool Eval_jac_c(const Vector& x, Matrix& jac_c) override;

    virtual bool Eval_d(const Vector& x, Vector& d) override;

    virtual bool Eval_jac_d(const Vector& x, Matrix& jac_d) override;

    virtual bool Eval_h(const Vector& x,
                        Number obj_factor,
                        const Vector& yc,
                        const Vector& yd,
                        SymMatrix& h) override;

    virtual void GetScalingParameters(
      const SmartPtr<const VectorSpace> x_space,
      const SmartPtr<const VectorSpace> c_space,
      const SmartPtr<const VectorSpace> d_space,
      Number& obj_scaling,
      SmartPtr<Vector>& x_scaling,
      SmartPtr<Vector>& c_scaling,
      SmartPtr<Vector>& d_scaling) const override;
    //@}

    /** @name Solution Reporting Methods */
    //@{
    virtual void FinalizeSolution(SolverReturn status,
                                  const Vector& x, const Vector& z_L, const Vector& z_U,
                                  const Vector& c, const Vector& d,
                                  const Vector& y_c, const Vector& y_d,
                                  Number obj_value) override;
    virtual bool IntermediateCallBack(AlgorithmMode mode,
                                      Index iter, Number obj_value,
                                      Number inf_pr, Number inf_du,
                                      Number mu, Number d_norm,
                                      Number regularization_size,
                                      Number alpha_du, Number alpha_pr,
                                      Index ls_trials,
                                      const IpoptData* ip_data,
                                      IpoptCalculatedQuantities* ip_cq) override;
    //@}

    /** Method returning information on quasi-Newton approximation. */
    virtual void
    GetQuasiNewtonApproximationSpaces(SmartPtr<VectorSpace>& approx_space,
                                      SmartPtr<Matrix>& P_approx) override;

    /** Enum for treatment of fixed variables option */
    enum FixedVariableTreatmentEnum {
      MAKE_PARAMETER=0,
      MAKE_CONSTRAINT
    };

    /** Enum for specifying which derivative test is to be performed. */
    enum DerivativeTestEnum {
      NO_TEST=0,
      FIRST_ORDER_TEST,
      SECOND_ORDER_TEST
    };

    /** Method for performing the derivative test */
    bool CheckDerivatives(DerivativeTestEnum deriv_test);

    /** @name Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

    /** Accessor method for the underlying TNLP. */
    SmartPtr<TNLP> tnlp() const
    {
      return tnlp_;
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
    /** Copy Constructor */
    TNLPAdapter(const TNLPAdapter&);

    /** Overloaded Equals Operator */
    void operator=(const TNLPAdapter&);
    //@}

    /** @name Methods for translating data for IpoptNLP into the TNLP
     *  data.  These methods are used to obtain the current (or
     *  final) data for the TNLP formulation from the IpoptNLP
     *  structure. */
    //@{
    /** Sort the primal variables, and add the fixed values in x */
    void ResortX(const Vector& x, Number* x_orig);
    void ResortG(const Vector& c, const Vector& d, Number *g_orig);
    void ResortBnds(const Vector& x_L, Number* x_L_orig,
                    const Vector& x_U, Number* x_U_orig);
    //@}

    /** Pointer to the TNLP class (class specific to Number* vectors and
     *  harwell triplet matrices) */
    SmartPtr<TNLP> tnlp_;

    /** Journalist */
    SmartPtr<const Journalist> jnlst_;

    /**@name Algorithmic parameters */
    //@{
    /** Value for a lower bound that denotes -infinity */
    Number nlp_lower_bound_inf_;
    /** Value for a upper bound that denotes infinity */
    Number nlp_upper_bound_inf_;
    /** Flag indicating how fixed variables should be handled */
    FixedVariableTreatmentEnum fixed_variable_treatment_;
    /* Maximal slack for one-sidedly bounded variables.  If a
     *  variable has only one bound, say a lower bound xL, then an
     *  upper bound xL + max_onesided_bound_slack_.  If this value is
     *  zero, no upper bound is added. */
    /* Took this out:  Number max_onesided_bound_slack_; */
    /** Enum indicating whether and which derivative test should be
     *  performed at starting point. */
    DerivativeTestEnum derivative_test_;
    /** Size of the perturbation for the derivative test */
    Number derivative_test_perturbation_;
    /** Relative threshold for marking deviation from finite
     *  difference test */
    Number derivative_test_tol_;
    /** Flag indicating if all test values should be printed, or only
     *  those violating the threshold. */
    bool derivative_test_print_all_;
    /** Flag indicating whether the TNLP with identical structure has
     *  already been solved before. */
    bool warm_start_same_structure_;
    /** Flag indicating what Hessian information is to be used. */
    HessianApproximationType hessian_approximation_;
    //@}

    /**@name Problem Size Data */
    //@{
    /** full dimension of x (fixed + non-fixed) */
    Index n_full_x_;
    /** full dimension of g (c + d) */
    Index n_full_g_;
    /** non-zeros of the jacobian of c */
    Index nz_jac_c_;
    /** non-zeros of the jacobian of c without added constraints for
     *  fixed variables. */
    Index nz_jac_c_no_fixed_;
    /** non-zeros of the jacobian of d */
    Index nz_jac_d_;
    /** number of non-zeros in full-size Jacobian of g */
    Index nz_full_jac_g_;
    /** number of non-zeros in full-size Hessian */
    Index nz_full_h_;
    /** number of non-zeros in the non-fixed-size Hessian */
    Index nz_h_;
    /** Number of fixed variables */
    Index n_x_fixed_;
    //@}

    /** Numbering style of variables and constraints */
    TNLP::IndexStyleEnum index_style_;

    /** @name Local copy of spaces (for warm start) */
    //@{
    SmartPtr<const VectorSpace> x_space_;
    SmartPtr<const VectorSpace> c_space_;
    SmartPtr<const VectorSpace> d_space_;
    SmartPtr<const VectorSpace> x_l_space_;
    SmartPtr<const MatrixSpace> px_l_space_;
    SmartPtr<const VectorSpace> x_u_space_;
    SmartPtr<const MatrixSpace> px_u_space_;
    SmartPtr<const VectorSpace> d_l_space_;
    SmartPtr<const MatrixSpace> pd_l_space_;
    SmartPtr<const VectorSpace> d_u_space_;
    SmartPtr<const MatrixSpace> pd_u_space_;
    SmartPtr<const MatrixSpace> Jac_c_space_;
    SmartPtr<const MatrixSpace> Jac_d_space_;
    SmartPtr<const SymMatrixSpace> Hess_lagrangian_space_;
    //@}

    /**@name Local Copy of the Data */
    //@{
    Number* full_x_; /** copy of the full x vector (fixed & non-fixed) */
    Number* full_lambda_; /** copy of lambda (yc & yd) */
    Number* full_g_; /** copy of g (c & d) */
    Number* jac_g_; /** the values for the full jacobian of g */
    Number* c_rhs_; /** the rhs values of c */
    //@}

    /**@name Tags for deciding when to update internal copies of vectors */
    //@{
    TaggedObject::Tag x_tag_for_iterates_;
    TaggedObject::Tag y_c_tag_for_iterates_;
    TaggedObject::Tag y_d_tag_for_iterates_;
    TaggedObject::Tag x_tag_for_g_;
    TaggedObject::Tag x_tag_for_jac_g_;
    //@}

    /**@name Methods to update the values in the local copies of vectors */
    //@{
    bool update_local_x(const Vector& x);
    bool update_local_lambda(const Vector& y_c, const Vector& y_d);
    //@}

    /**@name Internal routines for evaluating g and jac_g (values stored since
     * they are used in both c and d routins */
    //@{
    bool internal_eval_g(bool new_x);
    bool internal_eval_jac_g(bool new_x);
    //@}

    /**@name Internal Permutation Spaces and matrices
     */
    //@{
    /** Expansion from fixed x (ipopt) to full x */
    SmartPtr<ExpansionMatrix> P_x_full_x_;
    SmartPtr<ExpansionMatrixSpace> P_x_full_x_space_;

    /** Expansion from fixed x_L (ipopt) to full x */
    SmartPtr<ExpansionMatrix> P_x_x_L_;
    SmartPtr<ExpansionMatrixSpace> P_x_x_L_space_;

    /** Expansion from fixed x_U (ipopt) to full x */
    SmartPtr<ExpansionMatrix> P_x_x_U_;
    SmartPtr<ExpansionMatrixSpace> P_x_x_U_space_;

    /** Expansion from c only (ipopt) to full ampl c */
    SmartPtr<ExpansionMatrixSpace> P_c_g_space_;
    SmartPtr<ExpansionMatrix> P_c_g_;

    /** Expansion from d only (ipopt) to full ampl d */
    SmartPtr<ExpansionMatrixSpace> P_d_g_space_;
    SmartPtr<ExpansionMatrix> P_d_g_;

    Index* jac_idx_map_;
    Index* h_idx_map_;

    /** Position of fixed variables. This is required for a warm start */
    Index* x_fixed_map_;
    //@}
  };

} // namespace Ipopt

#endif
