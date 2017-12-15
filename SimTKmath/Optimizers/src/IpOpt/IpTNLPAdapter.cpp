// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpTNLPAdapter.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpTNLPAdapter.hpp"
#include "IpBlas.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpDenseVector.hpp"
#include "IpExpansionMatrix.hpp"
#include "IpGenTMatrix.hpp"
#include "IpSymTMatrix.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

#include <cstdio>

// Keeps MS VC++ 8 quiet about sprintf, strcpy, etc.
#ifdef _MSC_VER
#pragma warning(disable:4996)
#endif



namespace SimTKIpopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  TNLPAdapter::TNLPAdapter(const SmartPtr<TNLP> tnlp,
                           const SmartPtr<const Journalist> jnlst /* = NULL */)
      :
      tnlp_(tnlp),
      jnlst_(jnlst),
      n_full_x_(-1),
      n_full_g_(-1),
      nz_jac_c_(-1),
      nz_jac_c_no_fixed_(-1),
      nz_jac_d_(-1),
      nz_full_jac_g_(-1),
      nz_full_h_(-1),
      nz_h_(-1),
      full_x_(NULL),
      full_lambda_(NULL),
      full_g_(NULL),
      jac_g_(NULL),
      c_rhs_(NULL),
      x_tag_for_iterates_(0),
      y_c_tag_for_iterates_(0),
      y_d_tag_for_iterates_(0),
      x_tag_for_g_(0),
      x_tag_for_jac_g_(0),
      jac_idx_map_(NULL),
      h_idx_map_(NULL),
      x_fixed_map_(NULL)
  {
    ASSERT_EXCEPTION(IsValid(tnlp_), INVALID_TNLP,
                     "The TNLP passed to TNLPAdapter is NULL. This MUST be a valid TNLP!");
  }

  TNLPAdapter::~TNLPAdapter()
  {
    delete [] full_x_;
    delete [] full_lambda_;
    delete [] full_g_;
    delete [] jac_g_;
    delete [] c_rhs_;
    delete [] jac_idx_map_;
    delete [] h_idx_map_;
    delete [] x_fixed_map_;
  }

  void TNLPAdapter::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("NLP");
    roptions->AddNumberOption(
      "nlp_lower_bound_inf",
      "any bound less or equal this value will be considered -inf (i.e. not lower bounded).",
      -1e19);
    roptions->AddNumberOption(
      "nlp_upper_bound_inf",
      "any bound greater or this value will be considered +inf (i.e. not upper bounded).",
      1e19);
    roptions->AddStringOption2(
      "fixed_variable_treatment",
      "Determines how fixed variables should be handled.",
      "make_parameter",
      "make_parameter", "Remove fixed variable from optimization variables",
      "make_constraint", "Add equality constraints fixing variables",
      "The main difference between those two options is that the starting "
      "point in the \"make_constraint\" case still has the fixed variables at "
      "their given values, whereas in the other case the functions are always "
      "evaluated with the fixed values for those variables.  Also, for "
      "\"make_constraints\", bound multipliers are computed for the fixed "
      "variables.");
    roptions->SetRegisteringCategory("Derivative Checker");
    roptions->AddStringOption3(
      "derivative_test",
      "Enable derivative checker",
      "none",
      "none", "do not perform derivative test",
      "first-order", "perform test of first derivatives at starting point",
      "second-order", "perform test of first and second derivatives at starting point",
      "If this option is enabled, a (slow) derivative test will be performed "
      "before the optimization.  The test is performed at the user provided "
      "starting point and marks derivative values that seem suspicious");
    roptions->AddLowerBoundedNumberOption(
      "derivative_test_perturbation",
      "Size of the finite difference perturbation in derivative test.",
      0., true,
      1e-8,
      "This determines the relative perturbation of the variable entries.");
    roptions->AddLowerBoundedNumberOption(
      "derivative_test_tol",
      "Threshold for indicating wrong derivative.",
      0., true,
      1e-4,
      "If the relative deviation of the estimated derivative from the given "
      "one is larger than this value, the corresponding derivative is marked "
      "as wrong.");
    roptions->AddStringOption2(
      "derivative_test_print_all",
      "Indicates whether information for all estimated derivatives should be printed.",
      "no",
      "no", "Print only suspect derivatives",
      "yes", "Print all derivatives",
      "Determines verbosity of derivative checker.");
  }

  bool TNLPAdapter::ProcessOptions(const OptionsList& options,
                                   const std::string& prefix)
  {
    DBG_START_METH("TNLPAdapter::ProcessOptions", dbg_verbosity);
    options.GetNumericValue("nlp_lower_bound_inf", nlp_lower_bound_inf_, prefix);
    options.GetNumericValue("nlp_upper_bound_inf", nlp_upper_bound_inf_, prefix);

    ASSERT_EXCEPTION(nlp_lower_bound_inf_ < nlp_upper_bound_inf_,
                     OPTION_INVALID,
                     "Option \"nlp_lower_bound_inf\" must be smaller than \"nlp_upper_bound_inf\".");

    Index enum_int;
    options.GetEnumValue("fixed_variable_treatment", enum_int, prefix);
    fixed_variable_treatment_ = FixedVariableTreatmentEnum(enum_int);
    options.GetEnumValue("derivative_test", enum_int, prefix);
    derivative_test_ = DerivativeTestEnum(enum_int);
    options.GetNumericValue("derivative_test_perturbation",
                            derivative_test_perturbation_, prefix);
    options.GetNumericValue("derivative_test_tol",
                            derivative_test_tol_, prefix);
    options.GetBoolValue("derivative_test_print_all",
                         derivative_test_print_all_, prefix);

    // The option warm_start_same_structure is registered by OrigIpoptNLP
    options.GetBoolValue("warm_start_same_structure",
                         warm_start_same_structure_, prefix);
    // The following is registered in OrigIpoptNLP
    options.GetEnumValue("hessian_approximation", enum_int, prefix);
    hessian_approximation_ = HessianApproximationType(enum_int);

    return true;
  }

  bool TNLPAdapter::GetSpaces(SmartPtr<const VectorSpace>& x_space,
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
                              SmartPtr<const SymMatrixSpace>& Hess_lagrangian_space)
  {
    DBG_START_METH("TNLPAdapter::GetSpaces", dbg_verbosity);

    // First, if required, perform derivative test
    if (derivative_test_ != NO_TEST) {
      bool retval = CheckDerivatives(derivative_test_);
      if (!retval) {
        return retval;
      }
    }

    if (warm_start_same_structure_) {
      ASSERT_EXCEPTION(full_x_, INVALID_WARMSTART,
                       "warm_start_same_structure chosen, but TNLPAdapter is called for the first time.");
      if (IsValid(jnlst_)) {
        jnlst_->Printf(J_DETAILED, J_INITIALIZATION,
                       "Reusing previous information for warm start in TNLPAdapter.\n");
      }
    }
    else {
      // In case the Adapter has been used before, but this is not a
      // warm start, make sure we delete all previously allocated
      // memory
      delete [] full_x_;
      delete [] full_lambda_;
      delete [] full_g_;
      delete [] jac_g_;
      delete [] c_rhs_;
      delete [] jac_idx_map_;
      delete [] h_idx_map_;
      delete [] x_fixed_map_;
    }

    // Get the full dimensions of the problem
    Index n_full_x, n_full_g, nz_full_jac_g, nz_full_h;
    tnlp_->get_nlp_info(n_full_x, n_full_g, nz_full_jac_g,
                        nz_full_h, index_style_);
    ASSERT_EXCEPTION(!warm_start_same_structure_ ||
                     (n_full_x == n_full_x_ &&
                      n_full_g == n_full_g_ &&
                      nz_full_jac_g == nz_full_jac_g_ &&
                      nz_full_h == nz_full_h_),
                     INVALID_WARMSTART,
                     "warm_start_same_structure chosen, but problem dimensions are different.");
    n_full_x_ = n_full_x;
    n_full_g_ = n_full_g;
    nz_full_jac_g_ = nz_full_jac_g;
    nz_full_h_ = nz_full_h;

    if (!warm_start_same_structure_) {
      // create space to store vectors that are the full length of x
      full_x_ = new Number[n_full_x_];

      // create space to store vectors that area the full length of lambda
      full_lambda_ = new Number[n_full_g_];

      // create space to store vectors that are the full length of g
      full_g_ = new Number[n_full_g_];

      // allocate internal space to store the full jacobian
      jac_g_ = new Number[nz_full_jac_g_];

      //*********************************************************
      // Create the spaces and permutation spaces
      //*********************************************************
      /* Spaces for x, x_L, and x_U. We need to remove the fixed variables
       * and find out which bounds do not exist. */
      Number* x_l = new Number[n_full_x_];
      Number* x_u = new Number[n_full_x_];
      Number* g_l = new Number[n_full_g_];
      Number* g_u = new Number[n_full_g_];
      tnlp_->get_bounds_info(n_full_x_, x_l, x_u, n_full_g_, g_l, g_u);

      Index n_x_var = 0;
      Index n_x_l = 0;
      Index n_x_u = 0;
      Index* x_not_fixed_map = new Index[n_full_x_];
      Index* x_fixed_map_tmp = new Index[n_full_x_];
      Index* x_l_map = new Index[n_full_x_];
      Index* x_u_map = new Index[n_full_x_];
      n_x_fixed_ = 0;

      for (Index i=0; i<n_full_x_; i++) {
        Number lower_bound = x_l[i];
        Number upper_bound = x_u[i];
        if (lower_bound == upper_bound) {
          switch (fixed_variable_treatment_) {
            case MAKE_PARAMETER:
            // Variable is fixed, remove it from the problem
            full_x_[i] = lower_bound;
            x_fixed_map_tmp[n_x_fixed_] = i;
            n_x_fixed_++;
            break;
            case MAKE_CONSTRAINT:
            x_fixed_map_tmp[n_x_fixed_] = i; // don't really need this
            // array then
            n_x_fixed_++;
            x_not_fixed_map[n_x_var] = i;
            n_x_var++;
            break;
            default:
            DBG_ASSERT(false && "invalid fixed_variable_treatment_");
          }
        }
        else if (lower_bound > upper_bound) {
          char string[128];
          sprintf(string, "There are inconsistent bounds on variable %d: lower = %25.16e and upper = %25.16e.", i, lower_bound, upper_bound);
          delete [] x_not_fixed_map;
          delete [] x_fixed_map_tmp;
          delete [] x_l_map;
          delete [] x_u_map;
          THROW_EXCEPTION(INVALID_TNLP, string);
        }
        else {
          x_not_fixed_map[n_x_var] = i;
          if (lower_bound > nlp_lower_bound_inf_) {
            x_l_map[n_x_l] = n_x_var;
            n_x_l++;
          }

          if (upper_bound < nlp_upper_bound_inf_) {
            x_u_map[n_x_u] = n_x_var;
            n_x_u++;
          }
          n_x_var++;
        }
      }

      // If there are fixed variables, we keep their position around
      // for a possible warm start later or if fixed variables are
      // treated by added equality constraints
      if (n_x_fixed_>0) {
        x_fixed_map_ = new Index[n_x_fixed_];
        for (Index i=0; i<n_x_fixed_; i++) {
          x_fixed_map_[i] = x_fixed_map_tmp[i];
        }
      }
      else {
        x_fixed_map_ = NULL;
      }
      delete [] x_fixed_map_tmp;

      x_space_ = new DenseVectorSpace(n_x_var);
      x_l_space_ = new DenseVectorSpace(n_x_l);
      x_u_space_ = new DenseVectorSpace(n_x_u);

      if (n_x_fixed_>0 && fixed_variable_treatment_==MAKE_PARAMETER) {
        P_x_full_x_space_ = new ExpansionMatrixSpace(n_full_x_, n_x_var, x_not_fixed_map);
        P_x_full_x_ = P_x_full_x_space_->MakeNewExpansionMatrix();
      }
      else {
        P_x_full_x_space_ = NULL;
        P_x_full_x_ = NULL;
      }

      P_x_x_L_space_ = new ExpansionMatrixSpace(n_x_var, n_x_l, x_l_map);
      px_l_space_ = GetRawPtr(P_x_x_L_space_);
      P_x_x_L_ = P_x_x_L_space_->MakeNewExpansionMatrix();
      P_x_x_U_space_ = new ExpansionMatrixSpace(n_x_var, n_x_u, x_u_map);
      px_u_space_ = GetRawPtr(P_x_x_U_space_);
      P_x_x_U_ = P_x_x_U_space_->MakeNewExpansionMatrix();

      delete [] x_l;
      x_l = NULL;
      delete [] x_u;
      x_u = NULL;
      delete [] x_not_fixed_map;
      x_not_fixed_map = NULL;
      delete [] x_l_map;
      x_l_map = NULL;
      delete [] x_u_map;
      x_u_map = NULL;

      // Create the spaces for c and d
      // - includes the internal permutation matrices for
      //  full_g to c and d
      // - includes the permutation matrices for d_l and d_u
      // c(x) = (P_c)T * g(x)
      // d(x) = (P_d)T * g(x)
      // d_L = (P_d_L)T * (P_d)T * g_l
      // d_U = (P_d_U)T * (P_d)T * g_u
      Index n_c = 0;
      Index n_d = 0;
      Index n_d_l = 0;
      Index n_d_u = 0;

      Index* c_map = new Index[n_full_g_]; // we do not know n_c yet!
      Index* d_map = new Index[n_full_g_]; // we do not know n_d yet!
      Index* d_l_map = new Index[n_full_g_]; // "
      Index* d_u_map = new Index[n_full_g_]; // "
      for (Index i=0; i<n_full_g_; i++) {
        Number lower_bound = g_l[i];
        Number upper_bound = g_u[i];
        if (lower_bound == upper_bound) {
          // equality constraint
          c_map[n_c] = i;
          n_c++;
        }
        else if (lower_bound > upper_bound) {
          char string[128];
          sprintf(string, "There are inconsistent bounds on constraint %d: lower = %25.16e and upper = %25.16e.", i, lower_bound, upper_bound);
          delete [] d_map;
          delete [] c_map;
          delete [] d_u_map;
          delete [] d_l_map;
          THROW_EXCEPTION(INVALID_TNLP, string);
        }
        else {
          // inequality constraint
          d_map[n_d] = i;
          if (lower_bound > nlp_lower_bound_inf_) {
            d_l_map[n_d_l] = n_d;
            n_d_l++;
          }
          if (upper_bound < nlp_upper_bound_inf_) {
            d_u_map[n_d_u] = n_d;
            n_d_u++;
          }
          n_d++;
        }
      }

      // create the required c_space

      if (n_x_fixed_==0 || fixed_variable_treatment_==MAKE_PARAMETER) {
        SmartPtr<DenseVectorSpace> dc_space = new DenseVectorSpace(n_c);
        c_space_ = GetRawPtr(dc_space);
        c_rhs_ = new Number[n_c];
      }
      else {
        SmartPtr<DenseVectorSpace> dc_space =
          new DenseVectorSpace(n_c+n_x_fixed_);
        c_space_ = GetRawPtr(dc_space);
        c_rhs_ = new Number[n_c+n_x_fixed_];
      }
      // create the internal expansion matrix for c to g
      P_c_g_space_ = new ExpansionMatrixSpace(n_full_g_, n_c, c_map);
      P_c_g_ = P_c_g_space_->MakeNewExpansionMatrix();
      delete [] c_map;
      c_map = NULL;

      // create the required d_space
      d_space_ = new DenseVectorSpace(n_d);
      // create the internal expansion matrix for d to g
      P_d_g_space_ = new ExpansionMatrixSpace(n_full_g_, n_d, d_map);
      P_d_g_ = P_d_g_space_->MakeNewExpansionMatrix();
      delete [] d_map;
      d_map = NULL;

      // create the required d_l space
      d_l_space_ = new DenseVectorSpace(n_d_l);
      // create the required expansion matrix for d_L to d_L_exp
      pd_l_space_ = new ExpansionMatrixSpace(n_d, n_d_l, d_l_map);
      delete [] d_l_map;
      d_l_map = NULL;

      // create the required d_u space
      d_u_space_ = new DenseVectorSpace(n_d_u);
      // create the required expansion matrix for d_U to d_U_exp
      pd_u_space_ = new ExpansionMatrixSpace(n_d, n_d_u, d_u_map);
      delete [] d_u_map;
      d_u_map = NULL;

      delete [] g_l;
      g_l = NULL;
      delete [] g_u;
      g_u = NULL;

      /** Create the matrix space for the jacobians
       */
      // Get the non zero structure
      Index* g_iRow = new Index[nz_full_jac_g_];
      Index* g_jCol = new Index[nz_full_jac_g_];
      tnlp_->eval_jac_g(n_full_x_, NULL, false, n_full_g_, nz_full_jac_g_,
                        g_iRow, g_jCol, NULL);

      if (index_style_ != TNLP::FORTRAN_STYLE) {
        for (Index i=0; i<nz_full_jac_g_; i++) {
          g_iRow[i] += 1;
          g_jCol[i] += 1;
        }
      }

      // ... build the non-zero structure for jac_c
      // ... (the permutation from rows in jac_g to jac_c is
      // ...  the same as P_c_g_)
      Index nz_jac_all;
      if (fixed_variable_treatment_==MAKE_PARAMETER) {
        nz_jac_all = nz_full_jac_g_;
      }
      else {
        nz_jac_all = nz_full_jac_g_ + n_x_fixed_;
      }
      jac_idx_map_ = new Index[nz_jac_all];
      Index* jac_c_iRow = new Index[nz_jac_all];
      Index* jac_c_jCol = new Index[nz_jac_all];
      Index current_nz = 0;
      const Index* c_row_pos = P_c_g_->CompressedPosIndices();
      if (IsValid(P_x_full_x_)) {
        // there are missing variables x
        const Index* c_col_pos = P_x_full_x_->CompressedPosIndices();
        for (Index i=0; i<nz_full_jac_g_; i++) {
          const Index c_row = c_row_pos[g_iRow[i]-1];
          const Index c_col = c_col_pos[g_jCol[i]-1];
          if (c_col != -1 && c_row != -1) {
            jac_idx_map_[current_nz] = i;
            jac_c_iRow[current_nz] = c_row + 1;
            jac_c_jCol[current_nz] = c_col + 1;
            current_nz++;
          }
        }
      }
      else {
        for (Index i=0; i<nz_full_jac_g_; i++) {
          const Index c_row = c_row_pos[g_iRow[i]-1];
          const Index c_col = g_jCol[i]-1;
          if (c_row != -1) {
            jac_idx_map_[current_nz] = i;
            jac_c_iRow[current_nz] = c_row + 1;
            jac_c_jCol[current_nz] = c_col + 1;
            current_nz++;
          }
        }
      }
      nz_jac_c_no_fixed_ = current_nz;
      Index n_added_constr;
      if (fixed_variable_treatment_==MAKE_PARAMETER) {
        nz_jac_c_ = nz_jac_c_no_fixed_;
        n_added_constr = 0;
      }
      else {
        nz_jac_c_ = nz_jac_c_no_fixed_ + n_x_fixed_;
        for (Index i=0; i<n_x_fixed_; i++) {
          jac_c_iRow[current_nz] = n_c + i + 1;
          jac_c_jCol[current_nz] = x_fixed_map_[i]+1;
          current_nz++;
        }
        n_added_constr = n_x_fixed_;
      }

      Jac_c_space_ = new GenTMatrixSpace(n_c+n_added_constr, n_x_var,
                                         nz_jac_c_, jac_c_iRow, jac_c_jCol);
      delete [] jac_c_iRow;
      jac_c_iRow = NULL;
      delete [] jac_c_jCol;
      jac_c_jCol = NULL;

      // ... build the nonzero structure for jac_d
      // ... (the permuation from rows in jac_g to jac_c is the
      // ...  the same as P_d_g_)
      Index* jac_d_iRow = new Index[nz_full_jac_g_];
      Index* jac_d_jCol = new Index[nz_full_jac_g_];
      current_nz = 0;
      const Index* d_row_pos = P_d_g_->CompressedPosIndices();
      if (IsValid(P_x_full_x_)) {
        const Index* d_col_pos = P_x_full_x_->CompressedPosIndices();
        for (Index i=0; i<nz_full_jac_g_; i++) {
          const Index d_row = d_row_pos[g_iRow[i]-1];
          const Index d_col = d_col_pos[g_jCol[i]-1];
          if (d_col != -1 && d_row != -1) {
            jac_idx_map_[current_nz + nz_jac_c_no_fixed_] = i;
            jac_d_iRow[current_nz] = d_row + 1;
            jac_d_jCol[current_nz] = d_col + 1;
            current_nz++;
          }
        }
      }
      else {
        for (Index i=0; i<nz_full_jac_g_; i++) {
          const Index d_row = d_row_pos[g_iRow[i]-1];
          const Index d_col = g_jCol[i]-1;
          if (d_row != -1) {
            jac_idx_map_[current_nz + nz_jac_c_no_fixed_] = i;
            jac_d_iRow[current_nz] = d_row + 1;
            jac_d_jCol[current_nz] = d_col + 1;
            current_nz++;
          }
        }
      }
      nz_jac_d_ = current_nz;
      Jac_d_space_ = new GenTMatrixSpace(n_d, n_x_var, nz_jac_d_, jac_d_iRow, jac_d_jCol);
      delete [] jac_d_iRow;
      jac_d_iRow = NULL;
      delete [] jac_d_jCol;
      jac_d_jCol = NULL;

      delete [] g_iRow;
      g_iRow = NULL;
      delete [] g_jCol;
      g_jCol = NULL;

      if (hessian_approximation_==EXACT) {
        /** Create the matrix space for the hessian of the lagrangian */
        Index* full_h_iRow = new Index[nz_full_h_];
        Index* full_h_jCol = new Index[nz_full_h_];
        Index* h_iRow = new Index[nz_full_h_];
        Index* h_jCol = new Index[nz_full_h_];
        bool retval =tnlp_->eval_h(n_full_x_, NULL, false, 0, n_full_g_,
                                   NULL, false,
                                   nz_full_h_, full_h_iRow, full_h_jCol, NULL);
        if (!retval) {
          jnlst_->Printf(J_ERROR, J_INITIALIZATION,
                         "Option hessian_information is not chosen as limited_memory, but eval_h returns false.\n");
          delete [] h_iRow;
          delete [] h_jCol;
          THROW_EXCEPTION(OPTION_INVALID, "eval_h is called but has not been implemented");
        }

        if (index_style_ != TNLP::FORTRAN_STYLE) {
          for (Index i=0; i<nz_full_h_; i++) {
            full_h_iRow[i] += 1;
            full_h_jCol[i] += 1;
          }
        }

        current_nz = 0;
        if (IsValid(P_x_full_x_)) {
          h_idx_map_ = new Index[nz_full_h_];
          const Index* h_pos = P_x_full_x_->CompressedPosIndices();
          for (Index i=0; i<nz_full_h_; i++) {
            const Index h_row = h_pos[full_h_iRow[i]-1];
            const Index h_col = h_pos[full_h_jCol[i]-1];
            if (h_row != -1 && h_col != -1) {
              h_idx_map_[current_nz] = i;
              h_iRow[current_nz] = h_row + 1;
              h_jCol[current_nz] = h_col + 1;
              current_nz++;
            }
          }
        }
        else {
          h_idx_map_ = NULL;
          for (Index i=0; i<nz_full_h_; i++) {
            const Index h_row = full_h_iRow[i]-1;
            const Index h_col = full_h_jCol[i]-1;
            h_iRow[i] = h_row + 1;
            h_jCol[i] = h_col + 1;
            current_nz++;
          }
          current_nz = nz_full_h_;
        }
        nz_h_ = current_nz;
        Hess_lagrangian_space_ = new SymTMatrixSpace(n_x_var, nz_h_, h_iRow, h_jCol);
        delete [] full_h_iRow;
        full_h_iRow = NULL;
        delete [] full_h_jCol;
        full_h_jCol = NULL;
        delete [] h_iRow;
        h_iRow = NULL;
        delete [] h_jCol;
        h_jCol = NULL;
      }
      else {
        nz_h_ = 0;
        Hess_lagrangian_space_ = NULL;
      }
    } /* if (warm_start_same_structure_) { */

    // Assign the spaces to the returned pointers
    x_space = x_space_;
    c_space = c_space_;
    d_space = d_space_;
    x_l_space = x_l_space_;
    px_l_space = px_l_space_;
    x_u_space = x_u_space_;
    px_u_space = px_u_space_;
    d_l_space = d_l_space_;
    pd_l_space = pd_l_space_;
    d_u_space = d_u_space_;
    pd_u_space = pd_u_space_;
    Jac_c_space = Jac_c_space_;
    Jac_d_space = Jac_d_space_;
    Hess_lagrangian_space = Hess_lagrangian_space_;

    if (IsValid(jnlst_)) {
      jnlst_->Printf(J_ITERSUMMARY, J_STATISTICS,
                     "Number of nonzeros in equality constraint Jacobian...:%9d\n", nz_jac_c_);
      jnlst_->Printf(J_ITERSUMMARY, J_STATISTICS,
                     "Number of nonzeros in inequality constraint Jacobian.:%9d\n", nz_jac_d_);
      jnlst_->Printf(J_ITERSUMMARY, J_STATISTICS,
                     "Number of nonzeros in Lagrangian Hessian.............:%9d\n\n", nz_h_);
    }

    return true;
  }

  bool TNLPAdapter::GetBoundsInformation(const Matrix& Px_L,
                                         Vector& x_L,
                                         const Matrix& Px_U,
                                         Vector& x_U,
                                         const Matrix& Pd_L,
                                         Vector& d_L,
                                         const Matrix& Pd_U,
                                         Vector& d_U)
  {
    // This could be done more efficiently, I have already called this method
    // once to setup the structure for the problem, I could store the values
    // and use them here ?
    // Actually, this is better for a warm start
    Number* x_l = new Number[n_full_x_];
    Number* x_u = new Number[n_full_x_];
    Number* g_l = new Number[n_full_g_];
    Number* g_u = new Number[n_full_g_];
    tnlp_->get_bounds_info(n_full_x_, x_l, x_u, n_full_g_, g_l, g_u);

    if (fixed_variable_treatment_==MAKE_PARAMETER) {
      // Set the values of fixed variables
      for (Index i=0; i<n_x_fixed_; i++) {
        DBG_ASSERT(x_l[x_fixed_map_[i]] == x_u[x_fixed_map_[i]]);
        full_x_[x_fixed_map_[i]] = x_l[x_fixed_map_[i]];
      }
    }

    // Set the bounds values for x
    DenseVector* dx_L = dynamic_cast<DenseVector*>(&x_L);
    DBG_ASSERT(dx_L);
    Number* values = dx_L->Values();
    const ExpansionMatrix* em_Px_L =
      dynamic_cast<const ExpansionMatrix*>(&Px_L);
    DBG_ASSERT(em_Px_L);
    if (IsValid(P_x_full_x_)) {
      for (Index i=0; i<Px_L.NCols(); i++) {
        const Index ipopt_idx = em_Px_L->ExpandedPosIndices()[i];
        const Index full_idx = P_x_full_x_->ExpandedPosIndices()[ipopt_idx];
        const Number lower_bound = x_l[full_idx];
        values[i] = lower_bound;
      }
    }
    else {
      for (Index i=0; i<Px_L.NCols(); i++) {
        const Index ipopt_idx = em_Px_L->ExpandedPosIndices()[i];
        const Number lower_bound = x_l[ipopt_idx];
        values[i] = lower_bound;
      }
    }

    DenseVector* dx_U = dynamic_cast<DenseVector*>(&x_U);
    DBG_ASSERT(dx_U);
    values = dx_U->Values();
    const ExpansionMatrix* em_Px_U =
      dynamic_cast<const ExpansionMatrix*>(&Px_U);
    DBG_ASSERT(em_Px_U);
    if (IsValid(P_x_full_x_)) {
      for (Index i=0; i<Px_U.NCols(); i++) {
        const Index ipopt_idx = em_Px_U->ExpandedPosIndices()[i];
        const Index full_idx = P_x_full_x_->ExpandedPosIndices()[ipopt_idx];
        const Number upper_bound = x_u[full_idx];
        values[i] = upper_bound;
      }
    }
    else {
      for (Index i=0; i<Px_U.NCols(); i++) {
        const Index ipopt_idx = em_Px_U->ExpandedPosIndices()[i];
        const Number upper_bound = x_u[ipopt_idx];
        values[i] = upper_bound;
      }
    }

    // get the bounds values (rhs values to subtract) for c
    // i.e. if gL == gU, then we actually have g(x) = gL = gU,
    // since we solve c(x) = 0, we actually need c(x) - gL = 0
    for (Index i=0; i<P_c_g_->NCols(); i++) {
      Index full_idx = P_c_g_->ExpandedPosIndices()[i];
      Number rhs = g_l[full_idx];
      c_rhs_[i] = rhs;
    }
    // similarly, if we have fixed variables, consider them here
    if (fixed_variable_treatment_==MAKE_CONSTRAINT) {
      Index n_c_no_fixed = P_c_g_->NCols();
      for (Index i=0; i<n_x_fixed_; i++) {
        DBG_ASSERT(x_l[x_fixed_map_[i]]==x_u[x_fixed_map_[i]]);
        c_rhs_[n_c_no_fixed+i] = x_l[x_fixed_map_[i]];
      }
    }

    // get the bounds values for d
    DenseVector* dd_L = dynamic_cast<DenseVector*>(&d_L);
    DBG_ASSERT(dd_L);
    values = dd_L->Values();
    const ExpansionMatrix* em_Pd_L =
      dynamic_cast<const ExpansionMatrix*>(&Pd_L);
    DBG_ASSERT(em_Pd_L);
    for (Index i=0; i<Pd_L.NCols(); i++) {
      Index d_exp_idx = em_Pd_L->ExpandedPosIndices()[i];
      Index full_idx = P_d_g_->ExpandedPosIndices()[d_exp_idx];
      Number lower_bound = g_l[full_idx];
      values[i] = lower_bound;
    }

    DenseVector* dd_U = dynamic_cast<DenseVector*>(&d_U);
    DBG_ASSERT(dd_U);
    values = dd_U->Values();
    const ExpansionMatrix* em_Pd_U =
      dynamic_cast<const ExpansionMatrix*>(&Pd_U);
    DBG_ASSERT(em_Pd_U);
    for (Index i=0; i<Pd_U.NCols(); i++) {
      Index d_exp_idx = em_Pd_U->ExpandedPosIndices()[i];
      Index full_idx = P_d_g_->ExpandedPosIndices()[d_exp_idx];
      Number upper_bound = g_u[full_idx];
      values[i] = upper_bound;
    }

    delete [] x_l;
    x_l = NULL;
    delete [] x_u;
    x_u = NULL;
    delete [] g_l;
    g_l = NULL;
    delete [] g_u;
    g_u = NULL;

    return true;
  }

  bool TNLPAdapter::GetStartingPoint(SmartPtr<Vector> x,
                                     bool need_x,
                                     SmartPtr<Vector> y_c,
                                     bool need_y_c,
                                     SmartPtr<Vector> y_d,
                                     bool need_y_d,
                                     SmartPtr<Vector> z_L,
                                     bool need_z_L,
                                     SmartPtr<Vector> z_U,
                                     bool need_z_U
                                    )
  {
    Number* full_x = new Number[n_full_x_];
    Number* full_z_l = new Number[n_full_x_];
    Number* full_z_u = new Number[n_full_x_];
    Number* full_lambda = new Number[n_full_g_];
    bool init_x = need_x;
    bool init_z = need_z_L && need_z_U;
    bool init_lambda = need_y_c && need_y_d;

    bool retvalue =
      tnlp_->get_starting_point(n_full_x_, init_x, full_x, init_z, full_z_l, full_z_u, n_full_g_, init_lambda, full_lambda);

    if (!retvalue) {
      return false;
    }

    if (need_x) {
      DenseVector* dx = dynamic_cast<DenseVector*>(GetRawPtr(x));
      DBG_ASSERT(dx);
      Number* values = dx->Values();
      if (IsValid(P_x_full_x_)) {
        const Index* x_pos = P_x_full_x_->ExpandedPosIndices();
        for (Index i=0; i<x->Dim(); i++) {
          values[i] = full_x[x_pos[i]];
        }
      }
      else {
        IpBlasDcopy(x->Dim(), full_x, 1, values, 1);
      }
    }

    if (need_y_c) {
      DenseVector* dy_c = dynamic_cast<DenseVector*>(GetRawPtr(y_c));
      DBG_ASSERT(dy_c);
      Number* values = dy_c->Values();
      const Index* y_c_pos = P_c_g_->ExpandedPosIndices();
      for (Index i=0; i<P_c_g_->NCols(); i++) {
        values[i] = full_lambda[y_c_pos[i]];
      }
      if (fixed_variable_treatment_==MAKE_CONSTRAINT) {
        // ToDo maybe use info from z_L and Z_U here?
        const Number zero = 0.;
        IpBlasDcopy(n_x_fixed_, &zero, 0, &values[P_c_g_->NCols()], 1);
      }
    }

    if (need_y_d) {
      DenseVector* dy_d = dynamic_cast<DenseVector*>(GetRawPtr(y_d));
      DBG_ASSERT(dy_d);
      Number* values = dy_d->Values();
      const Index* y_d_pos = P_d_g_->ExpandedPosIndices();
      for (Index i=0; i<y_d->Dim(); i++) {
        values[i] = full_lambda[y_d_pos[i]];
      }
    }

    if (need_z_L) {
      DenseVector* dz_l = dynamic_cast<DenseVector*>(GetRawPtr(z_L));
      DBG_ASSERT(dz_l);
      Number* values = dz_l->Values();
      const Index* z_l_pos = P_x_x_L_->ExpandedPosIndices();
      if (IsValid(P_x_full_x_)) {
        const Index* x_pos = P_x_full_x_->ExpandedPosIndices();
        for (Index i=0; i<z_L->Dim(); i++) {
          Index idx = z_l_pos[i]; // convert from x_L to x (ipopt)
          idx = x_pos[idx]; // convert from x (ipopt) to x_full
          values[i] = full_z_l[idx];
        }
      }
      else {
        for (Index i=0; i<z_L->Dim(); i++) {
          Index idx = z_l_pos[i]; // convert from x_L to x (ipopt)
          values[i] = full_z_l[idx];
        }
      }
    }

    if (need_z_U) {
      DenseVector* dz_u = dynamic_cast<DenseVector*>(GetRawPtr(z_U));
      DBG_ASSERT(dz_u);
      Number* values = dz_u->Values();
      const Index* z_u_pos = P_x_x_U_->ExpandedPosIndices();
      if (IsValid(P_x_full_x_)) {
        const Index* x_pos = P_x_full_x_->ExpandedPosIndices();
        for (Index i=0; i<z_U->Dim(); i++) {
          Index idx = z_u_pos[i]; // convert from x_u to x (ipopt)
          idx = x_pos[idx]; // convert from x (ipopt) to x_full
          values[i] = full_z_u[idx];
        }
      }
      else {
        for (Index i=0; i<z_U->Dim(); i++) {
          Index idx = z_u_pos[i]; // convert from x_u to x (ipopt)
          values[i] = full_z_u[idx];
        }
      }
    }

    delete [] full_x;
    full_x = NULL;
    delete [] full_z_l;
    full_z_l = NULL;
    delete [] full_z_u;
    full_z_u = NULL;
    delete [] full_lambda;
    full_lambda = NULL;

    return true;
  }

  bool TNLPAdapter::GetWarmStartIterate(IteratesVector& warm_start_iterate)
  {
    return tnlp_->get_warm_start_iterate(warm_start_iterate);
  }

  bool TNLPAdapter::Eval_f(const Vector& x, Number& f)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }
    return tnlp_->eval_f(n_full_x_, full_x_, new_x, f);
  }

  bool TNLPAdapter::Eval_grad_f(const Vector& x, Vector& g_f)
  {
    bool retvalue = false;
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    if (IsValid(P_x_full_x_)) {
      Number* full_grad_f = new Number[n_full_x_];
      if (tnlp_->eval_grad_f(n_full_x_, full_x_, new_x, full_grad_f)) {
        const Index* x_pos = P_x_full_x_->ExpandedPosIndices();
        DenseVector* dg_f = dynamic_cast<DenseVector*>(&g_f);
        DBG_ASSERT(dg_f);
        Number* values = dg_f->Values();

        for (Index i=0; i<g_f.Dim(); i++) {
          values[i] = full_grad_f[x_pos[i]];
        }
        retvalue = true;
      }
      delete [] full_grad_f;
    }
    else {
      DenseVector* dg_f = dynamic_cast<DenseVector*>(&g_f);
      DBG_ASSERT(dg_f);
      Number* values = dg_f->Values();
      retvalue = tnlp_->eval_grad_f(n_full_x_, full_x_, new_x, values);
    }

    return retvalue;
  }

  bool TNLPAdapter::Eval_c(const Vector& x, Vector& c)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    DenseVector* dc = dynamic_cast<DenseVector*>(&c);
    DBG_ASSERT(dc);
    Number* values = dc->Values();
    if (internal_eval_g(new_x)) {
      const Index* c_pos = P_c_g_->ExpandedPosIndices();
      Index n_c_no_fixed = P_c_g_->NCols();
      for (Index i=0; i<n_c_no_fixed; i++) {
        values[i] = full_g_[c_pos[i]];
        values[i] -= c_rhs_[i];
      }
      if (fixed_variable_treatment_==MAKE_CONSTRAINT) {
        for (Index i=0; i<n_x_fixed_; i++) {
          values[n_c_no_fixed+i] =
            full_x_[x_fixed_map_[i]] - c_rhs_[n_c_no_fixed+i];
        }
      }
      return true;
    }

    return false;
  }

  bool TNLPAdapter::Eval_jac_c(const Vector& x, Matrix& jac_c)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    if (internal_eval_jac_g(new_x)) {
      GenTMatrix* gt_jac_c = dynamic_cast<GenTMatrix*>(&jac_c);
      DBG_ASSERT(gt_jac_c);
      Number* values = gt_jac_c->Values();

      for (Index i=0; i<nz_jac_c_no_fixed_; i++) {
        // Assume the same structure as initially given
        values[i] = jac_g_[jac_idx_map_[i]];
      }
      if (fixed_variable_treatment_==MAKE_CONSTRAINT) {
        const Number one = 1.;
        IpBlasDcopy(n_x_fixed_, &one, 0, &values[nz_jac_c_no_fixed_], 1);
      }
      return true;
    }
    return false;
  }

  bool TNLPAdapter::Eval_d(const Vector& x, Vector& d)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    DenseVector* dd = dynamic_cast<DenseVector*>(&d);
    DBG_ASSERT(dd);
    Number* values = dd->Values();
    if (internal_eval_g(new_x)) {
      const Index* d_pos = P_d_g_->ExpandedPosIndices();
      for (Index i=0; i<d.Dim(); i++) {
        values[i] = full_g_[d_pos[i]];
      }
      return true;
    }

    return false;
  }

  bool TNLPAdapter::Eval_jac_d(const Vector& x, Matrix& jac_d)
  {
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }

    if (internal_eval_jac_g(new_x)) {
      GenTMatrix* gt_jac_d = dynamic_cast<GenTMatrix*>(&jac_d);
      DBG_ASSERT(gt_jac_d);
      Number* values = gt_jac_d->Values();

      for (Index i=0; i<nz_jac_d_; i++) {
        // Assume the same structure as initially given
        values[i] = jac_g_[jac_idx_map_[nz_jac_c_no_fixed_ + i]];
      }
      return true;
    }
    return false;
  }

  bool TNLPAdapter::Eval_h(const Vector& x,
                           Number obj_factor,
                           const Vector& yc,
                           const Vector& yd,
                           SymMatrix& h)
  {
    // First see if all weights are set to zero (for example, when
    // computing the least square multiplier estimates, this is what
    // we do).  In that case, there is no need to compute values, just
    // set them to zero.
    if (obj_factor==0. && yc.Asum()==0. && yd.Asum()==0.) {
      SymTMatrix* st_h = dynamic_cast<SymTMatrix*>(&h);
      DBG_ASSERT(st_h);
      Number* values = st_h->Values();
      for (Index i=0; i<nz_h_; i++) {
        values[i] = 0.;
      }
      return true;
    }

    bool retval = false;
    bool new_x = false;
    if (update_local_x(x)) {
      new_x = true;
    }
    bool new_y = false;
    if (update_local_lambda(yc, yd)) {
      new_y = true;
    }

    SymTMatrix* st_h = dynamic_cast<SymTMatrix*>(&h);
    DBG_ASSERT(st_h);
    Number* values = st_h->Values();

    if (h_idx_map_) {
      Number* full_h = new Number[nz_full_h_];

      if (tnlp_->eval_h(n_full_x_, full_x_, new_x, obj_factor, n_full_g_,
                        full_lambda_, new_y, nz_full_h_, NULL, NULL, full_h)) {
        for (Index i=0; i<nz_h_; i++) {
          values[i] = full_h[h_idx_map_[i]];
        }
        retval = true;
      }
      delete [] full_h;
    }
    else {
      retval = tnlp_->eval_h(n_full_x_, full_x_, new_x, obj_factor, n_full_g_,
                             full_lambda_, new_y, nz_full_h_, NULL, NULL,
                             values);
    }

    return retval;
  }

  void TNLPAdapter::GetScalingParameters(
    const SmartPtr<const VectorSpace> x_space,
    const SmartPtr<const VectorSpace> c_space,
    const SmartPtr<const VectorSpace> d_space,
    Number& obj_scaling,
    SmartPtr<Vector>& x_scaling,
    SmartPtr<Vector>& c_scaling,
    SmartPtr<Vector>& d_scaling) const
  {
    x_scaling = x_space->MakeNew();
    c_scaling = c_space->MakeNew();
    d_scaling = d_space->MakeNew();
    DBG_ASSERT((c_scaling->Dim()+d_scaling->Dim()) == n_full_g_);

    DenseVector* dx = dynamic_cast<DenseVector*>(GetRawPtr(x_scaling));
    DenseVector* dc = dynamic_cast<DenseVector*>(GetRawPtr(c_scaling));
    DenseVector* dd = dynamic_cast<DenseVector*>(GetRawPtr(d_scaling));
    DBG_ASSERT(dx && dc && dd);
    Number* dx_values = dx->Values();
    Number* dc_values = dc->Values();
    Number* dd_values = dd->Values();

    Number* full_g_scaling = new Number[n_full_g_];
    bool use_x_scaling=true;
    bool use_g_scaling=true;

    if (IsValid(P_x_full_x_)) {
      Number* full_x_scaling = new Number[n_full_x_];
      bool retval = tnlp_->get_scaling_parameters(obj_scaling,
                    use_x_scaling, n_full_x_,
                    full_x_scaling,
                    use_g_scaling, n_full_g_,
                    full_g_scaling);
      if (!retval) {
        jnlst_->Printf(J_ERROR, J_INITIALIZATION,
                       "Option nlp_scaling_method selected as user-scaling, but no user-scaling available, or it cannot be computed.\n");
        THROW_EXCEPTION(OPTION_INVALID,
                        "User scaling chosen, but get_scaling_parameters returned false.");
      }

      if (use_x_scaling) {
        const Index* x_pos = P_x_full_x_->ExpandedPosIndices();
        for (Index i=0; i<dx->Dim(); i++) {
          dx_values[i] = full_x_scaling[x_pos[i]];
        }
      }
      delete [] full_x_scaling;
    }
    else {
      bool retval = tnlp_->get_scaling_parameters(obj_scaling,
                    use_x_scaling, n_full_x_,
                    dx_values,
                    use_g_scaling, n_full_g_,
                    full_g_scaling);
      if (!retval) {
        jnlst_->Printf(J_ERROR, J_INITIALIZATION,
                       "Option nlp_scaling_method selected as user-scaling, but no user-scaling available, or it cannot be computed.\n");
        THROW_EXCEPTION(OPTION_INVALID,
                        "User scaling chosen, but get_scaling_parameters returned false.");
      }
    }

    if (!use_x_scaling) {
      x_scaling = NULL;
    }

    if (use_g_scaling) {
      const Index* c_pos = P_c_g_->ExpandedPosIndices();
      for (Index i=0; i<P_c_g_->NCols(); i++) {
        dc_values[i] = full_g_scaling[c_pos[i]];
      }
      if (fixed_variable_treatment_==MAKE_CONSTRAINT) {
        const Number one = 1.;
        IpBlasDcopy(n_x_fixed_, &one, 0, &dc_values[P_c_g_->NCols()], 1);
      }

      const Index* d_pos = P_d_g_->ExpandedPosIndices();
      for (Index i=0; i<dd->Dim(); i++) {
        dd_values[i] = full_g_scaling[d_pos[i]];
      }
    }
    else {
      c_scaling = NULL;
      d_scaling = NULL;
    }

    delete [] full_g_scaling;
  }


  void TNLPAdapter::FinalizeSolution(SolverReturn status,
                                     const Vector& x, const Vector& z_L, const Vector& z_U,
                                     const Vector& c, const Vector& d,
                                     const Vector& y_c, const Vector& y_d,
                                     Number obj_value)
  {
    DBG_START_METH("TNLPAdapter::FinalizeSolution", dbg_verbosity);

    update_local_x(x);
    update_local_lambda(y_c, y_d);

    ResortX(x, full_x_);
    ResortG(y_c, y_d, full_lambda_);

    Number* full_g = new Number[n_full_g_];
    ResortG(c, d, full_g);

    Number* full_z_L = new Number[n_full_x_];
    Number* full_z_U = new Number[n_full_x_];
    for (int i=0; i<n_full_x_; i++) {
      full_z_L[i] = nlp_lower_bound_inf_;
      full_z_U[i] = nlp_upper_bound_inf_;
    }
    ResortBnds(z_L, full_z_L, z_U, full_z_U);

    // Hopefully the following is correct to recover the bound
    // multipliers for fixed variables (sign ok?)
    if (fixed_variable_treatment_==MAKE_CONSTRAINT && n_x_fixed_>0) {
      const DenseVector* dy_c = dynamic_cast<const DenseVector*>(&y_c);
      DBG_ASSERT(dy_c);
      DBG_ASSERT(!dy_c->IsHomogeneous());
      const Number* values = dy_c->Values();
      Index n_c_no_fixed = y_c.Dim() - n_x_fixed_;
      for (Index i=0; i<n_x_fixed_; i++) {
        full_z_L[x_fixed_map_[i]] = Max(0., -values[n_c_no_fixed+i]);
        full_z_U[x_fixed_map_[i]] = Max(0., values[n_c_no_fixed+i]);
      }
    }

    tnlp_->finalize_solution(status,
                             n_full_x_, full_x_, full_z_L, full_z_U,
                             n_full_g_, full_g, full_lambda_,
                             obj_value);

    delete [] full_z_L;
    full_z_L = NULL;
    delete [] full_z_U;
    full_z_U = NULL;
    delete [] full_g;
    full_g = NULL;
  }

  bool TNLPAdapter::
  IntermediateCallBack(AlgorithmMode mode,
                       Index iter, Number obj_value,
                       Number inf_pr, Number inf_du,
                       Number mu, Number d_norm,
                       Number regularization_size,
                       Number alpha_du, Number alpha_pr,
                       Index ls_trials,
                       const IpoptData* ip_data,
                       IpoptCalculatedQuantities* ip_cq)
  {
    return tnlp_->intermediate_callback(mode, iter, obj_value, inf_pr, inf_du,
                                        mu, d_norm, regularization_size,
                                        alpha_du, alpha_pr, ls_trials,
                                        ip_data, ip_cq);
  }


  void TNLPAdapter::
  GetQuasiNewtonApproximationSpaces(SmartPtr<VectorSpace>& approx_space,
                                    SmartPtr<Matrix>& P_approx)
  {
    Index num_nonlin_vars =
      tnlp_->get_number_of_nonlinear_variables();

    if (num_nonlin_vars<0) {
      approx_space = NULL;
      P_approx = NULL;
    }
    else {
      Index* pos_nonlin_vars = new Index[num_nonlin_vars];
      if (num_nonlin_vars>0) {
        bool retval = tnlp_->get_list_of_nonlinear_variables(num_nonlin_vars,
                      pos_nonlin_vars);
        if (!retval) {
          jnlst_->Printf(J_ERROR, J_INITIALIZATION,
                         "TNLP's get_number_of_nonlinear_variables returns non-negative number, but get_list_of_nonlinear_variables returns false.\n");
          THROW_EXCEPTION(INVALID_TNLP, "get_list_of_nonlinear_variables has not been overwritten");
        }
      }

      // Correct indices in case user starts counting variables at 1
      // and not 0
      if (index_style_ == TNLP::FORTRAN_STYLE) {
        for (Index i=0; i<num_nonlin_vars; i++) {
          pos_nonlin_vars[i]--;
        }
      }

      if (IsNull(P_x_full_x_)) {
        if (num_nonlin_vars == n_full_x_) {
          approx_space = NULL;
          P_approx = NULL;
        }
        else {
          SmartPtr<ExpansionMatrixSpace> ex_sp =
            new ExpansionMatrixSpace(n_full_x_, num_nonlin_vars,
                                     pos_nonlin_vars);
          P_approx = ex_sp->MakeNew();
          approx_space = new DenseVectorSpace(num_nonlin_vars);
        }
      }
      else {
        const Index* compr_pos = P_x_full_x_->CompressedPosIndices();
        Index* nonfixed_pos_nonlin_vars = new Index[num_nonlin_vars];

        Index nonfixed_nonlin_vars = 0;
        for (Index i=0; i<num_nonlin_vars; i++) {
          Index full_pos = pos_nonlin_vars[i];
          Index nonfixed_pos = compr_pos[full_pos];
          if (nonfixed_pos>=0) {
            nonfixed_pos_nonlin_vars[nonfixed_nonlin_vars] = nonfixed_pos;
            nonfixed_nonlin_vars++;
          }
        }

        const Index n_x_free = n_full_x_ - n_x_fixed_;
        if (nonfixed_nonlin_vars == n_x_free) {
          approx_space = NULL;
          P_approx = NULL;
        }
        else {
          SmartPtr<ExpansionMatrixSpace> ex_sp =
            new ExpansionMatrixSpace(n_x_free, nonfixed_nonlin_vars,
                                     nonfixed_pos_nonlin_vars);
          P_approx = ex_sp->MakeNew();
          approx_space = new DenseVectorSpace(nonfixed_nonlin_vars);
        }

        delete [] nonfixed_pos_nonlin_vars;
      }
      delete [] pos_nonlin_vars;
    }
  }

  void TNLPAdapter::ResortX(const Vector& x, Number* x_orig)
  {
    const DenseVector* dx = dynamic_cast<const DenseVector*>(&x);
    DBG_ASSERT(dx);

    if (IsValid(P_x_full_x_)) {
      const Index* x_pos = P_x_full_x_->CompressedPosIndices();

      if (dx->IsHomogeneous()) {
        Number scalar = dx->Scalar();
        for (Index i=0; i<n_full_x_; i++) {
          Index idx = x_pos[i];
          if (idx != -1) {
            x_orig[i] = scalar;
          }
          else {
            x_orig[i] = full_x_[i];
          }
        }
      }
      else {
        const Number* x_values = dx->Values();
        for (Index i=0; i<n_full_x_; i++) {
          Index idx = x_pos[i];
          if (idx != -1) {
            x_orig[i] = x_values[idx];
          }
          else {
            x_orig[i] = full_x_[i];
          }
        }
      }
    }
    else {
      if (dx->IsHomogeneous()) {
        Number scalar = dx->Scalar();
        IpBlasDcopy(n_full_x_, &scalar, 0, x_orig, 1);
      }
      else {
        IpBlasDcopy(n_full_x_, dx->Values(), 1, x_orig, 1);
      }
    }
  }

  void TNLPAdapter::ResortG(const Vector& c, const Vector& d, Number* g_orig)
  {
    const DenseVector* dc = dynamic_cast<const DenseVector*>(&c);
    DBG_ASSERT(dc);

    const Index* c_pos = P_c_g_->ExpandedPosIndices();
    if (dc->IsHomogeneous()) {
      Number scalar = dc->Scalar();
      for (Index i=0; i<P_c_g_->NCols(); i++) {
        g_orig[c_pos[i]] = scalar;
      }
    }
    else {
      const Number* c_values = dc->Values();
      for (Index i=0; i<P_c_g_->NCols(); i++) {
        g_orig[c_pos[i]] = c_values[i];
      }
    }

    const DenseVector* dd = dynamic_cast<const DenseVector*>(&d);
    DBG_ASSERT(dd);

    const Index* d_pos = P_d_g_->ExpandedPosIndices();
    if (dd->IsHomogeneous()) {
      Number scalar = dd->Scalar();
      for (Index i=0; i<d.Dim(); i++) {
        g_orig[d_pos[i]] = scalar;
      }
    }
    else {
      const Number* d_values = dd->Values();
      for (Index i=0; i<d.Dim(); i++) {
        g_orig[d_pos[i]] = d_values[i];
      }
    }
  }

  void TNLPAdapter::ResortBnds(const Vector& x_L, Number* x_L_orig,
                               const Vector& x_U, Number* x_U_orig)
  {
    if (x_L_orig) {
      const DenseVector* dx_L = dynamic_cast<const DenseVector*>(&x_L);
      DBG_ASSERT(dx_L);

      const Index* bnds_pos_not_fixed = P_x_x_L_->ExpandedPosIndices();

      if (IsValid(P_x_full_x_)) {
        const Index* bnds_pos_full = P_x_full_x_->ExpandedPosIndices();
        if (dx_L->IsHomogeneous()) {
          Number scalar = dx_L->Scalar();
          for (Index i=0; i<x_L.Dim(); i++) {
            int idx = bnds_pos_not_fixed[i];
            idx = bnds_pos_full[idx];
            x_L_orig[idx] = scalar;
          }
        }
        else {
          const Number* x_L_values = dx_L->Values();
          for (Index i=0; i<x_L.Dim(); i++) {
            int idx = bnds_pos_not_fixed[i];
            idx = bnds_pos_full[idx];
            x_L_orig[idx] = x_L_values[i];
          }
        }
      }
      else {
        if (dx_L->IsHomogeneous()) {
          Number scalar = dx_L->Scalar();
          for (Index i=0; i<x_L.Dim(); i++) {
            int idx = bnds_pos_not_fixed[i];
            x_L_orig[idx] = scalar;
          }
        }
        else {
          const Number* x_L_values = dx_L->Values();
          for (Index i=0; i<x_L.Dim(); i++) {
            int idx = bnds_pos_not_fixed[i];
            x_L_orig[idx] = x_L_values[i];
          }
        }
      }
    }

    if (x_U_orig) {
      const DenseVector* dx_U = dynamic_cast<const DenseVector*>(&x_U);
      DBG_ASSERT(dx_U);

      const Index* bnds_pos_not_fixed = P_x_x_U_->ExpandedPosIndices();

      if (IsValid(P_x_full_x_)) {
        const Index* bnds_pos_full = P_x_full_x_->ExpandedPosIndices();
        if (dx_U->IsHomogeneous()) {
          Number scalar = dx_U->Scalar();
          for (Index i=0; i<x_U.Dim(); i++) {
            int idx = bnds_pos_not_fixed[i];
            idx = bnds_pos_full[idx];
            x_U_orig[idx] = scalar;
          }
        }
        else {
          const Number* x_U_values = dx_U->Values();
          for (Index i=0; i<x_U.Dim(); i++) {
            int idx = bnds_pos_not_fixed[i];
            idx = bnds_pos_full[idx];
            x_U_orig[idx] = x_U_values[i];
          }
        }
      }
      else {
        if (dx_U->IsHomogeneous()) {
          Number scalar = dx_U->Scalar();
          for (Index i=0; i<x_U.Dim(); i++) {
            int idx = bnds_pos_not_fixed[i];
            x_U_orig[idx] = scalar;
          }
        }
        else {
          const Number* x_U_values = dx_U->Values();
          for (Index i=0; i<x_U.Dim(); i++) {
            int idx = bnds_pos_not_fixed[i];
            x_U_orig[idx] = x_U_values[i];
          }
        }
      }
    }
  }

  bool TNLPAdapter::update_local_x(const Vector& x)
  {
    if (x.GetTag() == x_tag_for_iterates_) {
      return false;
    }

    ResortX(x, full_x_);

    x_tag_for_iterates_ = x.GetTag();

    return true;
  }

  bool TNLPAdapter::update_local_lambda(const Vector& y_c, const Vector& y_d)
  {
    if (y_c.GetTag() == y_c_tag_for_iterates_
        && y_d.GetTag() == y_d_tag_for_iterates_) {
      return false;
    }

    ResortG(y_c, y_d, full_lambda_);

    y_c_tag_for_iterates_ = y_c.GetTag();
    y_d_tag_for_iterates_ = y_d.GetTag();

    return true;
  }

  bool TNLPAdapter::internal_eval_g(bool new_x)
  {
    if (x_tag_for_g_ == x_tag_for_iterates_) {
      // already calculated!
      return true;
    }

    x_tag_for_g_ = x_tag_for_iterates_;

    bool retval = tnlp_->eval_g(n_full_x_, full_x_, new_x, n_full_g_, full_g_);

    if (!retval) {
      x_tag_for_jac_g_ = 0;
    }

    return retval;
  }

  bool TNLPAdapter::internal_eval_jac_g(bool new_x)
  {
    if (x_tag_for_jac_g_ == x_tag_for_iterates_) {
      // already calculated!
      return true;
    }

    x_tag_for_jac_g_ = x_tag_for_iterates_;

    bool retval = tnlp_->eval_jac_g(n_full_x_, full_x_, new_x, n_full_g_,
                                    nz_full_jac_g_, NULL, NULL, jac_g_);

    if (!retval) {
      x_tag_for_jac_g_ = 0;
    }

    return retval;
  }

  bool TNLPAdapter::CheckDerivatives(TNLPAdapter::DerivativeTestEnum deriv_test)
  {
    if (deriv_test == NO_TEST) {
      return true;
    }

    Index nerrors = 0;

    ASSERT_EXCEPTION(IsValid(jnlst_), ERROR_IN_TNLP_DERIVATIVE_TEST,
                     "No Journalist given to TNLPAdapter.  Need Journalist, otherwise can't produce any output in DerivativeChecker!");

    jnlst_->Printf(J_SUMMARY, J_NLP,
                   "\nStarting derivative checker.\n\n");

    bool retval = true;
    // Since this method should be independent of all other internal
    // data (so that it can be called indpenendent of GetSpace etc),
    // we are not using any internal fields

    // Obtain the problem size
    Index nx; // number of variables
    Index ng; // number of constraints
    Index nz_jac_g; // number of nonzeros in constraint Jacobian
    Index nz_hess_lag; // number of nonzeros in Lagrangian Hessian
    TNLP::IndexStyleEnum index_style;
    tnlp_->get_nlp_info(nx, ng, nz_jac_g, nz_hess_lag, index_style);

    // Obtain starting point as reference point at which derivative
    // test should be performed
    Number* xref = new Number[nx];
    tnlp_->get_starting_point(nx, true, xref, false, NULL, NULL, ng, false, NULL);

    // Obtain value of objective and constraints at reference point
    bool new_x = true;
    Number fref;
    Number* gref = NULL;
    if (ng>0) {
      gref = new Number[ng];
    }
    retval = tnlp_->eval_f(nx, xref, new_x, fref);
    ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                     "In TNLP derivative test: f could not be evaluated at reference point.");
    new_x = false;
    if (ng>0) {
      retval = tnlp_->eval_g(nx, xref, new_x, ng, gref);
      ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                       "In TNLP derivative test: g could not be evaluated at reference point.");
    }

    // Obtain gradient of objective function at reference pont
    Number* grad_f = new Number[nx];
    retval = tnlp_->eval_grad_f(nx, xref, true, grad_f);
    ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                     "In TNLP derivative test: grad_f could not be evaluated at reference point.");

    Index* g_iRow = NULL;
    Index* g_jCol = NULL;
    Number* jac_g = NULL;
    if (ng>0) {
      // Obtain constraint Jacobian at reference point (including structure)
      g_iRow = new Index[nz_jac_g];
      g_jCol = new Index[nz_jac_g];
      retval = tnlp_->eval_jac_g(nx, NULL, false, ng, nz_jac_g,
                                 g_iRow, g_jCol, NULL);
      ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                       "In TNLP derivative test: Jacobian structure could not be evaluated.");
      // Correct counting if required to C-style
      if (index_style == TNLP::FORTRAN_STYLE) {
        for (Index i=0; i<nz_jac_g; i++) {
          g_iRow[i] -= 1;
          g_jCol[i] -= 1;
        }
      }
      // Obtain values at reference pont
      jac_g = new Number[nz_jac_g];
      retval = tnlp_->eval_jac_g(nx, xref, new_x, ng,
                                 nz_jac_g, NULL, NULL, jac_g);
      ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                       "In TNLP derivative test: Jacobian values could not be evaluated at reference point.");
    }

    // Space for the perturbed point
    Number* xpert = new Number[nx];
    IpBlasDcopy(nx, xref, 1, xpert, 1);

    // Space for constraints at perturbed point
    Number* gpert = NULL;
    if (ng>0) {
      gpert = new Number[ng];
    }

    Index index_correction = 0;
    if (index_style == TNLP::FORTRAN_STYLE) {
      index_correction = 1;
    }

    // Now go through all variables and check the partial derivatives
    for (Index ivar=0; ivar<nx; ivar++) {
      Number this_perturbation =
        derivative_test_perturbation_*Max(1.,fabs(xref[ivar]));
      xpert[ivar] = xref[ivar] + this_perturbation;

      Number fpert;
      new_x = true;
      retval = tnlp_->eval_f(nx, xpert, new_x, fpert);
      ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                       "In TNLP derivative test: f could not be evaluated at perturbed point.");
      new_x = false;

      Number deriv_approx = (fpert - fref)/this_perturbation;
      Number deriv_exact = grad_f[ivar];
      Number rel_error =
        fabs(deriv_approx-deriv_exact)/Max(fabs(deriv_approx),1.);
      char cflag=' ';
      if (rel_error >= derivative_test_tol_) {
        cflag='*';
        nerrors++;
      }
      if (cflag != ' ' || derivative_test_print_all_) {
        jnlst_->Printf(J_WARNING, J_NLP,
                       "%c grad_f[      %5d] = %23.16e    ~ %23.16e  [%10.3e]\n",
                       cflag, ivar+index_correction,
                       deriv_exact, deriv_approx, rel_error);
      }

      if (ng>0) {
        retval = tnlp_->eval_g(nx, xpert, new_x, ng, gpert);
        ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                         "In TNLP derivative test: g could not be evaluated at reference point.");
        for (Index icon=0; icon<ng; icon++) {
          deriv_approx = (gpert[icon] - gref[icon])/this_perturbation;
          deriv_exact = 0.;
          bool found = false;
          for (Index i=0; i<nz_jac_g; i++) {
            if (g_iRow[i]==icon && g_jCol[i]==ivar) {
              found = true;
              deriv_exact += jac_g[i];
            }
          }

          rel_error = fabs(deriv_approx-deriv_exact)/Max(fabs(deriv_approx),1.);
          cflag=' ';
          if (rel_error >= derivative_test_tol_) {
            cflag='*';
            nerrors++;
          }
          char sflag=' ';
          if (found) {
            sflag = 'v';
          }
          if (cflag != ' ' || derivative_test_print_all_) {
            jnlst_->Printf(J_WARNING, J_NLP,
                           "%c jac_g [%5d,%5d] = %23.16e %c  ~ %23.16e  [%10.3e]\n",
                           cflag, icon+index_correction, ivar+index_correction,
                           deriv_exact, sflag, deriv_approx, rel_error);
          }
        }
      }

      xpert[ivar] = xref[ivar];
    }

    const Number zero = 0.;
    if (deriv_test == SECOND_ORDER_TEST) {
      // Get sparsity structure of Hessian
      Index* h_iRow = new Index[nz_hess_lag];
      Index* h_jCol = new Index[nz_hess_lag];
      retval = tnlp_->eval_h(nx, NULL, false, 0., ng, NULL, false,
                             nz_hess_lag, h_iRow, h_jCol, NULL);
      ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                       "In TNLP derivative test: Hessian structure could not be evaluated.");

      if (index_style == TNLP::FORTRAN_STYLE) {
        for (Index i=0; i<nz_hess_lag; i++) {
          h_iRow[i] -= 1;
          h_jCol[i] -= 1;
        }
      }
      Number* h_values = new Number[nz_hess_lag];

      Number* lambda = NULL;
      if (ng>0) {
        lambda = new Number[ng];
        IpBlasDcopy(ng, &zero, 0, lambda, 1);
      }
      Number* gradref = new Number[nx]; // gradient of objective or constraint at reference point
      Number* gradpert = new Number[nx]; // gradient of objective or constraint at perturbed point
      Number* jacpert = new Number[nz_jac_g];

      // Check all Hessians
      for (Index icon=-1; icon<ng; icon++) {
        Number objfact = 0.;
        if (icon == -1) {
          objfact = 1.;
          IpBlasDcopy(nx, grad_f, 1, gradref, 1);
        }
        else {
          lambda[icon] = 1.;
          IpBlasDcopy(nx, &zero, 0, gradref, 1);
          for (Index i=0; i<nz_jac_g; i++) {
            if (g_iRow[i]==icon) {
              gradref[g_jCol[i]] += jac_g[i];
            }
          }
        }
        // Hessian at reference point
        new_x = true;
        bool new_y = true;
        retval = tnlp_->eval_h(nx, xref, new_x, objfact, ng, lambda, new_y,
                               nz_hess_lag, NULL, NULL, h_values);
        new_x = false;
        new_y = false;
        ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                         "In TNLP derivative test: Hessian could not be evaluated at reference point.");

        for (Index ivar=0; ivar<nx; ivar++) {
          Number this_perturbation =
            derivative_test_perturbation_*Max(1.,fabs(xref[ivar]));
          xpert[ivar] = xref[ivar] + this_perturbation;

          new_x = true;
          if (icon==-1) {
            // we are looking at the objective function
            retval = tnlp_->eval_grad_f(nx, xpert, new_x, gradpert);
            ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                             "In TNLP derivative test: grad_f could not be evaluated at perturbed point.");
          }
          else {
            // this is the icon-th constraint
            retval = tnlp_->eval_jac_g(nx, xpert, new_x, ng,
                                       nz_jac_g, NULL, NULL, jacpert);
            ASSERT_EXCEPTION(retval, ERROR_IN_TNLP_DERIVATIVE_TEST,
                             "In TNLP derivative test: Jacobian values could not be evaluated at reference point.");
            // ok, now we need to filter the gradient of the icon-th constraint
            IpBlasDcopy(nx, &zero, 0, gradpert, 1);
            IpBlasDcopy(nx, &zero, 0, gradref, 1);
            for (Index i=0; i<nz_jac_g; i++) {
              if (g_iRow[i]==icon) {
                gradpert[g_jCol[i]] += jacpert[i];
                gradref[g_jCol[i]] += jac_g[i];
              }
            }
          }
          new_x = false;

          for (Index ivar2=0; ivar2<nx; ivar2++) {
            Number deriv_approx = (gradpert[ivar2] - gradref[ivar2])/this_perturbation;
            Number deriv_exact = 0.;
            bool found = false;
            for (Index i=0; i<nz_hess_lag; i++) {
              if ( (h_iRow[i]==ivar && h_jCol[i]==ivar2) ||
                   (h_jCol[i]==ivar && h_iRow[i]==ivar2) ) {
                deriv_exact += h_values[i];
                found = true;
              }
            }
            Number rel_error =
              fabs(deriv_approx-deriv_exact)/Max(fabs(deriv_approx),1.);
            char cflag=' ';
            if (rel_error >= derivative_test_tol_) {
              cflag='*';
              nerrors++;
            }
            char sflag = ' ';
            if (found) {
              sflag = 'v';
            }
            if (cflag != ' ' || derivative_test_print_all_) {
              if (icon==-1) {
                jnlst_->Printf(J_WARNING, J_NLP,
                               "%c             obj_hess[%5d,%5d] = %23.16e %c  ~ %23.16e  [%10.3e]\n",
                               cflag, ivar+index_correction, ivar2+index_correction,
                               deriv_exact, sflag, deriv_approx, rel_error);
              }
              else {
                jnlst_->Printf(J_WARNING, J_NLP,
                               "%c %5d-th constr_hess[%5d,%5d] = %23.16e %c  ~ %23.16e  [%10.3e]\n",
                               cflag, icon+index_correction, ivar+index_correction, ivar2+index_correction,
                               deriv_exact, sflag, deriv_approx, rel_error);
              }
            }

          }

          xpert[ivar] = xref[ivar];
        }

        if (icon>=0) {
          lambda[icon] = 0.;
        }
      }

      delete [] h_iRow;
      delete [] h_jCol;
      delete [] h_values;
      delete [] lambda;
      delete [] gradref;
      delete [] gradpert;
      delete [] jacpert;
    }

    delete [] xref;
    delete [] gref;
    delete [] grad_f;
    delete [] xpert;
    delete [] g_iRow;
    delete [] g_jCol;
    delete [] jac_g;
    delete [] gpert;

    if (nerrors==0) {
      jnlst_->Printf(J_SUMMARY, J_NLP,
                     "\nNo errors detected by derivative checker.\n\n");
    }
    else {
      jnlst_->Printf(J_WARNING, J_NLP,
                     "\nDerivative checker detected %d error(s).\n\n", nerrors);
    }

    return retval;
  }

} // namespace Ipopt
