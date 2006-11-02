// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpIpoptApplication.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-02

#include "IpoptConfig.h"
#include "IpIpoptApplication.hpp"
#include "IpTNLPAdapter.hpp"
#include "IpIpoptAlg.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpAlgBuilder.hpp"
#include "IpSolveStatistics.hpp"
#include "IpLinearSolversRegOp.hpp"
#include "IpInterfacesRegOp.hpp"
#include "IpAlgorithmRegOp.hpp"

#ifdef HAVE_CMATH
# include <cmath>
#else
# ifdef HAVE_MATH_H
#  include <math.h>
# else
#  error "don't have header file for math"
# endif
#endif

#include <fstream>

namespace Ipopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  IpoptApplication::IpoptApplication(bool create_console_out)
      :
      jnlst_(new Journalist()),
      options_(new OptionsList()),
      statistics_(NULL),
      alg_(NULL),
      nlp_adapter_(NULL)
  {
    try {
# ifdef IP_DEBUG
      DebugJournalistWrapper::SetJournalist(GetRawPtr(jnlst_));
      SmartPtr<Journal> debug_jrnl = jnlst_->AddFileJournal("Debug", "debug.out", J_ITERSUMMARY);
      debug_jrnl->SetPrintLevel(J_DBG, J_ALL);
# endif

      DBG_START_METH("IpoptApplication::IpoptApplication()",
                     dbg_verbosity);

      if (create_console_out) {
        SmartPtr<Journal> stdout_jrnl =
          jnlst_->AddFileJournal("console", "stdout", J_ITERSUMMARY);
        stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);
      }

      // Register the valid options
      reg_options_ = new RegisteredOptions();
      RegisterAllOptions(reg_options_);

      options_->SetJournalist(jnlst_);
      options_->SetRegisteredOptions(reg_options_);
    }
    catch(IpoptException& exc) {
      exc.ReportException(*jnlst_);
      THROW_EXCEPTION(IPOPT_APPLICATION_ERROR,
                      "Caught unknown Ipopt exception");
    }
    catch(std::bad_alloc) {
      jnlst_->Printf(J_ERROR, J_MAIN, "\nEXIT: Not enough memory.\n");
      THROW_EXCEPTION(IPOPT_APPLICATION_ERROR,
                      "Not enough memory");
    }
    catch(...) {
      IpoptException exc("Unknown Exception caught in ipopt", "Unknown File", -1);
      exc.ReportException(*jnlst_);
      THROW_EXCEPTION(IPOPT_APPLICATION_ERROR,
                      "Caught unknown exception");
    }
  }

  void IpoptApplication::Initialize(std::string params_file /*= "ipopt.opt"*/)
  {
    std::ifstream is;
    if (params_file != "") {
      try {
        is.open(params_file.c_str());
      }
      catch(std::bad_alloc) {
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Not enough memory.\n");
        exit(-1);
      }
      catch(...) {
        IpoptException exc("Unknown Exception caught in ipopt", "Unknown File", -1);
        exc.ReportException(*jnlst_);
        exit(-1);
      }
    }
    Initialize(is);
    if (is) {
      is.close();
    }
  }

  void IpoptApplication::Initialize(std::istream& is)
  {
    try {
      // Get the options
      if (is.good()) {
        // stream exists, read the content
        options_->ReadFromStream(*jnlst_, is);
      }

      Index ivalue;
      options_->GetIntegerValue("print_level", ivalue, "");
      EJournalLevel print_level = (EJournalLevel)ivalue;
      SmartPtr<Journal> stdout_jrnl = jnlst_->GetJournal("console");
      if (IsValid(stdout_jrnl)) {
        // Set printlevel for stdout
        stdout_jrnl->SetAllPrintLevels(print_level);
        stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);
      }

      bool option_set;

#ifdef IP_DEBUG
      // Set printlevel for debug
      option_set = options_->GetIntegerValue("debug_print_level",
                                             ivalue, "");
      EJournalLevel debug_print_level;
      if (option_set) {
        debug_print_level = (EJournalLevel)ivalue;
      }
      else {
        debug_print_level = print_level;
      }
      SmartPtr<Journal> debug_jrnl = jnlst_->GetJournal("Debug");
      debug_jrnl->SetAllPrintLevels(debug_print_level);
      debug_jrnl->SetPrintLevel(J_DBG, J_ALL);
#endif

      // Open an output file if required
      std::string output_filename;
      options_->GetStringValue("output_file", output_filename, "");
      if (output_filename != "") {
        EJournalLevel file_print_level;
        option_set = options_->GetIntegerValue("file_print_level", ivalue, "");
        if (option_set) {
          file_print_level = (EJournalLevel)ivalue;
        }
        else {
          file_print_level = print_level;
        }
        OpenOutputFile(output_filename, file_print_level);
      }

      // output a description of all the options
      bool print_options_documentation;
      options_->GetBoolValue("print_options_documentation",
                             print_options_documentation, "");
      if (print_options_documentation) {
        bool latex;
        options_->GetBoolValue("print_options_latex_mode", latex, "");
        if (latex) {
          std::list<std::string> options_to_print;
          // Output
          options_to_print.push_back("print_level");
          options_to_print.push_back("print_user_options");
          options_to_print.push_back("print_options_documentation");
          options_to_print.push_back("output_file");
          options_to_print.push_back("file_print_level");
          // Termination
          options_to_print.push_back("tol");
          options_to_print.push_back("max_iter");
          options_to_print.push_back("compl_inf_tol");
          options_to_print.push_back("dual_inf_tol");
          options_to_print.push_back("constr_viol_tol");
          options_to_print.push_back("acceptable_tol");
          options_to_print.push_back("acceptable_compl_inf_tol");
          options_to_print.push_back("acceptable_constr_viol_tol");
          options_to_print.push_back("acceptable_dual_inf_tol");
          options_to_print.push_back("diverging_iterates_tol");
          options_to_print.push_back("barrier_tol_factor");
          // NLP scaling
          options_to_print.push_back("obj_scaling_factor");
          options_to_print.push_back("nlp_scaling_method");
          options_to_print.push_back("nlp_scaling_max_gradient");
          // NLP corrections
          options_to_print.push_back("bound_relax_factor");
          options_to_print.push_back("honor_original_bounds");
          options_to_print.push_back("check_derivatives_for_naninf");
          // Barrier parameter
          options_to_print.push_back("mu_strategy");
          options_to_print.push_back("mu_oracle");
          options_to_print.push_back("quality_function_max_section_steps");
          options_to_print.push_back("fixed_mu_oracle");
          options_to_print.push_back("mu_init");
          options_to_print.push_back("mu_max_fact");
          options_to_print.push_back("mu_max");
          options_to_print.push_back("mu_min");
          options_to_print.push_back("mu_linear_decrease_factor");
          options_to_print.push_back("mu_superlinear_decrease_power");
          //options_to_print.push_back("corrector_type");
          // Initialization
          options_to_print.push_back("bound_frac");
          options_to_print.push_back("bound_push");
          options_to_print.push_back("bound_mult_init_val");
          options_to_print.push_back("constr_mult_init_max");
          options_to_print.push_back("bound_mult_init_val");
          // Warm start
          options_to_print.push_back("warm_start_init_point");
          options_to_print.push_back("warm_start_bound_push");
          options_to_print.push_back("warm_start_bound_frac");
          options_to_print.push_back("warm_start_mult_bound_push");
          options_to_print.push_back("warm_start_mult_init_max");
          // Multiplier updates
          options_to_print.push_back("alpha_for_y");
          options_to_print.push_back("recalc_y");
          options_to_print.push_back("recalc_y_feas_tol");
          // Line search
          options_to_print.push_back("max_soc");
          options_to_print.push_back("watchdog_shortened_iter_trigger");
          options_to_print.push_back("watchdog_trial_iter_max");
          // Restoration phase
          options_to_print.push_back("expect_infeasible_problem");
          options_to_print.push_back("expect_infeasible_problem_ctol");
          options_to_print.push_back("start_with_resto");
          options_to_print.push_back("soft_resto_pderror_reduction_factor");
          options_to_print.push_back("required_infeasibility_reduction");
          options_to_print.push_back("bound_mult_reset_threshold");
          options_to_print.push_back("constr_mult_reset_threshold");
          options_to_print.push_back("evaluate_orig_obj_at_resto_trial");
          // Linear solver
          options_to_print.push_back("linear_solver");
          options_to_print.push_back("linear_system_scaling");
          options_to_print.push_back("linear_scaling_on_demand");
          options_to_print.push_back("max_refinement_steps");
          options_to_print.push_back("min_refinement_steps");
          // Hessian perturbation
          options_to_print.push_back("max_hessian_perturbation");
          options_to_print.push_back("min_hessian_perturbation");
          options_to_print.push_back("first_hessian_perturbation");
          options_to_print.push_back("perturb_inc_fact_first");
          options_to_print.push_back("perturb_inc_fact");
          options_to_print.push_back("perturb_dec_fact");
          options_to_print.push_back("jacobian_regularization_value");
          // Quasi-Newton
          options_to_print.push_back("hessian_approximation");
          options_to_print.push_back("limited_memory_max_history");
          options_to_print.push_back("limited_memory_max_skipping");
          // Derivative test
          options_to_print.push_back("derivative_test");
          options_to_print.push_back("derivative_test_perturbation");
          options_to_print.push_back("derivative_test_tol");
          options_to_print.push_back("derivative_test_print_all");

          // Special linear solver
#ifdef HAVE_MA27

          options_to_print.push_back("ma27_pivtol");
          options_to_print.push_back("ma27_pivtolmax");
          options_to_print.push_back("ma27_liw_init_factor");
          options_to_print.push_back("ma27_la_init_factor");
          options_to_print.push_back("ma27_meminc_factor");
#endif

#ifdef HAVE_MA57

          options_to_print.push_back("ma57_pivtol");
          options_to_print.push_back("ma57_pivtolmax");
          options_to_print.push_back("ma57_pre_alloc");
#endif

#ifdef HAVE_PARDISO

          options_to_print.push_back("pardiso_matching_strategy");
          options_to_print.push_back("pardiso_out_of_core_power");
#endif

#ifdef HAVE_WSMP

          options_to_print.push_back("wsmp_num_threads");
          options_to_print.push_back("wsmp_ordering_option");
          options_to_print.push_back("wsmp_pivtol");
          options_to_print.push_back("wsmp_pivtolmax");
          options_to_print.push_back("wsmp_scaling");
#endif

          reg_options_->OutputLatexOptionDocumentation(*jnlst_, options_to_print);
        }
        else {
          std::list<std::string> categories;
          categories.push_back("Output");
          /*categories.push_back("Main Algorithm");*/
          categories.push_back("Convergence");
          categories.push_back("NLP Scaling");
          categories.push_back("Barrier Parameter Update");
          categories.push_back("Line Search");
          categories.push_back("Initialization");
          categories.push_back("Warm Start");
          categories.push_back("Linear Solver");
          categories.push_back("Step Calculation");
          categories.push_back("Restoration Phase");
          categories.push_back("NLP");
          categories.push_back("Derivative Checker");
          categories.push_back("Hessian Approximation");
#ifdef HAVE_MA27

          categories.push_back("MA27 Linear Solver");
#endif
#ifdef HAVE_MA57

          categories.push_back("MA57 Linear Solver");
#endif
#ifdef HAVE_PARDISO

          categories.push_back("Pardiso Linear Solver");
#endif
#ifdef HAVE_WSMP

          categories.push_back("WSMP Linear Solver");
#endif

          categories.push_back("Uncategorized");
          //categories.push_back("Undocumented Options");
          reg_options_->OutputOptionDocumentation(*jnlst_, categories);
        }
      }

    }
    catch(OPTION_INVALID& exc) {
      exc.ReportException(*jnlst_, J_ERROR);
      exit(-1);
    }
    catch(IpoptException& exc) {
      exc.ReportException(*jnlst_, J_ERROR);
      exit(-1);
    }
    catch(std::bad_alloc) {
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Not enough memory.\n");
      exit(-1);
    }
    catch(...) {
      IpoptException exc("Unknown Exception caught in ipopt", "Unknown File", -1);
      exc.ReportException(*jnlst_);
      exit(-1);
    }
  }

  IpoptApplication::~IpoptApplication()
  {
    DBG_START_METH("IpoptApplication::~IpoptApplication()",
                   dbg_verbosity);
  }

  void IpoptApplication::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Output");
    roptions->AddBoundedIntegerOption(
      "print_level",
      "Output verbosity level.",
      0, J_LAST_LEVEL-1, J_ITERSUMMARY,
      "Sets the default verbosity level for console output. The "
      "larger this value the more detailed is the output.");

    roptions->AddStringOption1(
      "output_file",
      "File name of desired output file (leave unset for no file output).",
      "",
      "*", "Any acceptable standard file name",
      "NOTE: This option only works when read from the ipopt.opt options file! "
      "An output file with this name will be written (leave unset for no "
      "file output).  The verbosity level is by default set to \"print_level\", "
      "but can be overridden with \"file_print_level\".  The file name is "
      "changed to use only small letters.");
    roptions->AddBoundedIntegerOption(
      "file_print_level",
      "Verbosity level for output file.",
      0, J_LAST_LEVEL-1, J_ITERSUMMARY,
      "NOTE: This option only works when read from the ipopt.opt options file! "
      "Determines the verbosity level for the file specified by "
      "\"output_file\".  By default it is the same as \"print_level\".");
    roptions->AddStringOption2(
      "print_user_options",
      "Print all options set by the user.",
      "no",
      "no", "don't print options",
      "yes", "print options",
      "If selected, the algorithm will print the list of all options set by "
      "the user including their values and whether they have been used.");
    roptions->AddStringOption2(
      "print_options_documentation",
      "Switch to print all algorithmic options.",
      "no",
      "no", "don't print list",
      "yes", "print list",
      "If selected, the algorithm will print the list of all available "
      "algorithmic options with some documentation before solving the "
      "optimization problem.");

#ifdef IP_DEBUG

    roptions->AddBoundedIntegerOption(
      "debug_print_level",
      "Verbosity level for debug file.",
      0, J_LAST_LEVEL-1, J_ITERSUMMARY,
      "This Ipopt library has been compiled in debug mode, and a file "
      "\"debug.out\" is produced for every run.  This option determines "
      "the verbosity level for this file.  By default it is the same as "
      "\"print_level\".");

#endif

    roptions->AddStringOption2(
      "print_timing_statistics",
      "Switch to print timing statistics.",
      "no",
      "no", "don't print statistics",
      "yes", "print all timing statistics",
      "If selected, the program will print the CPU usage (user time) for "
      "selected tasks.");

    roptions->SetRegisteringCategory("Undocumented");
    roptions->AddStringOption2(
      "print_options_latex_mode",
      "Undocumented", "no",
      "no", "Undocumented",
      "yes", "Undocumented",
      "Undocumented"
    );
  }

  ApplicationReturnStatus
  IpoptApplication::OptimizeTNLP(const SmartPtr<TNLP>& tnlp)
  {
    nlp_adapter_ = new TNLPAdapter(GetRawPtr(tnlp), ConstPtr(jnlst_));

    return OptimizeNLP(nlp_adapter_);
  }

  ApplicationReturnStatus
  IpoptApplication::ReOptimizeTNLP(const SmartPtr<TNLP>& tnlp)
  {
    ASSERT_EXCEPTION(IsValid(nlp_adapter_), INVALID_WARMSTART,
                     "ReOptimizeTNLP called before OptimizeTNLP.");
    TNLPAdapter* adapter =
      dynamic_cast<TNLPAdapter*> (GetRawPtr(nlp_adapter_));
    DBG_ASSERT(adapter);
    ASSERT_EXCEPTION(adapter->tnlp()==tnlp, INVALID_WARMSTART,
                     "ReOptimizeTNLP called for different TNLP.")

    return ReOptimizeNLP(nlp_adapter_);
  }

  ApplicationReturnStatus
  IpoptApplication::OptimizeNLP(const SmartPtr<NLP>& nlp, SmartPtr<AlgorithmBuilder> alg_builder)
  {
    ApplicationReturnStatus retValue = Internal_Error;

    // Prepare internal data structures of the algorithm
    try {

      if (IsNull(alg_builder)) {
        alg_builder = new AlgorithmBuilder();
      }

      alg_builder->BuildIpoptObjects(*jnlst_, *options_, "", nlp,
                                     ip_nlp_, ip_data_, ip_cq_);

      alg_ = GetRawPtr(alg_builder->BuildBasicAlgorithm(*jnlst_, *options_, ""));

      // finally call the optimization
      retValue = call_optimize();
    }
    catch(OPTION_INVALID& exc) {
      exc.ReportException(*jnlst_, J_ERROR);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Invalid option encountered.\n");
      retValue = Invalid_Option;
    }
    catch(IpoptException& exc) {
      exc.ReportException(*jnlst_, J_ERROR);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Some uncaught Ipopt exception encountered.\n");
      retValue = Unrecoverable_Exception;
    }
    catch(std::bad_alloc) {
      retValue = Insufficient_Memory;
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Not enough memory.\n");
    }
    catch(...) {
      IpoptException exc("Unknown Exception caught in Ipopt", "Unknown File", -1);
      exc.ReportException(*jnlst_, J_ERROR);
      retValue = NonIpopt_Exception_Thrown;
    }

    jnlst_->FlushBuffer();

    return retValue;
  }

  ApplicationReturnStatus
  IpoptApplication::ReOptimizeNLP(const SmartPtr<NLP>& nlp)
  {
    ASSERT_EXCEPTION(IsValid(alg_), INVALID_WARMSTART,
                     "ReOptimizeNLP called before OptimizeNLP.");
    OrigIpoptNLP* orig_nlp =
      dynamic_cast<OrigIpoptNLP*> (GetRawPtr(ip_nlp_));
    DBG_ASSERT(orig_nlp);
    ASSERT_EXCEPTION(orig_nlp->nlp()==nlp, INVALID_WARMSTART,
                     "ReOptimizeTNLP called for different NLP.")

    return call_optimize();
  }


  ApplicationReturnStatus IpoptApplication::call_optimize()
  {
    // Reset the print-level for the screen output
    Index ivalue;
    options_->GetIntegerValue("print_level", ivalue, "");
    EJournalLevel print_level = (EJournalLevel)ivalue;
    SmartPtr<Journal> stdout_jrnl = jnlst_->GetJournal("console");
    if (IsValid(stdout_jrnl)) {
      // Set printlevel for stdout
      stdout_jrnl->SetAllPrintLevels(print_level);
      stdout_jrnl->SetPrintLevel(J_DBG, J_NONE);
    }

    statistics_ = NULL; /* delete old statistics */
    // Reset Timing statistics
    ip_data_->TimingStats().ResetTimes();

    ApplicationReturnStatus retValue = Internal_Error;
    try {
      // Get the pointers to the real objects (need to do it that
      // awkwardly since otherwise we would have to include so many
      // things in IpoptApplication, which a user would have to
      // include, too (SmartPtr doesn't work otherwise on AIX)
      IpoptAlgorithm* p2alg = dynamic_cast<IpoptAlgorithm*> (GetRawPtr(alg_));
      IpoptData* p2ip_data = dynamic_cast<IpoptData*> (GetRawPtr(ip_data_));
      OrigIpoptNLP* p2ip_nlp = dynamic_cast<OrigIpoptNLP*> (GetRawPtr(ip_nlp_));
      IpoptCalculatedQuantities* p2ip_cq = dynamic_cast<IpoptCalculatedQuantities*> (GetRawPtr(ip_cq_));
      // Set up the algorithm
      p2alg->Initialize(*jnlst_, *p2ip_nlp, *p2ip_data, *p2ip_cq,
                        *options_, "");

      // Process the options used below
      bool print_timing_statistics;
      options_->GetBoolValue("print_timing_statistics",
                             print_timing_statistics, "");

      // If selected, print the user options
      bool print_user_options;
      options_->GetBoolValue("print_user_options", print_user_options, "");
      if (print_user_options) {
        std::string liststr;
        options_->PrintUserOptions(liststr);
        jnlst_->Printf(J_ERROR, J_MAIN, "\nList of user-set options:\n\n%s", liststr.c_str());
      }

      if( jnlst_->ProduceOutput(J_DETAILED, J_MAIN) ) {
        // Print out the options (including the number of times they were used
        std::string liststr;
        options_->PrintList(liststr);
        jnlst_->Printf(J_DETAILED, J_MAIN, "\nList of options:\n\n%s", liststr.c_str());
      }

      // Run the algorithm
      SolverReturn status = p2alg->Optimize();

      // Since all the output below doesn't make any sense in this
      // case, we rethrow the TOO_FEW_DOF exception here
      ASSERT_EXCEPTION(status != TOO_FEW_DEGREES_OF_FREEDOM, TOO_FEW_DOF,
                       "Too few degrees of freedom (rethrown)!");

      jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                     "\nNumber of Iterations....: %d\n",
                     p2ip_data->iter_count());

      if (status!=INVALID_NUMBER_DETECTED) {
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "\n                                   (scaled)                 (unscaled)\n");
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "Objective...............: %24.16e  %24.16e\n",
                       p2ip_cq->curr_f(),
                       p2ip_cq->unscaled_curr_f());
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "Dual infeasibility......: %24.16e  %24.16e\n",
                       p2ip_cq->curr_dual_infeasibility(NORM_MAX),
                       p2ip_cq->unscaled_curr_dual_infeasibility(NORM_MAX));
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "Constraint violation....: %24.16e  %24.16e\n",
                       p2ip_cq->curr_nlp_constraint_violation(NORM_MAX),
                       p2ip_cq->unscaled_curr_nlp_constraint_violation(NORM_MAX));
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "Complementarity.........: %24.16e  %24.16e\n",
                       p2ip_cq->curr_complementarity(0., NORM_MAX),
                       p2ip_cq->unscaled_curr_complementarity(0., NORM_MAX));
        jnlst_->Printf(J_SUMMARY, J_SOLUTION,
                       "Overall NLP error.......: %24.16e  %24.16e\n\n",
                       p2ip_cq->curr_nlp_error(),
                       p2ip_cq->unscaled_curr_nlp_error());
      }

      p2ip_data->curr()->x()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "x");
      p2ip_data->curr()->y_c()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "y_c");
      p2ip_data->curr()->y_d()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "y_d");
      p2ip_data->curr()->z_L()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "z_L");
      p2ip_data->curr()->z_U()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "z_U");
      p2ip_data->curr()->v_L()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "v_L");
      p2ip_data->curr()->v_U()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "v_U");

      if (status==LOCAL_INFEASIBILITY) {
        p2ip_cq->curr_c()->Print(*jnlst_, J_VECTOR, J_SOLUTION, "curr_c");
        p2ip_cq->curr_d_minus_s()->Print(*jnlst_, J_VECTOR, J_SOLUTION,
                                         "curr_d_minus_s");
      }

      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "\nNumber of objective function evaluations             = %d\n",
                     p2ip_nlp->f_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of objective gradient evaluations             = %d\n",
                     p2ip_nlp->grad_f_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of equality constraint evaluations            = %d\n",
                     p2ip_nlp->c_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of inequality constraint evaluations          = %d\n",
                     p2ip_nlp->d_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of equality constraint Jacobian evaluations   = %d\n",
                     p2ip_nlp->jac_c_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of inequality constraint Jacobian evaluations = %d\n",
                     p2ip_nlp->jac_d_evals());
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Number of Lagrangian Hessian evaluations             = %d\n",
                     p2ip_nlp->h_evals());
      Number time_overall_alg = p2ip_data->TimingStats().OverallAlgorithm().TotalTime();
      Number time_funcs = p2ip_nlp->TotalFunctionEvaluationCPUTime();
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Total CPU secs in IPOPT (w/o function evaluations)   = %10.3f\n",
                     time_overall_alg-time_funcs);
      jnlst_->Printf(J_SUMMARY, J_STATISTICS,
                     "Total CPU secs in NLP function evaluations           = %10.3f\n",
                     time_funcs);

      // Write timing statistics information
      if (print_timing_statistics) {
        jnlst_->Printf(J_SUMMARY, J_TIMING_STATISTICS,
                       "\n\nTiming Statistics:\n\n");
        p2ip_data->TimingStats().PrintAllTimingStatistics(*jnlst_, J_SUMMARY,
            J_TIMING_STATISTICS);
        p2ip_nlp->PrintTimingStatistics(*jnlst_, J_SUMMARY,
                                        J_TIMING_STATISTICS);
      }

      // Write EXIT message
      if (status == SUCCESS) {
        retValue = Solve_Succeeded;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Optimal Solution Found.\n");
      }
      else if (status == MAXITER_EXCEEDED) {
        retValue = Maximum_Iterations_Exceeded;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Maximum Number of Iterations Exceeded.\n");
      }
      else if (status == STOP_AT_TINY_STEP) {
        retValue = Search_Direction_Becomes_Too_Small;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Search Direction is becoming Too Small.\n");
      }
      else if (status == STOP_AT_ACCEPTABLE_POINT) {
        retValue = Solved_To_Acceptable_Level;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Solved To Acceptable Level.\n");
      }
      else if (status == DIVERGING_ITERATES) {
        retValue = Diverging_Iterates;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Iterates divering; problem might be unbounded.\n");
      }
      else if (status == RESTORATION_FAILURE) {
        retValue = Restoration_Failed;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Restoration Failed!\n");
      }
      else if (status == ERROR_IN_STEP_COMPUTATION) {
        retValue = Error_In_Step_Computation;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Error in step computation (regularization becomes too large?)!\n");
      }
      else if (status == LOCAL_INFEASIBILITY) {
        retValue = Infeasible_Problem_Detected;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Converged to a point of local infeasibility. Problem may be infeasible.\n");
      }
      else if (status == USER_REQUESTED_STOP) {
        retValue = User_Requested_Stop;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Stopping optimization at current point as requested by user.\n");
      }
      else if (status == INVALID_NUMBER_DETECTED) {
        retValue = Invalid_Number_Detected;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Invalid number in NLP function or derivative detected.\n");
      }
      else {
        retValue = Internal_Error;
        jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: INTERNAL ERROR: Unknown SolverReturn value - Notify IPOPT Authors.\n");
        return retValue;
      }
      p2ip_nlp->FinalizeSolution(status,
                                 *p2ip_data->curr()->x(),
                                 *p2ip_data->curr()->z_L(),
                                 *p2ip_data->curr()->z_U(),
                                 *p2ip_cq->curr_c(), *p2ip_cq->curr_d(),
                                 *p2ip_data->curr()->y_c(),
                                 *p2ip_data->curr()->y_d(),
                                 p2ip_cq->curr_f());

      if (status!=INVALID_NUMBER_DETECTED) {
        // Create a SolveStatistics object
        statistics_ = new SolveStatistics(p2ip_nlp, p2ip_data, p2ip_cq);
      }
    }
    catch(TOO_FEW_DOF& exc) {
      exc.ReportException(*jnlst_, J_MOREDETAILED);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Problem has too few degrees of freedom.\n");
      retValue = Not_Enough_Degrees_Of_Freedom;
    }
    catch(OPTION_INVALID& exc) {
      exc.ReportException(*jnlst_, J_MOREDETAILED);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Invalid option encountered.\n");
      retValue = Invalid_Option;
    }
    catch(IpoptException& exc) {
      exc.ReportException(*jnlst_, J_ERROR);
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Some uncaught Ipopt exception encountered.\n");
      retValue = Unrecoverable_Exception;
    }
    catch(std::bad_alloc) {
      retValue = Insufficient_Memory;
      jnlst_->Printf(J_SUMMARY, J_MAIN, "\nEXIT: Not enough memory.\n");
    }
    catch(...) {
      IpoptException exc("Unknown Exception caught in Ipopt", "Unknown File", -1);
      exc.ReportException(*jnlst_, J_MOREDETAILED);
      retValue = NonIpopt_Exception_Thrown;
    }

    jnlst_->FlushBuffer();

    return retValue;
  }

  bool IpoptApplication::OpenOutputFile(std::string file_name,
                                        EJournalLevel print_level)
  {
    SmartPtr<Journal> file_jrnl = jnlst_->GetJournal("OutputFile:"+file_name);

    if (IsNull(file_jrnl)) {
      file_jrnl = jnlst_->AddFileJournal("OutputFile:"+file_name,
                                         file_name.c_str(),
                                         print_level);
    }

    file_jrnl->SetPrintLevel(J_DBG, J_NONE);

    return true;
  }

  void IpoptApplication::RegisterAllOptions(const SmartPtr<RegisteredOptions>& roptions)
  {
    RegisterOptions_Interfaces(roptions);
    RegisterOptions_Algorithm(roptions);
    RegisterOptions_LinearSolvers(roptions);
  }

  SmartPtr<SolveStatistics> IpoptApplication::Statistics()
  {
    return statistics_;
  }

} // namespace Ipopt
