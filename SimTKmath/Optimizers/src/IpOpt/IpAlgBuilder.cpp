// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpAlgBuilder.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-09-29

#include "IpoptConfig.h"
#include "IpAlgBuilder.hpp"

#include "IpOrigIpoptNLP.hpp"
#include "IpIpoptData.hpp"
#include "IpIpoptCalculatedQuantities.hpp"

#include "IpStdAugSystemSolver.hpp"
#include "IpAugRestoSystemSolver.hpp"
#include "IpPDFullSpaceSolver.hpp"
#include "IpPDPerturbationHandler.hpp"
#include "IpOptErrorConvCheck.hpp"
#include "IpBacktrackingLineSearch.hpp"
#include "IpFilterLSAcceptor.hpp"
#include "IpMonotoneMuUpdate.hpp"
#include "IpAdaptiveMuUpdate.hpp"
#include "IpLoqoMuOracle.hpp"
#include "IpProbingMuOracle.hpp"
#include "IpQualityFunctionMuOracle.hpp"
#include "IpRestoMinC_1Nrm.hpp"
#include "IpLeastSquareMults.hpp"
#include "IpDefaultIterateInitializer.hpp"
#include "IpWarmStartIterateInitializer.hpp"
#include "IpOrigIterationOutput.hpp"
#include "IpLimMemQuasiNewtonUpdater.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpLowRankAugSystemSolver.hpp"
#include "IpRestoIterationOutput.hpp"
#include "IpRestoFilterConvCheck.hpp"
#include "IpRestoIterateInitializer.hpp"
#include "IpRestoRestoPhase.hpp"
#include "IpTSymLinearSolver.hpp"
#include "IpUserScaling.hpp"
#include "IpGradientScaling.hpp"
#include "IpExactHessianUpdater.hpp"

# include "IpLapackSolverInterface.hpp"

namespace SimTKIpopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  void AlgorithmBuilder::BuildIpoptObjects(const Journalist& jnlst,
      const OptionsList& options,
      const std::string& prefix,
      const SmartPtr<NLP>& nlp,
      SmartPtr<IpoptNLP>& ip_nlp,
      SmartPtr<IpoptData>& ip_data,
      SmartPtr<IpoptCalculatedQuantities>& ip_cq)
  {
    DBG_ASSERT(prefix == "");

    SmartPtr<NLPScalingObject> nlp_scaling ;
    std::string nlp_scaling_method;
    options.GetStringValue("nlp_scaling_method", nlp_scaling_method, "");
    if (nlp_scaling_method == "user-scaling") {
      nlp_scaling = new UserScaling(ConstPtr(nlp));
    }
    else if (nlp_scaling_method == "gradient-based") {
      nlp_scaling = new GradientScaling(nlp);
    }
    else {
      nlp_scaling = new NoNLPScalingObject();
    }

    ip_nlp = new OrigIpoptNLP(&jnlst, GetRawPtr(nlp), nlp_scaling);

    // Create the IpoptData
    ip_data = new IpoptData();

    // Create the IpoptCalculators
    ip_cq = new IpoptCalculatedQuantities(ip_nlp, ip_data);
  }

  void AlgorithmBuilder::RegisterOptions(SmartPtr<RegisteredOptions> roptions)
  {
    roptions->SetRegisteringCategory("Linear Solver");
    roptions->AddStringOption6(
      "linear_solver",
      "Linear solver used for step computations.",
      "lapack",
      "ma27", "use the Harwell routine MA27",
      "ma57", "use the Harwell routine MA57",
      "pardiso", "use the Pardiso package",
      "taucs", "use TAUCS package (not yet working)",
      "mumps", "use MUMPS package (not yet working)",
      "lapack", "use LAPACK package",
      "Determines which linear algebra package is to be used for the "
      "solution of the augmented linear system (for obtaining the search "
      "directions). "
      "Note, the code must have been compiled with the linear solver you want "
      "to choose. Depending on your Ipopt installation, not all options are "
      "available.");
    roptions->SetRegisteringCategory("Linear Solver");
    roptions->AddStringOption2(
      "linear_system_scaling",
      "Method for scaling the linear system.",
#ifdef HAVE_MC19
      "mc19",
#else
      "none",
#endif
      "none", "no scaling will be performed",
      "mc19", "use the Harwell routine MC19",
      "Determines the method used to compute symmetric scaling "
      "factors for the augmented system (see also the "
      "\"linear_scaling_on_demand\" option).  This scaling is independent"
      "of the NLP problem scaling.  By default, MC19 is only used if MA27 or "
      "MA57 are selected as linear solvers. This option is only available if "
      "Ipopt has been compiled with MC19.");

    roptions->SetRegisteringCategory("NLP Scaling");
    roptions->AddStringOption3(
      "nlp_scaling_method",
      "Select the technique used for scaling the NLP.",
      "gradient-based",
      "none", "no problem scaling will be performed",
      "user-scaling", "scaling parameters will come from the user",
      "gradient-based", "scale the problem so the maximum gradient at the starting point is scaling_max_gradient",
      "Selects the technique used for scaling the problem internally before it is solved."
      " For user-scaling, the parameters come from the NLP. If you are using "
      "AMPL, they can be specified through suffixes (\"scaling_factor\")");

    roptions->SetRegisteringCategory("Barrier Parameter Update");
    roptions->AddStringOption2(
      "mu_strategy",
      "Update strategy for barrier parameter.",
      "monotone",
      "monotone", "use the monotone (Fiacco-McCormick) strategy",
      "adaptive", "use the adaptive update strategy",
      "Determines which barrier parameter update strategy is to be used.");
    roptions->AddStringOption3(
      "mu_oracle",
      "Oracle for a new barrier parameter in the adaptive strategy.",
      "quality-function",
      "probing", "Mehrotra's probing heuristic",
      "loqo", "LOQO's centrality rule",
      "quality-function", "minimize a quality function",
      "Determines how a new barrier parameter is computed in each "
      "\"free-mode\" iteration of the adaptive barrier parameter "
      "strategy. (Only considered if \"adaptive\" is selected for "
      "option \"mu_strategy\").");
    roptions->AddStringOption4(
      "fixed_mu_oracle",
      "Oracle for the barrier parameter when switching to fixed mode.",
      "average_compl",
      "probing", "Mehrotra's probing heuristic",
      "loqo", "LOQO's centrality rule",
      "quality-function", "minimize a quality function",
      "average_compl", "base on current average complementarity",
      "Determines how the first value of the barrier parameter should be "
      "computed when switching to the \"monotone mode\" in the adaptive "
      "strategy. (Only considered if \"adaptive\" is selected for option "
      "\"mu_strategy\".)");
  }

  SmartPtr<IpoptAlgorithm>
  AlgorithmBuilder::BuildBasicAlgorithm(const Journalist& jnlst,
                                        const OptionsList& options,
                                        const std::string& prefix)
  {
    DBG_START_FUN("AlgorithmBuilder::BuildBasicAlgorithm",
                  dbg_verbosity);
    // Create the convergence check
    SmartPtr<ConvergenceCheck> convCheck =
      new OptimalityErrorConvergenceCheck();

    // Create the solvers that will be used by the main algorithm

    SmartPtr<SparseSymLinearSolverInterface> SolverInterface;
    std::string linear_solver;
    options.GetStringValue("linear_solver", linear_solver, prefix);
    if (linear_solver=="ma27") {

      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear solver MA27 not available.");

    }
    else if (linear_solver=="ma57") {

      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear solver MA57 not available.");

    }
    else if (linear_solver=="pardiso") {
      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear solver Pardiso not available.");

    }
    else if (linear_solver=="taucs") {
      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear solver TAUCS not available.");

    }
    else if (linear_solver=="wsmp") {

      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear solver WSMP not available.");

    }
    else if (linear_solver=="mumps") {

      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear solver MUMPS not available.");

    }
    else if (linear_solver=="lapack") {
      SolverInterface = new LapackSolverInterface();

    }

    SmartPtr<TSymScalingMethod> ScalingMethod;
    std::string linear_system_scaling;
    if (!options.GetStringValue("linear_system_scaling",
                                linear_system_scaling, prefix)) {
      // By default, don't use mc19 for non-HSL solvers
      if (linear_solver!="ma27" && linear_solver!="ma57") {
        linear_system_scaling="none";
      }
    }
    if (linear_system_scaling=="mc19") {

      THROW_EXCEPTION(OPTION_INVALID,
                      "Selected linear system scaling method MC19 not available.");

    }

    SmartPtr<SymLinearSolver> ScaledSolver =
      new TSymLinearSolver(SolverInterface, ScalingMethod);

    SmartPtr<AugSystemSolver> AugSolver =
      new StdAugSystemSolver(*ScaledSolver);
    Index enum_int;
    options.GetEnumValue("hessian_approximation", enum_int, prefix);
    HessianApproximationType hessian_approximation =
      HessianApproximationType(enum_int);
    if (hessian_approximation==LIMITED_MEMORY) {
      SmartPtr<AugSystemSolver> tmp =
        new LowRankAugSystemSolver(*AugSolver);
      AugSolver = tmp;
    }

    SmartPtr<PDPerturbationHandler> pertHandler =
      new PDPerturbationHandler();
    SmartPtr<PDSystemSolver> PDSolver =
      new PDFullSpaceSolver(*AugSolver, *pertHandler);

    // Create the object for initializing the iterates Initialization
    // object.  We include both the warm start and the default
    // initializer, so that the warm start options can be activated
    // without having to rebuild the algorithm
    SmartPtr<EqMultiplierCalculator> EqMultCalculator =
      new LeastSquareMultipliers(*AugSolver);
    SmartPtr<IterateInitializer> WarmStartInitializer =
      new WarmStartIterateInitializer();
    SmartPtr<IterateInitializer> IterInitializer =
      new DefaultIterateInitializer(EqMultCalculator, WarmStartInitializer);

    // Solver for the restoration phase
    SmartPtr<AugSystemSolver> resto_AugSolver =
      new AugRestoSystemSolver(*AugSolver);
    SmartPtr<PDPerturbationHandler> resto_pertHandler =
      new PDPerturbationHandler();
    SmartPtr<PDSystemSolver> resto_PDSolver =
      new PDFullSpaceSolver(*resto_AugSolver, *resto_pertHandler);

    // Convergence check in the restoration phase
    SmartPtr<RestoFilterConvergenceCheck> resto_convCheck =
      new RestoFilterConvergenceCheck();

    // Line search method for the restoration phase
    SmartPtr<RestoRestorationPhase> resto_resto =
      new RestoRestorationPhase();
    SmartPtr<FilterLSAcceptor> resto_filterLSacceptor =
      new FilterLSAcceptor(GetRawPtr(resto_PDSolver));
    SmartPtr<LineSearch> resto_LineSearch =
      new BacktrackingLineSearch(GetRawPtr(resto_filterLSacceptor),
                                 GetRawPtr(resto_resto), GetRawPtr(resto_convCheck));

    // Create the mu update that will be used by the restoration phase
    // algorithm
    SmartPtr<MuUpdate> resto_MuUpdate;
    std::string resto_smuupdate;
    if (!options.GetStringValue("mu_strategy", resto_smuupdate, "resto."+prefix)) {
      // Change default for quasi-Newton option (then we use adaptive)
      Index enum_int;
      if (options.GetEnumValue("hessian_approximation", enum_int, prefix)) {
        HessianApproximationType hessian_approximation =
          HessianApproximationType(enum_int);
        if (hessian_approximation==LIMITED_MEMORY) {
          resto_smuupdate = "adaptive";
        }
      }
    }

    std::string resto_smuoracle;
    std::string resto_sfixmuoracle;
    if (resto_smuupdate=="adaptive" ) {
      options.GetStringValue("mu_oracle", resto_smuoracle, "resto."+prefix);
      options.GetStringValue("fixed_mu_oracle", resto_sfixmuoracle, "resto."+prefix);
    }

    if (resto_smuupdate=="monotone" ) {
      resto_MuUpdate = new MonotoneMuUpdate(GetRawPtr(resto_LineSearch));
    }
    else if (resto_smuupdate=="adaptive") {
      SmartPtr<MuOracle> resto_MuOracle;
      if (resto_smuoracle=="loqo") {
        resto_MuOracle = new LoqoMuOracle();
      }
      else if (resto_smuoracle=="probing") {
        resto_MuOracle = new ProbingMuOracle(resto_PDSolver);
      }
      else if (resto_smuoracle=="quality-function") {
        resto_MuOracle = new QualityFunctionMuOracle(resto_PDSolver);
      }
      SmartPtr<MuOracle> resto_FixMuOracle;
      if (resto_sfixmuoracle=="loqo") {
        resto_FixMuOracle = new LoqoMuOracle();
      }
      else if (resto_sfixmuoracle=="probing") {
        resto_FixMuOracle = new ProbingMuOracle(resto_PDSolver);
      }
      else if (resto_sfixmuoracle=="quality-function") {
        resto_FixMuOracle = new QualityFunctionMuOracle(resto_PDSolver);
      }
      else {
        resto_FixMuOracle = NULL;
      }
      resto_MuUpdate =
        new AdaptiveMuUpdate(GetRawPtr(resto_LineSearch),
                             resto_MuOracle, resto_FixMuOracle);
    }

    // Initialization of the iterates for the restoration phase
    SmartPtr<EqMultiplierCalculator> resto_EqMultCalculator =
      new LeastSquareMultipliers(*resto_AugSolver);
    SmartPtr<IterateInitializer> resto_IterInitializer =
      new RestoIterateInitializer(resto_EqMultCalculator);

    // Create the object for the iteration output during restoration
    SmartPtr<OrigIterationOutput> resto_OrigIterOutput = NULL;
    //   new OrigIterationOutput();
    SmartPtr<IterationOutput> resto_IterOutput =
      new RestoIterationOutput(resto_OrigIterOutput);

    // Get the Hessian updater for the restoration phase
    SmartPtr<HessianUpdater> resto_HessUpdater;
    switch(hessian_approximation) {
      case EXACT:
      resto_HessUpdater = new ExactHessianUpdater();
      break;
      case LIMITED_MEMORY:
      // ToDo This needs to be replaced!
      resto_HessUpdater  = new LimMemQuasiNewtonUpdater(true);
      break;
    }

    // Put together the overall restoration phase IP algorithm
    SmartPtr<IpoptAlgorithm> resto_alg =
      new IpoptAlgorithm(resto_PDSolver,
                         GetRawPtr(resto_LineSearch),
                         GetRawPtr(resto_MuUpdate),
                         GetRawPtr(resto_convCheck),
                         resto_IterInitializer,
                         resto_IterOutput,
                         resto_HessUpdater,
                         resto_EqMultCalculator);

    // Set the restoration phase
    SmartPtr<RestorationPhase> resto_phase =
      new MinC_1NrmRestorationPhase(*resto_alg, EqMultCalculator);

    // Create the line search to be used by the main algorithm
    SmartPtr<FilterLSAcceptor> filterLSacceptor =
      new FilterLSAcceptor(GetRawPtr(PDSolver));
    SmartPtr<LineSearch> lineSearch =
      new BacktrackingLineSearch(GetRawPtr(filterLSacceptor),
                                 GetRawPtr(resto_phase), convCheck);

    // The following cross reference is not good: We have to store a
    // pointer to the lineSearch object in resto_convCheck as a
    // non-SmartPtr to make sure that things are properly deleted when
    // the IpoptAlgorithm return by the Builder is destructed.
    resto_convCheck->SetOrigFilterLSAcceptor(*filterLSacceptor);

    // Create the mu update that will be used by the main algorithm
    SmartPtr<MuUpdate> MuUpdate;
    std::string smuupdate;
    if (!options.GetStringValue("mu_strategy", smuupdate, prefix)) {
      // Change default for quasi-Newton option (then we use adaptive)
      Index enum_int;
      if (options.GetEnumValue("hessian_approximation", enum_int, prefix)) {
        HessianApproximationType hessian_approximation =
          HessianApproximationType(enum_int);
        if (hessian_approximation==LIMITED_MEMORY) {
          smuupdate = "adaptive";
        }
      }
    }
    std::string smuoracle;
    std::string sfixmuoracle;
    if (smuupdate=="adaptive" ) {
      options.GetStringValue("mu_oracle", smuoracle, prefix);
      options.GetStringValue("fixed_mu_oracle", sfixmuoracle, prefix);
    }

    if (smuupdate=="monotone" ) {
      MuUpdate = new MonotoneMuUpdate(GetRawPtr(lineSearch));
    }
    else if (smuupdate=="adaptive") {
      SmartPtr<MuOracle> muOracle;
      if (smuoracle=="loqo") {
        muOracle = new LoqoMuOracle();
      }
      else if (smuoracle=="probing") {
        muOracle = new ProbingMuOracle(PDSolver);
      }
      else if (smuoracle=="quality-function") {
        muOracle = new QualityFunctionMuOracle(PDSolver);
      }
      SmartPtr<MuOracle> FixMuOracle;
      if (sfixmuoracle=="loqo") {
        FixMuOracle = new LoqoMuOracle();
      }
      else if (sfixmuoracle=="probing") {
        FixMuOracle = new ProbingMuOracle(PDSolver);
      }
      else if (sfixmuoracle=="quality-function") {
        FixMuOracle = new QualityFunctionMuOracle(PDSolver);
      }
      else {
        FixMuOracle = NULL;
      }
      MuUpdate = new AdaptiveMuUpdate(GetRawPtr(lineSearch),
                                      muOracle, FixMuOracle);
    }

    // Create the object for the iteration output
    SmartPtr<IterationOutput> IterOutput =
      new OrigIterationOutput();

    // Get the Hessian updater for the main algorithm
    SmartPtr<HessianUpdater> HessUpdater;
    switch(hessian_approximation) {
      case EXACT:
      HessUpdater = new ExactHessianUpdater();
      break;
      case LIMITED_MEMORY:
      // ToDo This needs to be replaced!
      HessUpdater  = new LimMemQuasiNewtonUpdater(false);
      break;
    }

    // Create the main algorithm
    SmartPtr<IpoptAlgorithm> alg =
      new IpoptAlgorithm(PDSolver,
                         GetRawPtr(lineSearch), MuUpdate,
                         convCheck, IterInitializer, IterOutput,
                         HessUpdater, EqMultCalculator);

    return alg;
  }

} // namespace
