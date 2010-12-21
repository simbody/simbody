// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpAlgorithmRegOp.cpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-08-16

#include "IpAlgorithmRegOp.hpp"
#include "IpRegOptions.hpp"

#include "IpAdaptiveMuUpdate.hpp"
#include "IpAlgBuilder.hpp"
#include "IpDefaultIterateInitializer.hpp"
#include "IpBacktrackingLineSearch.hpp"
#include "IpFilterLSAcceptor.hpp"
#include "IpGradientScaling.hpp"
#include "IpIpoptAlg.hpp"
#include "IpIpoptCalculatedQuantities.hpp"
#include "IpIpoptData.hpp"
#include "IpMonotoneMuUpdate.hpp"
#include "IpNLPScaling.hpp"
#include "IpOptErrorConvCheck.hpp"
#include "IpOrigIpoptNLP.hpp"
#include "IpOrigIterationOutput.hpp"
#include "IpLimMemQuasiNewtonUpdater.hpp"
#include "IpPDFullSpaceSolver.hpp"
#include "IpPDPerturbationHandler.hpp"
#include "IpProbingMuOracle.hpp"
#include "IpQualityFunctionMuOracle.hpp"
#include "IpRestoFilterConvCheck.hpp"
#include "IpRestoIpoptNLP.hpp"
#include "IpRestoMinC_1Nrm.hpp"
#include "IpWarmStartIterateInitializer.hpp"


namespace Ipopt
{

  void RegisterOptions_Algorithm(const SmartPtr<RegisteredOptions>& roptions)
  {
    roptions->SetRegisteringCategory("Barrier Parameter Update");
    AdaptiveMuUpdate::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Initialization");
    DefaultIterateInitializer::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Main Algorithm");
    AlgorithmBuilder::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Line Search");
    BacktrackingLineSearch::RegisterOptions(roptions);
    FilterLSAcceptor::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("NLP Scaling");
    StandardScalingBase::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("NLP Scaling");
    GradientScaling::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Uncategorized");
    IpoptAlgorithm::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Uncategorized");
    IpoptData::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Uncategorized");
    IpoptCalculatedQuantities::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Hessian Approximation");
    LimMemQuasiNewtonUpdater::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Barrier Parameter Update");
    MonotoneMuUpdate::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Convergence");
    OptimalityErrorConvergenceCheck::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("NLP");
    OrigIpoptNLP::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Output");
    OrigIterationOutput::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Step Calculation");
    PDFullSpaceSolver::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Step Calculation");
    PDPerturbationHandler::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Barrier Parameter Update");
    ProbingMuOracle::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Barrier Parameter Update");
    QualityFunctionMuOracle::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Restoration Phase");
    RestoFilterConvergenceCheck::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Restoration Phase");
    RestoIpoptNLP::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Uncategorized");
    roptions->SetRegisteringCategory("Restoration Phase");
    MinC_1NrmRestorationPhase::RegisterOptions(roptions);
    roptions->SetRegisteringCategory("Warm Start");
    WarmStartIterateInitializer::RegisterOptions(roptions);
  }

} // namespace Ipopt
