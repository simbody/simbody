/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-13 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors: Michael Sherman                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "SimTKmath.h"

#include "LBFGSOptimizer.h"
#include "LBFGSBOptimizer.h"
#include "InteriorPointOptimizer.h"
#include "CFSQPOptimizer.h"
#include "CMAESOptimizer.h"
#include <string>

namespace SimTK {

Optimizer::~Optimizer() {
   delete (OptimizerRep *)rep;
}

bool Optimizer::isAlgorithmAvailable(OptimizerAlgorithm algorithm) {
    switch(algorithm) {
        case InteriorPoint: return InteriorPointOptimizer::isAvailable();
        case LBFGS:         return LBFGSOptimizer::isAvailable();
        case LBFGSB:        return LBFGSBOptimizer::isAvailable();
#if SimTK_DEFAULT_PRECISION==2 // double only
        case CFSQP:         return CFSQPOptimizer::isAvailable();
#endif
        case CMAES:         return CMAESOptimizer::isAvailable();
        default:            return false;
    }
}

Optimizer::Optimizer( const OptimizerSystem& sys) : rep(0) {
    rep = constructOptimizerRep(sys, BestAvailable );
}
Optimizer::Optimizer( const OptimizerSystem& sys, OptimizerAlgorithm algorithm) : rep(0) {
    rep = constructOptimizerRep(sys, algorithm);
}

Optimizer::Optimizer() : rep(0) {
    rep = (OptimizerRep *) new DefaultOptimizer();
}

void Optimizer::setOptimizerSystem( const OptimizerSystem& sys ) {
    delete rep;
    rep = constructOptimizerRep( sys, BestAvailable );
}
void Optimizer::setOptimizerSystem( const OptimizerSystem& sys, OptimizerAlgorithm algorithm ) {
    delete rep;
    rep = constructOptimizerRep( sys, algorithm );
}

const OptimizerSystem&
Optimizer::getOptimizerSystem() const {
    assert(rep);
    return rep->getOptimizerSystem();
}

// copy constructor
Optimizer::Optimizer( const Optimizer& c ) : rep(0) {
    if (c.rep)
        rep = c.rep->clone();
}

// copy assignment operator
Optimizer& Optimizer::operator=(const Optimizer& rhs) {
    if (&rhs != this) {
        delete rep; rep = 0;
        if (rhs.rep)
            rep = rhs.rep->clone();
    }
    return *this;
}

Optimizer::OptimizerRep*
Optimizer::constructOptimizerRep( const OptimizerSystem& sys, OptimizerAlgorithm algorithm ) {
    OptimizerRep* newRep = 0;

    // if constructor specifies which algorithm, use it else select based on
    // problem paramters
    if ( algorithm == InteriorPoint ) {
        newRep = (OptimizerRep *) new InteriorPointOptimizer( sys  );
    } else if( algorithm == LBFGSB ) {
        newRep = (OptimizerRep *) new LBFGSBOptimizer( sys  );
    } else if( algorithm == LBFGS ) {
        newRep = (OptimizerRep *) new LBFGSOptimizer( sys  );
    }
#if SimTK_DEFAULT_PRECISION==2 // double only
    else if( algorithm == CFSQP ) {
        try {
            newRep = (OptimizerRep *) new CFSQPOptimizer( sys  );
        } catch (const SimTK::Exception::Base &ex) {
            std::cout << ex.getMessage() << std::endl;
            std::cout << "Failed to load CFSQP library. "
                "Will fall back to default optimizer choice." << std::endl;
            newRep = 0;
        }
    }
#endif
    else if( algorithm == CMAES ) {
        newRep = (OptimizerRep *) new CMAESOptimizer( sys  );
    }

    SimTK_APIARGCHECK_ALWAYS(
            algorithm != UnknownOptimizerAlgorithm &&
            algorithm != UserSuppliedOptimizerAlgorithm,
            "Optimizer", "constructOptimizerRep",
            "UnknownOptimizerAlgorithm and UserSuppliedOptimizerAlgorithm "
            "do not specify specific algorithms.");

    if(!newRep) {
        if( sys.getNumConstraints() > 0)   {
            newRep = (OptimizerRep *) new InteriorPointOptimizer( sys  );
        } else if( sys.getHasLimits() ) {
            newRep = (OptimizerRep *) new LBFGSBOptimizer( sys  );
        } else {
            newRep = (OptimizerRep *) new LBFGSOptimizer( sys  );
        }
    }
    newRep->setMyHandle(*this);
    newRep->sysp = &sys;

    return newRep;
}

void Optimizer::useNumericalGradient(bool flag, Real objEstAccuracy) {
    updRep().useNumericalGradient(flag, objEstAccuracy);
}
void Optimizer::useNumericalJacobian(bool flag, Real consEstAccuracy) {
     updRep().useNumericalJacobian(flag, consEstAccuracy);
}

void Optimizer::setDifferentiatorMethod( Differentiator::Method method) {
     updRep().setDifferentiatorMethod(method);
}

void Optimizer::setConvergenceTolerance( Real accuracy ) {
     updRep().setConvergenceTolerance(accuracy);
}

void Optimizer::setConstraintTolerance( Real tolerance ) {
     updRep().setConstraintTolerance(tolerance);
}

void Optimizer::setMaxIterations( int iter ) {
     updRep().setMaxIterations(iter);
}


void Optimizer::setLimitedMemoryHistory( int history ) {
     updRep().setLimitedMemoryHistory(history);
}

void Optimizer::setDiagnosticsLevel( int  level ) {
     updRep().setDiagnosticsLevel(level);
}

bool Optimizer::setAdvancedStrOption( const char *option, const char *value ) {
    return updRep().setAdvancedStrOption( option, value);
}

bool Optimizer::setAdvancedRealOption( const char *option, const Real value ) {
    return updRep().setAdvancedRealOption( option, value);
}

bool Optimizer::setAdvancedIntOption( const char *option, const int value ) {
    return updRep().setAdvancedIntOption( option, value);
}

bool Optimizer::setAdvancedBoolOption( const char *option, const bool value ) {
    return updRep().setAdvancedBoolOption( option, value);
}

Real Optimizer::optimize(SimTK::Vector   &results) {
    return updRep().optimize(results);
}

bool Optimizer::isUsingNumericalGradient() const {
    return getRep().isUsingNumericalGradient();
}
bool Optimizer::isUsingNumericalJacobian() const {
    return getRep().isUsingNumericalJacobian();
}

Differentiator::Method Optimizer::getDifferentiatorMethod() const {
    return getRep().getDifferentiatorMethod();
}

OptimizerAlgorithm Optimizer::getAlgorithm() const {
    return getRep().getAlgorithm();
}

Real Optimizer::getEstimatedAccuracyOfObjective() const {
    return getRep().getEstimatedAccuracyOfObjective();
}

Real Optimizer::getEstimatedAccuracyOfConstraints() const {
    return getRep().getEstimatedAccuracyOfConstraints();
}

} // namespace SimTK
