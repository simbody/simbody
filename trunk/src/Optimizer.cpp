
/* Portions copyright (c) 2006 Stanford University and Jack Middleton.
 * Contributors:
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
#include "SimTKmath.h"

#include "LBFGSOptimizer.h"
#include "LBFGSBOptimizer.h"
#include "InteriorPointOptimizer.h"
#include "CFSQPOptimizer.h"
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
        case CFSQP:         return CFSQPOptimizer::isAvailable();
        default:            return false;
    }
}

Optimizer::Optimizer( const OptimizerSystem& sys) : rep(0) {
    rep = constructOptimizerRep(sys, BestAvailiable );
}
Optimizer::Optimizer( const OptimizerSystem& sys, OptimizerAlgorithm algorithm) : rep(0) {
    rep = constructOptimizerRep(sys, algorithm);
}

Optimizer::Optimizer() : rep(0) {
    rep = (OptimizerRep *) new DefaultOptimizer();
}

void Optimizer::setOptimizerSystem( const OptimizerSystem& sys ) {
    delete rep;
    rep = constructOptimizerRep( sys, BestAvailiable );
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

    // if constructor specifies which algorithm, use it else select base on problem paramters 
    if ( algorithm == InteriorPoint ) {
        newRep = (OptimizerRep *) new InteriorPointOptimizer( sys  );
    } else if( algorithm == LBFGSB ) {
        newRep = (OptimizerRep *) new LBFGSBOptimizer( sys  );
    } else if( algorithm == LBFGS ) {
        newRep = (OptimizerRep *) new LBFGSOptimizer( sys  );
    } else if( algorithm == CFSQP ) {
        try {
            newRep = (OptimizerRep *) new CFSQPOptimizer( sys  );
        } catch (const SimTK::Exception::Base &ex) {
            std::cout << ex.getMessage() << std::endl;
            std::cout << "Failed to load CFSQP library.  Will fall back to default optimizer choice." << std::endl;
            newRep = 0;
        }
    }

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

void Optimizer::useNumericalGradient( const bool flag ) {
    updRep().useNumericalGradient(flag);
}
void Optimizer::useNumericalJacobian( const bool flag ) {
     updRep().useNumericalJacobian(flag);
}

void Optimizer::setConvergenceTolerance( const Real tolerance ) {
     updRep().setConvergenceTolerance(tolerance);
}

void Optimizer::setMaxIterations( const int iter ) {
     updRep().setMaxIterations(iter);
}

void Optimizer::setDifferentiatorMethod( Differentiator::Method method) {
     updRep().setDifferentiatorMethod(method);
}

void Optimizer::setLimitedMemoryHistory( const int history ) {
     updRep().setLimitedMemoryHistory(history);
}

void Optimizer::setDiagnosticsLevel( const int  level ) {
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


} // namespace SimTK
