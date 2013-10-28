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
#include "simmath/internal/OptimizerRep.h"

namespace SimTK {

//////////////////////
// DefaultOptimizer //
//////////////////////

Optimizer::OptimizerRep* DefaultOptimizer::clone() const {
    return( new DefaultOptimizer(*this));
}
Real DefaultOptimizer::optimize(  Vector &results ) {
    SimTK_APIARGCHECK_ALWAYS(false,"Optimizer","optimize",
    "the OptimizerSystem has not been set \n");
}
OptimizerAlgorithm DefaultOptimizer::getAlgorithm() const {
    SimTK_APIARGCHECK_ALWAYS(false,"Optimizer","getAlgorithm",
    "the OptimizerSystem has not been set \n");
}

//////////////////
// OptimizerRep //
//////////////////

Optimizer::OptimizerRep::~OptimizerRep() {
    delete jacDiff;
    delete gradDiff;
    delete cf;
    delete of;
}

void Optimizer::OptimizerRep::setConvergenceTolerance(Real accuracy ) {
   convergenceTolerance = accuracy;
}
void Optimizer::OptimizerRep::setConstraintTolerance(Real tolerance ) {
   constraintTolerance = tolerance;
}
void Optimizer::OptimizerRep::setMaxIterations( const int iter ) {
   maxIterations = iter;
}
void Optimizer::OptimizerRep::setLimitedMemoryHistory( const int history ) {
   limitedMemoryHistory = history;
}

void Optimizer::OptimizerRep::setDiagnosticsLevel( const int  level ) {
   diagnosticsLevel = level;
}

// These helper methods are private to this compilation unit.
template<class T> 
static bool setAdvancedOptionHelper(std::map<std::string,T> &optionMap, const std::string &option, const T &value) {
    optionMap[option] = value;
    return true;
}

template<class T> 
static bool getAdvancedOptionHelper(const std::map<std::string,T> &optionMap, const std::string &option, T &value) {
    typename std::map<std::string,T>::const_iterator iter = optionMap.find(option);
    if(iter != optionMap.end()) {
        value = iter->second;
        return true;
    } else
        return false;
}

bool Optimizer::OptimizerRep::setAdvancedStrOption( const std::string &option, const std::string &value ) {
    return setAdvancedOptionHelper(advancedStrOptions, option, value);
}
bool Optimizer::OptimizerRep::setAdvancedRealOption( const std::string &option, const Real value ) {
    return setAdvancedOptionHelper(advancedRealOptions, option, value);
}
bool Optimizer::OptimizerRep::setAdvancedIntOption( const std::string &option, const int value ) {
    return setAdvancedOptionHelper(advancedIntOptions, option, value);
}
bool Optimizer::OptimizerRep::setAdvancedBoolOption( const std::string &option, const bool value ) {
    return setAdvancedOptionHelper(advancedBoolOptions, option, value);
}

bool Optimizer::OptimizerRep::getAdvancedStrOption( const std::string &option, std::string &value ) const {
    return getAdvancedOptionHelper(advancedStrOptions, option, value);
}
bool Optimizer::OptimizerRep::getAdvancedRealOption( const std::string &option, Real &value ) const {
    return getAdvancedOptionHelper(advancedRealOptions, option, value);
}
bool Optimizer::OptimizerRep::getAdvancedIntOption( const std::string &option, int &value ) const {
    return getAdvancedOptionHelper(advancedIntOptions, option, value);
}
bool Optimizer::OptimizerRep::getAdvancedBoolOption( const std::string &option, bool &value ) const {
    return getAdvancedOptionHelper(advancedBoolOptions, option, value);
}

// TODO: this only works if called *prior* to the routines below.
void Optimizer::OptimizerRep::
setDifferentiatorMethod(Differentiator::Method method) {
     diffMethod = method;
}

void Optimizer::OptimizerRep::
useNumericalGradient(bool flag, Real objEstAccuracy) {
    objectiveEstimatedAccuracy = 
        objEstAccuracy > 0 ? objEstAccuracy : SignificantReal;
    delete gradDiff; gradDiff=0; delete of; of=0;
    if (flag) {     // turn on numerical jacbobian
        of = new SysObjectiveFunc(sysp->getNumParameters(), sysp);
        of->setEstimatedAccuracy(objectiveEstimatedAccuracy);
        gradDiff = new Differentiator(*of, diffMethod);
    }
    numericalGradient = flag;
}
void Optimizer::OptimizerRep::
useNumericalJacobian(bool flag, Real consEstAccuracy) {
    constraintsEstimatedAccuracy = 
        consEstAccuracy > 0 ? consEstAccuracy : SignificantReal;
    delete jacDiff; jacDiff=0; delete cf; cf=0;
    if (flag) {     // turn on numerical gradients
        cf = new SysConstraintFunc(sysp->getNumConstraints(), 
                                   sysp->getNumParameters(), sysp);
        cf->setEstimatedAccuracy(constraintsEstimatedAccuracy);
        jacDiff = new Differentiator(*cf, diffMethod); 
    }
    numericalJacobian = flag;
}


int Optimizer::OptimizerRep::objectiveFuncWrapper
   (int n, const Real* x, int newX, Real* f, void* vrep)
{
    assert(vrep);
    const OptimizerRep* rep = reinterpret_cast<const OptimizerRep*>(vrep);

    // This Vector is just a reference to existing space.
    const Vector    params(n, x, true); 
    const bool      isNewParam  = (newX==1);
    Real&           frep        = *f;

    return (rep->getOptimizerSystem().objectiveFunc(params, isNewParam, frep)==0) 
            ? 1 : 0;
}

int Optimizer::OptimizerRep::gradientFuncWrapper
   (int n, const Real* x, int newX, Real* gradient, void* vrep)
{
    assert(vrep);
    const OptimizerRep* rep = reinterpret_cast<const OptimizerRep*>(vrep);

    // These Vectors are just references to existing space.
    const Vector    params(n, x, true); 
    Vector          grad_vec(n,gradient,true);
    Real            fy0;
    const bool      isNewParam = (newX==1);

    const OptimizerSystem&  osys = rep->getOptimizerSystem();

    if( rep->isUsingNumericalGradient() ) {
        osys.objectiveFunc(params, true, fy0);
        rep->getGradientDifferentiator().calcGradient(params, fy0, grad_vec);
        return 1;
    }

    return (osys.gradientFunc(params, isNewParam, grad_vec)==0)
            ? 1 : 0;
}

int Optimizer::OptimizerRep::constraintFuncWrapper
   (int n, const Real* x, int newX, int m, Real* g, void* vrep)
{
    assert(vrep);
    const OptimizerRep* rep = reinterpret_cast<const OptimizerRep*>(vrep);

    // These Vectors are just references to existing space.
    const Vector    params( n, x, true);
    Vector          constraints(m, g, true);
    const bool      isNewParam = (newX==1);

    return (rep->getOptimizerSystem().constraintFunc(params, isNewParam, constraints)==0)
            ? 1 : 0;
}

int Optimizer::OptimizerRep::constraintJacobianWrapper
   (int n, const Real* x, int newX, int m, int nele_jac,
    int* iRow, int* jCol, Real* values, void* vrep)
{
    assert(vrep);
    const OptimizerRep* rep = reinterpret_cast<const OptimizerRep*>(vrep);

    if(m==0) return 1; // m==0 case occurs if you run IPOPT with no constraints

    const bool isNewParam = (newX==1);

    if (values == NULL) {
        // always assume  the jacobian is dense
        int index = 0;
        for(int j=0; j<m; ++j)
            for(int i=0; i<n; ++i) {
                iRow[index]     = j;
                jCol[index++]   = i;
            }
        return 1;   // success
    }

    // Calculate the Jacobian of the constraints.

    const Vector    params(n,x,true);   // This Vector refers to existing space 
    Matrix          jac(m,n);           // This is a new local temporary. TODO: get rid of this

    int status = -1;
    if( rep->isUsingNumericalJacobian() ) {
        Vector sfy0(m);            
        status = rep->getOptimizerSystem().constraintFunc(params, true, sfy0);
        rep->getJacobianDifferentiator().calcJacobian( params, sfy0, jac);
    } else {
        status = rep->getOptimizerSystem().constraintJacobian(params, isNewParam, jac);
    }

    // Transpose the jacobian because Ipopt indexes in Row major format.
    Real *ptr = values;
    for(int j=0; j<m; ++j)
        for(int i=0; i<n; ++i,++ptr)
            *ptr = jac(j,i);

    return (status==0) ? 1 : 0;
}

// TODO finish hessianWrapper
int Optimizer::OptimizerRep::hessianWrapper
   (int n, const Real* x, int newX, Real obj_factor,
    int m, Real* lambda, int new_lambda,
    int nele_hess, int* iRow, int* jCol,
    Real* values, void* vrep)
{
    assert(vrep);
    const OptimizerRep* rep = reinterpret_cast<const OptimizerRep*>(vrep);

    // These Vectors refer to existing space.
    const Vector coeff(n,x,true); 
    Vector hess(n*n,values,true); 
    const bool isNewParam = (newX==1);

    return rep->getOptimizerSystem().hessian(coeff, isNewParam, hess)==0
            ? 1 : 0;
}

} // namespace SimTK
