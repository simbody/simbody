
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
void Optimizer::OptimizerRep::setDifferentiatorMethod( Differentiator::Method method) {
     diffMethod = method;
}


void Optimizer::OptimizerRep::useNumericalGradient( const bool flag ) {
   if( flag ) {     // turn on numerical gradients
       if( !numericalGradient ) initNumericalGrad();       
       numericalGradient = true;
   } else {        // turn off numerical graidents
       numericalGradient = false;
   }
}
void Optimizer::OptimizerRep::useNumericalJacobian( const bool flag ) {
   if( flag ) {     // turn on numerical jacbobian
       if( !numericalJacobian )  initNumericalJac();       
       numericalJacobian = true;
   } else {
       numericalJacobian = false;
   }
}

void Optimizer::OptimizerRep::initNumericalJac() {  // instaniates a jacobian Differentiator
    cf      = new SysConstraintFunc(sysp->getNumConstraints(), sysp->getNumParameters(), sysp );
    jacDiff = new Differentiator(*cf, diffMethod);  // construct Differentiator
}
void Optimizer::OptimizerRep::initNumericalGrad() {  // instaniates a gradient Differentiator
    of       = new SysObjectiveFunc( sysp->getNumParameters(), sysp );
    gradDiff = new Differentiator(*of, diffMethod);  // construct Differentiator
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
