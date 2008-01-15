
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
#include "simmath/Optimizer.h"
#include "LBFGSOptimizer.h"
#include "LBFGSBOptimizer.h"
#include "InteriorPointOptimizer.h"
#include "CFSQPOptimizer.h"
#include <string>

namespace SimTK {

Optimizer::~Optimizer() {
   delete( (OptimizerRep *)rep );
}

bool Optimizer::isAlgorithmAvailable(OptimizerAlgorithm algorithm) {
    switch(algorithm) {
        case InteriorPoint: return InteriorPointOptimizer::isAvailable();
        case LBFGS: return LBFGSOptimizer::isAvailable();
        case LBFGSB: return LBFGSBOptimizer::isAvailable();
        case CFSQP: return CFSQPOptimizer::isAvailable();
        default: return false;
    }
}
Optimizer::Optimizer( const OptimizerSystem& sys) {
    // Perform construction of the OptimizerRep on the library side.
    librarySideOptimizerConstructor(sys, BestAvailiable );
    // But fill in function pointers from the client side.
    clientSideOptimizerConstructor();
}
Optimizer::Optimizer( const OptimizerSystem& sys, OptimizerAlgorithm algorithm) {
    // Perform construction of the OptimizerRep on the library side.
    librarySideOptimizerConstructor(sys, algorithm);
    // But fill in function pointers from the client side.
    clientSideOptimizerConstructor();
}

Optimizer::Optimizer() {
    rep = (OptimizerRep *) new DefaultOptimizer();
}
void Optimizer::setOptimizerSystem( const OptimizerSystem& sys ) {
    delete rep;
    librarySideOptimizerConstructor( sys, BestAvailiable );
    clientSideOptimizerConstructor();
}
void Optimizer::setOptimizerSystem( const OptimizerSystem& sys, OptimizerAlgorithm algorithm ) {
    delete rep;
    librarySideOptimizerConstructor( sys, algorithm );
    clientSideOptimizerConstructor();
}
void Optimizer::librarySideOptimizerConstructor( const OptimizerSystem& sys, OptimizerAlgorithm algorithm ) {
 
    rep = 0;

    // if constructor specifies which algorithm, use it else select base on problem paramters 
    if ( algorithm == InteriorPoint ) {
         rep = (OptimizerRep *) new InteriorPointOptimizer( sys  );
    } else if( algorithm == LBFGSB ) {
         rep = (OptimizerRep *) new LBFGSBOptimizer( sys  );
    } else if( algorithm == LBFGS ) {
         rep = (OptimizerRep *) new LBFGSOptimizer( sys  );
    } else if( algorithm == CFSQP ) {
         try {
             rep = (OptimizerRep *) new CFSQPOptimizer( sys  );
         } catch (const SimTK::Exception::Base &ex) {
             std::cout << ex.getMessage() << std::endl;
            std::cout << "Failed to load CFSQP library.  Will fall back to default optimizer choice." << std::endl;
            rep = 0;
        }
    }

    if(!rep) { 
        if( sys.getNumConstraints() > 0)   {
            rep = (OptimizerRep *) new InteriorPointOptimizer( sys  );
        } else if( sys.getHasLimits() ) {
            rep = (OptimizerRep *) new LBFGSBOptimizer( sys  );
        } else {
            rep = (OptimizerRep *) new LBFGSOptimizer( sys  );
        }
    } 
    rep->setMyHandle(*this);
    updRep().sysp = &sys;
}

void Optimizer::useNumericalGradient( const bool flag ) {

    ((OptimizerRep *)rep)->useNumericalGradient(flag);
    return;
}
void Optimizer::useNumericalJacobian( const bool flag )  {

    ((OptimizerRep *)rep)->useNumericalJacobian(flag);
    return;
}

void Optimizer::setConvergenceTolerance( const Real tolerance ){

    ((OptimizerRep *)rep)->setConvergenceTolerance(tolerance);
    return;
}

void Optimizer::setMaxIterations( const int iter ){

    ((OptimizerRep *)rep)->setMaxIterations(iter);
    return;
}

void Optimizer::setDifferentiatorMethod( Differentiator::Method method){

    ((OptimizerRep *)rep)->setDifferentiatorMethod(method);
    return;
}

void Optimizer::setLimitedMemoryHistory( const int history ){

    ((OptimizerRep *)rep)->setLimitedMemoryHistory(history);
    return;
}

void Optimizer::setDiagnosticsLevel( const int  level ){

    ((OptimizerRep *)rep)->setDiagnosticsLevel(level);
    return;
}

bool Optimizer::setAdvancedStrOption( const char *option, const char *value ){

    return( ((OptimizerRep *)rep)->setAdvancedStrOption( option, value) );
}

bool Optimizer::setAdvancedRealOption( const char *option, const Real value ){

    return( ((OptimizerRep *)rep)->setAdvancedRealOption( option, value) );
}

bool Optimizer::setAdvancedIntOption( const char *option, const int value ){

    return( ((OptimizerRep *)rep)->setAdvancedIntOption( option, value) );
}

bool Optimizer::setAdvancedBoolOption( const char *option, const bool value ){

    return( ((OptimizerRep *)rep)->setAdvancedBoolOption( option, value) );
}

Real Optimizer::optimize(SimTK::Vector   &results) {
    return( ((OptimizerRep *)rep)->optimize(results) );
}
   
int objectiveFuncWrapper( int n, Real *x, int newX,  Real *f, void* userData) {
    Vector parameters( n, x, true);
    Real& frep = *f;
    bool newParam;

    if( newX == 1 )
        newParam = true;
    else
         newParam = false;

    const OptimizerRep& rep = *reinterpret_cast<const OptimizerRep*>(userData);

    if( 0 == rep.objectiveFunc( rep.getOptimizerSystem(), parameters, newParam, frep )) {
        return(1);
    } else {
        return(0);
    }
}
int gradientFuncWrapper( int n, Real *x, int newX, Real *gradient, void* userData) {

    Real fy0;
    Real& sfy0 = fy0;
    Vector params( n, x, true);
    Vector grad_vec(n,gradient,true);
    const OptimizerRep& rep = *reinterpret_cast<const OptimizerRep*>(userData);
    bool newParam;

    if( newX == 1 )
         newParam = true;
    else
         newParam = false;

      if( rep.getNumericalGradient() ) {
          rep.getOptimizerSystem().objectiveFunc( params, true, sfy0 );
          rep.gradDiff->calcGradient( params, sfy0, grad_vec);
          return(1);
      } else {
          if( 0 == rep.gradientFunc( rep.getOptimizerSystem(), params, newParam, grad_vec )) {
            return(1);
          } else {
            return(0);
          }
      }

}
int constraintFuncWrapper( int n, Real *x, int newX, int m, Real *g,  void*userData) {
      Vector parameters( n, x, true);
      Vector constraints(m, g, true);
      const OptimizerRep& rep = *reinterpret_cast<const OptimizerRep*>(userData);
      bool newParam;

      if( newX == 1 )
          newParam = true;
      else
          newParam = false;
      if( 0 == rep.constraintFunc( rep.getOptimizerSystem(), parameters, newParam, constraints )) {
         return 1;
      } else { 
         return 0;
      }
}
int constraintJacobianWrapper(int n, Real *x, int newX, int m, Index nele_jac,
                int *iRow, int *jCol, Real *values, void *userData)
{
    if(m==0) return 1; // m==0 case occurs if you run IPOPT with no constraints

    int i,j,index;
    bool newParam;
    if( newX == 1 )
        newParam = true;
    else
        newParam = false;

    if (values == NULL) {

         /* always assume  the jacobian is dense */
         index = 0;
         for(j=0;j<m;j++) {
              for(i=0;i<n;i++) {
                   iRow[index] = j;
                   jCol[index++] = i;
               }
          }
    } else {
         /* jacobian of the constraints */
    
          Vector params(n,x,true); 
          const OptimizerRep& rep = *reinterpret_cast<const OptimizerRep*>(userData);

          Matrix jac(m,n);      

          if( rep.getNumericalJacobian() ) {
               Vector sfy0(m);            
               rep.getOptimizerSystem().constraintFunc( params, true, sfy0 );
               rep.jacDiff->calcJacobian( params, sfy0, jac);
          } else {
               rep.constraintJacobian( rep.getOptimizerSystem(), params, newParam, jac );
          }

         /* transpose the jacobian because Ipopt indexes in Row major format */
         Real *ptr = values;
         for(j=0;j<m;j++) {
             for(i=0;i<n;i++,ptr++) {
                 *ptr = jac(j,i);
             }
         }
    } 
    return( 1 );
}

// TODO finish hessianWrapper
int hessianWrapper(int n, Real *x, int newX, Real obj_factor,
    int m, Number *lambda, int new_lambda,
    int nele_hess, int *iRow, int *jCol,
    Real *values, void *userData) {


    Vector coeff(n,x,true); 
    Vector hess(n*n,values,true); 
    const OptimizerRep& rep = *reinterpret_cast<const OptimizerRep*>(userData);
    bool newParam;

    if( newX == 1 )
        newParam = true;
    else
        newParam = false;

    return( rep.hessian( rep.getOptimizerSystem(), coeff, newParam, hess ));
}


void Optimizer::registerObjectiveFunc(Optimizer::ObjectiveFunc f) {
    updRep().objectiveFunc = f;
}
void Optimizer::registerGradientFunc(Optimizer::GradientFunc f) {
    updRep().gradientFunc = f;
}
void Optimizer::registerConstraintFunc(Optimizer::ConstraintFunc f) {
    updRep().constraintFunc = f;
}
void Optimizer::registerConstraintJacobian(Optimizer::ConstraintJacobian f) {
    updRep().constraintJacobian = f;
}
void Optimizer::registerHessian(Optimizer::Hessian f) {
    updRep().hessian = f;
}


} // namespace SimTK
