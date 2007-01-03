
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
#include "Simmath.h"
#include "Optimizer.h"
#include "LBFGSOptimizer.h"
#include "LBFGSBOptimizer.h"
#include "InteriorPointOptimizer.h"

namespace SimTK {

   Optimizer::~Optimizer() {
         delete( (OptimizerRep *)rep );
      }

void Optimizer::librarySideOptimizerConstructor( OptimizerSystem& sys ) {


      if( sys.getNumConstraints() > 0)   {
         rep = (OptimizerRep *) new InteriorPointOptimizer( sys  );
      }else if( sys.getHasLimits() ) {
         rep = (OptimizerRep *) new LBFGSBOptimizer( sys  );
      } else {
         rep = (OptimizerRep *) new LBFGSOptimizer( sys  );
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

void Optimizer::setDiagnosticsLevel( const int  level ){

      ((OptimizerRep *)rep)->setDiagnosticsLevel(level);
      return;
}

int Optimizer::setAdvancedOptions( const char *option, const Real *values ){

      return( ((OptimizerRep *)rep)->setAdvancedOptions( option, values) );
}

Real Optimizer::optimize(SimTK::Vector   &results) {
      return( ((OptimizerRep *)rep)->optimize(results) );
}
   
int objectiveFuncWrapper( int n, Real *x, int new_x,  Real *f, void* user_data) {
      Vector parameters( n, x, true);
      Real& frep = *f;
      const OptimizerRep& rep = *reinterpret_cast<const OptimizerRep*>(user_data);
      rep.objectiveFunc( rep.getOptimizerSystem(), parameters, new_x, frep );
      return( true );
}
int gradientFuncWrapper( int n, Real *x, int new_x, Real *gradient, void* user_data) {

      Real fy0;
      Real& sfy0 = fy0;
      Vector params( n, x, true);
      Vector grad_vec(n,gradient,true);
      const OptimizerRep& rep = *reinterpret_cast<const OptimizerRep*>(user_data);

      if( rep.getNumericalGradient() ) {
          rep.getOptimizerSystem().objectiveFunc( params, true, sfy0 );
          rep.gradDiff->calcGradient( params, sfy0, grad_vec);
      } else {
          rep.gradientFunc( rep.getOptimizerSystem(), params, new_x, grad_vec );
      }
//printf("gradf = ");for(int i=0;i<n;i++){printf("%f ",gradient[i]);} printf("\n");

      return( true );
}
int constraintFuncWrapper( int n, Real *x, int new_x, int m, Real *g,  void*user_data) {
      Vector parameters( n, x, true);
      Vector constraints(m, g, true);
      const OptimizerRep& rep = *reinterpret_cast<const OptimizerRep*>(user_data);
      rep.constraintFunc( rep.getOptimizerSystem(), parameters, new_x, constraints );
      return( true );
}
int constraintJacobianWrapper(int n, Real *x, int new_x, int m, Index nele_jac,
                int *iRow, int *jCol, Real *values, void *user_data)
{
  int i,j,index;
  Real *nx;
  if (values == NULL) {

    /* always assume  the jacobian is dense */
    index = 0;
    for(j=0;j<m;j++) {
      for(i=0;i<n;i++) {
          iRow[index] = j;
          jCol[index++] = i;
//printf("IROW=%d JCol=%d \n",iRow[index-1],jCol[index-1]);
       }
    }
  } else {
    /* jacobian of the constraints */
    
    Vector params(n,x,true); 
    const OptimizerRep& rep = *reinterpret_cast<const OptimizerRep*>(user_data);



    Matrix jac(m,n);      

    if( rep.getNumericalJacobian() ) {
        Vector sfy0(m);            
        rep.getOptimizerSystem().constraintFunc( params, true, sfy0 );
        rep.jacDiff->calcJacobian( params, sfy0, jac);

    } else {

        rep.constraintJacobian( rep.getOptimizerSystem(), params, new_x, jac );

    }
    /* transpose the jacobian because Ipopt indexes in Row major format */
    Real *ptr = values;
    for(j=0;j<m;j++) {
        for(i=0;i<n;i++,ptr++) {
            *ptr = jac(j,i);
        }
    }

  } 
  return( true );
}
// TODO finish hessianWrapper
int hessianWrapper(int n, Real *x, int new_x, Real obj_factor,
            int m, Number *lambda, int new_lambda,
            int nele_hess, int *iRow, int *jCol,
            Real *values, void *user_data) {


    Vector coeff(n,x,true); 
    Vector hess(n*n,values,true); 
    const OptimizerRep& rep = *reinterpret_cast<const OptimizerRep*>(user_data);
    rep.hessian( rep.getOptimizerSystem(), coeff, new_x, hess );
    return( true );
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
