#ifndef _SimTK_OPTIMIZER_REP_H_
#define _SimTK_OPTIMIZER_REP_H_

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
#include "SimTKcommon.h"
#include "SimTKcommon/internal/common.h"
#include "simmatrix/internal/BigMatrix.h"
#include "Optimizer.h"


namespace SimTK {
extern int objectiveFuncWrapper( int n, Real *x, int new_x,  Real *f, void*user_data);
extern int gradientFuncWrapper( int n,  Real *x, int new_x, Real *gradient, void*user_data);
extern int constraintFuncWrapper( int n, Real *x, int new_x, int m, Real *g,  void*user_data);
extern int constraintJacobianWrapper( int n, Real *x, int new_x,int m, int nele_jac,
                int *iRow, int *jCol, Real *values, void *user_data);
extern int hessianWrapper(int n, Real *x, int new_x, Real obj_factor,
            int m, Real *lambda, int new_lambda,
            int nele_hess, int *iRow, int *jCol,
            Real *values, void *user_data);

class OptimizerRep {
public:
    OptimizerRep(OptimizerSystem& sys) 
       : sysp(0), myHandle(0)
    {
       zeroFunctionPointers();
    }

//    virtual ~OptimizerRep(){};

    virtual void setOptimizerParameters(unsigned int parameter, double *values ){};
    virtual void getOptimizerParameters(unsigned int parameter, double *values ){};
    virtual double optimize(  Vector &results ) =  0;

    const OptimizerSystem& getOptimizerSystem() const {return *sysp;}

    // Client-side function pointers
    Optimizer::ObjectiveFunc       objectiveFunc; // points to objectiveFunc_static
    Optimizer::GradientFunc        gradientFunc;
    Optimizer::ConstraintFunc      constraintFunc;
    Optimizer::ConstraintJacobian  constraintJacobian;
    Optimizer::Hessian             hessian;

    void zeroFunctionPointers() {
        objectiveFunc      = 0;
        gradientFunc       = 0;
        constraintFunc     = 0;
        constraintJacobian = 0;
        hessian            = 0;
    }

    void  setMyHandle(Optimizer& cp) {myHandle = &cp;}
    const Optimizer& getMyHandle() const {assert(myHandle); return *myHandle;}
    void  clearMyHandle() {myHandle=0;} private:
    const OptimizerSystem* sysp;

    friend class Optimizer;
    Optimizer* myHandle;   // The owner handle of this Rep.
    
}; // end class OptimizerRep
} // namespace SimTK
#endif  //_SimTK_OPTIMIZER_REP_H_
