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
#include "SimTKcommon/internal/BigMatrix.h"
#include "Optimizer.h"
#include "Differentiator.h"


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

    /*  class for Diff jacobian */
    class SysObjectiveFunc  : public Differentiator::GradientFunction {
public:
    SysObjectiveFunc(int ny, const OptimizerSystem* sysPtr )
        : Differentiator::GradientFunction(ny) { sysp = sysPtr; }

    // Must provide this pure virtual function.
    int f(const Vector& y, Real& fy) const  {
         sysp->objectiveFunc(y, true, fy);   // class user's objectiveFunc
    }
    const OptimizerSystem* sysp;
};


    /*  class for Diff gradient */
    class SysConstraintFunc : public Differentiator::JacobianFunction {
       public:
    SysConstraintFunc(int nf, int ny, const OptimizerSystem* sysPtr)
        : Differentiator::JacobianFunction(nf,ny) { sysp = sysPtr; }

    // Must provide this pure virtual function.
    int f(const Vector& y, Vector& fy) const  {
       sysp->constraintFunc(y, true, fy);  // calls user's contraintFunc
    }
    const OptimizerSystem* sysp;
};


class OptimizerRep {
public:
    ~OptimizerRep();
    OptimizerRep(OptimizerSystem& sys) 
       : sysp(0), 
         myHandle(0), 
         cf(0),
         of(0),
         jacDiff(0),
         gradDiff(0),
         convergenceTolerance(1e-4),
         diagnosticsLevel(0),
         numericalGradient(false), 
         numericalJacobian(false)

    {
       zeroFunctionPointers();
    }

    virtual Real optimize(  Vector &results ) =  0;

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

    Differentiator *gradDiff;   
    Differentiator *jacDiff;   
    void setDiagnosticsLevel( const int  level );
    void setConvergenceTolerance( const Real tolerance );
    int setAdvancedOptions( const char *option, const Real *values );
    void  setMyHandle(Optimizer& cp) {myHandle = &cp;}
    const Optimizer& getMyHandle() const {assert(myHandle); return *myHandle;}
    void  clearMyHandle() {myHandle=0;} 
    void useNumericalGradient( const bool flag ); 
    void useNumericalJacobian( const bool flag );  
    bool getNumericalGradient() const { return( numericalGradient ); }
    bool getNumericalJacobian() const { return( numericalJacobian ); }
    static int numericalGradient_static( const OptimizerSystem&, const Vector & parameters,  const bool new_parameters,  Vector &gradient );
    static int numericalJacobian_static(const OptimizerSystem&,
                                   const Vector& parameters, const bool new_parameters, Matrix& jacobian );

    protected:
    int diagnosticsLevel;
    Real convergenceTolerance;

    private:
    OptimizerSystem* sysp;
    bool numericalGradient; // true if optimizer will compute an numerical gradient
    bool numericalJacobian; // true if optimizer will compute an numerical gradient

    SysObjectiveFunc  *of;   
    SysConstraintFunc *cf; 
    void initNumericalJac();
    void disableNumericalJac();
    void initNumericalGrad();
    void disableNumericalGrad();

    friend class Optimizer;
    Optimizer* myHandle;   // The owner handle of this Rep.
    
}; // end class OptimizerRep
} // namespace SimTK
#endif  //_SimTK_OPTIMIZER_REP_H_
