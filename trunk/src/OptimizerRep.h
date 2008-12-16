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
#include "Optimizer.h"
#include "Differentiator.h"
#include <map>

namespace SimTK {


/*  class for Diff jacobian */
class SysObjectiveFunc : public Differentiator::GradientFunction {
public:
    SysObjectiveFunc(int ny, const OptimizerSystem* sysPtr )
        : Differentiator::GradientFunction(ny) { sysp = sysPtr; }

    // Must provide this pure virtual function.
    int f(const Vector& y, Real& fy) const  {
         return(sysp->objectiveFunc(y, true, fy));   // class user's objectiveFunc
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
       return(sysp->constraintFunc(y, true, fy));  // calls user's contraintFunc
    }
    const OptimizerSystem* sysp;
};


class Optimizer::OptimizerRep {
public:
    virtual ~OptimizerRep();
    OptimizerRep(const OptimizerSystem& sys) 
       : sysp(0), 
         myHandle(0), 
         cf(0),
         of(0),
         jacDiff(0),
         gradDiff(0),
         convergenceTolerance(1e-4),
         maxIterations(100),
         limitedMemoryHistory(5),
         diagnosticsLevel(0),
         diffMethod(Differentiator::CentralDifference),
         numericalGradient(false), 
         numericalJacobian(false)

    {
    }
    OptimizerRep()
       : sysp(0), 
         myHandle(0), 
         cf(0),
         of(0),
         jacDiff(0),
         gradDiff(0),
         convergenceTolerance(1e-4),
         maxIterations(100),
         limitedMemoryHistory(5),
         diagnosticsLevel(0),
         diffMethod(Differentiator::CentralDifference),
         numericalGradient(false), 
         numericalJacobian(false)

    {
    }

    virtual OptimizerRep* clone() const { return 0; };
    static bool isAvailable() { return true; }

    virtual Real optimize(  Vector &results ) =  0;

    const OptimizerSystem& getOptimizerSystem() const {return *sysp;}

  
    void setDiagnosticsLevel( const int  level );
    void setConvergenceTolerance( const Real tolerance );
    void setMaxIterations( const int iter );
    void setLimitedMemoryHistory( const int history );
    void setDifferentiatorMethod( Differentiator::Method method);

    bool setAdvancedStrOption( const std::string &option, const std::string &value );
    bool setAdvancedRealOption( const std::string &option, const Real value );
    bool setAdvancedIntOption( const std::string &option, const int value );
    bool setAdvancedBoolOption( const std::string &option, const bool value );

    bool getAdvancedStrOption( const std::string &option, std::string &value ) const;
    bool getAdvancedRealOption( const std::string &option, Real &value ) const;
    bool getAdvancedIntOption( const std::string &option, int &value ) const;
    bool getAdvancedBoolOption( const std::string &option, bool &value ) const;

    void  setMyHandle(Optimizer& cp) {myHandle = &cp;}
    const Optimizer& getMyHandle() const {assert(myHandle); return *myHandle;}
    void  clearMyHandle() {myHandle=0;} 
    void useNumericalGradient( const bool flag ); 
    void useNumericalJacobian( const bool flag );  
    bool isUsingNumericalGradient() const { return numericalGradient; }
    bool isUsingNumericalJacobian() const { return numericalJacobian; }

    const Differentiator& getGradientDifferentiator() const {
        assert(gradDiff);
        return *gradDiff;
    }

    const Differentiator& getJacobianDifferentiator() const {
        assert(jacDiff);
        return *jacDiff;
    }

    static int numericalGradient_static( const OptimizerSystem&, const Vector & parameters,  const bool new_parameters,  Vector &gradient );
    static int numericalJacobian_static(const OptimizerSystem&,
                                   const Vector& parameters, const bool new_parameters, Matrix& jacobian );

protected:
    // These methods are to be called by derived classes as an interface
    // to the OptimizerSystem virtuals. The signature must match that required by
    // IpOpt's matching callbacks. We're using the "user data" argument to pass in
    // the current OptimizerRep, making these behave like non-static members.

    static int objectiveFuncWrapper ( int n, const Real* x, int new_x, Real* f, void* rep);
    static int gradientFuncWrapper  ( int n, const Real* x, int new_x, Real* gradient, void* rep);
    static int constraintFuncWrapper( int n, const Real* x, int new_x, int m, Real* g, void* rep);
    static int constraintJacobianWrapper( int n, const Real* x, int new_x,int m, int nele_jac,
                                          int* iRow, int* jCol, Real* values, void* rep);
    static int hessianWrapper(  int n, const Real* x, int new_x, Real obj_factor,
                                int m, Real* lambda, int new_lambda,
                                int nele_hess, int* iRow, int* jCol,
                                Real* values, void* rep);

    int diagnosticsLevel;
    Real convergenceTolerance;
    int maxIterations;
    int limitedMemoryHistory;
    Differentiator::Method   diffMethod;

private:
    const OptimizerSystem* sysp;
    bool numericalGradient; // true if optimizer will compute an numerical gradient
    bool numericalJacobian; // true if optimizer will compute an numerical Jacobian
    Differentiator *gradDiff;   
    Differentiator *jacDiff; 

    SysObjectiveFunc  *of;   
    SysConstraintFunc *cf; 
    void initNumericalJac();
    void disableNumericalJac();
    void initNumericalGrad();
    void disableNumericalGrad();

    std::map<std::string, std::string> advancedStrOptions;
    std::map<std::string, Real> advancedRealOptions;
    std::map<std::string, int> advancedIntOptions;
    std::map<std::string, bool> advancedBoolOptions;

    friend class Optimizer;
    Optimizer* myHandle;   // The owner handle of this Rep.
    
}; // end class OptimizerRep

class DefaultOptimizer: public Optimizer::OptimizerRep {
    Real optimize(  Vector &results );
    OptimizerRep* clone() const;
};

} // namespace SimTK


#endif  //_SimTK_OPTIMIZER_REP_H_
