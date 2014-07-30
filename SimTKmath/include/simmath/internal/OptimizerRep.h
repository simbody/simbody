#ifndef SimTK_SIMMATH_OPTIMIZER_REP_H_
#define SimTK_SIMMATH_OPTIMIZER_REP_H_

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

#include "SimTKcommon.h"
#include "simmath/Optimizer.h"
#include "simmath/Differentiator.h"
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


class SimTK_SIMMATH_EXPORT Optimizer::OptimizerRep {
public:
    virtual ~OptimizerRep();
    OptimizerRep(const OptimizerSystem& sys) 
       : sysp(&sys), 
         myHandle(0), 
         cf(0),
         of(0),
         jacDiff(0),
         gradDiff(0),
         convergenceTolerance(Real(1e-3)),
         constraintTolerance(Real(1e-4)),
         maxIterations(1000),
         limitedMemoryHistory(50),
         diagnosticsLevel(0),
         diffMethod(Differentiator::CentralDifference),
         objectiveEstimatedAccuracy(SignificantReal),
         constraintsEstimatedAccuracy(SignificantReal),
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
         convergenceTolerance(Real(1e-3)),
         constraintTolerance(Real(1e-4)),
         maxIterations(1000),
         limitedMemoryHistory(50),
         diagnosticsLevel(0),
         diffMethod(Differentiator::CentralDifference),
         objectiveEstimatedAccuracy(SignificantReal),
         constraintsEstimatedAccuracy(SignificantReal),
         numericalGradient(false), 
         numericalJacobian(false)
    {
    }

    virtual OptimizerRep* clone() const { return 0; };
    static bool isAvailable() { return true; }

    virtual Real optimize(  Vector &results ) =  0;

    const OptimizerSystem& getOptimizerSystem() const {return *sysp;}

  
    void setDiagnosticsLevel( const int  level );
    void setConvergenceTolerance( Real accuracy );
    void setConstraintTolerance( Real tolerance );
    void setMaxIterations( const int iter );
    void setLimitedMemoryHistory( const int history );

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

    void useNumericalGradient(bool flag, Real objEstAccuracy); 
    void useNumericalJacobian(bool flag, Real consEstAccuracy);  
    void setDifferentiatorMethod( Differentiator::Method method);

    bool isUsingNumericalGradient() const { return numericalGradient; }
    bool isUsingNumericalJacobian() const { return numericalJacobian; }
    Differentiator::Method getDifferentiatorMethod() const {return diffMethod;}
    Real getEstimatedAccuracyOfObjective() const 
    {   return objectiveEstimatedAccuracy; }
    Real getEstimatedAccuracyOfConstraints() const 
    {   return constraintsEstimatedAccuracy; }

    const Differentiator& getGradientDifferentiator() const {
        assert(gradDiff);
        return *gradDiff;
    }

    const Differentiator& getJacobianDifferentiator() const {
        assert(jacDiff);
        return *jacDiff;
    }

    virtual OptimizerAlgorithm getAlgorithm() const {
        return UnknownOptimizerAlgorithm;
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
    Real constraintTolerance;
    int maxIterations;
    int limitedMemoryHistory;
    Differentiator::Method   diffMethod;
    Real objectiveEstimatedAccuracy;
    Real constraintsEstimatedAccuracy;

private:
    const OptimizerSystem* sysp;
    bool numericalGradient; // true if optimizer will compute an numerical gradient
    bool numericalJacobian; // true if optimizer will compute an numerical Jacobian
    Differentiator *gradDiff;   
    Differentiator *jacDiff; 

    SysObjectiveFunc  *of;   
    SysConstraintFunc *cf; 

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
    OptimizerAlgorithm getAlgorithm() const;
};

} // namespace SimTK


#endif  // SimTK_SIMMATH_OPTIMIZER_REP_H_
