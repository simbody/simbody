#ifndef SimTK_SIMMATH_OPTIMIZER_H_
#define SimTK_SIMMATH_OPTIMIZER_H_

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
#include "simmath/internal/common.h"
#include "simmath/Differentiator.h"

namespace SimTK {

enum OptimizerAlgorithm {
     BestAvailable  = 0, // Simmath will select best Optimizer based on problem type
     InteriorPoint  = 1, // IPOPT interior point optimizer
     LBFGS          = 2, // LBFGS optimizer
     LBFGSB         = 3, // LBFGS optimizer with simple bounds
     CFSQP          = 4  // CFSQP sequential quadratic programming optimizer (requires external library)
};

/**
 * Abstract class which defines an objective/cost function which is optimized by
 * and Optimizer object. The OptimizerSystem also defines any constraints which 
 * must be satisfied. 
 */
class SimTK_SIMMATH_EXPORT OptimizerSystem {
public:
    OptimizerSystem() : numParameters(0),
                        numEqualityConstraints(0),
                        numInequalityConstraints(0),
                        numLinearEqualityConstraints(0),
                        numLinearInequalityConstraints(0),
                        useLimits( false ),
                        lowerLimits(0),
                        upperLimits(0) { 
    }

    explicit OptimizerSystem(int nParameters ) { 
        new (this) OptimizerSystem(); // call the above constructor
        setNumParameters(nParameters);
    }

    virtual ~OptimizerSystem() {
        if( useLimits ) {
            delete lowerLimits;
            delete upperLimits;
        }
    }

    /// Objective/cost function which is to be optimized; return 0 when successful.
    /// This method must be supplied by concrete class.
    virtual int objectiveFunc      ( const Vector& parameters, 
                                 bool new_parameters, Real& f ) const {
                                 SimTK_THROW2(SimTK::Exception::UnimplementedVirtualMethod , "OptimizerSystem", "objectiveFunc" );
                                 return -1; }
  
    /// Computes the gradient of the objective function; return 0 when successful.
    /// This method does not have to be supplied if a numerical gradient is used.
    virtual int gradientFunc       ( const Vector &parameters, 
                                 bool new_parameters, Vector &gradient ) const  {
                                 SimTK_THROW2(SimTK::Exception::UnimplementedVirtualMethod , "OptimizerSystem", "gradientFunc" );
                                 return -1; }
    /// Computes the value of the constraints; return 0 when successful.
    /// This method must be supplied if the objective function has constraints.
    virtual int constraintFunc     ( const Vector & parameters, 
                                 bool new_parameters, Vector & constraints ) const {
                                 SimTK_THROW2(SimTK::Exception::UnimplementedVirtualMethod , "OptimizerSystem", "constraintFunc" );
                                 return -1; }
    /// Computes Jacobian of the constraints; return 0 when successful.
    /// This method does not have to be supplied if a numerical jacobian is used.
    virtual int constraintJacobian ( const Vector& parameters, 
                                  bool new_parameters, Matrix& jac ) const {
                                 SimTK_THROW2(SimTK::Exception::UnimplementedVirtualMethod , "OptimizerSystem", "constraintJacobian" );
                                 return -1; }
    /// Computes Hessian of the objective function; return 0 when successful.
    /// This method does not have to be supplied if limited memory is used.
    virtual int hessian            (  const Vector &parameters, 
                                 bool new_parameters, Vector &gradient) const {
                                 SimTK_THROW2(SimTK::Exception::UnimplementedVirtualMethod , "OptimizerSystem", "hessian" );
                                 return -1; }

   /// Sets the number of parameters in the objective function.
   void setNumParameters( const int nParameters ) {
       if(   nParameters < 1 ) {
           const char* where = " OptimizerSystem  Constructor";
           const char* szName = "number of parameters";
           SimTK_THROW5(SimTK::Exception::ValueOutOfRange, szName, 1, nParameters, INT_MAX, where);
       } else {
           numParameters = nParameters;
       }
   }
   /// Sets the number of equality constraints.
   void setNumEqualityConstraints( const int n ) {
       if( n < 0 ) {
           const char* where = " OptimizerSystem  setNumEqualityConstraints";
           const char* szName = "number of equality constraints";
           SimTK_THROW3(SimTK::Exception::SizeWasNegative, szName, n, where);
       } else {
           numEqualityConstraints = n;
       }
   }
   /// Sets the number of inequality constraints.
   void setNumInequalityConstraints( const int n ) {
       if( n < 0 ) {
           const char* where = " OptimizerSystem  setNumInequalityConstraints";
           const char* szName = "number of inequality constraints";
           SimTK_THROW3(SimTK::Exception::SizeWasNegative, szName, n, where);
       } else {
           numInequalityConstraints = n;
       }
   }
   /// Sets the number of lineaer equality constraints. 
   void setNumLinearEqualityConstraints( const int n ) {
       if( n < 0 || n > numEqualityConstraints ) {
           const char* where = " OptimizerSystem  setNumLinearEqualityConstraints";
           const char* szName = "number of linear equality constraints";
           SimTK_THROW4(SimTK::Exception::SizeOutOfRange, szName, n, numEqualityConstraints, where);
       } else {
           numLinearEqualityConstraints = n;
       }
   }
   /// Sets the number of lineaer inequality constraints.
   void setNumLinearInequalityConstraints( const int n ) {
       if( n < 0 || n > numInequalityConstraints ) {
           const char* where = " OptimizerSystem  setNumLinearInequalityConstraints";
           const char* szName = "number of linear inequality constraints";
           SimTK_THROW4(SimTK::Exception::SizeOutOfRange, szName, n, numInequalityConstraints, where);
       } else {
           numLinearInequalityConstraints = n;
       }
   }
   /// Set the upper and lower bounds on the paramters.
   void setParameterLimits( const Vector& lower, const Vector& upper  ) {
       if(   upper.size() != numParameters  && upper.size() != 0) {
           const char* where = " OptimizerSystem  setParamtersLimits";
           const char* szName = "upper limits length";
           SimTK_THROW5(Exception::IncorrectArrayLength, szName, upper.size(), "numParameters", numParameters, where);
       }
       if(   lower.size() != numParameters  && lower.size() != 0 ) {
           const char* where = " OptimizerSystem  setParamtersLimits";
           const char* szName = "lower limits length";
           SimTK_THROW5(Exception::IncorrectArrayLength, szName, lower.size(), "numParameters", numParameters, where);
       } 

       // set the upper and lower limits
       if( useLimits ) {
           delete lowerLimits;
           delete upperLimits;
       }

       if( upper.size() == 0 ) {
          useLimits = false;
       } else {
          lowerLimits = new Vector( lower );
          upperLimits = new Vector( upper );
          useLimits = true;
       }
   }

   /// Returns the number of parameters, that is, the number of variables that
   /// the Optimizer may adjust while searching for a solution.
   int getNumParameters() const {return numParameters;}
   /// Returns the total number of constraints.
   int getNumConstraints() const {return numEqualityConstraints+numInequalityConstraints;}
   /// Returns the number of equality constraints.
   int getNumEqualityConstraints() const {return numEqualityConstraints;}
   /// Returns the number of inequality constraints.
   int getNumInequalityConstraints() const {return numInequalityConstraints;}
   /// Returns the number of linear equality constraints.
   int getNumLinearEqualityConstraints() const {return numLinearEqualityConstraints;}
   /// Returns the number of nonlinear equality constraints.
   int getNumNonlinearEqualityConstraints() const {return numEqualityConstraints-numLinearEqualityConstraints;}
   /// Returns the number of linear inequality constraints.
   int getNumLinearInequalityConstraints() const {return numLinearInequalityConstraints;}
   /// Returns the number of linear inequality constraints.
   int getNumNonlinearInequalityConstraints() const {return numInequalityConstraints-numLinearInequalityConstraints;}

   /// Returns true if there are limits on the parameters.
   bool getHasLimits() const { return useLimits; }
   /// Returns the limits on the allowed values of each parameter, as
   /// an array of lower bounds and an array of upper bounds, with
   /// assumed lengths matching the number of parameters.
   void getParameterLimits( Real **lower, Real **upper ) const {
        *lower = &(*lowerLimits)[0];
        *upper = &(*upperLimits)[0];
   }

private:
   int numParameters;
   int numEqualityConstraints;
   int numInequalityConstraints;
   int numLinearEqualityConstraints;
   int numLinearInequalityConstraints;
   bool useLimits;
   Vector* lowerLimits;
   Vector* upperLimits;

}; // class OptimizerSystem

/**
 * API for SimTK Simmath's optimizers.
 * An optimizer finds a local minimum to an objective function. The
 * optimizer can be constrained to search for a minimum within a feasible 
 * region. The feasible region can be defined by setting limits on the 
 * parameters of the objective function and/or supplying constraint 
 * functions that must be satisfied. 
 * The optimizer starts searching for a minimum beginning at a user supplied 
 * initial value for the set of parameters.
 *
 * The objective function and constraints are specified by supplying the
 * Optimizer with a concrete implemenation of an OptimizerSystem class.
 * The OptimizerSystem can be passed to the Optimizer either through the 
 * Optimizer constructor or by calling the setOptimizerSystem method.  
 * The Optimizer class will select the best optimization algorithm to solve the
 * problem based on the constraints supplied by the OptimizerSystem. 
 * A user can also override the optimization algorithm selected by the
 * Optimizer by specifying the optimization algorithm. 
 *  
 */
class SimTK_SIMMATH_EXPORT Optimizer {
public:
    Optimizer();
    Optimizer( const OptimizerSystem& sys);
    Optimizer( const OptimizerSystem& sys, OptimizerAlgorithm algorithm);
    ~Optimizer();

    static bool isAlgorithmAvailable(OptimizerAlgorithm algorithm);
   
    /// Sets the relative accuracy used determine if the problem has converged.
    void setConvergenceTolerance(Real accuracy );
    /// Sets the absolute tolerance used to determine whether constraint
    /// violation is acceptable.
    void setConstraintTolerance(Real tolerance);


    /// Set the maximum number of iterations allowed of the optimization method's outer
    /// stepping loop. Most optimizers also have an inner loop ("line search") which is
    /// also iterative but is not affected by this setting. Inner loop convergence is
    /// typically prescribed by theory, and failure there is often an indication of
    /// an ill-formed problem.
    void setMaxIterations( int iter );
    /// Set the maximum number of previous hessians used in a limitied memory hessian approximation.
    void setLimitedMemoryHistory( int history );
    /// Set the level of debugging info displayed.
    void setDiagnosticsLevel( int level ); 

    void setOptimizerSystem( const OptimizerSystem& sys  );
    void setOptimizerSystem( const OptimizerSystem& sys, OptimizerAlgorithm algorithm );

    /// Set the value of an advanced option specified by a string.
    bool setAdvancedStrOption( const char *option, const char *value );
    /// Set the value of an advanced option specified by a real value.
    bool setAdvancedRealOption( const char *option, const Real value );
    /// Set the value of an advanced option specified by an integer value.
    bool setAdvancedIntOption( const char *option, const int value );
    /// Set the value of an advanced option specified by an boolean value.
    bool setAdvancedBoolOption( const char *option, const bool value );

    /// Enable numerical calculation of gradient, with optional estimation of
    /// the accuracy to which the objective function is calculated. For example,
    /// if you are calculate about 6 significant digits, supply the estimated
    /// accuracy as 1e-6. Providing the estimated accuracy improves the quality 
    /// of the calculated derivative. If no accuracy is provided we'll assume 
    /// the objective is calculated to near machine precision. See 
    /// SimTK::Differentiator for more information.
    void useNumericalGradient(bool flag, 
        Real estimatedAccuracyOfObjective = SignificantReal);
    /// Enable numerical calculation of the constraint Jacobian, with optional 
    /// estimation of the accuracy to which the constraint functions are 
    /// calculated.  For example, if you are calculating about 6 significant
    /// digits, supply the estimated accuracy as 1e-6. Providing the estimated 
    /// accuracy improves the quality of the calculated derivative. If no 
    /// accuracy is provided we'll assume the constraints are calculated to near
    /// machine precision. See SimTK::Differentiator for more information.
    void useNumericalJacobian(bool flag, 
        Real estimatedAccuracyOfConstraints = SignificantReal);

    /// Set which numerical gradient algorithm is used when numerical gradient
    /// or Jacobian is being used. Choices are Differentiator::ForwardDifference 
    /// (first order) or Differentiator::CentralDifference (second order) with 
    /// central the default.
    void setDifferentiatorMethod(Differentiator::Method method);

    /// Compute optimization.
    Real optimize(Vector&);

    /// Return a reference to the OptimizerSystem currently associated with this Optimizer.
    const OptimizerSystem& getOptimizerSystem() const;

    /// Indicate whether the Optimizer is currently set to use a numerical gradient.
    bool isUsingNumericalGradient() const;
    /// Indicate whether the Optimizer is currently set to use a numerical Jacobian.
    bool isUsingNumericalJacobian() const;
    /// Return the differentiation method that will be used if numerical
    /// gradient or Jacobian are required.
    Differentiator::Method getDifferentiatorMethod() const;
    /// Return the estimated accuracy last specified in useNumericalGradient().
    Real getEstimatedAccuracyOfObjective() const;
    /// Return the estimated accuracy last specified in useNumericalJacobian().
    Real getEstimatedAccuracyOfConstraints() const;

    // This is a local class.
    class OptimizerRep;
private:
    Optimizer( const Optimizer& c );
    Optimizer& operator=(const Optimizer& rhs);

    OptimizerRep* constructOptimizerRep(const OptimizerSystem&, OptimizerAlgorithm);
    const OptimizerRep& getRep() const {assert(rep); return *rep;}
    OptimizerRep&       updRep()       {assert(rep); return *rep;}

    // Hidden implementation to preserve binary compatibility.
    OptimizerRep* rep;

friend class OptimizerRep;
}; // class Optimizer
 
} // namespace SimTK

#endif //SimTK_SIMMATH_OPTIMIZER_H_

