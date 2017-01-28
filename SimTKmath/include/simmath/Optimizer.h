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

/**
 * The available Optimizer algorithms.
 * Gradient descent algorithms seek to find a local minimum, and are not
 * guaranteed to find the global minimum. See the description of Optimizer for
 * specific information about how to use the algorithms.
 */
enum OptimizerAlgorithm {
     /// Simmath will select best Optimizer based on problem type.
     BestAvailable = 0,
     /// IpOpt algorithm (https://projects.coin-or.org/ipopt);
     /// gradient descent.
     InteriorPoint = 1,
     /// Limited-memory Broyden-Fletcher-Goldfarb-Shanno algorithm; 
     /// gradient descent.
     LBFGS         = 2,
     /// LBFGS with simple bound constraints;
     /// gradient descent.
     LBFGSB        = 3,
     /// C implementation of sequential quadratic programming
     /// (requires external library:
     /// ftp://frcatel.fri.uniza.sk/pub/soft/math/matprog/doc/fsqp.html);
     /// gradient descent.
     CFSQP         = 4,
     /// Covariance matrix adaptation, evolution strategy
     /// (https://github.com/cma-es/c-cmaes);
     /// this is a randomized algorithm that attempts to find a global minimum.
     CMAES         = 5,
     UnknownOptimizerAlgorithm = 6, // the default impl. of getAlgorithm.
     /// An algorithm that is implemented outside of Simmath.
     UserSuppliedOptimizerAlgorithm = 7
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
    /// The value of f upon entry into the function is undefined.
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
   /// Set the upper and lower bounds on the parameters.
   void setParameterLimits( const Vector& lower, const Vector& upper  ) {
       if(   upper.size() != numParameters  && upper.size() != 0) {
           const char* where = " OptimizerSystem  setParametersLimits";
           const char* szName = "upper limits length";
           SimTK_THROW5(Exception::IncorrectArrayLength, szName, upper.size(), "numParameters", numParameters, where);
       }
       if(   lower.size() != numParameters  && lower.size() != 0 ) {
           const char* where = " OptimizerSystem  setParametersLimits";
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
 * An optimizer finds a minimum to an objective function. Usually, this minimum
 * is a local minimum. Some algorithms, like CMAES, are designed to find the
 * global minumum. The optimizer can be constrained to search for a minimum
 * within a feasible region. The feasible region is defined in two ways: via
 * limits on the parameters of the objective function; and, for algorithms
 * other than CMAES, by supplying constraint functions that must be satisfied.
 * The optimizer starts searching for a minimum beginning at a user supplied
 * initial value for the set of parameters.
 *
 * The objective function and constraints are specified by supplying the
 * Optimizer with a concrete implemenation of an OptimizerSystem class.
 * The OptimizerSystem can be passed to the Optimizer either through the
 * Optimizer constructor or by calling the Optimizer::setOptimizerSystem
 * method.  The Optimizer class will select the best optimization algorithm to
 * solve the problem based on the constraints supplied by the OptimizerSystem.
 * A user can also override the optimization algorithm selected by the
 * Optimizer by specifying the optimization algorithm.
 *
 * <h3> Optimization algorithms and advanced options </h3>
 *
 * See OptimizerAlgorithm for a brief description of the available algorithms.
 * Most of these algorithms have options that are specific to the algorithm.
 * These options are set via methods like Optimizer::setAdvancedStrOption. If
 * you want to get going quickly, you can just use the default values of these
 * options and ignore this section. As an example, an int option
 * <b>popsize</b> would be set via:
 *
 * @code
 * opt.setAdvancedIntOption("popsize", 5);
 * @endcode
 *
 * For now, we only have detailed documentation for the CMAES algorithm.
 *
 * <h4> CMAES </h4>
 *
 * This is the c-cmaes algorithm written by Niko Hansen
 * (https://github.com/cma-es/c-cmaes).
 *
 * Some notes:
 * - This algorithm obeys parameter limits.
 * - This is a derivative-free optimization algorithm, so methods like the
 *   following have no effect:
 *      - Optimizer::useNumericalGradient
 *      - Optimizer::setDifferentiatorMethod
 *      - Optimizer::setLimitedMemoryHistory
 *      - OptimizerSystem::gradientFunc
 *      - OptimizerSystem::hessian
 * - This algorithm does not obey constraint functions, so methods like the
 *   following have no effect:
 *      - Optimizer::setConstraintTolerance
 *      - Optimizer::useNumericalJacobian
 *      - OptimizerSystem::constraintFunc
 *      - OptimizerSystem::constraintJacobian
 *      - OptimizerSystem::setNumEqualityConstraints
 *      - OptimizerSystem::setNumInequalityConstraints
 *      - OptimizerSystem::setNumLinearEqualityConstraints
 *      - OptimizerSystem::setNumLinearInequalityConstraints
 * - The effect of the diagnostics level is as follows:
 *      - 0: minimal output to console (warnings, errors), some files are
 *      written to the current directory (errcmaes.err error log).
 *      - 1: additional output to console.
 *      - 2: all files are written to the current directory.
 *      - 3: output to console, and all files are written to the current directory.
 *
 *
 * Encoding of Variables
 *
 * Inappropriate initialization of the algorithm may lead to resampling of the
 * parameter distribution. Since some parameters may be defined in a bounded
 * region [a, b]. The CMA algorithm involves sampling a random distribution for
 * the variables (parameters) in the optimization problem. If the sampled values
 * do not lie within the variables' bounds, CMA must resample the distribution
 * until the variables lie within the bounds. This resampling is undesirable,
 * and may prevent the algorithm from functioning properly. There are two ways
 * to avoid excessive resampling:To overcome this the problem the user can
 * formulate the problem as follows:
 *
 * 1) Either define the bounds of the parameter and choose the appropriate stddev
 * and initial value for each parameter (see init_stepsize)
 * 2) Or to reformulate the problem, by rescaling the parameter space so that
 * each parameter map in a region between e.g. [0, 1] and each parameter has an
 * initial value of 0.5 and stddev of 0.2 (see comments below).
 *
 * The specific formulation of a (real) optimization problem has a tremendous
 * impact on the optimization performance. In particular, a reasonable parameter
 * encoding is essential. All parameters should be rescaled such that they have
 * presumably similar sensitivity (this makes the identity as initial covariance
 * matrix the right choice). Usually, the best approach is to write a wrapper
 * around the objective function that transforms the parameters before the actual
 * function call. The wrapper scales, for example, in each parameter/coordinate
 * the value [0; 10] into the typical actual domain of the parameter/coordinate.
 *
 * The natural encoding of (some of) the parameters can also be "logarithmic".
 * That is, for a parameter that must always be positive, with a ratio between
 * typical upper and lower value being larger than 100, we might use 10x instead
 * of x to call the objective function. More specifically, to achieve the parameter
 * range [10^–4,10^–1], we use 10^–4×10^3x/10 with x in [0; 10]. Again, the idea
 * is to have similar sensitivity: this makes sense if we expect the change from
 * 10^–4 to 10^–3 to have an impact similar to the change from 10^–2 to 10^–1.
 * In order to avoid the problem that changes of very small values have too less
 * an impact, an alternative is to choose 10^–1 × (x/10)2 >= 0. In the case where
 * only a lower bound at zero is necessary, a simple and natural transformation is
 * x2 × default_x, such that x=1 represents the default (or initial) value and x
 * remains unbounded during optimization.
 *
 * In summary, to map the values [0;10] into [a;b] we have the alternative
 * transformations a + (b-a) × x/10 or a + (b-a) × (x/10)2 >= a or a × (b/a)x/10 >= 0.
 *
 * For more details see: https://www.lri.fr/~hansen/cmaes_inmatlab.html
 *
 * Advanced options:
 * 
 * The default values for options whose name begins with "stop" are specified
 * at https://github.com/CMA-ES/c-cmaes/blob/master/cmaes_initials.par
 *
 * - <b>popsize</b> (int; default: depends on number of parameters) The
 *   population size (also known as lambda).
 * - <b>init_stepsize</b> (Vector/Real; default: 0.3) Initial stddev; After the
 *   encoding of variables, the initial solution point x0 and the initial standard
 *   deviation (step_size) sigma0 must be chosen. In a practical application, one
 *   often wants to start by trying to improve a given solution locally. In this
 *   case we choose a rather small sigma0 (say in [0.001, 0.1], given the x-values
 *   "live" in [0,10]). Thereby we can also check whether the initial solution is
 *   possibly a local optimum. When a global optimum is sought-after on rugged or
 *   multimodal landscapes, sigma0 should be chosen such that the final desirable
 *   location (or at least some of its domain of attraction) is not far outside of
 *   x0 ± 2sigma0 in each coordinate. (Remark that in Rn, if each boundary domain is
 *   in distance sigma, then the boundary corner is sigma*sqrt(n) away, which poses
 *   a slight dilemma for larger n.)
 *
 *   A warning is emitted if this is not set and default value is used for each
 *   variable.
 *
 *   Example setting the init_stepsize:
 *
 *   Vector initStepSize(N, 0.3);
 *   opt.setAdvancedVectorOption("init_stepsize", initStepSize);
 *   "or"
 *   opt.setAdvancedRealOption("init_stepsize", 0.3);
 * - <b>seed</b> (int; default: 0, which uses clock time) Seed for the random
 *   number generator that is used to sample the population from a normal
 *   distribution. See note below.
 * - <b>maxTimeFractionForEigendecomposition</b> (real; default: 0.2)
 *   Controls the amount of time spent generating eigensystem
 *   decompositions.
 * - <b>stopMaxFunEvals</b> (int) Stop optimization after this
 *   number of evaluations of the objective function.
 * - <b>stopFitness</b> (real) Stop if function value is smaller than
 *   stopFitness.
 * - <b>stopTolFunHist</b> (real) Stop if function value differences of best
 *   values are smaller than stopTolFunHist.
 * - <b>stopTolX</b> (real) Stop if step sizes are smaller than stopTolX.
 * - <b>stopTolUpXFactor</b> (real) Stop if standard deviation increases
 *   by more than stopTolUpXFactor.
 * - <b>parallel</b> (str) To run the optimization with multiple threads, set
 *   this to "multithreading". Only use this if your OptimizerSystem is
 *   threadsafe: you can't reliably modify any mutable variables in your
 *   OptimizerSystem::objectiveFun().
 * - <b>nthreads</b> (int) If the <b>parallel</b> option is set to
 *   "multithreading", this is the number of threads to use (by default, this
 *   is the number of processors/threads on the machine).
 *
 * If you want to generate identical results with repeated optimizations,
 * you can set the <b>seed</b> option. In addition, you *must* set the
 * <b>maxTimeFractionForEigendecomposition</b> option to be greater than or
 * equal to 1.0.
 *
 * @code
 * opt.setAdvancedIntOption("seed", 42);
 * opt.setAdvancedRealOption("maxTimeFractionForEigendecomposition", 1);
 * @endcode
 *
 */
class SimTK_SIMMATH_EXPORT Optimizer {
public:
    Optimizer();
    Optimizer( const OptimizerSystem& sys);
    Optimizer( const OptimizerSystem& sys, OptimizerAlgorithm algorithm);
    ~Optimizer();

    /// BestAvailable, UnknownAlgorithm, and UserSuppliedAlgorithm
    /// are treated as never available.
    static bool isAlgorithmAvailable(OptimizerAlgorithm algorithm);
   
    /// Sets the relative accuracy used determine if the problem has converged.
    void setConvergenceTolerance(Real accuracy );
    /// Sets the absolute tolerance used to determine whether constraint
    /// violation is acceptable.
    void setConstraintTolerance(Real tolerance);


    /// Set the maximum number of iterations allowed of the optimization
    /// method's outer stepping loop. Most optimizers also have an inner loop
    /// ("line search") which is also iterative but is not affected by this
    /// setting. Inner loop convergence is typically prescribed by theory, and
    /// failure there is often an indication of an ill-formed problem.
    void setMaxIterations( int iter );
    /// Set the maximum number of previous hessians used in a limited memory
    /// hessian approximation.
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
    /// Set the value of an advanced option specified by an Vector value.
    bool setAdvancedVectorOption( const char *option, const Vector value );

    
    /// Set which numerical differentiation algorithm is to be used for the next
    /// useNumericalGradient() or useNumericalJacobian() call. Choices are 
    /// Differentiator::ForwardDifference (first order) or 
    /// Differentiator::CentralDifference (second order) with central the 
    /// default.
    /// @warning This has no effect if you have already called 
    /// useNumericalGradient() or useNumericalJacobian(). Those must be called
    /// \e after setDifferentiatorMethod().
    /// @see SimTK::Differentiator
    void setDifferentiatorMethod(Differentiator::Method method);
    /// Return the differentiation method last supplied in a call to
    /// setDifferentiatorMethod(), \e not necessarily the method currently
    /// in use. See setDifferentiatorMethod() for more information.
    /// @see SimTK::Differentiator
    Differentiator::Method getDifferentiatorMethod() const;

    /// Return the algorithm used for the optimization. You may be interested
    /// in this value if you didn't specify an algorithm, or specified for
    /// Simbody to choose the BestAvailable algorithm. This method won't return
    /// BestAvailable, even if it's the 'algorithm' that you chose.
    OptimizerAlgorithm getAlgorithm() const;

    /// Enable numerical calculation of gradient, with optional estimation of
    /// the accuracy to which the objective function is calculated. For example,
    /// if you are calculate about 6 significant digits, supply the estimated
    /// accuracy as 1e-6. Providing the estimated accuracy improves the quality 
    /// of the calculated derivative. If no accuracy is provided we'll assume 
    /// the objective is calculated to near machine precision. The method used
    /// for calculating the derivative will be whatever was \e previously 
    /// supplied in a call to setDifferentiatorMethod(), or the default which
    /// is to use central differencing (two function evaluations per 
    /// gradient entry). See SimTK::Differentiator for more information.
    /// @see setDifferentiatorMethod(), SimTK::Differentiator
    void useNumericalGradient(bool flag, 
        Real estimatedAccuracyOfObjective = SignificantReal);
    /// Enable numerical calculation of the constraint Jacobian, with optional 
    /// estimation of the accuracy to which the constraint functions are 
    /// calculated.  For example, if you are calculating about 6 significant
    /// digits, supply the estimated accuracy as 1e-6. Providing the estimated 
    /// accuracy improves the quality of the calculated derivative. If no 
    /// accuracy is provided we'll assume the constraints are calculated to near
    /// machine precision.  The method used for calculating the derivative will 
    /// be whatever was \e previously supplied in a call to 
    /// setDifferentiatorMethod(), or the default which is to use central 
    /// differencing (two function evaluations per Jacobian column. See 
    /// SimTK::Differentiator for more information.
    /// @see setDifferentiatorMethod(), SimTK::Differentiator
    void useNumericalJacobian(bool flag, 
        Real estimatedAccuracyOfConstraints = SignificantReal);

    /// Compute optimization.
    Real optimize(Vector&);

    /// Return a reference to the OptimizerSystem currently associated with this Optimizer.
    const OptimizerSystem& getOptimizerSystem() const;

    /// Indicate whether the Optimizer is currently set to use a numerical gradient.
    bool isUsingNumericalGradient() const;
    /// Indicate whether the Optimizer is currently set to use a numerical Jacobian.
    bool isUsingNumericalJacobian() const;
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

