
#ifndef _SimTK_OPTIMIZER_H
#define _SimTK_OPTIMIZER_H

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

#include <limits.h>
#include "Simmath.h"
#include "simmatrix/internal/BigMatrix.h"

namespace SimTK {

class OptimizerSystem {
public:
    OptimizerSystem(int nParameters ) : numConstraints(0),
                                      numEqualityConstraints(0),  
                                      lowerLimits(0),
                                      upperLimits(0),
                                      useLimits( false ) { 
       if(   nParameters < 1 ) {
           char *where = " OptimizerSystem  Constructor";
           char *szName= "number of parameters";
           SimTK_THROW5(SimTK::Exception::ValueOutOfRange, szName, 1,  INT_MAX, nParameters, where);
       }else {
            numParameters = nParameters;
       }
    }

    OptimizerSystem(int nParameters, int nConstraints) : numEqualityConstraints(0),
                                                         upperLimits(0),
                                                         lowerLimits(0),
                                                         useLimits( false ) { 

       if(   nParameters < 1 ) {
           char *where = " OptimizerSystem  Constructor";
           char *szName= "number of parameters";
           SimTK_THROW5(SimTK::Exception::ValueOutOfRange, szName, 1,  INT_MAX, nParameters, where);
       }else {
            numParameters = nParameters;
       }

       if(   nConstraints < 0 ) {
           char *where = " OptimizerSystem  Constructor";
           char *szName= "number of constraints";
           SimTK_THROW4(SimTK::Exception::SizeOutOfRange, szName,  INT_MAX, nConstraints, where);
       }else {
            numConstraints = nConstraints;
       }
    }

  /* this method must be supplied by concreate class */
  virtual int objectiveFunc      ( const Vector& parameters, 
                                 const bool new_parameters, Real& f ) const {
                                 SimTK_THROW2(SimTK::Exception::UnimplementedVirtualMethod , "OptimizerSystem", "objectiveFunc" );
                                 return -1; }


  virtual int gradientFunc       ( const Vector &parameters, 
                                 const bool new_parameters, Vector &gradient ) const  {
                                 SimTK_THROW2(SimTK::Exception::UnimplementedVirtualMethod , "OptimizerSystem", "gradientFunc" );
                                 return -1; }

  virtual int constraintFunc     ( const Vector & parameters, 
                                 const bool new_parameters, Vector & constraints ) const {
                                 SimTK_THROW2(SimTK::Exception::UnimplementedVirtualMethod , "OptimizerSystem", "constraintFunc" );
                                 return -1; }

  virtual int constraintJacobian ( const Vector& parameters, 
                                  const bool new_parameters, Vector& jac ) const {
                                 SimTK_THROW2(SimTK::Exception::UnimplementedVirtualMethod , "OptimizerSystem", "constraintJacobian" );
                                 return -1; }

  virtual int hessian            (  const Vector &parameters, 
                                 const bool new_parameters, Vector &gradient) const {
                                 SimTK_THROW2(SimTK::Exception::UnimplementedVirtualMethod , "OptimizerSystem", "hessian" );
                                 return -1; }

   void setNumEqualityConstraints( const int n ) {
 
        if( n < 0 || n > numConstraints ) {
           char *where = " OptimizerSystem  setNumberOfEqualityConstraints";
           char *szName= "Number of Constraints";
           SimTK_THROW4(SimTK::Exception::SizeOutOfRange, szName,  n, numConstraints, where);
        } else {
           numEqualityConstraints = n;
        }

        return;
   }
   void setParameterLimits( const Vector& lower, const Vector& upper  ) {
 
        if(   upper.size() != numParameters  && upper.size() != 0) {
           char *where = " OptimizerSystem  setParamtersLimits";
           char *szName= "upper limits length";
           SimTK_THROW5(Exception::IncorrectArrayLength, szName, upper.size(), "numParameters", numParameters, where);
        }
        if(   lower.size() != numParameters  && lower.size() != 0 ) {
           char *where = " OptimizerSystem  setParamtersLimits";
           char *szName= "lower limits length";
           SimTK_THROW5(Exception::IncorrectArrayLength, szName, lower.size(), "numParameters", numParameters, where);
        } 
       // set the upper and lower limits
       if( lowerLimits ) delete lowerLimits;
       if( upperLimits ) delete upperLimits;

       if( upper.size() == 0 ) {
          useLimits = false;
       } else {
          lowerLimits = new Vector( lower );
          upperLimits = new Vector( upper );
          useLimits = true;
       }

       return;
   }
   int getNumParameters() const {return numParameters;};
   int getNumConstraints() const {return numConstraints;};
   int getNumEqualityConstraints() const {return numEqualityConstraints;};
   bool getHasLimits() const { return useLimits; }
   void getParameterLimits( double **lower, double **upper ) const {
        *lower = &(*lowerLimits)[0];
        *upper = &(*upperLimits)[0];
   }


   int numParameters;
   int numConstraints;
   int numEqualityConstraints;
   bool useLimits;
   Vector* lowerLimits;
   Vector* upperLimits;

}; // class OptimizerSystem

// These static functions are private to the current (client-side) compilation
// unit. They are used to navigate the client-side OptimizerSystem virtual function
// table, which cannot be done on the library side. Note that these are defined
// in the SimTK namespace so don't need "SimTK" in their names.
static int objectiveFunc_static(const OptimizerSystem& sys,
                                 const Vector& parameters, 
                                const bool new_parameters, Real& f ) {

    return sys.objectiveFunc(parameters, new_parameters,  f);
}
static int gradientFunc_static(const OptimizerSystem& sys,
                                const Vector &parameters, 
                               const bool new_parameters, Vector &gradient ) {
    return sys.gradientFunc( parameters, new_parameters, gradient);
}
static int constraintFunc_static(const OptimizerSystem& sys,
                                 const Vector &parameters, 
                                 const bool new_parameters, Vector& constraints ) {
    return sys.constraintFunc( parameters, new_parameters, constraints);
}
static int constraintJacobian_static(const OptimizerSystem& sys,
                                 const Vector &parameters, 
                                 const bool new_parameters, Vector& jac ) {
    return sys.constraintJacobian( parameters, new_parameters, jac);
}
static int hessian_static(const OptimizerSystem& sys,
                                 const Vector &parameters,
                                 const bool new_parameters, Vector& gradient ) {
    return sys.hessian( parameters, new_parameters, gradient);
}

/*
** Class for API interface to Simmath's optimizers.
** The OptimizerSystem class describes the optimization by
** specifying the objective function and constraints. OptimizerFactory()
** instantiates the correct optimizer based on the objective function 
** and constraints specified in the OptimizerSystem object. 
** If the user calls the Optimizer constructor and 
** supplies the algorithm argument the OptimizerFactory() will ignore the 
** will create instatiate the Optimizer asked for.
**  
*/

class Optimizer  {

   public:
    Optimizer( OptimizerSystem& sys) {
        // Perform construction of the CPodesRep on the library side.
        librarySideOptimizerConstructor(sys);
        // But fill in function pointers from the client side.
        clientSideOptimizerConstructor();
    }

    ~Optimizer();

    void setOptimizerParameters(unsigned int param, double *values); 
    void getOptimizerParameters(unsigned int param, double *values);

    double optimize(Vector&);
    
private:
  typedef int (*ObjectiveFunc)      ( const OptimizerSystem&,
                                    const Vector& parameters,  const bool new_parameters,
                                    Real& f );

  typedef int (*GradientFunc)       ( const OptimizerSystem&,
                                    const Vector &parameters, const bool new_parameters,
                                    Vector &gradient );

  typedef int (*ConstraintFunc)     ( const OptimizerSystem&,
                                    const Vector & parameters, const bool new_parameters,
                                    Vector &constraints );

  typedef int (*ConstraintJacobian) ( const OptimizerSystem&,
                                    const Vector& parameters, const bool new_parameters,
                                    Vector& jac );

  typedef int (*Hessian)            ( const OptimizerSystem&,
                                    const Vector &parameters, const bool new_parameters,
                                    Vector &gradient);

  void registerObjectiveFunc( ObjectiveFunc );
  void registerGradientFunc( GradientFunc );
  void registerConstraintFunc( ConstraintFunc );
  void registerConstraintJacobian( ConstraintJacobian );
  void registerHessian( Hessian );

    // This is the library-side part of the CPodes constructor. This must
    // be done prior to the client side construction.
  void librarySideOptimizerConstructor(OptimizerSystem& sys);

  void clientSideOptimizerConstructor() {
       registerObjectiveFunc( objectiveFunc_static );
       registerGradientFunc( gradientFunc_static );
       registerConstraintFunc( constraintFunc_static );
       registerConstraintJacobian( constraintJacobian_static );
       registerHessian( hessian_static );
  }

    class OptimizerRep* rep;
    friend class OptimizerRep;

    const OptimizerRep& getRep() const {assert(rep); return *rep;}
    OptimizerRep&       updRep()       {assert(rep); return *rep;}

    // Suppress copy constructor and default assigment operator.
    Optimizer(const Optimizer&);
    Optimizer& operator=(const Optimizer&);

}; // Class Optimizer
 
} // namespace SimTK

#endif //_SimTK_OPTIMIZER_H

