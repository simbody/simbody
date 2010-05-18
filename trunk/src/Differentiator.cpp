/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-10 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/** @file
 * This is the private (library side) implementation of the Simmath
 * Differentiator family of classes.
 */

#include "SimTKcommon.h"
#include "simmath/Differentiator.h"

#include <exception>

namespace SimTK {

    ///////////////////////////////
    // DIFFERENTIATOR EXCEPTIONS //
    ///////////////////////////////

class Differentiator::OpNotAllowedForFunctionOfThisShape : public Exception::Base {
public:
    OpNotAllowedForFunctionOfThisShape(const char* fn, int ln, 
                                       String op, String reqShape, String funcKind,
                                       int nf, int ny) : Base(fn,ln) 
    {
        setMessage("Differentiator method " + op + "() requires a " + reqShape + 
            " function but was called with a " + funcKind + " that is " + 
            String(nf) + "x" + String(ny));
    }
};

class Differentiator::UserFunctionThrewAnException : public Exception::Base {
public:
    UserFunctionThrewAnException(const char* fn, int ln, 
                                 const char* msg) : Base(fn,ln) {
        setMessage("A user function threw an exception when invoked by Differentiator."
            "  The exception message was: " + String(msg));
    }
};

class Differentiator::UserFunctionReturnedNonzeroStatus : public Exception::Base {
public:
    UserFunctionReturnedNonzeroStatus(const char* fn, int ln, int status) : Base(fn,ln) {
        setMessage("A user function returned non-zero status " + String(status)
                   + " when invoked by Differentiator.");
    }
};

class Differentiator::UnknownMethodSpecified : public Exception::Base {
public:
    UnknownMethodSpecified(const char* fn, int ln, int enumValue,
                           const char* classname, const char* methodname) : Base(fn,ln) {
        setMessage("An unrecognized Differentiator::Method with enum value " + String(enumValue)
                   + " was passed to " + String(classname) + "::" + String(methodname) +"()");
    }
};



// We want to perturb y0 by h such that the perturbation is
// exactly representable in binary. hEst is our first guess
// at h, which we calculate h=(y0+hEst)-y0 taking care not
// to let a clever compiler optimize that away.
static Real cleanUpH(Real hEst, Real y0) {
    volatile Real temp = y0+hEst; 
    return temp-y0;
}

static void throwIfMethodInvalid(Differentiator::Method m, const char* op) {
    if (!Differentiator::isValidMethod(m))
        SimTK_THROW3(Differentiator::UnknownMethodSpecified, (int)m,
            "Differentiator", op);
}

// Operation "op" has been given the indicated method m. If m
// is a valid method, return it. If m is unspecified, return the
// indicated default. Otherwise, throw a "bad method" exception.
static Differentiator::Method getMethodOrThrow
   (Differentiator::Method m, Differentiator::Method def, const char* op) 
{
    throwIfMethodInvalid(m,op);
    return (m==Differentiator::UnspecifiedMethod ? def : m);
}

    ///////////////////////////////////////////////////////////
    // REP CLASS DECLARATIONS FOR DIFFERENTIATOR & FUNCTIONS //
    ///////////////////////////////////////////////////////////

// These classes are just local definitions in this compilation unit; they
// aren't visible anywhere else.
class ScalarFunctionRep;
class GradientFunctionRep;
class JacobianFunctionRep;

// This is used as a value for y in calculating the step size when
// the actual y is smaller.
static const Real YMin = Real(0.1);
class Differentiator::DifferentiatorRep {
public:
    DifferentiatorRep(Differentiator* handle,
                      const Differentiator::Function::FunctionRep&,
                      Differentiator::Method defaultMethod);
    // default destructor, no default constructor, no copy or copy assign

    // This constant is the algorithm we'll use by default.
    static const Differentiator::Method DefaultDefaultMethod 
        = Differentiator::ForwardDifference;

    void calcDerivative(const ScalarFunctionRep&, Differentiator::Method, 
                        Real y0, Real fy0, Real& dfdy) const;
    void calcGradient(const GradientFunctionRep&, Differentiator::Method, 
                      const Vector& y0, Real fy0, Vector& gf)   const;
    void calcJacobian(const JacobianFunctionRep&, Differentiator::Method, 
                      const Vector& y0, const Vector& fy0, Matrix& dfdy) const;

    const Real& getAccFac(int order) const {
        if (order==1) return AccFac1;
        if (order==2) return AccFac2;
        assert(!"Unrecognized Differentiator order");
        return NaN;
    }

    void resetAllStatistics() {
        nDifferentiations = nDifferentiationFailures = nCallsToUserFunction = 0;
    }

    // Statistics
    mutable int nDifferentiations; 
    mutable int nDifferentiationFailures; 
    mutable int nCallsToUserFunction;
private:
    Differentiator* myHandle;
    friend class Differentiator;

    // Problem dimensions. These are obtained from the Function at
    // Differentiator construction time and cannot change afterwards.

    const Differentiator::Function::FunctionRep& frep;
    const int  NFunctions, NParameters;
    const Real EstimatedAccuracy;

    // This is set on construction, but can be changed.
    Differentiator::Method defaultMethod;

    // These are pre-calculated accuracy factors for 1st order and
    // 2nd order step size estimates, derived from EstimatedAccuracy
    // upon construction.
    const Real AccFac1, AccFac2;

    // These temporaries are kept here so we can reuse their storage space
    // from call to call once they have been allocated at construction. 
    // The *values* do not persist across calls.
    mutable Vector ytmp;           // [NParameters]
    mutable Vector fyptmp, fymtmp; // [NFunctions]

    // suppress
    DifferentiatorRep(const DifferentiatorRep&);
    DifferentiatorRep& operator=(const DifferentiatorRep&);
};

Differentiator::~Differentiator() {
    delete rep;
}

class Differentiator::Function::FunctionRep {
    friend class Differentiator::Function;
public:
    FunctionRep(int nf, int np, Real acc)
      : nFunc(nf), nParam(np), estimatedAccuracy(acc)
    {
        if (estimatedAccuracy < 0) // use default
            estimatedAccuracy = SignificantReal; // ~1e-14 in double

        resetAllStatistics();
    }

    virtual ~FunctionRep() { }

    virtual String functionKind() const=0; // for error messages

    // Each kind of function can be differentiated with any of these calls, provided
    // the dimensions are reasonable. That is, a scalar function can produce a
    // scalar derivative or a 1x1 Jacobian. If an initial (unperturbed) value is
    // available, its address is passed in fy0. Otherwise, fy0 should be null and
    // we'll attempt to generate our own with an initial call to the user function.
    virtual void calcDerivative(const DifferentiatorRep&, Differentiator::Method,
                                Real y0, const Real* fy0p, Real& dfdy) const=0;
    virtual void calcGradient(const DifferentiatorRep&, Differentiator::Method,
                              const Vector& y0, const Real* fy0p, Vector& gf) const=0; 
    virtual void calcJacobian(const DifferentiatorRep&, Differentiator::Method,
                              const Vector& y0, const Vector* fy0p, Matrix& dfdy) const=0;

    int getNumFunctions()  const {assert(nFunc>=0);  return nFunc;}
    int getNumParameters() const {assert(nParam>=0); return nParam;}
    Real getEstimatedAccuracy() const {
        assert(estimatedAccuracy>0); 
        return estimatedAccuracy;
    }

    void resetAllStatistics() {
        nCalls = nFailures = 0;
    }

protected:
    // Stats
    mutable int nCalls;
    mutable int nFailures;

private:
    int  nFunc, nParam;
    Real estimatedAccuracy;

};

class ScalarFunctionRep : public Differentiator::Function::FunctionRep {
public:
    ScalarFunctionRep(const Differentiator::ScalarFunction& func, Real acc)
    :   FunctionRep(1,1,acc), sf(func) {}

    // Virtuals (from FunctionRep)
    String functionKind() const {return "ScalarFunction";}

    void calcDerivative(const Differentiator::DifferentiatorRep& diff, Differentiator::Method m,
                        Real y0, const Real* fy0p, Real& dfdy) const
    { 
        Real fy0;
        if (fy0p) fy0 = *fy0p;
        else {diff.nCallsToUserFunction++; call(y0, fy0);}
        diff.calcDerivative(*this,m,y0,fy0,dfdy);
    }

    void calcGradient(const Differentiator::DifferentiatorRep& diff, Differentiator::Method m, 
                      const Vector& y0, const Real* fy0p, Vector& gf) const 
    { 
        assert(y0.size()==1);
        gf.resize(1);

        Real fy0;
        if (fy0p) fy0 = *fy0p;
        else {diff.nCallsToUserFunction++; call(y0[0], fy0);}
        diff.calcDerivative(*this,m,y0[0],fy0,gf[0]);
    }

    void calcJacobian(const Differentiator::DifferentiatorRep& diff, Differentiator::Method m, 
                      const Vector& y0, const Vector* fy0p, Matrix& dfdy) const
    {
        assert(y0.size()==1 && (!fy0p || (*fy0p).size()==1));
        dfdy.resize(1,1);

        Real fy0;
        if (fy0p) fy0 = (*fy0p)[0];
        else {diff.nCallsToUserFunction++; call(y0[0], fy0);}
        diff.calcDerivative(*this, m,y0[0],fy0,dfdy(0,0));
    }

    // Local stuff
    void call(const Real& y, Real& fy) const {
        nCalls++;
        nFailures++; // assume failure unless proven otherwise

        int status;
        try 
          { status = sf.f(y,fy);  }
        catch (const std::exception& e)
          { SimTK_THROW1(Differentiator::UserFunctionThrewAnException, e.what()); }
        catch (...)
          { SimTK_THROW1(Differentiator::UserFunctionThrewAnException, 
                         "UNRECOGNIZED EXCEPTION TYPE"); }
        if (status != 0)
            SimTK_THROW1(Differentiator::UserFunctionReturnedNonzeroStatus, status);

        nFailures--;
    }

    const Differentiator::ScalarFunction&       sf;
};

class GradientFunctionRep : public Differentiator::Function::FunctionRep {
public:
    GradientFunctionRep(const Differentiator::GradientFunction& func, int np, Real acc)
    :   FunctionRep(1,np,acc), gf(func) {}

    // Virtuals (from FunctionRep)
    String functionKind() const {return "GradientFunction";}

    void calcDerivative(const Differentiator::DifferentiatorRep& diff, Differentiator::Method m,
                        Real y0, const Real* fy0p, Real& dfdy) const
    {
        if (getNumParameters() != 1)
            SimTK_THROW5(Differentiator::OpNotAllowedForFunctionOfThisShape,
                "calcDerivative", "1x1", "GradientFunction", 1, getNumParameters());

        const Vector y0v(1, y0);
        Real fy0;
        if (fy0p) fy0 = *fy0p; 
        else {diff.nCallsToUserFunction++; call(y0v, fy0);}
        Vector dfdyv(1, &dfdy); // refers to dfdy data directly
        diff.calcGradient(*this,m,y0v,fy0,dfdyv); 
    }

    void calcGradient(const Differentiator::DifferentiatorRep& diff, Differentiator::Method m,
                      const Vector& y0, const Real* fy0p, Vector& gf) const 
    { 
        Real fy0;
        if (fy0p) fy0 = *fy0p; 
        else {diff.nCallsToUserFunction++; call(y0, fy0);}
        diff.calcGradient(*this,m,y0,fy0,gf);
    }

    void calcJacobian(const Differentiator::DifferentiatorRep& diff, Differentiator::Method m,
                      const Vector& y0, const Vector* fy0p, Matrix& dfdy) const
    {
        assert(!fy0p || (*fy0p).size()==1);
        dfdy.resize(1,getNumParameters());
        Real fy0;
        if (fy0p) fy0 = (*fy0p)[0];
        else {diff.nCallsToUserFunction++; call(y0, fy0);}
        diff.calcGradient(*this,m,y0,fy0,~dfdy[0]);
    }

    // Local stuff
    void call(const Vector& y, Real& fy) const {
        nCalls++;
        nFailures++; // assume failure unless proven otherwise

        int status;
        try 
          { status = gf.f(y,fy); } 
        catch (const std::exception& e)
          { SimTK_THROW1(Differentiator::UserFunctionThrewAnException, e.what()); }
        catch (...)
          { SimTK_THROW1(Differentiator::UserFunctionThrewAnException, 
                         "UNRECOGNIZED EXCEPTION TYPE"); }

        if (status != 0)
            SimTK_THROW1(Differentiator::UserFunctionReturnedNonzeroStatus, status);

        --nFailures;
    }

    const Differentiator::GradientFunction&       gf;
};

class JacobianFunctionRep : public Differentiator::Function::FunctionRep {
public:
    JacobianFunctionRep(const Differentiator::JacobianFunction& func, int nf, int np, Real acc)
    :   FunctionRep(nf,np,acc), jf(func) {}

    // Virtuals (from FunctionRep)
    String functionKind() const {return "JacobianFunction";}

    void calcDerivative(const Differentiator::DifferentiatorRep& diff, Differentiator::Method m,
                        Real y0, const Real* fy0p, Real& dfdy) const
    {
        if (getNumFunctions() != 1 || getNumParameters() != 1)
            SimTK_THROW5(Differentiator::OpNotAllowedForFunctionOfThisShape,
                "calcDerivative", "1x1", "JacobianFunction", getNumFunctions(), getNumParameters());

        const Vector y0v(1, y0);
        Vector fy0v(1);
        if (fy0p) fy0v[0] = *fy0p;
        else {diff.nCallsToUserFunction++; call(y0v, fy0v);}
        Matrix dfdym(1,1, &dfdy); // refers to dfdy directly as data
        diff.calcJacobian(*this,m,y0v,fy0v,dfdym); 
    }

    void calcGradient(const Differentiator::DifferentiatorRep& diff, Differentiator::Method m,
                      const Vector& y0, const Real* fy0p, Vector& gf) const 
    { 
        if (getNumFunctions() != 1)
            SimTK_THROW5(Differentiator::OpNotAllowedForFunctionOfThisShape,
                "calcGradient", "1xn", "JacobianFunction", getNumFunctions(), getNumParameters());

        Vector fy0v(1);
        if (fy0p) fy0v[0] = *fy0p;
        else {diff.nCallsToUserFunction++; call(y0, fy0v);}
        diff.calcJacobian(*this,m,y0,fy0v,gf);
    }

    void calcJacobian(const Differentiator::DifferentiatorRep& diff, Differentiator::Method m,
                      const Vector& y0, const Vector* fy0p, Matrix& dfdy) const
    {
        if (fy0p)
            diff.calcJacobian(*this,m,y0,*fy0p,dfdy);
        else {
            Vector fy0(getNumFunctions());
            diff.nCallsToUserFunction++; call(y0, fy0);
            diff.calcJacobian(*this,m,y0,fy0,dfdy);
        }
    }

    void call(const Vector& y, Vector& fy) const {
        nCalls++;
        nFailures++; // assume failure unless proven otherwise

        int status;
        try 
          { status = jf.f(y,fy); } 
        catch (const std::exception& e)
          { SimTK_THROW1(Differentiator::UserFunctionThrewAnException, e.what()); }
        catch (...)
          { SimTK_THROW1(Differentiator::UserFunctionThrewAnException, 
                         "UNRECOGNIZED EXCEPTION TYPE"); }

        if (status != 0)
            SimTK_THROW1(Differentiator::UserFunctionReturnedNonzeroStatus, status);

        nFailures--;
    }

    const Differentiator::JacobianFunction&       jf;
};

    //////////////////////////////////////
    // IMPLEMENTATION OF DIFFERENTIATOR //
    //////////////////////////////////////

Differentiator::Differentiator(const Function& f, Differentiator::Method defaultMethod) {
    rep = new DifferentiatorRep(this, *f.rep, defaultMethod);
}



Differentiator& Differentiator::setDefaultMethod(Differentiator::Method m) {
    rep->defaultMethod = getMethodOrThrow(m, DifferentiatorRep::DefaultDefaultMethod, 
                                          "setDefaultMethod");
    return *this;
}

Differentiator::Method Differentiator::getDefaultMethod() const {
    return rep->defaultMethod;
}

void Differentiator::calcDerivative
   (Real y0, Real fy0, Real& dfdy, Differentiator::Method m) const 
{
    rep->nDifferentiations++;
    rep->nDifferentiationFailures++; // assume the worst

    rep->frep.calcDerivative(*rep,m,y0,&fy0,dfdy);

    rep->nDifferentiationFailures--;
}

// The slow version
Real Differentiator::calcDerivative
   (Real y0, Differentiator::Method m) const 
{
    rep->nDifferentiations++;
    rep->nDifferentiationFailures++; // assume the worst

    Real dfdy;
    rep->frep.calcDerivative(*rep,m,y0,0,dfdy);

    rep->nDifferentiationFailures--;
    return dfdy;
}

void Differentiator::calcGradient
   (const Vector& y0, Real fy0, Vector& gradf,
    Differentiator::Method m) const 
{
    rep->nDifferentiations++;
    rep->nDifferentiationFailures++; // assume the worst

    rep->frep.calcGradient(*rep,m,y0,&fy0,gradf);

    rep->nDifferentiationFailures--;
}

// The slow version
Vector Differentiator::calcGradient
   (const Vector& y0, Differentiator::Method m) const 
{
    rep->nDifferentiations++;
    rep->nDifferentiationFailures++; // assume the worst

    SimTK_APIARGCHECK2_ALWAYS(y0.size()==rep->NParameters, "Differentiator", "calcGradient",
        "Supplied number of parameters (state length) was %ld but should have been %d", 
        y0.size(), rep->NParameters);

    Vector grad(rep->NParameters); 
    rep->frep.calcGradient(*rep,m,y0,0,grad);

    rep->nDifferentiationFailures--;
    return grad; // TODO: use DeadVector to avoid copy
}


void Differentiator::calcJacobian
   (const Vector& y0, const Vector& fy0, Matrix& dfdy,
   Differentiator::Method m) const 
{
    rep->nDifferentiations++;
    rep->nDifferentiationFailures++; // assume the worst

    SimTK_APIARGCHECK2_ALWAYS(y0.size()==rep->NParameters, "Differentiator", "calcJacobian",
        "Expecting %d elements in the parameter (state) vector but got %d", 
        rep->NParameters, (int)y0.size());

    SimTK_APIARGCHECK2_ALWAYS(fy0.size()==rep->NFunctions, "Differentiator", "calcJacobian",
        "Expecting %d elements in the unperturbed function value but got %d", 
        rep->NFunctions, (int)fy0.size());

    rep->frep.calcJacobian(*rep,m,y0,&fy0,dfdy);

    rep->nDifferentiationFailures--;
}

// The slow version
Matrix Differentiator::calcJacobian
   (const Vector& y0, Differentiator::Method m) const 
{
    rep->nDifferentiations++;
    rep->nDifferentiationFailures++; // assume the worst

    SimTK_APIARGCHECK2_ALWAYS(y0.size()==rep->NParameters, "Differentiator", "calcJacobian",
        "Expecting %d elements in the parameter (state) vector but got %d", 
        rep->NParameters, (int)y0.size());

    Matrix dfdy(rep->NFunctions, rep->NParameters); 
    rep->frep.calcJacobian(*rep,m,y0,0,dfdy);

    rep->nDifferentiationFailures--;
    return dfdy; // TODO: use DeadMatrix to avoid copy
}

void Differentiator::resetAllStatistics() {
    rep->resetAllStatistics();
}

int Differentiator::getNumDifferentiations() const {
    return rep->nDifferentiations;
}

int Differentiator::getNumDifferentiationFailures() const {
    return rep->nDifferentiationFailures;
}

int Differentiator::getNumCallsToUserFunction() const {
    return rep->nCallsToUserFunction;
}


/*static*/ bool 
Differentiator::isValidMethod(Differentiator::Method m) {
    switch(m) {
    case UnspecifiedMethod:
    case ForwardDifference:
    case CentralDifference:
        return true;
    }
    return false;
}

/*static*/ const char* 
Differentiator::getMethodName(Differentiator::Method m) {
    throwIfMethodInvalid(m, "getMethodName");
    switch(m) {
    case UnspecifiedMethod: return "UnspecifiedMethod";
    case ForwardDifference: return "ForwardDifference";
    case CentralDifference: return "CentralDifference";
    }
    return 0; // can't happen
}

/*static*/ int 
Differentiator::getMethodOrder(Differentiator::Method m) {
    throwIfMethodInvalid(m, "getMethodOrder");
    if (m==UnspecifiedMethod) 
        m = DifferentiatorRep::DefaultDefaultMethod;
    switch(m) {
    case ForwardDifference: return 1;
    case CentralDifference: return 2;
    }
    assert(!"Shouldn't have gotten here since method was already checked");
    return -1;
}

    ////////////////////////////////////////////////
    // IMPLEMENTATION OF DIFFERENTIATOR::FUNCTION //
    ////////////////////////////////////////////////

// The derived classes ScalarFunction, GradientFunction, JacobianFunction
// are responsible for allocating the right kind of concrete FunctionRep.
Differentiator::Function::Function()
  : rep(0) { }

Differentiator::Function::~Function() {
    delete rep;
}

Differentiator::Function& 
Differentiator::Function::setNumFunctions(int nf) {
    SimTK_APIARGCHECK1_ALWAYS(nf>=0, "Differentiator::Function", "setNumFunctions",
        "The number of functions was %d but must be >= 0", nf);

    rep->nFunc = nf;
    return *this;
}
Differentiator::Function& 
Differentiator::Function::setNumParameters(int ny) {
    SimTK_APIARGCHECK1_ALWAYS(ny>=0, "Differentiator::Function", "setNumParameters",
        "The number of parameters was %d but must be >= 0", ny);

    rep->nParam = ny;
    return *this;
}
Differentiator::Function& 
Differentiator::Function::setEstimatedAccuracy(Real ea) {
    SimTK_APIARGCHECK1_ALWAYS(0<ea&&ea<1, "Differentiator::Function", "setNumParameters",
        "The estimated accuracy was %g but must be between 0 and 1 (noninclusive)", ea);

    rep->estimatedAccuracy = ea;
    return *this;
}

int Differentiator::Function::getNumFunctions() const {
    return rep->nFunc;
}
int Differentiator::Function::getNumParameters() const {
    return rep->nParam;
}
Real Differentiator::Function::getEstimatedAccuracy() const {
    return rep->estimatedAccuracy;
}

void Differentiator::Function::resetAllStatistics(){
    rep->resetAllStatistics();
}
int Differentiator::Function::getNumCalls()    const {
    return rep->nCalls;
}
int Differentiator::Function::getNumFailures() const {
    return rep->nFailures;
}

    ///////////////////////////////////////////////////////////////////
    // IMPLEMENTATION OF SCALAR, GRADIENT, JACOBIAN FUNCTION CLASSES //
    ///////////////////////////////////////////////////////////////////

Differentiator::ScalarFunction::ScalarFunction(Real acc) {
    rep = new ScalarFunctionRep(*this, acc);
}
Differentiator::GradientFunction::GradientFunction(int np, Real acc) {
    rep = new GradientFunctionRep(*this, np, acc);
}
Differentiator::JacobianFunction::JacobianFunction(int nf, int np, Real acc) {
    rep = new JacobianFunctionRep(*this, nf, np, acc);
}


    //////////////////////////////////////////
    // IMPLEMENTATION OF DIFFERENTIATOR REP //
    //////////////////////////////////////////

Differentiator::DifferentiatorRep::DifferentiatorRep
(   Differentiator*                                 handle,
    const Differentiator::Function::FunctionRep&    fr,
    Differentiator::Method                          defMthd
)
:   myHandle(handle), frep(fr),
    NParameters(fr.getNumParameters()), 
    NFunctions(fr.getNumFunctions()), 
    EstimatedAccuracy(fr.getEstimatedAccuracy()),
    defaultMethod(getMethodOrThrow(defMthd, DefaultDefaultMethod, "Differentiator")),
    AccFac1(std::sqrt(EstimatedAccuracy)),
    AccFac2(std::pow(EstimatedAccuracy, OneThird))
{
    //TODO
    assert(NParameters >= 0 && NFunctions >= 0 && EstimatedAccuracy > 0);

    resetAllStatistics();
    ytmp.resize(NParameters);
    fyptmp.resize(NFunctions);
    fymtmp.resize(NFunctions);
}

void Differentiator::DifferentiatorRep::calcDerivative
   (const ScalarFunctionRep& f, Differentiator::Method m, Real y0, Real fy0, Real& dfdy) const 
{
    // This won't return if the method is bad.
    const Differentiator::Method method = getMethodOrThrow(m, defaultMethod, "calcDerivative");

    //TODO
    assert(NParameters==1 && NFunctions==1);

    const int  order = Differentiator::getMethodOrder(method);
    const Real hEst  = getAccFac(order)*std::max(std::abs(y0), YMin);
    const Real h     = cleanUpH(hEst, y0);

    Real fyplus, fyminus;
    Real y = y0+h; 
    nCallsToUserFunction++; f.call(y, fyplus);
    if (order==1) {
        dfdy = (fyplus-fy0)/h;
    } else {
        y = y0-h; 
        nCallsToUserFunction++; f.call(y, fyminus);
        dfdy = (fyplus-fyminus)/(2*h);
    }
}

void Differentiator::DifferentiatorRep::calcGradient
   (const GradientFunctionRep& f, Differentiator::Method m, const Vector& y0, Real fy0, Vector& gradf) const 
{
    // This won't return if the method is bad.
    const Differentiator::Method method = getMethodOrThrow(m, defaultMethod, "calcGradient");

    //TODO
    assert(NFunctions==1);
    assert(ytmp.size()==NParameters);
    assert(y0.size() == NParameters);

    gradf.resize(NParameters);

    ytmp = y0;
    const int order = Differentiator::getMethodOrder(method);

    for (int i=0; i < f.getNumParameters(); ++i) {
        const Real hEst = getAccFac(order)*std::max(std::abs(y0[i]), YMin);
        const Real h = cleanUpH(hEst, y0[i]);
        Real fyplus, fyminus;
        ytmp[i] = y0[i]+h; 
        nCallsToUserFunction++; f.call(ytmp, fyplus);
        if (order==1) {
            gradf[i] = (fyplus-fy0)/h;
        } else {
            ytmp[i] = y0[i]-h; 
            nCallsToUserFunction++; f.call(ytmp, fyminus);
            gradf[i] = (fyplus-fyminus)/(2*h);
        }
        ytmp[i] = y0[i]; // restore
    }
}

void Differentiator::DifferentiatorRep::calcJacobian
   (const JacobianFunctionRep& f, Differentiator::Method m, 
    const Vector& y0, const Vector& fy0, Matrix& dfdy) const 
{
    // This won't return if the method is bad.
    const Differentiator::Method method = getMethodOrThrow(m, defaultMethod, "calcJacobian");

    //TODO
    assert(ytmp.size()==NParameters && fyptmp.size()==NFunctions && fymtmp.size()==NFunctions);
    assert(y0.size()  == NParameters);
    assert(fy0.size() == NFunctions);

    dfdy.resize(NFunctions,NParameters);

    const int order = Differentiator::getMethodOrder(method);

    ytmp = y0;
    for (int i=0; i < NParameters; ++i) {
        const Real hEst = getAccFac(order)*std::max(std::abs(y0[i]), YMin);
        const Real h = cleanUpH(hEst, y0[i]);
        ytmp[i] = y0[i]+h; 
        nCallsToUserFunction++; f.call(ytmp, fyptmp);
        if (order==1) {
            dfdy(i) = (fyptmp-fy0)/h;
        } else {
            ytmp[i] = y0[i]-h; 
            nCallsToUserFunction++; f.call(ytmp, fymtmp);
            dfdy(i) = (fyptmp-fymtmp)/(2*h);
        }
        ytmp[i] = y0[i]; // restore
    }
}

} // namespace SimTK


