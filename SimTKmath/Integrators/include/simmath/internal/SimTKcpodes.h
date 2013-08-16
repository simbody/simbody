#ifndef SimTK_SIMMATH_CPODES_H_
#define SimTK_SIMMATH_CPODES_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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

/** @file
 * This is the header file that user code should include to pick up the 
 * SimTK C++ interface to the Sundials CPODES coordinate-projection 
 * integrator.
 */

#include "SimTKcommon.h"
#include "simmath/internal/common.h"

#include <cstdio> // Needed for "FILE".

namespace SimTK {

/**
 * This abstract class defines the system to be integrated with SimTK
 * CPodes. Note that this defines a client-side virtual function table
 * which must be used only on the client side. Library-side access to
 * these virtual functions is done only through a set
 * of equivalent static functions (provided in this same header files)
 * whose addresses can be reliably "tossed over the fence" to the library
 * side without compromising binary compatibility.
 */
class SimTK_SIMMATH_EXPORT CPodesSystem {
public:
    virtual ~CPodesSystem() {}
    // The default implementations of these virtual functions
    // just throw an "unimplemented virtual function" exception. 

    // At least one of these two must be supplied by the concrete class.
    virtual int  explicitODE(Real t, const Vector& y, Vector& fout) const;
    virtual int  implicitODE(Real t, const Vector& y, const Vector& yp, 
                             Vector& fout) const;

    virtual int  constraint(Real t, const Vector& y, 
                            Vector& cout) const;
    virtual int  project(Real t, const Vector& ycur, Vector& corr,
                         Real epsProj, Vector& err) const; // err is in/out
    virtual int  quadrature(Real t, const Vector& y, 
                            Vector& qout) const;
    virtual int  root(Real t, const Vector& y, const Vector& yp,
                      Vector& gout) const;
    virtual int  weight(const Vector& y, Vector& weights) const;
    virtual void errorHandler(int error_code, const char* module,
                              const char* function, char* msg) const;

    //TODO: Jacobian functions
};


// These static functions are private to the current (client-side) compilation
// unit. They are used to navigate the client-side CPodesSystem virtual function
// table, which cannot be done on the library side. Note that these are defined
// in the SimTK namespace so don't need "SimTK" in their names.
static int explicitODE_static(const CPodesSystem& sys, 
                              Real t, const Vector& y, Vector& fout) 
  { return sys.explicitODE(t,y,fout); }

static int implicitODE_static(const CPodesSystem& sys, 
                              Real t, const Vector& y, const Vector& yp, Vector& fout)
  { return sys.implicitODE(t,y,yp,fout); }

static int constraint_static(const CPodesSystem& sys, 
                             Real t, const Vector& y, Vector& cout)
  { return sys.constraint(t,y,cout); }

static int project_static(const CPodesSystem& sys, 
                          Real t, const Vector& ycur, Vector& corr,
                          Real epsProj, Vector& err)
  { return sys.project(t,ycur,corr,epsProj,err); }

static int quadrature_static(const CPodesSystem& sys, 
                             Real t, const Vector& y, Vector& qout)
  { return sys.quadrature(t,y,qout); }

static int root_static(const CPodesSystem& sys, 
                       Real t, const Vector& y, const Vector& yp, Vector& gout)
  { return sys.root(t,y,yp,gout); }

static int weight_static(const CPodesSystem& sys, 
                         const Vector& y, Vector& weights)
  { return sys.weight(y,weights); }

static void errorHandler_static(const CPodesSystem& sys, 
                                int error_code, const char* module, 
                                const char* function, char* msg)
  { sys.errorHandler(error_code,module,function,msg); }

/**
 * This is a straightforward translation of the Sundials CPODES C 
 * interface into C++. The class CPodes represents a single instance
 * of a CPODES integrator, and handles the associated memory internally.
 * Methods here are identical to the corresponding CPODES functions (with
 * the "CPode" prefix removed) but are const-correct and use SimTK
 * Vector & Real rather than Sundials N_Vector and realtype.
 */
class SimTK_SIMMATH_EXPORT CPodes {
public:
    // no default constructor
    // copy constructor and default assignment are suppressed

    enum ODEType {
        UnspecifiedODEType=0,
        ExplicitODE,
        ImplicitODE
    };

    enum LinearMultistepMethod {
        UnspecifiedLinearMultistepMethod=0,
        BDF,
        Adams
    };

    enum NonlinearSystemIterationType {
        UnspecifiedNonlinearSystemIterationType=0,
        Newton,
        Functional
    };

    enum ToleranceType {
        UnspecifiedToleranceType=0,
        ScalarScalar,
        ScalarVector,
        WeightFunction
    };

    enum ProjectionNorm {
        UnspecifiedProjectionNorm=0,
        L2Norm,
        ErrorNorm
    };

    enum ConstraintLinearity {
        UnspecifiedConstraintLinearity=0,
        Linear,
        Nonlinear
    };

    enum ProjectionFactorizationType {
        UnspecifiedProjectionFactorizationType=0,
        ProjectWithLU,
        ProjectWithQR,
        ProjectWithSchurComplement,
        ProjectWithQRPivot  // for handling redundancy
    };

    enum StepMode {
        UnspecifiedStepMode=0,
        Normal,
        OneStep,
        NormalTstop,
        OneStepTstop
    };

    explicit CPodes
       (ODEType                      ode=UnspecifiedODEType, 
        LinearMultistepMethod        lmm=UnspecifiedLinearMultistepMethod, 
        NonlinearSystemIterationType nls=UnspecifiedNonlinearSystemIterationType)
    {
        // Perform construction of the CPodesRep on the library side.
        librarySideCPodesConstructor(ode, lmm, nls);
        // But fill in function pointers from the client side.
        clientSideCPodesConstructor();
    }

    // Values for 'flag' return values. These are just the "normal" return
    // values; there are many more which are all negative and represent
    // error conditions.
    static const int Success     = 0;
    static const int TstopReturn = 1;
    static const int RootReturn  = 2;
    static const int Warning     = 99;
    static const int TooMuchWork = -1;
    static const int TooClose    = -27;

    // These values should be used by user routines. "Success" is the
    // same as above. A positive return value means "recoverable error",
    // i.e., CPodes should cut the step size and try again, while a
    // negative number means "unrecoverable error" which will kill
    // CPODES altogether with a CP_xxx_FAIL error. The particular numerical
    // values here have no significance, just + vs. -.
    static const int RecoverableError = 9999;
    static const int UnrecoverableError = -9999;

    ~CPodes();

    // Depending on the setting of ode_type at construction, init()
    // and reInit() will tell CPodes to use either the explicitODE()
    // or implicitODE() function from the CPodesSystem, so the user
    // MUST have overridden at least one of those virtual methods.
    int init(CPodesSystem& sys,
             Real t0, const Vector& y0, const Vector& yp0,
             ToleranceType tt, Real reltol, void* abstol);

    int reInit(CPodesSystem& sys,
               Real t0, const Vector& y0, const Vector& yp0,
               ToleranceType tt, Real reltol, void* abstol);

    // This tells CPodes to make use of the user's constraint()
    // method from CPodesSystem, and perform projection internally.
    int projInit(ProjectionNorm, ConstraintLinearity,
                 const Vector& ctol);

    // This tells CPodes to make use of the user's project()
    // method from CPodesSystem.
    int projDefine();

    // These tell CPodes to make use of the user's quadrature()
    // method from CPodesSystem.
    int quadInit(const Vector& q0);
    int quadReInit(const Vector& q0);

    // This tells CPodes to make use of the user's root() method
    // from CPodesSystem.
    int rootInit(int nrtfn);

    // This tells CPodes to make use of the user's errorHandler() 
    // method from CPodesSystem.
    int setErrHandlerFn();

    // These tells CPodes to make use of the user's weight() 
    // method from CPodesSystem.
    int setEwtFn();

    // TODO: these routines should enable methods that are defined
    // in the CPodesSystem, but a proper interface to the Jacobian
    // routines hasn't been implemented yet.
    int dlsSetJacFn(void* jac, void* jac_data);
    int dlsProjSetJacFn(void* jacP, void* jacP_data);


    int step(Real tout, Real* tret, 
             Vector& y_inout, Vector& yp_inout, StepMode=Normal);

    int setErrFile(FILE* errfp);
    int setMaxOrd(int maxord);
    int setMaxNumSteps(int mxsteps);
    int setMaxHnilWarns(int mxhnil);
    int setStabLimDet(bool stldet) ;
    int setInitStep(Real hin);
    int setMinStep(Real hmin);
    int setMaxStep(Real hmax);
    int setStopTime(Real tstop);
    int setMaxErrTestFails(int maxnef);

    int setMaxNonlinIters(int maxcor);
    int setMaxConvFails(int maxncf);
    int setNonlinConvCoef(Real nlscoef);

    int setProjUpdateErrEst(bool proj_err);
    int setProjFrequency(int proj_freq);
    int setProjTestCnstr(bool test_cnstr);
    int setProjLsetupFreq(int proj_lset_freq);
    int setProjNonlinConvCoef(Real prjcoef);

    int setQuadErrCon(bool errconQ, 
                      int tol_typeQ, Real reltolQ, void* abstolQ);

    int setTolerances(int tol_type, Real reltol, void* abstol);
    
    int setRootDirection(Array_<int>& rootdir);

    int getDky(Real t, int k, Vector& dky);

    int getQuad(Real t, Vector& yQout);
    int getQuadDky(Real t, int k, Vector& dky);

    int getWorkSpace(int* lenrw, int* leniw);
    int getNumSteps(int* nsteps);
    int getNumFctEvals(int* nfevals);
    int getNumLinSolvSetups(int* nlinsetups);
    int getNumErrTestFails(int* netfails);
    int getLastOrder(int* qlast);
    int getCurrentOrder(int* qcur);
    int getNumStabLimOrderReds(int* nslred);
    int getActualInitStep(Real* hinused);
    int getLastStep(Real* hlast);
    int getCurrentStep(Real* hcur);
    int getCurrentTime(Real* tcur);
    int getTolScaleFactor(Real* tolsfac);
    int getErrWeights(Vector& eweight);
    int getEstLocalErrors(Vector& ele) ;
    int getNumGEvals(int* ngevals);
    int getRootInfo(int* rootsfound);
    int getRootWindow(Real* tLo, Real* tHi);
    int getIntegratorStats(int* nsteps,
                           int* nfevals, int* nlinsetups,
                           int* netfails, int* qlast,
                           int* qcur, Real* hinused, Real* hlast,
                           Real* hcur, Real* tcur);

    int getNumNonlinSolvIters(int* nniters);
    int getNumNonlinSolvConvFails(int* nncfails);
    int getNonlinSolvStats(int* nniters, int* nncfails);
    int getProjNumProj(int* nproj);
    int getProjNumCnstrEvals(int* nce);
    int getProjNumLinSolvSetups(int* nsetupsP);
    int getProjNumFailures(int* nprf) ;
    int getProjStats(int* nproj,
                     int* nce, int* nsetupsP,
                     int* nprf);
    int getQuadNumFunEvals(int* nqevals);
    int getQuadErrWeights(Vector& eQweight);
    char* getReturnFlagName(int flag);


    int dlsGetWorkSpace(int* lenrwLS, int* leniwLS);
    int dlsGetNumJacEvals(int* njevals);
    int dlsGetNumFctEvals(int* nfevalsLS);
    int dlsGetLastFlag(int* flag);
    char* dlsGetReturnFlagName(int flag);

    int dlsProjGetNumJacEvals(int* njPevals);
    int dlsProjGetNumFctEvals(int* ncevalsLS);

    int lapackDense(int N);
    int lapackBand(int N, int mupper, int mlower);
    int lapackDenseProj(int Nc, int Ny, ProjectionFactorizationType);

private:
    // This is how we get the client-side virtual functions to
    // be callable from library-side code while maintaining binary
    // compatibility.
    typedef int (*ExplicitODEFunc)(const CPodesSystem&, 
                                   Real t, const Vector& y, Vector& fout);
    typedef int (*ImplicitODEFunc)(const CPodesSystem&, 
                                   Real t, const Vector& y, const Vector& yp,
                                   Vector& fout);
    typedef int (*ConstraintFunc) (const CPodesSystem&, 
                                   Real t, const Vector& y, Vector& cout);
    typedef int (*ProjectFunc)    (const CPodesSystem&, 
                                   Real t, const Vector& ycur, Vector& corr,
                                   Real epsProj, Vector& err);
    typedef int (*QuadratureFunc) (const CPodesSystem&, 
                                   Real t, const Vector& y, Vector& qout);
    typedef int (*RootFunc)       (const CPodesSystem&, 
                                   Real t, const Vector& y, const Vector& yp, 
                                   Vector& gout);
    typedef int (*WeightFunc)     (const CPodesSystem&, 
                                   const Vector& y, Vector& weights);
    typedef void (*ErrorHandlerFunc)(const CPodesSystem&, 
                                     int error_code, const char* module, 
                                     const char* function, char* msg);

    // Note that these routines do not tell CPodes to use the supplied
    // functions. They merely provide the client-side addresses of functions
    // which understand how to find the user's virtual functions, should those
    // actually be provided. Control over whether to actually call any of these
    // is handled elsewhere, with user-visible methods. These private methods
    // are to be called only upon construction of the CPodes object here. They
    // are not even dependent on which user-supplied concrete CPodesSystem is
    // being used.
    void registerExplicitODEFunc(ExplicitODEFunc);
    void registerImplicitODEFunc(ImplicitODEFunc);
    void registerConstraintFunc(ConstraintFunc);
    void registerProjectFunc(ProjectFunc);
    void registerQuadratureFunc(QuadratureFunc);
    void registerRootFunc(RootFunc);
    void registerWeightFunc(WeightFunc);
    void registerErrorHandlerFunc(ErrorHandlerFunc);


    // This is the library-side part of the CPodes constructor. This must
    // be done prior to the client side construction.
    void librarySideCPodesConstructor(ODEType, LinearMultistepMethod, NonlinearSystemIterationType);

    // Note that this routine MUST be called from client-side code so that
    // it picks up exactly the static routines above which will agree with
    // the client about the layout of the CPodesSystem virtual function table.
    void clientSideCPodesConstructor() {
        registerExplicitODEFunc(explicitODE_static);
        registerImplicitODEFunc(implicitODE_static);
        registerConstraintFunc(constraint_static);
        registerProjectFunc(project_static);
        registerQuadratureFunc(quadrature_static);
        registerRootFunc(root_static);
        registerWeightFunc(weight_static);
        registerErrorHandlerFunc(errorHandler_static);
    }

    // FOR INTERNAL USE ONLY
private:
    class CPodesRep* rep;
    friend class CPodesRep;

    const CPodesRep& getRep() const {assert(rep); return *rep;}
    CPodesRep&       updRep()       {assert(rep); return *rep;}

    // Suppress copy constructor and default assigment operator.
    CPodes(const CPodes&);
    CPodes& operator=(const CPodes&);
};

} // namespace SimTK

#endif // SimTK_CPODES_H_
