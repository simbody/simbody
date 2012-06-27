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
 * This is the private (library side) implementation of the CPodes
 * handle class and its actual representation CPodesRep.
 */

#include "simmath/internal/SimTKcpodes.h"
#include "nvector_SimTK.h"
#include "cpodes/cpodes.h"
#include "cpodes/cpodes_dense.h"
#include "cpodes/cpodes_lapack_exports.h"

#include <limits>

namespace SimTK {

class CPodesRep {
public:
    CPodesRep()
      : useImplicitODEFunction(false), cpode_mem(0), sysp(0), myHandle(0) 
    { 
        zeroFunctionPointers();
    }

    CPodesRep(int ode_type, int lmm_type, int nls_type)
      : useImplicitODEFunction(ode_type == CP_IMPL), cpode_mem(0), sysp(0), myHandle(0)
    {
        cpode_mem = CPodeCreate(ode_type, lmm_type, nls_type);
    }

    ~CPodesRep() {
        if (cpode_mem)
            CPodeFree(&cpode_mem);
    }

    const CPodesSystem& getCPodesSystem() const {return *sysp;}

    // Client-side function pointers
    CPodes::ExplicitODEFunc     explicitODEFunc;
    CPodes::ImplicitODEFunc     implicitODEFunc;
    CPodes::ConstraintFunc      constraintFunc;
    CPodes::ProjectFunc         projectFunc;
    CPodes::QuadratureFunc      quadratureFunc;
    CPodes::RootFunc            rootFunc;
    CPodes::WeightFunc          weightFunc;
    CPodes::ErrorHandlerFunc    errorHandlerFunc;

    void zeroFunctionPointers() {
        explicitODEFunc  = 0;
        implicitODEFunc  = 0;
        constraintFunc   = 0;
        projectFunc      = 0;
        quadratureFunc   = 0;
        rootFunc         = 0;
        weightFunc       = 0;
        errorHandlerFunc = 0;
    }

    void setMyHandle(CPodes& cp) {myHandle = &cp;}
    const CPodes& getMyHandle() const {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}
private:
    bool  useImplicitODEFunction;
    void* cpode_mem;
    const CPodesSystem* sysp;

    friend class CPodes;
    CPodes* myHandle;   // The owner handle of this Rep.

    void* getSystemForSundials() const {
        return (void*)const_cast<CPodesSystem*>(sysp);
    }
};


//////////////////////
// CPODES FUNCTIONS //
//////////////////////

// These functions satisfy the C signatures required by Sundials CPODES. They provide
// translation between SimTK data types and those understood by Sundials and CPODES.
// These in turn call (indirectly) user-supplied C++ code that was used in
// defining concrete classes derived from the user-visible CPodesSystem abstract
// class. Note that we *never* access the client-side CPodesSystem virtual
// function table from library code. CPodesRep contains an alternate table which
// will be correct even if the client side is out of date with respect to the
// library.

static int explicitODEWrapper(realtype t, N_Vector nv_y, N_Vector nv_fout, void* f_data) {
    const Vector& y    = N_Vector_SimTK::getVector(nv_y);
    Vector&       fout = N_Vector_SimTK::updVector(nv_fout);
    const CPodesRep& rep = *reinterpret_cast<const CPodesRep*>(f_data);
    return rep.explicitODEFunc(rep.getCPodesSystem(), t, y, fout);
}

static int implicitODEWrapper(realtype t, N_Vector nv_y, N_Vector nv_yp, N_Vector nv_fout, void* f_data) {
    const Vector& y    = N_Vector_SimTK::getVector(nv_y);
    const Vector& yp   = N_Vector_SimTK::getVector(nv_yp);
    Vector&       fout = N_Vector_SimTK::updVector(nv_fout);
    const CPodesRep& rep = *reinterpret_cast<const CPodesRep*>(f_data);
    return rep.implicitODEFunc(rep.getCPodesSystem(), t, y, yp, fout);
}

static int constraintWrapper(realtype t, N_Vector nv_y, N_Vector nv_cout, void* c_data) {
    const Vector& y    = N_Vector_SimTK::getVector(nv_y);
    Vector&       cout = N_Vector_SimTK::updVector(nv_cout);
    const CPodesRep& rep = *reinterpret_cast<const CPodesRep*>(c_data);
    return rep.constraintFunc(rep.getCPodesSystem(), t, y, cout);
}

static int projectWrapper(realtype t, N_Vector nv_ycur,
                          N_Vector nv_corr,
                          realtype epsProj, N_Vector nv_err, void *pdata)
{
    const Vector& ycur   = N_Vector_SimTK::getVector(nv_ycur);
    Vector&       corr   = N_Vector_SimTK::updVector(nv_corr);
    Vector&       err    = N_Vector_SimTK::updVector(nv_err);
    const CPodesRep& rep = *reinterpret_cast<const CPodesRep*>(pdata);
    return rep.projectFunc(rep.getCPodesSystem(), t, ycur, corr, epsProj, err);
}

static int quadratureWrapper(realtype t, N_Vector nv_y, 
                             N_Vector nv_qout, void *q_data) 
{
    const Vector& y    = N_Vector_SimTK::getVector(nv_y);
    Vector&       qout = N_Vector_SimTK::updVector(nv_qout);
    const CPodesRep& rep = *reinterpret_cast<const CPodesRep*>(q_data);
    return rep.quadratureFunc(rep.getCPodesSystem(), t, y, qout);
}

static int rootWrapper(realtype t, N_Vector nv_y, N_Vector nv_yp,
                       realtype *goutp, void *g_data) 
{
    const Vector& y    = N_Vector_SimTK::getVector(nv_y);
    const Vector& yp   = N_Vector_SimTK::getVector(nv_yp);

    const CPodesRep& rep = *reinterpret_cast<const CPodesRep*>(g_data);

    // Create a temporary Vector to hold the result. This is awkward
    // but we don't know here how much space to allocate since for
    // some odd reason the output array here isn't an N_Vector.
    Vector gout_tmp; // the root routine will have to resize this
    const int flag = rep.rootFunc(rep.getCPodesSystem(), t, y, yp, gout_tmp);
    for (int i=0; i<gout_tmp.size(); ++i)
        goutp[i] = gout_tmp[i];
    return flag;
}

static int weightWrapper(N_Vector nv_y, N_Vector nv_ewt, void *e_data) {
    const Vector& y    = N_Vector_SimTK::getVector(nv_y);
    Vector&       ewt  = N_Vector_SimTK::updVector(nv_ewt);
    const CPodesRep& rep = *reinterpret_cast<const CPodesRep*>(e_data);
    return rep.weightFunc(rep.getCPodesSystem(), y, ewt);
}

static void errorHandlerWrapper(int error_code, 
                                const char *module, const char *function, 
                                char *msg, void *eh_data)
{
    const CPodesRep& rep = *reinterpret_cast<const CPodesRep*>(eh_data);
    return rep.errorHandlerFunc(rep.getCPodesSystem(), error_code,module,function,msg);
}

////////////////////////////////////////
// CLASS SimTK::CPodes IMPLEMENTATION //
////////////////////////////////////////

// These static functions map SimTK::CPodes enumerated types to C CPodes #defines.
static int mapODEType(CPodes::ODEType ode) {
    switch(ode) {
    case CPodes::ExplicitODE: return CP_EXPL;
    case CPodes::ImplicitODE: return CP_IMPL;
    default: return std::numeric_limits<int>::min();
    }
}

static int mapLinearMultistepMethod(CPodes::LinearMultistepMethod lmm) {
    switch(lmm) {
    case CPodes::BDF:   return CP_BDF;
    case CPodes::Adams: return CP_ADAMS;
    default: return std::numeric_limits<int>::min();
    }
}

static int mapNonlinearSystemIterationType(CPodes::NonlinearSystemIterationType nls) {
    switch(nls) {
    case CPodes::Functional: return CP_FUNCTIONAL;
    case CPodes::Newton:     return CP_NEWTON;
    default: return std::numeric_limits<int>::min();
    }
}

static int mapToleranceType(CPodes::ToleranceType type) {
    switch(type) {
    case CPodes::ScalarScalar:    return CP_SS;
    case CPodes::ScalarVector:    return CP_SV;
    case CPodes::WeightFunction:  return CP_WF;
    default: return std::numeric_limits<int>::min();
    }
}

static int mapProjectionNorm(CPodes::ProjectionNorm norm) {
    switch(norm) {
    case CPodes::L2Norm:    return CP_PROJ_L2NORM;
    case CPodes::ErrorNorm: return CP_PROJ_ERRNORM;
    default: return std::numeric_limits<int>::min();
    }
}

static int mapConstraintLinearity(CPodes::ConstraintLinearity lin) {
    switch(lin) {
    case CPodes::Linear:    return CP_CNSTR_LIN;
    case CPodes::Nonlinear: return CP_CNSTR_NONLIN;
    default: return std::numeric_limits<int>::min();
    }
}


static int mapProjectionFactorizationType(CPodes::ProjectionFactorizationType ft) {
    switch(ft) {
    case CPodes::ProjectWithLU:               return CPDIRECT_LU;
    case CPodes::ProjectWithQR:               return CPDIRECT_QR;
    case CPodes::ProjectWithSchurComplement:  return CPDIRECT_SC;
    case CPodes::ProjectWithQRPivot:          return CPDIRECT_QRP;
    default: return std::numeric_limits<int>::min();
    }
}

static int mapStepMode(CPodes::StepMode mode) {
    switch(mode) {
    case CPodes::Normal:       return CP_NORMAL;
    case CPodes::OneStep:      return CP_ONE_STEP;
    case CPodes::NormalTstop:  return CP_NORMAL_TSTOP;
    case CPodes::OneStepTstop: return CP_ONE_STEP_TSTOP;
    default: return std::numeric_limits<int>::min();
    }
}

// The actual constructor is defined on the client side but
// calls this library-side routine to do most of the work. (Registration of
// user functions has to be done on the client side.)
void CPodes::librarySideCPodesConstructor
   (ODEType ode, LinearMultistepMethod lmm, NonlinearSystemIterationType nls) 
{
    if (ode == UnspecifiedODEType) ode = ExplicitODE;
    if (lmm == UnspecifiedLinearMultistepMethod) lmm = BDF;
    if (nls == UnspecifiedNonlinearSystemIterationType) 
        nls = (ode==ExplicitODE ? Functional : Newton);

    rep = new CPodesRep(mapODEType(ode), 
                        mapLinearMultistepMethod(lmm), 
                        mapNonlinearSystemIterationType(nls));
    rep->setMyHandle(*this);
}

CPodes::~CPodes() {
    if (rep && rep->myHandle==this) // this is the owner handle
        delete rep;
    rep = 0;
}

int CPodes::init(CPodesSystem& sys,
                 Real t0, const Vector& y0, const Vector& yp0,
                 ToleranceType tol_type, Real reltol, void* abstol)
{
    updRep().sysp = &sys;

    N_Vector_SimTK nv_y0(y0);   // references, not copies
    N_Vector_SimTK nv_yp0(yp0);
    return CPodeInit(updRep().cpode_mem,
                     getRep().useImplicitODEFunction 
                        ? (void*)implicitODEWrapper
                        : (void*)explicitODEWrapper,
                     (void*)rep, 
                     t0, &nv_y0, &nv_yp0,
                     mapToleranceType(tol_type),reltol,abstol);
}

int CPodes::reInit(CPodesSystem& sys,
                   Real t0, const Vector& y0, const Vector& yp0,
                   ToleranceType tol_type, Real reltol, void* abstol)
{
    updRep().sysp = &sys;

    N_Vector_SimTK nv_y0(y0);   // references, not copies
    N_Vector_SimTK nv_yp0(yp0);
    return CPodeReInit(updRep().cpode_mem,
                       getRep().useImplicitODEFunction 
                          ? (void*)implicitODEWrapper
                          : (void*)explicitODEWrapper,
                       (void*)rep,
                       t0, &nv_y0, &nv_yp0,
                       mapToleranceType(tol_type),reltol,abstol);
}

int CPodes::projInit(ProjectionNorm norm, ConstraintLinearity lin,
                     const Vector& ctol)
{
    N_Vector_SimTK nv_ctol(ctol);
    return CPodeProjInit(updRep().cpode_mem, 
                         mapProjectionNorm(norm), 
                         mapConstraintLinearity(lin),
                         constraintWrapper, (void*)rep,
                         &nv_ctol);
}

int CPodes::projDefine() {
    return CPodeProjDefine(updRep().cpode_mem, projectWrapper, (void*)rep);
}

int CPodes::quadInit(const Vector& q0) {
    N_Vector_SimTK nv_q0(q0);
    return CPodeQuadInit(updRep().cpode_mem, quadratureWrapper, (void*)rep,
                         &nv_q0);
}

int CPodes::quadReInit(const Vector& q0) {
    N_Vector_SimTK nv_q0(q0);
    return CPodeQuadReInit(updRep().cpode_mem, quadratureWrapper, (void*)rep,
                           &nv_q0);
}

int CPodes::rootInit(int nrtfn) {
    return CPodeRootInit(updRep().cpode_mem, nrtfn, rootWrapper, (void*)rep);
}

int CPodes::setErrHandlerFn() {
    return CPodeSetErrHandlerFn(updRep().cpode_mem, errorHandlerWrapper, (void*)rep);
}
int CPodes::setErrFile(FILE* errfp) {
    return CPodeSetErrFile(updRep().cpode_mem,errfp);
}

int CPodes::setEwtFn() {
    return CPodeSetEwtFn(updRep().cpode_mem, weightWrapper, (void*)rep);
}
int CPodes::setMaxOrd(int maxord) {
    return CPodeSetMaxOrd(updRep().cpode_mem,maxord);
}
int CPodes::setMaxNumSteps(int mxsteps) {
    return CPodeSetMaxNumSteps(updRep().cpode_mem,(long)mxsteps);
}
int CPodes::setMaxHnilWarns(int mxhnil) {
    return CPodeSetMaxHnilWarns(updRep().cpode_mem,mxhnil);
}
int CPodes::setStabLimDet(bool stldet) {
    return CPodeSetStabLimDet(updRep().cpode_mem,(booleantype)stldet);
}
int CPodes::setInitStep(Real hin) {
    return CPodeSetInitStep(updRep().cpode_mem,hin);
}
int CPodes::setMinStep(Real hmin) {
    return CPodeSetMinStep(updRep().cpode_mem,hmin);
}
int CPodes::setMaxStep(Real hmax) {
    return CPodeSetMaxStep(updRep().cpode_mem,hmax);
}
int CPodes::setStopTime(Real tstop) {
    return CPodeSetStopTime(updRep().cpode_mem,tstop);
}
int CPodes::setMaxErrTestFails(int maxnef) {
    return CPodeSetMaxErrTestFails(updRep().cpode_mem,maxnef);
}

int CPodes::setMaxNonlinIters(int maxcor) {
    return CPodeSetMaxNonlinIters(updRep().cpode_mem,maxcor);
}
int CPodes::setMaxConvFails(int maxncf) {
    return CPodeSetMaxConvFails(updRep().cpode_mem,maxncf);
}
int CPodes::setNonlinConvCoef(Real nlscoef) {
    return CPodeSetNonlinConvCoef(updRep().cpode_mem,nlscoef);
}

int CPodes::setProjUpdateErrEst(bool proj_err) {
    return CPodeSetProjUpdateErrEst(updRep().cpode_mem,(booleantype)proj_err);
}
int CPodes::setProjFrequency(int proj_freq) {
    return CPodeSetProjFrequency(updRep().cpode_mem,(long)proj_freq);
}
int CPodes::setProjTestCnstr(bool test_cnstr) {
    return CPodeSetProjTestCnstr(updRep().cpode_mem,(booleantype)test_cnstr);
}
int CPodes::setProjLsetupFreq(int proj_lset_freq) {
    return CPodeSetProjLsetupFreq(updRep().cpode_mem,(long)proj_lset_freq);
}
int CPodes::setProjNonlinConvCoef(Real prjcoef) {
    return CPodeSetProjNonlinConvCoef(updRep().cpode_mem,prjcoef);
}

int CPodes::setQuadErrCon(bool errconQ, 
                  int tol_typeQ, Real reltolQ, void* abstolQ) {
    return CPodeSetQuadErrCon(updRep().cpode_mem,(booleantype)errconQ,
                              tol_typeQ,reltolQ,abstolQ);
}

int CPodes::setTolerances(int tol_type, Real reltol, void* abstol) {
    return CPodeSetTolerances(updRep().cpode_mem,tol_type,reltol,abstol);
}

int CPodes::setRootDirection(Array_<int>& rootdir) {
    int* array = new int[rootdir.size()];
    for (int i = 0; i < (int)rootdir.size(); ++i)
        array[i] = rootdir[i];
    int result = CPodeSetRootDirection(updRep().cpode_mem, array);
    delete[] array;
    return result;
}

int CPodes::step(Real tout, Real* tret, 
         Vector& yout, Vector& ypout, StepMode mode)
{
    // The pretty code commented out here heap allocates the N_VectorContent_SimTK
    // objects, which are then deleted in the N_Vector_SimTK destructor:
    //   N_Vector_SimTK nv_yout(yout);
    //   N_Vector_SimTK nv_ypout(ypout);
    // Instead, we'll use stack allocated Content objects, and
    // construct the N_Vectors by hand here avoiding the per-step
    // heap allocation/deallocation, which gets tiresome after a while.

    N_VectorContent_SimTK content_yout(yout);
    N_VectorContent_SimTK content_ypout(ypout);
    _generic_N_Vector nv_yout, nv_ypout;
    nv_yout.content = (void*)&content_yout;
    nv_ypout.content = (void*)&content_ypout;
    nv_yout.ops = nv_ypout.ops = N_Vector_Ops_SimTK::getSundialsOpsPtr();

    return CPode(updRep().cpode_mem,tout,tret,
                 &nv_yout, &nv_ypout,
                 mapStepMode(mode));
}

int CPodes::getDky(Real t, int k, Vector& dky) {
    N_Vector_SimTK nv_dky(dky);
    return CPodeGetDky(updRep().cpode_mem,t,k,&nv_dky);
}

int CPodes::getQuad(Real t, Vector& yQout) {
    N_Vector_SimTK nv_yQout(yQout);
    return CPodeGetQuad(updRep().cpode_mem,t,&nv_yQout);
}
int CPodes::getQuadDky(Real t, int k, Vector& dky) {
    N_Vector_SimTK nv_dky(dky);
    return CPodeGetQuadDky(updRep().cpode_mem,t,k,&nv_dky);
}

int CPodes::getWorkSpace(int* lenrw, int* leniw) {
    long llenrw, lleniw; 
    int stat = CPodeGetWorkSpace(updRep().cpode_mem,&llenrw,&lleniw);
    *lenrw = (int)llenrw; *leniw = (int)lleniw;
    return stat;
}
int CPodes::getNumSteps(int* nsteps) {
    long lnsteps;
    int stat = CPodeGetNumSteps(updRep().cpode_mem,&lnsteps);
    *nsteps = (int)lnsteps;
    return stat;
}
int CPodes::getNumFctEvals(int* nfevals) {
    long lnfevals;
    int stat = CPodeGetNumFctEvals(updRep().cpode_mem,&lnfevals);
    *nfevals = (int)lnfevals;
    return stat;
}
int CPodes::getNumLinSolvSetups(int* nlinsetups) {
    long lnlinsetups;
    int stat = CPodeGetNumLinSolvSetups(updRep().cpode_mem,&lnlinsetups);
    *nlinsetups = (int)lnlinsetups;
    return stat;
}
int CPodes::getNumErrTestFails(int* netfails) {
    long lnetfails;
    int stat = CPodeGetNumErrTestFails(updRep().cpode_mem,&lnetfails);
    *netfails = (int)lnetfails;
    return stat;
}
int CPodes::getLastOrder(int* qlast) {
    return CPodeGetLastOrder(updRep().cpode_mem,qlast);
}
int CPodes::getCurrentOrder(int* qcur) {
    return CPodeGetCurrentOrder(updRep().cpode_mem,qcur);
}
int CPodes::getNumStabLimOrderReds(int* nslred) {
    long lnslred;
    int stat = CPodeGetNumStabLimOrderReds(updRep().cpode_mem,&lnslred);
    *nslred = (int)lnslred;
    return stat;
}
int CPodes::getActualInitStep(Real* hinused) {
    return  CPodeGetActualInitStep(updRep().cpode_mem,hinused);
}
int CPodes::getLastStep(Real* hlast) {
    return CPodeGetLastStep(updRep().cpode_mem,hlast);
}
int CPodes::getCurrentStep(Real* hcur) {
    return CPodeGetCurrentStep(updRep().cpode_mem,hcur);
}
int CPodes::getCurrentTime(Real* tcur) {
    return CPodeGetCurrentTime(updRep().cpode_mem,tcur);
}
int CPodes::getTolScaleFactor(Real* tolsfac) {
    return CPodeGetTolScaleFactor(updRep().cpode_mem,tolsfac);
}
int CPodes::getErrWeights(Vector& eweight) {
    N_Vector_SimTK nv_eweight(eweight);
    return CPodeGetErrWeights(updRep().cpode_mem, &nv_eweight);
}
int CPodes::getEstLocalErrors(Vector& ele) {
    N_Vector_SimTK nv_ele(ele);
    return CPodeGetEstLocalErrors(updRep().cpode_mem, &nv_ele);
}
int CPodes::getNumGEvals(int* ngevals) {
    long lngevals;
    int stat = CPodeGetNumGEvals(updRep().cpode_mem,&lngevals);
    *ngevals = (int)lngevals;
    return stat;
}
int CPodes::getRootInfo(int* rootsfound) {
    return CPodeGetRootInfo(updRep().cpode_mem,rootsfound);
}
int CPodes::getRootWindow(Real* tLo, Real* tHi) {
    return CPodeGetRootWindow(updRep().cpode_mem,tLo,tHi);
}
int CPodes::getIntegratorStats(int* nsteps,
                          int* nfevals, int* nlinsetups,
                          int* netfails, int* qlast,
                          int* qcur, Real* hinused, Real* hlast,
                          Real* hcur, Real* tcur) 
{
    long lnsteps, lnfevals, lnlinsetups, lnetfails;
    int stat = CPodeGetIntegratorStats(updRep().cpode_mem,&lnsteps,
                          &lnfevals,&lnlinsetups,
                          &lnetfails,qlast,
                          qcur,hinused,hlast,
                          hcur,tcur);
    *nsteps = (int)lnsteps; *nfevals = (int)lnfevals;
    *nlinsetups = (int)lnlinsetups; *netfails = (int)lnetfails;
    return stat;
}

int CPodes::getNumNonlinSolvIters(int* nniters) {
    long lnniters;
    int stat = CPodeGetNumNonlinSolvIters(updRep().cpode_mem,&lnniters);
    *nniters = (int)lnniters;
    return stat;
}
int CPodes::getNumNonlinSolvConvFails(int* nncfails) {
    long lnncfails;
    int stat = CPodeGetNumNonlinSolvConvFails(updRep().cpode_mem,&lnncfails);
    *nncfails = (int)lnncfails;
    return stat;
}
int CPodes::getNonlinSolvStats(int* nniters, int* nncfails) 
{
    long lnniters, lnncfails;
    int stat = CPodeGetNonlinSolvStats(updRep().cpode_mem,&lnniters,&lnncfails);
    *nniters = (int)lnniters; *nncfails = (int)lnncfails;
    return stat;
}
int CPodes::getProjNumProj(int* nproj) {
    long lnproj;
    int stat = CPodeGetProjNumProj(updRep().cpode_mem,&lnproj);
    *nproj = (int)lnproj;
    return stat;
}
int CPodes::getProjNumCnstrEvals(int* nce) {
    long lnce;
    int stat = CPodeGetProjNumCnstrEvals(updRep().cpode_mem,&lnce);
    *nce = (int)lnce;
    return stat;
}
int CPodes::getProjNumLinSolvSetups(int* nsetupsP) {
    long lnsetupsP;
    int stat = CPodeGetProjNumLinSolvSetups(updRep().cpode_mem,&lnsetupsP);
    *nsetupsP = (int)lnsetupsP;
    return stat;
}
int CPodes::getProjNumFailures(int* nprf) {
    long lnprf;
    int stat = CPodeGetProjNumFailures(updRep().cpode_mem,&lnprf);
    *nprf = (int)lnprf;
    return stat;
}
int CPodes::getProjStats(int* nproj,
                 int* nce, int* nsetupsP,
                 int* nprf) 
{
    long lnproj,lnce,lnsetupsP,lnprf;
    int stat = CPodeGetProjStats(updRep().cpode_mem,&lnproj,&lnce,&lnsetupsP,&lnprf);
    *nproj=(int)lnproj; *nce=(int)lnce;
    *nsetupsP=(int)lnsetupsP; *nprf=(int)lnprf;
    return stat;
}
int CPodes::getQuadNumFunEvals(int* nqevals) {
    long lnqevals;
    int stat = CPodeGetQuadNumFunEvals(updRep().cpode_mem,&lnqevals);
    *nqevals = (int)lnqevals;
    return stat;
}
int CPodes::getQuadErrWeights(Vector& eQweight) {
    N_Vector_SimTK nv_eQweight(eQweight);
    return CPodeGetQuadErrWeights(updRep().cpode_mem, &nv_eQweight);
}
char* CPodes::getReturnFlagName(int flag) {
    return CPodeGetReturnFlagName(flag);
}

int CPodes::dlsSetJacFn(void* jac, void* jac_data) {
    return CPDlsSetJacFn(updRep().cpode_mem,jac,jac_data);
}
int CPodes::dlsGetWorkSpace(int* lenrwLS, int* leniwLS) {
    long llenrwLS, lleniwLS;
    int stat = CPDlsGetWorkSpace(updRep().cpode_mem,&llenrwLS,&lleniwLS);
    *lenrwLS = (int)llenrwLS; *leniwLS = (int)lleniwLS;
    return stat;
}
int CPodes::dlsGetNumJacEvals(int* njevals) {
    long lnjevals;
    int stat = CPDlsGetNumJacEvals(updRep().cpode_mem,&lnjevals);
    *njevals = (int)lnjevals;
    return stat;
}
int CPodes::dlsGetNumFctEvals(int* nfevalsLS) {
    long lnfevalsLS;
    int stat = CPDlsGetNumFctEvals(updRep().cpode_mem,&lnfevalsLS);
    *nfevalsLS = (int)lnfevalsLS;
    return stat;
}
int CPodes::dlsGetLastFlag(int* flag) {
    return CPDlsGetLastFlag(updRep().cpode_mem,flag);
}
char* CPodes::dlsGetReturnFlagName(int flag) {
    return CPDlsGetReturnFlagName(flag);
}
int CPodes::dlsProjSetJacFn(void* jacP, void* jacP_data) {
    return CPDlsProjSetJacFn(updRep().cpode_mem,jacP,jacP_data);
}
int CPodes::dlsProjGetNumJacEvals(int* njPevals) {
    long lnjPevals;
    int stat = CPDlsProjGetNumJacEvals(updRep().cpode_mem,&lnjPevals);
    *njPevals = (int)lnjPevals;
    return stat;
}
int CPodes::dlsProjGetNumFctEvals(int* ncevalsLS) {
    long lncevalsLS;
    int stat = CPDlsProjGetNumFctEvals(updRep().cpode_mem,&lncevalsLS);
    *ncevalsLS = (int)lncevalsLS;
    return stat;
}

int CPodes::lapackDense(int N) {
    return CPLapackDense(updRep().cpode_mem,N);
}
int CPodes::lapackBand(int N, int mupper, int mlower) {
    return CPLapackBand(updRep().cpode_mem,N,mupper,mlower);
}
int CPodes::lapackDenseProj(int Nc, int Ny, ProjectionFactorizationType fact_type) {
    return CPLapackDenseProj(updRep().cpode_mem,Nc,Ny,
        mapProjectionFactorizationType(fact_type));
}



// Client-side function registration
void CPodes::registerExplicitODEFunc(CPodes::ExplicitODEFunc f) {
    updRep().explicitODEFunc = f;
}
void CPodes::registerImplicitODEFunc(CPodes::ImplicitODEFunc f) {
    updRep().implicitODEFunc = f;
}
void CPodes::registerConstraintFunc(CPodes::ConstraintFunc f) {
    updRep().constraintFunc = f;
}
void CPodes::registerProjectFunc(CPodes::ProjectFunc f) {
    updRep().projectFunc = f;
}
void CPodes::registerQuadratureFunc(CPodes::QuadratureFunc f) {
    updRep().quadratureFunc = f;
}
void CPodes::registerRootFunc(CPodes::RootFunc f) {
    updRep().rootFunc = f;
}
void CPodes::registerWeightFunc(CPodes::WeightFunc f) {
    updRep().weightFunc = f;
}
void CPodes::registerErrorHandlerFunc(CPodes::ErrorHandlerFunc f) {
    updRep().errorHandlerFunc = f;
}

/////////////////////////////////
// CPodesSystem IMPLEMENTATION //
/////////////////////////////////

// These are the default implementations for the CPodesSystem virtual functions.
// Note that this DOES NOT cause binary compatibility problems. The addresses of
// these functions will be supplied from the library side, but these addresses will
// get filled in to the default virtual function table on the *client* side which
// knows where to put each function by name.

int CPodesSystem::explicitODE(Real, const Vector&, Vector&) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "CPodesSystem", "explicitODE"); 
    return std::numeric_limits<int>::min();
}
int CPodesSystem::implicitODE(Real, const Vector&, const Vector&, Vector&) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "CPodesSystem", "implicitODE"); 
    return std::numeric_limits<int>::min();
}

int CPodesSystem::constraint(Real, const Vector&, Vector&) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "CPodesSystem", "constraint"); 
    return std::numeric_limits<int>::min();
}

int CPodesSystem::project(Real, const Vector&, Vector&, Real, Vector&) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "CPodesSystem", "project"); 
    return std::numeric_limits<int>::min();
}

int CPodesSystem::quadrature(Real, const Vector&, Vector&) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "CPodesSystem", "quadrature"); 
    return std::numeric_limits<int>::min();
}

int CPodesSystem::root(Real, const Vector&, const Vector&, Vector&) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "CPodesSystem", "root"); 
    return std::numeric_limits<int>::min();
}

int CPodesSystem::weight(const Vector&, Vector&) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "CPodesSystem", "weight"); 
    return std::numeric_limits<int>::min();
}

void CPodesSystem::errorHandler(int, const char*, const char*, char*) const {
    SimTK_THROW2(Exception::UnimplementedVirtualMethod, "CPodesSystem", "errorHandler"); 
}

} // namespace SimTK


