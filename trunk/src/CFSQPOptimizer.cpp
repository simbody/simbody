/* Portions copyright (c) 2006 Stanford University and Jack Middleton.
 * Contributors: Eran Guendelman, Frank C. Anderson
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
#include "CFSQPOptimizer.h"
#include <string>

// ----------------------------------------------------------------------
// Code to manage the library loading
// ----------------------------------------------------------------------

#ifndef _WIN32
#define SIMTK_PORTABLE_HMODULE void *
#define SIMTK_PORTABLE_HINSTANCE void *
#define WINAPI
#include <dlfcn.h>
#define GetProcAddress(handle, proc) dlsym(handle, proc)
#define LoadLibrary(name) dlopen(name, RTLD_LAZY | RTLD_GLOBAL)
const char *LoadLibraryErrorMessage() { return dlerror(); }
#else
#include <windows.h>
#define SIMTK_PORTABLE_HMODULE HMODULE
#define SIMTK_PORTABLE_HINSTANCE HINSTANCE
const char *LoadLibraryErrorMessage() { return 0; }
#endif

#ifndef _WIN32
static const char *CFSQP_LIBRARY_NAME = "libosimCFSQP.so";
#else
static const char *CFSQP_LIBRARY_NAME = "osimCFSQP.dll";
#endif

typedef void* (*CFSQP_FUNCTION)(int,int,int,int,int,int,int,int,int,int *,int,int,
                                int,int *,double,double,double,double,double *,
                                double *,double *,double *,double *,double *,
                                void (*)(int,int,double *,double *,void *),
                                void (*)(int,int,double *,double *,void *),
                                void (*)(int,int,double *,double *,
                                void (*)(int,int,double *,double *,void *),void *),
                                void (*)(int,int,double *,double *,
                                void (*)(int,int,double *,double *,void *),void *),
                                void *);
CFSQP_FUNCTION cfsqp_fptr = 0;

// ----------------------------------------------------------------------

#define USE_CONSTRAINT_CACHE

using namespace SimTK;

//=============================================================================
// CONSTRUCTION
//=============================================================================
//_____________________________________________________________________________
/**
 * Destructor.
 */
CFSQPOptimizer::
~CFSQPOptimizer()
{
    delete[] _mesh;
    delete[] _p;
    delete[] _c;
    delete[] _lambda;
}

OptimizerRep* CFSQPOptimizer::clone() const {
    return( new CFSQPOptimizer(*this) );
}

//_____________________________________________________________________________
CFSQPOptimizer::CFSQPOptimizer(const OptimizerSystem& sys)
    : OptimizerRep(sys)
{
    bindToCFSQPLibrary();

	_infinity = 1.0e10;
	_mode = 100;
	_epseqn = 1.0e-8;
	_inform = 0;
	_udelta = 0.0;

    int nx = sys.getNumParameters();

	// NEW MAX ITERATIONS
	maxIterations = 4 * nx;

	// SQP STUFF- NOT CURRENTLY SUPPORTED
	int i;
	_ncsrl = 0;
	_ncsrn = 0;
	_nfsr = 0;
	int lenmesh = _nfsr+_ncsrn+_ncsrl;  if(lenmesh<1) lenmesh=1;
	_mesh = new int[lenmesh];
	for(i=0;i<lenmesh;i++)  _mesh[i] = 0;

	// PERFORMANCE
	int lenp = 1 - _nfsr; // _target->getNumContacts() - _nfsr;
	for(int i=0;i<_nfsr;i++)  lenp += _mesh[i];
	if(lenp<1) lenp = 1;
    _p = new double[lenp];

	// CONSTRAINTS
	int lenc = sys.getNumConstraints() - _ncsrl - _ncsrn;
	for(int i=0;i<_ncsrn;i++)  lenc += _mesh[i+_nfsr];
	for(int i=0;i<_ncsrl;i++)  lenc += _mesh[i+_nfsr+_ncsrn];
	if(lenc<1) lenc = 1;
	_c = new double[lenc];

	// LAGRANGE MULTIPLIERS
	_lambda = new double[nx+lenp+lenc];
}

bool CFSQPOptimizer::isAvailable()
{
    try {
        bindToCFSQPLibrary();
    } catch (const SimTK::Exception::Base &ex) {
       std::cout << ex.getMessageText() << std::endl;
       return false;
    }
    return true;
}

void CFSQPOptimizer::bindToCFSQPLibrary()
{
    if(cfsqp_fptr) return;

    SIMTK_PORTABLE_HINSTANCE libraryHandle = LoadLibrary(CFSQP_LIBRARY_NAME);
    if(!libraryHandle) {
        const char *msg=LoadLibraryErrorMessage();
        if(msg)
            SimTK_THROW1(SimTK::Exception::Cant, std::string("Could not load CFSQP library '") + CFSQP_LIBRARY_NAME + "': " + msg); 
        else
            SimTK_THROW1(SimTK::Exception::Cant, std::string("Could not load CFSQP library '") + CFSQP_LIBRARY_NAME + "'"); 
    }

    cfsqp_fptr = (CFSQP_FUNCTION)GetProcAddress(libraryHandle, "cfsqp");
    if(!cfsqp_fptr) {
        SimTK_THROW1(SimTK::Exception::Cant, std::string("CFSQP library '") + CFSQP_LIBRARY_NAME + "' not valid (could not find 'cfsqp' symbol)"); 
    }

    std::cout << "Successfully linked to CFSQP library '" << CFSQP_LIBRARY_NAME << "'" << std::endl;
}

//=============================================================================
// COMPUTE OPTIMAL CONTROLS
//=============================================================================
//_____________________________________________________________________________
/**
 * Compute a set of optimal controls, given the input controls (xin) and
 * the current state of the optimization target.
 *
 * Whether or not the optimization terminates normally, the latest value of
 * the controls are copied to xout.  It is safe to use the same pointer
 * for xin and xout.  However, in all cases the calling routine must
 * allocate enough space for xin and xout.
 *
 * @param xin Values of the controls.
 * @param xout Optimal values of the controls.
 *
 * @return Parameter inform of cfsqp. -1 means a fatal error.
 */
Real CFSQPOptimizer::
optimize(Vector &results)
{
    const OptimizerSystem &sys = getOptimizerSystem();

    int nx = sys.getNumParameters();

	double *bl, *bu;
	if(sys.getHasLimits()) {
		sys.getParameterLimits(&bl,&bu);
	} else {
		bl = new double[nx];
		bu = new double[nx];
		for(int i=0; i<nx; i++) {
			bl[i] = -_infinity;
			bu[i] =  _infinity;
		}
	}

    double *x = &results[0];
    int numObjectiveFunctions = 1;

	// Clear cache before starting optimization (used to speed up constraint computations)
	clearCache();

	// OPTIMIZE
	cfsqp_fptr(getOptimizerSystem().getNumParameters(),numObjectiveFunctions,_nfsr,
        sys.getNumNonlinearInequalityConstraints(), sys.getNumInequalityConstraints(),
        sys.getNumNonlinearEqualityConstraints(), sys.getNumEqualityConstraints(),
		_ncsrl,_ncsrn,_mesh,
		_mode,diagnosticsLevel,maxIterations,&_inform,_infinity,convergenceTolerance,_epseqn,_udelta,
		bl,bu,x,_p,_c,_lambda,pFunc,cFunc,dpdxFunc,dcdxFunc,(void *)this);

	if(!getOptimizerSystem().getHasLimits()) {
		delete[] bl;
		delete[] bu;
	}

    if(diagnosticsLevel > 0) PrintInform(_inform,std::cout);
    
    if (_inform != 0 && _inform != 3 && _inform != 4 && _inform != 8) {
        char buf[1024];
        sprintf(buf, "CFSQP failed with status = %d",_inform);
        SimTK_THROW1(SimTK::Exception::OptimizerFailed, SimTK::String(buf));
    }

    return *_p;
}

//=============================================================================
// STATIC PERFORMANCE AND CONSTRAINT EVALUATIONS
//=============================================================================
//______________________________________________________________________________
/**
 * Compute the performance criterion.
 */
void CFSQPOptimizer::
pFunc(int nparam,int j,double *x,double *p,void *cd)
{
    CFSQPOptimizer *cfsqp = (CFSQPOptimizer *)cd;
    int nx=cfsqp->getOptimizerSystem().getNumParameters();
    cfsqp->objectiveFunc(cfsqp->getOptimizerSystem(),Vector(nx,x,true),true,*p);
}
//______________________________________________________________________________
/**
 * Compute the derivatives of the performance criterion.
 */
void CFSQPOptimizer::
dpdxFunc(int nparam,int j,double *x,double *dpdx,
	void (*dummy)(int,int,double *,double *,void *),void *cd)
{
    // TODO; support numerical gradients
    CFSQPOptimizer *cfsqp = (CFSQPOptimizer *)cd;
    int nx=cfsqp->getOptimizerSystem().getNumParameters();
    cfsqp->gradientFunc(cfsqp->getOptimizerSystem(),Vector(nx,x,true),true,Vector(nx,dpdx,true));
}

//______________________________________________________________________________
/**
 * Compute the constraints.
 */
void CFSQPOptimizer::
cFunc(int nparam,int j,double *x,double *c,void *cd)
{
    CFSQPOptimizer *cfsqp = (CFSQPOptimizer *)cd;
	int nx=cfsqp->getOptimizerSystem().getNumParameters();
    // special wrapper to deal with caching
	cfsqp->computeConstraint(Vector(nx,x,true),true,*c,j-1);
}
//______________________________________________________________________________
/**
 * Compute the gradient of the constraints.
 */
void CFSQPOptimizer::
dcdxFunc(int nparam,int j,double *x,double *dcdx,
	void (*dummy)(int,int,double *,double *,void *),
	void *cd)
{
    // TODO: support numerical gradients
    CFSQPOptimizer *cfsqp = (CFSQPOptimizer *)cd;
	int nx=cfsqp->getOptimizerSystem().getNumParameters();
	int nc=cfsqp->getOptimizerSystem().getNumConstraints();
    // special wrapper to deal with caching
	cfsqp->computeConstraintGradient(Vector(nx,x,true),true,Vector(nx,dcdx,true),j-1);
}

//=============================================================================
// PRINT
//=============================================================================
//______________________________________________________________________________
/**
 * Print the meaning of the value returned by computeOptimalControls().
 */
void CFSQPOptimizer::
PrintInform(int aInform,std::ostream &aOStream)
{
	switch(aInform) {
		case(0):
			aOStream<<"CFSQP(0): Normal termination.\n";
			break;
		case(1):
			aOStream<<"CFSQP(1): User-provided initial guess is infeasible ";
			aOStream<<"for linear constraints\n";
			aOStream<<"and CFSQP is unable to generate a point satisfying these ";
			aOStream<<"conditions.\n";
			break;
		case(2):
			aOStream<<"CFSQP(2): The user-provided initial guess is infeasible ";
			aOStream<<"for non-linear inequality constraints\n";
			aOStream<<"and linear constraints, and CFSQP is unable to generate ";
			aOStream<<"a point satisfying these constraints.\n";
			aOStream<<"This may be due to insucient accuracy of the QP solver.\n";
			break;
		case(3):
			aOStream<<"CFSQP(3): The maximum number of iterations has been ";
			aOStream<<"reached before a solution was obtained.\n";
			break;
		case(4):
			aOStream<<"CFSQP(4): The line search failed to find a new iterate.";
			aOStream<<" The step size was smaller than\n";
			aOStream<<"the machine precision.\n";
			break;
		case(5):
			aOStream<<"CFSQP(5): Failure of the QP solver in attempting to construct d0.\n";
			break;
		case(6):
			aOStream<<"CFSQP(6): Failure of the QP solver in attempting to construct d1.\n";
			break;
		case(7):
			aOStream<<"CFSQP(7): Input data are not consistent.  Set the print level";
			aOStream<<" greater than 0 for more information.\n";
			break;
		case(8):
			aOStream<<"CFSQP(8): The new iterate is numerically equivalent to ";
			aOStream<<"the previous iterate,\n";
			aOStream<<"though the stopping criterion is not yet satisfied. ";
			aOStream<<"Relaxing the stopping criterion\n";
			aOStream<<"shouldsolve this problem.\n";
			break;
		case(9):
			aOStream<<"CFSQP(9): One of the penalty parameters exceeded ";
			aOStream<<"the largest allowed bound.\n";
			aOStream<<"The algorithm is having trouble satisfying a non-linear ";
			aOStream<<"equality constraint.\n";
			break;
		default:
			aOStream<<"CFSQP("<<aInform<<"): Unrecognized inform value.\n";
	}
}
//=============================================================================
// CACHING
//=============================================================================
void CFSQPOptimizer::
clearCache()
{
#ifdef USE_CONSTRAINT_CACHE
	int nx=getOptimizerSystem().getNumParameters();
	int nc=getOptimizerSystem().getNumConstraints();
	_cachedConstraintJacobian.resize(nc,nx);
	_cachedConstraint.resize(nc);
	_cachedConstraintJacobianParameters.resize(0);
	_cachedConstraintParameters.resize(0);
#endif
}
int CFSQPOptimizer::
computeConstraint(const SimTK::Vector &x, const bool new_coefficients, double &c, int ic) const
{
	int nx=getOptimizerSystem().getNumParameters();
	int nc=getOptimizerSystem().getNumConstraints();
	int status = 0;
#ifdef USE_CONSTRAINT_CACHE
	bool cached_value_available = false;
	if(_cachedConstraintParameters.size()) {
		cached_value_available = true;
		for(int i=0; i<nx; i++) 
			if(x[i] != _cachedConstraintParameters[i]) {
				cached_value_available = false;
				break;
			}
	}

	if(!cached_value_available) {
		status = constraintFunc(getOptimizerSystem(),x,new_coefficients,_cachedConstraint);
		_cachedConstraintParameters.resize(nx);
		_cachedConstraintParameters = x;
	} 
	c = _cachedConstraint[ic];
#else
	SimTK::Vector allc(nc);
	status = constraintFunc(getOptimizerSystem(),x,new_coefficients,allc);
	c=allc[ic];
#endif
	return status;
}
int CFSQPOptimizer::
computeConstraintGradient(const SimTK::Vector &x, const bool new_coefficients, SimTK::Vector &dcdx, int ic) const
{
	int nx=getOptimizerSystem().getNumParameters();
	int nc=getOptimizerSystem().getNumConstraints();
	int status = 0;
#ifdef USE_CONSTRAINT_CACHE
	bool cached_value_available = false;
	if(_cachedConstraintJacobianParameters.size()) {
		cached_value_available = true;
		for(int i=0; i<nx; i++) 
			if(x[i] != _cachedConstraintJacobianParameters[i]) {
				cached_value_available = false;
				break;
			}
	}

	if(!cached_value_available) {
		status = constraintJacobian(getOptimizerSystem(),x,new_coefficients,_cachedConstraintJacobian);
		_cachedConstraintJacobianParameters.resize(nx);
		_cachedConstraintJacobianParameters = x;
	} 
	for(int col=0;col<nx;col++) dcdx[col]=_cachedConstraintJacobian(ic,col);
#else
	SimTK::Matrix jacobian(nc,nx);
	status = constraintJacobian(getOptimizerSystem(),x,new_coefficients,jacobian);
	for(int col=0;col<nx;col++) dcdx[col]=jacobian(ic,col);
#endif
	return status;
}
