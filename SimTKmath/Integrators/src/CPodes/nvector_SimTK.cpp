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
 * This is the implementation of N_Vector_SimTK anad related classes. It
 * is a concrete implementation of the generic interface defined by
 * Sundials for its abstract data structure N_Vector. The idea here is to
 * provide the needed operations using the SimTK Vector class.
 */

#include "simmath/internal/SimTKcpodes.h"
#include "nvector_SimTK.h"

using SimTK::Real;
using SimTK::Vector;


// These are the "virtual functions" required by the N_Vector interface.
// Each is an implementation of one of the function types defined in
// the Sundials C struct _generic_N_Vector_Ops. At startup we'll initialize
// a single static "ops" object and use it for all the N_Vector_SimTK's 
// that ever get created.

static N_Vector    nvclone_SimTK(N_Vector);
static N_Vector    nvcloneempty_SimTK(N_Vector);
static void        nvdestroy_SimTK(N_Vector);
static void        nvspace_SimTK(N_Vector, long int*, long int*);
static realtype*   nvgetarraypointer_SimTK(N_Vector);
static void        nvsetarraypointer_SimTK(realtype*, N_Vector);
static void        nvlinearsum_SimTK(realtype, N_Vector, realtype, N_Vector, N_Vector); 
static void        nvconst_SimTK(realtype, N_Vector);
static void        nvprod_SimTK(N_Vector, N_Vector, N_Vector);
static void        nvdiv_SimTK(N_Vector, N_Vector, N_Vector);
static void        nvscale_SimTK(realtype, N_Vector, N_Vector);
static void        nvabs_SimTK(N_Vector, N_Vector);
static void        nvinv_SimTK(N_Vector, N_Vector);
static void        nvaddconst_SimTK(N_Vector, realtype, N_Vector);
static realtype    nvdotprod_SimTK(N_Vector, N_Vector);
static realtype    nvmaxnorm_SimTK(N_Vector);
static realtype    nvwrmsnorm_SimTK(N_Vector, N_Vector);
static realtype    nvwrmsnormmask_SimTK(N_Vector, N_Vector, N_Vector);
static realtype    nvmin_SimTK(N_Vector);
static realtype    nvwl2norm_SimTK(N_Vector, N_Vector);
static realtype    nvl1norm_SimTK(N_Vector);
static void        nvcompare_SimTK(realtype, N_Vector, N_Vector);
static booleantype nvinvtest_SimTK(N_Vector, N_Vector);
static booleantype nvconstrmask_SimTK(N_Vector, N_Vector, N_Vector);
static realtype    nvminquotient_SimTK(N_Vector, N_Vector);

// The default constructor for this static const member takes care
// of initializing all the function pointers.
const N_Vector_Ops_SimTK N_Vector_Ops_SimTK::Ops;

// This is the default constructor. It is used exactly once to initialize
// the static const member variable above.
N_Vector_Ops_SimTK::N_Vector_Ops_SimTK() {
  nvclone           = nvclone_SimTK;
  nvcloneempty      = nvcloneempty_SimTK;
  nvdestroy         = nvdestroy_SimTK;
  nvspace           = nvspace_SimTK;
  nvgetarraypointer = nvgetarraypointer_SimTK;
  nvsetarraypointer = nvsetarraypointer_SimTK;
  nvlinearsum       = nvlinearsum_SimTK;
  nvconst           = nvconst_SimTK;
  nvprod            = nvprod_SimTK;
  nvdiv             = nvdiv_SimTK;
  nvscale           = nvscale_SimTK;
  nvabs             = nvabs_SimTK;
  nvinv             = nvinv_SimTK;
  nvaddconst        = nvaddconst_SimTK;
  nvdotprod         = nvdotprod_SimTK;
  nvmaxnorm         = nvmaxnorm_SimTK;
  nvwrmsnorm        = nvwrmsnorm_SimTK;
  nvwrmsnormmask    = nvwrmsnormmask_SimTK;
  nvmin             = nvmin_SimTK;
  nvwl2norm         = nvwl2norm_SimTK;
  nvl1norm          = nvl1norm_SimTK;
  nvcompare         = nvcompare_SimTK;
  nvinvtest         = nvinvtest_SimTK;
  nvconstrmask      = nvconstrmask_SimTK;
  nvminquotient     = nvminquotient_SimTK;
}

//////////////////////////////
// DEFINITIONS OF OPERATORS //
//////////////////////////////

// N_VClone
// Allocate a new N_Vector of the same size & type as
// existing vector w, but do not copy the data.
static N_Vector    
nvclone_SimTK(N_Vector w) {
    return N_Vector_SimTK::nvclone(w);
}

// N_VCloneEmpty
// Allocate a new N_Vector of the same type as
// existing vector w, but allocate no space.
static N_Vector    
nvcloneempty_SimTK(N_Vector w) {
    assert(N_Vector_SimTK::isA(w));
    return new N_Vector_SimTK(); 
}

// N_VDestroy
// Free all space associated with a given N_Vector. This
// may or may not destroy the underlying data, depending
// on whether the given N_Vector is the owner.
static void        
nvdestroy_SimTK(N_Vector v) {
    delete N_Vector_SimTK::updDowncast(v);
}

// N_VSpace
// Returns the storage requirements for one N_Vector, specified
// as a number of 'realtype' objects and a number of 'int'
// objects.
// sherm 061128: I discussed this with Radu and he confirmed
// that the routine does not make sense and is never used by
// Sundials for anything significant. It may show up in some
// kind of unused storage estimate. So I'm just returning
// the length as the number of reals and "1" as the number
// of ints.
static void        
nvspace_SimTK(N_Vector v, long int* lrw, long int* liw) {
    *lrw = N_Vector_SimTK::getVector(v).size();
    *liw = 1; // doesn't mean anything
}

// N_VGetArrayPointer
// This returns the address of the first entry in the N_Vector's
// data, with the built-in assumption that the data is stored
// contiguously as an array of reals.
static realtype*   
nvgetarraypointer_SimTK(N_Vector nvx) {
    Vector& x = N_Vector_SimTK::updVector(nvx);
    return x.updContiguousScalarData();
}

// N_VSetArrayPointer
// Replace the data portion of an N_Vector with a new one. This
// assumes that the existing N_Vector uses contiguous data
// and is the owner of its data. Note: we are NOT going to
// delete the old data here; we're assuming that the N_Vector
// user is holding onto a pointer to it obtained with
// N_VGetArrayPointer above. If that's not right then there
// is going to be a leak here.
static void        
nvsetarraypointer_SimTK(realtype* vdata, N_Vector nvz) {
    Vector& z = N_Vector_SimTK::updVector(nvz);
    assert(z.hasContiguousData());
    realtype* oldData;
    z.swapOwnedContiguousScalarData(vdata, z.size(), oldData);
    // don't do anything with oldData here
}

// N_VLinearSum
// z = ax + by
// This will blow up if x and y aren't the same size.
// z will get resized if necessary.
static void        
nvlinearsum_SimTK(realtype a, N_Vector nvx, realtype b, N_Vector nvy, N_Vector nvz) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);
    const Vector& y = N_Vector_SimTK::getVector(nvy);
    Vector&       z = N_Vector_SimTK::updVector(nvz);

    assert(y.size() == x.size() && z.size() == x.size());

    z = a*x + b*y;
}

// N_VConst
// Set all components of vector z to value c.
static void        
nvconst_SimTK(realtype c, N_Vector nvz) {
    Vector& z = N_Vector_SimTK::updVector(nvz);
    z = c;
}

// N_VProd
// z = x .* y (componentwise multiply)
static void        
nvprod_SimTK(N_Vector nvx, N_Vector nvy, N_Vector nvz) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);
    const Vector& y = N_Vector_SimTK::getVector(nvy);
    Vector&       z = N_Vector_SimTK::updVector(nvz);

    const int sz = x.size();
    assert(y.size() == sz && z.size() == sz);

    // TODO: should be a built-in Vector operation for this
    const Real* xp = x.getContiguousScalarData();
    const Real* yp = y.getContiguousScalarData();
    Real*       zp = z.updContiguousScalarData();

    for (int i=0; i<sz; ++i)
        zp[i] = xp[i]*yp[i];
}

// N_VDiv
// z = x ./ y (componentwise divide)
static void        
nvdiv_SimTK(N_Vector nvx, N_Vector nvy, N_Vector nvz) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);
    const Vector& y = N_Vector_SimTK::getVector(nvy);
    Vector&       z = N_Vector_SimTK::updVector(nvz);

    const int sz = x.size();
    assert(y.size() == sz && z.size() == sz);

    // TODO: should be a built-in Vector operation for this
    const Real* xp = x.getContiguousScalarData();
    const Real* yp = y.getContiguousScalarData();
    Real*       zp = z.updContiguousScalarData();

    for (int i=0; i<sz; ++i)
        zp[i] = xp[i]/yp[i];
}

// N_VScale
// z = c*x
static void        
nvscale_SimTK(realtype c, N_Vector nvx, N_Vector nvz) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);
    Vector&       z = N_Vector_SimTK::updVector(nvz);

    assert(z.size() == x.size());
    z = c*x;
}

// N_VAbs
// Set zi = |xi|, that is, componentwise absolute value.
static void        
nvabs_SimTK(N_Vector nvx, N_Vector nvz) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);
    Vector&       z = N_Vector_SimTK::updVector(nvz);

    assert(z.size() == x.size());
    z = x.abs();
}

// N_VInv
// zi = 1/xi, that is, componentwise inversion.
static void        
nvinv_SimTK(N_Vector nvx, N_Vector nvz) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);
    Vector&       z = N_Vector_SimTK::updVector(nvz);

    const int sz = x.size();
    assert(z.size() == sz);

    // TODO: should be a built-in Vector operation for this
    const Real* xp = x.getContiguousScalarData();
    Real*       zp = z.updContiguousScalarData();

    for (int i=0; i<sz; ++i)
        zp[i] = 1/xp[i];
}

// N_VAddConst
// zi = xi + b, that is, componentwise scalar addition.
static void        
nvaddconst_SimTK(N_Vector nvx, realtype b, N_Vector nvz) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);
    Vector&       z = N_Vector_SimTK::updVector(nvz);

    const int sz = x.size();
    assert(z.size() == sz);

    // TODO: should be a built-in Vector operation for this
    const Real* xp = x.getContiguousScalarData();
    Real*       zp = z.updContiguousScalarData();
    for (int i=0; i<sz; ++i)
        zp[i] = xp[i] + b;
}

// N_VDotProd
// result = dot(x,y)
static realtype    
nvdotprod_SimTK(N_Vector nvx, N_Vector nvy) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);
    const Vector& y = N_Vector_SimTK::getVector(nvy);

    assert(y.size() == x.size());

    return ~x * y;
}

// N_VMaxNorm
// result = max_i |xi|, that is, return the absolute value
// of the element whose absolute value is largest.
static realtype    
nvmaxnorm_SimTK(N_Vector nvx) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);

    const int sz = x.size();
    assert(sz > 0);

    // TODO: should be a built-in Vector operation for this
    const Real* xp = x.getContiguousScalarData();
    Real maxabs = std::abs(xp[0]);
    for (int i=1; i<sz; ++i) {
        const Real absval = std::abs(xp[i]);
        if (absval > maxabs) maxabs=absval;
    }
    return maxabs;
}

// N_VWrmsNorm
// result = sqrt( sum_i((xi*wi)^2)/n )
// that is, weighted RMS norm of x.
static realtype    
nvwrmsnorm_SimTK(N_Vector nvx, N_Vector nvw) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);
    const Vector& w = N_Vector_SimTK::getVector(nvw);

    const int sz = x.size();
    assert(sz > 0 && w.size() == sz);

    // TODO: should be a built-in Vector operation for this
    const Real* xp = x.getContiguousScalarData();
    const Real* wp = w.getContiguousScalarData();

    Real sumsq = 0;
    for (int i=0; i<sz; ++i) {
        const Real xw = xp[i]*wp[i];
        sumsq += xw*xw;
    }
    return std::sqrt( sumsq / sz );
}

// N_VWrmsNormMask
// result = sqrt( sum_id[i]!=0((xi*wi)^2)/n );
// that is, weighted RMS norm of x including only those terms
// where id[i] != 0.
static realtype    
nvwrmsnormmask_SimTK(N_Vector nvx, N_Vector nvw, N_Vector nvid) {
    const Vector& x  = N_Vector_SimTK::getVector(nvx);
    const Vector& w  = N_Vector_SimTK::getVector(nvw);
    const Vector& id = N_Vector_SimTK::getVector(nvid);

    const int sz = x.size();
    assert(sz > 0 && w.size()==sz && id.size()==sz);

    const Real* xp  = x.getContiguousScalarData();
    const Real* wp  = w.getContiguousScalarData();
    const Real* idp = id.getContiguousScalarData();

    Real sumsq = 0;
    for (int i=0; i<sz; ++i)
        if (idp[i] != 0) {
            const Real xw = xp[i]*wp[i];
            sumsq += xw*xw;
        }
    return std::sqrt( sumsq / sz );
}

// N_VMin
// Return the smallest element of x.
static realtype    
nvmin_SimTK(N_Vector nvx) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);

    const int sz = x.size();
    assert(sz > 0);

    // TODO: should be a built-in Vector operation for this
    const Real* xp  = x.getContiguousScalarData();

    Real minelt = xp[0];
    for (int i=1; i<sz; ++i)
        if (xp[i] < minelt) minelt=xp[i];
    return minelt;
}

// N_VWL2Norm
//   result = sqrt( sum_i( (xi*wi)^2 ) )
// Return the weighted Euclidean L2 norm of x.
static realtype    
nvwl2norm_SimTK(N_Vector nvx, N_Vector nvw) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);
    const Vector& w = N_Vector_SimTK::getVector(nvw);

    const int sz = x.size();
    assert(w.size() == sz);

    // TODO: should be a built-in Vector operation for this
    const Real* xp  = x.getContiguousScalarData();
    const Real* wp  = w.getContiguousScalarData();

    Real sumsq = 0;
    for (int i=0; i<sz; ++i) {
        const Real xw = xp[i]*wp[i];
        sumsq += xw*xw;
    }
    return std::sqrt(sumsq);
}

// N_VL1Norm
//   result = sum_i |xi|
// Return the L1 norm of x (sum of absolute values).
static realtype    
nvl1norm_SimTK(N_Vector nvx) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);

    const int sz = x.size();

    // TODO: should be a built-in Vector operation for this
    const Real* xp = x.getContiguousScalarData();

    Real sumabs = 0;
    for (int i=0; i<sz; ++i) {
        sumabs += std::abs(xp[i]);
    }
    return sumabs;
}

// N_VCompare
// Compare components of x to scalar c and return
// z such that zi=1 if |xi|>=c, else 0.
static void        
nvcompare_SimTK(realtype c, N_Vector nvx, N_Vector nvz) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);
    Vector&       z = N_Vector_SimTK::updVector(nvz);

    const int sz = x.size();
    assert(z.size() == sz);

    const Real* xp = x.getContiguousScalarData();
    Real*       zp = z.updContiguousScalarData();

    for (int i=0; i<sz; ++i)
        zp[i] = (std::abs(xp[i]) >= c ? 1 : 0);
}

// N_VInvTest
// Set z[i] = 1/x[i] if ALL x[i] != 0 and return TRUE.
// If ANY x[i]==0 return FALSE with z undefined.
static booleantype 
nvinvtest_SimTK(N_Vector nvx, N_Vector nvz) {
    const Vector& x = N_Vector_SimTK::getVector(nvx);
    Vector&       z = N_Vector_SimTK::updVector(nvz);

    const int sz = x.size();
    assert(z.size() == sz);

    const Real* xp = x.getContiguousScalarData();
    Real*       zp = z.updContiguousScalarData();

    for (int i=0; i<sz; ++i) {
        if (xp[i] == 0) return FALSE;
        zp[i] = 1/xp[i];
    }
    return TRUE;
}

// N_VConstrMask
// Perform constraint tests on xi based on ci:
//    ci ==  2  =>  xi >  0
//    ci ==  1  =>  xi >= 0
//    ci ==  0  =>  xi is unconstrained
//    ci == -1  =>  xi <= 0
//    ci == -2  =>  xi <  0
// We set mask entry mi to 1 if xi failed its constraint, 0 otherwise.
// Return TRUE if all passed (=> m==0), else FALSE.
static booleantype 
nvconstrmask_SimTK(N_Vector nvc, N_Vector nvx, N_Vector nvm) {
    const Vector& c = N_Vector_SimTK::getVector(nvc);
    const Vector& x = N_Vector_SimTK::getVector(nvx);
    Vector&       m = N_Vector_SimTK::updVector(nvm);

    const int sz = x.size();
    assert(c.size()==sz && m.size()==sz);

    m = 0; // assume success
    booleantype allGood = TRUE;

    const Real* cp = c.getContiguousScalarData();
    const Real* xp = x.getContiguousScalarData();
    Real*       mp = m.updContiguousScalarData();

    for (int i=0; i<sz; ++i) {
        if (cp[i]==0) continue;
        if (cp[i]== 2 && xp[i] >  0) continue;
        if (cp[i]== 1 && xp[i] >= 0) continue;
        if (cp[i]==-1 && xp[i] <= 0) continue;
        if (cp[i]==-2 && xp[i] <  0) continue;
        mp[i] = 1;
        allGood = FALSE;
    }

    return allGood;
}

// N_VMinQuotient
// Return the minimum of the quotients num_i/denom_i skipping
// any elements where denom_i==0. If all the denom_i are zero,
// return BIG_REAL.
static realtype    
nvminquotient_SimTK(N_Vector nvnum, N_Vector nvdenom) {
    const Vector& num   = N_Vector_SimTK::getVector(nvnum);
    const Vector& denom = N_Vector_SimTK::getVector(nvdenom);

    const int sz = num.size();
    assert(denom.size() == sz);

    const Real* nump   = num.getContiguousScalarData();
    const Real* denomp = denom.getContiguousScalarData();

    Real result = BIG_REAL;
    for (int i=0; i<sz; ++i) {
        if (denomp[i] == 0) continue;
        const Real quot = nump[i]/denomp[i];
        if (quot < result) result=quot;
    }

    return result;
}
