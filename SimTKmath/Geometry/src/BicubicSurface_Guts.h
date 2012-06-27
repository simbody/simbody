#ifndef SimTK_BICUBIC_SURFACE_GUTS_H_
#define SimTK_BICUBIC_SURFACE_GUTS_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-12 Stanford University and the Authors.        *
 * Authors: Matthew Millard, Michael Sherman                                  *
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
This file defines the internal BicubicSurface::Guts class, which is the private
implementation class of BicubicSurface, and BicubicSurface::PatchHint::Guts
which is used to improve performance for repeated access to the same patch or
nearby patches. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_BicubicHermitePatch.h"
#include "simmath/internal/Geo_BicubicBezierPatch.h"
#include "simmath/internal/BicubicSurface.h"

#include <cassert>

namespace SimTK { 
    
//==============================================================================
//                CLASS BICUBIC SURFACE :: PATCH HINT :: GUTS
//==============================================================================
/** This object is used to hold precalculated data about the most recently
accessed patch to accelerate the common case of repeated access to the same
patch or to nearby patches. **/
class BicubicSurface::PatchHint::Guts {
public:
    Guts() {clear();}
    // Default copy constructor, copy assignment, destructor

    bool isEmpty() const {return x0 < 0;}
    void clear() {
        x0=y0=-1;
        level = -2; // means "everything is invalid"
        #ifndef NDEBUG
            xS=ooxS=ooxS2=ooxS3=yS=ooyS=ooyS2=ooyS3 = NaN;
            f=fx=fy=fxy=fxx=fyy=fxxx=fxxy=fyyy=fxyy = NaN; 
        #endif
    }


    // Patch information. 
    // This is valid whenever level >= -1 and does not change
    // for repeated access anywhere within the same patch.
    int x0, y0; // Indices of the lower-left corner of the patch.
    // These are the precalculated patch dimensions and their reciprocals. 
    // xScale=x[x0+1]-x[x0], yScale=y[y0+1]-y[y0].
    Real xS, ooxS, ooxS2, ooxS3;
    Real yS, ooyS, ooyS2, ooyS3;
    // These are the bicubic patch coefficients a00-a33 for this patch, in
    // column order (a00 a10 a20 a30 a01 a11 ...).
    Vec<16> a;
    // These are the scaled function values at the corners of this patch, in
    // the order f00,f10,f01,f11,
    //           fx00,fx10,fx01,fx11,
    //           fy00,fy10,fy01,fy11,
    //           fxy00,fxy10,fxy01,fxy11
    Vec<16> fV;

    // Point information. 
    // This information records the results of the last 
    // access to the above patch, which was point (x,y,f(x,y)). This will 
    // change whenever a new point is accessed, even if it is on the same patch.

    // -2: nothing, 
    // -1: patch info only, 
    //  0: patch info and function value only, 
    //  1: patch info, function value, and 1st derivs,
    //  2: ... plus 2nd derivatives, 
    //  3: ... plus 3rd derivatives.
    int level;
    Vec2 xy;
    Real f; // f(x,y); valid if level >= 0.

    // These are valid only if level has reached the value indicated.
    Real fx, fy;                    // level >= 1
    Real fxy, fxx, fyy;             // level >= 2 (fyx==fxy)
    Real fxxx, fxxy, fyyy, fxyy;    // level == 3 (fyxx==fxyx==fxxy,
                                    //             fyyx==fyxy==fxyy)
};



//==============================================================================
//                  CLASS BICUBIC SURFACE :: GUTS
//==============================================================================
// Private implementation of BicubicSurface.
class SimTK_SIMMATH_EXPORT BicubicSurface::Guts {
friend class BicubicSurface;
public:
    // Construct an uninitialized BicubicSurface object. This can be filled
    // in later by assignment. Reference count is zero after construction.   
    Guts() {construct();}
    // This class should not be destructed until the reference count drops
    //to zero.
    ~Guts() {assert(referenceCount==0);}

    
    // Implementation of the specified-sample points constructor.
    Guts(const Vector& x, const Vector& y, const Matrix& f, 
         Real smoothness);

    // Implementation the regularly-sampled constructor.
    Guts(const Vec2& XY, const Vec2& spacing, 
         const Matrix& f, Real smoothness);

    // Implementation of the advanced constructor that gives all derivatives
    // at each grid point.
    Guts(const Vector& x, const Vector& y, const Matrix& f, 
         const Matrix& fx, const Matrix& fy, const Matrix& fxy);
    Guts(const Vec2& XY, const Vec2& spacing, const Matrix& f, 
         const Matrix& fx, const Matrix& fy, const Matrix& fxy);

    // Calculate the value of the surface at a particular XY coordinate.
    Real calcValue(const Vec2& XY, PatchHint& hint) const;
    
    // Calculate a partial derivative of this function at a particular XY
    // coordinate.
    Real calcDerivative(const Array_<int>& derivComponents, 
                        const Vec2& XY, PatchHint& hint) const;

    // Calculate a paraboloid and the max/min principal curvatures at a contact
    // point at XY.
    void calcParaboloid
       (const Vec2& XY, PatchHint& hint, Transform& X_SP, Vec2& k) const;
   
    void getNumPatches(int& nx, int &ny) const {
        nx = _ff.nrow()-1;
        ny = _ff.ncol()-1;
    }

    // Need to match out-of-order algebraic coefficients to the required
    // matrix form:
    //      [ a33 a32 a31 a30 ] 
    //  A = [ a23 a22 a21 a20 ] 
    //      [ a13 a12 a11 a10 ]
    //      [ a03 a02 a01 a00 ]
    // in which each entry is a 3-vector,
    // given scalar algebraic coefficients for the z coordinates in this order:
    //   0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15
    //  a00 a10 a20 a30 a01 a11 a21 a31 a02 a12 a22 a32 a03 a13 a23 a33
    //                                        
    // All x,y coefficients will be zero except for the constant and linear
    // terms a00, a01, a10.
    Geo::BicubicHermitePatch
    calcHermitePatch(int i, int j, BicubicSurface::PatchHint& hint) const {
        int nx, ny; getNumPatches(nx,ny);
        SimTK_ERRCHK4_ALWAYS(0<=i && i<nx && 0<=j && j<ny,
            "BicubicSurface::calcHermitePatch()",
            "Patch coordinates (%d,%d) out of range [0..%d,0..%d].",
            i,j,nx-1,ny-1);
        BicubicSurface::PatchHint::Guts& h = hint.updGuts();
        getPatchInfoIfNeeded(i,j,h); // need patch info
        const Vec<16>& a = h.a;
        const Vec3 
        a33(0,0,a[15]), a32(0,0,a[11]), a31(0,  0, a[7]), a30(   0,   0,  a[3]),
        a23(0,0,a[14]), a22(0,0,a[10]), a21(0,  0, a[6]), a20(   0,   0,  a[2]),
        a13(0,0,a[13]), a12(0,0,a[9] ), a11(0,  0, a[5]), a10( h.xS,  0,  a[1]),
        a03(0,0,a[12]), a02(0,0,a[8] ), a01(0,h.yS,a[4]), a00(_x[i],_y[j],a[0]);

        const Mat<4,4,Vec3> A(a33, a32, a31, a30, 
                              a23, a22, a21, a20, 
                              a13, a12, a11, a10,
                              a03, a02, a01, a00);
        return Geo::BicubicHermitePatch(A);
    }

    // Return as a Bezier patch. Currently we just create an algebraic patch
    // which is cheap (~110 flops), and then do an expensive algebraic-to-Bezier 
    // conversion (~220 flops). I don't know of a cheaper approach.
    Geo::BicubicBezierPatch
    calcBezierPatch(int i, int j, BicubicSurface::PatchHint& hint) const {
        const Geo::BicubicHermitePatch hpatch = calcHermitePatch(i,j,hint);
        const Mat<4,4,Vec3>& A = hpatch.getAlgebraicCoefficients();
        const Mat<4,4,Vec3>  B = Geo::BicubicBezierPatch::calcBFromA(A);
        return Geo::BicubicBezierPatch(B);
    }

    // Determine if a point is within the defined surface.
    bool isSurfaceDefined(const Vec2& XY) const;

    // Put a visualization mesh into the given PolygonalMesh argument.
    void createPolygonalMesh(Real resolution, PolygonalMesh& mesh) const;

    int getReferenceCount() const {return referenceCount;}
    void incrReferenceCount() const {++referenceCount;}
    int decrReferenceCount() const 
    {assert(referenceCount);return --referenceCount;}

    // return The total number of elements in the matrix f
    int getFnelt() const {return _ff.nelt();}
    // return The number of elements in the vector X
    int getXnelt() const {return _x.nelt();}
    // return  The number of elements in the vector Y
    int getYnelt() const {return _y.nelt();}
    
    // Methods for debugging and testing

    // Sets a flag that prints a lot of useful debugging data to the screen
    // when this class is used.
    void setDebug(bool aDebug) {_debug = aDebug;}
     
    // return  the matrix ff, which defines the grid points that the surface
    // actually passes through, and derivatives fx,fy,fxy.
    const Matrix_<Vec4>& getff() const {return _ff;}

    // return  the vector x, which defines one side of the mesh grid for 
    // which f(x,y) is defined.
    const Vector& getx() const {return _x;}
    
    // return  the vector y, which defines one side of the mesh grid for 
    // which f(x,y) is defined.
    const Vector& gety() const {return _y;}

    /** Return the function values for the patch containing a particular point.
    @param  XY an (X,Y) location within the grid
    @return  The 16x1 vector that defines the values of f,fx,fy, and fxy at the 
             corners. The vector is given in the following order

             f(0,0)   f(1,0)   f(0,1)   f(1,1),
             fx(0,0)  fx(1,0)  fx(0,1)  fx(1,1),
             fy(0,0)  fy(1,0)  fy(0,1)  fy(1,1),
             fxy(0,0) fxy(1,0) fxy(0,1) fxy(1,1)
    **/
    Vec<16> getPatchFunctionVector(const Vec2& XY) const {
        BicubicSurface::PatchHint hint;
        getFdF(XY,-1,hint); // just need patch info
        const BicubicSurface::PatchHint::Guts& h = hint.getGuts();
        return h.fV;

    }
    
    /** Return the patch coefficients for the patch containing a particular
    point.
    @param  XY an (X,Y) location within the grid
    @return  The 16x1 vector that defines the values of the coefficients aij of 
             the bicubic surface interpolation in this patch. The coefficient 
             vector is given in the following order
             
             a00,a10,a20,a30,a01,a11,a21,a31,a02,a12,a22,a32,a03,a13,a23,a33
    */
    Vec<16> getPatchBicubicCoefficients(const Vec2& XY) const {
        BicubicSurface::PatchHint hint;
        getFdF(XY,-1,hint); // just need patch info
        const BicubicSurface::PatchHint::Guts& h = hint.getGuts();
        return h.a;
    }


private:
    int calcLowerBoundIndex(const Vector& vecV, Real value, int pIdx,
                            int& howResolved) const;
    void getCoefficients(const Vec<16>& f, Vec<16>& aV) const;
    void getFdF(const Vec2& aXY, int wantLevel,
                BicubicSurface::PatchHint& hint) const;
    void getPatchInfoIfNeeded(int x0, int y0, 
                              BicubicSurface::PatchHint::Guts& h) const;

    // This is called from each constructor to initialize this object.
    void construct() {
        referenceCount = 0;
        resetStatistics();
        _hasRegularSpacing = false;
        _debug = false;
    }

    // Return true if the entries in this vector are monotonically increasing
    // (no duplicates allowed).
    static bool isMonotonicallyIncreasing(const Vector& v) {
        for (int i=1; i < v.size(); ++i)
            if (v[i] <= v[i-1]) return false;
        return true;
    }

    // Shared by the irregular and regular-spaced constructors.
    void constructFromSplines
       (const Matrix& f, Real smoothness);
    void constructFromKnownFunction
       (const Matrix& f, const Matrix& fx, const Matrix& fy,
        const Matrix& fxy);

//=============================================================================
// MEMBER VARIABLES
//=============================================================================
private:
    // This class is reference counted and is not destructed until the 
    // reference count goes to zero.
    mutable int referenceCount;

    // Interesting statistics about the use of this surface.
    mutable int numAccesses; 
    mutable int numAccessesSamePoint;
    mutable int numAccessesSamePatch;
    mutable int numAccessesNearbyPatch;
    void resetStatistics() const
    {   numAccesses = numAccessesSamePoint = 0;
    numAccessesSamePatch = numAccessesNearbyPatch = 0; }

    // PROPERTIES
    // Array of values for the independent variables (i.e., the spline knot
    // sequence). Each array must be monotonically increasing, of size nx and 
    // ny elements. These are filled in even if the grid is regularly 
    // spaced, although we'll use the spacing to avoid a search. The values
    // here define the actual grid points; if roundoff causes a calculated
    // grid location to be different, these take priority.
    Vector _x, _y;

    // If the grid has regular spacing we remember the specified spacing
    // here and use it instead of searching when we need to move to a new patch.
    // Note that the grid values are _x[0]+i*spacing[0] and _y[0]+j*spacing[1],
    // but you still have to check in _x and _y to avoid roundoff problems.
    bool _hasRegularSpacing;
    Vec2 _spacing;

    // 2D nx X ny z values that correspond to the values at the grid defined
    // by x and y, and the partial differentials at those grid points. The
    // entries at each grid point are ordered f, fx, fy, fxy.
    enum {F=0, Fx=1, Fy=2, Fxy=3}; 
    Matrix_<Vec4> _ff;

    //A private debugging flag - if set to true, a lot of useful debugging
    //data will be printed tot the screen
    bool _debug;


//=============================================================================
};    // END class BicubicSurface::Guts



}; //namespace
//=============================================================================
//=============================================================================

#endif  // SimTK_BICUBIC_SURFACE_GUTS_H_
