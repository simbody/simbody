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

/**@file
This file contains the library-side implementations of classes
BicubicSurface, BicubicSurface::Guts, and BicubicSurface::PatchHint. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/BicubicSurface.h"
#include "simmath/internal/SplineFitter.h"
#include "simmath/internal/ContactGeometry.h"

#include "BicubicSurface_Guts.h"

#include <algorithm>

using namespace SimTK;
using namespace std;

//==============================================================================
//                              BICUBIC SURFACE
//==============================================================================

// This is just a handle for BicubicSurface::Guts, which is the shared,
// underlying surface representation. Here we just manage the reference
// count and otherwise forward all requests to the Guts class.

// Default constructor is inline; just sets guts pointer to null.

// Destructor: delete guts if this is the last referencing handle.
BicubicSurface::~BicubicSurface() {clear();}

// Copy constructor is shallow.
BicubicSurface::BicubicSurface(const BicubicSurface& source) {
    guts = source.guts; // copy pointer
    if (guts) guts->incrReferenceCount();
}

// Copy assignment is shallow.
BicubicSurface& BicubicSurface::operator=(const BicubicSurface& source) {
    if (guts != source.guts) {
        clear();
        guts = source.guts; // copy pointer
        if (guts) guts->incrReferenceCount();
    }
    return *this;
}

// Return handle to its default-constructed state, deleting guts if this was
// the last reference.
void BicubicSurface::clear() {
    if (guts && guts->decrReferenceCount()==0) 
        delete guts;
    guts = 0;
}

// Constructor for irregularly spaced samples.
BicubicSurface::BicubicSurface
   (const Vector& x, const Vector& y, const Matrix& f, Real smoothness) 
:   guts(0) {
    guts = new BicubicSurface::Guts(x,y,f,smoothness);
    guts->incrReferenceCount(); // will be 1
}

// Constructor for regularly spaced samples.
BicubicSurface::BicubicSurface
   (const Vec2& XY, const Vec2& spacing, const Matrix& f, Real smoothness)
:   guts(0) {
    guts = new BicubicSurface::Guts(XY,spacing,f,smoothness);
    guts->incrReferenceCount(); // will be 1
}

// Constructor from known patch derivatives.
BicubicSurface::BicubicSurface
   (const Vector& x, const Vector& y, const Matrix& f, 
    const Matrix& fx, const Matrix& fy, const Matrix& fxy)
:   guts(0) {
    guts = new BicubicSurface::Guts(x,y,f,fx,fy,fxy);
    guts->incrReferenceCount(); // will be 1
}

// Constructor from known patch derivatives with regular spacing.
BicubicSurface::BicubicSurface
   (const Vec2& XY, const Vec2& spacing, const Matrix& f, 
    const Matrix& fx, const Matrix& fy, const Matrix& fxy)
:   guts(0) {
    guts = new BicubicSurface::Guts(XY,spacing,f,fx,fy,fxy);
    guts->incrReferenceCount(); // will be 1
}

// calcValue(), the fast version.
Real BicubicSurface::calcValue(const Vec2& XY, PatchHint& hint) const {
    SimTK_ERRCHK_ALWAYS(!isEmpty(), "BicubicSurface::calcValue()",
        "This method can't be called on an empty handle.");
    return guts->calcValue(XY, hint); 
}

// calcValue(), the slow version.
Real BicubicSurface::calcValue(const Vec2& XY) const {
    PatchHint hint; // create an empty hint
    return calcValue(XY, hint); 
}

// calcDerivative(), the fast version.
Real BicubicSurface::calcDerivative
   (const Array_<int>& components, const Vec2& XY, PatchHint& hint) const {
    SimTK_ERRCHK_ALWAYS(!isEmpty(), "BicubicSurface::calcDerivative()",
        "This method can't be called on an empty handle.");
    return guts->calcDerivative(components, XY, hint); 
}

// calcDerivative(), the slow version.
Real BicubicSurface::calcDerivative
   (const Array_<int>& components, const Vec2& XY) const {
    PatchHint hint; // create an empty hint
    return calcDerivative(components, XY, hint); 
}

UnitVec3 BicubicSurface::calcUnitNormal(const Vec2& XY, PatchHint& hint) const {
    const Array_<int> dx(1,0), dy(1,1);
    const Real fx = calcDerivative(dx,XY,hint);
    const Real fy = calcDerivative(dy,XY,hint);
    return UnitVec3(-fx, -fy, 1); // (1,0,fx) X (0,1,fy)
}
UnitVec3 BicubicSurface::calcUnitNormal(const Vec2& XY) const {
    PatchHint hint;
    return calcUnitNormal(XY,hint);
}

void BicubicSurface::calcParaboloid
   (const Vec2& XY, PatchHint& hint, Transform& X_SP, Vec2& k) const
{
    SimTK_ERRCHK_ALWAYS(!isEmpty(), "BicubicSurface::calcParaboloid()",
        "This method can't be called on an empty handle.");
    guts->calcParaboloid(XY, hint, X_SP, k);
}

void BicubicSurface::calcParaboloid
   (const Vec2& XY, Transform& X_SP, Vec2& k) const
{
    PatchHint hint;
    calcParaboloid(XY, hint, X_SP, k);
}

void BicubicSurface::getNumPatches(int& nx, int& ny) const 
{   guts->getNumPatches(nx,ny); }

Geo::BicubicHermitePatch BicubicSurface::
calcHermitePatch(int x, int y) const
{   PatchHint hint; return guts->calcHermitePatch(x,y, hint); }

Geo::BicubicBezierPatch BicubicSurface::
calcBezierPatch(int x, int y) const
{   PatchHint hint; return guts->calcBezierPatch(x,y, hint); }

bool BicubicSurface::isSurfaceDefined(const Vec2& XY) const 
{   return getGuts().isSurfaceDefined(XY); }

Vec2 BicubicSurface::getMinXY() const {
    const BicubicSurface::Guts& guts = getGuts();
    return Vec2(guts._x[0],guts._y[0]);
}

Vec2 BicubicSurface::getMaxXY() const {
    const BicubicSurface::Guts& guts = getGuts();
    return Vec2(guts._x[guts._x.size()-1],guts._y[guts._y.size()-1]);
}

PolygonalMesh BicubicSurface::createPolygonalMesh(Real resolution) const {
    PolygonalMesh mesh;
    getGuts().createPolygonalMesh(resolution, mesh);
    return mesh;
}

int BicubicSurface::getNumAccesses() const
{   return getGuts().numAccesses; }

int BicubicSurface::getNumAccessesSamePoint() const
{   return getGuts().numAccessesSamePoint; }

int BicubicSurface::getNumAccessesSamePatch() const
{   return getGuts().numAccessesSamePatch; }

int BicubicSurface::getNumAccessesNearbyPatch() const
{   return getGuts().numAccessesNearbyPatch; }

void BicubicSurface::resetStatistics() const
{   return getGuts().resetStatistics(); }



//==============================================================================
//                          BICUBIC SURFACE :: GUTS
//==============================================================================

// This is the constructor for irregularly-spaced samples.
BicubicSurface::Guts::Guts(const Vector& aX, const Vector& aY, 
                           const Matrix& af, Real smoothness)
{
    construct();

    _x.resize(aX.size()); _y.resize(aY.size());
    _x = aX; _y = aY;

    _hasRegularSpacing = false;

    constructFromSplines(af, smoothness);
}


// This is the constructor for regularly spaced samples.
BicubicSurface::Guts::Guts
   (const Vec2& XY, const Vec2& spacing, const Matrix& af, Real smoothness)
{
    construct();
    
    // Check for reasonable spacing.
    SimTK_ERRCHK2_ALWAYS(spacing > 0,
        "BicubicSurface::BicubicSurface(XY,spacing,af,smoothness)", 
        "A BicubicSurface requires positive spacing in both x and y"
        " but spacing was %g and %g.", spacing[0], spacing[1]);

    const int nx = af.nrow(), ny = af.ncol();
    _x.resize(nx); _y.resize(ny);

    for (int i=0; i < nx; ++i)
        _x[i] = XY[0] + i*spacing[0];
    for (int j=0; j < ny; ++j)
        _y[j] = XY[1] + j*spacing[1];

    _hasRegularSpacing = true;
    _spacing = spacing;

    constructFromSplines(af, smoothness);
}

// This implementation is shared by the regular and irregular constructors.
// We expect _x and _y already to have been filled in with the grid sample
// locations (either as supplied or as generated from regular spacing).
void BicubicSurface::Guts::
constructFromSplines(const Matrix& af, Real smoothness)
{
    const int nx = af.nrow(), ny = af.ncol();

    // Check for sufficient sample size.
    SimTK_ERRCHK2_ALWAYS(af.nrow() >= 4 && af.ncol() >= 4,
        "BicubicSurface::BicubicSurface()", 
        "A BicubicSurface requires at least 4 sample in x and y directions"
        " but grid dimensions were %d and %d.", nx, ny);

    SimTK_ERRCHK4_ALWAYS(_x.size() == nx && _y.size() == ny,
        "BicubicSurface::BicubicSurface()", 
        "Number of samples must match the grid dimension (%d X %d) but"
        "the number of supplied sample points was %d X %d.",
        nx, ny, _x.size(), _y.size());

    // We're checking this even for generated regularly-spaced samples
    // to catch the rare case that the spacing was so small that it 
    // produced two identical sample locations.

    SimTK_ERRCHK_ALWAYS(   isMonotonicallyIncreasing(_x)
                        && isMonotonicallyIncreasing(_y),
        "BicubicSurface::BicubicSurface()",
        "Sample vectors must each be monotonically increasing.");

    // The grid data stores a Vec4 at each grid point with (possibly smoothed)
    // function value f, and derivatives fx, fy, fxy in that order.
    _ff.resize(nx, ny);

    _debug = false;

    // These temporaries are needed for indexing into the Spline functions.
    Vector            coord(1);
    const Array_<int> deriv1(1,0); // just one zero to pick 1st derivative

    // This temporary will hold either the original function values or the
    // smoothed ones.
    Matrix fSmooth(nx,ny);

    // Smoothing strategy: we want something that is symmetric in x and y
    // so that you will get the same surface if you rotate the grid 90 degrees
    // to exchange the meaning of x and y. To accomplish that, we smooth
    // the grid separately along the rows and columns, and then average the
    // results to produce a new grid to which we fit the surface.

    if (smoothness > 0) {
        // This temp holds the function values as they look after smoothing
        // in the x direction (that is, down the columns of constant y).
        Matrix xf(nx,ny);
    
        // Smooth position data along lines of constant y.
        for(int j=0; j < ny; ++j){       
            Spline_<Real> xspline = SplineFitter<Real>::fitForSmoothingParameter
                                            (3,_x,af(j),smoothness).getSpline();
            for(int i=0; i < nx; ++i){    
                coord[0] = _x[i];
                xf(i,j) = xspline.calcValue(coord);
            }
        }

        // Smooth position data along lines of constant x, then average the
        // result with the corresponding value from smoothing the other way. 
        for(int i=0; i < nx; ++i){       
            Spline_<Real> yspline = SplineFitter<Real>::fitForSmoothingParameter
                                           (3,_y,~af[i],smoothness).getSpline();
            for(int j=0; j < ny; ++j){    
                coord[0] = _y[j];
                const Real yfij = yspline.calcValue(coord);
                fSmooth(i,j) = (xf(i,j) + yfij) / 2; // average xf,xy
            }
        }
    } else {
        // Not smoothing.
        fSmooth = af;
    }

    // Now fill in the f and fx entries in our internal grid by exactly
    // fitting a spline to the already-smoothed data.
    for(int j=0; j < ny; ++j){       
        Spline_<Real> xspline = SplineFitter<Real>::fitForSmoothingParameter
                                        (3,_x,fSmooth(j),0).getSpline();
        for(int i=0; i < nx; ++i){    
            coord[0] = _x[i];
            Vec4& fij = _ff(i,j);
            fij[F]  = fSmooth(i,j);
            fij[Fx] = xspline.calcDerivative(deriv1,coord);
        }
    }

    // Compute fy and fxy by fitting splines along the rows.
    // Note that we are using the already-smoothed value of f here.
    Vector tmpRow(ny);
    for(int i=0; i < nx; ++i){
        // Fit splines along rows of constant x to go exactly through the 
        // already-smoothed sample points in order to get fy=Df/Dy.
        Spline_<Real> yspline = SplineFitter<Real>::fitForSmoothingParameter
                                               (3,_y,~fSmooth[i],0).getSpline();

        // Fit splines along rows of constant x to interpolate fx in the y
        // direction to give fxy=Dfx/Dy.
        for (int j=0; j<_y.size(); ++j) tmpRow[j] = _ff(i,j)[Fx];
        Spline_<Real> ydxspline = SplineFitter<Real>::fitForSmoothingParameter
                                                 (3,_y,tmpRow,0).getSpline();

        for(int j=0; j < ny; ++j){    
            coord[0] = _y[j];
            Vec4& fij = _ff(i,j);
            fij[Fy]  = yspline.calcDerivative(deriv1,coord);            
            fij[Fxy] = ydxspline.calcDerivative(deriv1,coord);
        }
    }
}

// This is the advanced constructor where everything is known already.
BicubicSurface::Guts::Guts
   (const Vector& aX, const Vector& aY, const Matrix& af, 
    const Matrix& afx, const Matrix& afy, const Matrix& afxy)
{
    construct();

    _x.resize(aX.size()); _y.resize(aY.size());
    _x = aX; _y = aY;

    _hasRegularSpacing = false;

    constructFromKnownFunction(af, afx, afy, afxy);
}

BicubicSurface::Guts::Guts
   (const Vec2& XY, const Vec2& spacing, const Matrix& af, 
    const Matrix& afx, const Matrix& afy, const Matrix& afxy)
{
    construct();

    // Check for reasonable spacing.
    SimTK_ERRCHK2_ALWAYS(spacing > 0,
        "BicubicSurface::BicubicSurface(XY,spacing,af,smoothness)", 
        "A BicubicSurface requires positive spacing in both x and y"
        " but spacing was %g and %g.", spacing[0], spacing[1]);

    const int nx = af.nrow(), ny = af.ncol();
    _x.resize(nx); _y.resize(ny);

    for (int i=0; i < nx; ++i)
        _x[i] = XY[0] + i*spacing[0];
    for (int j=0; j < ny; ++j)
        _y[j] = XY[1] + j*spacing[1];

    _hasRegularSpacing = true;
    _spacing = spacing;

    constructFromKnownFunction(af, afx, afy, afxy);
}

// Expects _x and _y already to be filled in.
void BicubicSurface::Guts::
constructFromKnownFunction
   (const Matrix& af, const Matrix& afx, const Matrix& afy,
    const Matrix& afxy)
{ 
    const int nx = af.nrow(), ny = af.ncol();

    // Check for sufficient sample size.
    SimTK_ERRCHK2_ALWAYS(nx >= 2 && ny >= 2,
        "BicubicSurface::BicubicSurface(f,fx,fy,fxy)", 
        "A BicubicSurface requires at least 2 sample in x and y directions"
        " but grid dimensions were %d and %d.", nx, ny);

    SimTK_ERRCHK4_ALWAYS(_x.size() == nx && _y.size() == ny,
        "BicubicSurface::BicubicSurface(f,fx,fy,fxy)", 
        "Number of samples must match the grid dimension (%d X %d) but"
        "the number of supplied sample points was %d X %d.",
        nx, ny, _x.size(), _y.size());

    SimTK_ERRCHK2_ALWAYS(   afx.nrow()  == nx && afx.ncol()  == ny
                         && afy.nrow()  == nx && afy.ncol()  == ny
                         && afxy.nrow() == nx && afxy.ncol() == ny,
        "BicubicSurface::BicubicSurface(f,fx,fy,fxy)", 
        "All the derivative sample matrices must match the grid dimension"
        " (%d X %d).", nx, ny);

    SimTK_ERRCHK_ALWAYS(   isMonotonicallyIncreasing(_x)
                        && isMonotonicallyIncreasing(_y),
        "BicubicSurface::BicubicSurface(f,fx,fy,fxy)",
        "Sample vectors must each be monotonically increasing.");

    _ff.resize(nx,ny);

    _debug = false;

    //Copy the data into our packed data structure.
    for(int i=0; i < nx; ++i) {
        for(int j=0; j < ny; ++j){
            Vec4& fij = _ff(i,j);
            fij[F]    = af(i,j);      
            fij[Fx]   = afx(i,j);
            fij[Fy]   = afy(i,j);
            fij[Fxy]  = afxy(i,j);
        }
    }
}

//_____________________________________________________________________________

Real BicubicSurface::Guts::calcValue(const Vec2& aXY, PatchHint& hint) const
{    
    getFdF(aXY,0,hint); // just function value
    const PatchHint::Guts& h = hint.getGuts();
    assert(h.xy == aXY && h.level >= 0);
    return h.f;
}


Real BicubicSurface::Guts::calcDerivative
   (const Array_<int>& aDerivComponents, const Vec2& aXY, PatchHint& hint) const
{
    const int wantLevel = (int)aDerivComponents.size();

    if (wantLevel == 0)
        return calcValue(aXY, hint);  // "0th" deriv is the function value

    for (int i=0; i < wantLevel; ++i) {
        SimTK_ERRCHK2_ALWAYS(aDerivComponents[i]==0 || aDerivComponents[i]==1,
            "BicubicSurface::calcDerivative()",
            "Component %d was %d but must be 0 or 1 for x or y.",
            i, aDerivComponents[i]);
    }

    if (wantLevel > 3)
        return 0;   // 4th and higher derivatives are all zero

    getFdF(aXY, wantLevel, hint);
    const PatchHint::Guts& h = hint.getGuts();
    assert(h.xy == aXY && h.level >= wantLevel);

    if (aDerivComponents.size() == 1)
        return aDerivComponents[0]==0 ? h.fx : h.fy;            // fx : fy

    if (aDerivComponents.size() == 2) {
        if (aDerivComponents[0]==0) //x
            return aDerivComponents[1]==0 ? h.fxx : h.fxy;      // fxx:fxy
        else //y (fyx==fxy)
            return aDerivComponents[1]==0 ? h.fxy : h.fyy;      // fyx:fyy
    }

    // Third derivative.
    if (aDerivComponents[0]==0) { //x
        if (aDerivComponents[1]==0) // xx
            return aDerivComponents[2]==0 ? h.fxxx : h.fxxy;    // fxxx:fxxy
        else // xy (fxyx==fxxy)
            return aDerivComponents[2]==0 ? h.fxxy : h.fxyy;    // fxyx:fxyy
    } else { //y (fyx==fxy)
        if (aDerivComponents[1]==0) // yx (fyxx==fxxy, fyxy==fxyy)
            return aDerivComponents[2]==0 ? h.fxxy : h.fxyy;    // fyxx:fyxy
        else // yy (fyyx==fxyy)
            return aDerivComponents[2]==0 ? h.fxyy : h.fyyy;    // fyyx:fyyy
    }
}

// Cost is patch evaluation + about 200 flops.
void BicubicSurface::Guts::calcParaboloid
   (const Vec2& aXY, PatchHint& hint, Transform& X_SP, Vec2& k) const
{
    getFdF(aXY, 2, hint); // calculate through 2nd derivatives
    const PatchHint::Guts& h = hint.getGuts();
    assert(h.xy == aXY && h.level >= 2);

    const Vec3 P = aXY.append1(h.f);
    const Vec3 dPdx(1,0,h.fx);
    const Vec3 dPdy(0,1,h.fy);
    const UnitVec3 nn(-h.fx, -h.fy, 1); // dPdx X dPdy (normalizing, ~35 flops)
    const Vec3 d2Pdx2(0,0,h.fxx);
    const Vec3 d2Pdy2(0,0,h.fyy);
    const Vec3 d2Pdxdy(0,0,h.fxy);

    // TODO: could save a little time here by taking advantage of the known
    // sparsity of these vectors (probably not worth the trouble).
    k = ContactGeometry::evalParametricCurvature
                                (P,nn,dPdx,dPdy,d2Pdx2,d2Pdy2,d2Pdxdy,X_SP);
}

bool BicubicSurface::Guts::isSurfaceDefined(const Vec2& XYval) const
{
    const bool valueDefined = 
            (_x[0] <= XYval[0] &&  XYval[0] <= _x[_x.size()-1])
        &&  (_y[0] <= XYval[1] &&  XYval[1] <= _y[_y.size()-1]);

    return valueDefined;
}

/* This function ensures that the given hint contains correct information for
the patch indexed (x0,y0). */
void BicubicSurface::Guts::
getPatchInfoIfNeeded(int x0, int y0, BicubicSurface::PatchHint::Guts& h) const {
    const int x1 = x0+1, y1 = y0+1;

    // Compute Bicubic coefficients only if we're in a new patch
    // else use the old ones, because this is an expensive step!
    if( !(h.x0 == x0 && h.y0 == y0) ) {
        // The hint is no good at all since it is for the wrong patch.
        h.clear();
        h.x0 = x0; h.y0 = y0;

        // Compute the scaling of the new patch. Note that neither patch 
        // dimension can be zero since we don't allow duplicates in x or y.
        h.xS = _x(x1)-_x(x0);
        h.yS = _y(y1)-_y(y0);
        h.ooxS = 1/h.xS; h.ooxS2 = h.ooxS*h.ooxS; h.ooxS3=h.ooxS*h.ooxS2;
        h.ooyS = 1/h.yS; h.ooyS2 = h.ooyS*h.ooyS; h.ooyS3=h.ooyS*h.ooyS2;

        // Form the vector f and multiply Ainv*f to form coefficient vector a.

        const Vec4& f00 = _ff(x0,y0);
        const Vec4& f01 = _ff(x0,y1);
        const Vec4& f10 = _ff(x1,y0);
        const Vec4& f11 = _ff(x1,y1);

        h.fV[0] = f00[F];
        h.fV[1] = f10[F];
        h.fV[2] = f01[F];
        h.fV[3] = f11[F];

        // Can't precalculate these scaled values because the same grid point
        // is used for up to four different patches, each scaled differently.
        h.fV[4] = f00[Fx]*h.xS;
        h.fV[5] = f10[Fx]*h.xS;
        h.fV[6] = f01[Fx]*h.xS;
        h.fV[7] = f11[Fx]*h.xS;
    
        h.fV[8]  = f00[Fy]*h.yS;
        h.fV[9]  = f10[Fy]*h.yS;
        h.fV[10] = f01[Fy]*h.yS;
        h.fV[11] = f11[Fy]*h.yS;

        h.fV[12]  = f00[Fxy]*h.xS*h.yS;
        h.fV[13]  = f10[Fxy]*h.xS*h.yS;
        h.fV[14]  = f01[Fxy]*h.xS*h.yS;
        h.fV[15]  = f11[Fxy]*h.xS*h.yS;

        getCoefficients(h.fV,h.a);
    }
}


/**
This function computes the surface value and all derivatives (because it is
cheap to compute these derivatives once the bicubic interpolation 
coefficients have been solved for). These values are stored in mutable
function members so that if a user repeatedly asks for information about
the same location the stored data is supplied. Additionally, if the user 
is staying within a single patch, the the coefficients required for this 
patch are only computed once, and they are re-used until a location outside
the patch is requested.

@param aXY the X,Y location of interest
@param fV   an empty vector that is set to the values of f,fx,fy and fxy of
            the patch corners
                f(0,0)   f(1,0)   f(0,1)   f(1,1),
                fx(0,0)  fx(1,0)  fx(0,1)  fx(1,1),
                fy(0,0)  fy(1,0)  fy(0,1)  fy(1,1),
                fxy(0,0) fxy(1,0) fxy(0,1) fxy(1,1)
@param aijV an empty vector that is set to the coefficient values of the 
            bicubic surface polynomials for the patch that the current
            XY location resides in. The coefficients are returned in this 
            order

            a00,a10,a20,a30,a01,a11,a21,a31,a02,a12,a22,a32,a03,a13,a23,a33

@param aFdF an empty vector that is set to the 10 values of f and its partial
derivatives at the point XY. These values are stored in the following order:
            f(x,y) fx(x,y)  fy(x,y)  fxy(x,y)  fxx(x,y)  fyy(x,y) 
            fxxx(x,y) fxxy(x,y) fyyy(x,y) fxyy(x,y)
*/
void BicubicSurface::Guts::
getFdF(const Vec2& aXY, int wantLevel, PatchHint& hint) const {
    ++numAccesses; // All surface accesses come through here.

    //0. Check if the surface is defined for the XY value given.
    SimTK_ERRCHK6_ALWAYS(isSurfaceDefined(aXY), 
        "BicubicSurface::getFdF (private fcn)", 
        "BicubicSurface is not defined at requested location (%g,%g)."
        " The surface is valid from x[%g %g], y[%g %g].", aXY(0), aXY(1),
        _x[0], _x[_x.size()-1], _y[0], _y[_y.size()-1]);

    // -1 means just do the patch
    // 0 means patch and function value
    // 1 means add 1st deriv, 2 is 2nd, 3 is 3rd
    assert(-1 <= wantLevel && wantLevel <= 3);

    BicubicSurface::PatchHint::Guts& h = hint.updGuts();

    //1. Check to see if we have already computed values for the requested point.
    if(h.level >= wantLevel && aXY == h.xy){
        ++numAccessesSamePoint;
        return;    
    }

    // Nope. We're at least changing points.
    h.xy = aXY;
    h.level = -1; // we don't know anything about this point

    // We're going to feed calcLowerBoundIndex() our best guess as to the
    // patch this point is on. For regularly-spaced grid points we can find
    // it exactly, except for some possible roundoff. Otherwise the best we
    // can do is supply the current index from the hint.
    int pXidx = h.x0, pYidx = h.y0;
    if (_hasRegularSpacing) {
        pXidx = clamp(0, (int)std::floor((h.xy[0]-_x[0])/_spacing[0]),
                      _x.size()-2); // can't be last index
        pYidx = clamp(0, (int)std::floor((h.xy[1]-_y[0])/_spacing[1]),
                      _y.size()-2);
    }

    // Compute the indices that define the patch containing this value.
    int howResolvedX, howResolvedY;
    const int x0 = calcLowerBoundIndex(_x,aXY[0],pXidx,howResolvedX);
    const int x1 = x0+1;
    const int y0 = calcLowerBoundIndex(_y,aXY[1],pYidx,howResolvedY);
    const int y1 = y0+1;

    // 0->same patch, 1->nearby patch, 2->had to search
    const int howResolved = std::max(howResolvedX, howResolvedY);
    if      (howResolved == 0) ++numAccessesSamePatch;
    else if (howResolved == 1) ++numAccessesNearbyPatch;
 
    // Compute Bicubic coefficients only if we're in a new patch
    // else use the old ones, because this is an expensive step!
    getPatchInfoIfNeeded(x0,y0,h);

    // At this point we know that the hint contains valid patch information,
    // but it contains no valid point information.

    if (wantLevel == -1)
        return; // caller just wanted patch info

    const Vec<16>& a = h.a; // abbreviate for convenience


    // Compute where in the patch we are. This has to
    // be done everytime the location within the patch changes
    // ... which has happened if this code gets executed.

    //--------------------------------------------------------------------------
    // Evaluate function value f (38 flops).

    // 8 flops
    const Real xpt = (aXY(0)-_x(x0))*h.ooxS, xpt2=xpt*xpt, xpt3=xpt*xpt2;
    const Real ypt = (aXY(1)-_y(y0))*h.ooyS, ypt2=ypt*ypt, ypt3=ypt*ypt2;
    // 12 flops
    const Mat44 mx(a[ 0],   a[ 1]*xpt,   a[ 2]*xpt2,   a[ 3]*xpt3,
                   a[ 4],   a[ 5]*xpt,   a[ 6]*xpt2,   a[ 7]*xpt3,
                   a[ 8],   a[ 9]*xpt,   a[10]*xpt2,   a[11]*xpt3,
                   a[12],   a[13]*xpt,   a[14]*xpt2,   a[15]*xpt3);
    // 12 flops
    const Vec4 xsum = mx.rowSum();
    // 6 flops
    h.f = xsum[0] + ypt*xsum[1] + ypt2*xsum[2] + ypt3*xsum[3];
    h.level = 0; // function value is ready
    if (wantLevel == 0)
        return;

    //--------------------------------------------------------------------------
    // Evaluate first derivatives fx, fy (43 flops).

    // fy is 9 flops
    const Real dypt = h.ooyS, dypt2= 2*ypt*h.ooyS, dypt3= 3*ypt2*h.ooyS;
    h.fy = dypt*xsum[1] + dypt2*xsum[2] + dypt3*xsum[3];

    // fx is 34 flops
    const Real dxpt=h.ooxS, dxpt2=2*xpt*h.ooxS, dxpt3=3*xpt2*h.ooxS;
    const Mat43 mdx(a[ 1]*dxpt,    a[ 2]*dxpt2,    a[ 3]*dxpt3,
                    a[ 5]*dxpt,    a[ 6]*dxpt2,    a[ 7]*dxpt3,
                    a[ 9]*dxpt,    a[10]*dxpt2,    a[11]*dxpt3,
                    a[13]*dxpt,    a[14]*dxpt2,    a[15]*dxpt3);
    const Vec4 dxsum = mdx.rowSum();
    h.fx   = dxsum[0] + ypt*dxsum[1] + ypt2*dxsum[2] + ypt3*dxsum[3];
    h.level = 1; // first derivatives are ready
    if (wantLevel == 1)
        return;


    //--------------------------------------------------------------------------
    // Evaluate second derivatives fxy, fxx, fyy (40 flops).

    // fxy, fyy are 11 flops
    h.fxy = dypt*dxsum[1] + dypt2*dxsum[2] + dypt3*dxsum[3];
    const Real dyypt2=2*h.ooyS2, dyypt3=6*ypt*h.ooyS2;
    h.fyy = dyypt2*xsum[2] + dyypt3*xsum[3];

    // fxx is 29 flops
    const Real dxxpt2=2*h.ooxS2, dxxpt3=6*xpt*h.ooxS2;
    const Mat42 mdxx(a[ 2]*dxxpt2,    a[ 3]*dxxpt3,
                     a[ 6]*dxxpt2,    a[ 7]*dxxpt3,
                     a[10]*dxxpt2,    a[11]*dxxpt3,
                     a[14]*dxxpt2,    a[15]*dxxpt3);
    const Vec4 dxxsum = mdxx.rowSum();
    h.fxx  = dxxsum[0] + ypt*dxxsum[1] + ypt2*dxxsum[2] + ypt3*dxxsum[3];
    h.level = 2; // second derivatives are ready
    if (wantLevel == 2)
        return;

    //--------------------------------------------------------------------------
    // Evaluate third derivatives fxxx, fxxy, fyyy, fxyy (21 flops).

    // 10 flops
    const Real dyyypt3=6*h.ooyS3;
    h.fyyy = dyyypt3*xsum[3];
    h.fxyy = dyypt2*dxsum[2] + dyypt3*dxsum[3];
    h.fxxy = dypt*dxxsum[1] + dypt2*dxxsum[2] + dypt3*dxxsum[3];

    // 11 flops
    const Real dxxxpt3=6*h.ooxS3;
    const Vec4 mdxxx(a[ 3]*dxxxpt3,
                     a[ 7]*dxxxpt3,
                     a[11]*dxxxpt3,
                     a[15]*dxxxpt3);
    h.fxxx = mdxxx[0] + ypt* mdxxx[1] + ypt2*mdxxx[2] + ypt3*mdxxx[3];
    h.level = 3; // third derivatives are ready

    if(_debug == true){
        cout<<" getFdF" << endl;
        cout <<"Member variables" <<endl;
        cout << "_x" << _x << endl;
        cout << "\n"<<endl;
        cout << "_y" << _y << endl;
        cout << "\n"<<endl;
        cout << "_ff" << _ff << endl;
        cout << "\n"<<endl;

        cout <<" Intermediate variables " << endl;
        cout <<"XY: " << aXY << endl;
        printf("[x0 x1], [y0 y1]: [%d %d],[%d %d]\n",x0,x1,y0,y1);
        printf("[x0V x1V], [y0V y1V]: [%f %f],[%f %f]\n",_x(x0),_x(x1),_y(y0),_y(y1));
        cout <<" xS " << h.xS << " yS " << h.yS << endl;
        printf("(xp,yp): (%f %f)\n",xpt,ypt);
        cout << "\n\n\n"<<endl;
    }
}


/** Given a search value aVal and an n-vector aVec (n>=2) containing 
monotonically increasing values (no duplicates), this method finds the unique
pair of indices (i,i+1) such that either 
    - aVec[i] <= aVal < aVec[i+1], or
    - aVec[n-2] < aVal == aVec[n-1].

The following preconditions must be satisfied by the arguments to this
method:
    - aVec.size() >= 2
    - aVal >= aVec[0]
    - aVal <= aVec[n-1]


aVec = [0.1 0.2 0.3 0.4 0.5];
aVal = 0.125
idxLB = calcLowerBoundIndex(aVec,aVal,-1);
 
Then idxLB should be 0. Some effort has been put into making this code 
efficient, as  it is expected that very large vectors could be used in this 
function. If the data is evenly spaced, the index is computed. If not data 
near a previous index (set by the user) is searched. If the data is still 
not found a binary search is performed

@params aVec: A SimTK vector containing monotonically increasing data
@params aVal: A value bounded by the numbers in the vector
@params pIdx: A hint index provided by the user (of this function)
@returns int: The index of the entry in the vector that is the as close to
            aVal without exceeding it.

*/

// Return true if aVal is on the patch (really the line) between 
// aVec[indxL] and aVec[indxL+1]. By "on the patch" we mean that
// aVec[indxL] <= aVal < aVec[indxL+1] unless this is the last patch in
// which case we allow aVec[indxL] < aVal <= aVec[indxL+1].
static bool isOnPatch(const Vector& aVec, int indxL, Real aVal) {
    assert(aVec.size() >= 2);
    const int maxLB = aVec.size() - 2;
    assert(0 <= indxL && indxL <= maxLB);

    const Real low = aVec[indxL], high = aVec[indxL+1];

    if (aVal < low || aVal > high)
        return false;

    // Here we know low <= aVal <= high.
    if (aVal < high)
        return true;

    // Here we know low < aVal == high. This is only allowed on the last patch.
    return indxL == maxLB;
}

// howResolved: 0->same patch, 1->nearby patch, 2->search
int BicubicSurface::Guts::
calcLowerBoundIndex(const Vector& aVec, Real aVal, int pIdx,
                    int& howResolved) const {
    int idxLB = -1;
    bool idxComputed = false;

    assert(aVec.size() >= 2);

    // Because we're trying to find the lower index, it can't be the very
    // last knot.
    const int maxLB = aVec.size() - 2;
        
    // 1. Do a local search around the previous index, if one is given.
    if(0 <= pIdx && pIdx <= maxLB) {    
        // Are we still on the same patch? Caution -- can't be equal to the
        // upper knot unless it is the last one.
        if (isOnPatch(aVec, pIdx, aVal)) {
            howResolved = 0;
            return pIdx;
        }

        // Not on the same patch, how about adjacent patches?
        if (aVal < aVec[pIdx]) {
            // Value moved below the current patch.
            if (pIdx > 0 && isOnPatch(aVec, pIdx-1, aVal)) {
                howResolved = 1;
                return pIdx-1; // found it next door!
            }
        } else if (aVal >= aVec[pIdx+1]) {
            // Value moved above the current patch.
            if (pIdx < maxLB && isOnPatch(aVec, pIdx+1, aVal)) {
                howResolved = 1;
                return pIdx+1; // found it next door!
            }
        }    
    }
    
    // Either we didn't have a previous index to try, or it didn't help.

    // 2. Check the end points. We'll count these as "nearby patches" since
    // they are about the same amount of work.
    if (aVal <= aVec[0]) {
        howResolved = 1;
        return 0;
    }
    if (aVal >= aVec[maxLB+1])  {
        howResolved = 1;
        return maxLB;
    }
        
    // 3. If still not found, use binary search to find the appropriate index.
    
    // std::upper_bound returns the index of the first element that
    // is strictly greater than aVal (one after the last element
    // if aVal is exactly equal to the last knot).
    const Real* upper = 
        std::upper_bound(&aVec[0], &aVec[0] + aVec.size(), aVal);
    const int upperIx = clamp(0, (int)(upper-&aVec[1]), maxLB);

    howResolved = 2;
    return upperIx;
}

/**
This function will compute the coefficients aij for a bicubic patch that has
the values of f, fx, fy and fxy defined at its four corners. This function is
relatively expensive, as it involves multiplying a 16x16 matrix by a 16x1,
so it is called only when absolutely necessary.

@param f  Vector f defining the values of f, fx, fy and fxy at the four 
            corners of the patch. Vector f is in this order:
            f(0,0)   f(1,0)   f(0,1)   f(1,1),
            fx(0,0)  fx(1,0)  fx(0,1)  fx(1,1),
            fy(0,0)  fy(1,0)  fy(0,1)  fy(1,1),
            fxy(0,0) fxy(1,0) fxy(0,1) fxy(1,1)
@param aV  An empty vector for which the coefficients aij are written into
                in the following order:

                a00,a10,a20,a30,a01,a11,a21,a31,a02,a12,a22,a32,a03,a13,a23,a33
*/
void BicubicSurface::Guts::
getCoefficients(const Vec<16>& fV, Vec<16>& aV) const {
    // This is what the full matrix inverse looks like (copied here from
    // Wikipedia). It is very sparse and contains only a few unique values
    // so the matrix-vector product can be done very cheaply if worked out
    // in painstaking detail (with some help from Maple). The full 
    // matrix-vector product takes 496 flops; we can do it in 80.
    /*
    const Real Ainv[] = { 
        1, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0, 
        0, 0, 0, 0,     1, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0, 
       -3, 3, 0, 0,    -2,-1, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0, 
        2,-2, 0, 0,     1, 1, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0, 
        0, 0, 0, 0,     0, 0, 0, 0,     1, 0, 0, 0,     0, 0, 0, 0, 
        0, 0, 0, 0,     0, 0, 0, 0,     0, 0, 0, 0,     1, 0, 0, 0, 
        0, 0, 0, 0,     0, 0, 0, 0,    -3, 3, 0, 0,    -2,-1, 0, 0, 
        0, 0, 0, 0,     0, 0, 0, 0,     2,-2, 0, 0,     1, 1, 0, 0, 
       -3, 0, 3, 0,     0, 0, 0, 0,    -2, 0,-1, 0,     0, 0, 0, 0, 
        0, 0, 0, 0,    -3, 0, 3, 0,     0, 0, 0, 0,    -2, 0,-1, 0, 
        9,-9,-9, 9,     6, 3,-6,-3,     6,-6, 3,-3,     4, 2, 2, 1, 
       -6, 6, 6,-6,    -3,-3, 3, 3,    -4, 4,-2, 2,    -2,-2,-1,-1, 
        2, 0,-2, 0,     0, 0, 0, 0,     1, 0, 1, 0,     0, 0, 0, 0, 
        0, 0, 0, 0,     2, 0,-2, 0,     0, 0, 0, 0,     1, 0, 1, 0, 
       -6, 6, 6,-6,    -4,-2, 4, 2,    -3, 3,-3, 3,    -2,-1,-2,-1, 
        4,-4,-4, 4,     2, 2,-2,-2,     2,-2, 2,-2,     1, 1, 1, 1 
    };
    Mat<16,16> AinvM(Ainv);
    aV = AinvM*fV; // So cool that I can do this in C++! Go Sherm!
    */


    // Matt's masterful Maple work, plus a little manual hacking by Sherm
    // produced this version:

    // Caution: note index change.
    Real f1=fV[0], f2=fV[1], f3=fV[2], f4=fV[3], f5=fV[4], f6=fV[5],
         f7=fV[6], f8=fV[7], f9=fV[8], f10=fV[9], f11=fV[10], f12=fV[11],
         f13=fV[12], f14=fV[13], f15=fV[14], f16=fV[15];

    // 54 add/subtract
    // 26 multiplies
    Real f86 = f8-f6, f1211=f12-f11, f109=f10-f9, f75=f7-f5;
    Real f1234 = f1-f2-f3+f4;
    Real t312 = 2*f86;
    Real t311 = 2*f1211;
    Real t310 = 3*f86 - 2*f14;
    Real t309 = 3*f1211 - 2*f15;
    Real t10 = 2*f9;
    Real t308 = f14 - 2*f10 + t10;
    Real t307 = 3*f109 - f14;
    Real t306 = 3*f75 - f15;
    Real t15 = 2*f5;
    Real t305 = t15 - 2*f7 + f13 + f15;
    Real t289 = -2*f13;
    Real t304 = t289 - 6*f1234 - f16;
    Real t302 = -3*f1;
    Real t301 = 2*f1;
    aV(0) = f1;
    aV(1) = f5;
    aV(2) = t302 + 3*f2 - t15 - f6;
    aV(3) = t301 - 2*f2 + f5 + f6;
    aV(4) = f9;
    aV(5) = f13;
    aV(6) = t289 + t307;
    aV(7) = f13 + t308;
    aV(8) = t302 + 3*f3 - t10 - f11;
    aV(9) = t289 + t306;
    aV(10) = 9*f1234 + 6*(f5-f7+f9-f10) + 4*f13 + f16 - t309 - t310;
    aV(11) = 4*f109 + t304 + t306 + t310 + t311;
    aV(12) = t301 - 2*f3 + f9 + f11;
    aV(13) = t305;
    aV(14) = 4*f75  + t304 + t307 + t309 + t312;
    aV(15) = 4*f1234 + f16 + t305 + t308 - t311 - t312;

    if(_debug == true){
        cout << "getCoefficients" << endl;
        cout << "\tfV" << fV << endl;
        cout << "\taijV" << aV << endl;
    }
}

// We're going to generate quads in strips like this:
//    *  <- *  <-  *  <-   *
//    v    ^ v    ^ v      ^
//    * ->  *  ->  *  ->   *
// Each patch will generate the same number of quads so there will be more
// quads where the surface is denser.
void BicubicSurface::Guts::
createPolygonalMesh(Real resolution, PolygonalMesh& mesh) const {
    PatchHint hint;
    // Number of patches in x and y direction.
    const int nxpatch = _x.size()-1;
    const int nypatch = _y.size()-1;
    // n is the number of subdivisions per patch
    const int n = std::max(1 +  (int)(resolution+.5), 1); // round
    const int nx = nxpatch*n + 1; // number of vertices along each row
    const int ny = nypatch*n + 1;
    // These will alternately serve as previous and current.
    Array_<int> row1(ny), row2(ny);
    Array_<int>* prevRowVerts = &row1;
    Array_<int>* curRowVerts  = &row2;
    Vector xVals(nx), yVals(ny);

    // Calculate the x and y sample values.
    int nxt = 0;
    for (int px=0; px < nxpatch; ++px) {
        const Real xlo = _x[px], xhi = _x[px+1];
        const Real xs = xhi - xlo, width = xs/n;
        // For each quad within patch px
        for (int qx=0; qx < n; ++qx)
            xVals[nxt++] = xlo + qx*width;
    }
    xVals[nxt++] = _x[_x.size()-1];
    assert(nxt == nx);

    nxt = 0;
    for (int py=0; py < nypatch; ++py) {
        const Real ylo = _y[py], yhi = _y[py+1];
        const Real ys = yhi - ylo, width = ys/n;
        // For each quad within patch py
        for (int qx=0; qx < n; ++qx)
            yVals[nxt++] = ylo + qx*width;
    }
    yVals[nxt++] = _y[_y.size()-1];
    assert(nxt == ny);

    // Fill in the zeroth row.
    Vec2 pt; pt[0] = xVals[0];
    for (int j=0; j < ny; ++j) {
        pt[1] = yVals[j];
        const Real z = calcValue(pt, hint);
        (*prevRowVerts)[j] = mesh.addVertex(Vec3(pt[0],pt[1],z));
    }

    // For each remaining row, generate vertices and then a strip of faces.
    Array_<int> face(4);
    for (int i=1; i < nx; ++i) {
        Vec2 pt; pt[0] = xVals[i];
        for (int j=0; j < ny; ++j) {
            pt[1] = yVals[j];
            const Real z = calcValue(pt, hint);
            (*curRowVerts)[j] = mesh.addVertex(Vec3(pt[0],pt[1],z));
        }
        for (int j=1; j < ny; ++j) {
            face[0] = (*curRowVerts)[j-1]; // counterclockwise
            face[1] = (*curRowVerts)[j];
            face[2] = (*prevRowVerts)[j];
            face[3] = (*prevRowVerts)[j-1];
            mesh.addFace(face);
        }
        std::swap(prevRowVerts, curRowVerts);
    }
}

//==============================================================================
//                     BICUBIC SURFACE :: PATCH HINT
//==============================================================================

// Here we ensure that there is always a valid PatchHint::Guts object present
// and forward all operations to it.

BicubicSurface::PatchHint::PatchHint() : guts(new PatchHint::Guts) {}
BicubicSurface::PatchHint::~PatchHint() {delete guts;}

BicubicSurface::PatchHint::PatchHint(const PatchHint& src)
:   guts(new PatchHint::Guts(*src.guts)) {}

BicubicSurface::PatchHint&
BicubicSurface::PatchHint::operator=(const PatchHint& src) 
{   *guts = *src.guts; return *this; }

bool BicubicSurface::PatchHint::isEmpty() const {return guts->isEmpty();}
void BicubicSurface::PatchHint::clear()         {guts->clear();}




