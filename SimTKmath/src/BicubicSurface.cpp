/* -------------------------------------------------------------------------- *
 *                        SimTK Simbody: SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2011 Stanford University and the Authors.           *
 * Authors: Matthew Millard                                                   *
 * Contributors: Michael Sherman                                              *
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

/**@file
This file contains the library-side implementations of classes
BicubicSurface, BicubicSurface::Guts, and BicubicSurface::PatchHint. **/

#include "SimTKCommon.h" 
#include "simmath/internal/common.h"
#include "simmath/internal/BicubicSurface.h"
#include "simmath/internal/SplineFitter.h"

#include "BicubicSurface_Guts.h"

#include <algorithm>

using namespace SimTK;
using namespace std;

//==============================================================================
//                              BICUBIC SURFACE
//==============================================================================

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
   (Real xSpacing, Real ySpacing, const Matrix& f, Real smoothness) 
:   guts(0) {
    guts = new BicubicSurface::Guts(xSpacing,ySpacing,f,smoothness);
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

// calcValue(), the fast version.
Real BicubicSurface::calcValue(const Vec2& XY, PatchHint& hint) const {
    SimTK_ERRCHK_ALWAYS(!isEmpty(), "BicubicSurface::calcValue()",
        "This method can't be called on an empty handle.");
    return guts->calcValue(XY, hint); 
}

// calcValue(), the slow version.
Real BicubicSurface::calcValue(const Vec2& XY) const {
    SimTK_ERRCHK_ALWAYS(!isEmpty(), "BicubicSurface::calcValue()",
        "This method can't be called on an empty handle.");
    PatchHint hint; // create an empty hint
    return guts->calcValue(XY, hint); 
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
    SimTK_ERRCHK_ALWAYS(!isEmpty(), "BicubicSurface::calcDerivative()",
        "This method can't be called on an empty handle.");
    PatchHint hint; // create an empty hint
    return guts->calcDerivative(components, XY, hint); 
}



//==============================================================================
//                          BICUBIC SURFACE :: GUTS
//==============================================================================

BicubicSurface::Guts::Guts(const Vector& aX, const Vector& aY, 
                           const Matrix& af, Real smoothness)
{
    construct();
    
    // CHECK NUMBER OF DATA POINTS
    SimTK_ERRCHK_ALWAYS( ((aX.size() >= 4 && aY.size() >= 4)
                      && ( af.ncol() >= 4 && af.nrow() >= 4) ), 
        "BicubicSurface::BicubicSurface", 
        "A BicubicSurface requires aX and aY to be of length 4,"
        " and af to be 4x4");

    // CHECK DIMENSIONS OF AF
    SimTK_ERRCHK_ALWAYS((af.ncol() == aX.size()
                      && af.nrow() == aY.size()), 
        "BicubicSurface::BicubicSurface", 
        "Matrix f(x,y) must have a row dimension that matches the size of "
        "\n vector X, and a column dimension that matches the size of "
        "\n vector Y");


    // INDEPENDENT VALUES (KNOT SEQUENCE)
    _x.resize(aX.size());
    _y.resize(aY.size());
    _f.resize(af.nrow(),af.ncol());

    _fx.resize(af.nrow(),af.ncol());
    _fy.resize(af.nrow(),af.ncol());
    _fxy.resize(af.nrow(),af.ncol());

    _flagXEvenlySpaced = true;
    _flagYEvenlySpaced = true;
    _debug = false;

    Matrix tmpf(af.nrow(),af.ncol());
    Vector tmpfcol(_x.size());
    Vector tmpVal(1);

    Real xsp = aX(1)-aX(0);
    Real ysp = aY(1)-aY(0);
    Real xspi=0.0;
    Real yspj=0.0;
    
    //Interpolate position data according to the smoothness factor.
    for(int j=0; j<_y.size();j++){
        tmpfcol = af(j);        
        Spline_<Real> xspline = SplineFitter<Real>::fitForSmoothingParameter(3,aX,tmpfcol,smoothness).getSpline();
        for(int i=0; i<_x.size(); i++){    
            tmpVal(0) = aX(i);
            tmpf(i,j) = xspline.calcValue(tmpVal);
        }
    }

    //Copy the data over, check for even spacing.
    for(int i=0; i<_x.size();i++)
    {
        _x(i) = aX(i);
        if(i > 0){ 
            xspi = _x(0) + i*xsp;//aX[i]-aX[i-1];
            if(std::abs(xspi-_x(i))/xsp > 1e-6) 
                _flagXEvenlySpaced = false;
        }

        for(int j=0; j<_y.size();j++)
            _f(i,j) = tmpf(i,j);                        
    }

    for(int j=0; j<_y.size();j++)
    {
        _y(j) = aY(j);
        if(j > 0){ 
            yspj = _y(0) + j*ysp;//yspj = aY[j]-aY[j-1];
            if(std::abs(yspj-_y(j))/ysp > 1e-6) 
                _flagYEvenlySpaced = false;
        }
    }

    Vector tmpVX(_x.size());
    RowVector tmpVY(_y.size());
    RowVector tmpVXY(_y.size());
    Vector tmpdFV;    
    Vector tmpdFXYV;
    Array_<int, unsigned int> deriv1(1);
    deriv1[0]=0;

    //Compute fx using NaturalCubicSplines
    //Future upgrade: make the type of spline user selectable
    //tmpdFV.resize(_x.size());
    
    tmpVal.resize(1);

    for(int j=0; j<_y.size();j++){
        tmpVX = _f(j);        
        Spline_<Real> xspline = 
            SplineFitter<Real>::fitForSmoothingParameter(3,_x,tmpVX,0.0).getSpline();
        for(int i=0; i<_x.size(); i++){    
            tmpVal(0) = _x(i);
            _fx(i,j) = xspline.calcDerivative(deriv1,tmpVal);
        }
    }

    //Compute fy and fxy using NaturalCubicSplines
    //Future upgrade: make the type of spline user selectable
    tmpdFV.resize(_y.size());
    tmpdFXYV.resize(_y.size());

    for(int i=0; i<_x.size();i++){
        tmpVY = _f[i]; //square brackets grabs a row
        tmpVXY = _fx[i];
        Spline_<Real> yspline   = 
            SplineFitter<Real>::fitForSmoothingParameter(3,_y,~tmpVY,0.0).getSpline();        
        Spline_<Real> ydxspline = 
            SplineFitter<Real>::fitForSmoothingParameter(3,_y,~tmpVXY,0.0).getSpline();


        for(int j=0; j<_y.size(); j++){    
            tmpVal(0) = _y(j);            
            _fy(i,j)    =   yspline.calcDerivative(deriv1,tmpVal);            
            _fxy(i,j)    = ydxspline.calcDerivative(deriv1,tmpVal);
        }
    }

    //These are mutable so I can set them in the const functions that the user
    //uses to call calcValue, and calcDerivative
    _pXYVal     = Vector(2);
    _pXYIdx     = Array_<int>(4);
    _pFdF       = Vector(10);
    _fV         = Vector(16);
    _aV         = Vector(16);

    //Initialize the record data with values that will not show up in any
    //user's mesh grid.
    _pXYVal = NaN;
    _pXYIdx[0] = -1;
    _pXYIdx[1] = -1;
    _pXYIdx[2] = -1;
    _pXYIdx[3] = -1;
    _pFdF = NaN;
    _fV = NaN;
    _aV =  NaN;
}

BicubicSurface::Guts::Guts
   (const Vector& aX, const Vector& aY, const Matrix& af, 
    const Matrix& afx, const Matrix& afy, const Matrix& afxy)
{
    construct();

    // CHECK NUMBER OF DATA POINTS
    SimTK_ERRCHK_ALWAYS((aX.size() >= 4     && aY.size() >= 4) 
                     && (af.ncol() >= 4     && af.nrow() >= 4)
                     && (afx.ncol() >= 4    && afx.nrow() >= 4)
                     && (afy.ncol() >= 4    && afy.nrow() >= 4)
                     && (afxy.ncol() >= 4   && afxy.nrow() >= 4), 
        "BicubicSurface::BicubicSurface", 
        "A BicubicSurface requires aX and aY to be of length 4, and af to"
        " be 4x4");


    // CHECK DIMENSIONS OF AF
    SimTK_ERRCHK_ALWAYS((( af.nrow() == aX.size() &&   af.ncol() == aY.size())
                     &&  (afx.nrow() == aX.size() &&  afx.ncol() == aY.size())
                     &&  (afy.nrow() == aX.size() &&  afy.ncol() == aY.size())
                     && (afxy.nrow() == aX.size() && afxy.ncol() == aY.size())), 
        "BicubicSurface::BicubicSurface", 
        "Matrix f,fx,fy,fxy must have row dimensions that match the length"
        "\nof vector X, and column dimensions that match the length of "
        "\nof vector Y.");

    // INDEPENDENT VALUES (KNOT SEQUENCE)
    _x.resize(aX.size());
    _y.resize(aY.size());
    _f.resize(af.nrow(),af.ncol());

    _fx.resize(af.nrow(),af.ncol());
    _fy.resize(af.nrow(),af.ncol());
    _fxy.resize(af.nrow(),af.ncol());

    _flagXEvenlySpaced = true;
    _flagYEvenlySpaced = true;
    _debug = false;

    Matrix tmpf(af.nrow(),af.ncol());
    Vector tmpfcol(_x.size());
    Vector tmpVal(1);

    Real xsp = aX(1)-aX(0);
    Real ysp = aY(1)-aY(0);
    Real xspi=0.0;
    Real yspj=0.0;

    //Copy the data over, check for even spacing.
    for(int i=0; i<_x.size();i++)
    {
        _x(i) = aX(i);
        if(i > 0){ 
            xspi = _x(0) + i*xsp;
            if(std::abs(xspi-_x(i))/xsp > 1e-6) 
                _flagXEvenlySpaced = false;
        }

        for(int j=0; j<_y.size();j++){
            _f(i,j)     = af(i,j);      
            _fx(i,j)    = afx(i,j);
            _fy(i,j)    = afy(i,j);
            _fxy(i,j)   = afxy(i,j);
        }
    }

    for(int j=0; j<_y.size();j++)
    {
        _y(j) = aY(j);
        if(j > 0){ 
            yspj = _y(0) + j*ysp;
            if(std::abs(yspj-_y(j))/ysp > 1e-6) 
                _flagYEvenlySpaced = false;
        }
    }

    //These are mutable so I can set them in the const functions that the user
    //uses to call calcValue, and calcDerivative
    _pXYVal     = Vector(2);
    _pXYIdx     = Array_<int>(4);
    _pFdF       = Vector(10);
    _fV         = Vector(16);
    _aV         = Vector(16);

    //Initialize the record data with values that will not show up in any
    //user's mesh grid.
    _pXYVal = NaN;
    _pXYIdx[0] = -1;
    _pXYIdx[1] = -1;
    _pXYIdx[2] = -1;
    _pXYIdx[3] = -1;
    _pFdF = NaN;
    _fV = NaN;
    _aV =  NaN;
}
//_____________________________________________________________________________

Real BicubicSurface::Guts::calcValue(const Vec2& aXY, PatchHint& hint) const
{    
    Vector fV, aV, aFdF; 
    getFdF(aXY,fV,aV,aFdF, hint);
    return aFdF(0);
}


Real BicubicSurface::Guts::calcDerivative
   (const Array_<int>& aDerivComponents, const Vec2& aXY, PatchHint& hint) const
{
    if (aDerivComponents.empty())
        return calcValue(aXY, hint);  // "0th" deriv is the function value

    for (int i=0; i < (int)aDerivComponents.size(); ++i) {
        SimTK_ERRCHK2_ALWAYS(aDerivComponents[i]==0 || aDerivComponents[i]==1,
            "BicubicSurface::calcDerivative()",
            "Component %d was %d but must be 0 or 1 for x or y.",
            i, aDerivComponents[i]);
    }

    if (aDerivComponents.size() > 3)
        return 0;   // 4th and higher derivatives are all zero

    Vector fV, aV, aFdF; 
    getFdF(aXY,fV,aV,aFdF, hint);
    // 0=f, 1=fx, 2=fy, 3=fxy, 4=fxx, 5=fyy, 6=fxxx, 7=fxxy, 8=fyyy, 9=fxyy
    //                   =fyx                         =fyxx           =fyyx
    //                                                =fxyx           =fyxy

    if (aDerivComponents.size() == 1)
        return aDerivComponents[0]==0 ? aFdF[1] : aFdF[2];      // fx : fy

    if (aDerivComponents.size() == 2)
        if (aDerivComponents[0]==0) //x
            return aDerivComponents[1]==0 ? aFdF[4] : aFdF[3];  // fxx:fxy
        else //y (fyx==fxy)
            return aDerivComponents[1]==0 ? aFdF[3] : aFdF[5];  // fyx:fyy

    // Third derivative.
    if (aDerivComponents[0]==0) { //x
        if (aDerivComponents[1]==0) // xx
            return aDerivComponents[2]==0 ? aFdF[6] : aFdF[7];  // fxxx:fxxy
        else // xy (fxyx==fxxy)
            return aDerivComponents[2]==0 ? aFdF[7] : aFdF[9];  // fxyx:fxyy
    } else { //y (fyx==fxy)
        if (aDerivComponents[1]==0) // yx (fyxx==fxxy, fyxy==fxyy)
            return aDerivComponents[2]==0 ? aFdF[7] : aFdF[9];  // fyxx:fyxy
        else // yy (fyyx==fxyy)
            return aDerivComponents[2]==0 ? aFdF[9] : aFdF[8];  // fyyx:fyyy
    }
}

bool BicubicSurface::Guts::isSurfaceDefined(const Vec2& XYval) const
{
    const bool valueDefined = 
            (_x[0] <= XYval[0] &&  XYval[0] <= _x[_x.size()-1])
        &&  (_y[0] <= XYval[1] &&  XYval[1] <= _y[_y.size()-1]);

    return valueDefined;
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
getFdF(const Vec2& aXY, Vector& fV, Vector& aijV, Vector& aFdF,
       PatchHint& hint) const
{
    //0. Check if the surface is defined for the XY value given
    //   Check if desired point is inside the grid, else throw an exception
    SimTK_ERRCHK6_ALWAYS(isSurfaceDefined(aXY), 
        "BicubicSurface::getFdF (private fcn)", 
        "BicubicSurface is not defined at requested location (%g,%g)."
        " The surface is valid from x[%g %g], y[%g %g].", aXY(0), aXY(1),
        _x[0], _x[_x.size()-1], _y[0], _y[_y.size()-1]);

    Real val = 0;
    int pXidx = -1;
    int pYidx = -1;
    aFdF.resize(10);

    //1. Check to see if we have already computed values for the requested point.
    if(aXY(0) == _pXYVal(0) && aXY(1) == _pXYVal(1)){
        aFdF    = _pFdF;
        fV      = _fV;
        aijV    = _aV;
        return;    
    }

    // Nope.
    //1. Compute the indices that define the patch that the value is in    
    int x0 = calcLowerBoundIndex(_x,aXY(0),pXidx,_flagXEvenlySpaced);
    int x1 = x0+1;
    int y0 = calcLowerBoundIndex(_y,aXY(1),pYidx,_flagYEvenlySpaced);
    int y1 = y0+1;

    //2. Form the vector f

    // Compute the scaling of the local patch. Note that neither patch 
    // dimension can be zero since we don't allow duplicates in x or y.
    const Real xS = _x(x1)-_x(x0);
    const Real yS = _y(y1)-_y(y0);
    const Real ooxS = 1/xS, ooxS2 = ooxS*ooxS, ooxS3=ooxS*ooxS2;
    const Real ooyS = 1/yS, ooyS2 = ooyS*ooyS, ooyS3=ooyS*ooyS2;

    //3. Multiply by Ainv to form coefficient vector a
    fV.resize(16);                       
    aijV.resize(16);
            
    //Compute Bicubic coefficients only if we're in a new patch
    //else use the old ones, because this is an expensive step!
    if( !(_pXYIdx[0] == x0 && _pXYIdx[2] == y0) ){
        fV(0) =  _f(x0,y0);
        fV(1) =  _f(x1,y0);
        fV(2) =  _f(x0,y1);
        fV(3) =  _f(x1,y1);

        fV(4) = _fx(x0,y0)*xS;
        fV(5) = _fx(x1,y0)*xS;
        fV(6) = _fx(x0,y1)*xS;
        fV(7) = _fx(x1,y1)*xS;
    
        fV(8)  = _fy(x0,y0)*yS;
        fV(9)  = _fy(x1,y0)*yS;
        fV(10) = _fy(x0,y1)*yS;
        fV(11) = _fy(x1,y1)*yS;

        fV(12)  = _fxy(x0,y0)*xS*yS;
        fV(13)  = _fxy(x1,y0)*xS*yS;
        fV(14)  = _fxy(x0,y1)*xS*yS;
        fV(15)  = _fxy(x1,y1)*xS*yS;

        getCoefficients(fV,aijV);
    }else{
        fV      = _fV;
        aijV    = _aV;
    }

    const Vec<16>& a = Vec<16>::getAs(&aijV[0]);


    //Compute where in the patch we are. This has to
    //be done everytime the location within the patch changes
    // ... which has happened if this code gets executed

    //--------------------------------------------------------------------------
    // Evaluate function value f (38 flops).

    // 8 flops
    const Real xpt = (aXY(0)-_x(x0))*ooxS, xpt2=xpt*xpt, xpt3=xpt*xpt2;
    const Real ypt = (aXY(1)-_y(y0))*ooyS, ypt2=ypt*ypt, ypt3=ypt*ypt2;
    // 12 flops
    const Mat44 mx(a[ 0],   a[ 1]*xpt,   a[ 2]*xpt2,   a[ 3]*xpt3,
                   a[ 4],   a[ 5]*xpt,   a[ 6]*xpt2,   a[ 7]*xpt3,
                   a[ 8],   a[ 9]*xpt,   a[10]*xpt2,   a[11]*xpt3,
                   a[12],   a[13]*xpt,   a[14]*xpt2,   a[15]*xpt3);
    // 12 flops
    const Vec4 xsum = mx.rowSum();
    // 6 flops
    const Real f = xsum[0] + ypt*xsum[1] + ypt2*xsum[2] + ypt3*xsum[3];

    //--------------------------------------------------------------------------
    // Evaluate first derivatives fx, fy (43 flops).

    // fy is 9 flops
    const Real dypt = ooyS, dypt2= 2*ypt*ooyS, dypt3= 3*ypt2*ooyS;
    const Real fy = dypt*xsum[1] + dypt2*xsum[2] + dypt3*xsum[3];

    // fx is 34 flops
    const Real dxpt=ooxS, dxpt2=2*xpt*ooxS, dxpt3=3*xpt2*ooxS;
    const Mat43 mdx(a[ 1]*dxpt,    a[ 2]*dxpt2,    a[ 3]*dxpt3,
                    a[ 5]*dxpt,    a[ 6]*dxpt2,    a[ 7]*dxpt3,
                    a[ 9]*dxpt,    a[10]*dxpt2,    a[11]*dxpt3,
                    a[13]*dxpt,    a[14]*dxpt2,    a[15]*dxpt3);
    const Vec4 dxsum = mdx.rowSum();
    const Real fx   = dxsum[0] + ypt*dxsum[1] + ypt2*dxsum[2] + ypt3*dxsum[3];

    //--------------------------------------------------------------------------
    // Evaluate second derivatives fxy, fxx, fyy (40 flops).

    // fxy, fyy are 11 flops
    const Real fxy = dypt*dxsum[1] + dypt2*dxsum[2] + dypt3*dxsum[3];
    const Real dyypt2=2*ooyS2, dyypt3=6*ypt*ooyS2;
    const Real fyy = dyypt2*xsum[2] + dyypt3*xsum[3];

    // fxx is 29 flops
    const Real dxxpt2=2*ooxS2, dxxpt3=6*xpt*ooxS2;
    const Mat42 mdxx(a[ 2]*dxxpt2,    a[ 3]*dxxpt3,
                     a[ 6]*dxxpt2,    a[ 7]*dxxpt3,
                     a[10]*dxxpt2,    a[11]*dxxpt3,
                     a[14]*dxxpt2,    a[15]*dxxpt3);
    const Vec4 dxxsum = mdxx.rowSum();
    const Real fxx  = dxxsum[0] + ypt*dxxsum[1] + ypt2*dxxsum[2] + ypt3*dxxsum[3];

    //--------------------------------------------------------------------------
    // Evaluate third derivatives fxxx, fxxy, fyyy, fxyy (21 flops).

    // 10 flops
    const Real dyyypt3=6*ooyS3;
    const Real fyyy = dyyypt3*xsum[3];
    const Real fxyy = dyypt2*dxsum[2] + dyypt3*dxsum[3];
    const Real fxxy = dypt*dxxsum[1] + dypt2*dxxsum[2] + dypt3*dxxsum[3];

    // 11 flops
    const Real dxxxpt3=6*ooxS3;
    const Vec4 mdxxx(a[ 3]*dxxxpt3,
                     a[ 7]*dxxxpt3,
                     a[11]*dxxxpt3,
                     a[15]*dxxxpt3);
    const Real fxxx = mdxxx[0] + ypt* mdxxx[1] + ypt2*mdxxx[2] + ypt3*mdxxx[3];

    //Populate the output vector
    aFdF(0) = f;
    aFdF(1) = fx;
    aFdF(2) = fy;
    aFdF(3) = fxy;
    aFdF(4) = fxx;
    aFdF(5) = fyy;
    aFdF(6) = fxxx;
    aFdF(7) = fxxy;
    aFdF(8) = fyyy;
    aFdF(9) = fxyy;

    if(_debug == true){
        cout<<" getFdF" << endl;
        cout <<"Member variables" <<endl;
        cout << "_x" << _x << endl;
        cout << "\n"<<endl;
        cout << "_y" << _y << endl;
        cout << "\n"<<endl;
        cout << "_f" << _f << endl;
        cout << "\n"<<endl;
        cout << "_fx" << _fx << endl;
        cout << "\n"<<endl;
        cout << "_fy" << _fy << endl;
        cout << "\n"<<endl;
        cout << "_fxy" << _fxy << endl;
        cout << "\n\n\n"<<endl;

        cout <<" Intermediate variables " << endl;
        cout <<"XY: " << aXY << endl;
        printf("[x0 x1], [y0 y1]: [%d %d],[%d %d]\n",x0,x1,y0,y1);
        printf("[x0V x1V], [y0V y1V]: [%f %f],[%f %f]\n",_x(x0),_x(x1),_y(y0),_y(y1));
        cout <<" xS " << xS << " yS " << yS << endl;
        printf("(xp,yp): (%f %f)\n",xpt,ypt);
        //cout << " axpt : " << axVec << endl;
        //cout << " adxpt : " << adxVec << endl;
        //cout << " ypt : " << yVec << endl;
        //cout << " dypt : " << dyVec << endl;
        cout << "\n\n\n"<<endl;

        cout <<" Final Output Vector " << endl;
        cout << "f,fx,fy,fxy" << aFdF << endl;
    }

    //Update the previous value records
            
    _pXYVal(0)    = aXY(0);
    _pXYVal(1)    = aXY(1);
    _pXYIdx[0]  = x0;
    _pXYIdx[1]  = x1;
    _pXYIdx[2]  = y0;
    _pXYIdx[3]  = y1;
            
    _pFdF       = aFdF;
    _fV         = fV;
    _aV         = aijV;
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
int BicubicSurface::Guts::
calcLowerBoundIndex(const Vector& aVec, Real aVal, int pIdx, bool evenlySpaced) const
{
    int idxLB = -1;
    bool idxComputed = false;
    bool scalarVector = false;

    if(aVec.size()<=1) scalarVector = true;

    
    if(scalarVector == false){
        //Compute index if the data is evenly spaced    
        if(evenlySpaced==true){
            Real sp = aVec(1)-aVec(0);
            Real minVal = aVec(0);
            idxLB = (int)floor((aVal-minVal)/sp);
            idxComputed = true;
        }else{ 
            //We have to search for the index
        
            //1. Do a local search around the previous index, if one is given
            if(pIdx >= 0 && pIdx <= aVec.size()){                
                int idxL = std::max(pIdx-1,0);
                int idxU = std::min(pIdx+1,aVec.size());
                if(aVal < aVec(idxU)  && aVal > aVec(idxL)){
                    if(aVal > aVec(pIdx))
                        idxLB = pIdx;    
                    else
                        idxLB = idxL;                        
                    idxComputed = true;
                }    
            }
    
            
            if(idxComputed == false){
                //2. If still not found check the end points
                if(aVal <= aVec(0)){ 
                    idxLB = 0;
                    idxComputed = true;
                }else if(aVal >= aVec(aVec.size()-1)){
                    idxLB = aVec.size()-1;
                    idxComputed = true;
                }
        
                //If still not found, use bisection method to find the appropriate index
                if(idxComputed == false){
                    int idxL = 0;
                    int idxU = aVec.size()-1;
                    int idxM = 0;
                    while (idxU-idxL > 1) {
                        idxM = (idxU+idxL)/2;
                        if (aVec(idxM) > aVal)
                            idxU = idxM;
                        else
                            idxL = idxM;
                    }
                    if(aVal > aVec(idxL) && aVal < aVec(idxU)){
                        idxLB = idxL;
                        idxComputed = true;
                    }


                }
            }
        }
    }
    //Check to ensure that idxLB is within the bounds of the vector. 
    //This index could be outside the bounds if the data is evenly spaced, and aVal
    //is outside the range of the vector.
    if(idxLB < 0)
        idxLB = 0;
    if(idxLB >= (aVec.size()-1))
        idxLB = aVec.size()-2;

    if (!evenlySpaced) {
        const Real* lower = std::lower_bound(&aVec[0], &aVec[0] + aVec.size(), aVal);
        int lowerIx = std::max((int)(lower-&aVec[1]), 0);

        //printf("-----> Matt computed %d, --lower=%d\n",
        //    idxLB, lowerIx);

        if (idxLB != lowerIx) {
            printf("*** DISAGREED WITH LOWER FOR aVal=%g Matt=%g --lower=%g\n", 
            aVal, aVal-aVec[idxLB], aVal-aVec[lowerIx]);
            assert(false);
        }
    }

    return idxLB;
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
getCoefficients(const Vector& fV, Vector& aV) const {
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
    Matrix AinvM(16,16,Ainv);
    aV = AinvM*fV; // So cool that I can do this in C++! Go Sherm!
    */

    // Matt's masterful Maple work, plus a little manual hacking by Sherm
    // produced this version:
    aV.resize(16);

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




