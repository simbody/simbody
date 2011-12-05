#ifndef SimTK_BICUBIC_SURFACE_GUTS_H_
#define SimTK_BICUBIC_SURFACE_GUTS_H_

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

/** @file
This file defines the internal BicubicSurface::Guts class, which is the private
implementation class of BicubicSurface, and BicubicSurface::PatchHint::Guts
which is used to improve performance for repeated access to the same patch or
nearby patches. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
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
    Guts() : level(-1) {}
    // Default copy constructor, copy assignment, destructor

    bool isEmpty() const {return x0 < 0;}
    void clear() {x0=y0=level = -1;}


    // Patch information. 
    // This is valid whenever level >= 0 and does not change
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

    // -1: nothing, 0: function value only, 1: function value and 1st derivs,
    // 2: function value, 1st & 2nd derivatives, 3: plus 3rd derivatives.
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
/** Private implementation of BicubicSurface. **/
class SimTK_SIMMATH_EXPORT BicubicSurface::Guts {
friend class BicubicSurface;
public:
    /** Construct an uninitialized BicubicSurface object. This can be filled
    in later by assignment. Reference count is zero after construction. **/    
    Guts() {construct();}
    ~Guts() {assert(referenceCount==0);}

    
    /** Construct a bicubic surface that approximates f(x,y)
    with the spacing between each grid point in f defined by the vectors x and 
    y. The smoothness paramter controls how closely the surface approaches the 
    grid points specified in matrix f, with the default being that the surface
    will pass exactly through those points.

    @param x    Vector of sample locations along the X axis (minimum 4 values).
                Must be monotonically increasing (no duplicates).
    @param y    Vector of sample locations along the Y axis (minimum 4 values).
                Must be monotonically increasing (no duplicates).            
    @param f    
        Matrix of function values (or surface heights) evaluated at the grid 
        points formed by x and y (dimension x.size() X y.size()), such that 
        f(i,j) is F(x[i],y[j]) where F is the function being approximated here.
    @param smoothness 
        A value of 0 will force surface to pass through all of the 
        points in f(x,y). As smoothness tends to 1, the surface will 
        become smoother and smoother, but will not pass through the knot 
        points stored in matrix \a f.

    If your sample points are regularly spaced, use the other constructor. **/
    Guts(const Vector& x, const Vector& y, const Matrix& f, 
                   Real smoothness);

    /** Construct a bicubic surface that approximates f(x,y)
    over a grid with regular spacing in both the x and y directions. The 
    smoothness parameter controls how closely the surface approaches the 
    grid points specified in matrix f, with the default being that the surface
    will pass exactly through those points.

    @param xSpacing     Spacing for x-axis sampling; must be greater than 0.
    @param ySpacing     Spacing for y-axis sampling; must be greater than 0.
    @param f            
        Matrix of function values (or surface heights) evaluated at points of 
        the x-y plane regularly sampled using the upplied spacings. Can be 
        rectangular but must have minimum dimension 4x4. Here 
        f(i,j)=F(i*xSpacing,j*ySpacing) where F is the function being 
        approximated.
    @param smoothness 
        A value of 0 will force surface to pass through all of the 
        points in f(x,y). As smoothness tends to 1, the surface will 
        become smoother and smoother, but will not pass through the knot 
        points stored in matrix \a f.

    If your sample points are not regularly spaced, use the other constructor.
    **/
    Guts(Real xSpacing, Real ySpacing, const Matrix& f, Real smoothness)
    {   SimTK_THROW2(Exception::UnimplementedVirtualMethod, 
        "BicubicSurface::Guts", "ctor(xSpacing,ySpacing,f,smoothness)"); }

    /** Calculate the value of the surface at a particular XY coordinate. Note
    that XY must be a vector with only 2 elements in it (because this is a
    2-argument function), anything else will throw an exception. This is the
    required implementation of the Function base class pure virtual.
     
    @param XY the 2-Vector of input arguments X and Y. 
    @return The interpolated value of the function at point (X,Y). **/
    Real calcValue(const Vec2& XY, PatchHint& hint) const;
    
    /** Calculate a partial derivative of this function at a particular point.  
    Which derivative to take is specified by listing the input components
    (0==x, 1==y) with which to take it. For example, if derivComponents=={0}, 
    that indicates a first derivative with respective to argument x.  
    If derivComponents=={0, 0, 0}, that indicates a third derivative with
    respective to argument x.  If derivComponents=={0, 1}, that indicates 
    a partial second derivative with respect to x and y, that is Df(x,y)/DxDy.
    (We use capital D to indicate partial derivative.)
     
    @param derivComponents  
        The input components with respect to which the derivative should be 
        taken. Its size must be less than or equal to the  value returned by 
        getMaxDerivativeOrder().      
    @param XY    
        The vector of two input arguments that define the XY location on the 
        surface. 
    @return The interpolated value of the selected function partial derivative
    for arguments (X,Y). **/
    Real calcDerivative(const Array_<int>& derivComponents, 
                        const Vec2& XY, PatchHint& hint) const;
    
    /** The surface interpolation only works within the grid defined by the 
    vectors x and y used in the constructor. This function check to see if an 
    XYval is within the defined bounds of this particular BicubicSurface.
     
    @param XY   The vector of exactly 2 input arguments that define the XY 
                location on the surface.
    @return \c true if the point is in range, \c false otherwise. 
    
    An attempt to invoke calcValue() or calcDerivative() on an out-of-range
    point will raise an exception; use this method to check first if you 
    are not sure. **/
    bool isSurfaceDefined(const Vec2& XY) const;

    int getReferenceCount() const {return referenceCount;}
    void incrReferenceCount() const {++referenceCount;}
    int decrReferenceCount() const 
    {assert(referenceCount);return --referenceCount;}

    /** @return The total number of elements in the matrix f **/
    int getFnelt() const {return _ff.nelt();}
    /** @return The number of elements in the vector X **/
    int getXnelt() const {return _x.nelt();}
    /** @return  The number of elements in the vector Y **/
    int getYnelt() const {return _y.nelt();}
    
    /** @name           Methods for debugging and testing
    Don't call these methods unless you are sure you know what you're doing.
    There is no guarantee that these will remain in the API from release to
    release. **/
    /**@{**/

    /** Sets a flag that prints a lot of useful debugging data to the screen
    when this class is used.
    @param aDebug   setting this value to true will cause a lot of data to be
                    printed to the screen **/
    void setDebug(bool aDebug) {_debug = aDebug;}
     
    /** @return  the matrix ff, which defines the grid points that the surface
    actually passes through, and derivatives fx,fy,fxy. **/
    const Matrix_<Vec4>& getff() const {return _ff;}

    ///** @return  the matrix fx, which defines the values of fx(x,y) at the 
    //grid pts **/
    //const Matrix& getfx() const {return _fx;}
    
    ///** @return  the matrix fy, which defines the values of fy(x,y) at the 
    //grid pts **/
    //const Matrix& getfy() const {return _fy;}

    ///** @return  the matrix fxy, which defines the values of fxy(x,y) at the 
    //grid pts **/
    //const Matrix& getfxy() const {return _fxy;}

    /** @return  the vector x, which defines one side of the mesh grid for 
    which f(x,y) is defined. **/
    const Vector& getx() const {return _x;}
    
    /** @return  the vector y, which defines one side of the mesh grid for 
    which f(x,y) is defined. **/
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
        Vec<16> fV, aV;
        Vec<10> aFdF;
        BicubicSurface::PatchHint hint;
        getFdF(XY,-1,fV,aV,aFdF,hint); // just need patch info
        return fV;

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
        Vec<16> fV, aV;
        Vec<10> aFdF;
        BicubicSurface::PatchHint hint;
        getFdF(XY,-1,fV,aV,aFdF,hint); // just need patch info
        return aV;
    }

    /** A constructor for a bicubic surface that sets the partial derivatives 
    of the surface to the values specified by fx, fy, and fxy.

    @param x: vector of X grid points (minimum 4 values)
    @param y: vector of Y grid points (minimum 4 values)
    @param f:   matrix of the surface heights evaluated at the grid formed 
                by x and y (minumum 4x4)
    @param fx:  matrix of the partial derivative of f w.r.t to x (minumum 4x4)
    @param fy:  matrix of the partial derivative of f w.r.t to y (minumum 4x4)
    @param fxy: matrix of the partial derivative of f w.r.t to x,y (minumum 4x4)
    */
    Guts(const Vector& x, const Vector& y, const Matrix& f, 
         const Matrix& fx, const Matrix& fy, const Matrix& fxy);
    /**@}**/

private:
    int calcLowerBoundIndex(const Vector& vecV, Real value, int pIdx,
                                                    bool evenlySpaced) const;
    void getCoefficients(const Vec<16>& f, Vec<16>& aV) const;
    void getFdF(const Vec2& aXY, int wantLevel,
                Vec<16>& fV, Vec<16>& aijV, Vec<10>& aFdF,
                BicubicSurface::PatchHint& hint) const;

    // This is called from each constructor to initialize this object.
    void construct() {
        referenceCount = 0;
        resetStatistics();
        _flagXEvenlySpaced = _flagYEvenlySpaced = false;
        _debug = false;
    }

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
    // ny elements.
    Vector _x, _y;
    bool _flagXEvenlySpaced, _flagYEvenlySpaced;

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
