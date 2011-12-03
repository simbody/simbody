#ifndef SimTK_BICUBIC_SURFACE_H_
#define SimTK_BICUBIC_SURFACE_H_

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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"

#include <limits>

namespace SimTK { 

/** This function will create a smooth surface that approximates a two-argument
function F(X,Y) from a given set of samples of that function on a rectangular
grid. A bicubic surface interpolation is used to approximate the function
between the sample points. That is desirable for simulation use because it is 
continuous up to the second derivative, providing smoothly varying first 
derivatives, and a very smooth surface. The third derivatives will be 
discontinuous between grid boundaries; all higher derivatives are zero.

The user only need provide two vectors x and y defining the sample points,
and a matrix f that defines the value of the function at each sample (you can
think of that as the height Z of the surface over the X-Y plane). If the 
samples along both axes are regularly spaced, x and y can be defined just by
giving the spacing; otherwise, the sample locations are given explicitly.

Graphically if these vectors and matrices were laid next to each other 
consistently with how the surface is computed the diagram would look like this:
<pre>
             y(0)       y(1)    ...   y(ny-1)
            ------     ------         --------
    x(0)  |  f(0,0)     f(0,1)  ...   f(0,ny-1)
    x(1)  |  f(1,0)     f(1,1)  ...   f(1,ny-1)
     .    |    .          .              .
     .    |    .          .              .
     .    |    .          .              .
  x(nx-1) | f(nx-1,0)  f(nx-1,1)    f(nx-1,ny-1)
</pre>
such that f(i,j)=F(x(i),y(j)).

Note that the each XY location can only have a unique value associated with it
-- cave-like structures cannot be represented using this interpolation method.
    
Technically a bicubic surface interpolation requires the partial derivatives 
fx, fy and fxy at each of the grid points. To take this burden from the user, 
these partial derivative matricies are computed using only the the supplied 
points for X, Y and F. For the interested reader, these partial derivatives are 
computed by fitting splines through the points provided, and then taking 
derivatives of splines. 

These splines will pass through the points exactly when the smoothness 
parameter of the surface is set to 0, and will be interpolated using natural 
cubic splines. When the smoothness paramter is between 0 and 1, the surface 
will be 'relaxed' using the algorithm used in SplineFitter, and will not
exactly pass through the points given, but will smoothly come close to the 
points. The smoothness parameter can thus
be used to generate a surface that smoothly interpolates noisy surface data.

Here is the Wikipedia entry from which we implemented this method:
http://en.wikipedia.org/wiki/Bicubic_interpolation

@see SplineFitter for implementation notes regarding smoothing. **/
class SimTK_SIMMATH_EXPORT BicubicSurface : public Function_<Real> {
public:
    //--------------------------------------------------------------------------
    // CONSTRUCTION
    //--------------------------------------------------------------------------
    
    /** Construct an uninitialized BicubicSurface object. This can be filled
    in later by assignment. **/    
    BicubicSurface();

    // default destructor, copy constructor, copy assignment

    
    /** Construct a bicubic surface that approaches the grid points in f(x,y)
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
        points stored in matrix \a f. (Optional, default is 0.)

    If your sample points are regularly spaced, use the other constructor. **/
    BicubicSurface(const Vector& x, const Vector& y, const Matrix& f, 
                   Real smoothness=0);

    /** Construct a bicubic surface that approaches the grid points in f(x,y)
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
        points stored in matrix \a f. (Optional, default is 0.)

    If your sample points are not regularly spaced, use the other constructor.
    **/
    BicubicSurface(Real xSpacing, Real ySpacing, Matrix& f, Real smoothness=0);

    /** Calculate the value of the surface at a particular XY coordinate. Note
    that XY must be a vector with only 2 elements in it (because this is a
    2-argument function), anything else will throw an exception. This is the
    required implementation of the Function base class pure virtual.
     
    @param XY the 2-Vector of input arguments X and Y. 
    @return The interpolated value of the function at point (X,Y). **/
    virtual Real calcValue(const Vector& XY) const;
    
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
    virtual Real calcDerivative(const Array_<int>& derivComponents, 
                                const Vector& XY) const;
    
    /** The surface interpolation only works within the grid defined by the 
    vectors x and y used in the constructor. This function check to see if an 
    XYval is within the defined bounds of this particular BicubicSurface.
     
    @param XY   The vector of exactly 2 input arguments that define the XY 
                location on the surface.
    @return \c true if the point is in range, \c false otherwise. 
    
    An attempt to invoke calcValue() or calcDerivative() on an out-of-range
    point will raise an exception; use this method to check first if you 
    are not sure. **/
    bool isSurfaceDefined(const Vector& XY) const;

    /** This implements the Function base class pure virtual; here it
    always returns 2 (X and Y). **/
    virtual int getArgumentSize() const {return 2;}

    /** This implements the Function base class pure virtual specifying how
    many derivatives can be taken of this function; here it is unlimited.
    However, note that a bicubic surface is continuous up to the second 
    derivative, discontinuous at the third, and zero for any derivatives equal 
    to or higher than the fourth. **/
    virtual int getMaxDerivativeOrder() const 
    {   return std::numeric_limits<int>::max(); }


    /**
     
     @return The total number of elements in the matrix f 
    */
    int getFnelt() const;
    /**
    
    @return The number of elements in the vector X
    */
    int getXnelt() const;
    /**
    
    @return  The number of elements in the vector Y
    */
    int getYnelt() const;
    
    /** @name           Methods for debugging and testing
    Don't call these methods unless you are sure you know what you're doing.
    There is no guarantee that these will remain in the API from release to
    release. **/
    /**@{**/

    /**
    DEBUGGING CODE ONLY. DO NOT USE.    
    Sets a flag that prints a lot of useful debugging data to the screen when
    this class is used

    @param aDebug   setting this value to true will cause a lot of data to be
                    printed to the screen
    */
    void setDebug(bool aDebug);
    
    
    /**
    DEBUGGING CODE ONLY. DO NOT USE.    
    
    @return  the matrix f, which defines the grid points that the surface
    actually passes through
    */
    Matrix getf();

    /**
    DEBUGGING CODE ONLY. DO NOT USE.   
    
    @return  the matrix fx, which defines the values of fx(x,y) at the grid pts
    */
    Matrix getfx();
    
    /**
    DEBUGGING CODE ONLY. DO NOT USE.   
    
    @return  the matrix fy, which defines the values of fy(x,y) at the grid pts
    */
    Matrix getfy();

    /**
    DEBUGGING CODE ONLY. DO NOT USE.   
    
    @return  the matrix fxy, which defines the values of fxy(x,y) at the grid pts
    */
    Matrix getfxy();

    /**
    DEBUGGING CODE ONLY. DO NOT USE.   
    
    @return  the vector x, which defines one side of the mesh grid for which
             f(x,y) is defined.
    */
    Vector getx();
    
    /**
    DEBUGGING CODE ONLY. DO NOT USE.   
    
    @return  the vector y, which defines one side of the mesh grid for which
             f(x,y) is defined.
    */
    Vector gety();

    /**
    DEBUGGING CODE ONLY. DO NOT USE.  

    @param  XY an (X,Y) location within the grid
    @return  The 16x1 vector that defines the values of f,fx,fy, and fxy at the 
             corners. The vector is given in the following order

             f(0,0)   f(1,0)   f(0,1)   f(1,1),
             fx(0,0)  fx(1,0)  fx(0,1)  fx(1,1),
             fy(0,0)  fy(1,0)  fy(0,1)  fy(1,1),
             fxy(0,0) fxy(1,0) fxy(0,1) fxy(1,1)
    */
    Vector getPatchFunctionVector(Vector XY);
    
    /**
    DEBUGGING CODE ONLY. DO NOT USE. 

    @param  XY an (X,Y) location within the grid
    @return  The 16x1 vector that defines the values of the coefficients aij of 
             the bicubic surface interpolation in this patch. The coefficient 
             vector is given in the following order
             
             a00,a10,a20,a30,a01,a11,a21,a31,a02,a12,a22,a32,a03,a13,a23,a33
    */
    Vector getPatchBicubicCoefficients(Vector XY);
    /**
    DEBUGGING CODE ONLY. DO NOT USE.  

    A constructor for a bicubic surface that sets the partial derivatives of the
    surface to the values specified by fx, fy, and fxy.

    @param x: vector of X grid points (minimum 4 values)
    @param y: vector of Y grid points (minimum 4 values)
    @param f:   matrix of the surface heights evaluated at the grid formed 
                by x and y (minumum 4x4)
    @param fx:  matrix of the partial derivative of f w.r.t to x (minumum 4x4)
    @param fy:  matrix of the partial derivative of f w.r.t to y (minumum 4x4)
    @param fxy: matrix of the partial derivative of f w.r.t to x,y (minumum 4x4)
    */
    BicubicSurface(Vector& x, Vector& y, Matrix& f, 
               Matrix& fx, Matrix& fy, Matrix& fxy);
    /**@}**/

private:
    void setNull();
    void setEqual(const BicubicSurface &aSpline);
    int calcLowerBoundIndex(const Vector& vecV, Real value, int pIdx,
                                                    bool evenlySpaced) const;
    void getCoefficients(const Vector& f, Vector& aV) const;
    void getFdF(const Vector& aXY, 
                Vector& fV, Vector& aijV, Vector& aFdF) const;

//=============================================================================
// MEMBER VARIABLES
//=============================================================================
private:
    // PROPERTIES
    /* Array of values for the independent variables (i.e., the spline knot
    sequence).  Each array must be monotonically increasing, of size nx and 
    ny elements */
    Vector _x;
    Vector _y;

    /*  2D nx X ny z values that correspond to the values at the grid defined
        by x and y. */  
    Matrix _f;
    bool _flagXEvenlySpaced;
    bool _flagYEvenlySpaced;

    //This is all 'hint' data that is stored between calls to this
    //object to reduce the time required to find the patch the user
    //is in.
    mutable Vector      _pXYVal;
    mutable Array_<int> _pXYIdx;
    mutable Vector      _pFdF;
    mutable Vector      _fV;
    mutable Vector      _aV;

    //The matrices of partial differentials of the surface w.r.t. x, y and xy
    Matrix _fx;
    Matrix _fy;
    Matrix _fxy;

    //A private debugging flag - if set to true, a lot of useful debugging
    //data will be printed tot the screen
    bool _debug;


//=============================================================================
};    // END class BicubicSurface


}; //namespace
//=============================================================================
//=============================================================================

#endif  // SimTK_BICUBIC_SURFACE_H_
