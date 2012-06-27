/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2011-12 Stanford University and the Authors.        *
 * Authors: Matthew Millard                                                   *
 * Contributors: Michael Sherman                                              *
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

#include "SimTKmath.h"

// Include the private implementation class declaration for testing purposes;
// this is not part of the API.
#include "../src/BicubicSurface_Guts.h"


#include <cstdio>
#include <iostream>
#include <fstream>

using namespace SimTK;
using namespace std;

/**
This function computes a standard central difference dy/dx. 
If extrap_endpoints is set to 1, then the derivative at the end points 
is estimated by linearly extrapolating the dy/dx values beside the  end points

 @param x domain vector
 @param y range vector
 @param extrap_endpoints:(false) Endpoints of the returned vector will be zero, 
                                  because a central difference is undefined at 
                                  these endpoints
                          (true) Endpoints are computed by linearly 
                                 extrapolating using a first difference from 
                                 the neighboring 2 points

 @returns dy/dx computed using central differences
*/
Vector getCentralDifference(Vector x, Vector y, bool extrap_endpoints) {
    Vector dy(x.size());
    Real dx1,dx2;
    Real dy1,dy2;
    int size = x.size();
    for(int i=1; i<x.size()-1; i++){
        dx1 = x(i)-x(i-1);
        dx2 = x(i+1)-x(i);
        dy1 = y(i)-y(i-1);
        dy2 = y(i+1)-y(i);
        dy(i)= 0.5*dy1/dx1 + 0.5*dy2/dx2;
    }

    if(extrap_endpoints == true){
        dy1 = dy(2)-dy(1);
        dx1 = x(2)-x(1);
        dy(0) = dy(1) + (dy1/dx1)*(x(0)-x(1));

        dy2 = dy(size-2)-dy(size-3);
        dx2 = x(size-2)-x(size-3);
        dy(size-1) = dy(size-2) + (dy2/dx2)*(x(size-1)-x(size-2));
    }
    return dy;
}

/**
 Return the value, and set of first partial deriviatives of a 2D function f
 that defines a surface f(x,y). There are 4 functions of choice, which
 can be specified by the paramter fcnType

 @param x : x argument of the function f(x,y)
 @param y : y argument of the function f(x,y)
 @param fcnType [0,1,2,3,4]. Chooses one of the following functions for
                             f(x,y):
 
         fcnType = 0 :f(x,y) = 0;
         fcnType = 1 :f(x,y) = 2*x + y
         fcnType = 2 :f(x,y) = xy
         fcnType = 3 :f(x,y) = cos( (3x^2+y^2)^0.5 )
         fcnType = 4 : f(x,y) = 3x^2 + y^2
*/
Vector getAnalyticFunction(Real x, Real y, int fcnType){
    Vector fdF(4);
    fdF = -1;

    switch(fcnType){
        case 0:    //f(x,y) = 0;
            fdF = 0;
            break;
        case 1: //f(x,y) = 2*x + y
            fdF(0) = 2*x + y;    //f
            fdF(1) = 2;            //fx
            fdF(2) = 1;            //fy
            fdF(3) = 0;            //fxy
            break;
        case 2: //f(x,y) = xy
            fdF(0) = x*y;        //f
            fdF(1) = y;            //fx
            fdF(2) = x;            //fy
            fdF(3) = 1;            //fxy
            break;
        case 3:    //f(x,y) = cos( (3x^2+y^2)^0.5 );
            //f
            fdF(0) = cos( sqrt((3*x*x + y*y)) );        
            //fx - exported from Maple (didn't trust myself not to make a typing mistake
            fdF(1) = -0.3e1*x*sin( sqrt( (3*x*x + y*y)) ) * pow( (3*x*x + y*y + 1e-6), -0.1e1 / 0.2e1) ;            
            //fy
            fdF(2) = -sin(sqrt( (3 * x * x + y * y))) * pow((3 * x * x + y * y + 1e-6), -0.1e1 / 0.2e1) * y;            
            //fxy
            fdF(3) =  -0.3e1 * cos(sqrt((3 * x * x + y * y)))/(3 * x * x + y * y + 1e-6) *  y *  x + 0.3e1 * sin(sqrt( (3 * x * x + y * y))) * pow( (3 * x * x + y * y + 1e-6), -0.3e1 / 0.2e1) * x *  y;            

            break;
        case 4: //f(x,y) = 3x^2 + y^2
            fdF(0) = 3*x*x + y*y;
            fdF(1) = 6*x;
            fdF(2) = 2*y;
            fdF(3) = 0;
            break;
        default:
            cout << "Invalid fcnType in testBicubicSurface.cpp: getAnayticFunction";
    }


    return fdF;
}


/**
 This function will generate a rectangular grid that spans from xmin to xmax 
 in size number of steps, and also from ymin to ymax in size number of steps.
 Although the grid spacing can be different in the x and y dimensions, within 
 these dimensions the grids are equally spaced (by xDelta and yDelta).

 An analytic function (chosen using the fcnType variable) is used to generate
 f(x,y) values at each grid point. These values are used to initialize a 
 bicubic surface using the advanced test constructor that sets the partial
 derivatives fx, fy, and fxy directly.

 The values of the bicubic surface are evaluated at the grid points, we'll 
 call them knot points, and are asserted to be equal to the analytic function
 at these values. Additionally, every grid is evaluated at its center, and 
 the value of the bicubic surface is asserted to be equal to the analytic 
 function at the mid point to within a tolerance. This tolerance is a function
 of the grid size. This tolerance has been determined hueristically, so if you 
 try a new function and the test fails, look closely at the values to see if 
 its really failing or if the tolerance is just too tight.

 @params xmin: the minimum value of the x,y grid in the x dimension
 @params xmax: the maximum value of the x,y grid in the x dimension
 @params ymin: the minimum value of the y grid in the y dimension
 @params ymax: the maximum value of the y grid in the y dimension
 @params size: the number of steps to take to go from xmin to xmax, and
               ymin to ymax
 @params fcnType: An integer value [0-4] that picks an analytical function
                  to use for comparision purposes.
 @params flag_verbosePrint: false: print only the maximum error at the 
                                   knot points, mid points and the 
                                   tolerance used at the assertions

                            true: Additionally print the values of f, fx,
                                  fy, and fxy at the knots and the mid
                                  points if there are less than 10 steps 
@params flag_matlabcompre: true:    Will print
@returns nothing
*/
void testBicubicAgainstAnalyticFcn(Real xmin, Real xmax, Real ymin, 
                Real ymax, int size, int fcnType, bool flag_verbosePrint,
                                                   bool flag_matlabcompare){
        
    Real deltaX,deltaY;
    deltaX = (xmax-xmin)/(size-1);
    deltaY = (ymax-ymin)/(size-1);
        
    //Generate initialization data
    // two constant spaced vectors & height matrix & first derivatives to initialize the grid
    Vector x(size), y(size);
    Matrix z(size,size),zx(size,size),zy(size,size),zxy(size,size);

    //Generate test data to evaluate the error of the surface interpolation at the mid
    //point of each grid square. The `M' stands for mid-point
    Vector xM(size-1), yM(size-1);
    Matrix zM(size-1,size-1),zMx(size-1,size-1),zMy(size-1,size-1),zMxy(size-1,size-1);

    for (int i = 0; i < size; i++) {
        x(i) = xmin + ((Real)i)*deltaX;
        y(i) = ymin + ((Real)i)*deltaY;
        if(i<size-1){
            xM(i) = xmin + deltaX/(Real)2 + ((Real)i)*deltaX;
            yM(i) = ymin + deltaY/(Real)2 + ((Real)i)*deltaY;
        }
    }


    switch(fcnType){
        case 0:
            cout << "Testing bicubic surface against: f(x,y) = 0" <<endl;
            break;
        case 1:
            cout << "Testing bicubic surface against: f(x,y) = 2*x+y" <<endl;
            break;
        case 2:
            cout << "Testing bicubic surface against: f(x,y) = x*y" <<endl;
            break;
        case 3:
            cout << "Testing bicubic surface against: f(x,y) = cos( (3*x^2 + y^2)^0.5 ) " <<endl;
            break;
        case 4:
            cout << "Testing bicubic surface against: f(x,y) = 3*x^2 + y^2 " <<endl;
            break;
    }


    Vector fdF(4);
    Vector fdFM(4);
    fdF = 0;
    fdFM= 0;

    for(int i=0; i<size;i++){
        for(int j=0; j<size; j++){
            fdF = getAnalyticFunction(x(i),y(j),fcnType);
            //printf("i:%d, j:%d, x:%f, y:%f, f:%f, fx:%f, fy:%f, fxy:%f\n",i,j,x(i),y(j),fdF(0),fdF(1),fdF(2),fdF(3));

            z(i,j)         = fdF(0);
            zx(i,j)        = fdF(1);
            zy(i,j)        = fdF(2);
            zxy(i,j)       = fdF(3);

            if( i < size-1 && j < size-1){
                fdFM = getAnalyticFunction(xM(i),yM(j),fcnType);
                zM(i,j)         = fdFM(0);
                zMx(i,j)        = fdFM(1);
                zMy(i,j)        = fdFM(2);
                zMxy(i,j)       = fdFM(3);
            }
        }        
    }


    //Initialize the Bicubic Surface
    Real smoothness = 0.0;
    BicubicSurface bcs(x, y, z, zx, zy, zxy);
    const BicubicSurface::Guts& bcsg = bcs.getGuts();
    BicubicFunction bcsf(bcs);

    if(flag_verbosePrint == true && size <= 10){
        cout << "\n\nx:\n" << bcsg.getx() << endl;
        cout << "\n\ny:\n" << bcsg.gety() << endl;
        cout << "\n\nff:\n" << bcsg.getff() << endl;
    }

    //Test it at the knot points, mid grid and compute the error
    Vector errV(4); //Knot point error vector: f,fx,fy,fxy error
    Vector errVM(4);//Mid grid error vector:    f,fx,fy,fxy error

    Vector bcsV(4);    //Spline surface values at the knots
    Vector bcsMV(4);    //Spline surface values at the midpoints
        
    Vector XY(2); //XY value at the knot points;
    Vector XYM(2); //XY value at mid grid;

    const int ifxy[] = {1,0};
    const int ifxx[] = {0,0};
    const int ifyy[] = {1,1};
    const int ifxxx[] = {0,0,0};
    const int ifxxy[] = {0,0,1};
    const int ifyyy[] = {1,1,1};
    const int ifxyy[] = {1,1,0};

    Array_<int> fx(1); //Arguments required to get the correct derivative 
    Array_<int> fy(1); // from the calcDerivatie interface
    Array_<int> fxy(ifxy,ifxy+2);
    Array_<int> fxx(ifxx,ifxx+2);
    Array_<int> fyy(ifyy,ifyy+2);
    Array_<int> fxxx(ifxxx,ifxxx+3);
    Array_<int> fyyy(ifyyy,ifyyy+3);
    Array_<int> fxxy(ifxxy,ifxxy+3);
    Array_<int> fxyy(ifxyy,ifxyy+3);                

    fx[0]   =0;
    fy[0]   =1;

    errV = 0;
    errVM = 0;

    Matrix fk(size,size),fxk(size,size),fyk(size,size);
    Matrix fxyk(size,size),fxxk(size,size),fyyk(size,size);
    Matrix fxxyk(size,size),fxyyk(size,size);
    Matrix fxxxk(size,size),fyyyk(size,size);

    Matrix fMk(size-1,size-1),fxMk(size-1,size-1),fyMk(size-1,size-1);
    Matrix fxyMk(size-1,size-1),fxxMk(size-1,size-1),fyyMk(size-1,size-1);
    Matrix fxxyMk(size-1,size-1),fxyyMk(size-1,size-1);
    Matrix fxxxMk(size-1,size-1),fyyyMk(size-1,size-1);

    for(int i=0; i<size; i++){
        for(int j=0; j<size; j++){
            XY(0)=x(i);
            XY(1)=y(j);

            fk(i,j)     = bcsf.calcValue(XY);
            fxk(i,j)    = bcsf.calcDerivative(fx,XY);
            fyk(i,j)    = bcsf.calcDerivative(fy,XY);
            fxyk(i,j)   = bcsf.calcDerivative(fxy,XY);
            fxxk(i,j)   = bcsf.calcDerivative(fxx,XY);
            fyyk(i,j)   = bcsf.calcDerivative(fyy,XY);

            fxxxk(i,j)   = bcsf.calcDerivative(fxxx,XY);
            fyyyk(i,j)   = bcsf.calcDerivative(fyyy,XY);
            fxxyk(i,j)   = bcsf.calcDerivative(fxxy,XY);
            fxyyk(i,j)   = bcsf.calcDerivative(fxyy,XY);

            if( errV(0) < abs(fk(i,j) - z(i,j)) )
                errV(0) = abs(fk(i,j) - z(i,j));
            if( errV(1) < abs(fxk(i,j) - zx(i,j)) )
                errV(1) = abs(fxk(i,j) - zx(i,j));
            if( errV(2) < abs(fyk(i,j) - zy(i,j)) )
                errV(2) = abs(fyk(i,j) - zy(i,j));
            if( errV(3) < abs(fxyk(i,j) - zxy(i,j)) )
                errV(3) = abs(fxyk(i,j) - zxy(i,j));
                        
                /*if(abs(errV(0)) > 1e-4 ){
                    printf("Analytic (x,y),f,fx,fy,fxy: (%g,%g),%g, %g, %g, %g\n",
                        x(i),y(j),z(i,j),zx(i,j),zy(i,j),zxy(i,j));
                    printf("Approx.  (x,y),f,fx,fy,fxy: (%g,%g),%g, %g, %g, %g\n\n",
                        x(i),y(j),fk(i,j),fxk(i,j),fyk(i,j),fxyk(i,j));
                    bcs.setDebug(true);
                }*/

            if(i<size-1 && j<size-1){
                XYM(0)=xM(i);
                XYM(1)=yM(j);
                fMk(i,j)    = bcsf.calcValue(XYM);                            
                fxMk(i,j)   = bcsf.calcDerivative(fx,XYM);
                fyMk(i,j)   = bcsf.calcDerivative(fy,XYM);
                fxyMk(i,j)  = bcsf.calcDerivative(fxy,XYM);    
                fxxMk(i,j)  = bcsf.calcDerivative(fxx,XYM);
                fyyMk(i,j)  = bcsf.calcDerivative(fyy,XYM);

                fxxxMk(i,j)   = bcsf.calcDerivative(fxxx,XYM);
                fyyyMk(i,j)   = bcsf.calcDerivative(fyyy,XYM);
                fxxyMk(i,j)   = bcsf.calcDerivative(fxxy,XYM);
                fxyyMk(i,j)   = bcsf.calcDerivative(fxyy,XYM);

                if( errVM(0) < abs(fMk(i,j) - zM(i,j)) )
                    errVM(0) = abs(fMk(i,j) - zM(i,j));
                if( errVM(1) < abs(fxMk(i,j) - zMx(i,j)) )
                    errVM(1) = abs(fxMk(i,j) - zMx(i,j));
                if( errVM(2) < abs(fyMk(i,j) - zMy(i,j)) )
                    errVM(2) = abs(fyMk(i,j) - zMy(i,j));
                if( errVM(3) < abs(fxyMk(i,j) - zMxy(i,j)) )
                    errVM(3) = abs(fxyMk(i,j) - zMxy(i,j));
                        
            }
        }
    }


    if(flag_verbosePrint == true && size <= 10){

        cout << "\n\n Err f (@knot, calc):\n" << fk-z << endl;
        cout << "\n\n Err fx (@knot, calc):\n" << fxk -zx << endl;
        cout << "\n\n Err fy (@knot, calc):\n" << fyk -zy << endl;
        cout << "\n\n Err fxy (@knot, calc):\n" << fxyk -zxy << endl;
                    
        if(flag_matlabcompare == true){
            cout << "\n\n    x (@knot):\n" << x << endl;
            cout << "\n\n    x (@knot):\n" << y << endl;
            cout << "\n\n    f (@knot, calc):\n" << fk << endl;
            cout << "\n\n    fx (@knot, calc):\n" << fxk << endl;
            cout << "\n\n    fy (@knot, calc):\n" << fyk << endl;
            cout << "\n\n    fxy (@knot, calc):\n" << fxyk << endl;
            cout << "\n\n    fxx (@knot, calc):\n" << fxxk << endl;
            cout << "\n\n    fyy (@knot, calc):\n" << fyyk << endl;
            cout << "\n\n    fxxx (@knot, calc):\n" << fxxxk << endl;
            cout << "\n\n    fyyy (@knot, calc):\n" << fyyyk << endl;
            cout << "\n\n    fxxy (@knot, calc):\n" << fxxyk << endl;
            cout << "\n\n    fxyy (@knot, calc):\n" << fxyyk << endl;
                    

            cout << "\n\n    x (@mid):\n" << xM << endl;
            cout << "\n\n    y (@mid):\n" << yM << endl;
            cout << "\n\n    f (@mid, calc):\n" << fMk << endl;
            cout << "\n\n    fx (@mid, calc):\n" << fxMk << endl;
            cout << "\n\n    fy (@mid, calc):\n" << fyMk << endl;
            cout << "\n\n    fxy (@mid, calc):\n" << fxyMk << endl;
            cout << "\n\n    fxx (@mid, calc):\n" << fxxMk << endl;
            cout << "\n\n    fyy (@mid, calc):\n" << fyyMk << endl;
            cout << "\n\n    fxxx (@mid, calc):\n" << fxxxMk << endl;
            cout << "\n\n    fyyy (@mid, calc):\n" << fyyyMk << endl;
            cout << "\n\n    fxxy (@mid, calc):\n" << fxxyMk << endl;
            cout << "\n\n    fxyy (@mid, calc):\n" << fxyyMk << endl;
        }

    }

                
    Real mid_tol = (1e-1)*(deltaX/2+deltaY/2);
    Real knot_tol = 0;

    if(flag_verbosePrint == true){
        printf("    Smoothness set to : %f\n", smoothness);
        printf("    f:  err@knots %f, err@mid %f\n",errV(0),errVM(0));
        printf("    fx: err@knots %f, err@mid %f\n",errV(1),errVM(1));
        printf("    fy: err@knots %f, err@mid %f\n",errV(2),errVM(2));
        printf("    fxy:err@knots %f, err@mid %f\n\n",errV(3),errVM(3));                       
        printf("    Test tolerance for f(x,y) @knots : %f, @mid: %f\n\n", 
                                                        knot_tol,mid_tol);
        cout<< "    First derivatives are not tested because these derivatives" << endl;
        cout<< "    shouldn't match: the bicubic interpolation estimates" << endl;
        cout<< "    these derivatives using the derivative of a natural" << endl;
        cout<< "    cubic spline" << endl; 
    }

    //See if the values for f, fx, fy and fxy match the knot points
    //exactly
    SimTK_TEST_EQ(fk,z);
    SimTK_TEST_EQ(fxk,zx);
    SimTK_TEST_EQ(fyk,zy);
    SimTK_TEST_EQ_TOL(fxyk,zxy,1e-10);
    //See if the maximum error at the mid points are acceptable
    SimTK_TEST_EQ_TOL(errVM(0),0,mid_tol);
    
}

/**
This function will construct a single bicubic surface patch that goes from xmin,ymin
to xmax, ymax. A series of points within this patch will be computed using the bicubic
interpolation method, and the coefficients will be checked to ensure that the 
relationship between the 16 coefficients, aV, and the 16 corner conditions, fV, are
related to eachother through the endpoint conditions that define a bicubic surface 
interpolation (http://en.wikipedia.org/wiki/Bicubic_interpolation)

fV = A*aV

aV: [a00,   a10     a20     a30,    
     a01    a11     a21     a31, 
     a02    a12     a22     a32, 
     a03    a13     a23     a33]^T

fV:[f(0,0)   f(1,0)   f(0,1)   f(1,1)  
   fx(0,0)  fx(1,0)  fx(0,1)  fx(1,1) 
   fy(0,0)  fy(1,0)  fy(0,1)  fy(1,1)
  fxy(0,0) fxy(1,0) fxy(0,1) fxy(1,1)]

A is a 16x16 matrix that defines the relationship between the polynomial that enforces
the conditions that the polynomial has the same values and partial derivatives as the
function at the corners. To see this matrix in detail refer to the wikipedia page,
or to the code below. Note that A^(-1) is the one that is shown in the wikipedia page,
where as the one in the test code is a hand derived version of A.

 @params xmin: the minimum value of the x,y grid in the x dimension
 @params xmax: the maximum value of the x,y grid in the x dimension
 @params ymin: the minimum value of the y grid in the y dimension
 @params ymax: the maximum value of the y grid in the y dimension
 @params fcnType: An integer value [0-4] that picks an analytical function
                  to use for comparision purposes.
 @params smoothness: A value of 0 will make sure the patch goes through the 
                     desired points exactly. A value between 0 and 1 will
                     relax the surface.
 @returns nothing

*/
void testBicubicCoefficients(Real xmin,Real xmax,Real ymin, Real ymax, 
                                              int fcnType, Real smoothness){
    int size = 4;

    const Real A[] = {1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
                        1, 1, 1, 1,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
                        1, 0, 0, 0,  1, 0, 0, 0,  1, 0, 0, 0,  1, 0, 0, 0,
                        1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,  1, 1, 1, 1,
                        0, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
                        0, 1, 2, 3,  0, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
                        0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,  0, 1, 0, 0,
                        0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,  0, 1, 2, 3,
                        0, 0, 0, 0,  1, 0, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
                        0, 0, 0, 0,  1, 1, 1, 1,  0, 0, 0, 0,  0, 0, 0, 0,
                        0, 0, 0, 0,  1, 0, 0, 0,  2, 0, 0, 0,  3, 0, 0, 0,
                        0, 0, 0, 0,  1, 1, 1, 1,  2, 2, 2, 2,  3, 3, 3, 3,                  
                        0, 0, 0, 0,  0, 1, 0, 0,  0, 0, 0, 0,  0, 0, 0, 0,
                        0, 0, 0, 0,  0, 1, 2, 3,  0, 0, 0, 0,  0, 0, 0, 0,
                        0, 0, 0, 0,  0, 1, 0, 0,  0, 2, 0, 0,  0, 3, 0, 0,
                        0, 0, 0, 0,  0, 1, 2, 3,  0, 2, 4, 6,  0, 3, 6, 9};

    /*Ok we need at least a 4x4 grid to use the default bicubic surface 
    interpolation because the constructor forms the partial derivatives 
    using natural cubic splines. Natural cubic splines require at least 
    4 knot points to be defined.
    */
    Vector xV(size), yV(size), xeV(2*size-1), yeV(2*size-1);
    Matrix zM(size,size);
    Vector tmpV(size);

    Vec<16> fT, aV, fV, fVerr;
    Mat<16,16> AM(A), ATest;

    //Initialize the grid
    for(int i=0; i<size; i++){
        xV(i) = xmin + i*(xmax-xmin)/((Real)size-1.0);
        yV(i) = xmin + i*(ymax-ymin)/((Real)size-1.0);        
    }
    for(int i=0; i<size;i++){
        for(int j=0; j<size; j++){
            tmpV = getAnalyticFunction(xV(i),yV(j),fcnType);
            zM(i,j) = tmpV(0);
        }
    }

    //Create the bicubic surface using the regular constuctor
    BicubicSurface bcs(xV, yV, zM, smoothness);
    const BicubicSurface::Guts& bcsg = bcs.getGuts();

    //Initialize the grid to evaluate the surface at the knots and at 
    //the midpoints
    for(int i=0; i<(2*size-1); i++){
        xeV(i) = xmin + i*(xmax-xmin)/((Real)(2*size)-1.0);
        yeV(i) = ymin + i*(ymax-ymin)/((Real)(2*size)-1.0);
    }

    //Evaluate the surface at the knot points, and at the 
    //mid grid points and test if fV = A*aV holds
    Vec2 aXY;    
    for(int i=0; i<(2*size-1); i++){
        for(int j=0; j<(2*size-1); j++){
            aXY = Vec2(xeV(i), yeV(j));
            fV = bcsg.getPatchFunctionVector(aXY);
            aV = bcsg.getPatchBicubicCoefficients(aXY);
            fT = AM*aV;
            fVerr = fV-fT;
            
            //printf(" (%d,%d) ",i,j);
            //cout << fVerr.norm() << endl;

            //Due to the relatively large number of floating point
            //operations required, the tolerance needs to be 1e-12
            SimTK_TEST_EQ_TOL(fV,fT,1e-12);
        }
    }

}

/**
 This function will check that numerical derivatives of fx, fy, fxy, fxx,
 fyy, fxyy, fxxy, fxxx and fyyy match the values that the Bicubic surface
 function are returning. In addition, the surfaces that are defined by fx,
 fy, fxy, fxx, and fyy will be tested by continuity. Continuity is checked
 by moving a distance away from the knot point, computing the local derivative
 at the point along the direction towards the knot point, and then linearly
 extrapolating back to the knot point. If the linear extrapolation (of f, fx
 fy, fxy, fxx or fyy) matches the value of the function (f, fx, fy, fxy, fxx
 or fyy) at the knot point closely, then we can have some confidence that the
 surface is continuous. I say confidence rather than certaintity because for 
 certaintity we'd have to take the limit as that distance approches zero, and 
 that doesn't make sense in floating point.

 @params xmin: the minimum value of the x,y grid in the x dimension
 @params xmax: the maximum value of the x,y grid in the x dimension
 @params ymin: the minimum value of the y grid in the y dimension
 @params ymax: the maximum value of the y grid in the y dimension
 @params fcnType: An integer value [0-4] that picks an analytical function
                  to use for comparision purposes.
 @params smoothness: A value of 0 will make sure the patch goes through the 
                     desired points exactly. A value between 0 and 1 will
                     relax the surface.
 @params verbosePrint: true:  will print all of the detailed results for the
                              derivative comparisons, and the continuity 
                              checks
 @returns nothing
*/
void testBicubicConsistencyContinuity(Real xmin, Real xmax, Real ymin, 
             Real ymax, int fcnType, Real smoothness, bool verbosePrint){
    int size = 4;
    Real minstep = min((xmax-xmin),(ymax-ymin));
    Real dh = (minstep/(Real)size)/100.0;

    Vector xV(size), yV(size),dxV(4),dyV(4), tmpV(4), aXY(2);
    Matrix zM(size,size);


    //Initialize the 4x4 grid with a non-even grid spacing
    Real spacingX =1*(xmax-xmin)/((Real)size-1.0);
    Real spacingY =1*(ymax-ymin)/((Real)size-1.0);

    for(int i=0; i<size; i++){
        xV(i) = xmin + i*(xmax-xmin)/((Real)size-1.0);
        yV(i) = xmin + i*(ymax-ymin)/((Real)size-1.0);        
    }

    //Adjust the interior points a little bit to make
    //the spacing of the grid non-even. This will test
    //that BicubicSurface correctly handling the stretching 
    //of each individual patch correctly.
    for(int i=1; i<size-1;i++){
        xV(i) = xV(i) + 0.1*spacingX*pow(-1.0,i);
        yV(i) = yV(i) + 0.1*spacingY*pow(-1.0,i);
    }

    if(verbosePrint==true){
        cout << "X Spacing: " << xV << endl;
        cout << "Y Spacing: " << yV << endl;
    }

    for(int i=0; i<size;i++){
        for(int j=0; j<size; j++){
            tmpV = getAnalyticFunction(xV(i),yV(j),fcnType);
            zM(i,j) = tmpV(0);
        }
    }
    
    //Create the bicubic surface
    BicubicSurface bcs(xV, yV, zM, smoothness);
    BicubicFunction bcsf(bcs);

    //Initialize the vectors dxV and dyV to be near an interior knot
    //with the inner patch a distance h away from the knot point
    //and the second patch a distance h+dx away from the knot point

    int tsize = 17;
    Real tsizeh = floor((Real)tsize/2.0);    

    Matrix meshX(tsize,tsize), meshY(tsize,tsize);

    aXY(0) = xV(1);
    aXY(1) = yV(1);

    //Set up all of the partial derivative vectors required for the bench mark
    Array_<int> derivX(1);
    Array_<int> derivY(1);
    Array_<int> derivXY(2);
    Array_<int> derivXX(2);
    Array_<int> derivYY(2);
    Array_<int> derivXXY(3);
    Array_<int> derivXYY(3);
    Array_<int> derivXXX(3);
    Array_<int> derivYYY(3);
    Array_<int> deriv4X(4);
    Array_<int> deriv4Y(4);

    derivX[0] = 0;
        derivY[0] = 1;
    derivXY[0]= 0;
    derivXY[1]= 1;
        derivXX[0]= 0;
        derivXX[1]= 0;
    derivYY[0]= 1;
    derivYY[1]= 1;
        derivXXY[0]= 0;
        derivXXY[1]= 0;
        derivXXY[2]= 1;
    derivXYY[0]= 0;
    derivXYY[1]= 1;
    derivXYY[2]= 1;
        derivXXX[0]= 0;
        derivXXX[1]= 0;
        derivXXX[2]= 0;
    derivYYY[0]= 1;
    derivYYY[1]= 1;
    derivYYY[2]= 1;
    for(int i=0;i<4;i++){
        deriv4X[i]=0;
        deriv4Y[i]=1;
    }

    //Function computed derivatives
    Matrix bcsF(tsize,tsize),    bcsFx(tsize,tsize),     bcsFy(tsize,tsize);
    Matrix bcsFxy(tsize,tsize),  bcsFxx(tsize,tsize),    bcsFyy(tsize,tsize);
    Matrix bcsFxxy(tsize,tsize), bcsFxyy(tsize,tsize),   bcsFxxx(tsize,tsize);
    Matrix bcsFyyy(tsize,tsize), bcsF4x(tsize,tsize),    bcsF4y(tsize,tsize);

    //Numerically computed derivatives
    Matrix                       numFx(tsize,tsize),     numFy(tsize,tsize);
    Matrix numFxy(tsize,tsize),  numFxx(tsize,tsize),    numFyy(tsize,tsize);
    Matrix numFxxy(tsize,tsize), numFxyy(tsize,tsize),   numFxxx(tsize,tsize);
    Matrix numFyyy(tsize,tsize);

    aXY(0) = xV(1);
    aXY(1) = yV(1);

    //Sample the surface about aXY over a 17x17 grid
    for(int i=0;i<tsize;i++){
        for(int j=0;j<tsize;j++){
            meshX(i,j) = (xV(1) - tsizeh*dh) + dh*i;
            meshY(i,j) = (yV(1) - tsizeh*dh) + dh*j;
            aXY(0) = meshX(i,j);
            aXY(1) = meshY(i,j);
            
            bcsF(i,j) = bcsf.calcValue(aXY);
            
            bcsFx(i,j)= bcsf.calcDerivative(derivX,aXY);
            bcsFy(i,j)= bcsf.calcDerivative(derivY,aXY);
            
            bcsFxy(i,j)= bcsf.calcDerivative(derivXY,aXY);
            bcsFxx(i,j)= bcsf.calcDerivative(derivXX,aXY);
            bcsFyy(i,j)= bcsf.calcDerivative(derivYY,aXY);

            bcsFxxy(i,j)= bcsf.calcDerivative(derivXXY,aXY);
            bcsFxyy(i,j)= bcsf.calcDerivative(derivXYY,aXY);
            bcsFxxx(i,j)= bcsf.calcDerivative(derivXXX,aXY);
            bcsFyyy(i,j)= bcsf.calcDerivative(derivYYY,aXY);

            //Should be zero, just testing.
            bcsF4x(i,j) = bcsf.calcDerivative(deriv4X,aXY);
            bcsF4y(i,j) = bcsf.calcDerivative(deriv4Y,aXY);
        }
    }

    //Now compute the equivalent numerical derivatives using
    //central differences on the values in bcsF
    for(int i=0;i<tsize;i++){
        numFx(i)    = getCentralDifference(meshX(i),    bcsF(i),    true);
        numFxx(i)   = getCentralDifference(meshX(i),    numFx(i),   true);
        numFxxx(i)  = getCentralDifference(meshX(i),    numFxx(i),  true);

        numFy[i]    = ~getCentralDifference(~meshY[i],    ~bcsF[i],    true);
        numFyy[i]   = ~getCentralDifference(~meshY[i],    ~numFy[i],   true);
        numFyyy[i]  = ~getCentralDifference(~meshY[i],    ~numFyy[i],  true);
    }
    for(int i=0;i<tsize;i++){
        numFxy[i]   = ~getCentralDifference(~meshY[i],    ~numFx[i],    true);
        numFxxy[i]  = ~getCentralDifference(~meshY[i],    ~numFxx[i],   true);
    }
    for(int i=0;i<tsize;i++){
        numFxyy[i]  = ~getCentralDifference(~meshY[i],    ~numFxy[i],    true); 
    }       

    Real tol1 = dh;
    Real tol2 = dh*10;
    Real tol3 = dh*100;
    Vector dirXY(2);

    
    for(int i=3;i<tsize-3;i++){
        for(int j=3;j<tsize-3;j++){

            //1. Now compare the inner 10x10 numerical values 
            //   for each of the derivatives to the values computed 
            //   by the bicubic function
            if(verbosePrint==true){
                printf("\n\nCheck Derivatives (i,j): %d, %d", i,j);
                printf("\nbcs: fx:%f fy:%f" , bcsFx(i,j),bcsFy(i,j));
                printf("\nnum: fx:%f fy:%f" , numFx(i,j),numFy(i,j));
                printf("\n\n|bcs: fxy:%f fxx:%f fyy:%f ",
                    bcsFxy(i,j), bcsFxx(i,j), bcsFyy(i,j));
                printf("\nnum: fxy:%f fxx:%f fyy:%f ",
                    numFxy(i,j), numFxx(i,j), numFyy(i,j));
                printf("\n\n|bcs: fxxy:%f fxyy:%f fxxx:%f fyyy:%f",
                    bcsFxxy(i,j), bcsFxyy(i,j), bcsFxxx(i,j), bcsFyyy(i,j));
                printf("\nnum: fxxy:%f fxyy:%f fxxx:%f fyyy:%f",
                    numFxxy(i,j), numFxyy(i,j), numFxxx(i,j), numFyyy(i,j));
            }
           
            SimTK_TEST_EQ_TOL(bcsFx(i,j),numFx(i,j),tol1);
            SimTK_TEST_EQ_TOL(bcsFy(i,j),numFy(i,j),tol1);
            
            SimTK_TEST_EQ_TOL(bcsFxy(i,j),numFxy(i,j),tol2);
            SimTK_TEST_EQ_TOL(bcsFxx(i,j),numFxx(i,j),tol2);
            SimTK_TEST_EQ_TOL(bcsFyy(i,j),numFyy(i,j),tol2);

            /*The numerical 3rd derivatives will not match at the boundaries
            between patches. They are discontinuous in this region in the 
            formulation, and make the numerical derivatives around these 
            boundaries poorly estimated.*/

            //2. Continuity testing:
            //Test that a linear extrapolation from the current location
            //to the knot point matches the value of the knot point
            if(j != tsizeh || i != tsizeh){
                
                dirXY(0) = meshX(i,j)-meshX(8,8);
                dirXY(1) = meshY(i,j)-meshY(8,8);

                Real dist = pow(dirXY(0)*dirXY(0) + dirXY(1)*dirXY(1),0.5);

                //Test for surface continuity
                Real f0 = bcsF(i,j) -(bcsFx(i,j)*dirXY(0) 
                                      + bcsFy(i,j)*dirXY(1));
                Real err0 =f0-bcsF(8,8);
                Real errR0= abs(err0)/( abs(bcsF(8,8)) + 1e-10);
                
                //Test for fx derivative continuity
                Real f1x = bcsFx(i,j) -(bcsFxx(i,j)*dirXY(0));
                Real err1x =f1x-bcsFx(8,8);
                Real errR1x= abs(err1x)/( abs(bcsFx(8,8)) + 1e-10);

                //Test for fy derivative continity
                Real f1y = bcsFy(i,j) -(bcsFyy(i,j)*dirXY(1));
                Real err1y =f1y-bcsFy(8,8);
                Real errR1y= abs(err1y)/( abs(bcsFy(8,8)) + 1e-10);

                //Test for fxx derivative continuity
                Real f2x = bcsFxx(i,j) -(bcsFxxx(i,j)*dirXY(0));
                Real err2x =f2x-bcsFxx(8,8);
                Real errR2x= abs(err2x)/( abs(bcsFxx(8,8)) + 1e-10);

                //Test for fyy derivative continuity
                Real f2y = bcsFyy(i,j) -(bcsFyyy(i,j)*dirXY(1));
                Real err2y =f2y-bcsFyy(8,8);
                Real errR2y= abs(err2y)/( abs(bcsFyy(8,8)) + 1e-10);

                //Test for fxy derivative continuity
                Real fxy = bcsFxy(i,j) -(bcsFxxy(i,j)*dirXY(0) + bcsFxyy(i,j)*dirXY(1));
                Real errxy =fxy-bcsFxy(8,8);
                Real errRxy= abs(errxy)/( abs(bcsFxy(8,8)) + 1e-10);


                if(verbosePrint==true){
                    printf("\n\nCheck Continuity (i,j): %d, %d", i,j);
                    printf("\nf(x,y)  : %f num f  : %f  errR: %f" 
                                                , bcsF(8,8),  f0,  errR0);  
                    printf("\nfx(x,y) : %f num fx : %f  errR: %f" 
                                                , bcsFx(8,8), f1x, errR1x); 
                    printf("\nfy(x,y) : %f num fx : %f  errR: %f" 
                                                , bcsFy(8,8), f1y, errR1y); 
                    printf("\nfxx(x,y): %f num fxx: %f  errR: %f" 
                                                , bcsFxx(8,8),f2x, errR2x); 
                    printf("\nfyy(x,y): %f num fyy: %f  errR: %f" 
                                                , bcsFyy(8,8),f2y, errR2y); 
                    printf("\nfxy(x,y): %f num fxy: %f  errR: %f" 
                                                , bcsFxy(8,8),fxy, errRxy); 
                }
                   


                SimTK_TEST_EQ_TOL(errR0,0, dh);
                SimTK_TEST_EQ_TOL(errR1x,0,dh*5);
                SimTK_TEST_EQ_TOL(errR1y,0,dh*5);
                SimTK_TEST_EQ_TOL(errR2x,0,dh*5);
                SimTK_TEST_EQ_TOL(errR2y,0,dh*5);
                SimTK_TEST_EQ_TOL(errRxy,0,dh*10);
            }

        }
    }

}

/**
 This test function will create a bicubic surface and then test that 
 a version of this surface initialized using the copy constructor and
 the equal operator returns the same values over the surface as the original
*/
void testCopyConstEqOp(){
    int fcnType = 4;
    Real xmin = 0;
    Real xmax = 2*Pi;
    Real ymin = 0;
    Real ymax = Pi;
    Real smoothness = 0.1;

    int size = 4;
    Real minstep = min((xmax-xmin),(ymax-ymin));
    Real dh = (minstep/(Real)size)/100.0;

    Vector xV(size), yV(size),dxV(4),dyV(4), tmpV(4), aXY(2);
    Matrix zM(size,size);


    //Initialize the 4x4 grid with a non-even grid spacing
    Real spacingX =1*(xmax-xmin)/((Real)size-1.0);
    Real spacingY =1*(ymax-ymin)/((Real)size-1.0);

    for(int i=0; i<size; i++){
        xV(i) = xmin + i*(xmax-xmin)/((Real)size-1.0);
        yV(i) = xmin + i*(ymax-ymin)/((Real)size-1.0);        
    }

    //Adjust the interior points a little bit to make
    //the spacing of the grid non-even. This will test
    //that BicubicSurface correctly handling the stretching 
    //of each individual patch correctly.
    for(int i=1; i<size-1;i++){
        xV(i) = xV(i) + 0.1*spacingX*pow(-1.0,i);
        yV(i) = yV(i) + 0.1*spacingY*pow(-1.0,i);
    }

    /*if(verbosePrint==true){
        cout << "X Spacing: " << xV << endl;
        cout << "Y Spacing: " << yV << endl;
    }*/

    for(int i=0; i<size;i++){
        for(int j=0; j<size; j++){
            tmpV = getAnalyticFunction(xV(i),yV(j),fcnType);
            zM(i,j) = tmpV(0);
        }
    }
    
    //Create the bicubic surface
    BicubicSurface bcs(xV, yV, zM, smoothness);
    BicubicSurface bcsCC(bcs);
    BicubicSurface bcsEQOP;
    bcsEQOP = bcs;

    // Extract the implementation objects so we can look at the internals.
    const BicubicSurface::Guts& bcsg     = bcs.getGuts();
    const BicubicSurface::Guts& bcsCCg   = bcsCC.getGuts();
    const BicubicSurface::Guts& bcsEQOPg = bcsEQOP.getGuts();

    // These should all be the same underlying object, and the reference
    // count should be 3.
    SimTK_TEST(&bcsCCg == &bcsg);
    SimTK_TEST(&bcsEQOPg == &bcsg);
    SimTK_TEST(bcsg.getReferenceCount() == 3);

    // Create Function objects referencing the surface(s).
    BicubicFunction bcsf(bcs);
    BicubicFunction bcsCCf(bcs);
    BicubicFunction bcsEQOPf(bcs);

    // Reference count should now be 6.
    SimTK_TEST(bcsg.getReferenceCount() == 6);


    // These tests are meaningless now if the above ones succeed, since
    // obviously if they are the same object they will produce the same info!

    //Just to be extra sure, we'll actually check some values
    //computed from each of these different surfaces as well
    Real deltaX = (xmax-xmin)/15;
    Real deltaY = (ymax-ymin)/15;
    Array_<int> dX(1);
    Array_<int> dY(1);
    Array_<int> dXY(2);
    Array_<int> dXX(2);
    Array_<int> dYY(2);
    Array_<int> dXXY(3);
    Array_<int> dXYY(3);
    Array_<int> dXXX(3);
    Array_<int> dYYY(3);
    Array_<int> d4X(4);
    Array_<int> d4Y(4);

    dX[0] = 0;
        dY[0] = 1;
    dXY[0]= 0;
    dXY[1]= 1;
        dXX[0]= 0;
        dXX[1]= 0;
    dYY[0]= 1;
    dYY[1]= 1;
        dXXY[0]= 0;
        dXXY[1]= 0;
        dXXY[2]= 1;
    dXYY[0]= 0;
    dXYY[1]= 1;
    dXYY[2]= 1;
        dXXX[0]= 0;
        dXXX[1]= 0;
        dXXX[2]= 0;
    dYYY[0]= 1;
    dYYY[1]= 1;
    dYYY[2]= 1;

for(int i=0;i<16;i++){
    aXY(0) = xmin + i*deltaX;
    for(int j=0;j<16;j++){
        aXY(1) = ymin + j*deltaY;
        SimTK_TEST_EQ(bcsf.calcValue(aXY),  bcsCCf.calcValue(aXY));
        SimTK_TEST_EQ(bcsf.calcValue(aXY),bcsEQOPf.calcValue(aXY));

        SimTK_TEST_EQ(bcsf.calcDerivative(dX,aXY),  bcsCCf.calcDerivative(dX,aXY));
        SimTK_TEST_EQ(bcsf.calcDerivative(dX,aXY),bcsEQOPf.calcDerivative(dX,aXY));

        SimTK_TEST_EQ(bcsf.calcDerivative(dY,aXY),  bcsCCf.calcDerivative(dY,aXY));
        SimTK_TEST_EQ(bcsf.calcDerivative(dY,aXY),bcsEQOPf.calcDerivative(dY,aXY));

        SimTK_TEST_EQ(bcsf.calcDerivative(dXY,aXY),  bcsCCf.calcDerivative(dXY,aXY));
        SimTK_TEST_EQ(bcsf.calcDerivative(dXY,aXY),bcsEQOPf.calcDerivative(dXY,aXY));

        SimTK_TEST_EQ(bcsf.calcDerivative(dXXY,aXY),  bcsCCf.calcDerivative(dXXY,aXY));
        SimTK_TEST_EQ(bcsf.calcDerivative(dXXY,aXY),bcsEQOPf.calcDerivative(dXXY,aXY));

        SimTK_TEST_EQ(bcsf.calcDerivative(dXYY,aXY),  bcsCCf.calcDerivative(dXYY,aXY));
        SimTK_TEST_EQ(bcsf.calcDerivative(dXYY,aXY),bcsEQOPf.calcDerivative(dXYY,aXY));

        SimTK_TEST_EQ(bcsf.calcDerivative(dXXX,aXY),  bcsCCf.calcDerivative(dXXX,aXY));
        SimTK_TEST_EQ(bcsf.calcDerivative(dXXX,aXY),bcsEQOPf.calcDerivative(dXXX,aXY));

        SimTK_TEST_EQ(bcsf.calcDerivative(dYYY,aXY),  bcsCCf.calcDerivative(dYYY,aXY));
        SimTK_TEST_EQ(bcsf.calcDerivative(dYYY,aXY),bcsEQOPf.calcDerivative(dYYY,aXY));
    }
}

}

void testHint() {
    const Real xData[4] = { .1, 1, 2, 10 };
    const Real yData[5] = { -3, -2, 0, 1, 3 };
    const Real fData[] = { 1,   2,   3,   4,   5,
                           1.1, 2.1, 3.1, 4.1, 5.1,
                           1,   2,   3,   4,   5,
                           1.2, 2.2, 3.2, 4.2, 5.2 };
    const Vector x(4,   xData);
    const Vector y(5,   yData);
    const Matrix f(4,5, fData);
    BicubicSurface surf(x, y, f, 0); // not smoothed

    SimTK_TEST(surf.getNumAccesses() == 0);

    BicubicSurface::PatchHint hint;
    Real val = surf.calcValue(Vec2(.5, .5), hint);
    SimTK_TEST(surf.getNumAccesses() == 1);
    val = surf.calcValue(Vec2(.5, .5), hint); // should be free
    SimTK_TEST(surf.getNumAccesses() == 2);
    SimTK_TEST(surf.getNumAccessesSamePoint() == 1);

    val = surf.calcValue(Vec2(.50001, .50002), hint);
    SimTK_TEST(surf.getNumAccessesSamePatch() == 1);

    val = surf.calcValue(Vec2(1.5, -1), hint);
    SimTK_TEST(surf.getNumAccessesNearbyPatch() == 1);

    // This should report "same patch" rather than "same point" because
    // derivative info hasn't been calculated yet.
    Array_<int> deriv1(1, 1), deriv2(2, 0); // fy, fxx
    val = surf.calcDerivative(deriv2, Vec2(1.5, -1), hint);
    SimTK_TEST(surf.getNumAccessesSamePatch() == 2);

    // When 2nd deriv info is calculated we get 1st deriv also. So now
    // we should get "same point" even though we haven't asked for this yet.
    val = surf.calcDerivative(deriv1, Vec2(1.5, -1), hint);
    SimTK_TEST(surf.getNumAccessesSamePoint() == 2);

}

int main() {
    //Evaluate the bicubic surface interpolation against an analytical 
    //function. Throw an error if the values of the function are different
    //at the knot points, or different within tolerance at the mid grid points
    SimTK_START_TEST("Testing Bicubic Interpolation");
        SimTK_SUBTEST(testHint);

    cout << "\n---------------------------------------------"<< endl;
    cout<< "\n\nANALYTICAL FUNCTION COMPARISON:" << endl;
    testBicubicAgainstAnalyticFcn(0.0, 1.0, 0.0, 1.0,9,0,false,false);
    testBicubicAgainstAnalyticFcn(0.0, 1.0, 0.0, 1.0,9,1,false,false);
    testBicubicAgainstAnalyticFcn(0.0, 1.0, 0.0, 1.0,9,2,false,false);
    testBicubicAgainstAnalyticFcn(0.0, 1.0, 0.0, 1.0,9,3,false,false);
    testBicubicAgainstAnalyticFcn(0.0, 1.0, 0.0, 1.0,9,4,false,false);
    printf("\n\n*Test Passed*. Constructor with x,y,f,fx,fy,fxy specified,"
            " \n\tSmoothness parameter %f tested\n"
                "\tAdditional smoothness parameters not tested because"
                "\n\tsurface will not pass through the knot points",Real(0));
    cout << "\n---------------------------------------------"<< endl;

    cout << "\n---------------------------------------------"<< endl;
    cout << "\n\nBICUBIC COEFFICIENT VALIDATION:" << endl;
    cout << "  Testing that the bicubic interpolation coefficients" <<endl;
    cout << " are being solved correctly by asserting fV - A*aV = 0"<<endl;
    testBicubicCoefficients(      0.0, 1.0, 0.0, 1.0,  3, 0.0);          
    testBicubicCoefficients(      0.0, 1.0, 0.0, 1.0,  3, 0.5);
    printf("\n\n*Test Passed*. Constructor with x,y,f specified,"
           " \n\tSmoothness parameter %f and %f tested",(Real)0.0,(Real)0.5);
    cout << "\n---------------------------------------------"<< endl;

    cout << "\n\n---------------------------------------------"<< endl;
    cout << "\n\nBICUBIC DERIVATIVE & CONTINUITY TESTING:" <<endl;
    cout << " 1. Derivative are tested for consistency by ensuring that" << endl;
    cout << "    numerical derivatives of f(x,y) match values returned " << endl;
    cout << "    by the function." << endl;
    cout << "    Partial derivatives tested: fx,fy,fxy,fxx,fyy" << endl; 
    cout << "\n 2. Continuity is tested by asserting that a linear extrapolation" << endl;
    cout << "    from a a point near a knot is equal to the value of the surface" << endl;
    cout << "    of f(x,y) at the knot. Surfaces tested f, fx, fy, fxy, fxx, fyy." << endl;
    testBicubicConsistencyContinuity(        0.0, 1.0, 0.0, 1.0,  3, 0.0, false);
    testBicubicConsistencyContinuity(        0.0, 1.0, 0.0, 1.0,  3, 0.5, false);
    printf("\n\n*Test Passed*. Constructor with x,y,f specified,"
            " \n\tSmoothness parameter %f and %f tested",(Real)0.0,(Real)0.5);
    cout << "\n---------------------------------------------"<< endl;

    cout << "\n\n---------------------------------------------"<< endl;
    cout << "\n\nCOPY CONSTRUCTOR AND = OPERATOR TESTING:" <<endl;
    cout <<" Tested by using the copy contructor and equal operator" << endl;
    cout <<" and comparing the values of the internally stored matrices" << endl;
    cout <<" of x,y,f,fx,fy,fxy between the different surfaces, and then" <<endl;
    cout <<" comparing values of f,fx,fy,fxy,fxx,fyy,fxyy,fxxy,fxxx,fyyy" <<endl;
    cout <<" between the different surface objects across the patch"<<endl;
    testCopyConstEqOp();
    cout <<" *Test Passed*." << endl;
    cout << "\n---------------------------------------------"<< endl;

    SimTK_END_TEST();
}

