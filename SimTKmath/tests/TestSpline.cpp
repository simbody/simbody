/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Matthew Millard (the testNaturalCubicSpline code)            *
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
#include <vector> // use some std::vectors to test Array_ interoperability
#include <iostream>
#include <fstream>
#include <cstdio>

using namespace SimTK;
using namespace std;

const Real TESTTOL = 1e-9;

void testSpline() {
    Vector_<Vec3> coeff(5);
    coeff[0] = Vec3(0, 1, 2);
    coeff[1] = Vec3(1, 4, 1);
    coeff[2] = Vec3(2, 2, 20);
    coeff[3] = Vec3(1, -1, 2);
    coeff[4] = Vec3(0, 0, 1);
    Vector x(Vec5(0, 1, 2, 5, 10));
    
    // Create a linear spline, and verify that it interpolates linearly between
    // the control points.
    
    Spline_<Vec3> spline(1, x, coeff);
    for (int i = 0; i < x.size(); ++i)
        SimTK_TEST_EQ(coeff[i], spline.calcValue(Vector(1, x[i])));
    std::vector<int> deriv;
    deriv.push_back(0);
    for (int i = 0; i < x.size()-1; ++i) {
        for (int j = 0; j < 10; ++j) {
            Real fract = (i+1.0)/12.0;
            Real t = x[i]+fract*(x[i+1]-x[i]);
            SimTK_TEST_EQ_TOL(spline.calcValue(Vector(1, t)), 
                              coeff[i]+fract*(coeff[i+1]-coeff[i]), TESTTOL);
            SimTK_TEST_EQ_TOL(spline.calcDerivative(deriv, Vector(1, t)), 
                              (coeff[i+1]-coeff[i])/(x[i+1]-x[i]),TESTTOL);
        }
    }
    
    // Create a cubic spline and verify the derivative calculations.
    
    spline = Spline_<Vec3>(3, x, coeff);
    Real delta = 1e-10;
    for (int i = 0; i < x.size()-1; ++i) {
        for (int j = 0; j < 10; ++j) {
            Real fract = (i+1.0)/12.0;
            Real t = x[i]+fract*(x[i+1]-x[i]);
            Vec3 value1 = spline.calcValue(Vector(1, t-delta));
            Vec3 value2 = spline.calcValue(Vector(1, t+delta));
            SimTK_TEST_EQ_TOL(spline.calcDerivative(deriv, Vector(1, t)), 
                              (value2-value1)/(2*delta), 1e-4);
        }
    }
}

void testSplineFitter() {
    Real stddev = 0.5;
    int n = 100;
    Random::Gaussian random(0.0, stddev);
    Vector x(n);
    Vector_<Vec3> truey(n);
    Vector_<Vec3> y(n);
    for (int i = 0; i < x.size(); ++i) {
        x[i] = i*0.1;
        truey[i] = Vec3(sin(x[i]), 3.0*sin(2*x[i]), cos(x[i]));
        y[i] = truey[i] + Vec3(random.getValue(),random.getValue(),random.getValue());
    }
    SplineFitter<Vec3> fitter = SplineFitter<Vec3>::fitFromGCV(3, x, y);
    Spline_<Vec3> spline1 = fitter.getSpline();
    
    // The fitting should have reduced the error.
    
    Vec3 originalError = mean(abs(y-truey));
    Vec3 fittedError = mean(abs(spline1.getControlPointValues()-truey));
    SimTK_TEST(fittedError[0] < originalError[0]);
    SimTK_TEST(fittedError[1] < originalError[1]);
    SimTK_TEST(fittedError[2] < originalError[2]);
    
    // If we perform the fitting again, explicitly specifying the same value for
    // the smoothing parameter, it should produce identical results.   
    SimTK_TEST_EQ_TOL(SplineFitter<Vec3>::fitForSmoothingParameter
                                    (3, x, y, fitter.getSmoothingParameter())
                                        .getSpline().getControlPointValues(), 
                      spline1.getControlPointValues(), 
                      TESTTOL);
    
    // Likewise, specifying the same number of degrees of freedom should produce
    // identical results.
    SimTK_TEST_EQ_TOL(SplineFitter<Vec3>::fitFromDOF
                                    (3, x, y, fitter.getDegreesOfFreedom())
                                        .getSpline().getControlPointValues(), 
                      spline1.getControlPointValues(),
                      TESTTOL);
    
    // If we specify a smoothing parameter of 0, it should exactly reproduce 
    // the original data.
    Spline_<Vec3> nosmoothing = SplineFitter<Vec3>::fitForSmoothingParameter
                                                    (3, x, y, 0.0).getSpline();
    for (int i = 0; i < x.size(); ++i)
        SimTK_TEST_EQ_TOL(y[i], nosmoothing.calcValue(Vector(1, x[i])),TESTTOL);
}

void testRealSpline() {
    Vector coeff(5);
    coeff[0] = 0;
    coeff[1] = 1;
    coeff[2] = 2;
    coeff[3] = 1;
    coeff[4] = 0;
    Vector x(Vec5(0, 1, 2, 5, 10));

    // Create a linear spline, and verify that it interpolates linearly between
    // the control points.

    Spline spline(1, x, coeff);
    for (int i = 0; i < x.size(); ++i)
        SimTK_TEST_EQ_TOL(coeff[i], spline.calcValue(Vector(1, x[i])), TESTTOL);
    Array_<int> deriv;
    deriv.push_back(0);
    for (int i = 0; i < x.size()-1; ++i) {
        for (int j = 0; j < 10; ++j) {
            Real fract = (i+1.0)/12.0;
            Real t = x[i]+fract*(x[i+1]-x[i]);
            SimTK_TEST_EQ_TOL(spline.calcValue(Vector(1, t)), 
                              coeff[i]+fract*(coeff[i+1]-coeff[i]), TESTTOL);
            SimTK_TEST_EQ_TOL(spline.calcDerivative(deriv, Vector(1, t)), 
                              (coeff[i+1]-coeff[i])/(x[i+1]-x[i]), TESTTOL);
        }
    }
    SimTK_TEST_EQ_TOL(1, spline.getControlPointValues()[1], TESTTOL);

    // Try using a SplineFitter.

    SplineFitter<Real> fitter = SplineFitter<Real>::fitFromGCV(3, x, coeff);
    Spline spline2 = fitter.getSpline();
    SimTK_TEST_EQ_TOL(3, spline2.getSplineDegree(),TESTTOL);
}

//MM bits added to test the numerical accuracy of the natural cubic splines.
/**
* This function computes a standard central difference dy/dx. 
* If extrap_endpoints is set to 1, then the derivative at the 
* end points is estimated by linearly extrapolating the dy/dx 
* values beside the end points
*
* @param x domain vector
* @param y range vector
& @param extrap_endpoints: 
*   (false)   Endpoints of the returned vector will be zero, 
*             because a central difference is undefined at 
*             these endpoints
*    (true)   Endpoints are computed by linearly extrapolating 
*             using a first difference from the neighboring 2
*             points
*                                    
* @returns dy/dx computed using central differences
*/
Vector getCentralDifference(Vector x, Vector y, 
                                   bool extrap_endpoints){
    Vector dy(x.size());
    double dx1,dx2;
    double dy1,dy2;
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
* This function will print cvs file of the column vector 
* col0 and the matrix data
*
* @params col0: A vector that must have the same number of rows 
*               as the data matrix. This column vector is 
*               printed as the first column
* @params data: A matrix of data
* @params filename: The name of the file to print
*/
void printMatrixToFile(Vector col0,Matrix data, string filename){
    ofstream datafile;
    datafile.open(filename.c_str());

    for(int i = 0; i < data.nrow(); i++){
        datafile << col0(i) << ",";
        for(int j = 0; j < data.ncol(); j++){
            if(j<data.ncol()-1)
                datafile << data(i,j) << ",";
            else
                datafile << data(i,j) << "\n";
        }    
    }
    datafile.close();
} 

/**
* This function will compute the value and first two 
* derivatives of an analytic function at the point x. 
*
* @params x       the input value
* @params fcnType the function to compute. There are 
*                 currently 5 choices (see below)
* @returns Vector: a 3x1 vector of the value, 
*                 first derivative and second derivative
*/
Vector getAnalyticFunction(double x,int fcnType){
    Vector fdF(3);
    fdF = -1;

    switch(fcnType){
        case 0:     //f(x) = 0;
            fdF = 0;
            break;
        case 1:     //f(x) = 2*x
            fdF(0) = 2*x;   //f
            fdF(1) = 2;     //fx
            fdF(2) = 0;
            break;
        case 2:     //f(x) = x^2
            fdF(0) = x*x;   //f
            fdF(1) = 2*x;   //fx
            fdF(2) = 2;
            break;
        case 3:     //f(x) = 2*x + x*x;
            //f
            fdF(0) = 2*x + x*x;
            fdF(1) = 2 + 2*x;
            fdF(2) = 2;
            break;
        case 4:     //f(x)  =2*x + x*x + 5*x*x*x
            fdF(0) = 2*x + x*x + 5*x*x*x;
            fdF(1) = 2 + 2*x + 15*x*x;
            fdF(2) = 2 + 30*x;
            break;
        case 5:     //fx(x) = sin(x)
            fdF(0) = sin(x);
            fdF(1) = cos(x);
            fdF(2) = -sin(x);
            break;
        default:
            cout << "Invalid fcnType in testBicubicSurface.cpp: getAnayticFunction";
    }


    return fdF;
}
/**
* This function tests the accuracy of the natural cubic spline sp. 
* The accuracy of the spline is tested in the following manner:
*
*    a.    Spline must pass through the knots given
*             -Error between spline and input data at the knots 
*             (should be zero)
*    b.   The first derivatives are continuous at the knot points
*             -Error between the value of the first derivative at 
*             the knot point, and what a linear extrapolation would 
*             predict just to the left and right ofthe knot point. 
*             (should be zero, within a tolerace affected by the
*             step size in xD)
*    c.   The second derivatives are continuous at the knots points
*             -Error between the value of the numerically calculated 
*             derivative at the knot point, and what a linear
*             extrapolation would predict just to the left and 
*             right of the knot point. (should be zero, within a 
*             tolerace affected by the step size in xD)
*    d.  The second derivative is zero at the end points.
*             -Numerically calculated extrapolation of the 2nd 
*             derivative should be zero at the end points within
*             some tolerance
*
*/
Vector benchmarkNaturalCubicSpline
   (Function* sp, Vector xK, Vector yK, Vector xM,Vector xD,
    string name, bool print)
{
    int size = xK.size();
    int sizeD= xD.size();
    int sizeDK = xD.size()/(xK.size()-1);
    double deltaD = (xK(xK.size()-1)-xK(0))/xD.size();

    Matrix ysp_K(size,2),ysp_M(size-1,2),ysp_D(sizeD,4);
    Vector errVec(4);
    errVec = 1;
    ysp_K = 0;
    ysp_M = 0;
    ysp_D = 0;

    vector<int> derOrder(1);
    derOrder[0] = 0;
    
    

    ///////////////////////////////////////////
    //1. Evaluate the spline at the knots, the mid points and then a dense sample
    ///////////////////////////////////////////
        Vector tmpV1(1);
        double xVal=0;
        for(int i=0;i<size;i++){
            xVal = xK(i);
            tmpV1(0)=xK(i);
            ysp_K(i,0) = sp->calcValue(tmpV1);
            ysp_K(i,1) = sp->calcDerivative(derOrder,tmpV1);
        }            
        for(int i=0;i<size-1;i++){
            xVal = xM(i);
            tmpV1(0) = xM(i);
            ysp_M(i,0) = sp->calcValue(tmpV1);
            ysp_M(i,1) = sp->calcDerivative(derOrder,tmpV1);
        }
        for(int i=0;i<sizeD;i++){
            xVal = xD(i);
            tmpV1(0) = xD(i);
            ysp_D(i,0) = sp->calcValue(tmpV1);
            ysp_D(i,1) = sp->calcDerivative(derOrder,tmpV1);
        }

    //////////////////////////////////////
    //2.    Compute the second derivative of the spline (using central 
    //differences), and linearly  interpolate to get the end points. 
    //The end points should go to exactly zero because the second
    // derivative is linear in a cubic spline, as is the linear 
    // extrapolation
    //
    // Also compute the 3rd derivative using the same method. The 3rd 
    // derivative is required in the test to determine if the second 
    // derivative is continuous at the knots or not.
    //////////////////////////////////////

        ysp_D(2)    = getCentralDifference(xD, ysp_D(1),    true);
        ysp_D(3)    = getCentralDifference(xD, ysp_D(2),    true);

    //////////////////////////////////////
    //3. Now check to see if the splines meet the conditions of a 
    //natural cubic spline:
    //////////////////////////////////////

        Vector tmpK(size,size),tmpM(size-1,size-1);

        //* a.    Spline passes through all knot points given
                    tmpK = yK-ysp_K(0);
                errVec(0) = tmpK.norm();
            
        //    b. The first derivative is continuous at the knot points. 
        //    Apply a continuity test to the data points that define 
        //    the second derivative
        //
        //        Continuity test:    a linear extrapolation of first 
        //                            derivative of the curve in interest
        //                            on either side of the point in 
        //                            interest should equal the point in 
        //                            interest;

            double ykL,ykR,y0L,dydxL,y0R,dydxR = 0;
            for(int i=1; i<size-1; i++){
                y0L = ysp_D(i*sizeDK-1,1);
                y0R = ysp_D(i*sizeDK+1,1);
                dydxL = ysp_D(i*sizeDK-1,2); //Found using central differences
                dydxR = ysp_D(i*sizeDK+1,2); //Found using central differences
                ykL = y0L + dydxL*deltaD;
                ykR = y0R - dydxR*deltaD;
                errVec(1) = (ysp_D(i*sizeDK,1)-ykL)+(ysp_D(i*sizeDK,1)-ykR);
            }
            
            
        //    c. The second derivative is continuous at the knot points. 
        //    Apply a continuity test to the data points that define 
        //    the second derivative. This also tests if the first 
        //    derivative is smooth.
        //
        //        Continuity test:    a linear extrapolation of first
        //                            derivative of the curve in interest
        //                            on either side of the point in 
        //                            interest should equal the point in 
        //                            interest;
            for(int i=1; i<size-1; i++){
                y0L = ysp_D(i*sizeDK-1,2);
                y0R = ysp_D(i*sizeDK+1,2);
                dydxL = ysp_D(i*sizeDK-1,3); //Found using central differences
                dydxR = ysp_D(i*sizeDK+1,3); //Found using central differences
                ykL = y0L + dydxL*deltaD;
                ykR = y0R - dydxR*deltaD;
                errVec(2) = (ysp_D(i*sizeDK,2)-ykL)+(ysp_D(i*sizeDK,2)-ykR);
            }    

        //////////////////////////////////////
        //* d.    The second derivative is zero at the end points
        //////////////////////////////////////

        errVec(3) = abs(ysp_D(0,2)) + abs(ysp_D(sizeD-1,2));

        

        //////////////////////////////////////
        //print the data for analysis
        //////////////////////////////////////

        if(print==true){
            string fname = name;
            fname.append("_K.csv");
            printMatrixToFile(xK,    ysp_K,    fname);

            fname = name;
            fname.append("_M.csv");
            printMatrixToFile(xM,    ysp_M,    fname);

            fname = name;
            fname.append("_D.csv");
            printMatrixToFile(xD,    ysp_D,    fname);
        }

        return errVec;

}

/**
* The test works by seeing if the tested splines have the properties of 
* a natural cubic spline. To do so, this test file has several steps
*
* User Steps: Configure the script
*        a. Choose the function to be interpolated
*        b. Choose the location and number of knot points
*        c. Choose the density of a high resolution interpolation
*
* Test Script Steps:
* 0. Initialize the input vectors xK, xM, and xD for the knot locations, 
*    mid knot location and high resolution step locations respectively
* 1. Initialize the analytically computed output yK, yM and yD
* 2. Create each of the spline objects. 
* 3. Evaluate the numerical accuracy of the splines by calling 
*    testNaturalCubicSpline
*/
void testNaturalCubicSpline() {
    /////////////////////////////
    //Configuration Variables
    ////////////////////////////
        bool printToTerminal = false; //Setting this to true will print some 
                                      //useful data to the terminal
        bool printData = false;       //Set to true to print the knot, 
                                      //mid knot, and 
                                      //dense vector values, first derivatives, 
                                      //and second derivatives (for the splines) 
                                      //for analysis outside of this script.
        int fcnType = 5;    //Chooses what kind of analytical test function to 
                            //use to initialize and test the various 
                            //spline classes
        const int size =6;             //Number of knot points
        const int sizeDK = 100;        //Number of points per knot in 
                                       //the densely sampled vector
        int sizeD=sizeDK*(size-1);     //Number of points in a densly sampled 
                                       //interpolation

        //Domain vector variables
        double xmin,xmax,deltaX,deltaD;
        xmin = Pi/4;            //Value of first knot
        xmax = Pi/2;            //Value of the final knot

        deltaX = (xmax-xmin)/(size-1);    
        deltaD = (xmax-xmin)/(sizeD-1);
        double etime = 0;
    /////////////////////////////
    //Test Code body
    ////////////////////////////

        Matrix testResults(4,1); //This matrix stores the results of the 
                                        //5 tests in each row entry, for each 
                                        //of the 2 spline classes tested.
                                        // SimTK SplineFitter results are stored
                                        // in column 0
                                        //  OpenSim::NaturalCubicSpline results
                                        //  are stored in column 1
        //testResults.elementwiseAssign(0.0);
        testResults = -1;
        Vector tmpV1(1);

        //Generate initialization knot points (denoted by a 'K') 
        //        and the mid points (denoted by a 'M')        
        //        and for the densely sampled interpolation vector 
        //        (denoted by a 'D')
        Vector xK(size), xM(size-1), xD(sizeD);
        Matrix yK(size,3), yM(size-1,3), yD(sizeD,3);



    ///////////////////////////////////////////
    //0. Initialize the input vectors xK, xM and xD
    ///////////////////////////////////////////
        for (int i = 0; i < size; i++) {
            xK(i) = xmin + ((double)i)*deltaX;
            if(i<size-1){
                xM(i) = xmin + deltaX/(double)2 + ((double)i)*deltaX;
            }
        }
        for(int i = 0; i < sizeD; i++)
            xD(i) = xmin + deltaD*(double)i;

    ///////////////////////////////////////////        
    //1.    Initialize the analytic function vector data to interpolate
    //        Let the user know which function is being used
    ///////////////////////////////////////////
        if(printToTerminal==true){
            switch(fcnType){
                case 0:
                    cout << "f(x) = 0" <<endl;
                    break;
                case 1:
                    cout << "f(x) = 2*x" <<endl;
                    break;
                case 2:
                    cout << "f(x) = x^2" <<endl;
                    break;
                case 3:
                    cout << "f(x) = 2*x + x^2 " <<endl;
                    break;
                case 4:
                    cout << "f(x) = 2*x + x^2 + 5x^3 " <<endl;
                    break;
            }
        }


        //Get the function values at the knot points
        Vector tmp(3);
        tmp = 0;
        for(int i=0; i<size;i++){
            tmp = getAnalyticFunction(xK(i),fcnType);
            for(int k=0;k<3;k++)
                yK(i,k) = tmp(k);

        }
         Vector yKVal = yK(0);

        //Get the function y, dy, ddy at the mid points
        for(int i=0; i<size-1;i++){
            tmp = getAnalyticFunction(xM(i),fcnType);
            for(int k=0;k<3;k++)
                yM(i,k) = tmp(k);
        }
        //Get the function y, dy, ddy at the dense points
        for(int i=0; i<sizeD;i++){
            tmp = getAnalyticFunction(xD(i),fcnType);
            for(int k=0;k<3;k++)
                yD(i,k) = tmp(k);
        }

    ///////////////////////////////////////////
    //2. Create each of the splines
    ///////////////////////////////////////////

        //SplineFitter
        Vector sfDerivs1(xK.size());
        

        
        Spline_<Real> sTK = SplineFitter<Real>::fitForSmoothingParameter(3,xK,yKVal,0.0).getSpline();
            
    ///////////////////////////////////////////
    //3. Test the splines
    ///////////////////////////////////////////

        testResults(0) = benchmarkNaturalCubicSpline(&sTK, xK, yK(0), xM, xD, "simtk_splinefitter",true);
        if(printToTerminal==true){
            cout << "Test Result Matrix: 0 or small numbers pass" <<endl;
            cout << "  column 0: SplineFitter, column 1: OpenSim::NaturalCubicSpline" <<endl;
            cout << "  row 0: Passes through knots (tol 1e-14)" << endl;
            cout << "  row 1: First derivative is continuous and smooth (tol " << deltaD << ")" << endl;
            cout << "  row 2: Second derivative is continuous (tol " << 10*deltaD << ")" << endl;
            cout << "  row 3: Second derivative is zero at endpoints (tol "<< deltaD/10 <<")" << endl;
            cout << testResults << endl;
        }

    //////////////////////////////////////
    //4. Run numerical assertions on each test
    //////////////////////////////////////

        double tol = 0;
            
        for(int k=0;k<testResults.ncol();k++){
            for(int i=0;i<testResults.nrow();i++){                    
                switch(i){
                    case 0: //Equal at knots
                        tol = 1e-14;
                        break;
                    case 1: //Continuous 1st derivative
                        tol = deltaD;
                        break;
                    case 2: //Continuous 2nd derivative
                        tol = 10*deltaD;
                        break;    
                    case 3: //2nd derivative zero at end points
                        tol = deltaD/10;
                        break;
                    default:
                        cout << "testNCSpline: Invalid error type selected" << endl;
                }
                //cout << "Testing (i,k) " << i << " " << k << " tol " << tol << " \tval " << testResults(i,k) << endl;
                SimTK_TEST_EQ_TOL(testResults(i,k),0,tol);
            }
        }
        // getchar();

}

int main () {
    SimTK_START_TEST("TestSpline");
        SimTK_SUBTEST(testSpline);
        SimTK_SUBTEST(testSplineFitter);
        SimTK_SUBTEST(testRealSpline);
        SimTK_SUBTEST(testNaturalCubicSpline);
    SimTK_END_TEST();
}
