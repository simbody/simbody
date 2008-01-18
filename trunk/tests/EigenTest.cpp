/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Jack Middleton                                                   *
 * Contributors:                                                              *
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
 * This is a test program which uses the Eigen  class to compute 
 * eigen values and eigen vectors
 */

/*
The data for this test is from an example FORTRAN  program from the
Numerical Algorithms Group (NAG)
URL:http://www.nag.com/lapack-ex/lapack-ex.html


Solves for the eigen valus and vectors for the 
following system

   Ax = 0 where :



     0.35   0.45   -0.14   -0.17 
     0.09   0.07   -0.54    0.35 
A = -0.44  -0.33   -0.03    0.17 
     0.25  -0.32   -0.13    0.11 



SOLUTION = 
reciprocal condition number =  9.9E-01
 Error bound                 =  1.3E-16

 Eigenvector( 1)
 -6.5509E-01
 -5.2363E-01
  5.3622E-01
 -9.5607E-02

 Reciprocal condition number =  8.2E-01
 Error bound                 =  1.6E-16

 Eigenvalue( 2) = (-9.9412E-02, 4.0079E-01)

 Reciprocal condition number =  7.0E-01
 Error bound                 =  1.8E-16

 Eigenvector( 2)
 (-1.9330E-01, 2.5463E-01)
 ( 2.5186E-01,-5.2240E-01)
 ( 9.7182E-02,-3.0838E-01)
 ( 6.7595E-01, 0.0000E+00)

 Reciprocal condition number =  4.0E-01
 Error bound                 =  3.3E-16

 Eigenvalue( 3) = (-9.9412E-02,-4.0079E-01)

 Reciprocal condition number =  7.0E-01
 Error bound                 =  1.8E-16

 Eigenvector( 3)
 (-1.9330E-01,-2.5463E-01)
 ( 2.5186E-01, 5.2240E-01)
 ( 9.7182E-02, 3.0838E-01)
 ( 6.7595E-01,-0.0000E+00)

 Reciprocal condition number =  4.0E-01
 Error bound                 =  3.3E-16

 Eigenvalue( 4) = -1.0066E-01

 Reciprocal condition number =  5.7E-01
 Error bound                 =  2.3E-16

 Eigenvector( 4)
  1.2533E-01
  3.3202E-01
  5.9384E-01
  7.2209E-01

 Reciprocal condition number =  3.1E-01
 Error bound                 =  4.2E-16


estimated rank = 4

*/



//#define SimTK_USE_STATIC_LIBRARIES

#include "SimTKcommon.h"
#include "SimTKmath.h"

#include <cstdio>
#include <cassert>
#include <iostream>

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

static const double EPS = 0.00001;

using namespace SimTK;

using std::printf;
using std::cout;
using std::endl;

Real A[16] = {  0.35,  0.45,  -0.14,  -0.17,
                0.09,  0.07,  -0.54,   0.35,
               -0.44, -0.33,  -0.03,   0.17,
                0.25, -0.32,  -0.13,   0.11 };

typedef std::complex<double> cd;
cd expEigen[4] = { cd(0.79948,   0.0), 
                                      cd(-0.099412,  0.40079), 
                                      cd(-0.099412, -0.40079),
                                      cd(-0.10066,   0.0) };

cd expVectors[16] = { cd( -.65509, 0.0), cd( -.52363, 0.0), cd(  .53622, 0.0), cd( -.095607, 0.0),
                      cd(-.1933001,  .25463), cd( .2518601, -.52240), cd( .09718202,-.30838), cd( .67595,   0.000),
                      cd(-.1933001, -.25463), cd( .2518601,  .52240), cd( .09718202, .30838), cd( .67595,   -0.000),
                      cd( .12533, 0.0), cd( .33202, 0.0), cd( .59384, 0.0), cd( .72209, 0.0) };
 template <typename T> 
T absNormComplex( Vector_<std::complex<T> >& values, Vector_<std::complex<T> >& expected) {
   T norm = 0;
 
   for(int i=0;i<values.size(); i++ ) {
       norm += (fabs(values(i).real()) - fabs(expected(i).real())) * (fabs(values(i).real()) - fabs(expected(i).real())) + 
             (fabs(values(i).imag()) - fabs(expected(i).imag())) * (fabs(values(i).imag()) - fabs(expected(i).imag()));  
   } 

   return( sqrt(norm) );
}
template <typename T> 
T absNorm( Vector_<T>& values, Vector_<T>& expected) {
   T norm = 0;

   for(int i=0;i<values.size(); i++ ) {
       norm += (fabs(values(i)) - fabs(expected(i))) * (fabs(values(i)) - fabs(expected(i))) + 
              (fabs(values(i)) - fabs(expected(i))) * (fabs(values(i)) - fabs(expected(i)));  
   }
   return( sqrt(norm) );
}

int main () {
    double errnorm;
    try { 
           // Default precision (Real, normally double) test.

        Matrix a(4,4, A);
        Vector_<std::complex<double> > expectedValues(4);
        for(int i=0;i<4;i++) expectedValues[i] = expEigen[i];
        Matrix_<std::complex<double> > expectedVectors(4,4);
        for(int i=0;i<4;i++) for(int j=0;j<4;j++) expectedVectors(i,j) = expVectors[j*4+i];
        Vector_<std::complex<double> > values; // should get sized automatically to 4 by getAllEigenValuesAndVectors()
        Vector_<std::complex<double> > expectVec(4);
        Vector_<std::complex<double> > computeVec(4);
        Matrix_<std::complex<double> > vectors; // should get sized automatically to 4x4 by getAllEigenValuesAndVectors()

        Eigen  es(a);   // setup the eigen system 

        es.getAllEigenValuesAndVectors( values, vectors );  // solve for the eigenvalues and eigenvectors of the system 

        printf("\n *****  NON-SYMMETRIC:  ***** \n\n"  );
        cout << " Real SOLUTION: " << values << "  errnorm=" << absNormComplex(values,expectedValues) << endl;
        ASSERT(absNormComplex(values,expectedValues) < 0.001);

        cout << "Vectors = "  << endl;
        for(int i=0;i<4;i++) {
            computeVec = vectors(i); 
            expectVec = expectedVectors(i);

            errnorm =  absNormComplex( computeVec, expectVec );
            cout << computeVec << "  errnorm=" << errnorm << endl;
            ASSERT( errnorm < 0.00001 );

        }
  
        cout << endl << endl;

        Vector_<std::complex<float> > expectedValuesf(4);
        Matrix_<std::complex<float> > expectedVectorsf(4,4);
        Vector_<std::complex<float> > expectVecf(4);
        Vector_<std::complex<float> > computeVecf(4);
        Matrix_<std::complex<float> > vectorsf; // should get sized automatically to 4x4 by getAllEigenValuesAndVectors()
        Vector_<std::complex<float> > valuesf; // should get sized automatically to 4 by getAllEigenValuesAndVectors()
        Matrix_<float> af(4,4); for (int i=0; i<4; ++i) for (int j=0; j<4; ++j) af(i,j)=(float)a(i,j); 

        for(int i=0;i<4;i++) expectedValuesf[i] = (std::complex<float>)expEigen[i];
        for(int i=0;i<4;i++) for(int j=0;j<4;j++) expectedVectorsf(i,j) = (std::complex<float>)expVectors[j*4+i];

        Eigen  esf(af);   // setup the eigen system

        esf.getAllEigenValuesAndVectors( valuesf, vectorsf);   // solve for the eigenvalues and vectors of the system

        cout << " float SOLUTION: " << valuesf << "  errnorm=" << absNormComplex(valuesf,expectedValuesf) << endl;
        ASSERT(absNormComplex(valuesf,expectedValuesf) < 0.001);

        cout << "Vectors = " << endl;
        for(int i=0;i<4;i++) {
            computeVecf = vectorsf(i); 
            expectVecf = expectedVectorsf(i);

            errnorm = absNormComplex( computeVecf, expectVecf );
            cout << computeVecf << "  errnorm=" << errnorm << endl;
            ASSERT( errnorm < 0.0001 );
        }
/*
//
//         SYMMETRIC TESTS:
//
        Real S[16] = {     1.0,  2.0,  3.0,  4.0,
                           0.0,  2.0,  3.0,  4.0,
                           0.0,  0.0,  3.0,  4.0,
                           0.0,  0.0,  0.0,  4.0 };

        Real expectSymVals[4] = { -2.0531, -0.5146, -0.2943, 12.8621 };
        Real expectSymVecs[16] = {  -0.7003, -0.5144,  0.2767, -0.4103,
                                    -0.3592,  0.4851, -0.6634, -0.4422,
                                     0.1569,  0.5420,  0.6504, -0.5085,
                                     0.5965, -0.4543, -0.2457, -0.6144 };



        Matrix as(4,4, S);
        as.setMatrixStructure( MatrixStructures::Symmetric ) ;

        Eigen  esym(as);   // setup the eigen system
        Vector symValues;
        Vector expectedSymValues(4);
        for(int i=0;i<4;i++) expectedSymValues[i] = expectSymVals[i];
        Matrix symVectors;
        Vector symVector;
        Vector expectedSymVector; 
        Matrix expectedSymVectors(4,4,expectSymVecs );
        esym.getAllEigenValuesAndVectors( symValues, symVectors );  // solve for the eigenvalues and eigenvectors of the system 
        printf("\n *****  SYMMETRIC:  ***** \n\n"  );

        cout << " Real SOLUTION:  All Values and Vectors" << symValues << "  errnorm=" << absNorm(symValues,expectedSymValues) << endl;

        cout << " Eigen Vectors = " << endl;
        for(int i=0;i<4;i++) {
            symVector = symVectors(i);
            expectedSymVector = expectedSymVectors(i);
            errnorm = absNorm(symVector,expectedSymVector);
            cout << symVectors(i) << "  errnorm=" << errnorm << endl;
//            cout << expectedSymVectors(i) <<  endl;
        }

        Eigen  efewi(as);   // setup the eigen system
        efewi.getFewEigenValuesAndVectors( symValues, symVectors, 2, 3 );  // solve for few eigenvalues and eigenvectors of the system 
        cout << " Real SOLUTION:  Values/Vectors between indices 2 and 3" << symValues << endl;
        for(int j=0;j<symVectors.ncol();j++)  cout << symVectors(j) <<  endl;

        Eigen  efewr(as);   // setup the eigen system

        // solve for few eigenvalues/vectors between -1, 1 
        efewr.getFewEigenValuesAndVectors( symValues, symVectors, -1.0, 1.0 );  
        cout << " Real SOLUTION:  Values/Vectors  between -1, 1 " << symValues << endl;
        for(int j=0;j<symVectors.ncol();j++)  cout << symVectors(j) <<  endl;
*/       

        return 0;
    } 
    catch (std::exception& e) {
        std::printf("FAILED: %s\n", e.what());
        return 1;
    }
}


