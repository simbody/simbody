/* -------------------------------------------------------------------------- *
 *                          Simbody(tm): SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
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
 * This is a test program which uses the FactorSVD  class to compute
 * eigen values and eigen vectors
 */

/*
The data for this test is from an example FORTRAN  program from the
Numerical Algorithms Group (NAG)
URL:http://www.nag.com/lapack-ex/lapack-ex.html


Solves for the singular valus and vectors for the
following matrix



A = 2.27   0.28  -0.48   1.07  -2.35   0.62
   -1.54  -1.67  -3.09   1.22   2.93  -7.39
    1.15   0.94   0.99   0.79  -1.45   1.03
   -1.94  -0.78  -0.21   0.63   2.30  -2.57



SOLUTION =
Singular values
     9.9966  3.6831  1.3569  0.5000
 Left singular vectors
          1       2       3       4
 1  -0.1921  0.8030  0.0041 -0.5642
 2   0.8794  0.3926 -0.0752  0.2587
 3  -0.2140  0.2980  0.7827  0.5027
 4   0.3795 -0.3351  0.6178 -0.6017

 Right singular vectors by row (first m rows of V**T)
          1       2       3       4       5       6
 1  -0.2774 -0.2020 -0.2918  0.0938  0.4213 -0.7816
 2   0.6003  0.0301 -0.3348  0.3699 -0.5266 -0.3353
 3  -0.1277  0.2805  0.6453  0.6781  0.0413 -0.1645
 4   0.1323  0.7034  0.1906 -0.5399 -0.0575 -0.3957

 Error estimate for the singular values
        1.1E-15

 Error estimates for the left singular vectors
        1.8E-16    4.8E-16    1.3E-15    1.3E-15

 Error estimates for the right singular vectors
        1.8E-16    4.8E-16    1.3E-15    2.2E-15


*/

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

Real A[24] = {    2.27,   0.28,  -0.48,   1.07,  -2.35,   0.62,
                 -1.54,  -1.67,  -3.09,   1.22,   2.93,  -7.39,
                  1.15,   0.94,   0.99,   0.79,  -1.45,   1.03,
                 -1.94,  -0.78,  -0.21,   0.63,   2.30,  -2.57 };
Real X[4] =  { 9.9966,  3.6831,  1.3569,  0.5000 };


int main () {

    try {
           // Default precision (Real, normally double) test.

        Matrix a(4,6, A);
        Vector singularValues( 4 );
        Vector expectedValues( 4, X );
        Matrix rightVectors;
        Matrix leftVectors;

        FactorSVD  svd(a, 0.01);   // setup the eigen system

        svd.getSingularValues( singularValues );  // solve for the singular values
        cout << " SingularValues rcond = 0.01 : " << singularValues << endl;

        svd.factor(a );   // setup the eigen system
        cout << " SingularValues rcond = default : " << singularValues << "  errnorm=" << (singularValues-expectedValues).norm() << endl;
        svd.getSingularValuesAndVectors( singularValues, leftVectors, rightVectors );  // solve for the singular values
        ASSERT((singularValues-expectedValues).norm() < 0.001);


         printf("Left Vectors = \n");
         for(int i=0;i<leftVectors.ncol();i++) {
             for(int j=0;j<leftVectors.nrow();j++)  printf("%f  ",leftVectors(i,j) );
             printf("\n");
         }


         printf("Right Vectors = \n");
         for(int i=0;i<rightVectors.ncol();i++) {
             for(int j=0;j<rightVectors.nrow();j++)  printf("%f  ",rightVectors(i,j) );
             printf("\n");
         }

       Real C[4] = { 1.0,   2.0,
              1.0,   3.0  };

        Matrix c(2,2, C);
        FactorSVD csvd(c);
        Matrix invSVD;
        csvd.inverse(invSVD);
        cout << " FactorSVD.inverse : " << endl;
        cout << invSVD[0] << endl;
        cout << invSVD[1] << endl;

        Real Z[4] = { 0.0,   0.0,
                     0.0,   0.0  };

        Matrix z(2,2, Z);
        FactorSVD zsvd(z);
        Vector_<double> xz;
        Vector_<double> bz(2);
        bz(1) = bz(0) = 0.0;
        zsvd.solve( bz, xz );
        cout << " solve with mat all zeros : " << endl;
        for(int i=0;i<xz.size();i++) printf("%f ", xz(i) );  printf("\n");

        Matrix_<double> z0;
        FactorSVD z0svd(z0);
        Vector_<double> bz0(0);
        z0svd.solve( bz0, xz );
        cout << " solve with mat(0,0) : " << endl;
        for(int i=0;i<xz.size();i++) printf("%f ", xz(i) );  printf("\n");


        cout << " SVD factorization with mat(0,0) : " << endl;
        FactorSVD z0fsvd(z0);
        z0fsvd.getSingularValuesAndVectors( singularValues, leftVectors, rightVectors );  // solve for the singular values
        cout << " Real SOLUTION: " << singularValues <<  endl;

         printf("Left Vectors = \n");
         for(int i=0;i<leftVectors.ncol();i++) {
             for(int j=0;j<leftVectors.nrow();j++)  printf("%f  ",leftVectors(i,j) );
             printf("\n");
         }


         printf("Right Vectors = \n");
         for(int i=0;i<rightVectors.ncol();i++) {
             for(int j=0;j<rightVectors.nrow();j++)  printf("%f  ",rightVectors(i,j) );
             printf("\n");
         }

        return 0;
    }
    catch (std::exception& e) {
        std::printf("FAILED: %s\n", e.what());
        return 1;
    }
}

