/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKcommon.h"

#include <iostream>

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using std::cout;
using std::endl;
using namespace SimTK;
using namespace std;

static Real sqrt(int a) {return std::sqrt(Real(a));}
static Real sin(int a) {return std::sin(Real(a));}
static Real cos(int a) {return std::cos(Real(a));}
static Real tan(int a) {return std::tan(Real(a));}
static Real asin(int a) {return std::asin(Real(a));}
static Real acos(int a) {return std::acos(Real(a));}
static Real atan(int a) {return std::atan(Real(a));}
static Real exp(int a) {return std::exp(Real(a));}
static Real log(int a) {return std::log(Real(a));}
static Real cosh(int a) {return std::cosh(Real(a));}
static Real sinh(int a) {return std::sinh(Real(a));}
static Real tanh(int a) {return std::tanh(Real(a));}

template <class T, int N>
void testVector(const T& value, const Vec<N>& expected) {
    ASSERT(value.size() == N);
    for (int i = 0; i < N; ++i) {
        if (isNaN(expected[i])) {
            ASSERT(isNaN(value[i]));
        }
        else {
            ASSERT(value[i] == expected[i]);
        }
    }
}

template <class T, int M, int N>
void testMatrix(const T& value, const Mat<M, N>& expected) {
    ASSERT(value.nrow() == M);
    ASSERT(value.ncol() == N);
    for (int i = 0; i < M; ++i)
        for (int j = 0; j < N; ++j) {
            if (isNaN(expected(i, j))) {
                ASSERT(isNaN(value(i, j)));
            }
            else {
                ASSERT(value(i, j) == expected(i, j));
            }
        }
}

template <class T, int N>
void testSymMat(const T& value, const SymMat<N>& expected) {
    ASSERT(value.nrow() == N);
    for (int i = 0; i < N; ++i)
        for (int j = 0; j <= i; ++j) {
            if (isNaN(expected(i, j))) {
                ASSERT(isNaN(value(i, j)));
            }
            else {
                ASSERT(value(i, j) == expected(i, j));
            }
        }
}

int main() {
    try {
        // Create a bunch of vectors and matrices that will be used for testing.
        
        Vector vector(5);
        vector[0] = -1;
        vector[1] = 2;
        vector[2] = -3;
        vector[3] = 4;
        vector[4] = -5;
        RowVector rowvector = ~vector;
        Vec5 vec = Vec5(-1, 2, -3, 4, -5);
        Row5 row = Row5(-1, 2, -3, 4, -5);
        Matrix matrix(2, 3);
        matrix[0][0] = -1;
        matrix[0][1] = 2;
        matrix[0][2] = -3;
        matrix[1][0] = 4;
        matrix[1][1] = -5;
        matrix[1][2] = 6;
        Mat<2,3> mat = Mat<2,3>(-1,  2, -3, 
                                 4, -5,  6);
        Real temp[] = {-1, 2, -3};
        SymMat<2> symmat = SymMat<2>(temp);
        
        // Test the abs function.
        
        Vec5 expectedVec = Vec5(1, 2, 3, 4, 5);
        Mat<2,3> expectedMat = Mat<2,3>(1, 2, 3, 
                                        4, 5, 6);
        testVector(abs(vector), expectedVec);
        testVector(abs(rowvector), expectedVec);
        testVector(abs(vec), expectedVec);
        testVector(abs(row), expectedVec);
        testMatrix<Matrix,2,3>(abs(matrix), expectedMat);
        testMatrix<Mat<2,3>,2,3>(abs(mat), expectedMat);
        testSymMat(abs(symmat), SymMat<2>(Mat22(1, 2, 
                                                2, 3)));
        
        // Test the exp function.
        
        expectedVec = Vec5(exp(-1), exp(2), exp(-3), exp(4), exp(-5));
        expectedMat = Mat<2,3>(exp(-1), exp( 2), exp(-3), 
                               exp( 4), exp(-5), exp( 6));
        testVector(exp(vector), expectedVec);
        testVector(exp(rowvector), expectedVec);
        testVector(exp(vec), expectedVec);
        testVector(exp(row), expectedVec);
        testMatrix<Matrix,2,3>(exp(matrix), expectedMat);
        testMatrix<Mat<2,3>,2,3>(exp(mat), expectedMat);
        testSymMat(exp(symmat), SymMat<2>(Mat22(exp(-1), exp( 2), 
                                                exp( 2), exp(-3))));
        
        // Test the log function.
        
        expectedVec = Vec5(log(-1), log(2), log(-3), log(4), log(-5));
        expectedMat = Mat<2,3>(log(-1), log(2), log(-3), log(4), log(-5), log(6));
        testVector(log(vector), expectedVec);
        testVector(log(rowvector), expectedVec);
        testVector(log(vec), expectedVec);
        testVector(log(row), expectedVec);
        testMatrix<Matrix,2,3>(log(matrix), expectedMat);
        testMatrix<Mat<2,3>,2,3>(log(mat), expectedMat);
        testSymMat(log(symmat), SymMat<2>(Mat22(log(-1), log(2), log(2), log(-3))));
        
        // Test the sqrt function.
        
        expectedVec = Vec5(sqrt(-1), sqrt(2), sqrt(-3), sqrt(4), sqrt(-5));
        expectedMat = Mat<2,3>(sqrt(-1), sqrt(2), sqrt(-3), sqrt(4), sqrt(-5), sqrt(6));
        testVector(sqrt(vector), expectedVec);
        testVector(sqrt(rowvector), expectedVec);
        testVector(sqrt(vec), expectedVec);
        testVector(sqrt(row), expectedVec);
        testMatrix<Matrix,2,3>(sqrt(matrix), expectedMat);
        testMatrix<Mat<2,3>,2,3>(sqrt(mat), expectedMat);
        testSymMat(sqrt(symmat), SymMat<2>(Mat22(sqrt(-1), sqrt(2), sqrt(2), sqrt(-3))));
        
        // Test the sin function.
        
        expectedVec = Vec5(sin(-1), sin(2), sin(-3), sin(4), sin(-5));
        expectedMat = Mat<2,3>(sin(-1), sin(2), sin(-3), sin(4), sin(-5), sin(6));
        testVector(sin(vector), expectedVec);
        testVector(sin(rowvector), expectedVec);
        testVector(sin(vec), expectedVec);
        testVector(sin(row), expectedVec);
        testMatrix<Matrix,2,3>(sin(matrix), expectedMat);
        testMatrix<Mat<2,3>,2,3>(sin(mat), expectedMat);
        testSymMat(sin(symmat), SymMat<2>(Mat22(sin(-1), sin(2), sin(2), sin(-3))));
        
        // Test the cos function.
        
        expectedVec = Vec5(cos(-1), cos(2), cos(-3), cos(4), cos(-5));
        expectedMat = Mat<2,3>(cos(-1), cos(2), cos(-3), cos(4), cos(-5), cos(6));
        testVector(cos(vector), expectedVec);
        testVector(cos(rowvector), expectedVec);
        testVector(cos(vec), expectedVec);
        testVector(cos(row), expectedVec);
        testMatrix<Matrix,2,3>(cos(matrix), expectedMat);
        testMatrix<Mat<2,3>,2,3>(cos(mat), expectedMat);
        testSymMat(cos(symmat), SymMat<2>(Mat22(cos(-1), cos(2), cos(2), cos(-3))));
        
        // Test the tan function.
        
        expectedVec = Vec5(tan(-1), tan(2), tan(-3), tan(4), tan(-5));
        expectedMat = Mat<2,3>(tan(-1), tan(2), tan(-3), tan(4), tan(-5), tan(6));
        testVector(tan(vector), expectedVec);
        testVector(tan(rowvector), expectedVec);
        testVector(tan(vec), expectedVec);
        testVector(tan(row), expectedVec);
        testMatrix<Matrix,2,3>(tan(matrix), expectedMat);
        testMatrix<Mat<2,3>,2,3>(tan(mat), expectedMat);
        testSymMat(tan(symmat), SymMat<2>(Mat22(tan(-1), tan(2), tan(2), tan(-3))));
        
        // Test the asin function.
        
        expectedVec = Vec5(asin(-1), asin(2), asin(-3), asin(4), asin(-5));
        expectedMat = Mat<2,3>(asin(-1), asin(2), asin(-3), asin(4), asin(-5), asin(6));
        testVector(asin(vector), expectedVec);
        testVector(asin(rowvector), expectedVec);
        testVector(asin(vec), expectedVec);
        testVector(asin(row), expectedVec);
        testMatrix<Matrix,2,3>(asin(matrix), expectedMat);
        testMatrix<Mat<2,3>,2,3>(asin(mat), expectedMat);
        testSymMat(asin(symmat), SymMat<2>(Mat22(asin(-1), asin(2), asin(2), asin(-3))));
        
        // Test the asin function.
        
        expectedVec = Vec5(acos(-1), acos(2), acos(-3), acos(4), acos(-5));
        expectedMat = Mat<2,3>(acos(-1), acos(2), acos(-3), acos(4), acos(-5), acos(6));
        testVector(acos(vector), expectedVec);
        testVector(acos(rowvector), expectedVec);
        testVector(acos(vec), expectedVec);
        testVector(acos(row), expectedVec);
        testMatrix<Matrix,2,3>(acos(matrix), expectedMat);
        testMatrix<Mat<2,3>,2,3>(acos(mat), expectedMat);
        testSymMat(acos(symmat), SymMat<2>(Mat22(acos(-1), acos(2), acos(2), acos(-3))));
        
        // Test the atan function.
        
        expectedVec = Vec5(atan(-1), atan(2), atan(-3), atan(4), atan(-5));
        expectedMat = Mat<2,3>(atan(-1), atan(2), atan(-3), atan(4), atan(-5), atan(6));
        testVector(atan(vector), expectedVec);
        testVector(atan(rowvector), expectedVec);
        testVector(atan(vec), expectedVec);
        testVector(atan(row), expectedVec);
        testMatrix<Matrix,2,3>(atan(matrix), expectedMat);
        testMatrix<Mat<2,3>,2,3>(atan(mat), expectedMat);
        testSymMat(atan(symmat), SymMat<2>(Mat22(atan(-1), atan(2), atan(2), atan(-3))));
        
        // Test the sinh function.
        
        expectedVec = Vec5(sinh(-1), sinh(2), sinh(-3), sinh(4), sinh(-5));
        expectedMat = Mat<2,3>(sinh(-1), sinh(2), sinh(-3), sinh(4), sinh(-5), sinh(6));
        testVector(sinh(vector), expectedVec);
        testVector(sinh(rowvector), expectedVec);
        testVector(sinh(vec), expectedVec);
        testVector(sinh(row), expectedVec);
        testMatrix<Matrix,2,3>(sinh(matrix), expectedMat);
        testMatrix<Mat<2,3>,2,3>(sinh(mat), expectedMat);
        testSymMat(sinh(symmat), SymMat<2>(Mat22(sinh(-1), sinh(2), sinh(2), sinh(-3))));
        
        // Test the cosh function.
        
        expectedVec = Vec5(cosh(-1), cosh(2), cosh(-3), cosh(4), cosh(-5));
        expectedMat = Mat<2,3>(cosh(-1), cosh(2), cosh(-3), cosh(4), cosh(-5), cosh(6));
        testVector(cosh(vector), expectedVec);
        testVector(cosh(rowvector), expectedVec);
        testVector(cosh(vec), expectedVec);
        testVector(cosh(row), expectedVec);
        testMatrix<Matrix,2,3>(cosh(matrix), expectedMat);
        testMatrix<Mat<2,3>,2,3>(cosh(mat), expectedMat);
        testSymMat(cosh(symmat), SymMat<2>(Mat22(cosh(-1), cosh(2), cosh(2), cosh(-3))));
        
        // Test the tanh function.
        
        expectedVec = Vec5(tanh(-1), tanh(2), tanh(-3), tanh(4), tanh(-5));
        expectedMat = Mat<2,3>(tanh(-1), tanh(2), tanh(-3), tanh(4), tanh(-5), tanh(6));
        testVector(tanh(vector), expectedVec);
        testVector(tanh(rowvector), expectedVec);
        testVector(tanh(vec), expectedVec);
        testVector(tanh(row), expectedVec);
        testMatrix<Matrix,2,3>(tanh(matrix), expectedMat);
        testMatrix<Mat<2,3>,2,3>(tanh(mat), expectedMat);
        testSymMat(tanh(symmat), SymMat<2>(Mat22(tanh(-1), tanh(2), tanh(2), tanh(-3))));
        
        // Test the sum function.
        
        ASSERT(sum(vector) == -3);
        ASSERT(sum(rowvector) == -3);
        ASSERT(sum(vec) == -3);
        ASSERT(sum(row) == -3);
        testVector(sum(matrix), Vec3(3, -3, 3));
        testVector(sum(mat), Vec3(3, -3, 3));
        testVector(sum(symmat), Vec2(1, -1));
        
        // Test the min function.
        
        ASSERT(min(vector) == -5);
        ASSERT(min(rowvector) == -5);
        ASSERT(min(vec) == -5);
        ASSERT(min(row) == -5);
        testVector(min(matrix), Vec3(-1, -5, -3));
        testVector(min(mat), Vec3(-1, -5, -3));
        testVector(min(symmat), Vec2(-1, -3));
        
        // Test the max function.
        
        ASSERT(max(vector) == 4);
        ASSERT(max(rowvector) == 4);
        ASSERT(max(vec) == 4);
        ASSERT(max(row) == 4);
        testVector(max(matrix), Vec3(4, 2, 6));
        testVector(max(mat), Vec3(4, 2, 6));
        testVector(max(symmat), Vec2(2, 2));
        
        // Test the mean function.
        
        ASSERT(mean(vector) == -0.6);
        ASSERT(mean(rowvector) == -0.6);
        ASSERT(mean(vec) == -0.6);
        ASSERT(mean(row) == -0.6);
        testVector(mean(matrix), Vec3(1.5, -1.5, 1.5));
        testVector(mean(mat), Vec3(1.5, -1.5, 1.5));
        testVector(mean(symmat), Vec2(0.5, -0.5));
        
        // Test the sort function.

        expectedVec = Vec5(-5, -3, -1, 2, 4);
        expectedMat = Mat<2,3>(-1, -5, -3, 4, 2, 6);
        testVector(sort(vector), expectedVec);
        testVector(sort(rowvector), expectedVec);
        testVector(sort(vec), expectedVec);
        testVector(sort(row), expectedVec);
        testMatrix<Matrix,2,3>(sort(matrix), expectedMat);
        testMatrix<Mat<2,3>,2,3>(sort(mat), expectedMat);
        testMatrix<Mat<2,2>,2,2>(sort(symmat), Mat22(-1, -3, 
                                                      2,  2));
        
        // Test the median function.
        
        ASSERT(median(vector) == -1);
        ASSERT(median(rowvector) == -1);
        ASSERT(median(vec) == -1);
        ASSERT(median(row) == -1);
        testVector(median(matrix), Vec3(1.5, -1.5, 1.5));
        testVector(median(mat), Vec3(1.5, -1.5, 1.5));
        testVector(median(symmat), Vec2(0.5, -0.5));
} catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}

