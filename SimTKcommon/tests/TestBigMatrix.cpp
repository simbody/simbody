/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

#include "SimTKcommon.h"
#include "SimTKcommon/Testing.h"

#include <iostream>

#define ASSERT(cond) {SimTK_ASSERT_ALWAYS(cond, "Assertion failed");}

using std::cout;
using std::endl;
using namespace SimTK;
using namespace std;

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


void testMatDivision() {
    Mat22 m1( 4, 0,
              0, 1);
    Mat22 oom1( .25, 0,
                  0, 1 );
    Mat<2,2,Mat22> m2(Mat22( 2, 0,
                             0, 3 ));
    Mat<2,2,Mat22> oom2(Mat22( .5, 0,
                               0, OneThird));

    SimTK_TEST_EQ(1/m1, oom1); 
    SimTK_TEST_EQ(1/m2, oom2);
}

static void f(const Vec3& v) {
    cout << "f(v)=" << v << endl;
}

void testTransform() {
    Transform X;
    Rotation R;
    Mat33 m;
    Vec3 v;
    Vec<3,Real,6> vs(1,2,3); // funny stride
    Vec<4,Real,9> vs2(1,2,3,0); // funny stride
    f(vs);

    SimTK_TEST(X*vs == -(X*-vs));
    SimTK_TEST(X*vs2 == -(X*-vs2));
    
    SimTK_TEST(R*vs == -(R*-vs));
    SimTK_TEST(~vs*R == -(-~vs*R));
}

// Make sure we can instantiate all of these successfully.
template class MatrixBase<double>;
template class VectorBase<double>;
template class RowVectorBase<double>;
template class MatrixView_<double>;
template class VectorView_<double>;
template class RowVectorView_<double>;
template class Matrix_<double>;
template class Vector_<double>;
template class RowVector_<double>;

template class MatrixBase<negator<double> >;
template class VectorBase<negator<double> >;
template class RowVectorBase<negator<double> >;
template class MatrixView_<negator<double> >;
template class VectorView_<negator<double> >;
template class RowVectorView_<negator<double> >;
template class Matrix_<negator<double> >;
template class Vector_<negator<double> >;
template class RowVector_<negator<double> >;

int main() {
    try {
        // Currently, this only tests a small number of operations that were recently added.
        // It should be expanded into a more comprehensive test of the big matrix classes.

        testMatDivision();
        testTransform();
        
        Matrix m(Mat22(1, 2, 3, 4));
        testMatrix<Matrix,2,2>(m, Mat22(1, 2, 3, 4));
        m += 3;
        testMatrix<Matrix,2,2>(m, Mat22(4, 2, 3, 7));
        m -= 3;
        testMatrix<Matrix,2,2>(m, Mat22(1, 2, 3, 4));
        testMatrix<Matrix,2,2>(m-1, Mat22(0, 2, 3, 3));
        testMatrix<Matrix,2,2>(m+1, Mat22(2, 2, 3, 5));
        testMatrix<Matrix,2,2>(1-m, Mat22(0, -2, -3, -3));
        testMatrix<Matrix,2,2>(1+m, Mat22(2, 2, 3, 5));
        Vector v(Vec3(1, 2, 3));
        testVector(v, Vec3(1, 2, 3));
        v += 2;
        testVector(v, Vec3(3, 4, 5));
        v -= 2;
        testVector(v, Vec3(1, 2, 3));
        testVector(v-1, Vec3(0, 1, 2));
        testVector(v+1, Vec3(2, 3, 4));
        testVector(1-v, Vec3(0, -1, -2));
        testVector(1+v, Vec3(2, 3, 4));
        RowVector r(Row3(1, 2, 3));
        testVector(r, Vec3(1, 2, 3));
        r += 2;
        testVector(r, Vec3(3, 4, 5));
        r -= 2;
        testVector(r, Vec3(1, 2, 3));
        testVector(r-1, Vec3(0, 1, 2));
        testVector(r+1, Vec3(2, 3, 4));
        testVector(1-r, Vec3(0, -1, -2));
        testVector(1+r, Vec3(2, 3, 4));

        Matrix mm( Mat23( 1, 2, 3,
                          7, 8, 9 ) );
        testMatrix<Matrix,2,3>(mm, Mat23(1,2,3,7,8,9));

            // Test copying a column or row of a Matrix into
            // a Vector or RowVector.

        // Test assignment constructor
        Vector vv = mm(1); testVector(vv, Vec2(2,8));
        // Test copy assignment
        vv = mm(0); testVector(vv, Vec2(1,7));
        // Test assignment constructor
        RowVector rr = mm[1]; testVector(rr, Vec3(7,8,9));
        // Test copy assignment
        rr = mm[0]; testVector(rr, Vec3(1,2,3));

            // Test copying a row into a Vector and column into RowVector.

        // Test assignment (copy) constructor
        RowVector rrr = ~mm(1); 
        testVector(rrr, Vec2(2,8));
        // Test copy assignment
        rrr = ~mm(0); testVector(rrr, Vec2(1,7));

        // Test assignment (copy) constructor
        Vector vvv = ~mm[1]; 
        testVector(vvv, Vec3(7,8,9));
        // Test copy assignment
        vvv = ~mm[0]; testVector(vvv, Vec3(1,2,3));

            // Test creating a Matrix that shares space with an Array

        // Easy case: sizeof(element) == sizeof(scalar)
        Array_<Real> rarrmat;
        rarrmat.push_back(1.1); rarrmat.push_back(2.2); // col(0)
        rarrmat.push_back(3.3); rarrmat.push_back(4.4); // col(1)
        Matrix rmatrix(2,2, 2/*lda*/, &rarrmat[0]);
        testMatrix<Matrix,2,2>(rmatrix, Mat22(1.1, 3.3,
                                              2.2, 4.4));

        // Here sizeof(element) != sizeof(scalar)
        Array_<SpatialVec> svarrmat;                                 
        svarrmat.push_back(SpatialVec(Vec3(1,2,3),Vec3(4,5,6)));
        svarrmat.push_back(SpatialVec(Vec3(1.1,2.1,3.1),Vec3(4.1,5.1,6.1)));
        svarrmat.push_back(SpatialVec(Vec3(1.2,2.2,3.2),Vec3(4.2,5.2,6.2)));
        svarrmat.push_back(SpatialVec(Vec3(1.3,2.3,3.3),Vec3(4.3,5.3,6.3)));
        const int szInScalars = sizeof(SpatialVec)/sizeof(Real);
        Matrix_<SpatialVec> svmatrix(2,2, 2*szInScalars/*lda*/, 
                                    (Real*)&svarrmat[0]); 
        Matrix_<SpatialVec> svmatans(2,2);
        svmatans(0,0) = svarrmat[0]; svmatans(1,0)=svarrmat[1];
        svmatans(0,1) = svarrmat[2]; svmatans(1,1)=svarrmat[3];
        SimTK_TEST_EQ_TOL(svmatrix, svmatans, 1e-16); // should be exact

            // Test creating a Vector that shares space with an Array

        // Easy case: sizeof(element) == sizeof(scalar)
        Array_<Real> rarray;
        rarray.push_back(1.1); rarray.push_back(2.2); rarray.push_back(3.3);
        Vector rvector(3, &rarray[0], true);
        testVector(rarray, Vec3(1.1,2.2,3.3));

        // Here sizeof(element) != sizeof(scalar)
        Array_<SpatialVec> svarray;
        svarray.push_back(SpatialVec(Vec3(1,2,3),Vec3(4,5,6)));
        svarray.push_back(SpatialVec(Vec3(1.1,2.1,3.1),Vec3(4.1,5.1,6.1)));
        svarray.push_back(SpatialVec(Vec3(1.2,2.2,3.2),Vec3(4.2,5.2,6.2)));
        Vector_<SpatialVec> svvector(3, (Real*)&svarray[0], true);
        Vector_<SpatialVec> svanswer(3); 
        svanswer[0]=svarray[0];svanswer[1]=svarray[1];svanswer[2]=svarray[2];
        SimTK_TEST_EQ_TOL(svvector, svanswer, 1e-16); // should be exact

        // Create 0-width slices of Matrix that has general shape,
        // vector shape, and row vector shape. This caused trouble before
        // because vector and row shapes use 1d matrix storage; when making
        // a 0-width slice of those they have to go back to general shape.
        // Note that you are allowed to index off the bottom and right if
        // you make a zero-width slice.

        Matrix general(3, 4);
        MatrixView gslice1 = general(1,1,0,2); // middle
        SimTK_TEST(gslice1.nrow()==0 && gslice1.ncol()==2);
        MatrixView gslice2 = general(1,1,1,0); // middle
        SimTK_TEST(gslice2.nrow()==1 && gslice2.ncol()==0);
        MatrixView gslice3 = general(0,0,3,0); // left side
        SimTK_TEST(gslice3.nrow()==3 && gslice3.ncol()==0);
        MatrixView gslice4 = general(0,0,0,4); // top
        SimTK_TEST(gslice4.nrow()==0 && gslice4.ncol()==4);
        MatrixView gslice5 = general(3,0,0,4); // off the bottom
        SimTK_TEST(gslice5.nrow()==0 && gslice5.ncol()==4);
        MatrixView gslice6 = general(0,4,3,0); // off the right side
        SimTK_TEST(gslice6.nrow()==3 && gslice6.ncol()==0);
        MatrixView gslice7 = general(0,0,0,0);
        SimTK_TEST(gslice7.nrow()==0 && gslice7.ncol()==0);
        MatrixView gslice8 = general(1,2,0,0);
        SimTK_TEST(gslice8.nrow()==0 && gslice8.ncol()==0);
        MatrixView gslice9 = general(2,3,0,0);
        SimTK_TEST(gslice9.nrow()==0 && gslice9.ncol()==0);

        MatrixView vector = general(0,1,3,1);
        SimTK_TEST(vector.nrow()==3 && vector.ncol()==1);
        MatrixView vslice1 = vector(0,0,3,0);
        SimTK_TEST(vslice1.nrow()==3 && vslice1.ncol()==0);
        MatrixView vslice2 = vector(0,0,0,0);
        SimTK_TEST(vslice2.nrow()==0 && vslice2.ncol()==0);
        MatrixView vslice3 = vector(2,0,1,0);
        SimTK_TEST(vslice3.nrow()==1 && vslice3.ncol()==0);
        MatrixView vslice4 = vector(3,0,0,1); // off the bottom
        SimTK_TEST(vslice4.nrow()==0 && vslice4.ncol()==1);
        MatrixView vslice5 = vector(0,1,3,0); // off the right
        SimTK_TEST(vslice5.nrow()==3 && vslice5.ncol()==0);
        vslice5 = Matrix(3,0);


    } catch(const std::exception& e) {
        cout << "exception: " << e.what() << endl;
        return 1;
    }
    cout << "Done" << endl;
    return 0;
}
