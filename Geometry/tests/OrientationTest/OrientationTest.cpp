/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 * Tests for the classes defined in Orientation.h.
 */

#include "simtk/SimTK.h"
#include "simmatrix/Scalar.h"
#include "simmatrix/SmallMatrix.h"
#include "simmatrix/BigMatrix.h"
#include "simbody/internal/Geometry.h"

#include <string>
#include <typeinfo>
#include <iostream>
using std::cout;
using std::endl;

using namespace simtk;

int main() {
try {
    UnitVec3 u;
    cout << "should be NaN: u()=" << u << endl;
    cout << "typename(u)=" << typeid(u).name() << " ~u=" << typeid(~u).name() << endl;

    u = UnitVec3(1,2,3);
    cout << "u=UnitVec(1,2,3): " << u << " transposed: " << ~u << endl;
    cout << "u[2]=" << u[2] << " u(2)=" << u(2) << endl;
    cout << "(~u)[2]=" << (~u)[2] << " (~u)(2)=" << (~u)(2) << endl;

    cout << "u.perp()=" << u.perp() << " ~u*u.perp()=" << ~u*u.perp() << endl;
    cout << "(~u).perp()=" << (~u).perp() << " (~u).perp()*u=" << (~u).perp()*u << endl;

    Vec3 v(u);
    Vec3 w = u;
    cout << "Vec3 v(u)=" << v << "  Vec3 w=u, w=" << w << endl;

    // these lines shouldn't compile
    //u = Vec3(1,2,3);
    //u += UnitVec3(1,2,3);
    //u[2] = 5.;
    //u(2) = 5.;
    //u=0.;
    // end of lines which shouldn't compile


    RotationMat R;
    cout << "should be identity: R()=" << R;
    cout << "typename(R)=" << typeid(R).name() << " ~R=" << typeid(~R).name() << endl;

    TransformMat X;
    cout << "should be identity: X()=" << X;
    cout << "typename(X)=" << typeid(X).name() << " ~X=" << typeid(~X).name() << endl;


    R = RotationMat(u);
    cout << "RotationMat with u as z axis (norm=" << R.norm() << "): " << R; 
    cout << "~R: " << ~R;
    cout << "typename(R[1])=" << typeid(R[1]).name() << " R(1)=" << typeid(R(1)).name() << endl;

    X = TransformMat(R, Vec3(-1, 2, 20));
    cout << "Transform X=" << X;
    cout << "   X.R=" << X.R() << "  X.T=" << X.T() << "  X.TInv=" << X.TInv() << endl;
    cout << "   X.x=" << X.x() << "  y=" << X.y() << "  z=" << X.z() << endl;
    cout << "Transform ~X=" << ~X;
    cout << "   (~X).R=" << (~X).R() << "  (~X).T=" << (~X).T() << "  (~X).TInv=" << (~X).TInv() << endl;
    cout << "   (~X).x=" << (~X).x() << "  y=" << (~X).y() << "  z=" << (~X).z() << endl;    
    cout << "should be identity: X*~X=" << X*~X;
    cout << "should be identity: ~X*X=" << ~X*X;

    return 0;
}
catch(const Exception::Base& e) {
    std::cout << e.getMessage() << std::endl;
}
}
