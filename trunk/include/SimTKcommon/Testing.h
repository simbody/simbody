#ifndef SimTK_SimTKCOMMON_TESTING_H_
#define SimTK_SimTKCOMMON_TESTING_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"

#include <cmath>     
#include <ctime>
#include <algorithm> 
#include <iostream>

/** @file
 * This file defines a SimTK::Test class and some related macros which
 * provide functionality useful in regression tests. This header file is
 * <em>not</em> automatically included with SimTKcommon.h; you have to
 * ask for it explicilty. Here's how you use this facility:
 * <pre>
 * \#include "SimTKcommon/Testing.h"
 * void myFirstSubtest() {...}
 * void myNextSubtest() {...}
 * int main() {
 *      SimTK_START_TEST("OverallTestName");
 *          SimTK_SUBTEST(myFirstSubtest);
 *          SimTK_SUBTEST(myNextSubtest);
 *      SimTK_END_TEST();
 * }
 * </pre>
 * This will result in nice output including execution times for
 * the overall test and the individual subtests, and arrange for 
 * any exceptions raised in the tests to be caught, properly reported,
 * and cause a non-zero return from main(). If everything runs
 * successfully, main() will return 0. Here is an example of the
 * output produced:
 * <pre>
 * Starting test TestScalar ...
 *   testIsNaN            ... done. testIsNaN            time: 0s.
 *   testIsInf            ... done. testIsInf            time: 0s.
 *   testIsFinite         ... done. testIsFinite         time: 0s.
 *   testSignBit          ... done. testSignBit          time: 0s.
 *   testSign             ... done. testSign             time: 0s.
 *   testSquareAndCube    ... done. testSquareAndCube    time: 0s.
 * Done. TestScalar time: 0s.
 * </pre>
 * (Admittedly the timings aren't much use in that example!)
 *
 * Within your subtests, several useful macros and static functions
 * are available. By using these macros, the resulting message will 
 * include the actual line number at which the test failure occurred.
 * <pre>
 *      SimTK_TEST(cond)       -- this is like assert(cond)
 *      SimTK_TEST_EQ(a,b)     -- like assert(a==b)
 *      SimTK_TEST_NUMEQ(a,b)  -- like assert(numericallyEqual(a,b))
 * </pre>
 * The NUMEQ macro tests scalar and composite numerical values for
 * equality to within a numerical tolerance, using both relative
 * and absolute tolerances. TODO: currently there is no way to 
 * specify your own tolerance.
 *
 * The SimTK::Test class has a number of static methods that are useful
 * in tests. Currently these are all for generating numerical objects
 * filled with random numbers (all uniform between -1 and 1). These are:
 * <pre>
 *      randReal()      randFloat()     randDouble()
 *      randComplex()   randConjugate()
 *      randVec<M>()    randRow<N>()    randMat<M,N>()  randSymMat<N>()
 *      randVector(m)   randMatrix(m,n)
 *      randVec3()      randMat33()
 *      randSpatialVec() randSpatialMat()
 * </pre>
 * These are invoked Test::randReal() etc.
 */

namespace SimTK {

/// This is the main class to support testing. Objects of this type are
/// created by the SimTK_START_TEST macro; don't allocate them directly.
/// The class name appears directly in tests only for access to its 
/// static members like Test::randMatrix().
class Test {
public:
    class Subtest;
    Test(const std::string& name) : testName(name)
    {
        std::clog << "Starting test " << testName << " ...\n";
        startTime = std::clock(); 
    }
    ~Test() {
        std::clog << "Done. " << testName << " time: " 
                  << (std::clock()-startTime)/CLOCKS_PER_SEC << "s.\n";
    }

    // Scale by the magnitude of the quantities being compared, so that we don't
    // ask for unreasonable precision. For magnitudes near zero, we'll be satisfied
    // if both are very small without demanding that they must also be relatively
    // close. That is, we use a relative tolerance for big numbers and an absolute
    // tolerance for small ones.
    static bool numericallyEqual(float v1, float v2) {
        const float scale = std::max(std::max(std::abs(v1), std::abs(v2)), 0.1f);
        return std::abs(v1-v2) < scale*NTraits<float>::getSignificant();
    }
    static bool numericallyEqual(double v1, double v2) {
        const double scale = std::max(std::max(std::abs(v1), std::abs(v2)), 0.1);
        return std::abs(v1-v2) < scale*NTraits<double>::getSignificant();
    }
    template <class P>
    static bool numericallyEqual(const std::complex<P>& v1, const std::complex<P>& v2) {
        return numericallyEqual(v1.real(), v2.real())
                && numericallyEqual(v1.imag(), v2.imag());
    }
    template <class P>
    static bool numericallyEqual(const conjugate<P>& v1, const conjugate<P>& v2) {
        return numericallyEqual(v1.real(), v2.real())
                && numericallyEqual(v1.imag(), v2.imag());
    }
    template <class P>
    static bool numericallyEqual(const std::complex<P>& v1, const conjugate<P>& v2) {
        return numericallyEqual(v1.real(), v2.real())
                && numericallyEqual(v1.imag(), v2.imag());
    }
    template <class P>
    static bool numericallyEqual(const conjugate<P>& v1, const std::complex<P>& v2) {
        return numericallyEqual(v1.real(), v2.real())
                && numericallyEqual(v1.imag(), v2.imag());
    }
    template <class P>
    static bool numericallyEqual(const negator<P>& v1, const negator<P>& v2) {
        return numericallyEqual(-v1, -v2);  // P, P
    }
    template <class P>
    static bool numericallyEqual(const P& v1, const negator<P>& v2) {
        return numericallyEqual(-v1, -v2);  // P, P
    }
    template <class P>
    static bool numericallyEqual(const negator<P>& v1, const P& v2) {
        return numericallyEqual(-v1, -v2);  // P, P
    }
    template <class P>
    static bool numericallyEqual(const negator<std::complex<P> >& v1, const conjugate<P>& v2) {
        return numericallyEqual(-v1, -v2);  // complex, conjugate
    }
    template <class P>
    static bool numericallyEqual(const negator<conjugate<P> >& v1, const std::complex<P>& v2) {
        return numericallyEqual(-v1, -v2);  // conjugate, complex
    }
    template <class P>
    static bool numericallyEqual(const std::complex<P>& v1, const negator<conjugate<P> >& v2) {
        return numericallyEqual(-v1, -v2); // complex, conjugate
    }
    template <class P>
    static bool numericallyEqual(const conjugate<P>& v1, const negator<std::complex<P> >& v2) {
        return numericallyEqual(-v1, -v2); // conjugate, complex
    }
    template <int M, class E1, int S1, class E2, int S2>
    static bool numericallyEqual(const Vec<M,E1,S1>& v1, const Vec<M,E2,S2>& v2) {
        for (int i=0; i<M; ++i) if (!numericallyEqual(v1[i],v2[i])) return false;
        return true;
    }
    template <int N, class E1, int S1, class E2, int S2>
    static bool numericallyEqual(const Row<N,E1,S1>& v1, const Row<N,E2,S2>& v2) {
        for (int j=0; j<N; ++j) if (!numericallyEqual(v1[j],v2[j])) return false;
        return true;
    }
    template <int M, int N, class E1, int CS1, int RS1, class E2, int CS2, int RS2>
    static bool numericallyEqual(const Mat<N,M,E1,CS1,RS1>& v1, const Mat<N,M,E2,CS2,RS2>& v2) {
        for (int j=0; j<N; ++j) if (!numericallyEqual(v1(j),v2(j))) return false;
        return true;
    }
    template <int N, class E1, int S1, class E2, int S2>
    static bool numericallyEqual(const SymMat<N,E1,S1>& v1, const SymMat<N,E2,S2>& v2) {
        return numericallyEqual(v1.getAsVec(), v2.getAsVec());
    }
    template <class E>
    static bool numericallyEqual(const Vector_<E>& v1, const Vector_<E>& v2) {
        if (v1.size() != v2.size()) return false;
        for (int i=0; i < v1.size(); ++i)
            if (!numericallyEqual(v1[i], v2[i])) return false;
        return true;
    }
    template <class E>
    static bool numericallyEqual(const RowVector_<E>& v1, const RowVector_<E>& v2) {
        if (v1.size() != v2.size()) return false;
        for (int i=0; i < v1.size(); ++i)
            if (!numericallyEqual(v1[i], v2[i])) return false;
        return true;
    }
    template <class E>
    static bool numericallyEqual(const Matrix_<E>& v1, const Matrix_<E>& v2) {
        if (v1.nrow() != v2.nrow() || v1.ncol() != v2.ncol()) return false;
        for (int j=0; j < v1.ncol(); ++j)
            if (!numericallyEqual(v1(j), v2(j))) return false;
        return true;
    }

    // Random numbers
    static Real randReal() {
        static Random::Uniform rand(-1,1);
        return rand.getValue();
    }
    static Complex randComplex() {return Complex(randReal(),randReal());}
    static Conjugate randConjugate() {return Conjugate(randReal(),randReal());}
    static float randFloat() {return (float)randReal();}
    static double randDouble() {return (double)randReal();}

    template <int M> static Vec<M> randVec() 
    {   Vec<M> v; for (int i=0; i<M; ++i) v[i]=randReal(); return v;}
    template <int N> static Row<N> randRow() {return ~randVec<N>();}
    template <int M, int N> static Mat<M,N> randMat()
    {   Mat<M,N> m; for (int j=0; j<N; ++j) m(j)=randVec<M>(); return m;}
    template <int N> static SymMat<N> randSymMat() {SymMat<N> s; s.updAsVec() = randVec<N*(N+1)/2>(); return s;}

    static Vector randVector(int m)
    {   Vector v(m); for (int i=0; i<m; ++i) v[i]=randReal(); return v;}
    static Matrix randMatrix(int m, int n)
    {   Matrix M(m,n); for (int j=0; j<n; ++j) M(j)=randVector(m); return M;}

    static Vec3 randVec3() {return randVec<3>();}
    static Mat33 randMat33() {return randMat<3,3>();}
    static SpatialVec randSpatialVec() {
        return SpatialVec(randVec3(), randVec3());
    }
    static SpatialMat randSpatialMat() {
        return SpatialMat(randMat33(), randMat33(),
                          randMat33(), randMat33());
    }
private:
    std::clock_t startTime;
    std::string  testName;
};

/// Internal utility class for generating test messages for subtests.
class Test::Subtest {
public:
    Subtest(const std::string& name) : subtestName(name)
    {
        char padded[128];
        sprintf(padded, "%-20s", name.c_str());
        paddedName = std::string(padded);
        std::clog << "  " << paddedName << " ... " << std::flush;
        startTime = std::clock(); 
    }
    ~Subtest() {
        std::clog << "done. " << paddedName << " time: " 
                  << (std::clock()-startTime)/CLOCKS_PER_SEC << "s.\n";
    }
private:
    std::clock_t startTime;
    std::string  subtestName;
    std::string  paddedName; // name plus some blanks
};

} // namespace SimTK

/// Invoke this macro before anything else in your test's main().
#define SimTK_START_TEST(testName)      \
    SimTK::Test simtk_test_(testName);  \
    try {

/// Invoke this macro as the last thing in your test's main().
#define SimTK_END_TEST() \
    } catch(const std::exception& e) {                  \
        std::cerr << "Test failed due to exception: "   \
                  << e.what() << std::endl;             \
        return 1;                                       \
    } catch(...) {                                      \
        std::cerr << "Test failed due to unrecognized exception.\n";    \
        return 1;                                       \
    }                                                   \
    return 0;

/// Invoke a subtest in the form of a no-argument function, arranging for some 
/// friendly output and timing information.
#define SimTK_SUBTEST(testFunction) \
    do {SimTK::Test::Subtest sub(#testFunction); (testFunction)();} while(false)
/// Invoke a subtest in the form of a 1-argument function, arranging for some 
/// friendly output and timing information.
#define SimTK_SUBTEST1(testFunction,arg1) \
    do {SimTK::Test::Subtest sub(#testFunction); (testFunction)(arg1);} while(false)
/// Invoke a subtest in the form of a 2-argument function, arranging for some 
/// friendly output and timing information.
#define SimTK_SUBTEST2(testFunction,arg1,arg2) \
    do {SimTK::Test::Subtest sub(#testFunction); (testFunction)(arg1,arg2);} while(false)

/// Test that some condition holds and complain if it doesn't.
#define SimTK_TEST(cond) {SimTK_ASSERT_ALWAYS((cond), "Test condition failed.");}

/// Test that two values are exactly equal, using whatever operator==() is defined
/// for their comparison. If these are numerical values this is an extremely 
/// stringent test; consider using SimTK_TEST_NUMEQ() instead to compare them
/// to a tolerance.
#define SimTK_TEST_EQ(v1,v2)    \
    {SimTK_ASSERT_ALWAYS((v1)==(v2)),   \
     "Test values should have been exactly equal.");}

/// Test that two numerical values are equal to within a reasonable numerical
/// error tolerance, using a relative and absolute error tolerance. In the
/// case of composite types, the test is performed elementwise.
#define SimTK_TEST_NUMEQ(v1,v2)    \
    {SimTK_ASSERT_ALWAYS(SimTK::Test::numericallyEqual((v1),(v2)),   \
     "Test values should have been numerically equivalent.");}

#endif // SimTK_SimTKCOMMON_TESTING_H_
