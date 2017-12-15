#ifndef SimTK_SimTKCOMMON_TESTING_H_
#define SimTK_SimTKCOMMON_TESTING_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2009-15 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/Random.h"
#include "SimTKcommon/internal/Timing.h"

#include <cmath>     
#include <algorithm> 
#include <iostream>

/** @file
 * This file defines a SimTK::Test class and some related macros which
 * provide functionality useful in regression tests.
 */

namespace SimTK {

/**@defgroup RegressionTesting    SimTK Regression Testing
 *
 * SimTK defines some utilities to facilitate the creation of regression
 * tests for SimTK facilities. These utilities consist of a SimTK::Test
 * class and related support macros.
 *
 * Features include:
 *      - uniform, readable output including execution times
 *      - identical comparison tests for all numerical types, scalar and composite
 *      - careful treatment of numerical tolerances using relative and absolute
 *          comparisons, with provision for size-dependent, reduced accuracy
 *          expectations for matrix operations
 *      - default tolerance varies with precision (caller can override)
 *      - convenient generation of random test data
 *      - convenient testing of required argument checking (i.e., test
 *          fails unless an exception is thrown)
 *
 * Here's how you use this facility:
 * <pre>
 * \#include "SimTKcommon.h"
 * void myFirstSubtest() {...}
 * void myNextSubtest() {...}
 * int main() {
 *      SimTK_START_TEST("OverallTestName");
 *          SimTK_SUBTEST(myFirstSubtest);
 *          SimTK_SUBTEST(myNextSubtest);
 *      SimTK_END_TEST();
 * }
 * </pre>
 * The arguments to SimTK_SUBTEST are function names and will be
 * called with "()" appended. If your subtest functions have arguments,
 * use SimTK_SUBTEST1(name,arg) or SimTK_SUBTEST2(name,arg1,arg2) which
 * will call name(arg) or name(arg1,arg2) as appropriate.
 *
 * This will result in nice output including execution times for
 * the overall test and the individual subtests, and arrange for 
 * any exceptions raised in the tests to be caught, properly reported,
 * and cause a non-zero return from main(). If everything runs
 * successfully, main() will return 0. Here is an example of the
 * output produced:
 * <pre>
 * Starting test TestScalar ...
 *   testIsNaN            ... done. testIsNaN            time: 0ms.
 *   testIsInf            ... done. testIsInf            time: 0ms.
 *   testIsFinite         ... done. testIsFinite         time: 0ms.
 *   testSignBit          ... done. testSignBit          time: 0ms.
 *   testSign             ... done. testSign             time: 0ms.
 *   testSquareAndCube    ... done. testSquareAndCube    time: 0ms.
 * Done. TestScalar time: 15ms.
 * </pre>
 * (Admittedly the timings aren't much use in that example!)
 *
 * Within your subtests, several useful macros and static functions
 * are available. By using these macros, the resulting message will 
 * include the actual line number at which the test failure occurred.
 * <pre>
 *      SimTK_TEST(cond)             -- this is like assert(cond)
 *      SimTK_TEST_FAILED("message") -- like assert(!"message")
 *
 *      SimTK_TEST_EQ(a,b)           -- equal to within a default tolerance
 *      SimTK_TEST_NOTEQ(a,b)        -- not equal to within a default tolerance
 *
 *      SimTK_TEST_EQ_SIZE(a,b,n)    -- equal to within n * default tolerance
 *      SimTK_TEST_NOTEQ_SIZE(a,b,n) -- not equal to within n * default tolerance
 *
 *      SimTK_TEST_EQ_TOL(a,b,tol)   -- same as above with specified tolerance
 *      SimTK_TEST_NOTEQ_TOL(a,b,tol)
 *
 *      SimTK_TEST_MUST_THROW(statement)        -- we expect the statement to throw some exception
 *      SimTK_TEST_MUST_THROW_EXC(statement, exception) -- we expect a particular exception type
 *      SimTK_TEST_MUST_THROW_DEBUG(statement)  -- same as above but only checked in Debug builds
 *      SimTK_TEST_MUST_THROW_EXC_DEBUG(statement, exception) -- ditto
 *
 *      \#define SimTK_TEST_SUPPRESS_MUST_THROW 
 *          -- Define this temporarily at top of test programs to make it 
 *             easier to locate an *unexpected* exception while debugging. It
 *             simply disables the "MUST_THROW" macros so you don't have to wade
 *             through them in the debugger to get to the actual problem.
 * </pre>
 * The SimTK_TEST_EQ macros test scalar and composite numerical values for
 * equality to within a numerical tolerance, using both relative
 * and absolute tolerances. The default is the value of SignificantReal
 * for the underlying numerical type. For composite types the equality test is done
 * elementwise; that is, we apply it strictly to each pair of elements not
 * to an overall norm.
 *
 * The SimTK_TEST_EQ_SIZE macros allows you to specify a multiple of default
 * tolerance to be used. This is necessary for most Matrix operations since 
 * attainable accuracy falls off with the size of the matrix. Typically, if 
 * the smallest dimension of the Matrix is n, then the tolerance you should allow
 * is n*scalarTol where scalarTol is the default tolerance for a scalar
 * operation. Note that you still need to specify size when comparing 
 * Vector or scalar values if those values were produced using a matrix
 * computation.
 *
 * The SimTK_TEST_EQ_TOL macros take a user-specified tolerance value for 
 * the elementwise tests, overriding the default.
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
 *      randRotation()  randTransform()
 * </pre>
 * These are invoked Test::randReal() etc.
 *
 * @{
 */


/// This is the main class to support testing. Objects of this type are
/// created by the SimTK_START_TEST macro; don't allocate them directly.
/// The class name appears directly in tests only for access to its 
/// static members like Test::randMatrix().
class Test {
public:
    class Subtest;
    Test(const std::string& name)
    :   startCpuTime(SimTK::cpuTime()),
        startRealTime(SimTK::realTime()),
        testName(name)
    {
        std::clog << "Starting test " << testName << " ...\n";
    }
    ~Test() {
        const double finalRealTime=SimTK::realTime();
        const double finalCpuTime=SimTK::cpuTime();
        std::ostringstream fmt;
        fmt << std::fixed << std::setprecision(1);
        fmt << "\n" << testName << " done." 
            << " real/CPU ms: " << (finalRealTime-startRealTime)*1000
            << " / "  << (finalCpuTime-startCpuTime)*1000 <<std::endl;
        std::clog << fmt.str();
    }

    template <class T>
    static double defTol() {return (double)NTraits<typename CNT<T>::Precision>::getSignificant();}

    // For dissimilar types, the default tolerance is the narrowest of the two.
    template <class T1, class T2>
    static double defTol2() {return std::max(defTol<T1>(), defTol<T2>());}

    // Scale by the magnitude of the quantities being compared, so that we don't
    // ask for unreasonable precision. For magnitudes near zero, we'll be satisfied
    // if both are very small without demanding that they must also be relatively
    // close. That is, we use a relative tolerance for big numbers and an absolute
    // tolerance for small ones.
    static bool numericallyEqual(float v1, float v2, int n, double tol=defTol<float>()) {
        const float scale = n*std::max(std::max(std::abs(v1), std::abs(v2)), 1.0f);
        return std::abs(v1-v2) < scale*(float)tol;
    }
    static bool numericallyEqual(double v1, double v2, int n, double tol=defTol<double>()) {
        const double scale = n*std::max(std::max(std::abs(v1), std::abs(v2)), 1.0);
        return std::abs(v1-v2) < scale*(double)tol;
    }

    // For integers we ignore tolerance.
    static bool numericallyEqual(int i1, int i2, int n, double tol=0) {return i1==i2;}
    static bool numericallyEqual(unsigned u1, unsigned u2, int n, double tol=0) {return u1==u2;}

    // Mixed floating types use default tolerance for the narrower type.
    static bool numericallyEqual(float v1, double v2, int n, double tol=defTol<float>())
    {   return numericallyEqual((double)v1, v2, n, tol); }
    static bool numericallyEqual(double v1, float v2, int n, double tol=defTol<float>())
    {   return numericallyEqual(v1, (double)v2, n, tol); }

    // Mixed int/floating just upgrades int to floating type.
    static bool numericallyEqual(int i1, float f2, int n, double tol=defTol<float>())
    {   return numericallyEqual((float)i1,f2,n,tol); }
    static bool numericallyEqual(float f1, int i2, int n, double tol=defTol<float>())
    {   return numericallyEqual(f1,(float)i2,n,tol); }
    static bool numericallyEqual(unsigned i1, float f2, int n, double tol=defTol<float>())
    {   return numericallyEqual((float)i1,f2,n,tol); }
    static bool numericallyEqual(float f1, unsigned i2, int n, double tol=defTol<float>())
    {   return numericallyEqual(f1,(float)i2,n,tol); }
    static bool numericallyEqual(int i1, double f2, int n, double tol=defTol<double>())
    {   return numericallyEqual((double)i1,f2,n,tol); }
    static bool numericallyEqual(double f1, int i2, int n, double tol=defTol<double>())
    {   return numericallyEqual(f1,(double)i2,n,tol); }
    static bool numericallyEqual(unsigned i1, double f2, int n, double tol=defTol<double>())
    {   return numericallyEqual((double)i1,f2,n,tol); }
    static bool numericallyEqual(double f1, unsigned i2, int n, double tol=defTol<double>())
    {   return numericallyEqual(f1,(double)i2,n,tol); }

    template <class P>
    static bool numericallyEqual(const std::complex<P>& v1, const std::complex<P>& v2, int n, double tol=defTol<P>()) {
        return numericallyEqual(v1.real(), v2.real(), n, tol)
            && numericallyEqual(v1.imag(), v2.imag(), n, tol);
    }
    template <class P>
    static bool numericallyEqual(const conjugate<P>& v1, const conjugate<P>& v2, int n, double tol=defTol<P>()) {
        return numericallyEqual(v1.real(), v2.real(), n, tol)
            && numericallyEqual(v1.imag(), v2.imag(), n, tol);
    }
    template <class P>
    static bool numericallyEqual(const std::complex<P>& v1, const conjugate<P>& v2, int n, double tol=defTol<P>()) {
        return numericallyEqual(v1.real(), v2.real(), n, tol)
            && numericallyEqual(v1.imag(), v2.imag(), n, tol);
    }
    template <class P>
    static bool numericallyEqual(const conjugate<P>& v1, const std::complex<P>& v2, int n, double tol=defTol<P>()) {
        return numericallyEqual(v1.real(), v2.real(), n, tol)
            && numericallyEqual(v1.imag(), v2.imag(), n, tol);
    }
    template <class P>
    static bool numericallyEqual(const negator<P>& v1, const negator<P>& v2, int n, double tol=defTol<P>()) {
        return numericallyEqual(-v1, -v2, n, tol);  // P, P
    }
    template <class P>
    static bool numericallyEqual(const P& v1, const negator<P>& v2, int n, double tol=defTol<P>()) {
        return numericallyEqual(-v1, -v2, n, tol);  // P, P
    }
    template <class P>
    static bool numericallyEqual(const negator<P>& v1, const P& v2, int n, double tol=defTol<P>()) {
        return numericallyEqual(-v1, -v2, n, tol);  // P, P
    }
    template <class P>
    static bool numericallyEqual(const negator<std::complex<P> >& v1, const conjugate<P>& v2, int n, double tol=defTol<P>()) {
        return numericallyEqual(-v1, -v2, n, tol);  // complex, conjugate
    }
    template <class P>
    static bool numericallyEqual(const negator<conjugate<P> >& v1, const std::complex<P>& v2, int n, double tol=defTol<P>()) {
        return numericallyEqual(-v1, -v2, n, tol);  // conjugate, complex
    }
    template <class P>
    static bool numericallyEqual(const std::complex<P>& v1, const negator<conjugate<P> >& v2, int n, double tol=defTol<P>()) {
        return numericallyEqual(-v1, -v2, n, tol); // complex, conjugate
    }
    template <class P>
    static bool numericallyEqual(const conjugate<P>& v1, const negator<std::complex<P> >& v2, int n, double tol=defTol<P>()) {
        return numericallyEqual(-v1, -v2, n, tol); // conjugate, complex
    }
    template <int M, class E1, int S1, class E2, int S2>
    static bool numericallyEqual(const Vec<M,E1,S1>& v1, const Vec<M,E2,S2>& v2, int n, double tol=(defTol2<E1,E2>())) {
        for (int i=0; i<M; ++i) if (!numericallyEqual(v1[i],v2[i], n, tol)) return false;
        return true;
    }
    template <int N, class E1, int S1, class E2, int S2>
    static bool numericallyEqual(const Row<N,E1,S1>& v1, const Row<N,E2,S2>& v2, int n, double tol=(defTol2<E1,E2>())) {
        for (int j=0; j<N; ++j) if (!numericallyEqual(v1[j],v2[j], n, tol)) return false;
        return true;
    }
    template <int M, int N, class E1, int CS1, int RS1, class E2, int CS2, int RS2>
    static bool numericallyEqual(const Mat<M,N,E1,CS1,RS1>& v1, const Mat<M,N,E2,CS2,RS2>& v2, int n, double tol=(defTol2<E1,E2>())) {
        for (int j=0; j<N; ++j) if (!numericallyEqual(v1(j),v2(j), n, tol)) return false;
        return true;
    }
    template <int N, class E1, int S1, class E2, int S2>
    static bool numericallyEqual(const SymMat<N,E1,S1>& v1, const SymMat<N,E2,S2>& v2, int n, double tol=(defTol2<E1,E2>())) {
        return numericallyEqual(v1.getAsVec(), v2.getAsVec(), n, tol);
    }
    template <class E1, class E2>
    static bool numericallyEqual(const VectorView_<E1>& v1, const VectorView_<E2>& v2, int n, double tol=(defTol2<E1,E2>())) {
        if (v1.size() != v2.size()) return false;
        for (int i=0; i < v1.size(); ++i)
            if (!numericallyEqual(v1[i], v2[i], n, tol)) return false;
        return true;
    }
    template <class E1, class E2>
    static bool numericallyEqual(const Vector_<E1>& v1, const Vector_<E2>& v2, int n, double tol=(defTol2<E1,E2>()))
    {   return numericallyEqual((const VectorView_<E1>&)v1, (const VectorView_<E2>&)v2, n, tol); }
    template <class E1, class E2>
    static bool numericallyEqual(const Vector_<E1>& v1, const VectorView_<E2>& v2, int n, double tol=(defTol2<E1,E2>()))
    {   return numericallyEqual((const VectorView_<E1>&)v1, (const VectorView_<E2>&)v2, n, tol); }
    template <class E1, class E2>
    static bool numericallyEqual(const VectorView_<E1>& v1, const Vector_<E2>& v2, int n, double tol=(defTol2<E1,E2>()))
    {   return numericallyEqual((const VectorView_<E1>&)v1, (const VectorView_<E2>&)v2, n, tol); }

    template <class E1, class E2>
    static bool numericallyEqual(const RowVectorView_<E1>& v1, const RowVectorView_<E2>& v2, int n, double tol=(defTol2<E1,E2>())) {
        if (v1.size() != v2.size()) return false;
        for (int i=0; i < v1.size(); ++i)
            if (!numericallyEqual(v1[i], v2[i], n, tol)) return false;
        return true;
    }
    template <class E1, class E2>
    static bool numericallyEqual(const RowVector_<E1>& v1, const RowVector_<E2>& v2, int n, double tol=(defTol2<E1,E2>()))
    {   return numericallyEqual((const RowVectorView_<E1>&)v1, (const RowVectorView_<E2>&)v2, n, tol); }
    template <class E1, class E2>
    static bool numericallyEqual(const RowVector_<E1>& v1, const RowVectorView_<E2>& v2, int n, double tol=(defTol2<E1,E2>()))
    {   return numericallyEqual((const RowVectorView_<E1>&)v1, (const RowVectorView_<E2>&)v2, n, tol); }
    template <class E1, class E2>
    static bool numericallyEqual(const RowVectorView_<E1>& v1, const RowVector_<E2>& v2, int n, double tol=(defTol2<E1,E2>()))
    {   return numericallyEqual((const RowVectorView_<E1>&)v1, (const RowVectorView_<E2>&)v2, n, tol); }

    template <class E1, class E2>
    static bool numericallyEqual(const MatrixView_<E1>& v1, const MatrixView_<E2>& v2, int n, double tol=(defTol2<E1,E2>())) {
        if (v1.nrow() != v2.nrow() || v1.ncol() != v2.ncol()) return false;
        for (int j=0; j < v1.ncol(); ++j)
            if (!numericallyEqual(v1(j), v2(j), n, tol)) return false;
        return true;
    }
    template <class E1, class E2>
    static bool numericallyEqual(const Matrix_<E1>& m1, const Matrix_<E2>& m2, int n, double tol=(defTol2<E1,E2>()))
    {   return numericallyEqual((const MatrixView_<E1>&)m1, (const MatrixView_<E2>&)m2, n, tol); }
    template <class E1, class E2>
    static bool numericallyEqual(const Matrix_<E1>& m1, const MatrixView_<E2>& m2, int n, double tol=(defTol2<E1,E2>()))
    {   return numericallyEqual((const MatrixView_<E1>&)m1, (const MatrixView_<E2>&)m2, n, tol); }
    template <class E1, class E2>
    static bool numericallyEqual(const MatrixView_<E1>& m1, const Matrix_<E2>& m2, int n, double tol=(defTol2<E1,E2>()))
    {   return numericallyEqual((const MatrixView_<E1>&)m1, (const MatrixView_<E2>&)m2, n, tol); }

    template <class P>
    static bool numericallyEqual(const Rotation_<P>& R1, const Rotation_<P>& R2, int n, double tol=defTol<P>()) {
        return R1.isSameRotationToWithinAngle(R2, (Real)(n*tol));
    }

    template <class P>
    static bool numericallyEqual(const Transform_<P>& T1, const Transform_<P>& T2, int n, double tol=defTol<P>()) {
        return numericallyEqual(T1.R(), T2.R(), n, tol)
            && numericallyEqual(T1.p(), T2.p(), n, tol);
    }

    template <class P>
    static bool numericallyEqual(const UnitInertia_<P>& G1, const UnitInertia_<P>& G2, int n, double tol=defTol<P>()) {
        return numericallyEqual(G1.asSymMat33(),G2.asSymMat33(), n, tol);
    }

    template <class P>
    static bool numericallyEqual(const Inertia_<P>& I1, const Inertia_<P>& I2, int n, double tol=defTol<P>()) {
        return numericallyEqual(I1.asSymMat33(),I2.asSymMat33(), n, tol);
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
    template <int N> static SymMat<N> randSymMat() 
    {   SymMat<N> s; s.updAsVec() = randVec<N*(N+1)/2>(); return s; }

    static Vector randVector(int m)
    {   Vector v(m); for (int i=0; i<m; ++i) v[i]=randReal(); return v;}
    static Matrix randMatrix(int m, int n)
    {   Matrix M(m,n); for (int j=0; j<n; ++j) M(j)=randVector(m); return M;}

    static Vec3 randVec3() {return randVec<3>();}
    static Mat33 randMat33() {return randMat<3,3>();}
    static SymMat33 randSymMat33() {return randSymMat<3>();}
    static SpatialVec randSpatialVec() {
        return SpatialVec(randVec3(), randVec3());
    }
    static SpatialMat randSpatialMat() {
        return SpatialMat(randMat33(), randMat33(),
                          randMat33(), randMat33());
    }
    static Rotation randRotation() {
        // Generate random angle and random axis to rotate around.
        return Rotation((Pi/2)*randReal(), randVec3());
    }
    static Transform randTransform() {
        return Transform(randRotation(), randVec3());
    }
private:
    const double startCpuTime;
    const double startRealTime;
    std::string  testName;
};

/// Internal utility class for generating test messages for subtests.
class Test::Subtest {
public:
    Subtest(const std::string& name) 
    :   startCpuTime(SimTK::cpuTime()),
        startRealTime(SimTK::realTime()),
        subtestName(name)
    {
        std::clog << "  " << subtestName << " ...\n" << std::flush;
    }
    ~Subtest() {
        const double finalRealTime=SimTK::realTime();
        const double finalCpuTime=SimTK::cpuTime();
        std::ostringstream fmt;
        fmt << std::fixed << std::setprecision(1);
        fmt << "  " << subtestName << " done."
            << " real/CPU ms: " << (finalRealTime-startRealTime)*1000
            << " / "  << (finalCpuTime-startCpuTime)*1000 <<std::endl;
        std::clog << fmt.str();
    }
private:
    const double startCpuTime;
    const double startRealTime;
    std::string  subtestName;
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
/// Invoke a subtest in the form of a 3-argument function, arranging for some 
/// friendly output and timing information.
#define SimTK_SUBTEST3(testFunction,arg1,arg2,arg3) \
    do {SimTK::Test::Subtest sub(#testFunction); (testFunction)(arg1,arg2,arg3);} while(false)
/// Invoke a subtest in the form of a 4-argument function, arranging for some 
/// friendly output and timing information.
#define SimTK_SUBTEST4(testFunction,arg1,arg2,arg3,arg4) \
    do {SimTK::Test::Subtest sub(#testFunction); (testFunction)(arg1,arg2,arg3,arg4);} while(false)

/// Test that some condition holds and complain if it doesn't.
#define SimTK_TEST(cond) {SimTK_ASSERT_ALWAYS((cond), "Test condition failed.");}

/// Call this if you have determined that a test case has failed and just need
/// to report it and die. Pass the message as a string in quotes.
#define SimTK_TEST_FAILED(msg) {SimTK_ASSERT_ALWAYS(!"Test case failed.", msg);}

/// Call this if you have determined that a test case has failed and just need
/// to report it and die. The message is a printf format string in quotes; here
/// with one argument expected.
#define SimTK_TEST_FAILED1(fmt,a1) {SimTK_ASSERT1_ALWAYS(!"Test case failed.",fmt,a1);}

/// Call this if you have determined that a test case has failed and just need
/// to report it and die. The message is a printf format string in quotes; here
/// with two arguments expected.
#define SimTK_TEST_FAILED2(fmt,a1,a2) {SimTK_ASSERT2_ALWAYS(!"Test case failed.",fmt,a1,a2);}

/// Test that two numerical values are equal to within a reasonable numerical
/// error tolerance, using a relative and absolute error tolerance. In the
/// case of composite types, the test is performed elementwise.
#define SimTK_TEST_EQ(v1,v2)    \
    {SimTK_ASSERT_ALWAYS(SimTK::Test::numericallyEqual((v1),(v2),1),   \
     "Test values should have been numerically equivalent at default tolerance.");}

/// Test that two numerical values are equal to within a specified multiple of the
/// default error tolerance.
#define SimTK_TEST_EQ_SIZE(v1,v2,n)    \
    {SimTK_ASSERT1_ALWAYS(SimTK::Test::numericallyEqual((v1),(v2),(n)),   \
     "Test values should have been numerically equivalent at size=%d times default tolerance.",(n));}

/// Test that two numerical values are equal to within a specified numerical
/// error tolerance, using a relative and absolute error tolerance. In the
/// case of composite types, the test is performed elementwise.
#define SimTK_TEST_EQ_TOL(v1,v2,tol)    \
    {SimTK_ASSERT1_ALWAYS(SimTK::Test::numericallyEqual((v1),(v2),1,(tol)),   \
     "Test values should have been numerically equivalent at tolerance=%g.",(tol));}

/// Test that two numerical values are NOT equal to within a reasonable numerical
/// error tolerance, using a relative and absolute error tolerance. In the
/// case of composite types, the equality test is performed elementwise.
#define SimTK_TEST_NOTEQ(v1,v2)    \
    {SimTK_ASSERT_ALWAYS(!SimTK::Test::numericallyEqual((v1),(v2),1),   \
     "Test values should NOT have been numerically equivalent (at default tolerance).");}

/// Test that two numerical values are NOT equal to within a specified multiple of
/// the default error tolerance, using a relative and absolute error tolerance. In the
/// case of composite types, the equality test is performed elementwise.
#define SimTK_TEST_NOTEQ_SIZE(v1,v2,n)    \
    {SimTK_ASSERT1_ALWAYS(!SimTK::Test::numericallyEqual((v1),(v2),(n)),   \
     "Test values should NOT have been numerically equivalent at size=%d times default tolerance.",(n));}

/// Test that two numerical values are NOT equal to within a specified numerical
/// error tolerance, using a relative and absolute error tolerance. In the
/// case of composite types, the equality test is performed elementwise.
#define SimTK_TEST_NOTEQ_TOL(v1,v2,tol)    \
    {SimTK_ASSERT1_ALWAYS(!SimTK::Test::numericallyEqual((v1),(v2),1,(tol)),   \
     "Test values should NOT have been numerically equivalent at tolerance=%g.",(tol));}

#ifndef SimTK_TEST_SUPPRESS_EXPECTED_THROW

/// Test that the supplied statement throws an std::exception of some kind.
#define SimTK_TEST_MUST_THROW(stmt)             \
    do {int threw=0; try {stmt;}                \
        catch(const std::exception&){threw=1;}  \
        catch(...){threw=2;}                    \
        if (threw==0) SimTK_TEST_FAILED1("Expected statement\n----\n%s\n----\n  to throw an exception but it did not.",#stmt); \
        if (threw==2) SimTK_TEST_FAILED1("Expected statement\n%s\n  to throw an std::exception but it threw something else.",#stmt); \
    }while(false)

/// Test that the supplied statement throws an std::exception of some kind, and
/// show what message got thrown.
#define SimTK_TEST_MUST_THROW_SHOW(stmt)        \
    do {int threw=0; try {stmt;}                \
        catch(const std::exception& e) {threw=1; \
            std::cout << "(OK) Threw: " << e.what() << std::endl;}  \
        catch(...){threw=2;}                    \
        if (threw==0) SimTK_TEST_FAILED1("Expected statement\n----\n%s\n----\n  to throw an exception but it did not.",#stmt); \
        if (threw==2) SimTK_TEST_FAILED1("Expected statement\n%s\n  to throw an std::exception but it threw something else.",#stmt); \
    }while(false)

/// Test that the supplied statement throws a particular exception.
#define SimTK_TEST_MUST_THROW_EXC(stmt,exc)     \
    do {int threw=0; try {stmt;}                \
        catch(const exc&){threw=1;}             \
        catch(...){threw=2;}                    \
        if (threw==0) SimTK_TEST_FAILED1("Expected statement\n----\n%s\n----\n  to throw an exception but it did not.",#stmt); \
        if (threw==2) SimTK_TEST_FAILED2("Expected statement\n----\n%s\n----\n  to throw exception type %s but it threw something else.",#stmt,#exc); \
    }while(false)

/// Allow the supplied statement to throw any std::exception without failing.
#define SimTK_TEST_MAY_THROW(stmt)             \
    do {int threw=0; try {stmt;}                \
        catch(const std::exception&){threw=1;}  \
        catch(...){threw=2;}                    \
        if (threw==2) SimTK_TEST_FAILED1("Expected statement\n%s\n  to throw an std::exception but it threw something else.",#stmt); \
    }while(false)

/// Allow the supplied statement to throw a particular exception without failing.
#define SimTK_TEST_MAY_THROW_EXC(stmt,exc)     \
    do {int threw=0; try {stmt;}                \
        catch(const exc&){threw=1;}             \
        catch(...){threw=2;}                    \
        if (threw==2) SimTK_TEST_FAILED2("Expected statement\n----\n%s\n----\n  to throw exception type %s but it threw something else.",#stmt,#exc); \
    }while(false)

// When we're only required to throw in Debug, we have to suppress the
// test case altogether in Release because it may cause damage. 
#ifdef NDEBUG
    /// Include a bad statement when in Debug and insist that it get caught,
    /// but don't include the statement at all in Release.
    #define SimTK_TEST_MUST_THROW_DEBUG(stmt)
    /// Include a bad statement when in Debug and insist that it get caught,
    /// but don't include the statement at all in Release.
    #define SimTK_TEST_MUST_THROW_EXC_DEBUG(stmt,exc)
#else
    /// Include a bad statement when in Debug and insist that it get caught,
    /// but don't include the statement at all in Release.
    #define SimTK_TEST_MUST_THROW_DEBUG(stmt) SimTK_TEST_MUST_THROW(stmt)
    /// Include a bad statement when in Debug and insist that it get caught,
    /// but don't include the statement at all in Release.
    #define SimTK_TEST_MUST_THROW_EXC_DEBUG(stmt,exc) \
                SimTK_TEST_MUST_THROW_EXC(stmt,exc)
#endif

#else // expected throws are suppressed
#define SimTK_TEST_MUST_THROW(stmt)
#define SimTK_TEST_MUST_THROW_SHOW(stmt)
#define SimTK_TEST_MUST_THROW_EXC(stmt,exc)
#define SimTK_TEST_MAY_THROW(stmt)
#define SimTK_TEST_MAY_THROW_EXC(stmt,exc)
#define SimTK_TEST_MUST_THROW_DEBUG(stmt)
#define SimTK_TEST_MUST_THROW_EXC_DEBUG(stmt,exc)
#endif


//  End of Regression testing group.
/// @}

#endif // SimTK_SimTKCOMMON_TESTING_H_
