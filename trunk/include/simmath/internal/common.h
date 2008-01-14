#ifndef SIMMATH_COMMON_H_
#define SIMMATH_COMMON_H_ 

/* Portions copyright (c) 2006 Stanford University and Jack Middleton.
 * Contributors:
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
 * IN NO EVENT SHALL THE AUTHORS, COPYRIGHT HOLDERS, OR CONTRIBUTORS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/** @file
 * This is the header file that user code should include to pick up the 
 * SimTK Simmath numerical differentiation tools.
 */

#include "SimTKcommon.h"

/* Shared libraries are messy in Visual Studio. We have to distinguish three
 * cases:
 *   (1) this header is being used to build the simmath shared library (dllexport)
 *   (2) this header is being used by a *client* of the simmath shared
 *       library (dllimport)
 *   (3) we are building the simmath static library, or the client is
 *       being compiled with the expectation of linking with the
 *       simmath static library (nothing special needed)
 * In the CMake script for building this library, we define one of the symbols
 *     SIMMATH_BUILDING_{SHARED|STATIC}_LIBRARY
 * Client code normally has no special symbol defined, in which case we'll
 * assume it wants to use the shared library. However, if the client defines
 * the symbol SimTK_USE_STATIC_LIBRARIES we'll suppress the dllimport so
 * that the client code can be linked with static libraries. Note that
 * the client symbol is not library dependent, while the library symbols
 * affect only the simmath library, meaning that other libraries can
 * be clients of this one. However, we are assuming all-static or all-shared.
*/

#ifdef WIN32
    #if defined(SimTK_SIMMATH_BUILDING_SHARED_LIBRARY)
        #define SimTK_SIMMATH_EXPORT __declspec(dllexport)
    #elif defined(SimTK_SIMMATH_BUILDING_STATIC_LIBRARY) || defined(SimTK_USE_STATIC_LIBRARIES)
        #define SimTK_SIMMATH_EXPORT
    #else
        /* i.e., a client of a shared library */
        #define SimTK_SIMMATH_EXPORT __declspec(dllimport)
    #endif
#else
    /* Linux, Mac */
    #define SimTK_SIMMATH_EXPORT
#endif


// Every SimTK Core library must provide these two routines, with the library
// name appearing after the "version_" and "about_".
extern "C" {
    SimTK_SIMMATH_EXPORT void SimTK_version_simmath(int* major, int* minor, int* build);
    SimTK_SIMMATH_EXPORT void SimTK_about_simmath(const char* key, int maxlen, char* value);
}



const static double POSITIVE_INF =  2e19;
const static double NEGATIVE_INF = -2e19;

namespace SimTK {
template <class P>
bool calcEigenValuesRightEigenVectors( Matrix_<P> &m, Vector_< std::complex<P> > &eigenValues, Matrix_< std::complex<P> > &eigenVectors ); 


namespace Exception {

class OptimizerFailed : public Base {
public:
        OptimizerFailed( const char * fn, int ln, String msg) : Base(fn, ln)
        {
            setMessage("Optmizer failed: " + msg );
        }
private:
};

class UnrecognizedParameter : public Base {
public:
        UnrecognizedParameter( const char * fn, int ln, String msg) : Base(fn, ln)
        {
            setMessage("Unrecognized Parameter: " + msg );
        }
private:
};

class IllegalLapackArg : public Base {
public:
        IllegalLapackArg( const char *fn, int ln, const char *lapackRoutine, 
                  int argNum ) : Base(fn, ln)
        {
        char buf[1024];

        sprintf(buf, "SimTK internal error: %s called with an illegal value to arguement #%d \n Please report this to SimTK",
            lapackRoutine, argNum );
        setMessage(String(buf));

        }
private:
};
class IncorrectArrayLength : public Base {
public:
        IncorrectArrayLength( const char *fn, int ln, const char *valueName, int length,  
                              const char *paramName, int paramValue, const char *where) : Base(fn, ln)
        {
        char buf[1024];

        sprintf(buf, "Incorrect array length in %s : %s is %d and must equal %s which is %d",
            where, valueName, length, paramName, paramValue );
        setMessage(String(buf));

        }
private:
};

class SingularMatrix : public Base {
public:
        SingularMatrix( const char *fn, int ln, int index,  
                               const char *where) : Base(fn, ln)
        {
        char buf[1024];

        sprintf(buf, "%s failed because index %d in matrix was singular and factorization failed",
            where, index );
        setMessage(String(buf));

        }
private:
};

class ConvergedFailed : public Base {
public:
        ConvergedFailed( const char *fn, int ln, const char *algorithm,  
                               const char *where) : Base(fn, ln)
        {
        char buf[1024];

        sprintf(buf, "%s failed because %s failed to converge", where, algorithm );
        setMessage(String(buf));

        }
private:
};

class NotPositiveDefinite : public Base {
public:
        NotPositiveDefinite( const char *fn, int ln, int index,  
                               const char *where) : Base(fn, ln)
        {
        char buf[1024];

        sprintf(buf, "%s failed because index %d in matrix was not positive definite and factorization failed ",
            where, index );
        setMessage(String(buf));

        }
private:
};
} // namespace Exception

} //  namespace SimTK

#endif // SIMMATH_COMMON_H_
