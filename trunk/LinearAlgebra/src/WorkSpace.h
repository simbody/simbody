#ifndef _WORK_SPACE_H
#define _WORK_SPACE_H
/* Portions copyright (c) 2007 Stanford University and Jack Middleton.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include <cstdio>
#include <complex>
#include <cassert>
#include <iostream>
#include "SimTKcommon.h"

namespace SimTK {

/*  the WorkSpace class was ment to be used as a resource manager 
    class for generic  blocks of memory the TypeWorkSpace class
    is easier to use. The code for WorkSpace and the associated traits
    should be deleted once we are sure we do not need them.


struct WorkSpaceTypes {
     enum eltTypes {
         IntWA           = 1,
         FloatWA         = 2,
         DoubleWA        = 3,
         ComplexDoubleWA = 4,
         ComplexFloatWA  = 5,
     };

};
template <typename T> struct WSTraits;
template <> struct WSTraits<int>{     
    static const WorkSpaceTypes::eltTypes dType = WorkSpaceTypes::IntWA; 
    static int const size = sizeof(int);
};
template <> struct WSTraits<float> {  
    static const WorkSpaceTypes::eltTypes dType = WorkSpaceTypes::FloatWA; 
    static int const size = sizeof(float);
};
template <> struct WSTraits<double> { 
    static const WorkSpaceTypes::eltTypes dType = WorkSpaceTypes::DoubleWA; 
    static int const size = sizeof(double);
};
template <> struct WSTraits<std::complex<double> >{ 
    static const WorkSpaceTypes::eltTypes dType = WorkSpaceTypes::ComplexDoubleWA;
    static int const size = sizeof(std::complex<double>);
};
template <> struct WSTraits<std::complex<float> >{ 
    static const WorkSpaceTypes::eltTypes dType = WorkSpaceTypes::ComplexFloatWA; 
    static int const size = sizeof(std::complex<float>);
};

class WorkSpace {

    public:
    WorkSpace(long size) {

   //     data = (void*)malloc(size);
        data = (void*)new char[size];
// ***  Throw exception if malloc was unsuccessful
        lwork = size;
    }
    template <typename T>
    static long WorkStorage( Matrix_<T> mat) {
        return( WSTraits<T>::size*mat.size() );
    }
    template <typename T>
    static long PivotStorage( Matrix_<T> mat) {
        return( WSTraits<int>::size*mat.ncol() );
    }
    template <typename T>
    static long WorkStorage( Vector_<T> vec) {
        return( WSTraits<T>::size*vec.size() );
    }
    template <typename T>
    WorkSpace(Matrix_<T> mat) {
        lwork =  WorkStorage( mat );
        data = (void*)malloc(lwork );
    }
    template <typename T>
    WorkSpace(Vector_<T> vec) {
        lwork =  WorkStorage( vec );
        data = (void*)malloc(lwork );
    }

    ~WorkSpace() {
//        free(data);
          delete [] data;
    }

    private:
   
    int lwork;
    void*  data;
};
*/
template <typename T>
class TypedWorkSpace {
    public:

    // copy constructor
    TypedWorkSpace( const TypedWorkSpace& c ) {
        size = c.size;

        if( size == 0 ) {
             data = 0;
         } else {
             data = new T[size];
             for(int i=0;i<size;i++) data[i] = c.data[i];
         }
       
    }
    TypedWorkSpace& operator=(const TypedWorkSpace& rhs) {
        if (&rhs == this)
            return *this;

        delete [] data;
        data = 0;
        size = rhs.size;

        if( size > 0) {
             data = new T[size];
             for(int i=0;i<size;i++) data[i] = rhs.data[i];
        }
        return *this;
    }

    explicit TypedWorkSpace( int n ) {
        size = n;
        data = (n==0 ? 0 : new T[n]);
    }

    TypedWorkSpace() : size(0), data(0) { }

    ~TypedWorkSpace() {
        delete [] data;
    }
    
    void resize( int n ) {
        delete [] data;
        size = n;
        data = (n==0 ? 0 : new T[n]);
    }

    int size;
    T* data; 
};

} // namespace SimTK
#endif   //  _FACTOR_REP_H_
