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

        data = (void*)malloc(size);
// TODO  Throw exception if malloc was unsuccessful
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
        free(data);
    }

    template<typename T> T* getData() { assert(false); }

    private:
   
    int lwork;
    void*  data;
};

    template<> int* WorkSpace::getData(){ return( static_cast<int*>(data)); }
    template<> double* WorkSpace::getData(){ return( static_cast<double*>(data)); }
    template<> float* WorkSpace::getData(){ return( static_cast<float*>(data)); }
    template<> std::complex<double>* WorkSpace::getData(){ return( static_cast<std::complex<double>* >(data)); }
    template<> std::complex<float>* WorkSpace::getData(){ return( static_cast<std::complex<float>* >(data)); }

template <typename T>
class TypedWorkSpace {
    public:

    TypedWorkSpace( long n ) {
        size = n;
        data = new T[n];
    }

    ~TypedWorkSpace() {
        delete [] data;
    }

    long size;
    T* data; 
};

} // namespace SimTK
#endif   //  _FACTOR_REP_H_
