#ifndef SimTK_SIMMATH_WORK_SPACE_H_
#define SimTK_SIMMATH_WORK_SPACE_H_
/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmath(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-10 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
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
#include <cstdio>
#include <complex>
#include <cassert>
#include <iostream>


namespace SimTK {

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
#endif // SimTK_SIMMATH_WORK_SPACE_H_

