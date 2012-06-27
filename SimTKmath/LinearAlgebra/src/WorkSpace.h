#ifndef SimTK_SIMMATH_WORK_SPACE_H_
#define SimTK_SIMMATH_WORK_SPACE_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
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

