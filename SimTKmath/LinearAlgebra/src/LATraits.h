#ifndef SimTK_SIMMATH_LA_TRAITS_H_
#define SimTK_SIMMATH_LA_TRAITS_H_

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

#include <complex>
#include "SimTKcommon.h"

namespace SimTK {

template <typename T> struct LANT;
        
template <> struct LANT<float>{
     static const int sign =  1;
     static const bool conjugate = false;
};
template <> struct LANT<double>{
     static const int sign     = 1;
     static const bool conjugate = false;
};
template <> struct LANT<std::complex<float> >{
     static const int sign    = 1;
     static const bool conjugate = false;
};
template <> struct LANT<std::complex<double> >{
     static const int sign     = 1;
     static const bool conjgate = false;
};
template <> struct LANT<negator<float> >{
     static const int sign    = -1;
     static const bool conjugate = false;
};
template <> struct LANT<negator<double> >{
     static const int sign     = -1;
     static const bool conjugate = false;
};
template <> struct LANT<negator<std::complex<float> > >{
     static const int sign    = -1;
     static const bool conjugate = false;
};
template <> struct LANT<negator<std::complex<double> > >{
     static const int sign     = -1;
     static const bool conjugate = false;
};
template <> struct LANT<conjugate<float> >{
     static const int sign    = 1;
     static const bool conjugate = true;
};
template <> struct LANT<conjugate<double> >{
     static const int sign     = 1;
     static const bool conjugate = true;
};
template <> struct LANT<negator<conjugate<float> > >{
     static const int sign    = -1;
     static const bool conjugate = true;
};
template <> struct LANT<negator<conjugate<double> > >{
     static const int sign     = -1;
     static const bool conjugate = true;
};

} // namespace SimTK
#endif   // SimTK_SIMMATH_LA_TRAITS_H_
