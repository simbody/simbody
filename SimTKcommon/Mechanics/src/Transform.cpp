/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Paul Mitiguy                                                 *
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

/**@file
 * Implementations of non-inline methods of classes dealing with Transforms.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/internal/Transform.h"

//-------------------------------------------------------------------
namespace SimTK {

// Instantiate Transforms and InverseTransforms for float and double. This
// will catch bugs now rather than when some poor user is the first to
// instantiate some broken method.
template class Transform_<float>;
template class Transform_<double>;
template class InverseTransform_<float>;
template class InverseTransform_<double>;


// Define the stream output operators and instantiate them for float and
// double Transforms.
template <class P> std::ostream&
operator<<(std::ostream& o, const Transform_<P>& x )
{   return o << x.asMat34() << Row<4,P>(0,0,0,1) << std::endl;}
template <class P> std::ostream&
operator<<(std::ostream& o, const InverseTransform_<P>& x )
{   return o << Transform_<P>(x);}

template SimTK_SimTKCOMMON_EXPORT std::ostream&
operator<<(std::ostream& o, const Transform_<float>& x );
template SimTK_SimTKCOMMON_EXPORT std::ostream&
operator<<(std::ostream& o, const Transform_<double>& x );

template SimTK_SimTKCOMMON_EXPORT std::ostream&
operator<<(std::ostream& o, const InverseTransform_<float>& x );
template SimTK_SimTKCOMMON_EXPORT std::ostream&
operator<<(std::ostream& o, const InverseTransform_<double>& x );

//------------------------------------------------------------------------------
}  // End of namespace SimTK


