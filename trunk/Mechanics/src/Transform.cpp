//-----------------------------------------------------------------------------
// File:     Transform.cpp
// Class:    Transform and InverseTransform
// Parent:   None:  Data contains Rotation (for orientation) and Vec3 (translation)
// Purpose:  Transform (orientation and translation) relating two right-handed orthogonal bases
//-----------------------------------------------------------------------------

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simmatrix(tm)                       *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors: Paul Mitiguy                                                 *
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


