#ifndef SimTK_SIMMATRIX_H_
#define SimTK_SIMMATRIX_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-13 Stanford University and the Authors.        *
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

/** @file
 * This is the header which should be included in user programs that would
 * like to make use of all the Simmatrix facilities, but none of the other
 * parts of SimTKcommon.
 */

// This is so Doxygen can locate the symbols we mention.
namespace SimTK {

/** @defgroup MatVecUtilities   Matrix and Vector Utilities

@brief Simbody contains an extensive library for manipulating Matrix and Vector
objects, modeled after Matlab's similar features.

Simbody's matrix library contains two separate but related sets of classes for
vector and matrix objects, one set for small, fixed-size objects and the other
for larger, run-time allocated objects.

First, there are classes to represent small, fixed size vectors
and matrices with zero runtime overhead: Vec for column vectors, and Mat for
matrices. There is also a Row type that does not normally appear in user
programs. These classes are templatized based on size and element type. Synonyms
(typedefs) are defined for common combinations; for example, 
@ref SimTK::Vec3 "Vec3" is a synonym for \c Vec<3,Real>, while 
@ref SimTK::Mat22 "Mat22" is a synonym for \c Mat<2,2,Real>. (Typedef
@ref SimTK::Real "Real" is synonymous with C++ \c double unless Simbody was 
compiled with \c float as the default precision.) You can also
create other combinations, such as \c Mat<2,10,Real> or
\c Vec<4,std::complex<Real>>. However, the size must always be determinable at
compile time. The in-memory representation of these small objects is minimal:
only the data elements are stored.

Second, there are classes to represent large vectors and matrices whose sizes
are determined at runtime: Vector_ for column vectors and Matrix_ for matrices.
There is also a RowVector_ type that rarely appears in user programs. These 
classes are templatized based on element type. In user code, it is most common
to see typedefs @ref SimTK::Vector "Vector" and @ref SimTK::Matrix "Matrix" 
which are synonyms for \c Vector_<Real> and \c Matrix_<Real>. As for small matrices,
you can use other element types. In fact, the element type can even be one of 
the fixed-size vector or matrix objects. For example, \c Vector_<Vec3> is a 
variable-length vector, where each element is itself a fixed-size, 3-component 
vector. The type \c Vec<2,Vec3>, called a spatial vector (@ref SimTK::SpatialVec
"SpatialVec"), is useful for 
combining rotational and translational quantities into a single object 
representing a spatial velocity or spatial force, for example. However, it is not permissible to use the variable-size 
Vector_ or Matrix_ objects as element types. The in-memory representation of 
these objects includes, in addition to the data, an opaque descriptor containing
the length and information on how the data is laid out; the declared objects 
actually consist only of a pointer (essentially a \c void*) to the descriptors. 
This has many advantages for implementation and binary compatibility, but makes 
it difficult to look through these objects in a debugger as you can with the 
small Vec and Mat classes.

<h2>Implementation</h2>
The intent of Simbody's matrix library, which we call \e Simmatrix, is similar
to that of the %Eigen library (http://eigen.tuxfamily.org). If you want to know
more about the design goals and implementation of Simmatrix, see the design
document here: https://simtk.org/home/simbody, Documents tab.
**/

} // namespace SimTK

// Each of these is independently user-includable, with later ones including
// former ones.
#include "SimTKcommon/Scalar.h"         // self-contained
#include "SimTKcommon/SmallMatrix.h"    // includes Scalar.h
#include "SimTKcommon/Orientation.h"    // includes SmallMatrix.h
#include "SimTKcommon/Mechanics.h"      // includes Orientation.h

// Here we add the missing pieces that provide large matrix functionality,
// and some additional small matrix functionality that depends on having
// access to large matrix capabilities.
#include "SimTKcommon/internal/BigMatrix.h"
#include "SimTKcommon/internal/SmallDefsThatNeedBig.h"
#include "SimTKcommon/internal/VectorMath.h"

#endif // SimTK_SIMMATRIX_H_
