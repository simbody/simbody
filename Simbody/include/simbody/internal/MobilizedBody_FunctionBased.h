#ifndef SimTK_SIMBODY_MOBILIZED_BODY_FUNCTIONBASED_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_FUNCTIONBASED_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-13 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Ajay Seth                                                    *
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
Declares the MobilizedBody::FunctionBased class. **/

#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/MobilizedBody_Custom.h"

namespace SimTK {

/** This is a subclass of MobilizedBody::Custom which uses a set of Function 
objects to define the behavior of the %MobilizedBody. When you create it, you 
specify the number of generalized coordinates, and six Function objects which
calculate the spatial rotations and translations based on those coordinates.  
It assumes there is a one to one correspondence between generalized coordinates 
and generalized speeds, so qdot == u.

Each of the Function objects must take some subset of the generalized 
coordinates as inputs, and produce a single number as its output. It also must 
support derivatives up to second order.  Taken together, the six functions
define a SpatialVec giving the body's mobilizer transform. **/
class SimTK_SIMBODY_EXPORT MobilizedBody::FunctionBased 
:   public MobilizedBody::Custom {
public:
    /** Default constructor provides an empty handle that can be assigned to
    reference any %MobilizedBody::FunctionBased. **/
    FunctionBased() {}

    /** Create a FunctionBased MobilizedBody.
    
    @param parent         the MobilizedBody's parent body
    @param bodyInfo       describes this MobilizedBody's physical properties
    @param nmobilities    the number of generalized coordinates belonging to this MobilizedBody
    @param functions      the Functions describing how the body moves based on its generalized coordinates.
                          This must be of length 6.  The elements correspond to, in order, x rotation, y rotation, z rotation,
                          x translation, y translation, and z translation.  The MobilizedBody takes over ownership of the functions,
                          and automatically deletes them when the MobilizedBody is deleted.
    @param coordIndices   the indices of the generalized coordinates that are inputs to each function.  For example, if coordIndices[2] = {0, 1},
                          that means that functions[2] takes two input arguments, and q[0] and q[1] respectively should be passed as those arguments.
    @param direction      whether you want the coordinates defined as though parent & child were swapped
    @see MobilizedBody for a diagram and explanation of terminology. **/
    FunctionBased(MobilizedBody& parent, const Body& bodyInfo, 
                  int nmobilities, const Array_<const Function*>& functions,
                  const Array_<Array_<int> >& coordIndices,
                  Direction direction=Forward);

    /** For compatibility with std::vector. **/
    FunctionBased(MobilizedBody& parent, const Body& bodyInfo, 
                  int nmobilities, const std::vector<const Function*>& functions,
                  const std::vector<std::vector<int> >& coordIndices,
                  Direction direction=Forward) 
    {
        Array_< Array_<int> > coordCopy(coordIndices); // sorry, must copy
        // Use the above constructor.
        new(this) FunctionBased(parent,bodyInfo,nmobilities,
                                ArrayViewConst_<const Function*>(functions), 
                                coordCopy, direction);
    }

    /** Create a FunctionBased MobilizedBody.
    
    @param parent         the MobilizedBody's parent body
    @param X_PF           the default inboard frame
    @param bodyInfo       describes this MobilizedBody's physical properties
    @param X_BM           the default outboard frame
    @param nmobilities    the number of generalized coordinates belonging to this MobilizedBody
    @param functions      the Functions describing how the body moves based on its generalized coordinates.
                          This must be of length 6.  The elements correspond to, in order, x rotation, y rotation, z rotation,
                          x translation, y translation, and z translation.  The MobilizedBody takes over ownership of the functions,
                          and automatically deletes them when the MobilizedBody is deleted.
    @param coordIndices   the indices of the generalized coordinates that are inputs to each function.  For example, if coordIndices[2] = {0, 1},
                          that means that functions[2] takes two input arguments, and q[0] and q[1] respectively should be passed as those arguments.
    @param direction      whether you want the coordinates defined as though parent & child were swapped
    @see MobilizedBody for a diagram and explanation of terminology. **/
    FunctionBased(MobilizedBody& parent, const Transform& X_PF, 
                  const Body& bodyInfo, const Transform& X_BM, 
                  int nmobilities, const Array_<const Function*>& functions,
                  const Array_<Array_<int> >& coordIndices,
                  Direction direction=Forward);

    /** For compatibility with std::vector. **/
    FunctionBased(MobilizedBody& parent, const Transform& X_PF, 
                  const Body& bodyInfo, const Transform& X_BM, 
                  int nmobilities, const std::vector<const Function*>& functions,
                  const std::vector<std::vector<int> >& coordIndices,
                  Direction direction=Forward)
    {
        Array_< Array_<int> > coordCopy(coordIndices); // sorry, must copy
        // Use the above constructor.
        new(this) FunctionBased(parent,X_PF,bodyInfo,X_BM,
                                nmobilities, ArrayViewConst_<const Function*>(functions), 
                                coordCopy, direction);
    }

    /** Create a FunctionBased MobilizedBody.
    
    @param parent         the MobilizedBody's parent body
    @param bodyInfo       describes this MobilizedBody's physical properties
    @param nmobilities    the number of generalized coordinates belonging to this MobilizedBody
    @param functions      the Functions describing how the body moves based on its generalized coordinates.
                          This must be of length 6.  The elements correspond to, in order, x rotation, y rotation, z rotation,
                          x translation, y translation, and z translation.  The MobilizedBody takes over ownership of the functions,
                          and automatically deletes them when the MobilizedBody is deleted.
    @param coordIndices   the indices of the generalized coordinates that are inputs to each function.  For example, if coordIndices[2] = {0, 1},
                          that means that functions[2] takes two input arguments, and q[0] and q[1] respectively should be passed as those arguments.
    @param axes           the axes directions (as Vec3's) for each spatial coordinate, which each function describes, and is therefore length 6.
                          First 3 and last 3 axes must be linearly independent, otherwise there will be redundant speeds for the same motion.
    @param direction      whether you want the coordinates defined as though parent & child were swapped
    @see MobilizedBody for a diagram and explanation of terminology. **/
    FunctionBased(MobilizedBody& parent, const Body& bodyInfo, 
                  int nmobilities, const Array_<const Function*>& functions,
                  const Array_<Array_<int> >& coordIndices, const Array_<Vec3>& axes,
                  Direction direction=Forward);

    /** For compatibility with std::vector. **/
    FunctionBased(MobilizedBody& parent, const Body& bodyInfo, 
                  int nmobilities, const std::vector<const Function*>& functions,
                  const std::vector<std::vector<int> >& coordIndices, const std::vector<Vec3>& axes,
                  Direction direction=Forward)
    {
        Array_< Array_<int> > coordCopy(coordIndices); // sorry, must copy
        // Use the above constructor.
        new(this) FunctionBased(parent,bodyInfo,
                                nmobilities, ArrayViewConst_<const Function*>(functions), 
                                coordCopy, ArrayViewConst_<Vec3>(axes), 
                                direction);
    }

    /** Create a FunctionBased MobilizedBody.
    
    @param parent         the MobilizedBody's parent body
    @param X_PF           the default inboard frame
    @param bodyInfo       describes this MobilizedBody's physical properties
    @param X_BM           the default outboard frame
    @param nmobilities    the number of generalized coordinates belonging to this MobilizedBody
    @param functions      the Functions describing how the body moves based on its generalized coordinates.
                          This must be of length 6.  The elements correspond to, in order, x rotation, y rotation, z rotation,
                          x translation, y translation, and z translation.  The MobilizedBody takes over ownership of the functions,
                          and automatically deletes them when the MobilizedBody is deleted.
    @param coordIndices   the indices of the generalized coordinates that are inputs to each function.  For example, if coordIndices[2] = {0, 1},
                          that means that functions[2] takes two input arguments, and q[0] and q[1] respectively should be passed as those arguments.
    @param axes           the axes directions (as Vec3's) for each spatial coordinate, which each function describes, and is therefore length 6.
                          First 3 and last 3 axes must be linearly independent, otherwise there will be redundant speeds for the same motion.
    @param direction      whether you want the coordinates defined as though parent & child were swapped
    @see MobilizedBody for a diagram and explanation of terminology. **/
    FunctionBased(MobilizedBody& parent, const Transform& X_PF, 
                  const Body& bodyInfo, const Transform& X_BM, 
                  int nmobilities, const Array_<const Function*>& functions,
                  const Array_<Array_<int> >& coordIndices, const Array_<Vec3>& axes,
                  Direction direction=Forward);

    /** For compatibility with std::vector. **/
    FunctionBased(MobilizedBody& parent, const Transform& X_PF, 
                  const Body& bodyInfo, const Transform& X_BM,
                  int nmobilities, const std::vector<const Function*>& functions,
                  const std::vector<std::vector<int> >& coordIndices, const std::vector<Vec3>& axes,
                  Direction direction=Forward)
    {
        Array_< Array_<int> > coordCopy(coordIndices); // sorry, must copy
        // Use the above constructor.
        new(this) FunctionBased(parent,X_PF,bodyInfo,X_BM,
                                nmobilities, ArrayViewConst_<const Function*>(functions), 
                                coordCopy, ArrayViewConst_<Vec3>(axes), 
                                direction);
    }
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_FUNCTIONBASED_H_



