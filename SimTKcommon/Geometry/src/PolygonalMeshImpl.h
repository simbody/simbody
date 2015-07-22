#ifndef SimTK_SimTKCOMMON_POLYGONAL_MESH_IMPL_H_
#define SimTK_SimTKCOMMON_POLYGONAL_MESH_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-13 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

#include "SimTKcommon/internal/PolygonalMesh.h"
#include "SimTKcommon/internal/Array.h"

namespace SimTK {

/**
 * This is the internal implementation of PolygonalMesh.
 */
class SimTK_SimTKCOMMON_EXPORT PolygonalMeshImpl
:   public PIMPLImplementation<PolygonalMesh, PolygonalMeshImpl> {
public:
    PolygonalMeshImpl() {faceVertexStart.push_back(0);}
    ~PolygonalMeshImpl() {}
    PolygonalMeshImpl* clone() const{return new PolygonalMeshImpl(*this);}
    void clear() {
        vertices.clear(); faceVertexIndex.clear(); faceVertexStart.clear();
        faceVertexStart.push_back(0);
    }
    Array_<Vec3>    vertices;
    Array_<int>     faceVertexIndex;
    Array_<int>     faceVertexStart;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_POLYGONAL_MESH_IMPL_H_
