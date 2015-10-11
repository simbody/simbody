/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
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
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"
#include "simmath/internal/ContactGeometry.h"

#include "ContactGeometryImpl.h"

#include <iostream>
#include <cmath>
#include <map>
#include <set>

using namespace SimTK;
using std::map;
using std::pair;
using std::set;
using std::string;
using std::cout; using std::endl;


//==============================================================================
//                 CONTACT GEOMETRY :: HALF SPACE & IMPL
//==============================================================================
ContactGeometry::HalfSpace::HalfSpace()
:   ContactGeometry(new HalfSpace::Impl()) {}

// TODO: currently the only possible normal is -XAxis.
UnitVec3 ContactGeometry::HalfSpace::getNormal() const {
    return UnitVec3(-XAxis);
}

/*static*/ ContactGeometryTypeId ContactGeometry::HalfSpace::classTypeId()
{   return ContactGeometry::HalfSpace::Impl::classTypeId(); }

DecorativeGeometry ContactGeometry::HalfSpace::Impl::createDecorativeGeometry() const {
    return DecorativeBrick(Vec3(Real(0.01),1,1));
}

// Point position is given in the half space frame.
Vec3 ContactGeometry::HalfSpace::Impl::findNearestPoint
   (const Vec3& position, bool& inside, UnitVec3& normal) const {
    inside = (position[0] >= 0);
    normal = -UnitVec3(XAxis); // this does not require normalization
    return Vec3(0, position[1], position[2]);
}

bool ContactGeometry::HalfSpace::Impl::intersectsRay
   (const Vec3& origin, const UnitVec3& direction,
    Real& distance, UnitVec3& normal) const
{
    if (std::abs(direction[0]) < SignificantReal)
        return false; // ray is parallel to halfspace surface

    const Real t = origin[0]/direction[0];
    if (t > 0)
        return false; // ray points away from surface

    distance = -t;
    normal = -UnitVec3(XAxis); // cheap; no normalization required
    return true;
}

void ContactGeometry::HalfSpace::Impl::getBoundingSphere
   (Vec3& center, Real& radius) const
{   center = Vec3(0);
    radius = Infinity; }

const ContactGeometry::HalfSpace::Impl& ContactGeometry::HalfSpace::
getImpl() const {
    assert(impl);
    return static_cast<const HalfSpace::Impl&>(*impl);
}

ContactGeometry::HalfSpace::Impl& ContactGeometry::HalfSpace::
updImpl() {
    assert(impl);
    return static_cast<HalfSpace::Impl&>(*impl);
}
