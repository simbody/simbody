/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2014 Stanford University and the Authors.           *
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

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/internal/Geo.h"
#include "simmath/internal/Geo_Point.h"
#include "simmath/internal/Geo_Sphere.h"
#include "simmath/internal/Geo_Box.h"
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
//                      CONTACT GEOMETRY :: BRICK
//==============================================================================

ContactGeometry::Brick::Brick(const Vec3& halfLengths)
:   ContactGeometry(new Brick::Impl(halfLengths)) {}

void ContactGeometry::Brick::setHalfLengths(const Vec3& halfLengths)
{   updImpl().setHalfLengths(halfLengths); }

const Vec3& ContactGeometry::Brick::getHalfLengths() const
{   return getImpl().getHalfLengths(); }

const Geo::Box& ContactGeometry::Brick::getGeoBox() const {
    return getImpl().getGeoBox();
}

/*static*/ ContactGeometryTypeId ContactGeometry::Brick::classTypeId()
{   return ContactGeometry::Brick::Impl::classTypeId(); }


//==============================================================================
//                               BRICK IMPL
//==============================================================================


const ContactGeometry::Brick::Impl& ContactGeometry::Brick::
getImpl() const {
    assert(impl);
    return static_cast<const Brick::Impl&>(*impl);
}

ContactGeometry::Brick::Impl& ContactGeometry::Brick::
updImpl() {
    assert(impl);
    return static_cast<Brick::Impl&>(*impl);
}

DecorativeGeometry ContactGeometry::Brick::Impl::
createDecorativeGeometry() const {
    return DecorativeBrick(getHalfLengths());
}

Vec3 ContactGeometry::Brick::Impl::
findNearestPoint(const Vec3& position, bool& inside, UnitVec3& normal) const {
    SimTK_ASSERT_ALWAYS(!"implemented",
                        "ContactGeometry::Brick::Impl::findNearestPoint()");
    return Vec3(NaN);

}

bool ContactGeometry::Brick::Impl::
intersectsRay(const Vec3& origin, const UnitVec3& direction,
              Real& distance, UnitVec3& normal) const {
    SimTK_ASSERT_ALWAYS(!"implemented",
                        "ContactGeometry::Brick::Impl::intersectsRay()");
    return false;
}

void ContactGeometry::Brick::Impl::
getBoundingSphere(Vec3& center, Real& radius) const {
    center = Vec3(0);
    radius = getHalfLengths().norm(); // length of diagonal
}


// Just an axis-aligned leaf box.
void ContactGeometry::Brick::Impl::createOBBTree() {
    OBBNode& root = obbTree.updRoot();
    root.box.setHalfLengths(getHalfLengths());
    root.normal = UnitVec3(XAxis);  // doesn't matter
    root.coneHalfAngle = Pi;        // has all possible normals
    root.pointOnSurface = Vec3(getHalfLengths()[0],0,0); // doesn't matter
    root.children.clear(); // This is a leaf

    // Leaf contents.
    root.centerUW = Vec2(NaN);
    root.dims     = Vec2(NaN);
}

