#ifndef SimTK_SIMMATH_ORIENTED_BOUNDING_BOX_H_
#define SimTK_SIMMATH_ORIENTED_BOUNDING_BOX_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
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

namespace SimTK {

/**
 * This class represents a rectangular box with arbitrary position and 
 * orientation.  It is used in collision detection as a bounding volume for 
 * geometry of various types.
 *
 * An OrientedBoundingBox is defined by a Transform that specifies its position
 * and orientation, and a Vec3 that specifies its size. In the reference frame
 * defined by the Transform, one corner is at the origin and the opposite 
 * corner is at the point returned by getSize().
 */
class SimTK_SIMMATH_EXPORT OrientedBoundingBox {
public:
    OrientedBoundingBox();
    /**
     * Create an OrientedBoundingBox.
     *
     * @param transform     specifies the position and orientation of the box
     * @param size          specifies the dimensions of the box
     */
    OrientedBoundingBox(const Transform& transform, const Vec3& size);
    /**
     * Create an OrientedBoundingBox which encloses a set of points.
     */
    explicit OrientedBoundingBox(const Vector_<Vec3>& points);
    /**
     * Get the position and orientation of the box.
     */
    const Transform& getTransform() const;
    /**
     * Get the dimensions of the box.
     */
    const Vec3& getSize() const;
    /**
     * Determine whether a point is inside the box.
     */
    bool containsPoint(const Vec3& point) const;
    /**
     * Determine whether this box intersects another bounding box at any point.
     */
    bool intersectsBox(const OrientedBoundingBox& box) const;
    /**
     * Determine whether a ray intersects this bounding box.
     *
     * @param origin     the position at which the ray begins
     * @param direction  the ray direction
     * @param distance   if an intersection is found, the distance from the ray
     *                   origin to the intersection point is stored in this.
     *                   Otherwise, it is left unchanged.
     * @return true if an intersection is found, false otherwise
     */
    bool intersectsRay(const Vec3& origin, const UnitVec3& direction, 
                       Real& distance) const;
    /**
     * Given a point in space, find the point inside the bounding box which is 
     * nearest to it.
     */
    Vec3 findNearestPoint(const Vec3& position) const;
    /**
     * Get the locations of the eight corners of the box.
     *
     * @param corners   the corner locations are stored in this array
     */
    void getCorners(Vec3 corners[8]) const;
private:
    Real calculateVolume(const Vector_<Vec3>& points, const Rotation& rotation);
    Transform transform;
    Vec3 size;
};

SimTK_SIMMATH_EXPORT OrientedBoundingBox 
operator*(const Transform& t, const OrientedBoundingBox& box);

} // namespace SimTK

#endif // SimTK_SIMMATH_ORIENTED_BOUNDING_BOX_H_
