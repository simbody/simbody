#ifndef SimTK_SIMMATH_ORIENTED_BOUNDING_BOX_H_
#define SimTK_SIMMATH_ORIENTED_BOUNDING_BOX_H_

/* -------------------------------------------------------------------------- *
 *                        SimTK Simbody: SimTKmath                            *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008-11 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
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
