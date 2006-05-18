#ifndef SimTK_SIMBODY_ANALYTIC_GEOMETRY_H_
#define SimTK_SIMBODY_ANALYTIC_GEOMETRY_H_

/* Copyright (c) 2006 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/** @file
 * This is the client-side interface to analytic geometry suitable for
 * using to generate mass properties, collisions, etc. Analytic geometry
 * types can also generate decorative geometry for visualization purposes.
 *
 * Each object has its own local coordinate system and is defined self-consistently
 * but independent of anything else. Clients can associate these with a Transform
 * to produce the effect of attaching AnalyticGeometry objects to bodies.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/DecorativeGeometry.h"

#include <cassert>

namespace SimTK {

/**
 * This abstract class represents a piece of high-quality geometry that
 * can be used for valid physical simulation. This is distinct from
 * DecorativeGeometry which is used for making animations. However, you
 * can use AnalyticGeometry to generate DecorativeGeometry although
 * not the reverse.
 */
class SimTK_SIMBODY_API AnalyticGeometry {
public:
    AnalyticGeometry() : rep(0) { }
    ~AnalyticGeometry();
    AnalyticGeometry(const AnalyticGeometry&);
    AnalyticGeometry& operator=(const AnalyticGeometry&);

    void setPlacement(const Transform& X_BG);
    const Transform& getPlacement() const;

    DecorativeGeometry generateDecorativeGeometry() const;

    /// Is this handle the owner of this rep? This is true if the
    /// handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // Internal use only
    explicit AnalyticGeometry(class AnalyticGeometryRep* r) : rep(r) { }
    bool hasRep() const {return rep!=0;}
    const AnalyticGeometryRep& getRep() const {assert(rep); return *rep;}
    AnalyticGeometryRep&       updRep()       {assert(rep); return *rep;}
protected:
    class AnalyticGeometryRep* rep;
};

class SimTK_SIMBODY_API AnalyticCurve : public AnalyticGeometry {
public:
    AnalyticCurve() { }
    Real calcArcLength() const;
    Vec3 calcPositionFromArcLengthParameter(const Real& s) const;

    bool isClosed() const;

    SimTK_PIMPL_DOWNCAST(AnalyticCurve, AnalyticGeometry);
};

class SimTK_SIMBODY_API AnalyticSurface : public AnalyticGeometry {
public:
    AnalyticSurface() { }

    Real calcArea() const;
    Vec3 calcPositionFromSurfaceParameters(const Vec2&) const;
    Vec3 calcNormalFromSurfaceParameters(const Vec2&) const;

    SimTK_PIMPL_DOWNCAST(AnalyticSurface, AnalyticGeometry);
};

class SimTK_SIMBODY_API AnalyticVolume : public AnalyticGeometry {
public:
    AnalyticVolume() { }
    Real calcVolume() const;

    /// In the local frame of this volume, is this point inside or out?
    bool isPointInside(const Vec3&) const;

    SimTK_PIMPL_DOWNCAST(AnalyticVolume, AnalyticGeometry);
};

/// An analytic line has only a length. The line's origin is at its
/// center, with the line running along the x axis. The arc length
/// goes from -length/2 to length/2 along x.
class SimTK_SIMBODY_API AnalyticLine : public AnalyticCurve {
public:
    AnalyticLine() { }
    AnalyticLine(Real length);

    SimTK_PIMPL_DOWNCAST(AnalyticLine, AnalyticGeometry);
};

/// An analytic circle has only a radius. The center of the circle
/// is the origin of its frame. The x axis points to the arc length 0
/// location (that is, arc length zero is point (r,0,0). The
/// y axis is the normal to the circle's plane. The arc length
/// increases counterclockwise looking down the normal at the
/// circle (right hand rule around y). z thus points at the 
/// arc length 3*pi*r point at (0,0,r).
class SimTK_SIMBODY_API AnalyticCircle : public AnalyticCurve {
public:
    AnalyticCircle() { }
    AnalyticCircle(Real radius);

    SimTK_PIMPL_DOWNCAST(AnalyticCircle, AnalyticGeometry);
};

class SimTK_SIMBODY_API AnalyticSphere : public AnalyticVolume {
public:
    AnalyticSphere() { }
    AnalyticSphere(Real radius);

    SimTK_PIMPL_DOWNCAST(AnalyticSphere, AnalyticGeometry);
};

/// The coordinate frame of the central cross section is the same as
/// for a circle; that is, x and z are radial and y points along
/// the cylinder's axis. This supports a cylindrical coordinate
/// system (h, theta), with height -halfLength <= h <= halfLength.
class SimTK_SIMBODY_API AnalyticCylinder : public AnalyticVolume {
public:
    AnalyticCylinder() { }
    AnalyticCylinder(Real radius, Real halfLength);

    SimTK_PIMPL_DOWNCAST(AnalyticCylinder, AnalyticGeometry);
};

/// This is a rectangular solid. It's local coordinate system
/// origin is at its center. Its dimensions are specified by
/// giving the half-length in x,y,z.
class SimTK_SIMBODY_API AnalyticBrick : public AnalyticVolume {
public:
    AnalyticBrick() { }
    AnalyticBrick(const Vec3& xyzHalfLengths);

    SimTK_PIMPL_DOWNCAST(AnalyticBrick, AnalyticGeometry);
};


} // namespace SimTK

#endif // SimTK_SIMBODY_ANALYTIC_GEOMETRY_H_
