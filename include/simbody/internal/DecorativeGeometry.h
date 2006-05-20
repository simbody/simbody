#ifndef SimTK_SIMBODY_DECORATIVE_GEOMETRY_H_
#define SimTK_SIMBODY_DECORATIVE_GEOMETRY_H_

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
 * This is the client-side interface to decorative geometry suitable for
 * visualization. DO NOT confuse this with AnalyticGeometry which can represent
 * physically meaningful objects. However, any AnalyticGeometry object can
 * generate some decorative geometry for visualization.
 *
 * Each object has its own local coordinate system and is defined self-consistently
 * but independent of anything else. Clients can associate these with a reference
 * frame (e.g. a body), and place the local frame of the geometry objects on the 
 * reference frame, or at a fixed transform from the reference frame. That places
 * the DecorativeGeometry objects in a scene. We support both 3D objects
 * which are attached to actors in the scene, and 2D "screen"
 * objects like titles which are attached to the display rather than the actors.
 * The classes here deal only with the local-frame definitions of the geometric
 * objects, not their placement in the scene.
 */

#include "simbody/internal/common.h"

#include <cassert>

class vtkPolyData;

namespace SimTK {

class AnalyticGeometry;


// Some common RGB values;
static const Vec3 Black   = Vec3( 0, 0, 0);
static const Vec3 Gray    = Vec3(.5,.5,.5);
static const Vec3 Red     = Vec3( 1, 0, 0);
static const Vec3 Green   = Vec3( 0, 1, 0);
static const Vec3 Blue    = Vec3( 0, 0, 1);
static const Vec3 Yellow  = Vec3( 1, 1, 0);
static const Vec3 Orange  = Vec3( 1,.5, 0);
static const Vec3 Magenta = Vec3( 1, 0, 1);
static const Vec3 Purple  = Vec3(.5, 0,.5);
static const Vec3 Cyan    = Vec3( 0, 1, 1);
static const Vec3 White   = Vec3( 1, 1, 1);

/**
 * This is an abstract handle class using the PIMPL design pattern to hide the private
 * implementation. This is effectively an abstract class although the virtual
 * function table is hidden in the private part.
 */
class SimTK_SIMBODY_API DecorativeGeometry {
public:
    DecorativeGeometry() : rep(0) { }
    ~DecorativeGeometry();
    DecorativeGeometry(const DecorativeGeometry&);
    DecorativeGeometry& operator=(const DecorativeGeometry&);

    /// Implicit conversion
    DecorativeGeometry(const AnalyticGeometry&);

    /// This returns a VTK pipeline suitable for generating the polygons
    /// needed for display. Although this is returned non-const for use with VTK,
    /// ownership is retained by the DecorativeGeometry object, which will
    /// call VTK's Delete() method on destruction. So the caller should not attempt
    /// to Delete() the PolyData. Also note that any changes to the generating
    /// parameters INVALIDATES the pipeline (TODO: should just update it).
    /// Generating parameters include resolution, transform, scale here, and any
    /// specific source parameters like radius for a sphere.
    vtkPolyData* updVTKPolyData();

    /// Each concrete DecorativeGeometry object is expected to have a default resolution
    /// that gets the point across but is cheap to draw and hence probably somewhat "chunky".
    /// The resolution parameter here scales that default up or down. A value less than
    /// or equal to zero here is interpreted as an instruction to "use the default".
    /// This value affects the generated polygonal data.
    DecorativeGeometry& setResolution(Real);

    /// This transform shifts the generated polygons with respect to this object's
    /// local frame. Subsequent calls with other transforms simply replace the earlier
    /// one; they do not accumulate. The default transform is identity and you can 
    /// call setPlacement(Transform()) to put the transform back into its original state.
    /// This value affects the generated polygonal data.
    DecorativeGeometry& setPlacement(const Transform& X_BG);

    /// Each concrete DecorativeGeometry object is expected to have a default size
    /// around "1", whatever that means for a particular object, and most objects also
    /// allow a user-specified size on construction. The scale factor here applies to
    /// the object as the user built it, or to the default if the user didn't specify
    /// a size. The default scaling is 1, and any value less than or equal to zero
    /// here is interpreted as a request to "use the default".
    /// This value affects the generated polygonal data.
    DecorativeGeometry& setScale(Real);

    /// Return the current setting of the "resolution" factor. A return value of -1
    /// means "use the default".
    Real getResolution() const;

    /// Return the current value of the object's transform. If none has been set this
    /// will be the identity transform. Note that this transform specifies how the
    /// polygons are placed with respect to the object's local frame.
    const Transform& getPlacement() const;

    /// Return the current setting of the "scale" factor. A return value of -1 
    /// means "use the default" (which is typically 1).
    Real getScale() const;

    /// Request a specific color for this DecorativeGeometry object. This does NOT
    /// affect the generated geometry here. The default is that the color is
    /// determined elsewhere.
    DecorativeGeometry& setColor(const Vec3& rgb); // 0-1 for each color; default is 0,0,0 (black)

    /// Request a level of transparency for this DecorativeGeometry. This does NOT
    /// affect the generated geometry here. The default is that opacity is 
    /// determined elsewhere.
    DecorativeGeometry& setOpacity(Real);          // 0-1; default is 1 (opaque)

    /// Request an adjustment to the default rendering of lines and curves. This 
    /// does NOT affect geometry generated here; it is a request passed on to the
    /// renderer which will probably pass it on to the hardware. A value less
    /// than or equal to zero here is interpreted as "use the default".
    DecorativeGeometry& setLineThickness(Real);

    const Vec3& getColor()      const;
    Real        getOpacity()    const;
    Real        getLineThickness() const;

    /// Request a particular rendering of this DecorativeGeometry object as a 
    /// set of points. The default is that the rendering choice is made elsewhere.
    DecorativeGeometry& setRepresentationToPoints();

    /// Request a particular rendering of this DecorativeGeometry object in
    /// wireframe. The default is that the rendering choice is made elsewhere.
    DecorativeGeometry& setRepresentationToWireframe();

    /// Request a particular rendering of this DecorativeGeometry object using
    /// shaded surfaces. The default is that the rendering choice is made elsewhere.
    DecorativeGeometry& setRepresentationToSurface();

    /// Specify that the representation for this DecorativeGeometry object should
    /// be chosen elsewhere. This is the default.
    DecorativeGeometry& setRepresentationToUseDefault();

    /// -1 means "use default"; otherwise we're not documenting the meaning here.
    int getRepresentation() const;

    // Bookkeeping below here -- internal use only. Don't look below or you will
    // turn into a pillar of salt.

    bool isOwnerHandle() const;
    bool isEmptyHandle() const;
    explicit DecorativeGeometry(class DecorativeGeometryRep* r) : rep(r) { }
    bool hasRep() const {return rep!=0;}
    const DecorativeGeometryRep& getRep() const {assert(rep); return *rep;}
    DecorativeGeometryRep&       updRep()       {assert(rep); return *rep;}
protected:
    DecorativeGeometryRep* rep;
};

/**
 * A line between two points. Note that the actual placement can be changed 
 * by the parent class transform & scale; here we are just generating the 
 * initial line in the geometry object's local frame. 
 *
 * There is a default constructor for this object but it is not much
 * use unless followed by endpoint specifications. By default we produce
 * a line going from (0,0,0) to (1,1,1) just so it will show up if you
 * forget to set it to something meaningful. Having a default constructor
 * allows us to have arrays of these objects.
 */
class SimTK_SIMBODY_API DecorativeLine : public DecorativeGeometry {
public:
    explicit DecorativeLine(const Vec3& p1=Vec3(0), const Vec3& p2=Vec3(1)); // line between p1 and p2

    void setPoint1(const Vec3& p1);
    void setPoint2(const Vec3& p2);
    void setEndpoints(const Vec3& p1, const Vec3& p2);

    const Vec3& getPoint1() const;
    const Vec3& getPoint2() const;

    SimTK_PIMPL_DOWNCAST(DecorativeLine, DecorativeGeometry);
};

/**
 * This defines a circle in the x-y plane, centered at the origin. The
 * default constructor creates a circle of diameter 1.
 */
class SimTK_SIMBODY_API DecorativeCircle : public DecorativeGeometry {
public:
    explicit DecorativeCircle(Real radius=0.5);

    void setRadius(Real);
    Real getRadius() const;

    SimTK_PIMPL_DOWNCAST(DecorativeCircle, DecorativeGeometry);
};

/**
 * This defines a sphere centered at the origin. The
 * default constructor creates a sphere of diameter 1.
 */
class SimTK_SIMBODY_API DecorativeSphere : public DecorativeGeometry {
public:
    explicit DecorativeSphere(Real radius=0.5);

    void setRadius(Real);
    Real getRadius() const;

    SimTK_PIMPL_DOWNCAST(DecorativeSphere, DecorativeGeometry);
};

/**
 * This defines a rectangular solid centered at the origin and
 * aligned with the local frame axes. The default constructor creates 
 * a cube of length 1 on each side.
 */
class SimTK_SIMBODY_API DecorativeBrick : public DecorativeGeometry {
public:
    explicit DecorativeBrick(const Vec3& halfLengths = Vec3(0.5));

    void setHalfLengths(const Vec3&);
    const Vec3& getHalfLengths() const;

    SimTK_PIMPL_DOWNCAST(DecorativeBrick, DecorativeGeometry);
};

/**
 * This defines a cylinder centered on the origin and aligned in the
 * y direction. The default constructor gives it a height of 1 and
 * the base circle a diameter of 1.
 */
class SimTK_SIMBODY_API DecorativeCylinder : public DecorativeGeometry {
public:
    explicit DecorativeCylinder(Real radius=0.5, Real halfHeight=0.5);

    void setRadius(Real);
    void setHalfHeight(Real);
    Real getRadius() const;
    Real getHalfHeight() const;

    SimTK_PIMPL_DOWNCAST(DecorativeCylinder, DecorativeGeometry);
};

/**
 * This defines geometry to represent a coordinate frame. The default
 * constructor makes three perpendicular lines beginning at the 
 * origin and extending in the +x, +y, and +z directions by 1 unit.
 */
class SimTK_SIMBODY_API DecorativeFrame : public DecorativeGeometry {
public:
    explicit DecorativeFrame(Real axisLength=1);

    void setAxisLength(Real);
    Real getAxisLength() const;

    SimTK_PIMPL_DOWNCAST(DecorativeFrame, DecorativeGeometry);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_DECORATIVE_GEOMETRY_H_
