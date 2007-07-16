#ifndef SimTK_DECORATIVE_GEOMETRY_H_
#define SimTK_DECORATIVE_GEOMETRY_H_

/* Portions copyright (c) 2005-7 Stanford University and Michael Sherman.
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
 * This is the client-side interface to an implementation-independent
 * representation of "Decorations" suitable for visualization, annotation,
 * logging, or debugging but which cannot have any effect on the behavior of
 * a System or the evolution of a Study. 
 * DO NOT confuse this with AnalyticGeometry which can represent
 * physically meaningful objects that may interact and change the behavior
 * of a System. However, an AnalyticGeometry object is likely to generate a
 * corresponding Decoration for visualization.
 *
 * Why is there a DecorativeGeometry facility at the System level at all, so far away from
 * any application program? That's because for crude visualization and debugging
 * purposes, the Subsystems themselves are best able to produce some illustrative
 * geometry. Otherwise, you need a special purpose visualization tool which understands
 * what's going on inside each subsystem. If you don't mind taking what you get, just
 * ask each subsystem to generate what it thinks would be helpful visualization. To do
 * that, the subysystems need a way to talk about geometry without knowing anything
 * about how that geometry will eventually get onto someone's screen. And that's why
 * we're here!
 *
 * Each DecorativeGeometry object has its own local coordinate system
 * and is defined self-consistently
 * but independent of anything else. Clients can associate these with a reference
 * frame (e.g. a body), and place the local frame of the geometry objects on the 
 * reference frame, or at a fixed transform from the reference frame. That places
 * the DecorativeGeometry objects in a scene. We support both 3D objects
 * which are attached to actors in the scene, and 2D "screen"
 * objects like titles which are attached to the display rather than the actors.
 * The classes here deal only with the local-frame definitions of the geometric
 * objects, not their placement in the scene.
 */

/* TODO: NEW TEXT -- NOT IMPLEMENTED YET (just thinking out loud)
 *
 * Each kind of Decoration defines a parameterized object, whose appearance
 * is specified once its parameters are known. Parameters will be scalars,
 * points, direction vectors, and transforms. These involve dependencies on
 * a fixed set of reference frames: the ground body frame, the moving body frames,
 * and the screen or viewing frame. The body frames are defined by the System;
 * while the viewing frame is defined by whatever application is communicating
 * to a user. The World frame of the viewer and the Ground frame of the System
 * are the same. At run time a System can always report the configuration of its
 * moving bodies with respect to Ground, and hence these can be placed in the
 * scene with respect to the World frame. A text label, for example, could be tethered
 * to a point on a moving body while defined to remain parallel to the viewing plane.
 *
 * In addition to their geometric representations, Decorations defined here can 
 * state preferences, if any, for a variety of viewing options such as color,
 * drawing style, and transparency. If left unspecified, the viewer will have to
 * determine how to set these options; if specified the viewer may want to use
 * the specified options but can override them if it has good reason to do so.
 * Many different implementations of Decorations are possible. The generic 
 * implementation is defined below as an abstract class; each viewer will use
 * a particular implementation that defines a concrete implementation for each of the virtual
 * methods of the abstraction.
 *
 * There are two Stages associated with a Decoration. When is the Decoration
 * known to exist (including its parameterization), and when are values for its 
 * parameters calculatable. For example, a permanent body-fixed sphere is known
 * to exist a the Topological stage (i.e., as soon as the body exists). Its 
 * configuration in the World frame, however, cannot be computed until we
 * know the body's configuration, that is, Stage::Position.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"

#include <cassert>


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

// Drawing representations

class DecorativeGeometryImplementation;

/**
 * This is an abstract handle class using the PIMPL design pattern to hide the private
 * implementation. This is effectively an abstract class although the virtual
 * function table is hidden in the private part.
 */
class SimTK_SimTKCOMMON_EXPORT DecorativeGeometry {
public:
    DecorativeGeometry() : rep(0) { }
    ~DecorativeGeometry();
    DecorativeGeometry(const DecorativeGeometry&);
    DecorativeGeometry& operator=(const DecorativeGeometry&);

    // Drawing modes.
    enum Representation {
        DrawPoints    =  1,
        DrawWireframe =  2,
        DrawSurface   =  3,

        DrawDefault   = -1
    };

    /// Implicit conversion
    DecorativeGeometry(const AnalyticGeometry&);

    /// By default the geometry will be placed on ground. If you want it attached to
    /// another reference frame (body), say so here. The geometry should be rendered
    /// with respect to the indicated body frame; however, the interpretation of this
    /// integer Id is left to the implementation.
    DecorativeGeometry& setBodyId(int);

    /// This transform shifts the generated polygons with respect to this object's
    /// local frame. Subsequent calls with other transforms simply replace the earlier
    /// one; they do not accumulate. The default transform is identity and you can 
    /// call setTransform(Transform()) to put the transform back into its original state.
    /// This value affects the generated polygonal data.
    DecorativeGeometry& setTransform(const Transform& X_BG);

    /// Each concrete DecorativeGeometry object is expected to have a default resolution
    /// that gets the point across but is cheap to draw and hence probably somewhat "chunky".
    /// The resolution parameter here scales that default up or down. A value less than
    /// or equal to zero here is interpreted as an instruction to "use the default".
    /// This value affects the generated polygonal data.
    DecorativeGeometry& setResolution(Real);

    /// Each concrete DecorativeGeometry object is expected to have a default size
    /// around "1", whatever that means for a particular object, and most objects also
    /// allow a user-specified size on construction. The scale factor here applies to
    /// the object as the user built it, or to the default if the user didn't specify
    /// a size. The default scaling is 1, and any value less than or equal to zero
    /// here is interpreted as a request to "use the default".
    /// This value affects the generated polygonal data.
    DecorativeGeometry& setScale(Real);

    /// Return the body to which this geometry is attached. The geometry's placement is
    /// interpreted relative to the body's frame.
    int getBodyId() const;

    /// Return the current setting of the "resolution" factor. A return value of -1
    /// means "use the default".
    Real getResolution() const;

    /// Return the current value of the object's transform. If none has been set this
    /// will be the identity transform. Note that this transform specifies how the
    /// polygons are placed with respect to the object's local frame.
    const Transform& getTransform() const;

    /// Return the current setting of the "scale" factor. A return value of -1 
    /// means "use the default" (which is typically 1).
    Real getScale() const;

    /// Request a specific color for this DecorativeGeometry object. This does NOT
    /// affect the generated geometry here. The default is that the color is
    /// determined elsewhere.
    DecorativeGeometry& setColor(const Vec3& rgb); // 0-1 for each color

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

    /// Request a particular rendering representation of this DecorativeGeometry
    /// object. The default is that the rendering representation choice is made elsewhere.
    DecorativeGeometry& setRepresentation(const Representation&);

    /// returns drawing mode: -1 means "use default"; see above for others
    Representation getRepresentation() const;

    void implementGeometry(DecorativeGeometryImplementation&) const;

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
class SimTK_SimTKCOMMON_EXPORT DecorativeLine : public DecorativeGeometry {
public:
    explicit DecorativeLine(const Vec3& p1=Vec3(0), const Vec3& p2=Vec3(1)); // line between p1 and p2

    // Retain the derived type when setting generic geometry options.
    DecorativeLine& setBodyId(int b)          {DecorativeGeometry::setBodyId(b);        return *this;}
    DecorativeLine& setTransform(const Transform& X_BD) {DecorativeGeometry::setTransform(X_BD); return *this;}
    DecorativeLine& setResolution(Real r)     {DecorativeGeometry::setResolution(r);    return *this;}
    DecorativeLine& setScale(Real s)          {DecorativeGeometry::setScale(s);         return *this;}
    DecorativeLine& setColor(const Vec3& rgb) {DecorativeGeometry::setColor(rgb);       return *this;}
    DecorativeLine& setOpacity(Real o)        {DecorativeGeometry::setOpacity(o);       return *this;}
    DecorativeLine& setLineThickness(Real t)  {DecorativeGeometry::setLineThickness(t); return *this;}
    DecorativeLine& setRepresentation(const Representation& r) {
        DecorativeGeometry::setRepresentation(r); return *this;
    }

    // These are specific to lines.
    DecorativeLine& setPoint1(const Vec3& p1);
    DecorativeLine& setPoint2(const Vec3& p2);
    DecorativeLine& setEndpoints(const Vec3& p1, const Vec3& p2);

    const Vec3& getPoint1() const;
    const Vec3& getPoint2() const;

    SimTK_PIMPL_DOWNCAST(DecorativeLine, DecorativeGeometry);
private:
    class DecorativeLineRep& updRep();
    const DecorativeLineRep& getRep() const;
};

/**
 * This defines a circle in the x-y plane, centered at the origin. The
 * default constructor creates a circle of diameter 1.
 */
class SimTK_SimTKCOMMON_EXPORT DecorativeCircle : public DecorativeGeometry {
public:
    explicit DecorativeCircle(Real radius=0.5);

    void setRadius(Real);
    Real getRadius() const;

    SimTK_PIMPL_DOWNCAST(DecorativeCircle, DecorativeGeometry);
private:
    class DecorativeCircleRep& updRep();
    const DecorativeCircleRep& getRep() const;
};

/**
 * This defines a sphere centered at the origin. The
 * default constructor creates a sphere of diameter 1.
 */
class SimTK_SimTKCOMMON_EXPORT DecorativeSphere : public DecorativeGeometry {
public:
    explicit DecorativeSphere(Real radius=0.5);

    void setRadius(Real);
    Real getRadius() const;

    SimTK_PIMPL_DOWNCAST(DecorativeSphere, DecorativeGeometry);
private:
    class DecorativeSphereRep& updRep();
    const DecorativeSphereRep& getRep() const;
};

/**
 * This defines an ellipsoidal solid centered at the origin and
 * aligned with the local frame axes. The default constructor creates 
 * an ellipsoid with radii (1/2, 1/3, 1/4) in x,y,z resp.
 */
class SimTK_SimTKCOMMON_EXPORT DecorativeEllipsoid : public DecorativeGeometry {
public:
    explicit DecorativeEllipsoid(const Vec3& radii = Vec3(0.5,1/3.,0.25));

    void setRadii(const Vec3&);
    const Vec3& getRadii() const;

    SimTK_PIMPL_DOWNCAST(DecorativeEllipsoid, DecorativeGeometry);
private:
    class DecorativeEllipsoidRep& updRep();
    const DecorativeEllipsoidRep& getRep() const;
};

/**
 * This defines a rectangular solid centered at the origin and
 * aligned with the local frame axes. The default constructor creates 
 * a cube of length 1 on each side.
 */
class SimTK_SimTKCOMMON_EXPORT DecorativeBrick : public DecorativeGeometry {
public:
    explicit DecorativeBrick(const Vec3& halfLengths = Vec3(0.5));

    void setHalfLengths(const Vec3&);
    const Vec3& getHalfLengths() const;

    SimTK_PIMPL_DOWNCAST(DecorativeBrick, DecorativeGeometry);
private:
    class DecorativeBrickRep& updRep();
    const DecorativeBrickRep& getRep() const;
};

/**
 * This defines a cylinder centered on the origin and aligned in the
 * y direction. The default constructor gives it a height of 1 and
 * the base circle a diameter of 1.
 */
class SimTK_SimTKCOMMON_EXPORT DecorativeCylinder : public DecorativeGeometry {
public:
    explicit DecorativeCylinder(Real radius=0.5, Real halfHeight=0.5);

    void setRadius(Real);
    void setHalfHeight(Real);
    Real getRadius() const;
    Real getHalfHeight() const;

    SimTK_PIMPL_DOWNCAST(DecorativeCylinder, DecorativeGeometry);
private:
    class DecorativeCylinderRep& updRep();
    const DecorativeCylinderRep& getRep() const;
};

/**
 * This defines geometry to represent a coordinate frame. The default
 * constructor makes three perpendicular lines beginning at the 
 * origin and extending in the +x, +y, and +z directions by 1 unit.
 */
class SimTK_SimTKCOMMON_EXPORT DecorativeFrame : public DecorativeGeometry {
public:
    explicit DecorativeFrame(Real axisLength=1);

    void setAxisLength(Real);
    Real getAxisLength() const;

    SimTK_PIMPL_DOWNCAST(DecorativeFrame, DecorativeGeometry);
private:
    class DecorativeFrameRep& updRep();
    const DecorativeFrameRep& getRep() const;
};


/**
 * Use this abstract class to connect your implementation of decorative geometry
 * to the implementation-independent classes above.
 */
class SimTK_SimTKCOMMON_EXPORT DecorativeGeometryImplementation {
public:
    virtual ~DecorativeGeometryImplementation() { }
    virtual void implementLineGeometry(     const DecorativeLine&)     = 0;
    virtual void implementBrickGeometry(    const DecorativeBrick&)    = 0;
    virtual void implementCylinderGeometry( const DecorativeCylinder&) = 0;
    virtual void implementCircleGeometry(   const DecorativeCircle&)   = 0; 
    virtual void implementSphereGeometry(   const DecorativeSphere&)   = 0;
    virtual void implementEllipsoidGeometry(const DecorativeEllipsoid&)= 0;
    virtual void implementFrameGeometry(    const DecorativeFrame&)    = 0;

    // TODO: wrappers for binary compatibility
};

} // namespace SimTK

#endif // SimTK_DECORATIVE_GEOMETRY_H_
