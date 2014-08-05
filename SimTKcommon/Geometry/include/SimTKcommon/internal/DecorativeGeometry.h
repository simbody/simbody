#ifndef SimTK_SimTKCOMMON_DECORATIVE_GEOMETRY_H_
#define SimTK_SimTKCOMMON_DECORATIVE_GEOMETRY_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-14 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Jack Middleton, Peter Eastman, Ayman Habib                   *
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
Declarations of DecorativeGeometry and related derived classes. **/

#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/PolygonalMesh.h"

#include <cassert>


namespace SimTK {

// Some common RGB values;
extern SimTK_SimTKCOMMON_EXPORT const Vec3 Black;   ///< RGB=( 0, 0, 0)
extern SimTK_SimTKCOMMON_EXPORT const Vec3 Gray;    ///< RGB=(.5,.5,.5)
extern SimTK_SimTKCOMMON_EXPORT const Vec3 Red;     ///< RGB=( 1, 0, 0)
extern SimTK_SimTKCOMMON_EXPORT const Vec3 Green;   ///< RGB=( 0, 1, 0)
extern SimTK_SimTKCOMMON_EXPORT const Vec3 Blue;    ///< RGB=( 0, 0, 1)
extern SimTK_SimTKCOMMON_EXPORT const Vec3 Yellow;  ///< RGB=( 1, 1, 0)
extern SimTK_SimTKCOMMON_EXPORT const Vec3 Orange;  ///< RGB=( 1,.5, 0)
extern SimTK_SimTKCOMMON_EXPORT const Vec3 Magenta; ///< RGB=( 1, 0, 1)
extern SimTK_SimTKCOMMON_EXPORT const Vec3 Purple;  ///< RGB=(.5, 0,.5)
extern SimTK_SimTKCOMMON_EXPORT const Vec3 Cyan;    ///< RGB=( 0, 1, 1)
extern SimTK_SimTKCOMMON_EXPORT const Vec3 White;   ///< RGB=( 1, 1, 1)

// Drawing representations

class DecorativeGeometryImplementation;

/** This is the client-side interface to an implementation-independent
representation of "Decorations" suitable for visualization, annotation,
logging, or debugging but which cannot have any effect on the behavior of
a System or the evolution of a Study. DO NOT confuse this with real geometry 
(like contact geometry) which can represent physically meaningful objects that 
may interact and change the behavior of a System. However, all geometry objects
can generate %DecorativeGeometry for their visualization.

Why is there a %DecorativeGeometry facility at the System level at all, so far 
away from any application program? That's because for crude visualization and 
debugging purposes, the Subsystems themselves are best able to produce some 
illustrative geometry. Otherwise, you need a special purpose visualization 
tool which understands what's going on inside each subsystem. If you don't mind
taking what you get, just ask each subsystem to generate what it thinks would 
be helpful visualization. To do that, the subysystems need a way to talk about 
geometry without knowing anything about how that geometry will eventually get 
onto someone's screen. And that's why we're here!

Each %DecorativeGeometry object has its own local coordinate system and is 
defined self-consistently but independent of anything else. Clients can 
associate these with a reference frame (e.g. a body), and place the local frame
of the geometry objects on the reference frame, or at a fixed transform from 
the reference frame. That places the %DecorativeGeometry objects in a scene. We
support both 3D objects which are attached to actors in the scene, and 2D 
"screen" objects like titles which are attached to the display rather than the 
actors. The classes here deal only with the local-frame definitions of the 
geometric objects, not their placement in the scene. 

This is an abstract handle class using the PIMPL design pattern to hide the
private implementation. This is effectively an abstract class although the 
virtual function table is hidden in the private part. **/
class SimTK_SimTKCOMMON_EXPORT DecorativeGeometry {
public:
/** Default constructor creates an empty handle. **/
DecorativeGeometry() : rep(0) { }
~DecorativeGeometry();
/** Copy construction is deep; the source object will be cloned to create an
independent copy. **/
DecorativeGeometry(const DecorativeGeometry& source);
/** Copy assignment is deep; the handle will be cleared if necessary and then
the source object will be cloned to create an independent copy. **/
DecorativeGeometry& operator=(const DecorativeGeometry& source);

/** Drawing modes. **/
enum Representation {
    DrawPoints    =  1, ///< Use a cloud of points.
    DrawWireframe =  2, ///< Use a line drawing.
    DrawSurface   =  3, ///< Use a shaded surface.

    DrawDefault   = -1  ///< Let someone else decide.
};

/** By default the geometry should be placed relative to the Ground frame. If 
you want it attached to another reference frame (body), say so here. The 
geometry should be rendered with respect to the indicated body frame; however, 
the interpretation of this integer Id is left to the implementation. If you 
don't set the \a bodyId yourself it will be zero. For use in Simbody, the 
\a bodyId is interpreted as a MobilizedBodyIndex that can be mapped to a 
MobilizedBody that carries a reference frame. 

The \a bodyId is copied if you copy the %DecorativeGeometry object. If you are
copying to a different body you'll need to change the bodyId afterwards. **/
DecorativeGeometry& setBodyId(int bodyId);

/** For selection or other purposes, you may want to use this method to store
an index that can identify this particular piece of geometry. As an alternative,
or addition, see setUserRef(). In any case the \a index is simply stored with 
the object and returned if you ask for it. If you don't set it the value 
is -1. The \a index is copied if you copy the %DecorativeGeometry object. Be
sure to change it afterwards if that is not the correct index for the copy. **/
DecorativeGeometry& setIndexOnBody(int index);

/** Use this method to store an arbitrary reference pointer with this 
%DecorativeGeometry object. This is particularly useful in selection operations
where the rendering of this object has been picked by a user and you want to
map it back to some meaningful object in your model. You can also use the
setIndexOnBody() method to store some additional identifying information. 
If you don't set this pointer it will be set to zero (nullptr). This value
is stored and returned only; no interpretation is done and the pointed-to
object will not be deleted when the %DecorativeGeometry object is deleted. 

@warning The value of the \a userRef pointer is copied if you make a copy of the 
%DecorativeGeometry object. That is likely to be incorrect in many 
circumstances, depending on how you are using this value. Be sure to clear or 
change the pointer if necessary after you make a copy. **/
DecorativeGeometry& setUserRef(void* userRef);

/** This transform shifts the generated polygons with respect to this object's
local frame. Subsequent calls with other transforms simply replace the earlier
one; they do not accumulate. The default transform is identity and you can call
setTransform(Transform()) to put the transform back into its original state.
This value affects the generated polygonal data. **/
DecorativeGeometry& setTransform(const Transform& X_BG);

/** Each concrete %DecorativeGeometry object is expected to have a default 
resolution that gets the point across but is cheap to draw and hence probably 
somewhat "chunky". The resolution parameter here scales that default up or 
down. The face density in the displayed representation is roughly 
proportional to this value. 1.0 means to use the default resolution. Values 
less than 1.0 are lower resolution, and values greater than 1.0 are higher 
resolution. A value less than or equal to zero here is interpreted as an 
instruction to "use the default". **/
DecorativeGeometry& setResolution(Real);

/** Each concrete DecorativeGeometry object is expected to have a default size
around "1", whatever that means for a particular object, and most objects also
allow a user-specified size on construction. The x,y,z scale factors here are
given in the object's coordinate frame, and apply to the object as the user 
built it, or to the default if the user didn't specify a size. The default 
scaling is 1,1,1 and any value less than or equal to zero here is interpreted 
as a request to "use the default" in that direction. **/
DecorativeGeometry& setScaleFactors(const Vec3& scale);

/** Convenience method to set all three scale factors to the same value. **/
DecorativeGeometry& setScale(Real scale) {return setScaleFactors(Vec3(scale));}

/** Return the \a bodyId that was supplied to the most recent setBodyId() 
call for this %DecorativeGeometry object, or zero if that method has not been
called. Copy construction and copy assignment copy the \a bodyId. This is
intended to identify the body frame to which this geometry's placement should be
relative; with the default zero value meaning the Ground or World frame. **/
int getBodyId() const;

/** Return the \a index that was supplied to the most recent setIndexOnBody() 
call for this %DecorativeGeometry object, or -1 if that method has not been
called. Copy construction and copy assignment copy the \a index. Interpretation 
of this integer is up to the caller. **/
int getIndexOnBody() const;

/** Return the pointer value that was supplied to the most recent setUserRef()
call for this %DecorativeGeometry object, or zero (nullptr) if that method has
not been called. Copy construction and copy assignment copy the pointer. 
Interpretation of this value is up to the caller. **/
void* getUserRef() const; 

/** Return the current setting of the "resolution" factor. A return value of
-1 means "use the default". **/
Real getResolution() const;

/** Return the current value of the object's transform. If none has been set 
this will be the identity transform. Note that this transform specifies how the
polygons are placed with respect to the object's local frame. **/
const Transform& getTransform() const;

/** Return the current setting of the "scale" factors. A return value of -1 
in one of the factors means to "use the default" (which is typically 1) in
that direction. **/
const Vec3& getScaleFactors() const;

/** Request a specific color for this DecorativeGeometry object. The default 
is that the color is determined elsewhere. To explicitly request the default,
set the color to Vec3(-1). The implementation will check the 0'th element
(that is, the "R" element) and if it is less than zero will ignore the other
two elements and use the default for all three. **/
DecorativeGeometry& setColor(const Vec3& rgb); // 0-1 for each color

/** Request a level of transparency for this DecorativeGeometry. This does NOT
affect the generated geometry here. The default is that opacity is 
determined elsewhere. **/
DecorativeGeometry& setOpacity(Real);          // 0-1; default is 1 (opaque)

/** Request an adjustment to the default rendering of lines and curves. This 
does NOT affect geometry generated here; it is a request passed on to the
renderer which will probably pass it on to the hardware. A value less
than or equal to zero here is interpreted as "use the default". **/
DecorativeGeometry& setLineThickness(Real);

/** Return the color specified for this object, if any, otherwise Vec3(-1)
indicating that the default color will be used. **/
const Vec3& getColor()      const;
/** Return the opacity specified for this object. **/
Real getOpacity()    const;
/** Return the line thickness specified for this object, if any, otherwise
return -1 to indicate that the default line thickness should be used. **/
Real getLineThickness() const;
    
/** Set whether the geometry acts as a billboard, always rotating to face the 
camera. The default is typically no except for text. If you want 3D text that
moves with your model, set this to false. Here 0 means false, 1 means true,
and -1 means "use default". **/
DecorativeGeometry& setFaceCamera(int shouldFace);
/** Get whether the geometry acts as a billboard, always rotating to face the 
camera. Returns 0 for no, 1 for yes, -1 for "is using default". **/
int getFaceCamera() const;

/** Request a particular rendering representation of this DecorativeGeometry
object. The default is that the rendering representation choice is made 
elsewhere. **/
DecorativeGeometry& setRepresentation(const Representation&);

/** Returns drawing mode: -1 means "use default"; see above for others. **/
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


/** A point of interest. Note that the point's location is given relative
to the DecorativeGeometry frame so it will move if the geometry is transformed
when attached somewhere or displayed. The default constructor will put the
point at (0,0,0). **/
class SimTK_SimTKCOMMON_EXPORT DecorativePoint : public DecorativeGeometry {
public:
    explicit DecorativePoint(const Vec3& p=Vec3(0));

    // These are specific to DecorativePoint.

    DecorativePoint& setPoint(const Vec3& p);
    const Vec3& getPoint() const;

    // Retain the derived type when setting generic geometry options.
    DecorativePoint& setBodyId(int b)          {DecorativeGeometry::setBodyId(b);        return *this;}
    DecorativePoint& setIndexOnBody(int x)     {DecorativeGeometry::setIndexOnBody(x);   return *this;}
    DecorativePoint& setUserRef(void* p)       {DecorativeGeometry::setUserRef(p);       return *this;}
    DecorativePoint& setTransform(const Transform& X_BD) {DecorativeGeometry::setTransform(X_BD); return *this;}
    DecorativePoint& setResolution(Real r)     {DecorativeGeometry::setResolution(r);    return *this;}
    DecorativePoint& setScaleFactors(const Vec3& s) {DecorativeGeometry::setScaleFactors(s); return *this;}
    DecorativePoint& setColor(const Vec3& rgb) {DecorativeGeometry::setColor(rgb);       return *this;}
    DecorativePoint& setOpacity(Real o)        {DecorativeGeometry::setOpacity(o);       return *this;}
    DecorativePoint& setLineThickness(Real t)  {DecorativeGeometry::setLineThickness(t); return *this;}
    DecorativePoint& setRepresentation(const Representation& r) 
    {   DecorativeGeometry::setRepresentation(r); return *this; }

    SimTK_PIMPL_DOWNCAST(DecorativePoint, DecorativeGeometry);
private:
    class DecorativePointRep& updRep();
    const DecorativePointRep& getRep() const;
};

/** A line between two points. Note that the actual placement can be changed 
by the parent class transform & scale; here we are just generating the 
initial line in the geometry object's local frame. 

There is a default constructor for this object but it is not much
use unless followed by endpoint specifications. By default we produce
a line going from (0,0,0) to (1,1,1) just so it will show up if you
forget to set it to something meaningful. Having a default constructor
allows us to have arrays of these objects. **/
class SimTK_SimTKCOMMON_EXPORT DecorativeLine : public DecorativeGeometry {
public:
    explicit DecorativeLine(const Vec3& p1=Vec3(0), const Vec3& p2=Vec3(1)); // line between p1 and p2

    // These are specific to lines.
    DecorativeLine& setPoint1(const Vec3& p1);
    DecorativeLine& setPoint2(const Vec3& p2);
    DecorativeLine& setEndpoints(const Vec3& p1, const Vec3& p2);

    // Retain the derived type when setting generic geometry options.
    DecorativeLine& setBodyId(int b)          {DecorativeGeometry::setBodyId(b);        return *this;}
    DecorativeLine& setIndexOnBody(int x)     {DecorativeGeometry::setIndexOnBody(x);   return *this;}
    DecorativeLine& setUserRef(void* p)       {DecorativeGeometry::setUserRef(p);       return *this;}
    DecorativeLine& setTransform(const Transform& X_BD) {DecorativeGeometry::setTransform(X_BD); return *this;}
    DecorativeLine& setResolution(Real r)     {DecorativeGeometry::setResolution(r);    return *this;}
    DecorativeLine& setScaleFactors(const Vec3& s) {DecorativeGeometry::setScaleFactors(s); return *this;}
    DecorativeLine& setColor(const Vec3& rgb) {DecorativeGeometry::setColor(rgb);       return *this;}
    DecorativeLine& setOpacity(Real o)        {DecorativeGeometry::setOpacity(o);       return *this;}
    DecorativeLine& setLineThickness(Real t)  {DecorativeGeometry::setLineThickness(t); return *this;}
    DecorativeLine& setRepresentation(const Representation& r) 
    {   DecorativeGeometry::setRepresentation(r); return *this; }

    const Vec3& getPoint1() const;
    const Vec3& getPoint2() const;

    SimTK_PIMPL_DOWNCAST(DecorativeLine, DecorativeGeometry);
private:
    class DecorativeLineRep& updRep();
    const DecorativeLineRep& getRep() const;
};

/** This defines a circle in the x-y plane, centered at the origin. The
default constructor creates a circle of diameter 1. **/
class SimTK_SimTKCOMMON_EXPORT DecorativeCircle : public DecorativeGeometry {
public:
    explicit DecorativeCircle(Real radius=0.5);

    DecorativeCircle& setRadius(Real);
    Real getRadius() const;

    // Retain the derived type when setting generic geometry options.
    DecorativeCircle& setBodyId(int b)          {DecorativeGeometry::setBodyId(b);        return *this;}
    DecorativeCircle& setIndexOnBody(int x)     {DecorativeGeometry::setIndexOnBody(x);   return *this;}
    DecorativeCircle& setUserRef(void* p)       {DecorativeGeometry::setUserRef(p);       return *this;}
    DecorativeCircle& setTransform(const Transform& X_BD) {DecorativeGeometry::setTransform(X_BD); return *this;}
    DecorativeCircle& setResolution(Real r)     {DecorativeGeometry::setResolution(r);    return *this;}
    DecorativeCircle& setScaleFactors(const Vec3& s) {DecorativeGeometry::setScaleFactors(s); return *this;}
    DecorativeCircle& setColor(const Vec3& rgb) {DecorativeGeometry::setColor(rgb);       return *this;}
    DecorativeCircle& setOpacity(Real o)        {DecorativeGeometry::setOpacity(o);       return *this;}
    DecorativeCircle& setLineThickness(Real t)  {DecorativeGeometry::setLineThickness(t); return *this;}
    DecorativeCircle& setRepresentation(const Representation& r) 
    {   DecorativeGeometry::setRepresentation(r); return *this; }

    SimTK_PIMPL_DOWNCAST(DecorativeCircle, DecorativeGeometry);
private:
    class DecorativeCircleRep& updRep();
    const DecorativeCircleRep& getRep() const;
};

/** This defines a sphere centered at the origin. The default constructor 
creates a sphere of diameter 1. **/
class SimTK_SimTKCOMMON_EXPORT DecorativeSphere : public DecorativeGeometry {
public:
    explicit DecorativeSphere(Real radius=0.5);

    DecorativeSphere& setRadius(Real);
    Real getRadius() const;

    // Retain the derived type when setting generic geometry options.
    DecorativeSphere& setBodyId(int b)          {DecorativeGeometry::setBodyId(b);        return *this;}
    DecorativeSphere& setIndexOnBody(int x)     {DecorativeGeometry::setIndexOnBody(x);   return *this;}
    DecorativeSphere& setUserRef(void* p)       {DecorativeGeometry::setUserRef(p);       return *this;}
    DecorativeSphere& setTransform(const Transform& X_BD) {DecorativeGeometry::setTransform(X_BD); return *this;}
    DecorativeSphere& setResolution(Real r)     {DecorativeGeometry::setResolution(r);    return *this;}
    DecorativeSphere& setScaleFactors(const Vec3& s) {DecorativeGeometry::setScaleFactors(s); return *this;}
    DecorativeSphere& setColor(const Vec3& rgb) {DecorativeGeometry::setColor(rgb);       return *this;}
    DecorativeSphere& setOpacity(Real o)        {DecorativeGeometry::setOpacity(o);       return *this;}
    DecorativeSphere& setLineThickness(Real t)  {DecorativeGeometry::setLineThickness(t); return *this;}
    DecorativeSphere& setRepresentation(const Representation& r) 
    {   DecorativeGeometry::setRepresentation(r); return *this; }

    SimTK_PIMPL_DOWNCAST(DecorativeSphere, DecorativeGeometry);
private:
    class DecorativeSphereRep& updRep();
    const DecorativeSphereRep& getRep() const;
};

/** This defines an ellipsoidal solid centered at the origin and aligned with
the local frame axes. The default constructor creates an ellipsoid with radii 
(1/2, 1/3, 1/4) in x,y,z resp. **/
class SimTK_SimTKCOMMON_EXPORT DecorativeEllipsoid : public DecorativeGeometry {
public:
    explicit DecorativeEllipsoid(const Vec3& radii = 
        Vec3(Real(0.5),Real(1/3.),Real(0.25)));

    DecorativeEllipsoid& setRadii(const Vec3&);
    const Vec3& getRadii() const;

    // Retain the derived type when setting generic geometry options.
    DecorativeEllipsoid& setBodyId(int b)          {DecorativeGeometry::setBodyId(b);        return *this;}
    DecorativeEllipsoid& setIndexOnBody(int x)     {DecorativeGeometry::setIndexOnBody(x);   return *this;}
    DecorativeEllipsoid& setUserRef(void* p)       {DecorativeGeometry::setUserRef(p);       return *this;}
    DecorativeEllipsoid& setTransform(const Transform& X_BD) {DecorativeGeometry::setTransform(X_BD); return *this;}
    DecorativeEllipsoid& setResolution(Real r)     {DecorativeGeometry::setResolution(r);    return *this;}
    DecorativeEllipsoid& setScaleFactors(const Vec3& s) {DecorativeGeometry::setScaleFactors(s); return *this;}
    DecorativeEllipsoid& setColor(const Vec3& rgb) {DecorativeGeometry::setColor(rgb);       return *this;}
    DecorativeEllipsoid& setOpacity(Real o)        {DecorativeGeometry::setOpacity(o);       return *this;}
    DecorativeEllipsoid& setLineThickness(Real t)  {DecorativeGeometry::setLineThickness(t); return *this;}
    DecorativeEllipsoid& setRepresentation(const Representation& r) 
    {   DecorativeGeometry::setRepresentation(r); return *this; }

    SimTK_PIMPL_DOWNCAST(DecorativeEllipsoid, DecorativeGeometry);
private:
    class DecorativeEllipsoidRep& updRep();
    const DecorativeEllipsoidRep& getRep() const;
};

/** This defines a rectangular solid centered at the origin and aligned with 
the local frame axes. The default constructor creates a cube of length 1 on 
each side. **/
class SimTK_SimTKCOMMON_EXPORT DecorativeBrick : public DecorativeGeometry {
public:
    explicit DecorativeBrick(const Vec3& halfLengths = Vec3(Real(0.5)));

    DecorativeBrick& setHalfLengths(const Vec3&);
    const Vec3& getHalfLengths() const;

    // Retain the derived type when setting generic geometry options.
    DecorativeBrick& setBodyId(int b)          {DecorativeGeometry::setBodyId(b);        return *this;}
    DecorativeBrick& setIndexOnBody(int x)     {DecorativeGeometry::setIndexOnBody(x);   return *this;}
    DecorativeBrick& setUserRef(void* p)       {DecorativeGeometry::setUserRef(p);       return *this;}
    DecorativeBrick& setTransform(const Transform& X_BD) {DecorativeGeometry::setTransform(X_BD); return *this;}
    DecorativeBrick& setResolution(Real r)     {DecorativeGeometry::setResolution(r);    return *this;}
    DecorativeBrick& setScaleFactors(const Vec3& s) {DecorativeGeometry::setScaleFactors(s); return *this;}
    DecorativeBrick& setColor(const Vec3& rgb) {DecorativeGeometry::setColor(rgb);       return *this;}
    DecorativeBrick& setOpacity(Real o)        {DecorativeGeometry::setOpacity(o);       return *this;}
    DecorativeBrick& setLineThickness(Real t)  {DecorativeGeometry::setLineThickness(t); return *this;}
    DecorativeBrick& setRepresentation(const Representation& r) 
    {   DecorativeGeometry::setRepresentation(r); return *this; }

    SimTK_PIMPL_DOWNCAST(DecorativeBrick, DecorativeGeometry);
private:
    class DecorativeBrickRep& updRep();
    const DecorativeBrickRep& getRep() const;
};

/** This defines a cylinder centered on the origin and aligned in the y 
direction. The default constructor gives it a height of 1 and the base circle 
a diameter of 1. **/
class SimTK_SimTKCOMMON_EXPORT DecorativeCylinder : public DecorativeGeometry {
public:
    explicit DecorativeCylinder(Real radius=0.5, Real halfHeight=0.5);

    DecorativeCylinder& setRadius(Real);
    DecorativeCylinder& setHalfHeight(Real);
    Real getRadius() const;
    Real getHalfHeight() const;

    // Retain the derived type when setting generic geometry options.
    DecorativeCylinder& setBodyId(int b)          {DecorativeGeometry::setBodyId(b);        return *this;}
    DecorativeCylinder& setIndexOnBody(int x)     {DecorativeGeometry::setIndexOnBody(x);   return *this;}
    DecorativeCylinder& setUserRef(void* p)       {DecorativeGeometry::setUserRef(p);       return *this;}
    DecorativeCylinder& setTransform(const Transform& X_BD) {DecorativeGeometry::setTransform(X_BD); return *this;}
    DecorativeCylinder& setResolution(Real r)     {DecorativeGeometry::setResolution(r);    return *this;}
    DecorativeCylinder& setScaleFactors(const Vec3& s) {DecorativeGeometry::setScaleFactors(s); return *this;}
    DecorativeCylinder& setColor(const Vec3& rgb) {DecorativeGeometry::setColor(rgb);       return *this;}
    DecorativeCylinder& setOpacity(Real o)        {DecorativeGeometry::setOpacity(o);       return *this;}
    DecorativeCylinder& setLineThickness(Real t)  {DecorativeGeometry::setLineThickness(t); return *this;}
    DecorativeCylinder& setRepresentation(const Representation& r) 
    {   DecorativeGeometry::setRepresentation(r); return *this; }

    SimTK_PIMPL_DOWNCAST(DecorativeCylinder, DecorativeGeometry);
private:
    class DecorativeCylinderRep& updRep();
    const DecorativeCylinderRep& getRep() const;
};

/** This defines geometry to represent a coordinate frame. The default 
constructor makes three perpendicular lines beginning at the origin and 
extending in the +x, +y, and +z directions by 1 unit. **/
class SimTK_SimTKCOMMON_EXPORT DecorativeFrame : public DecorativeGeometry {
public:
    explicit DecorativeFrame(Real axisLength=1);

    DecorativeFrame& setAxisLength(Real);
    Real getAxisLength() const;

    // Retain the derived type when setting generic geometry options.
    DecorativeFrame& setBodyId(int b)          {DecorativeGeometry::setBodyId(b);        return *this;}
    DecorativeFrame& setIndexOnBody(int x)     {DecorativeGeometry::setIndexOnBody(x);   return *this;}
    DecorativeFrame& setUserRef(void* p)       {DecorativeGeometry::setUserRef(p);       return *this;}
    DecorativeFrame& setTransform(const Transform& X_BD) {DecorativeGeometry::setTransform(X_BD); return *this;}
    DecorativeFrame& setResolution(Real r)     {DecorativeGeometry::setResolution(r);    return *this;}
    DecorativeFrame& setScaleFactors(const Vec3& s) {DecorativeGeometry::setScaleFactors(s); return *this;}
    DecorativeFrame& setColor(const Vec3& rgb) {DecorativeGeometry::setColor(rgb);       return *this;}
    DecorativeFrame& setOpacity(Real o)        {DecorativeGeometry::setOpacity(o);       return *this;}
    DecorativeFrame& setLineThickness(Real t)  {DecorativeGeometry::setLineThickness(t); return *this;}
    DecorativeFrame& setRepresentation(const Representation& r) 
    {   DecorativeGeometry::setRepresentation(r); return *this; }

    SimTK_PIMPL_DOWNCAST(DecorativeFrame, DecorativeGeometry);
private:
    class DecorativeFrameRep& updRep();
    const DecorativeFrameRep& getRep() const;
};

/** This defines a text label with its base at the origin. The default 
constructor creates a blank label. **/
class SimTK_SimTKCOMMON_EXPORT DecorativeText : public DecorativeGeometry {
public:
    explicit DecorativeText(const std::string& label="");

    DecorativeText& setText(const std::string& label);
    const std::string& getText() const;

    /** By default the text is part of the scene; set this flag if you want
    it to just show up in a fixed spot on the screen instead. **/
    DecorativeText& setIsScreenText(bool isScreen);
    bool getIsScreenText() const;

    // Retain the derived type when setting generic geometry options.
    DecorativeText& setBodyId(int b)          {DecorativeGeometry::setBodyId(b);        return *this;}
    DecorativeText& setIndexOnBody(int x)     {DecorativeGeometry::setIndexOnBody(x);   return *this;}
    DecorativeText& setUserRef(void* p)       {DecorativeGeometry::setUserRef(p);       return *this;}
    DecorativeText& setTransform(const Transform& X_BD) {DecorativeGeometry::setTransform(X_BD); return *this;}
    DecorativeText& setResolution(Real r)     {DecorativeGeometry::setResolution(r);    return *this;}
    DecorativeText& setScaleFactors(const Vec3& s) {DecorativeGeometry::setScaleFactors(s); return *this;}
    DecorativeText& setColor(const Vec3& rgb) {DecorativeGeometry::setColor(rgb);       return *this;}
    DecorativeText& setOpacity(Real o)        {DecorativeGeometry::setOpacity(o);       return *this;}
    DecorativeText& setLineThickness(Real t)  {DecorativeGeometry::setLineThickness(t); return *this;}
    DecorativeText& setFaceCamera(int yn)     {DecorativeGeometry::setFaceCamera(yn);   return *this;}
    DecorativeText& setRepresentation(const Representation& r) 
    {   DecorativeGeometry::setRepresentation(r); return *this; }

    SimTK_PIMPL_DOWNCAST(DecorativeText, DecorativeGeometry);
private:
    class DecorativeTextRep& updRep();
    const DecorativeTextRep& getRep() const;
};

/** This defines a displayable mesh by referencing an already-existing
PolygonalMesh object. **/
class SimTK_SimTKCOMMON_EXPORT DecorativeMesh : public DecorativeGeometry {
public:
    explicit DecorativeMesh(const PolygonalMesh& mesh);
    const PolygonalMesh& getMesh() const;

    // Retain the derived type when setting generic geometry options.
    DecorativeMesh& setBodyId(int b)          {DecorativeGeometry::setBodyId(b);        return *this;}
    DecorativeMesh& setIndexOnBody(int x)     {DecorativeGeometry::setIndexOnBody(x);   return *this;}
    DecorativeMesh& setUserRef(void* p)       {DecorativeGeometry::setUserRef(p);       return *this;}
    DecorativeMesh& setTransform(const Transform& X_BD) {DecorativeGeometry::setTransform(X_BD); return *this;}
    DecorativeMesh& setResolution(Real r)     {DecorativeGeometry::setResolution(r);    return *this;}
    DecorativeMesh& setScaleFactors(const Vec3& s) {DecorativeGeometry::setScaleFactors(s); return *this;}
    DecorativeMesh& setColor(const Vec3& rgb) {DecorativeGeometry::setColor(rgb);       return *this;}
    DecorativeMesh& setOpacity(Real o)        {DecorativeGeometry::setOpacity(o);       return *this;}
    DecorativeMesh& setLineThickness(Real t)  {DecorativeGeometry::setLineThickness(t); return *this;}
    DecorativeMesh& setRepresentation(const Representation& r) 
    {   DecorativeGeometry::setRepresentation(r); return *this; }

    SimTK_PIMPL_DOWNCAST(DecorativeMesh, DecorativeGeometry);
private:
    class DecorativeMeshRep& updRep();
    const DecorativeMeshRep& getRep() const;
};


/** This defines a displayable mesh by referencing a file name containing the mesh. If format is not supported
  by visualizer it will be ignored. . **/
class SimTK_SimTKCOMMON_EXPORT DecorativeMeshFile : public DecorativeGeometry {
public:
    explicit DecorativeMeshFile(const std::string& meshFileName);
    const std::string& getMeshFile() const;

    // Retain the derived type when setting generic geometry options.
    DecorativeMeshFile& setBodyId(int b)          {DecorativeGeometry::setBodyId(b);        return *this;}
    DecorativeMeshFile& setIndexOnBody(int x)     {DecorativeGeometry::setIndexOnBody(x);   return *this;}
    DecorativeMeshFile& setUserRef(void* p)       {DecorativeGeometry::setUserRef(p);       return *this;}
    DecorativeMeshFile& setTransform(const Transform& X_BD) {DecorativeGeometry::setTransform(X_BD); return *this;}
    DecorativeMeshFile& setResolution(Real r)     {DecorativeGeometry::setResolution(r);    return *this;}
    DecorativeMeshFile& setScaleFactors(const Vec3& s) {DecorativeGeometry::setScaleFactors(s); return *this;}
    DecorativeMeshFile& setColor(const Vec3& rgb) {DecorativeGeometry::setColor(rgb);       return *this;}
    DecorativeMeshFile& setOpacity(Real o)        {DecorativeGeometry::setOpacity(o);       return *this;}
    DecorativeMeshFile& setLineThickness(Real t)  {DecorativeGeometry::setLineThickness(t); return *this;}
    DecorativeMeshFile& setRepresentation(const Representation& r) 
    {   DecorativeGeometry::setRepresentation(r); return *this; }

    SimTK_PIMPL_DOWNCAST(DecorativeMeshFile, DecorativeGeometry);
private:
    class DecorativeMeshFileRep& updRep();
    const DecorativeMeshFileRep& getRep() const;
};


/** This defines a single DecorativeGeometry object that is composed of a
collection of other DecorativeGeometry objects. Parameters set for the
parent object serve as defaults for the contained objects, but those objects
can override or be composed with the default values. Body id, index, and
userRef pointer (if any) for the contained objects are overridden 
unconditionally from the %Decorations object so that implementations will see
these as identical for each piece of the composite object. Reference frames are
composed; scale factors, resolution, opacity, and line thickness are composed,
with unspecified values (-1) being treated as 1. **/
class SimTK_SimTKCOMMON_EXPORT Decorations : public DecorativeGeometry {
public:
    /** Construct an empty container for DecorativeGeometry objects. **/
    Decorations();
    /** Construct a Decorations container initially consting of just a single
    DecorativeGeometry object. **/
    explicit Decorations(const DecorativeGeometry& decoration);
    /** Add a DecorativeGeometry object to this collection. **/
    Decorations& addDecoration(const DecorativeGeometry& decoration);
    /** Add a DecorativeGeometry object to this collection and place it
    relative to the Decorations frame. **/
    Decorations& addDecoration(const Transform& placement,
                               const DecorativeGeometry& decoration);
    /** Determine how many DecorativeGeometry objects are in this 
    collection. **/
    int getNumDecorations() const;
    /** Get access to one of the DecorativeGeometry objects in this 
    collection. **/
    const DecorativeGeometry& getDecoration(int i) const;

    // Retain the derived type when setting generic geometry options.
    Decorations& setBodyId(int b)          {DecorativeGeometry::setBodyId(b);        return *this;}
    Decorations& setIndexOnBody(int x)     {DecorativeGeometry::setIndexOnBody(x);   return *this;}
    Decorations& setUserRef(void* p)       {DecorativeGeometry::setUserRef(p);       return *this;}
    Decorations& setTransform(const Transform& X_BD) {DecorativeGeometry::setTransform(X_BD); return *this;}
    Decorations& setResolution(Real r)     {DecorativeGeometry::setResolution(r);    return *this;}
    Decorations& setScaleFactors(const Vec3& s) {DecorativeGeometry::setScaleFactors(s); return *this;}
    Decorations& setColor(const Vec3& rgb) {DecorativeGeometry::setColor(rgb);       return *this;}
    Decorations& setOpacity(Real o)        {DecorativeGeometry::setOpacity(o);       return *this;}
    Decorations& setLineThickness(Real t)  {DecorativeGeometry::setLineThickness(t); return *this;}
    Decorations& setRepresentation(const Representation& r) 
    {   DecorativeGeometry::setRepresentation(r); return *this; }


    SimTK_PIMPL_DOWNCAST(Decorations, DecorativeGeometry);
private:
    class DecorationsRep& updRep();
    const DecorationsRep& getRep() const;
};

/** Use this abstract class to connect your implementation of decorative 
geometry to the implementation-independent classes above. **/
class SimTK_SimTKCOMMON_EXPORT DecorativeGeometryImplementation {
public:
    virtual ~DecorativeGeometryImplementation() { }
    virtual void implementPointGeometry(    const DecorativePoint&)    = 0;
    virtual void implementLineGeometry(     const DecorativeLine&)     = 0;
    virtual void implementBrickGeometry(    const DecorativeBrick&)    = 0;
    virtual void implementCylinderGeometry( const DecorativeCylinder&) = 0;
    virtual void implementCircleGeometry(   const DecorativeCircle&)   = 0; 
    virtual void implementSphereGeometry(   const DecorativeSphere&)   = 0;
    virtual void implementEllipsoidGeometry(const DecorativeEllipsoid&)= 0;
    virtual void implementFrameGeometry(    const DecorativeFrame&)    = 0;
    virtual void implementTextGeometry(     const DecorativeText&)     = 0;
    virtual void implementMeshGeometry(     const DecorativeMesh&)     = 0;
    virtual void implementMeshFileGeometry(     const DecorativeMeshFile&)    =0;
};

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_DECORATIVE_GEOMETRY_H_
