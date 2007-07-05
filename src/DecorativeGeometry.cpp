/* Portions copyright (c) 2006-7 Stanford University and Michael Sherman.
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

#include "simbody/internal/common.h"
#include "simbody/internal/AnalyticGeometry.h"
#include "simbody/internal/DecorativeGeometry.h"

#include "DecorativeGeometryRep.h"

#include <cmath>

namespace SimTK {


    ////////////////////////
    // DecorativeGeometry //
    ////////////////////////



// This is an owner handle if there is no rep or if the rep points back
// to this handle.
bool DecorativeGeometry::isOwnerHandle() const {
    return rep==0 || rep->myHandle == this;
}
bool DecorativeGeometry::isEmptyHandle() const {return rep==0;}

DecorativeGeometry::~DecorativeGeometry() {
    if (isOwnerHandle())
        delete rep;
    rep = 0;
}

DecorativeGeometry::DecorativeGeometry(const DecorativeGeometry& src) : rep(0) {
    if (src.rep) {
        rep = src.rep->clone();
        rep->setMyHandle(*this);
    }
}

DecorativeGeometry& DecorativeGeometry::operator=(const DecorativeGeometry& src) {
    if (&src != this) {
        if (isOwnerHandle()) delete rep;
        rep = 0;
        if (src.rep) {
            rep = src.rep->clone();
            rep->setMyHandle(*this);
        }
    }
    return *this;
}

DecorativeGeometry::DecorativeGeometry(const AnalyticGeometry& ag) : rep(0) {
    *this = ag.generateDecorativeGeometry(); // TODO: avoid copy of rep
}

DecorativeGeometry& DecorativeGeometry::setBodyId(MobilizedBodyId b) {updRep().setBodyId(b);return *this;}
MobilizedBodyId DecorativeGeometry::getBodyId() const {return getRep().getBodyId();}

DecorativeGeometry& DecorativeGeometry::setTransform(const Transform& X_BD) {updRep().setTransform(X_BD);return *this;}
const Transform& DecorativeGeometry::getTransform() const    {return getRep().getTransform();}

DecorativeGeometry& DecorativeGeometry::setResolution(Real r) {updRep().setResolution(r);return *this;}
Real DecorativeGeometry::getResolution() const {return getRep().getResolution();}

DecorativeGeometry& DecorativeGeometry::setScale(Real s) {updRep().setScale(s);return *this;}
Real DecorativeGeometry::getScale() const {return getRep().getScale();}

DecorativeGeometry& DecorativeGeometry::setColor(const Vec3& rgb) {updRep().setColor(rgb);return *this;}
const Vec3& DecorativeGeometry::getColor() const   {return getRep().getColor();}

DecorativeGeometry& DecorativeGeometry::setOpacity(Real o)  {updRep().setOpacity(o);return *this;}
Real DecorativeGeometry::getOpacity()  const {return getRep().getOpacity();}

DecorativeGeometry& DecorativeGeometry::setLineThickness(Real t) {updRep().setLineThickness(t);return *this;}
Real DecorativeGeometry::getLineThickness() const {return getRep().getLineThickness();}

DecorativeGeometry& DecorativeGeometry::setRepresentation(const DecorativeGeometry::Representation& r) {
    updRep().setRepresentation(r);return *this;
}

DecorativeGeometry::Representation
DecorativeGeometry::getRepresentation() const {return getRep().getRepresentation();}

void DecorativeGeometry::implementGeometry(DecorativeGeometryImplementation& geometry) const
{
    getRep().implementGeometry(geometry);
}


    ////////////////////
    // DecorativeLine //
    ////////////////////

/*static*/ bool 
DecorativeLine::isInstanceOf(const DecorativeGeometry& s) {
    return DecorativeLineRep::isA(s.getRep());
}
/*static*/ const DecorativeLine&
DecorativeLine::downcast(const DecorativeGeometry& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const DecorativeLine&>(s);
}
/*static*/ DecorativeLine&
DecorativeLine::updDowncast(DecorativeGeometry& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<DecorativeLine&>(s);
}

const DecorativeLineRep& 
DecorativeLine::getRep() const {
    return dynamic_cast<const DecorativeLineRep&>(*rep);
}
DecorativeLineRep&       
DecorativeLine::updRep() {
    return dynamic_cast<DecorativeLineRep&>(*rep);
}

DecorativeLine::DecorativeLine(const Vec3& p1, const Vec3& p2) {
    rep = new DecorativeLineRep(p1,p2);
    rep->setMyHandle(*this);
}
DecorativeLine& DecorativeLine::setPoint1(const Vec3& p1) {
    DecorativeLineRep::downcast(*rep).setPoint1(p1); return *this;
}
DecorativeLine& DecorativeLine::setPoint2(const Vec3& p2) {
    DecorativeLineRep::downcast(*rep).setPoint2(p2); return *this;
}
DecorativeLine& DecorativeLine::setEndpoints(const Vec3& p1, const Vec3& p2) {
    DecorativeLineRep::downcast(*rep).setEndpoints(p1,p2); return *this;
}
const Vec3& DecorativeLine::getPoint1() const {
    return DecorativeLineRep::downcast(*rep).getPoint1();
}
const Vec3& DecorativeLine::getPoint2() const {
    return DecorativeLineRep::downcast(*rep).getPoint2();
}


    //////////////////////
    // DecorativeCircle //
    //////////////////////

DecorativeCircle::DecorativeCircle(Real radius) {
    rep = new DecorativeCircleRep(radius);
    rep->setMyHandle(*this);
}
void DecorativeCircle::setRadius(Real r) {
    DecorativeCircleRep::downcast(*rep).setRadius(r);
}
Real DecorativeCircle::getRadius() const {
    return DecorativeCircleRep::downcast(*rep).getRadius();
}


    //////////////////////
    // DecorativeSphere //
    //////////////////////

DecorativeSphere::DecorativeSphere(Real radius) {
    rep = new DecorativeSphereRep(radius);
    rep->setMyHandle(*this);
}

void DecorativeSphere::setRadius(Real r) {
    DecorativeSphereRep::downcast(*rep).setRadius(r);
}
Real DecorativeSphere::getRadius() const {
    return DecorativeSphereRep::downcast(*rep).getRadius();
}

    /////////////////////
    // DecorativeBrick //
    /////////////////////

DecorativeBrick::DecorativeBrick(const Vec3& xyzHalfLengths) {
    rep = new DecorativeBrickRep(xyzHalfLengths);
    rep->setMyHandle(*this);
}

void DecorativeBrick::setHalfLengths(const Vec3& xyzHalfLengths) {
    DecorativeBrickRep::downcast(*rep).setHalfLengths(xyzHalfLengths);
}
const Vec3& DecorativeBrick::getHalfLengths() const {
    return DecorativeBrickRep::downcast(*rep).getHalfLengths();
}

    ////////////////////////
    // DecorativeCylinder //
    ////////////////////////

DecorativeCylinder::DecorativeCylinder(Real radius, Real halfHeight) {
    rep = new DecorativeCylinderRep(radius,halfHeight);
    rep->setMyHandle(*this);
}

void DecorativeCylinder::setRadius(Real r) {
    DecorativeCylinderRep::downcast(*rep).setRadius(r);
}
void DecorativeCylinder::setHalfHeight(Real r) {
    DecorativeCylinderRep::downcast(*rep).setHalfHeight(r);
}

Real DecorativeCylinder::getRadius() const {
    return DecorativeCylinderRep::downcast(*rep).getRadius();
}
Real DecorativeCylinder::getHalfHeight() const {
    return DecorativeCylinderRep::downcast(*rep).getHalfHeight();
}

    /////////////////////
    // DecorativeFrame //
    /////////////////////

DecorativeFrame::DecorativeFrame(Real axisLength) {
    rep = new DecorativeFrameRep(axisLength);
    rep->setMyHandle(*this);
}

void DecorativeFrame::setAxisLength(Real l) {
    DecorativeFrameRep::downcast(*rep).setAxisLength(l);
}
Real DecorativeFrame::getAxisLength() const {
    return DecorativeFrameRep::downcast(*rep).getAxisLength();
}

} // namespace SimTK

