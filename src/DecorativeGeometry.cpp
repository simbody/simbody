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
    if (src.rep)
        rep = src.rep->clone();
}

DecorativeGeometry& DecorativeGeometry::operator=(const DecorativeGeometry& src) {
    if (&src == this) return *this;
    if (isOwnerHandle()) delete rep;
    rep = src.rep ? src.rep->clone() : 0;
    return *this;
}

DecorativeGeometry::DecorativeGeometry(const AnalyticGeometry& ag) : rep(0) {
    *this = ag.generateDecorativeGeometry(); // TODO: avoid copy of rep
}

vtkPolyData* DecorativeGeometry::createVTKPolyData() const {
    return getRep().createVTKPolyData();
}

void DecorativeGeometry::setPlacement(const Transform& X_BG) {
    updRep().setPlacement(X_BG);
}
const Transform& DecorativeGeometry::getPlacement() const {
    return getRep().getPlacement();
}

    ////////////////////
    // DecorativeLine //
    ////////////////////

DecorativeLine::DecorativeLine(Real length) {
    rep = new DecorativeLineRep(length);
    rep->setMyHandle(*this);
}

    //////////////////////
    // DecorativeCircle //
    //////////////////////

DecorativeCircle::DecorativeCircle(Real radius) {
    rep = new DecorativeCircleRep(radius);
    rep->setMyHandle(*this);
}

    //////////////////////
    // DecorativeSphere //
    //////////////////////

DecorativeSphere::DecorativeSphere(Real radius) {
    rep = new DecorativeSphereRep(radius);
    rep->setMyHandle(*this);
}

    /////////////////////
    // DecorativeBrick //
    /////////////////////

DecorativeBrick::DecorativeBrick(const Vec3& xyzLengths) {
    rep = new DecorativeBrickRep(xyzLengths);
    rep->setMyHandle(*this);
}

    /////////////////////
    // DecorativeFrame //
    /////////////////////

DecorativeFrame::DecorativeFrame(Real axisLength) {
    rep = new DecorativeFrameRep(axisLength);
    rep->setMyHandle(*this);
}

    /////////////////////////////
    // VTK PolyData generation //
    /////////////////////////////

vtkPolyData* DecorativeLineRep::createVTKPolyData() const {
    assert(!"DecorativeLineRep::createVTKPolyData() NOT IMPLEMENTED YET");
    return 0;
}
vtkPolyData* DecorativeCircleRep::createVTKPolyData() const {
    assert(!"DecorativeCircleRep::createVTKPolyData() NOT IMPLEMENTED YET");
    return 0;
}
vtkPolyData* DecorativeSphereRep::createVTKPolyData() const {
    assert(!"DecorativeSphereRep::createVTKPolyData() NOT IMPLEMENTED YET");
    return 0;
}
vtkPolyData* DecorativeBrickRep::createVTKPolyData() const {
    assert(!"DecorativeBrickRep::createVTKPolyData() NOT IMPLEMENTED YET");
    return 0;
}
vtkPolyData* DecorativeFrameRep::createVTKPolyData() const {
    assert(!"DecorativeFrameRep::createVTKPolyData() NOT IMPLEMENTED YET");
    return 0;
}

} // namespace SimTK

