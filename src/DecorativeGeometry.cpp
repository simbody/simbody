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

#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkCylinderSource.h"
#include "vtkPolyData.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"

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
    if (&src != this) {
        if (isOwnerHandle()) delete rep;
        rep = src.rep ? src.rep->clone() : 0;
    }
    return *this;
}

DecorativeGeometry::DecorativeGeometry(const AnalyticGeometry& ag) : rep(0) {
    *this = ag.generateDecorativeGeometry(); // TODO: avoid copy of rep
}

vtkPolyData* DecorativeGeometry::getVTKPolyData() {
    return updRep().getVTKPolyData();
}

void DecorativeGeometry::setPlacement(const Transform& X_BG) {
    updRep().setPlacement(X_BG);
}
const Transform& DecorativeGeometry::getPlacement() const {
    return getRep().getPlacement();
}

void DecorativeGeometry::setColor(const Vec3& rgb) {
    return updRep().setColor(rgb);
}
void DecorativeGeometry::setOpacity(Real o) {
    return updRep().setOpacity(o);
}
void DecorativeGeometry::setResolution(Real r) {
    return updRep().setResolution(r);
}
void DecorativeGeometry::setScale(Real s) {
    return updRep().setScale(s);
}

const Vec3& DecorativeGeometry::getColor() const {
    return getRep().getColor();
}
Real DecorativeGeometry::getOpacity()  const {
    return getRep().getOpacity();
}
Real DecorativeGeometry::getResolution() const {
    return getRep().getResolution();
}
Real DecorativeGeometry::getScale() const {
    return getRep().getScale();
}

void DecorativeGeometry::setRepresentationToPoints() {
    return updRep().setRepresentationToPoints();
}
void DecorativeGeometry::setRepresentationToWireframe() {
    return updRep().setRepresentationToWireframe();
}
void DecorativeGeometry::setRepresentationToSurface() {
    return updRep().setRepresentationToSurface();
}
void DecorativeGeometry::setRepresentationToUseDefault() {
    return updRep().setRepresentationToUseDefault();
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

    ////////////////////////
    // DecorativeCylinder //
    ////////////////////////

DecorativeCylinder::DecorativeCylinder(Real radius, Real halfLength) {
    rep = new DecorativeCylinderRep(radius,halfLength);
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



/*static*/ vtkTransform* 
DecorativeGeometryRep::createVTKTransform(const Transform& X_BG, const Real& s) {
    static const Real Pi = std::acos(-1.), RadiansPerDegree = Pi/180;

    const Vec3 t = X_BG.T();
    const Vec4 r = X_BG.R().convertToAngleAxis();

    vtkTransform* xform = vtkTransform::New();  // starts out as identity
    rememberVTKObject(xform);

    xform->Translate(t[0],t[1],t[2]);
    xform->RotateWXYZ(r[0]/RadiansPerDegree,r[1],r[2],r[3]);     // angle, axis
    xform->Scale(s,s,s);
    return xform; // don't forget to Delete() this later!
}

/*static*/ vtkPolyData*
DecorativeGeometryRep::transformVTKPolyData(const Transform& X_BG, const Real& s, 
                                            vtkPolyData* in)
{
    vtkTransformPolyDataFilter* trf = vtkTransformPolyDataFilter::New();
    rememberVTKObject(trf);

    trf->SetTransform(createVTKTransform(X_BG, s)); 
    trf->SetInput(in);
    return trf->GetOutput();
}

vtkPolyData* DecorativeLineRep::createVTKPolyData() {
    assert(!"DecorativeLineRep::createVTKPolyData() NOT IMPLEMENTED YET");
    return 0;
}
vtkPolyData* DecorativeCircleRep::createVTKPolyData() {
    assert(!"DecorativeCircleRep::createVTKPolyData() NOT IMPLEMENTED YET");
    return 0;
}

vtkPolyData* DecorativeSphereRep::createVTKPolyData() {
    vtkSphereSource *sphere = vtkSphereSource::New();
    rememberVTKObject(sphere);

    sphere->SetRadius(getRadius());

    int res = DefaultResolution;
    if (getResolution() > 0.) 
        res = (int)(res*getResolution()+0.5);
    sphere->SetThetaResolution(res);
    sphere->SetPhiResolution(res);

    const Real scale = getScale() > 0. ? getScale() : 1.;

    vtkPolyData* data = transformVTKPolyData(getPlacement(), scale, 
                                             sphere->GetOutput());
    return data;
}

vtkPolyData* DecorativeBrickRep::createVTKPolyData() {
    vtkCubeSource* cube = vtkCubeSource::New();
    rememberVTKObject(cube);

    const Vec3& h = getXYZHalfLengths();

    cube->SetXLength(2*h[0]);
    cube->SetXLength(2*h[0]);
    cube->SetXLength(2*h[0]);

    // resolution is ignored -- our needs are few for rectangles!

    const Real scale = getScale() > 0. ? getScale() : 1.;

    vtkPolyData* data = transformVTKPolyData(getPlacement(), scale, 
                                             cube->GetOutput());
    return data;
}

vtkPolyData* DecorativeCylinderRep::createVTKPolyData() {
    vtkCylinderSource* cyl = vtkCylinderSource::New();
    rememberVTKObject(cyl);

    cyl->SetRadius(getRadius());
    cyl->SetHeight(2*getHalfLength());

    int res = DefaultResolution;
    if (getResolution() > 0.) 
        res = (int)(res*getResolution()+0.5);
    cyl->SetResolution(res);

    const Real scale = getScale() > 0. ? getScale() : 1.;

    vtkPolyData* data = transformVTKPolyData(getPlacement(), scale, 
                                             cyl->GetOutput());
    return data;
}

vtkPolyData* DecorativeFrameRep::createVTKPolyData() {
    assert(!"DecorativeFrameRep::createVTKPolyData() NOT IMPLEMENTED YET");
    return 0;
}

} // namespace SimTK

