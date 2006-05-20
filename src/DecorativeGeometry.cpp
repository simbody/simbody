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

#include "vtkLineSource.h"
#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkCylinderSource.h"
#include "vtkPolyData.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkAppendPolyData.h"
#include "vtkVectorText.h"

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

vtkPolyData* DecorativeGeometry::updVTKPolyData() {
    return updRep().updVTKPolyData();
}

DecorativeGeometry& DecorativeGeometry::setResolution(Real r) {updRep().setResolution(r);return *this;}
Real DecorativeGeometry::getResolution() const {return getRep().getResolution();}

DecorativeGeometry& DecorativeGeometry::setPlacement(const Transform& X_BG) {updRep().setPlacement(X_BG);return *this;}
const Transform& DecorativeGeometry::getPlacement() const    {return getRep().getPlacement();}

DecorativeGeometry& DecorativeGeometry::setScale(Real s) {updRep().setScale(s);return *this;}
Real DecorativeGeometry::getScale() const {return getRep().getScale();}

DecorativeGeometry& DecorativeGeometry::setColor(const Vec3& rgb) {updRep().setColor(rgb);return *this;}
const Vec3& DecorativeGeometry::getColor() const   {return getRep().getColor();}

DecorativeGeometry& DecorativeGeometry::setOpacity(Real o)  {updRep().setOpacity(o);return *this;}
Real DecorativeGeometry::getOpacity()  const {return getRep().getOpacity();}

DecorativeGeometry& DecorativeGeometry::setLineThickness(Real t) {updRep().setLineThickness(t);return *this;}
Real DecorativeGeometry::getLineThickness() const {return getRep().getLineThickness();}

DecorativeGeometry& DecorativeGeometry::setRepresentationToPoints()     {updRep().setRepresentationToPoints();return *this;}
DecorativeGeometry& DecorativeGeometry::setRepresentationToWireframe()  {updRep().setRepresentationToWireframe();return *this;}
DecorativeGeometry& DecorativeGeometry::setRepresentationToSurface()    {updRep().setRepresentationToSurface();return *this;}
DecorativeGeometry& DecorativeGeometry::setRepresentationToUseDefault() {updRep().setRepresentationToUseDefault();return *this;}
int  DecorativeGeometry::getRepresentation() const       {return getRep().getRepresentation();}

    ////////////////////
    // DecorativeLine //
    ////////////////////

DecorativeLine::DecorativeLine(const Vec3& p1, const Vec3& p2) {
    rep = new DecorativeLineRep(p1,p2);
    rep->setMyHandle(*this);
}
void DecorativeLine::setPoint1(const Vec3& p1) {
    DecorativeLineRep::downcast(*rep).setPoint1(p1);
}
void DecorativeLine::setPoint2(const Vec3& p2) {
    DecorativeLineRep::downcast(*rep).setPoint2(p2);
}
void DecorativeLine::setEndpoints(const Vec3& p1, const Vec3& p2) {
    DecorativeLineRep::downcast(*rep).setEndpoints(p1,p2);
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

    /////////////////////////////
    // VTK PolyData generation //
    /////////////////////////////



vtkTransform* 
DecorativeGeometryRep::createVTKTransform(const Transform& X_BG, const Real& s) {
    static const Real Pi = std::acos(-1.), RadiansPerDegree = Pi/180;

    const Vec3 t = X_BG.T();
    const Vec4 r = X_BG.R().convertToAngleAxis();

    vtkTransform* xform = vtkTransform::New();  // starts out as identity
    rememberVTKObject(xform);

    xform->Translate(t[0],t[1],t[2]);
    xform->RotateWXYZ(r[0]/RadiansPerDegree,r[1],r[2],r[3]);     // angle, axis
    xform->Scale(s,s,s);
    return xform;
}

vtkPolyData*
DecorativeGeometryRep::transformVTKPolyData(const Transform& X_BG, const Real& s, 
                                            vtkPolyData* in)
{
    // Careful -- the poly data filter has to be allocated last, so allocate
    // the transform first.
    vtkTransform* tr = createVTKTransform(X_BG, s);

    vtkTransformPolyDataFilter* trf = vtkTransformPolyDataFilter::New();
    rememberVTKObject(trf);

    trf->SetTransform(tr); 
    trf->SetInput(in);
    return trf->GetOutput();
}

void DecorativeLineRep::createVTKPolyData() {
    vtkLineSource* line = vtkLineSource::New(); 
    rememberVTKObject(line);

    const Vec3& p1 = getPoint1();
    const Vec3& p2 = getPoint2();

    line->SetPoint1(p1[0],p1[1],p1[2]); line->SetPoint2(p2[0],p2[1],p2[2]);

    const Real scale = getScale() > 0. ? getScale() : 1.;
    vtkPolyData* data = transformVTKPolyData(getPlacement(), scale, 
                                             line->GetOutput());
}
void DecorativeCircleRep::createVTKPolyData() {
    assert(!"DecorativeCircleRep::createVTKPolyData() NOT IMPLEMENTED YET");
}

void DecorativeSphereRep::createVTKPolyData() {
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
}

void DecorativeBrickRep::createVTKPolyData() {
    vtkCubeSource* cube = vtkCubeSource::New();
    rememberVTKObject(cube);

    const Vec3& h = getHalfLengths();

    cube->SetXLength(2*h[0]);
    cube->SetYLength(2*h[1]);
    cube->SetZLength(2*h[2]);

    // resolution is ignored -- our needs are few for rectangles!

    const Real scale = getScale() > 0. ? getScale() : 1.;

    vtkPolyData* data = transformVTKPolyData(getPlacement(), scale, 
                                             cube->GetOutput());
}

void DecorativeCylinderRep::createVTKPolyData() {
    vtkCylinderSource* cyl = vtkCylinderSource::New();
    rememberVTKObject(cyl);

    cyl->SetRadius(getRadius());
    cyl->SetHeight(2*getHalfHeight());

    int res = DefaultResolution;
    if (getResolution() > 0.) 
        res = (int)(res*getResolution()+0.5);
    cyl->SetResolution(res);

    const Real scale = getScale() > 0. ? getScale() : 1.;
    vtkPolyData* data = transformVTKPolyData(getPlacement(), scale, 
                                             cyl->GetOutput());
}

void DecorativeFrameRep::createVTKPolyData() {
    vtkAppendPolyData* app = vtkAppendPolyData::New();
    rememberVTKObject(app);

    const Real& length = getAxisLength();

    vtkLineSource* line = vtkLineSource::New(); 
    rememberVTKObject(line);
    line->SetPoint1(0,0,0); line->SetPoint2(length,0,0);
    app->AddInput(line->GetOutput());

    line = vtkLineSource::New(); 
    rememberVTKObject(line);    
    line->SetPoint1(0,0,0); line->SetPoint2(0,length,0);
    app->AddInput(line->GetOutput());

    line = vtkLineSource::New(); 
    rememberVTKObject(line); 
    line->SetPoint1(0,0,0); line->SetPoint2(0,0,length);
    app->AddInput(line->GetOutput());

    vtkVectorText* xtext = vtkVectorText::New();
    rememberVTKObject(xtext);
    xtext->SetText("x"); // default size is around 1
    vtkPolyData* label = transformVTKPolyData(Transform(Vec3(length,-0.05*length,0)), length*0.2, 
                                              xtext->GetOutput());
    app->AddInput(label);
    const Real scale = getScale() > 0. ? getScale() : 1.;
    vtkPolyData* data = transformVTKPolyData(getPlacement(), scale, 
                                             app->GetOutput());
}

} // namespace SimTK

