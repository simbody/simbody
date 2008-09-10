/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Jack Middleton                                                    *
 * Contributors: Michael Sherman                                              *
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

/** @file
 * Implementation of the VTKDecorativeGeometry class.
 */

#include "SimTKcommon.h"

#include "simbody/internal/common.h"

#include "VTKDecorativeGeometry.h"

#include "vtkCommand.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkActor.h"
#include "vtkFollower.h"

#include "vtkPolyDataMapper.h"
#include "vtkCaptionActor2D.h"
#include "vtkInteractorObserver.h"
#include "vtkInteractorStyleTrackballCamera.h"


#include "vtkProperty.h"
#include "vtkAssembly.h"
#include "vtkAppendPolyData.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkPolyData.h"
#include "vtkPolyDataAlgorithm.h"
#include "vtkLineSource.h"
#include "vtkDiskSource.h"
#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkCylinderSource.h"
#include "vtkVectorText.h"
#include "vtkCellArray.h"

#include "vtkObject.h"

#include <cmath>

using namespace SimTK;

static const Real RadiansPerDegree = Pi/180;


    ///////////////////////////
    // VTKDecorativeGeometry //
    ///////////////////////////

vtkPolyData* 
VTKDecorativeGeometry::getVTKPolyData() {
    assert(vtkObjects.size());
    return vtkPolyDataAlgorithm::SafeDownCast(vtkObjects.back())->GetOutput();
}

vtkTransform* 
VTKDecorativeGeometry::createVTKTransform(const Transform& X_BG, const Vec3& s) {

    const Vec3 t = X_BG.T();
    const Vec4 r = X_BG.R().convertRotationToAngleAxis();

    vtkTransform* xform = vtkTransform::New();  // starts out as identity
    rememberVTKObject(xform);

    xform->Translate(t[0],t[1],t[2]);
    xform->RotateWXYZ(r[0]/RadiansPerDegree,r[1],r[2],r[3]);     // angle, axis
    xform->Scale(s[0],s[1],s[2]); // scale and stretch
    return xform;
}

vtkPolyData*
VTKDecorativeGeometry::transformVTKPolyData(const Transform& X_BG, const Vec3& s, 
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


void VTKDecorativeGeometry::deleteVTKGeometry() { // Delete in reverse order of allocation
    for (int i=(int)vtkObjects.size()-1; i >= 0; --i) {
        vtkObject* obj = vtkObjects[i];
        obj->Delete();
        vtkObjects[i]=0;
    }
    vtkObjects.resize(0);
}

    // Implementations of pure virtuals for building specific geometric objects. //

void VTKDecorativeGeometry::implementLineGeometry(const DecorativeLine& dline) {
    const Vec3& p1 = dline.getPoint1();
    const Vec3& p2 = dline.getPoint2();
    vtkLineSource* vline = vtkLineSource::New(); 
    rememberVTKObject(vline);

    vline->SetPoint1(p1[0],p1[1],p1[2]); vline->SetPoint2(p2[0],p2[1],p2[2]);

    const Real scale = dline.getScale() > 0. ? dline.getScale() : 1.;
    vtkPolyData* data = transformVTKPolyData(dline.getTransform(), Vec3(scale), 
                                             vline->GetOutput());
    // Not using "data" here -- transform also appended to the vtkObjects list.
}

void VTKDecorativeGeometry::implementCircleGeometry(const DecorativeCircle& dcircle) {
    static const int DefaultResolution = 15;

    const Real r = dcircle.getRadius();
    vtkDiskSource *vcircle = vtkDiskSource::New();
    rememberVTKObject(vcircle);

    vcircle->SetOuterRadius(r);
    vcircle->SetInnerRadius(0.0);

    int res = DefaultResolution;
    if (dcircle.getResolution() > 0.) 
        res = (int)(res*dcircle.getResolution()+0.5);
    vcircle->SetCircumferentialResolution(res);

    const Real scale = dcircle.getScale() > 0. ? dcircle.getScale() : 1.;
    vtkPolyData* data = transformVTKPolyData(dcircle.getTransform(), Vec3(scale), 
                                             vcircle->GetOutput());
    // Not using "data" here -- transform also appended to the vtkObjects list.
}

void VTKDecorativeGeometry::implementSphereGeometry(const DecorativeSphere& dsphere) {
    static const int DefaultResolution = 15;

    const Real r = dsphere.getRadius();
    vtkSphereSource *vsphere = vtkSphereSource::New();
    rememberVTKObject(vsphere);

    vsphere->SetRadius(r);

    int res = DefaultResolution;
    if (dsphere.getResolution() > 0.) 
        res = (int)(res*dsphere.getResolution()+0.5);
    vsphere->SetThetaResolution(res);
    vsphere->SetPhiResolution(res);

    const Real scale = dsphere.getScale() > 0. ? dsphere.getScale() : 1.;
    vtkPolyData* data = transformVTKPolyData(dsphere.getTransform(), Vec3(scale), 
                                             vsphere->GetOutput());
    // Not using "data" here -- transform also appended to the vtkObjects list.
}



void VTKDecorativeGeometry::implementEllipsoidGeometry(const DecorativeEllipsoid& dellipsoid) {
    static const int DefaultResolution = 15;

    const Vec3& r = dellipsoid.getRadii();
    vtkSphereSource* vellipsoid = vtkSphereSource::New();
    rememberVTKObject(vellipsoid);

    // Make a sphere of unit radius, then we have to distort in x, y, and z.
    vellipsoid->SetRadius(1);

    int res = DefaultResolution;
    if (dellipsoid.getResolution() > 0.) 
        res = (int)(res*dellipsoid.getResolution()+0.5);
    vellipsoid->SetThetaResolution(res);
    vellipsoid->SetPhiResolution(res);

    const Real scale = dellipsoid.getScale() > 0. ? dellipsoid.getScale() : 1.;

    vtkPolyData* data = transformVTKPolyData(dellipsoid.getTransform(), scale*r, 
                                             vellipsoid->GetOutput());
    // Not using "data" here -- transform also appended to the vtkObjects list.
}

void VTKDecorativeGeometry::implementBrickGeometry(const DecorativeBrick& dbrick) {
    const Vec3& h = dbrick.getHalfLengths();
    vtkCubeSource* vcube = vtkCubeSource::New();
    rememberVTKObject(vcube);

    vcube->SetXLength(2*h[0]);
    vcube->SetYLength(2*h[1]);
    vcube->SetZLength(2*h[2]);

    // resolution is ignored -- our needs are few for rectangles!

    const Real scale = dbrick.getScale() > 0. ? dbrick.getScale() : 1.;

    vtkPolyData* data = transformVTKPolyData(dbrick.getTransform(), Vec3(scale), 
                                             vcube->GetOutput());
    // Not using "data" here -- transform also appended to the vtkObjects list.
}

void VTKDecorativeGeometry::implementCylinderGeometry(const DecorativeCylinder& dcyl) {
    static const int DefaultResolution = 10;

    const Real r          = dcyl.getRadius();
    const Real halfHeight = dcyl.getHalfHeight();

    vtkCylinderSource* vcyl = vtkCylinderSource::New();
    rememberVTKObject(vcyl);

    vcyl->SetRadius(r);
    vcyl->SetHeight(2*halfHeight);

    int res = DefaultResolution;
    if (dcyl.getResolution() > 0.) 
        res = (int)(res*dcyl.getResolution()+0.5);
    vcyl->SetResolution(res);

    const Real scale = dcyl.getScale() > 0. ? dcyl.getScale() : 1.;
    vtkPolyData* data = transformVTKPolyData(dcyl.getTransform(), Vec3(scale), 
                                             vcyl->GetOutput());
    // Not using "data" here -- transform also appended to the vtkObjects list.
}

void VTKDecorativeGeometry::implementFrameGeometry(const DecorativeFrame& dframe) {
    const Real length = dframe.getAxisLength();

    vtkAppendPolyData* app = vtkAppendPolyData::New();
    rememberVTKObject(app);

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
    vtkPolyData* label = transformVTKPolyData(Transform(Vec3(length,-0.05*length,0)), Vec3(length*0.2), 
                                              xtext->GetOutput());
    app->AddInput(label);
    const Real scale = dframe.getScale() > 0. ? dframe.getScale() : 1.;
    vtkPolyData* data = transformVTKPolyData(dframe.getTransform(), Vec3(scale), 
                                             app->GetOutput());
    // Not using "data" here -- transform also appended to the vtkObjects list.
}

void VTKDecorativeGeometry::implementTextGeometry(const DecorativeText& dtext) {
    vtkVectorText *vtext = vtkVectorText::New();
    rememberVTKObject(vtext);

    vtext->SetText(dtext.getText().c_str());

    const Real scale = dtext.getScale() > 0. ? dtext.getScale() : 1.;
    vtkPolyData* data = transformVTKPolyData(dtext.getTransform(), Vec3(scale), 
                                             vtext->GetOutput());
    // Not using "data" here -- transform also appended to the vtkObjects list.
}

void VTKDecorativeGeometry::implementMeshGeometry(const DecorativeMesh& dmesh) {
    vtkPolyData* vmesh = vtkPolyData::New();
    rememberVTKObject(vmesh);

    vtkPoints *points = vtkPoints::New();
    vtkCellArray *faces = vtkCellArray::New();
    const PolygonalMesh& mesh = dmesh.getMesh();
    for (int i = 0; i < mesh.getNumVertices(); i++) {
        const Vec3& pos = mesh.getVertexPosition(i);
        points->InsertNextPoint(pos[0], pos[1], pos[2]);
    }
    for (int i = 0; i < mesh.getNumFaces(); i++) {
        int numVerts = mesh.getNumVerticesForFace(i);
        faces->InsertNextCell(numVerts);
        for (int j = 0; j < numVerts; j++)
            faces->InsertCellPoint(mesh.getFaceVertex(i, j));
    }
    vmesh->SetPoints(points);
    points->Delete();
    vmesh->SetPolys(faces);
    faces->Delete();
    
    const Real scale = dmesh.getScale() > 0. ? dmesh.getScale() : 1.;
    vtkPolyData* data = transformVTKPolyData(dmesh.getTransform(), Vec3(scale), vmesh);
    // Not using "data" here -- transform also appended to the vtkObjects list.
}

