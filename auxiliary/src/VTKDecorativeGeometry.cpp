/* Portions copyright (c) 2006 Stanford University and Jack Middleton.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/AnalyticGeometry.h"
#include "simbody/internal/DecorativeGeometry.h"

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
#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkCylinderSource.h"
#include "vtkVectorText.h"

#include "vtkObject.h"
#include "VTKDecorativeGeometry.h"

#include <cmath>

namespace SimTK {
    static const Real Pi = std::acos(-1.), RadiansPerDegree = Pi/180;


    ///////////////////////////
    // VTKDecorativeGeometry //
    ///////////////////////////


/*  TODO JACKM  DO we need clone?
VTKDecorativeGeometry::VTKDecorativeGeometry(const VTKDecorativeGeometry& src)  {
             src.clone();
}

VTKDecorativeGeometry& VTKDecorativeGeometry::operator=(const VTKDecorativeGeometry& src) {
    if (&src != this) {
             src.clone();
    }
    return *this;
}
*/




vtkTransform* 
VTKDecorativeGeometry::createVTKTransform(const Transform& X_BG, const Real& s) {

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
VTKDecorativeGeometry::transformVTKPolyData(const Transform& X_BG, const Real& s, 
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

void VTKDecorativeLine::createVTKPolyData(const Vec3& p1, const Vec3& p2) {
    vtkLineSource* line = vtkLineSource::New(); 
    rememberVTKObject(line);


    line->SetPoint1(p1[0],p1[1],p1[2]); line->SetPoint2(p2[0],p2[1],p2[2]);

    const Real scale = myHandle->getScale() > 0. ? myHandle->getScale() : 1.;
    vtkPolyData* data = transformVTKPolyData(myHandle->getPlacement(), scale, 
                                             line->GetOutput());
}
void VTKDecorativeCircle::createVTKPolyData(Real r) {
    assert(!"DecorativeCircle::createVTKPolyData() NOT IMPLEMENTED YET");
}

void VTKDecorativeSphere::createVTKPolyData(Real r) {
    vtkSphereSource *sphere = vtkSphereSource::New();
    rememberVTKObject(sphere);

    sphere->SetRadius(r);

    int res = DefaultResolution;
    if (myHandle->getResolution() > 0.) 
        res = (int)(res*myHandle->getResolution()+0.5);
    sphere->SetThetaResolution(res);
    sphere->SetPhiResolution(res);

    const Real scale = myHandle->getScale() > 0. ? myHandle->getScale() : 1.;
    vtkPolyData* data = transformVTKPolyData(myHandle->getPlacement(), scale, 
                                             sphere->GetOutput());
}


void VTKDecorativeBrick::createVTKPolyData(const Vec3& h) {
    vtkCubeSource* cube = vtkCubeSource::New();
    rememberVTKObject(cube);


    cube->SetXLength(2*h[0]);
    cube->SetYLength(2*h[1]);
    cube->SetZLength(2*h[2]);

    // resolution is ignored -- our needs are few for rectangles!

    const Real scale = myHandle->getScale() > 0. ? myHandle->getScale() : 1.;

    vtkPolyData* data = transformVTKPolyData(myHandle->getPlacement(), scale, 
                                             cube->GetOutput());
}

void VTKDecorativeCylinder::createVTKPolyData(Real r, Real halfHeight) {
    vtkCylinderSource* cyl = vtkCylinderSource::New();
    rememberVTKObject(cyl);

    cyl->SetRadius(r);
    cyl->SetHeight(2*halfHeight);

    int res = DefaultResolution;
    if (myHandle->getResolution() > 0.) 
        res = (int)(res*myHandle->getResolution()+0.5);
    cyl->SetResolution(res);

    const Real scale = myHandle->getScale() > 0. ? myHandle->getScale() : 1.;
    vtkPolyData* data = transformVTKPolyData(myHandle->getPlacement(), scale, 
                                             cyl->GetOutput());
}

void VTKDecorativeFrame::createVTKPolyData(Real length) {
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
    vtkPolyData* label = transformVTKPolyData(Transform(Vec3(length,-0.05*length,0)), length*0.2, 
                                              xtext->GetOutput());
    app->AddInput(label);
    const Real scale = myHandle->getScale() > 0. ? myHandle->getScale() : 1.;
    vtkPolyData* data = transformVTKPolyData(myHandle->getPlacement(), scale, 
                                             app->GetOutput());
}

} // namespace SimTK

