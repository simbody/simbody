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
#include "simbody/internal/VTKReporter.h"

#include "vtkCommand.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkAssembly.h"
#include "vtkCamera.h"
#include "vtkActor.h"
#include "vtkFollower.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"

#include "vtkLineSource.h"
#include "vtkTextSource.h"
#include "vtkVectorText.h"
#include "vtkConeSource.h"
#include "vtkSphereSource.h"
#include "vtkCubeSource.h"
#include "vtkPointSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkCaptionActor2D.h"
#include "vtkAppendPolyData.h"
#include "vtkInteractorObserver.h"
#include "vtkInteractorStyleTrackballCamera.h"

#include "windows.h" // kludge

#include <cassert>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

namespace SimTK {
static const Real Pi = std::acos(-1.), RadiansPerDegree = Pi/180;
static const int  GroundBodyNum = 0; // ground is always body 0
    
typedef std::pair<vtkProp3D*, Transform> BodyActor;
typedef std::vector<BodyActor>           ActorList;

class VTKDecoration {
public:
    VTKDecoration() : data(0), actor(0) { }
    ~VTKDecoration() {clear();}
    VTKDecoration(const VTKDecoration& d) : data(0), actor(0) {
        copyFrom(d);
    }

    VTKDecoration& operator=(const VTKDecoration& d) {
        if (&d != this)
            copyFrom(d);
        return *this;
    }

    vtkActor* getActor() {return actor;} // not a copy; use ShallowCopy if you need one
    vtkPolyData* getPolyData() {return data;}

    void setScale  (Real s)        {assert(actor);actor->SetScale(s);}
    void setColor  (const Vec3& c) {assert(actor);actor->GetProperty()->SetColor(c[0],c[1],c[2]);}
    void setOpacity(Real o)        {assert(actor);actor->GetProperty()->SetOpacity(o);}

protected:
    void setPolyData(vtkPolyData* output) {
        data = vtkPolyData::New();
        data->DeepCopy(output);
        createActor();
        setDefaultColor();
    }

private:
    vtkPolyData* data;
    vtkActor*    actor;

    void clear() {
        if (actor) actor->Delete(); actor=0; 
        if (data)  data->Delete();  data=0;
    }

    void copyFrom(const VTKDecoration& d) {
        clear();
        if (d.data) {
            setPolyData(d.data);
            actor->SetProperty(d.actor->GetProperty());
            actor->SetScale(d.actor->GetScale());
        }
    }

    void setDefaultColor() {
        setColor(Vec3(0,0,1)); // blue
    }

    void createActor() {
        if (actor) actor->Delete();
        actor = vtkActor::New();
        vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
        mapper->SetInput(data);
        actor->SetMapper(mapper);
        mapper->Delete();
    }
};

class VTKSphere : public VTKDecoration {
public:
    VTKSphere(Real radius) {
        vtkSphereSource *sphere = vtkSphereSource::New();
        sphere->SetRadius(radius);
        sphere->SetThetaResolution(18);
        sphere->SetPhiResolution(18);
        sphere->Update();
        setPolyData(sphere->GetOutput());
        sphere->Delete();
    }
};

void VTKReporter::addDecoration(int bodyNum, VTKDecoration& d, Transform X_GD) {
    addActor(bodyNum, d.getActor(), X_GD);
}

void VTKReporter::addActor(int bodyNum, vtkActor* a, Transform X_GA) {
    ActorList& actors = bodies[bodyNum];
    vtkActor* acopy = vtkActor::New();
    acopy->ShallowCopy(a);
    actors.push_back(BodyActor(acopy,X_GA));
    renderer->AddActor(actors[actors.size()-1].first);
    if (bodyNum==0) setConfiguration(0, Transform()); // just do this once
}

VTKReporter::~VTKReporter() {deletePointers();}

VTKReporter::VTKReporter(const MultibodySystem& m) 
    : mbs(m) {
    zeroPointers();
    makeShapes();

    renWin = vtkRenderWindow::New();
    renWin->SetSize(600,600);
    
    // an interactor
    vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin); 
    vtkInteractorStyleTrackballCamera* style=vtkInteractorStyleTrackballCamera::New(); 
    iren->SetInteractorStyle(style);
    style->Delete();
    iren->Initialize(); // register interactor to pick up windows messages


    renderer = vtkRenderer::New();
    renderer->SetBackground(1,1,1); // white
    renderer->GetActiveCamera()->SetFocalPoint(0,0,0);
    renderer->GetActiveCamera()->SetPosition(0,1,5);


    renWin->AddRenderer(renderer);

    const MechanicalSubsystem& sbs = mbs.getMechanicalSubsystem();
    bodies.resize(sbs.getNBodies());

    for (int i=0; i<(int)bodies.size(); ++i) {
        vtkActor *aAxes = vtkActor::New();
        aAxes->SetMapper(axesMapper);
        aAxes->GetProperty()->SetColor(0,i==0?1:0,0); // green for ground, else black
        aAxes->GetProperty()->SetLineWidth(3);
        addActor(i, aAxes, Transform());
        aAxes->Delete();

        vtkActor *aSphere = vtkActor::New();
        //aSphere->SetMapper(sphereMapper);
        aSphere->SetMapper(cubeMapper);
        aSphere->GetProperty()->SetColor(0,0,1); // sphere color blue
        //aSphere->GetProperty()->SetRepresentationToWireframe();
        aSphere->GetProperty()->SetOpacity(0.3);
        aSphere->SetScale(0.1);

        addActor(i, aSphere, Transform()); // at the origin

        aSphere->Delete();
    }


    renWin->Render();
}

void VTKReporter::report(const State& s) {
    if (!renWin) return;

    const MechanicalSubsystem& mech = mbs.getMechanicalSubsystem();
    for (int i=1; i<mech.getNBodies(); ++i) {
        const Transform& config = mech.getBodyConfiguration(s, i);
        setConfiguration(i, config);
    }

    renWin->Render();

    // Process any window messages since last time
    //TODO: Win32 specific
    MSG msg;
    bool done = false;
    while (PeekMessage(&msg, NULL, 0, 0, PM_REMOVE)) {
        if (msg.message == WM_QUIT) {
            printf("quit!!\n");
            renWin->Delete();renWin=0;
            break;
        }
        TranslateMessage(&msg);
        DispatchMessage(&msg);
    }

}

void VTKReporter::zeroPointers() {
    renWin=0; renderer=0; sphereMapper=0; cubeMapper=0;
    lineMapper=0; cylinderMapper=0; axesMapper=0;
}
void VTKReporter::deletePointers() {
    for (int i=0; i<(int)bodies.size(); ++i) {
        ActorList& actors = bodies[i];
        for (int a=0; a<(int)actors.size(); ++a)
            actors[a].first->Delete(), actors[a].first=0;
    }
    if(axesMapper)axesMapper->Delete(); 
    if(cylinderMapper)cylinderMapper->Delete(); 
    if(lineMapper)lineMapper->Delete();
    if(cubeMapper)cubeMapper->Delete(); 
    if(sphereMapper)sphereMapper->Delete();
    if(renderer)renderer->Delete();
    if(renWin)renWin->Delete();
    zeroPointers();
}
void VTKReporter::makeShapes() {
    // sphere
    vtkSphereSource *sphere = vtkSphereSource::New();
    sphere->SetRadius(1.);
    sphere->SetThetaResolution(18);
    sphere->SetPhiResolution(18);
    sphereMapper = vtkPolyDataMapper::New();
    sphereMapper->SetInput(sphere->GetOutput());
    sphere->Delete();

    // cube
    vtkCubeSource *cube = vtkCubeSource::New();
    cube->SetXLength(1.); cube->SetYLength(1.); cube->SetZLength(1.);
    cubeMapper = vtkPolyDataMapper::New();
    cubeMapper->SetInput(cube->GetOutput());
    cube->Delete();    

    // line
    vtkLineSource *line = vtkLineSource::New();
    line->SetPoint1(0,0,0);
    line->SetPoint2(1,0,0);
    lineMapper = vtkPolyDataMapper::New();
    lineMapper->SetInput(line->GetOutput());
    line->Delete();

    // axes
    vtkAppendPolyData* app = vtkAppendPolyData::New();
    line = vtkLineSource::New(); line->SetPoint1(0,0,0); line->SetPoint2(1,0,0);
    app->AddInput(line->GetOutput());
    line->Delete();
    //cube = vtkCubeSource::New(); 
    //cube->SetXLength(1); cube->SetYLength(0.01); cube->SetZLength(0.01);
    //cube->SetCenter(.5,0,0);
    //app->AddInput(cube->GetOutput());
    //cube->Delete();
    line = vtkLineSource::New(); line->SetPoint1(0,0,0); line->SetPoint2(0,1,0);
    app->AddInput(line->GetOutput());
    line->Delete();
    line = vtkLineSource::New(); line->SetPoint1(0,0,0); line->SetPoint2(0,0,1);
    app->AddInput(line->GetOutput());
    line->Delete();

    vtkVectorText* xtext = vtkVectorText::New();
    xtext->SetText("x"); // default size is around 1
    vtkTransform* t = vtkTransform::New();
    t->Translate(1,-0.025,0); t->Scale(.1,.1,1); 
    vtkTransformPolyDataFilter* trf = vtkTransformPolyDataFilter::New();
    trf->SetInput(xtext->GetOutput());
    trf->SetTransform(t); 
    t->Delete();
    app->AddInput(trf->GetOutput());
    trf->Delete();

    xtext->Delete();


    axesMapper = vtkPolyDataMapper::New();
    axesMapper->SetInput(app->GetOutput());
    app->Delete();
}

void VTKReporter::setConfiguration(int bodyNum, const Transform& X_GB) {
    const ActorList& actors = bodies[bodyNum];
    for (int i=0; i < (int)actors.size(); ++i) {
        vtkProp3D*       actor = actors[i].first;
        const Transform& X_BA  = actors[i].second;
        const Transform  X_GA  = X_GB*X_BA;
        actor->SetPosition(X_GA.T()[0], X_GA.T()[1], X_GA.T()[2]);
        const Vec4 av = X_GA.R().convertToAngleAxis();
        actor->SetOrientation(0,0,0);
        actor->RotateWXYZ(av[0]/RadiansPerDegree, av[1], av[2], av[3]);
    }

}

} // namespace SimTK

