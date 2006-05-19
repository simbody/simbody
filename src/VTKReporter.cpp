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
static const Vec3 DefaultGroundBodyColor(0,1,0); // green
static const Vec3 DefaultBaseBodyColor(1,0,0);   // red
static const Vec3 DefaultBodyColor(0,0,0);       // black

class VTKReporterRep {
public:
    // no default constructor -- must have MultibodySystem always
    VTKReporterRep(const MultibodySystem& m);

    ~VTKReporterRep() {
        deletePointers();
    }

    // This will make a copy of the supplied DecorativeGeometry.
    void addDecoration(int bodyNum, const Transform& X_GD, const DecorativeGeometry&);

    void setDefaultBodyColor(int bodyNum, const Vec3& rgb) {
        bodies[bodyNum].defaultColorRGB = rgb;
    }
    const Vec3& getDefaultBodyColor(int body) const {return bodies[body].defaultColorRGB;}

    VTKReporterRep* clone() const {
        VTKReporterRep* dup = new VTKReporterRep(*this);
        dup->myHandle = 0;
        return dup;
    }

    void report(const State& s);

    void setMyHandle(VTKReporter& h) {myHandle = &h;}
    void clearMyHandle() {myHandle=0;}

private:
    friend class VTKReporter;
    VTKReporter* myHandle;     // the owner of this rep

    const MultibodySystem& mbs;

    struct PerBodyInfo {
        std::vector<vtkProp3D*>         aList;
        std::vector<DecorativeGeometry> gList; // one per actor (TODO)
        Vec3        defaultColorRGB;
    };
    std::vector<PerBodyInfo> bodies;

    vtkRenderWindow* renWin;
    vtkRenderer*     renderer;

    void zeroPointers();
    void deletePointers();
    void setConfiguration(int bodyNum, const Transform& X_GB);
};

    /////////////////
    // VTKReporter //
    /////////////////


bool VTKReporter::isOwnerHandle() const {
    return rep==0 || rep->myHandle==this;
}
bool VTKReporter::isEmptyHandle() const {return rep==0;}

VTKReporter::VTKReporter(const MultibodySystem& m) : rep(0) {
    rep = new VTKReporterRep(m);
    rep->setMyHandle(*this);
}

VTKReporter::VTKReporter(const VTKReporter& src) : rep(0) {
    if (src.rep) {
        rep = src.rep->clone();
        rep->setMyHandle(*this);
    }
}

VTKReporter& VTKReporter::operator=(const VTKReporter& src) {
    if (&src != this) {
        delete rep; rep=0;
        if (src.rep) {
            rep = src.rep->clone();
            rep->setMyHandle(*this);
        }
    }
    return *this;
}

VTKReporter::~VTKReporter() {
    delete rep; rep=0;
}

void VTKReporter::report(const State& s) {
    assert(rep);
    rep->report(s);
}

void VTKReporter::addDecoration(int body, const Transform& X_GD,
                                const DecorativeGeometry& g) 
{
    assert(rep);
    rep->addDecoration(body, X_GD, g);
}


    ////////////////////
    // VTKReporterRep //
    ////////////////////

void VTKReporterRep::addDecoration(int body, const Transform& X_GD,
                                   const DecorativeGeometry& g)
{
    // For now we create a unique actor for each piece of geometry
    vtkActor* actor = vtkActor::New();
    bodies[body].aList.push_back(actor);
    DecorativeGeometry tmp = g;
    tmp = g;
    bodies[body].gList.push_back(g);
    DecorativeGeometry& geom  = bodies[body].gList.back();

    // Apply the transformation.
    geom.setPlacement(X_GD*geom.getPlacement());
    vtkPolyData* poly = geom.updVTKPolyData();

    // Now apply the actor-level properties from the geometry.
    const Vec3 color = (geom.getColor()[0] != -1 ? geom.getColor() : getDefaultBodyColor(body)); 
    actor->GetProperty()->SetColor(color[0],color[1],color[2]);

    const Real opacity = (geom.getOpacity() != -1 ? geom.getOpacity() : Real(1));
    actor->GetProperty()->SetOpacity(opacity);

    const Real lineWidth = (geom.getLineThickness() != -1 ? geom.getLineThickness() : Real(1));
    actor->GetProperty()->SetLineWidth(lineWidth);

    const int representation = (geom.getRepresentation() != -1 ? geom.getRepresentation() : VTK_SURFACE);
    actor->GetProperty()->SetRepresentation(representation);

    // Set up the mapper & register actor with renderer
    vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
    mapper->SetInput(poly);
    actor->SetMapper(mapper);
    mapper->Delete(); mapper=0; // remove now-unneeded mapper reference
    renderer->AddActor(actor);

}

VTKReporterRep::VTKReporterRep(const MultibodySystem& m) 
    : myHandle(0), mbs(m) 
{
    zeroPointers();

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
    //renderer->GetActiveCamera()->SetFocalPoint(0,0,0);
    renderer->GetActiveCamera()->SetPosition(0,1,10);


    renWin->AddRenderer(renderer);

    const MechanicalSubsystem& sbs = mbs.getMechanicalSubsystem();
    bodies.resize(sbs.getNBodies());

    setDefaultBodyColor(GroundBodyNum, DefaultGroundBodyColor);
    for (int i=1; i<(int)bodies.size(); ++i) {
        if (sbs.getParent(i) == GroundBodyNum)
             setDefaultBodyColor(i, DefaultBaseBodyColor);
        else setDefaultBodyColor(i, DefaultBodyColor);
    }

    for (int i=0; i<(int)bodies.size(); ++i) {
        DecorativeFrame axes(1);
        axes.setLineThickness(3);
        addDecoration(i, Transform(), axes); // the body frame

        // Display the inboard joint frame (at half size), unless it is the
        // same as the body frame. Then find the corresponding frame on the
        // parent and display that in this body's color.
        if (i > 0) {
            const Transform& jInb = sbs.getJointFrame(State(), i);
            if (jInb.T() != Vec3(0) || jInb.R() != Mat33(1)) {
                addDecoration(i, jInb, DecorativeFrame(0.5));
                if (jInb.T() != Vec3(0))
                    addDecoration(i, Transform(), DecorativeLine(Vec3(0), jInb.T()));
            }
            const Transform& jParent = sbs.getJointFrameOnParent(State(), i);
            DecorativeFrame frameOnParent(0.5);
            frameOnParent.setColor(getDefaultBodyColor(i));
            addDecoration(sbs.getParent(i), jParent, frameOnParent);
            if (jParent.T() != Vec3(0))
                addDecoration(sbs.getParent(i), Transform(), DecorativeLine(Vec3(0),jParent.T()));
        }

        // Put a little black wireframe sphere at the COM, and add a line from 
        // body origin to the com.

        DecorativeSphere com(.1);
        com.setResolution(1./3); // chunky is fine
        com.setColor(Vec3(0,0,0));
        com.setRepresentationToWireframe();
        const Vec3& comPos = sbs.getBodyCenterOfMass(State(), i);
        addDecoration(i, Transform(comPos), com);
        if (comPos != Vec3(0))
            addDecoration(i, Transform(), DecorativeLine(Vec3(0), comPos));
    }

    renWin->Render();
}

void VTKReporterRep::report(const State& s) {
    if (!renWin) return;

    mbs.realize(s, Stage::Configured); // just in case

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

void VTKReporterRep::zeroPointers() {
    renWin=0; renderer=0;
}
void VTKReporterRep::deletePointers() {
    // Delete all the actors. The geometry gets deleted automatically
    // thanks to good design!
    for (int i=0; i<(int)bodies.size(); ++i) {
        std::vector<vtkProp3D*>& actors = bodies[i].aList;
        for (int a=0; a<(int)actors.size(); ++a)
            actors[a]->Delete(), actors[a]=0;
    }

    if(renderer)renderer->Delete();
    if(renWin)renWin->Delete();
    zeroPointers();
}

void VTKReporterRep::setConfiguration(int bodyNum, const Transform& X_GB) {
    const std::vector<vtkProp3D*>& actors = bodies[bodyNum].aList;
    for (int i=0; i < (int)actors.size(); ++i) {
        vtkProp3D*       actor = actors[i];
        actor->SetPosition(X_GB.T()[0], X_GB.T()[1], X_GB.T()[2]);
        const Vec4 av = X_GB.R().convertToAngleAxis();
        actor->SetOrientation(0,0,0);
        actor->RotateWXYZ(av[0]/RadiansPerDegree, av[1], av[2], av[3]);
    }
}

} // namespace SimTK

