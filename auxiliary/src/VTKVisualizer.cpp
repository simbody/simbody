/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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

#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/VTKVisualizer.h"

#include "VTKDecorativeGeometry.h"

#include "vtkCommand.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkProperty.h"
#include "vtkAssembly.h"
#include "vtkCamera.h"
#include "vtkLight.h"
#include "vtkActor.h"
#include "vtkFollower.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkToolkits.h"

#include "vtkPolyDataMapper.h"
#include "vtkCaptionActor2D.h"
#include "vtkAppendPolyData.h"
#include "vtkInteractorObserver.h"
#include "vtkInteractorStyleTrackballCamera.h"

#ifdef _WIN32
#include "windows.h" // kludge
#endif

#ifdef VTK_USE_CARBON
#include <Carbon/Carbon.h>
#endif

#ifdef VTK_USE_X
#include <X11/Intrinsic.h>
#include "vtkXRenderWindowInteractor.h"
#endif

#include <cassert>
#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

using namespace SimTK;

static const Real RadiansToDegrees = (Real)SimTK_RADIAN_TO_DEGREE;

static const Vec3 DefaultGroundBodyColor = Green;
static const Vec3 DefaultBaseBodyColor   = Red;
static const Vec3 DefaultBodyColor       = Gray;

class SimTK::VTKVisualizerRep {
public:
    // no default constructor -- must have MultibodySystem always
    VTKVisualizerRep(const MultibodySystem& m, VTKVisualizer *reporter);

    ~VTKVisualizerRep() {
        deletePointers();
    }

    void disableDefaultGeometry() { defaultBodyScaleForAutoGeometry=0.;}

    // This will make a copy of the supplied DecorativeGeometry.
    // These are topology-stage decorations which we can precalculate (at least in part)
    // since they will be present in every rendered frame.
    void addDecoration(MobilizedBodyIndex bodyNum, const Transform& X_GD, const DecorativeGeometry&);
    void addRubberBandLine(MobilizedBodyIndex b1, const Vec3& station1, MobilizedBodyIndex b2, const Vec3& station2,
                           const DecorativeLine&);

    // This geometry survives only until the next frame is rendered, then evaporates.
    void addEphemeralDecoration(const DecorativeGeometry&);

    // Make sure everything can be seen.
    void resetCamera() {cameraNeedsToBeReset=true;}
    
    VTKVisualizerRep* clone() const {
        VTKVisualizerRep* dup = new VTKVisualizerRep(*this);
        dup->myHandle = 0;
        return dup;
    }

    void report(const State& s);

    void setMyHandle(VTKVisualizer& h) {myHandle = &h;}
    void clearMyHandle() {myHandle=0;}

private:
    friend class VTKVisualizer;
    VTKVisualizer* myHandle;     // the owner of this rep

    bool cameraNeedsToBeReset; // report() checks and clears this

    Real defaultBodyScaleForAutoGeometry;

    const MultibodySystem& mbs;

    struct PerBodyInfo {
        PerBodyInfo() { }
        std::vector<vtkProp3D*>         aList;
        std::vector<DecorativeGeometry> gList; // one per actor (TODO)
    };
    std::vector<PerBodyInfo> bodies;

    struct PerDynamicGeomInfo {
        PerDynamicGeomInfo() : actor(0) { }
        vtkActor*      actor;
        DecorativeLine line;
        MobilizedBodyIndex  body1, body2;
        Vec3 station1, station2;
    };
    std::vector<PerDynamicGeomInfo> dynamicGeom;

    // This geometry gets displayed at the next frame render and then 
    // destroyed. We have to remember the actors we generate to do that so 
    // we can remove them from the renderer when we're done with the frame.
    std::vector<DecorativeGeometry> ephemeralGeometry;
    std::vector<vtkActor*>    ephemeralActors;

    vtkRenderWindow* renWin;
    vtkRenderer*     renderer;

    void initTopology();
    void createInstanceGeometry(const State& state);
    void zeroPointers();
    void deletePointers();
    void setConfiguration(MobilizedBodyIndex bodyNum, const Transform& X_GB);
    void setRubberBandLine(int dgeom, const Vec3& p1, const Vec3& p2);

    void displayEphemeralGeometry(const State& s);

    int convertToVTKRepresentation(DecorativeGeometry::Representation drawMode) {
        int vtkDrawMode = -1;
        switch(drawMode) {
            case DecorativeGeometry::DrawPoints:    vtkDrawMode=VTK_POINTS;    break;
            case DecorativeGeometry::DrawWireframe: vtkDrawMode=VTK_WIREFRAME; break; 
            case DecorativeGeometry::DrawSurface:   vtkDrawMode=VTK_SURFACE;   break; 
            default: assert(!"unrecognized drawing mode");
        }
        return vtkDrawMode;
    }

};

    ///////////////////
    // VTKVisualizer //
    ///////////////////


bool VTKVisualizer::isOwnerHandle() const {
    return rep==0 || rep->myHandle==this;
}
bool VTKVisualizer::isEmptyHandle() const {return rep==0;}

VTKVisualizer::VTKVisualizer(const MultibodySystem& m) : rep(0) {
    rep = new VTKVisualizerRep(m, this);
}

VTKVisualizer::VTKVisualizer(const VTKVisualizer& src) : rep(0) {
    if (src.rep) {
        rep = src.rep->clone();
        rep->setMyHandle(*this);
    }
}

VTKVisualizer& VTKVisualizer::operator=(const VTKVisualizer& src) {
    if (&src != this) {
        delete rep; rep=0;
        if (src.rep) {
            rep = src.rep->clone();
            rep->setMyHandle(*this);
        }
    }
    return *this;
}

VTKVisualizer::~VTKVisualizer() {
    delete rep; rep=0;
}

void VTKVisualizer::report(const State& s) const {
    assert(rep);
    rep->report(s);
}


void VTKVisualizer::setCameraLocation(const Vec3& p) {
    if (!rep || !rep->renderer) return;
    vtkCamera* camera = rep->renderer->GetActiveCamera();
    camera->SetPosition(p[0], p[1], p[2]);
    camera->ComputeViewPlaneNormal();
}

void VTKVisualizer::setCameraFocalPoint(const Vec3& p) {
    if (!rep || !rep->renderer) return;
    vtkCamera* camera = rep->renderer->GetActiveCamera();
    camera->SetFocalPoint(p[0], p[1], p[2]);
    camera->ComputeViewPlaneNormal();
}

void VTKVisualizer::setCameraUpDirection(const Vec3& d) {
    if (!rep || !rep->renderer) return;
    vtkCamera* camera = rep->renderer->GetActiveCamera();
    camera->SetViewUp(d[0], d[1], d[2]);
    camera->OrthogonalizeViewUp();
}

void VTKVisualizer::setCameraClippingRange(Real nearPlane, Real farPlane) {
    if (!rep || !rep->renderer) return;
    vtkCamera* camera = rep->renderer->GetActiveCamera();
    camera->SetClippingRange(nearPlane, farPlane);
}

void VTKVisualizer::zoomCameraToIncludeAllGeometry() {
    if (!rep || !rep->renderer) return;
    rep->renderer->ResetCamera();
}

void VTKVisualizer::zoomCamera(Real z) {
    if (!rep || !rep->renderer) return;
    vtkCamera* camera = rep->renderer->GetActiveCamera();
    camera->Zoom(z);
}

void VTKVisualizer::addDecoration(MobilizedBodyIndex body, const Transform& X_GD,
                                const DecorativeGeometry& g) 
{
    assert(rep);
    rep->addDecoration(body, X_GD, g);
}

void VTKVisualizer::addRubberBandLine(MobilizedBodyIndex b1, const Vec3& station1,
                                    MobilizedBodyIndex b2, const Vec3& station2,
                                    const DecorativeLine& g)
{
    assert(rep);
    rep->addRubberBandLine(b1,station1,b2,station2,g);
}

void VTKVisualizer::addEphemeralDecoration(const DecorativeGeometry& g) 
{
    assert(rep);
    rep->addEphemeralDecoration(g);
}

    //////////////////////
    // VTKVisualizerRep //
    //////////////////////

void VTKVisualizerRep::addDecoration(MobilizedBodyIndex body, const Transform& X_GD,
                                   const DecorativeGeometry& g)
{
    class DecorativeGeometryRep;

    // TODO: this should not be done here (sherm 071106)
    initTopology();

    // For now we create a unique actor for each piece of geometry
    vtkActor* actor;
    if (g.getFaceCamera()) {
        vtkFollower* follower = vtkFollower::New();
        follower->SetCamera(renderer->GetActiveCamera());
        actor = follower;
    }
    else
        actor = vtkActor::New();
    bodies[body].aList.push_back(actor);
    bodies[body].gList.push_back(g);
    DecorativeGeometry& dgeom  = bodies[body].gList.back();

    // Apply the transformation.
    dgeom.setTransform(X_GD*dgeom.getTransform());

    // Create the VTK geometry for this DecorativeGometry object. Calls VTKDecorativeGeometry's 
    // implementLine/Sphere/Brick/...Geometry routine as appropriate to generate the 
    // polygons/lines/points for the DecorativeGeometry object. 

    VTKDecorativeGeometry vgeom;
    dgeom.implementGeometry(vgeom);
    vtkPolyData* poly = vgeom.getVTKPolyData(); // retrieve the results

    // Now apply the actor-level properties from the geometry.
    const Vec3 color = (dgeom.getColor()[0] != -1 ? dgeom.getColor() : DefaultBodyColor); 
    actor->GetProperty()->SetColor(color[0],color[1],color[2]);

    const Real opacity = (dgeom.getOpacity() != -1 ? dgeom.getOpacity() : Real(1));
    actor->GetProperty()->SetOpacity(opacity);

    const Real lineWidth = (dgeom.getLineThickness() != -1 ? dgeom.getLineThickness() : Real(1));
    actor->GetProperty()->SetLineWidth(lineWidth);

    const DecorativeGeometry::Representation representation = 
       (dgeom.getRepresentation() != -1 ? dgeom.getRepresentation() 
                                        : DecorativeGeometry::DrawSurface);
    actor->GetProperty()->SetRepresentation( convertToVTKRepresentation(representation) );

    // Set up the mapper & register actor with renderer
    vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
    mapper->SetInput(poly);
    actor->SetMapper(mapper);
    mapper->Delete(); mapper=0; // remove now-unneeded mapper reference
    renderer->AddActor(actor);

    cameraNeedsToBeReset = true;
}

void VTKVisualizerRep::addRubberBandLine(MobilizedBodyIndex b1, const Vec3& station1,
                                       MobilizedBodyIndex b2, const Vec3& station2,
                                       const DecorativeLine& g)
{
        // TODO: this should not be done here (sherm 071106)
    initTopology();

    // Create a unique actor for each piece of geometry.
    int nxt = (int)dynamicGeom.size();
    dynamicGeom.resize(nxt+1);
    PerDynamicGeomInfo& info = dynamicGeom.back();

    info.actor = vtkActor::New();
    info.line  = g;
    info.body1 = b1; info.body2 = b2;
    info.station1 = station1; info.station2 = station2;

    // Now apply the actor-level properties from the geometry.
    const Vec3 color = (info.line.getColor()[0] != -1 ? info.line.getColor() : DefaultBodyColor); 
    info.actor->GetProperty()->SetColor(color[0],color[1],color[2]);

    const Real opacity = (info.line.getOpacity() != -1 ? info.line.getOpacity() : Real(1));
    info.actor->GetProperty()->SetOpacity(opacity);

    const Real lineWidth = (info.line.getLineThickness() != -1 ? info.line.getLineThickness() : Real(1));
    info.actor->GetProperty()->SetLineWidth(lineWidth);

    const DecorativeGeometry::Representation representation =
       (info.line.getRepresentation() != -1 ? info.line.getRepresentation() 
                                            : DecorativeGeometry::DrawSurface);
    info.actor->GetProperty()->SetRepresentation(convertToVTKRepresentation(representation));

    // Set up the mapper & register actor with renderer, but don't set up mapper's input yet.
    vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();
    info.actor->SetMapper(mapper);
    mapper->Delete(); mapper=0; // remove now-unneeded mapper reference
    renderer->AddActor(info.actor);

    cameraNeedsToBeReset = true;
}

void VTKVisualizerRep::addEphemeralDecoration(const DecorativeGeometry& g)
{
    // TODO: this should not be done here (sherm 071106)
    initTopology();

    ephemeralGeometry.push_back(g);
}

void VTKVisualizerRep::displayEphemeralGeometry(const State& s)
{
    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

    // Out with the old ...
    for (int i=0; i < (int)ephemeralActors.size(); ++i) {
        renderer->RemoveActor(ephemeralActors[i]);
        ephemeralActors[i]->Delete();
    }

    // And in with the new ...
    // Create a unique actor for each piece of geometry.
    // TODO: could probably do this with a single actor.
    ephemeralActors.resize(ephemeralGeometry.size());
    for (int i=0; i<(int)ephemeralGeometry.size(); ++i) {
        DecorativeGeometry& dgeom = ephemeralGeometry[i];
        ephemeralActors[i] = vtkActor::New();

        const MobilizedBodyIndex body = MobilizedBodyIndex(dgeom.getBodyId());
        const Transform& X_GB = matter.getMobilizedBody(body).getBodyTransform(s);

        // Apply the transformation.
        dgeom.setTransform(X_GB*dgeom.getTransform());

        // Now apply the actor-level properties from the geometry.
        const Vec3 color = (dgeom.getColor()[0] != -1 ? dgeom.getColor() : DefaultBodyColor); 
        ephemeralActors[i]->GetProperty()->SetColor(color[0],color[1],color[2]);

        const Real opacity = (dgeom.getOpacity() != -1 ? dgeom.getOpacity() : Real(1));
        ephemeralActors[i]->GetProperty()->SetOpacity(opacity);

        const Real lineWidth = (dgeom.getLineThickness() != -1 ? dgeom.getLineThickness() : Real(1));
        ephemeralActors[i]->GetProperty()->SetLineWidth(lineWidth);

        const DecorativeGeometry::Representation representation = 
           (dgeom.getRepresentation() != -1 ? dgeom.getRepresentation() 
                                            : DecorativeGeometry::DrawSurface);
        ephemeralActors[i]->GetProperty()->SetRepresentation(convertToVTKRepresentation(representation));

        // Set up the mapper and render the geometry into it.
        vtkPolyDataMapper* mapper = vtkPolyDataMapper::New();

        VTKDecorativeGeometry vgeom;
        dgeom.implementGeometry(vgeom);
        vtkPolyData* poly = vgeom.getVTKPolyData(); // retrieve the results

        mapper->SetInput(poly);

        // Associate the mapper with the actor, and add the actor to the renderer.
        ephemeralActors[i]->SetMapper(mapper);
        mapper->Delete(); mapper=0; // remove now-unneeded mapper reference
        renderer->AddActor(ephemeralActors[i]);
    }

    //if (ephemeralGeometry.size())
    //    cameraNeedsToBeReset = true; // TODO: is this really a good idea?

    ephemeralGeometry.clear();
}

VTKVisualizerRep::VTKVisualizerRep(const MultibodySystem& m, VTKVisualizer* reporter ) 
    :  mbs(m), cameraNeedsToBeReset(true)
{
    myHandle = reporter;
    const Real cameraScale = 1.0;
    zeroPointers();

    renWin = vtkRenderWindow::New();
    renWin->SetSize(1200,900);
    
    // an interactor
    vtkRenderWindowInteractor* iren = vtkRenderWindowInteractor::New();
    iren->SetRenderWindow(renWin); 
    vtkInteractorStyleTrackballCamera* style=vtkInteractorStyleTrackballCamera::New(); 
    iren->SetInteractorStyle(style);
    style->Delete();
    iren->Initialize(); // register interactor to pick up windows messages

    renderer = vtkRenderer::New();
    renderer->SetBackground(1,1,1); // white

    vtkCamera* camera = vtkCamera::New();
    camera->SetPosition(0, .1*cameraScale, cameraScale);
    camera->SetFocalPoint(0,0,0);
    camera->ComputeViewPlaneNormal();
    camera->SetViewUp(0,1,0);
    renderer->SetActiveCamera(camera);
    camera->Delete();

    // Use modern renderer base camera lights in preference to deprecated interactor version
    iren->LightFollowCameraOn();
    renderer->LightFollowCameraOn();

    vtkLight* light = vtkLight::New();
    light->SetLightTypeToCameraLight();
    light->SetPosition(-1,0,0);
    light->SetFocalPoint(0,0,0);
    light->SetColor(1,1,1);
    light->SetIntensity(.75);
    renderer->AddLight(light);
    light->Delete();

    light = vtkLight::New();
    light->SetLightTypeToCameraLight();
    light->SetPosition(1,0,0);
    light->SetFocalPoint(0,0,0);
    light->SetColor(1,1,1);
    light->SetIntensity(.75);
    renderer->AddLight(light);
    light->Delete();

    light = vtkLight::New();
    light->SetLightTypeToCameraLight();
    light->SetPosition(0,1,1);
    light->SetFocalPoint(0,0,0);
    light->SetColor(1,1,1);
    light->SetIntensity(.75);
    renderer->AddLight(light);
    light->Delete();


    renWin->AddRenderer(renderer);

    renderer->ResetCamera();
    renWin->Render();
}

void VTKVisualizerRep::initTopology() {
    static bool hasInitialized = false;
    if (hasInitialized)
        return;
    SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS(mbs.systemTopologyHasBeenRealized(), "MultibodySystem", mbs.getName(), "VTKVisualizerRep::initTopology()");
    hasInitialized = true;
    const SimbodyMatterSubsystem& sbs = mbs.getMatterSubsystem();
    bodies.resize(sbs.getNBodies());
}

void VTKVisualizerRep::createInstanceGeometry(const State& state) {
    static bool hasCreated = false;
    if (hasCreated)
        return;
    hasCreated = true;
    
    // Mine the system for any geometry it wants us to show.

    std::vector<DecorativeGeometry> sysGeom;
    for (Stage stage = Stage::Topology; stage < Stage::Time; stage++)
        mbs.calcDecorativeGeometryAndAppend(state, stage, sysGeom);
    for (int i=0; i<(int)sysGeom.size(); ++i)
        addDecoration(MobilizedBodyIndex(sysGeom[i].getBodyId()), Transform(), sysGeom[i]);
}

void VTKVisualizerRep::report(const State& s) {
    if (!renWin) return;
    
    initTopology();

    mbs.realize(s, Stage::Position); // just in case
    createInstanceGeometry(s);

    const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();
    for (MobilizedBodyIndex i(1); i<matter.getNBodies(); ++i) {
        const Transform& config = matter.getMobilizedBody(i).getBodyTransform(s);
        setConfiguration(i, config);
    }
    for (int i=0; i<(int)dynamicGeom.size(); ++i) {
        const PerDynamicGeomInfo& info = dynamicGeom[i];
        const Transform& X_GB1 = 
            matter.getMobilizedBody(info.body1).getBodyTransform(s);
        const Transform& X_GB2 = 
            matter.getMobilizedBody(info.body2).getBodyTransform(s);
        setRubberBandLine(i, X_GB1*info.station1, X_GB2*info.station2);
    }

    for (int stage=Stage::Time; stage <= s.getSystemStage(); ++stage)
        mbs.calcDecorativeGeometryAndAppend(s, Stage::getValue(stage), 
                                            ephemeralGeometry);

    displayEphemeralGeometry(s);

    if (cameraNeedsToBeReset) {
        renderer->ResetCamera();
        cameraNeedsToBeReset = false;
    }

    renWin->Render();

   // removeEphemeralGeometry();

    // Process any window messages since last time
#ifdef _WIN32
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
#endif
    
#ifdef VTK_USE_CARBON
    EventRef theEvent;
    EventTargetRef theTarget;
    theTarget = GetEventDispatcherTarget();
    while (ReceiveNextEvent(0, NULL, kEventDurationNoWait, true, &theEvent) == noErr) {
        SendEventToEventTarget(theEvent, theTarget);
        ReleaseEvent(theEvent);
    }
#endif

#ifdef VTK_USE_X
    vtkXRenderWindowInteractor* rwi = dynamic_cast<vtkXRenderWindowInteractor*>(renWin->GetInteractor());
    while (XtAppPending(rwi->GetApp()))
        XtAppProcessEvent(rwi->GetApp(), XtIMAll);
#endif

}

void VTKVisualizerRep::zeroPointers() {
    renWin=0; renderer=0;
}
void VTKVisualizerRep::deletePointers() {
    // Delete all the actors. The geometry gets deleted automatically
    // thanks to good design!
    for (int i=0; i<(int)bodies.size(); ++i) {
        std::vector<vtkProp3D*>& actors = bodies[i].aList;
        for (int a=0; a<(int)actors.size(); ++a)
            actors[a]->Delete(), actors[a]=0;
    }
    for (int i=0; i<(int)dynamicGeom.size(); ++i) {
        dynamicGeom[i].actor->Delete();
        dynamicGeom[i].actor = 0;
    }
    for (int i=0; i < (int)ephemeralActors.size(); ++i) {
        renderer->RemoveActor(ephemeralActors[i]);
        ephemeralActors[i]->Delete();
        ephemeralActors[i] = 0;
    }

    if(renderer)renderer->Delete();
    if(renWin)renWin->Delete();
    zeroPointers();
}

void VTKVisualizerRep::setConfiguration(MobilizedBodyIndex bodyNum, const Transform& X_GB) {
    const std::vector<vtkProp3D*>& actors = bodies[bodyNum].aList;
    for (int i=0; i < (int)actors.size(); ++i) {
        vtkProp3D*       actor = actors[i];
        actor->SetPosition(X_GB.T()[0], X_GB.T()[1], X_GB.T()[2]);
        const Vec4 av = X_GB.R().convertRotationToAngleAxis();
        actor->SetOrientation(0,0,0);
        if (!actor->IsA("vtkFollower"))
            actor->RotateWXYZ(av[0]*RadiansToDegrees, av[1], av[2], av[3]);
    }
}

// Provide two points in ground frame and generate the appropriate line between them.
void VTKVisualizerRep::setRubberBandLine(int dgeom, const Vec3& p1, const Vec3& p2) {
    vtkActor*       actor = dynamicGeom[dgeom].actor;
    DecorativeLine& line  = dynamicGeom[dgeom].line;
    line.setEndpoints(p1, p2);

    VTKDecorativeGeometry vgeom;
    line.implementGeometry(vgeom);
    vtkPolyData* poly = vgeom.getVTKPolyData();

    vtkPolyDataMapper::SafeDownCast(actor->GetMapper())->SetInput(poly);
}

const vtkRenderer* VTKVisualizer::getVtkRenderer() const {
  return getRep().renderer;
}

vtkRenderer* VTKVisualizer::updVtkRenderer() {
  return updRep().renderer;
}

const vtkRenderWindow* VTKVisualizer::getVtkRenderWindow() const {
  return getRep().renWin;
}

vtkRenderWindow* VTKVisualizer::updVtkRenderWindow() {
  return updRep().renWin;
}
