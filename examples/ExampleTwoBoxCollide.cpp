/* -------------------------------------------------------------------------- *
 *                   Simbody(tm) Example: Contact Playground                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-13 Stanford University and the Authors.        *
 * Authors: Kevin He, Michael Sherman                                         *
 * Contributors:                                                              *
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

/* Kevin He at Roblox created this example starting with
ExampleContactPlayground. It basically demonstrates why an Elastic Foundation
(EF) contact model is not a good way to handle coarsely-meshed simple objects,
like a box. EF uses the centroid of each face to generate an area-weighted
force, which can be a good physical representation with a dense mesh but
since the vertices and edges don't participate there is a lot of visible
penetration for a coarse mesh like the ones here.
*/

#include "Simbody.h"

#include <cstdio>
#include <exception>
#include <algorithm>
#include <iostream>
#include <fstream>
using std::cout; using std::endl;

using namespace SimTK;

Array_<State> saveEm;

static const Real TimeScale = 1;
static const Real FrameRate = 30;
static const Real ReportInterval = TimeScale/FrameRate;

// These are the item numbers for the entries on the Run menu.
static const int RunMenuId = 3, HelpMenuId = 7;
static const int GoItem = 1, ReplayItem=2, QuitItem=3;

// This is a periodic event handler that interrupts the simulation on a regular
// basis to poll the InputSilo for user input. If there has been some, process it.
// This one does nothing but look for the Run->Quit selection.
class UserInputHandler : public PeriodicEventHandler {
public:
    UserInputHandler(Visualizer::InputSilo& silo, Real interval)
    :   PeriodicEventHandler(interval), m_silo(silo) {}

    virtual void handleEvent(State& state, Real accuracy,
                             bool& shouldTerminate) const override
    {
        int menuId, item;
        if (m_silo.takeMenuPick(menuId, item) && menuId==RunMenuId && item==QuitItem)
            shouldTerminate = true;
    }

private:
    Visualizer::InputSilo& m_silo;
};


static void makeCube(Real h, PolygonalMesh& cube);
static void makeTetrahedron(Real r, PolygonalMesh& tet);
static void makePyramid(Real baseSideLength, PolygonalMesh& pyramid);
static void makeOctahedron(Real radius, PolygonalMesh& pyramid);


int main() {
  try
  { // Create the system.

    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::Gravity   gravity(forces, matter, UnitVec3(2,-10,0), 1);

    ContactTrackerSubsystem  tracker(system);
    CompliantContactSubsystem contactForces(system, tracker);
    contactForces.setTrackDissipatedEnergy(true);
    contactForces.setTransitionVelocity(1e-3);

    Vec3 halfSize(3,4,5);
    Vec3 halfSize2(200, 4, 200);
    ContactGeometry::TriangleMesh box(PolygonalMesh::createBrickMesh(halfSize));
    ContactGeometry::TriangleMesh box2(PolygonalMesh::createBrickMesh(halfSize2, 20));
    DecorativeMesh showBox(box.createPolygonalMesh());
    DecorativeMesh showBox2(box2.createPolygonalMesh());

    const Real boxMass = halfSize[0] * halfSize[1] * halfSize[2] * 8;
    const Real boxMass2 = halfSize2[0] * halfSize2[1] * halfSize2[2] * 8;
    Body::Rigid boxBody(MassProperties(boxMass, Vec3(0),
                                        boxMass * UnitInertia::brick(halfSize)));
    Body::Rigid boxBody2(MassProperties(boxMass2, Vec3(0),
                        boxMass2 * UnitInertia::brick(halfSize2)));
    boxBody.addDecoration(Transform(),
                            showBox.setColor(Red).setOpacity(1));
    boxBody.addDecoration(Transform(),
                            showBox.setColor(Gray).setRepresentation(DecorativeGeometry::DrawWireframe));
    boxBody2.addDecoration(Transform(),
                            showBox2.setColor(Cyan).setOpacity(.6));
    boxBody2.addDecoration(Transform(),
                            showBox2.setColor(Gray).setRepresentation(DecorativeGeometry::DrawWireframe));
//     boxBody.addDecoration(Transform(),
//                           DecorativeSphere(1).setColor(Gray).setOpacity(.1).setResolution(10));

    const Real fFac = 0.3;       // to turn off friction
    const Real fDis = 0.1;    // to turn off dissipation
    const Real fVis = 0.01;    // to turn off viscous friction
    const Real fK = 1e+8; // pascals

    boxBody.addContactSurface(Transform(),
                               ContactSurface(box,
                                               ContactMaterial(fK, fDis, fFac, fFac, fVis),
                                               .5 /*thickness*/)
                                               );
    boxBody2.addContactSurface(Transform(),
                                ContactSurface(box2,
                                                ContactMaterial(fK, fDis, fFac, fFac, fVis),
                                                .5 /*thickness*/)
                                                );
    MobilizedBody::Free boxMBody(matter.Ground(), Transform(Vec3(0)), boxBody, Transform(Vec3(0)));

    MobilizedBody::Weld boxMBody2(matter.Ground(), Transform(Vec3(0)), boxBody2, Transform(Vec3(0)));

    Visualizer viz(system);
    // viz.addDecorationGenerator(new ForceArrowGenerator(system,contactForces));
    viz.setMode(Visualizer::RealTime);
    viz.setDesiredBufferLengthInSec(1);
    viz.setDesiredFrameRate(FrameRate);
    viz.setGroundHeight(0);
    viz.setShowShadows(true);

    Visualizer::InputSilo* silo = new Visualizer::InputSilo();
    viz.addInputListener(silo);
    Array_<std::pair<String,int> > runMenuItems;
    runMenuItems.push_back(std::make_pair("Go", GoItem));
    runMenuItems.push_back(std::make_pair("Replay", ReplayItem));
    runMenuItems.push_back(std::make_pair("Quit", QuitItem));
    viz.addMenu("Run", RunMenuId, runMenuItems);

    Array_<std::pair<String,int> > helpMenuItems;
    helpMenuItems.push_back(std::make_pair("TBD - Sorry!", 1));
    viz.addMenu("Help", HelpMenuId, helpMenuItems);

    // system.addEventReporter(new MyReporter(system,contactForces,ReportInterval));
    system.addEventReporter(new Visualizer::Reporter(viz, ReportInterval));

    // Check for a Run->Quit menu pick every 1/4 second.
    system.addEventHandler(new UserInputHandler(*silo, .25));

    // Initialize the system and state.

    system.realizeTopology();
    State state = system.getDefaultState();
    //ball.setQToFitTransform(state, Transform(Rotation(Pi/2,XAxis),
    //                                         Vec3(0,-1.8,0)));
    boxMBody.setQToFitTransform(state, Transform(Vec3(0, 10, 0)));
    boxMBody2.setQToFitTransform(state, Transform(Vec3(0, 0, 0)));

    viz.report(state);

    cout << "\nChoose 'Go' from Run menu to simulate:\n";
    int menuId, item;
    do { silo->waitForMenuPick(menuId, item);
         if (menuId != RunMenuId || item != GoItem)
             cout << "\aDude ... follow instructions!\n";
    } while (menuId != RunMenuId || item != GoItem);

    //ball.setOneU(state, 2, -20);

   // ball.setOneU(state, 0, .05); // to break symmetry

    CPodesIntegrator integ(system,CPodes::BDF,CPodes::Newton);
    integ.setAccuracy(1e-3); // minimum for CPodes
    TimeStepper ts(system, integ);

    ts.initialize(state);

    double cpuStart = cpuTime();
    double realStart = realTime();

    ts.stepTo(100.0);

    const double timeInSec = realTime() - realStart;
    const int evals = integ.getNumRealizations();
    cout << "Done -- took " << integ.getNumStepsTaken() << " steps in " <<
        timeInSec << "s elapsed for " << ts.getTime() << "s sim (avg step="
        << (1000*ts.getTime())/integ.getNumStepsTaken() << "ms) "
        << (1000*ts.getTime())/evals << "ms/eval\n";
    cout << "  CPU time was " << cpuTime() - cpuStart << "s\n";

    printf("Using Integrator %s at accuracy %g:\n",
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), integ.getNumProjections());

    viz.dumpStats(std::cout);

    // Add as slider to control playback speed.
    viz.addSlider("Speed", 1, 0, 4, 1);
    viz.setMode(Visualizer::PassThrough);

    silo->clear(); // forget earlier input
    double speed = 1; // will change if slider moves
    while(true) {
        cout << "Choose Run/Replay to see that again ...\n";

        int menuId, item;
        silo->waitForMenuPick(menuId, item);


        if (menuId != RunMenuId) {
            cout << "\aUse the Run menu!\n";
            continue;
        }

        if (item == QuitItem)
            break;
        if (item != ReplayItem) {
            cout << "\aHuh? Try again.\n";
            continue;
        }

        for (double i=0; i < (int)saveEm.size(); i += speed ) {
            int slider; Real newValue;
            if (silo->takeSliderMove(slider,newValue)) {
                speed = newValue;
            }
            viz.report(saveEm[(int)i]);
        }
    }

  } catch (const std::exception& e) {
    std::printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);

  } catch (...) {
    std::printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

    return 0;
}


// Create a triangle mesh in the shape of a pyramid, with the
// square base in the x-z plane centered at 0,0,0 of given side length s.
// The base is split into two triangles. The apex will be at (0,s,0).
static void makePyramid(Real s, PolygonalMesh& pyramidMesh) {
    const Real h = s/2;
    Array_<Vec3> vertices;
    vertices.push_back(Vec3(-h, 0, -h));     // base
    vertices.push_back(Vec3( h, 0, -h));
    vertices.push_back(Vec3( h, 0,  h));
    vertices.push_back(Vec3(-h, 0,  h));
    vertices.push_back(Vec3( 0, s,  0)); // apex
    Array_<int> faceIndices;
    int faces[6][3] = {{0, 1, 2}, {0, 2, 3}, {1, 0, 4},
                       {2, 1, 4}, {3, 2, 4}, {0, 3, 4}};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 3; j++)
            faceIndices.push_back(faces[i][j]);

    for (unsigned i=0; i < vertices.size(); ++i)
        pyramidMesh.addVertex(vertices[i]);
    for (unsigned i=0; i < faceIndices.size(); i += 3) {
        const Array_<int> verts(&faceIndices[i], &faceIndices[i]+3);
        pyramidMesh.addFace(verts);
    }
}


// Create a triangle mesh in the shape of a tetrahedron with the
// points in the corners of a cube inscribed in a sphere of radius r.
static void makeTetrahedron(Real r, PolygonalMesh& tet) {
    const Real h = r/std::sqrt(Real(3)); // half-dim of cube
    Array_<Vec3> vertices;
    vertices.push_back(Vec3( h, h,  h));
    vertices.push_back(Vec3(-h,-h,  h));
    vertices.push_back(Vec3(-h, h, -h));
    vertices.push_back(Vec3( h,-h, -h));
    Array_<int> faceIndices;
    int faces[4][3] = {{0, 2, 1}, {1, 3, 0}, {0, 3, 2}, {2, 3, 1}};
    for (int i = 0; i < 4; i++)
        for (int j = 0; j < 3; j++)
            faceIndices.push_back(faces[i][j]);

    for (unsigned i=0; i < vertices.size(); ++i)
        tet.addVertex(vertices[i]);
    for (unsigned i=0; i < faceIndices.size(); i += 3) {
        const Array_<int> verts(&faceIndices[i], &faceIndices[i]+3);
        tet.addFace(verts);
    }
}

static void makeOctahedralMesh(const Vec3& r, Array_<Vec3>& vertices,
                               Array_<int>&  faceIndices) {
    vertices.push_back(Vec3( r[0],  0,  0));   //0
    vertices.push_back(Vec3(-r[0],  0,  0));   //1
    vertices.push_back(Vec3( 0,  r[1],  0));   //2
    vertices.push_back(Vec3( 0, -r[1],  0));   //3
    vertices.push_back(Vec3( 0,  0,  r[2]));   //4
    vertices.push_back(Vec3( 0,  0, -r[2]));   //5
    int faces[8][3] = {{0, 2, 4}, {4, 2, 1}, {1, 2, 5}, {5, 2, 0},
                       {4, 3, 0}, {1, 3, 4}, {5, 3, 1}, {0, 3, 5}};
    for (int i = 0; i < 8; i++)
        for (int j = 0; j < 3; j++)
            faceIndices.push_back(faces[i][j]);
}

// Create a triangle mesh in the shape of an octahedron (like two
// pyramids stacked base-to-base, with the square base in the x-z plane
// centered at 0,0,0 of given "radius" r.
// The apexes will be at (0,+/-r,0).
static void makeOctahedron(Real r, PolygonalMesh& mesh) {
    Array_<Vec3> vertices;
    Array_<int> faceIndices;
    makeOctahedralMesh(Vec3(r), vertices, faceIndices);

    for (unsigned i=0; i < vertices.size(); ++i)
        mesh.addVertex(vertices[i]);
    for (unsigned i=0; i < faceIndices.size(); i += 3) {
        const Array_<int> verts(&faceIndices[i], &faceIndices[i]+3);
        mesh.addFace(verts);
    }
}

static void makeCube(Real h, PolygonalMesh& cube) {
    Array_<Vec3> vertices;
    vertices.push_back(Vec3( h, h,  h));
    vertices.push_back(Vec3( h, h, -h));
    vertices.push_back(Vec3( h,-h,  h));
    vertices.push_back(Vec3( h,-h, -h));
    vertices.push_back(Vec3(-h, h,  h));
    vertices.push_back(Vec3(-h, h, -h));
    vertices.push_back(Vec3(-h,-h,  h));
    vertices.push_back(Vec3(-h,-h, -h));

    Array_<int> faceIndices;
    int faces[6][4] = {{0,2,3,1},{1,5,4,0},{0,4,6,2},
                       {2,6,7,3},{3,7,5,1},{4,5,7,6}};
    for (int i = 0; i < 6; i++)
        for (int j = 0; j < 4; j++)
            faceIndices.push_back(faces[i][j]);

    for (unsigned i=0; i < vertices.size(); ++i)
        cube.addVertex(vertices[i]);
    for (unsigned i=0; i < faceIndices.size(); i += 4) {
        const Array_<int> verts(&faceIndices[i], &faceIndices[i]+4);
        cube.addFace(verts);
    }
}


