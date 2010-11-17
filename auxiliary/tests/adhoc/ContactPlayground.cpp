/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2010 Stanford University and the Authors.           *
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

/**@file
 * Adhoc main program for playing with contact.
 */

#include "SimTKsimbody.h"
#include "SimTKsimbody_aux.h"   // requires VTK

#include <cstdio>
#include <exception>
#include <algorithm>
#include <iostream>
#include <ctime>
using std::cout; using std::endl;

using namespace SimTK;

Array_<State> saveEm;

static const Real ReportInterval = 0.01;
static const Real ForceScale = .25;
static const Real MomentScale = .5;

class MyReporter : public PeriodicEventReporter {
public:
    MyReporter(const MultibodySystem& system, 
               const CompliantContactSubsystem& complCont,
               VTKVisualizer& viz,
               Real reportInterval)
    :   PeriodicEventReporter(reportInterval), m_system(system),
        m_compliant(complCont), m_viz(viz) {}

    ~MyReporter() {}

    void generateForceArrows(const State& state) const {
        const Vec3 frcColors[] = {Red,Orange,Yellow};
        const Vec3 momColors[] = {Blue,Green,Purple};
        m_system.realize(state, Stage::Velocity);

        const int ncont = m_compliant.getNumContactForces(state);
        for (int i=0; i < ncont; ++i) {
            const ContactForce& force = m_compliant.getContactForce(state,i);
            const ContactId     id    = force.getContactId();
            const Vec3& frc = force.getForceOnSurface2()[1];
            const Vec3& mom = force.getForceOnSurface2()[0];
            Real  frcMag = frc.norm(), momMag=mom.norm();
            int frcThickness = 1, momThickness = 1;
            Real frcScale = ForceScale, momScale = ForceScale;
            while (frcMag > 10)
                frcThickness++, frcScale /= 10, frcMag /= 10;
            while (momMag > 10)
                momThickness++, momScale /= 10, momMag /= 10;
            DecorativeLine frcLine(force.getContactPoint(),
                force.getContactPoint() + frcScale*frc);
            DecorativeLine momLine(force.getContactPoint(),
                force.getContactPoint() + momScale*mom);
            frcLine.setColor(frcColors[id%3]);
            momLine.setColor(momColors[id%3]);
            frcLine.setLineThickness(2*frcThickness);
            momLine.setLineThickness(2*momThickness);
            m_viz.addEphemeralDecoration(frcLine);
            m_viz.addEphemeralDecoration(momLine);

            ContactPatch patch;
            const bool found = m_compliant.calcContactPatchDetailsById(state,id,patch);
            //cout << "patch for id" << id << " found=" << found << endl;
            //cout << "resultant=" << patch.getContactForce() << endl;
            //cout << "num details=" << patch.getNumDetails() << endl;
            for (int i=0; i < patch.getNumDetails(); ++i) {
                const ContactDetail& detail = patch.getContactDetail(i);
                const Real peakPressure = detail.getPeakPressure();
                // Make a black line from the element's contact point in the normal
                // direction, with length proportional to log(peak pressure)
                // on that element. 
                DecorativeLine normal(detail.getContactPoint(),
                    detail.getContactPoint()+ std::log10(peakPressure)
                                                * detail.getContactNormal());
                normal.setColor(Black);
                m_viz.addEphemeralDecoration(normal);
                // Make a red line that extends from the contact
                // point in the direction of the slip velocity, of length 3*slipvel.
                DecorativeLine slip(detail.getContactPoint(),
                    detail.getContactPoint()+3*detail.getSlipVelocity());
                slip.setColor(Red);
                m_viz.addEphemeralDecoration(slip);
            }
        }
    }

    void handleEvent(const State& state) const {
        m_system.realize(state, Stage::Dynamics);
        cout << state.getTime() << ": E = " << m_system.calcEnergy(state)
             << " Ediss=" << m_compliant.getDissipatedEnergy(state)
             << " E+Ediss=" << m_system.calcEnergy(state)
                               +m_compliant.getDissipatedEnergy(state)
             << endl;
        const int ncont = m_compliant.getNumContactForces(state);
        cout << "Num contacts: " << m_compliant.getNumContactForces(state) << endl;
        for (int i=0; i < ncont; ++i) {
            const ContactForce& force = m_compliant.getContactForce(state,i);
            //cout << force;
        }
        generateForceArrows(state);
        saveEm.push_back(state);
    }
private:
    const MultibodySystem&           m_system;
    const CompliantContactSubsystem& m_compliant;
    VTKVisualizer&                   m_viz;
};

static void makeCube(Real h, PolygonalMesh& cube);
static void makeTetrahedron(Real r, PolygonalMesh& tet);
static void makePyramid(Real baseSideLength, PolygonalMesh& pyramid);
static void makeOctahedron(Real radius, PolygonalMesh& pyramid);
static void makeSphere(Real radius, int level, PolygonalMesh& sphere);

int main() {
  try
  { // Create the system.
    
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::Gravity   gravity(forces, matter, UnitVec3(2,-10,0), 1);

    ContactTrackerSubsystem  tracker(system);
    CompliantContactSubsystem contactForces(system, tracker);
    contactForces.setTransitionVelocity(1e-3);

    GeneralContactSubsystem OLDcontact(system);
    const ContactSetIndex OLDcontactSet = OLDcontact.createContactSet();


    PolygonalMesh pyramidMesh;
    //makeCube(1, pyramidMesh);
    //makeTetrahedron(1, pyramidMesh);
    //pyramidMesh.transformMesh(Rotation(Pi/4, UnitVec3(-1,0,1)));
    //makePyramid(1, pyramidMesh);
    //makeOctahedron(1, pyramidMesh);
    makeSphere(1, 4, pyramidMesh);

    ContactGeometry::TriangleMesh pyramid(pyramidMesh);
    DecorativeMesh showPyramid(pyramid.createPolygonalMesh());
    Array_<DecorativeLine> normals;
    const Real NormalLength = .02;
    for (int fx=0; fx < pyramid.getNumFaces(); ++fx)
        normals.push_back(
        DecorativeLine(pyramid.findCentroid(fx),
                       pyramid.findCentroid(fx)
                           + NormalLength*pyramid.getFaceNormal(fx)));
    // not displaying mesh normals at the moment

    ContactCliqueId clique1 = ContactSurface::createNewContactClique();
    ContactCliqueId clique2 = ContactSurface::createNewContactClique();
    ContactCliqueId clique3 = ContactSurface::createNewContactClique();

    const Real fFac =1; // to turn off friction
    const Real fDis = .5*0.2; // to turn off dissipation
    const Real fVis =  .1*.1; // to turn off viscous friction
    const Real fK = 100*1e6; // pascals

    // Right hand wall
    matter.Ground().updBody().addDecoration(Vec3(.25+.01,0,0),
        DecorativeBrick(Vec3(.01,4,8)).setColor(Gray).setOpacity(.1));
    matter.Ground().updBody().addContactSurface(Vec3(.25,0,0),
        ContactSurface(ContactGeometry::HalfSpace(),
                       ContactMaterial(fK*.01,fDis*.9,fFac*.8,fFac*.7,fVis*10))
                       .joinClique(clique3));

    // Halfspace floor
    const Rotation R_xdown(-Pi/2,ZAxis);
    matter.Ground().updBody().addDecoration(
        Transform(R_xdown, Vec3(0,-3-.01,0)),
        DecorativeBrick(Vec3(.01,4,8)).setColor(Gray).setOpacity(.1));
    matter.Ground().updBody().addContactSurface(
        Transform(R_xdown, Vec3(0,-3,0)),
        ContactSurface(ContactGeometry::HalfSpace(),
                       ContactMaterial(fK*.1,fDis*.9,fFac*.8,fFac*.7,fVis*10))
                       //ContactMaterial(2e6,.01,.1,.05,.01))
                       .joinClique(clique1));

    //// Big Sphere floor
    //const Real FloorRadius = 10;
    //matter.Ground().updBody().addDecoration(
    //    Vec3(0,-FloorRadius-3,0),
    //    DecorativeSphere(FloorRadius).setColor(Green));
    //matter.Ground().updBody().addContactSurface(
    //    Vec3(0,-FloorRadius-3,0),
    //    ContactSurface(ContactGeometry::Sphere(FloorRadius),
    //                   ContactMaterial(1e6,fDis*.9,fFac*.8,fFac*.7,fVis*10))
    //                   //ContactMaterial(2e6,.01,.1,.05,.01))
    //                   .joinClique(clique1));

    const Real rad = .4;
    Body::Rigid pendulumBody1(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody1.addDecoration(Transform(), 
        DecorativeSphere(rad).setOpacity(.2));
    pendulumBody1.addContactSurface(Transform(),
        ContactSurface(ContactGeometry::Sphere(rad),
                       ContactMaterial(fK*.01,fDis*.9,fFac*.8,fFac*.7,fVis*10))
                       .joinClique(clique2));

    Body::Rigid pendulumBody2(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody2.addDecoration(Transform(), 
        DecorativeSphere(rad).setColor(Orange).setOpacity(.2));
    pendulumBody2.addContactSurface(Transform(),
        ContactSurface(ContactGeometry::Sphere(rad),
                       ContactMaterial(fK*.01,fDis*.9,fFac*.8,fFac*.7,fVis*10))
                       .joinClique(clique1));

    MobilizedBody::Pin pendulum(matter.Ground(), Transform(Vec3(0)), 
                                pendulumBody1,    Transform(Vec3(0, 1, 0)));

    MobilizedBody::Pin pendulum2(pendulum, Transform(Vec3(0)), 
                                 pendulumBody2, Transform(Vec3(0, 1, 0)));

    Body::Rigid pendulumBody3(MassProperties(100.0, Vec3(0), 100*Inertia(1)));
    PolygonalMesh body3contact;
    makeSphere(rad, 2, body3contact);
    ContactGeometry::TriangleMesh geo3(body3contact);
    const DecorativeMesh mesh3(geo3.createPolygonalMesh());
    pendulumBody3.addDecoration(Transform(), 
        DecorativeMesh(mesh3).setOpacity(.2));
    pendulumBody3.addDecoration(Transform(), 
        DecorativeMesh(mesh3).setColor(Gray)
                   .setRepresentation(DecorativeGeometry::DrawWireframe)
                   .setOpacity(.1));

    //ContactGeometry::Sphere geo3(rad);
    pendulumBody3.addContactSurface(Transform(),
        ContactSurface(geo3,
                       ContactMaterial(fK*.1,fDis*.9,fFac*.8,fFac*.7,fVis*10))
                       .joinClique(clique2));
    MobilizedBody::Pin pendulum3(matter.Ground(), Transform(Vec3(-2,0,0)), 
                                 pendulumBody3, Transform(Vec3(0, 2, 0)));

    Force::MobilityLinearSpring(forces, pendulum2, MobilizerUIndex(0),
        10, 0*(Pi/180));

    const Real ballMass = 200;
    Body::Rigid ballBody(MassProperties(ballMass, Vec3(0), 
                            ballMass*Gyration::sphere(1)));
    //ballBody.addDecoration(Transform(), DecorativeSphere(.3).setColor(Cyan));
    //ballBody.addContactSurface(Transform(),
    //    ContactSurface(ContactGeometry::Sphere(.3),
    //                   ContactMaterial(1e7,.05,fFac*.8,fFac*.7,fVis*10))
    //                   .joinClique(clique2));
    ballBody.addDecoration(Transform(), 
        showPyramid.setColor(Cyan).setOpacity(.2));
    ballBody.addDecoration(Transform(), 
        showPyramid.setColor(Gray)
                   .setRepresentation(DecorativeGeometry::DrawWireframe));
    //for (unsigned i=0; i < normals.size(); ++i)
    //    ballBody.addDecoration(Transform(),
    //        normals[i].setColor(Gray));
    ballBody.addDecoration(Transform(), DecorativeSphere(1).setColor(Gray)
                                             .setOpacity(.1).setResolution(10));
    ballBody.addContactSurface(Transform(),
        ContactSurface(pyramid,
                       ContactMaterial(fK*.1,fDis*.9,
                                       .1*fFac*.8,.1*fFac*.7,fVis*1))
                       //ContactMaterial(2e6,.01,.1,.05,.01))
                       //.joinClique(clique2)
                       );
    MobilizedBody::Free ball(matter.Ground(), Transform(Vec3(-2,0,0)),
        ballBody, Transform(Vec3(0)));

    //// The old way ...
    //OLDcontact.addBody(OLDcontactSet, ball,
    //    pyramid, Transform());

    //OLDcontact.addBody(OLDcontactSet, matter.updGround(),
    //    ContactGeometry::HalfSpace(), Transform(R_xdown, Vec3(0,-3,0)));
    //ElasticFoundationForce ef(forces, OLDcontact, OLDcontactSet);
    //Real stiffness = 1e6, dissipation = 0.01, us = 0.1, 
    //    ud = 0.05, uv = 0.01, vt = 0.01;
    ////Real stiffness = 1e6, dissipation = 0.1, us = 0.8, 
    ////    ud = 0.7, uv = 0.01, vt = 0.01;

    //ef.setBodyParameters(ContactSurfaceIndex(0), 
    //    stiffness, dissipation, us, ud, uv);
    //ef.setTransitionVelocity(vt);
    //// end of old way.


    VTKEventReporter& reporter = *new VTKEventReporter(system, ReportInterval);
    VTKVisualizer& viz = reporter.updVisualizer();
    
    MyReporter& myRep = *new MyReporter(system,contactForces,viz,ReportInterval);
    system.updDefaultSubsystem().addEventReporter(&myRep);
    system.updDefaultSubsystem().addEventReporter(&reporter);

    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    ball.setQToFitTransform(state, Transform(Rotation(Pi/2,XAxis),
                                             Vec3(0,-1.8,0)));

    pendulum.setOneQ(state, 0, -Pi/12);
    pendulum3.setOneQ(state, 0, -Pi/4);

    viz.report(state);
    printf("Default state\n");
    cout << "t=" << state.getTime() 
         << " q=" << pendulum.getQAsVector(state) << pendulum2.getQAsVector(state) 
         << " u=" << pendulum.getUAsVector(state) << pendulum2.getUAsVector(state) 
         << endl;
    cout << "Hit ENTER to simulate:\n";
    char c=getchar();


    pendulum.setOneU(state, 0, 5.0);
    ball.setOneU(state, 2, -20);

    ball.setOneU(state, 0, .05); // to break symmetry
    
    // Simulate it.
    const clock_t start = clock();


    //ExplicitEulerIntegrator integ(system);
    CPodesIntegrator integ(system,CPodes::BDF,CPodes::Newton);
    //RungeKuttaFeldbergIntegrator integ(system);
    //RungeKuttaMersonIntegrator integ(system);
    //RungeKutta3Integrator integ(system);
    //VerletIntegrator integ(system);
    //integ.setMaximumStepSize(1e-0001);
    integ.setAccuracy(1e-4);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(20.0);

    const double timeInSec = (double)(clock()-start)/CLOCKS_PER_SEC;
    const int evals = integ.getNumRealizations();
    cout << "Done -- took " << integ.getNumStepsTaken() << " steps in " <<
        timeInSec << "s for " << ts.getTime() << "s sim (avg step=" 
        << (1000*ts.getTime())/integ.getNumStepsTaken() << "ms) " 
        << (1000*ts.getTime())/evals << "ms/eval\n";

    printf("Using Integrator %s at accuracy %g:\n", 
        integ.getMethodName(), integ.getAccuracyInUse());
    printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), integ.getNumStepsAttempted());
    printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
    printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), integ.getNumProjections());


    while(true) {
        for (int i=0; i < (int)saveEm.size(); ++i) {
            //myRep.generateForceArrows(saveEm[i]);
            viz.report(saveEm[i]);
            //vtk.report(saveEm[i]); // half speed
        }
        cout << "Hit ENTER to animate again\n";
        getchar();
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

struct VertKey {
    VertKey(const Vec3& v) : v(v) {}
    Vec3 v;
    bool operator<(const VertKey& other) const {
        const Real tol = SignificantReal;
        const Vec3 diff = v - other.v;
        if (diff[0] < -tol) return true;
        if (diff[0] >  tol) return false;
        if (diff[1] < -tol) return true;
        if (diff[1] >  tol) return false;
        if (diff[2] < -tol) return true;
        if (diff[2] >  tol) return false;
        return false; // they are numerically equal
    }
};
typedef std::map<VertKey,int> VertMap;

/* Search a list of vertices for one close enough to this one and
return its index if found, otherwise add to the end. */
static int getVertex(const Vec3& v, VertMap& vmap, Array_<Vec3>& verts) {
    VertMap::const_iterator p = vmap.find(VertKey(v));
    if (p != vmap.end()) return p->second;
    const int ix = (int)verts.size();
    verts.push_back(v);
    vmap.insert(std::make_pair(VertKey(v),ix));
    return ix;
}

/* Each face comes in as below, with vertices 0,1,2 on the surface
of a sphere or radius r centered at the origin. We bisect the edges to get
points a',b',c', then move out from the center to make points a,b,c
on the sphere.
         1
        /\        
       /  \
    c /____\ b      Then construct new triangles
     /\    /\            [0,b,a]
    /  \  /  \           [a,b,c]
   /____\/____\          [c,2,a]
  2      a     0         [b,1,c]
*/
static void refineSphere(Real r, VertMap& vmap, 
                         Array_<Vec3>& verts, Array_<int>&  faces) {
    assert(faces.size() % 3 == 0);
    const int nVerts = faces.size(); // # face vertices on entry
    for (int i=0; i < nVerts; i+=3) {
        const int v0=faces[i], v1=faces[i+1], v2=faces[i+2];
        const Vec3 a = r*UnitVec3(verts[v0]+verts[v2]);
        const Vec3 b = r*UnitVec3(verts[v0]+verts[v1]);
        const Vec3 c = r*UnitVec3(verts[v1]+verts[v2]);
        const int va=getVertex(a,vmap,verts), 
                  vb=getVertex(b,vmap,verts), 
                  vc=getVertex(c,vmap,verts);
        // Replace the existing face with the 0ba triangle, then add the rest.
        // Refer to the above picture.
        faces[i+1] = vb; faces[i+2] = va;
        faces.push_back(va); faces.push_back(vb); faces.push_back(vc);//abc
        faces.push_back(vc); faces.push_back(v2); faces.push_back(va);//c2a
        faces.push_back(vb); faces.push_back(v1); faces.push_back(vc);//b1c
    }
}

// level  numfaces
//   0       8   <-- octahedron
//   1       32
//   2       128 <-- still lumpy
//   3       512 <-- very spherelike
//   n       2*4^(n+1)
static void makeSphere(Real radius, int level, PolygonalMesh& sphere) {
    Array_<Vec3> vertices;
    Array_<int> faceIndices;
    makeOctahedralMesh(Vec3(radius), vertices, faceIndices);

    VertMap vmap;
    for (unsigned i=0; i < vertices.size(); ++i)
        vmap[vertices[i]] = i;

    while (level > 0) {
        refineSphere(radius, vmap, vertices, faceIndices);
        --level;
    }

    for (unsigned i=0; i < vertices.size(); ++i)
        sphere.addVertex(vertices[i]);
    for (unsigned i=0; i < faceIndices.size(); i += 3) {
        const Array_<int> verts(&faceIndices[i], &faceIndices[i]+3);
        sphere.addFace(verts);
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
