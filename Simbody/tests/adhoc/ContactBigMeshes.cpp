/* -------------------------------------------------------------------------- *
 *                           SimTK Simbody(tm)                                *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2010-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
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

/**@file
 * Adhoc main program for testing contact behavior between two unreasonably
 * dense meshes.
 */

#include "SimTKsimbody.h"

#include <cstdio>
#include <exception>
#include <algorithm>
#include <iostream>
#include <fstream>
#include <ctime>
using std::cout; using std::endl;

using namespace SimTK;

Array_<State> saveEm;

static const Real ReportInterval = 0.01;
static const Real ForceScale = 10;
static const Real MomentScale = 20;

class ForceArrowGenerator : public DecorationGenerator {
public:
    ForceArrowGenerator(const MultibodySystem& system,
                        const CompliantContactSubsystem& complCont) 
    :   m_system(system), m_compliant(complCont) {}

    virtual void generateDecorations(const State& state, Array_<DecorativeGeometry>& geometry) {
        const Vec3 frcColors[] = {Red,Orange,Cyan};
        const Vec3 momColors[] = {Blue,Green,Purple};
        m_system.realize(state, Stage::Velocity);

        const int ncont = m_compliant.getNumContactForces(state);
        for (int i=0; i < ncont; ++i) {
            const ContactForce& force = m_compliant.getContactForce(state,i);
            const ContactId     id    = force.getContactId();
            printf("viz contact %d: id=%d\n", i, (int)id);
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
            geometry.push_back(frcLine);
            geometry.push_back(momLine);
        }
    }
private:
    const MultibodySystem&              m_system;
    const CompliantContactSubsystem&    m_compliant;
};

class MyReporter : public PeriodicEventReporter {
public:
    MyReporter(const MultibodySystem& system, 
               const CompliantContactSubsystem& complCont,
               Real reportInterval)
    :   PeriodicEventReporter(reportInterval), m_system(system),
        m_compliant(complCont) 
    {}

    ~MyReporter() {}

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

        saveEm.push_back(state);
    }
private:
    const MultibodySystem&           m_system;
    const CompliantContactSubsystem& m_compliant;
};

int main() {
  try
  { std::cout << "Current working directory: " 
              << Pathname::getCurrentWorkingDirectory() << std::endl;

    // Create the system.
    
    MultibodySystem         system; system.setUseUniformBackground(true);
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::UniformGravity   gravity(forces, matter, 0*Vec3(2, -9.8, 0));

    ContactTrackerSubsystem  tracker(system);
    CompliantContactSubsystem contactForces(system, tracker);
    contactForces.setTrackDissipatedEnergy(true);

    GeneralContactSubsystem OLDcontact(system);
    const ContactSetIndex OLDcontactSet = OLDcontact.createContactSet();

    contactForces.setTransitionVelocity(1e-3);

    std::ifstream meshFile1, meshFile2;
    PolygonalMesh femurMesh; 
    femurMesh.loadObjFile("ContactBigMeshes_Femur.obj"); 
    PolygonalMesh patellaMesh; 
    patellaMesh.loadObjFile("ContactBigMeshes_Patella.obj"); 

    ContactGeometry::TriangleMesh femurTri(femurMesh);
    ContactGeometry::TriangleMesh patellaTri(patellaMesh);

    DecorativeMesh showFemur(femurTri.createPolygonalMesh());
    Array_<DecorativeLine> femurNormals;
    const Real NormalLength = .02;
    //for (int fx=0; fx < femurTri.getNumFaces(); ++fx)
    //    femurNormals.push_back(
    //    DecorativeLine(femurTri.findCentroid(fx),
    //                   femurTri.findCentroid(fx)
    //                       + NormalLength*femurTri.getFaceNormal(fx)));

    DecorativeMesh showPatella(patellaTri.createPolygonalMesh());
    Array_<DecorativeLine> patellaNormals;
    //for (int fx=0; fx < patellaTri.getNumFaces(); ++fx)
    //    patellaNormals.push_back(
    //    DecorativeLine(patellaTri.findCentroid(fx),
    //                   patellaTri.findCentroid(fx)
    //                       + NormalLength*patellaTri.getFaceNormal(fx)));

    // This transform has the meshes close enough that their OBBs overlap
    // but in the end none of the faces are touching.
    const Transform X_FP(
        Rotation(Mat33( 0.97107625831404454, 0.23876955530133021, 0,
                       -0.23876955530133021, 0.97107625831404454, 0,
                        0,                   0,                   1), true),
        Vec3(0.057400580865008571, 0.43859170879135373, 
             -0.00016506240185135300)
        );


    const Real fFac =1; // to turn off friction
    const Real fDis = .5*0.2; // to turn off dissipation
    const Real fVis =  .1*.1; // to turn off viscous friction
    const Real fK = 100*1e6; // pascals

    // Put femur on ground at origin
    matter.Ground().updBody().addDecoration(Vec3(0,0,0),
        showFemur.setColor(Cyan).setOpacity(.2));
    matter.Ground().updBody().addContactSurface(Vec3(0,0,0),
        ContactSurface(femurTri,
            ContactMaterial(fK*.01,fDis*.9,fFac*.8,fFac*.7,fVis*10),
            .01 /*thickness*/));


    Body::Rigid patellaBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    patellaBody.addDecoration(Transform(), 
        showPatella.setColor(Red).setOpacity(.2));
    patellaBody.addContactSurface(Transform(),
        ContactSurface(patellaTri,
            ContactMaterial(fK*.001,fDis*.9,fFac*.8,fFac*.7,fVis*10),
            .01 /*thickness*/));

    MobilizedBody::Free patella(matter.Ground(), Transform(Vec3(0)), 
                                patellaBody,    Transform(Vec3(0)));


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

    Visualizer viz(system);
    Visualizer::Reporter& reporter = *new Visualizer::Reporter(viz, ReportInterval);
    viz.addDecorationGenerator(new ForceArrowGenerator(system,contactForces));
    MyReporter& myRep = *new MyReporter(system,contactForces,ReportInterval);

    system.adoptEventReporter(&myRep);
    system.adoptEventReporter(&reporter);

    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();

    viz.report(state);
    printf("Reference state -- hit ENTER\n");
    cout << "t=" << state.getTime() 
         << " q=" << patella.getQAsVector(state)  
         << " u=" << patella.getUAsVector(state)  
         << endl;
    char c=getchar();

    patella.setQToFitTransform(state, ~X_FP);
    viz.report(state);
    printf("Initial state -- hit ENTER\n");
    cout << "t=" << state.getTime() 
         << " q=" << patella.getQAsVector(state)  
         << " u=" << patella.getUAsVector(state)  
         << endl;
    c=getchar();
    
    // Simulate it.
    const clock_t start = clock();

    RungeKutta3Integrator integ(system);
    TimeStepper ts(integ);
    ts.initialize(state);
    ts.stepTo(2.0);

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
            viz.report(saveEm[i]);
        }
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

