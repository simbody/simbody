#include "SimTKsimbody.h"
#include "SimTKcommon/Testing.h"

#include <cstdio>
#include <iostream>
#include <exception>

using namespace SimTK;
using std::cout;using std::endl;


static const Real smallRad = 5*1.25/100;
static const Real largeRad = 10.5/100;
static const MobodIndex smallBodx = MobodIndex(1);
static const MobodIndex largeBodx = MobodIndex(0);
static const Vec3 smallPos(-.28,0,0);
static const Vec3 largePos(-.28,0,0);
static const Real hSmall = .03; // 6cm=100% strain
static const Real hLarge = .03;

class ForceReporter : public PeriodicEventReporter {
public:
    ForceReporter(const MultibodySystem& system, 
                  const CompliantContactSubsystem& complCont,
                  Real reportInterval)
    :   PeriodicEventReporter(reportInterval), m_system(system),
        m_compliant(complCont)
    {}

    ~ForceReporter() {}

    void handleEvent(const State& state) const override {
        m_system.realize(state, Stage::Dynamics);
        //cout << state.getTime() << ": E = " << m_system.calcEnergy(state)
        //     << " Ediss=" << m_compliant.getDissipatedEnergy(state)
        //     << " E+Ediss=" << m_system.calcEnergy(state)
        //                       +m_compliant.getDissipatedEnergy(state)
        //     << endl;
        const int ncont = m_compliant.getNumContactForces(state);
        //cout << "Num contacts: " << m_compliant.getNumContactForces(state) << endl;
        if (ncont==0)
            printf("%g 0.0 0.0 0.0\n", state.getTime());

        const SimbodyMatterSubsystem& matter = m_system.getMatterSubsystem();
        const MobilizedBody& smallBod=matter.getMobilizedBody(smallBodx);
        const MobilizedBody& largeBod=matter.getMobilizedBody(largeBodx);
        const Vec3 smallCtr = smallBod.findStationLocationInGround(state, smallPos);
        const Vec3 largeCtr = largeBod.findStationLocationInGround(state, largePos);
        const Real d = (smallRad+largeRad)-(smallCtr-largeCtr).norm();
        
        for (int i=0; i < ncont; ++i) {
            const ContactForce& force = m_compliant.getContactForce(state,i);
            const ContactId     id    = force.getContactId();
            cout << state.getTime() 
                 << " " << force.getForceOnSurface2()[1][1] // Normal
                 << " " << force.getForceOnSurface2()[1][0] // Tangential
                 << " " << d << "\n"; // penetration distance

            //ContactPatch patch;
            //const bool found = m_compliant.calcContactPatchDetailsById(state,id,patch);
            //cout << "patch for id" << id << " found=" << found << endl;
            //cout << "resultant=" << patch.getContactForce() << endl;
            //cout << "num details=" << patch.getNumDetails() << endl;
            //Real patchArea=0, maxDepth=0, maxf=0;
            //Vec3 ftot(0);
            //for (int i=0; i < patch.getNumDetails(); ++i) {
            //    const ContactDetail& detail = patch.getContactDetail(i);
            //    patchArea += detail.getPatchArea();
            //    maxDepth = std::max(maxDepth, detail.getDeformation());
            //    ftot += detail.getForceOnSurface2();
            //    maxf = std::max(maxf, detail.getForceOnSurface2().norm());
            //}
            //cout << "patchArea=" << patchArea 
            //     << " ftot=" << ftot << " maxDepth=" << maxDepth << "\n";
            //cout << "maxf=" << maxf << "\n";
        }
    }
private:
    const MultibodySystem&           m_system;
    const CompliantContactSubsystem& m_compliant;
};

// Nylon
static const Real nylon_density = 1100.;  // kg/m^3
static const Real nylon_young   = .05*2.5e9;  // pascals (N/m)
static const Real nylon_poisson = 0.4;    // ratio
static const Real nylon_planestrain = 
    ContactMaterial::calcPlaneStrainStiffness(nylon_young, nylon_poisson);
static const Real nylon_confined =
    ContactMaterial::calcConfinedCompressionStiffness(nylon_young, nylon_poisson);
static const Real nylon_dissipation = 10*0.005;

int main() {
  try
  { // Create the system.
    
    MultibodySystem         system; system.setUseUniformBackground(true);
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::Gravity          gravity(forces, matter, -YAxis, 10);
    //Force::GlobalDamper     damper(forces, matter, 1.0);

    ContactTrackerSubsystem  tracker(system);
    CompliantContactSubsystem contactForces(system, tracker);
    contactForces.setTrackDissipatedEnergy(true);
    contactForces.setTransitionVelocity(1e-3);

    // g=10, mass==.5 => weight = .5*10=5N.
    // body origin at the hinge on the right
    Real mass = 0.5;
    Vec3 com = Vec3(-.14,0,0);
    Vec3 halfDims(.14, .01, .01); // a rectangular solid
    Inertia inertia = mass*
        UnitInertia::brick(halfDims).shiftFromCentroid(Vec3(.14,0,0));

    Body::Rigid lidBody(MassProperties(mass,com,inertia));
    lidBody.addDecoration(Vec3(-.14,0,0), 
        DecorativeBrick(halfDims).setColor(Cyan).setOpacity(.1));

    PolygonalMesh smallMesh, largeMesh;

#define USE_VIZ
//#define SHOW_NORMALS
Array_<DecorativeLine> smallNormals;
Array_<DecorativeLine> largeNormals;
const Real NormalLength = .001;

#define USE_MESH_SMALL
#define USE_MESH_BIG
#ifdef USE_MESH_SMALL
    smallMesh = PolygonalMesh::createSphereMesh(smallRad, 5);
    std::cerr << "small mesh faces: " << smallMesh.getNumFaces() << "\n";
    ContactGeometry::TriangleMesh smallContact(smallMesh);
    DecorativeMesh smallArtwork(smallContact.createPolygonalMesh());
    smallArtwork.setColor(Red);
    #ifdef SHOW_NORMALS
    for (int fx=0; fx < smallContact.getNumFaces(); ++fx) {
        smallNormals.push_back(
        DecorativeLine(smallContact.findCentroid(fx),
                       smallContact.findCentroid(fx)
                           + NormalLength*smallContact.getFaceNormal(fx))
                           .setColor(Black).setLineThickness(2));
    }
    #endif
#else
    DecorativeSphere smallArtwork(smallRad);
    smallArtwork.setResolution(2).setRepresentation(DecorativeGeometry::DrawWireframe);
    ContactGeometry::Sphere smallContact(smallRad);
#endif
#ifdef USE_MESH_BIG
    largeMesh = PolygonalMesh::createSphereMesh(largeRad, 6);
    std::cerr << "large mesh faces: " << largeMesh.getNumFaces() << "\n";
    ContactGeometry::TriangleMesh largeContact(largeMesh);
    DecorativeMesh largeArtwork(largeContact.createPolygonalMesh());
    //DecorativeSphere largeArtwork(largeRad);
    largeArtwork.setColor(Green);
    #ifdef SHOW_NORMALS

    for (int fx=0; fx < largeContact.getNumFaces(); ++fx) {
        largeNormals.push_back(
        DecorativeLine(largeContact.findCentroid(fx),
                       largeContact.findCentroid(fx)
                           + NormalLength*largeContact.getFaceNormal(fx))
                           .setColor(Black).setLineThickness(2));
    }
    #endif
#else
    DecorativeSphere largeArtwork(largeRad);
    largeArtwork.setResolution(4).setRepresentation(DecorativeGeometry::DrawWireframe);
    ContactGeometry::Sphere largeContact(largeRad);
#endif


    ContactMaterial nylon(
        //nylon_planestrain, 
        nylon_confined, 
        nylon_dissipation, 
        .1*.9,.1*.8,.1*.6); // static, dynamic, viscous friction

    lidBody.addDecoration(smallPos, 
        smallArtwork.setOpacity(0.5));
    lidBody.addDecoration(smallPos, 
        smallArtwork.setRepresentation(DecorativeGeometry::DrawWireframe));

    for (unsigned i=0; i < smallNormals.size(); ++i)
        lidBody.addDecoration(smallPos, smallNormals[i]);

    lidBody.addContactSurface(smallPos, 
        ContactSurface(smallContact, nylon, hSmall));

    matter.Ground().updBody().addDecoration(largePos, 
        largeArtwork.setOpacity(0.5));
    matter.Ground().updBody().addDecoration(largePos, 
        largeArtwork.setRepresentation(DecorativeGeometry::DrawWireframe));

    for (unsigned i=0; i < largeNormals.size(); ++i)
        matter.Ground().updBody().addDecoration(largePos, 
        largeNormals[i]);
    matter.Ground().updBody().addContactSurface(largePos, 
        ContactSurface(largeContact, nylon, hLarge));

    MobilizedBody::Pin lid(matter.Ground(), Vec3(0.05,smallRad+largeRad,0), 
                           lidBody,         Vec3(0,0,0));

#ifdef USE_VIZ
    Visualizer viz(system);
    viz.setCameraClippingPlanes(.01,100);
#endif

    system.adoptEventReporter(new Visualizer::Reporter(viz, 1./1000));
    ForceReporter* frcReporter = 
        new ForceReporter(system, contactForces, 1./100000);
    system.adoptEventReporter(frcReporter);

    RungeKuttaMersonIntegrator integ(system);
    //CPodesIntegrator integ(system);
    //integ.setMaximumStepSize(.01);
    integ.setAccuracy(1e-5);
    TimeStepper ts(integ);

    system.realizeTopology();
    State startState = system.getDefaultState();
    lid.setOneQ(startState, 0, -Pi/10);
    //lid.setOneQ(startState, 0, 0.0610351);
    system.realize(startState, Stage::Position);
#ifdef USE_VIZ
    viz.report(startState);
#endif
    system.realize(startState);
    //frcReporter->handleEvent(startState);
    //printf("type to go ...\n");
   // getchar();
    ts.initialize(startState);
    ts.stepTo(0.5);

  } catch (const std::exception& e) {
    std::printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);

  } catch (...) {
    std::printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }

    return 0;
}

