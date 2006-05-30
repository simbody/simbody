/* Copyright (c) 2006 Stanford University and Michael Sherman.
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

/**@file
 * Test System-level functionality for Simbody.
 */

#include "SimTKcommon.h"
#include "Simmatrix.h"
#include "Simbody.h"

#include "simbody/internal/VTKReporter.h"

#include <string>
#include <iostream>
#include <exception>
#include <cmath>
using std::cout;
using std::endl;

using namespace SimTK;

static const Real Pi = std::acos(-1.);
static const Real RadiansToDegrees = 180/Pi;
static const Real DegreesToRadians = Pi/180;
static const int  GroundBody = 0;


class Measure {
public:
    Measure() : rep(0) { }
    ~Measure();
    Measure(const Measure&);
    Measure& operator=(const Measure&);

    void realize(const State&, Stage) const;

    /// Stage at which this Measure can be evaluated.
    Stage  getStage()        const;

    /// This is the index number within the owning Subsystem. It is 
    /// guaranteed to be a small integer index suitable for very
    /// fast access to the Measure.
    int    getIndex()    const;


    /// Is this handle the owner of this rep? This is true if the
    /// handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // Internal use only
    explicit Measure(class MeasureRep* r) : rep(r) { }
    bool                hasRep() const {return rep!=0;}
    const MeasureRep& getRep() const {assert(rep); return *rep;}
    MeasureRep&       updRep() const {assert(rep); return *rep;}
protected:
    class MeasureRep* rep;
};

class RealMeasure : public Measure {
public:
    RealMeasure();

    // virtual
    const Real& getValue(const State&) const;
};


int main() {
    try {


    // A System is a set of *interrelated* Subsystems. That is, the subsystems
    // contain objects whose definitions depend on objects in other subsystems.
    // To make that work, the subsystems have to be associated with the system
    // before they are fully constructed.

    // A SimbodySystem (a kind of MultibodySystem) pre-defines empty Subsystems
    // of particular types.
    SimbodySystem sys;

    // Get access to the subsystems within sys so that we can fill them with
    // goodies.

    SimbodyMechanics&          myPend   = sys.updSimbodyMechanicsSubsystem();
    SimbodyAnalyticGeometry&   geometry = sys.updAnalyticGeometrySubsystem();
    SimbodyMassElement&        masses   = sys.updMassPropertiesSubsystem();
    SimbodyForces&             forces   = sys.updSimbodyForcesSubsystem();
    SimbodyMeasures&           measures = sys.updSimbodyMeasuresSubsystem();
    SimbodyVisualization&      viz      = sys.updVisualizationSubsystem();

    // System-global variables.
    RealVariable& density = sys.addRealVariable(Stage::Parametrized).setDefault(1);

    int body1 = myPend.addRigidBody(GroundBody, Mobility::Pin);

    forces.addMobilitySpring(body1, 0, 10., 10*DegreesToRadians); // body, dof, k, l0
    forces.addGravity(Vec3(0,-9.8,0));





// TODO: I have a mobility that depends on a curve fixed to a body.
// Mobility is expressed as allowable movement of child 
// reference frame J with respect to parent frame Jb. The
// generalized speeds are defined via geometry fixed to J and Jb,
// via H(q) = partial(V_JbJ_G)/partial(u). The generalized coordinates
// are defined in some convenient way such that qdot=Qu for some matrix
// Q(q).
// Anything with a frame fixed to a body can be used to express mobility.
// So for example we can use an analytic curve on one body and the origin
// of a frame on another to express a point-on-curve 4 dof joint. Generalized
// speeds would be angular velocity and linear velocity of the follower point
// along the instantaneous curve tangent. The generalized coordinates would
// be quaternions and an arc-length-like scalar, although it wouldn't have
// to be uniform as long as we can calculate its derivative from the speed.
// A 2 dof bead-on-wire would have one speed like the above, and the other
// a rotation rate about the tangent as an axis.

    // Build topology.
    Body& ground = myPend.updGround();
    Body& body1  = myPend.addRigidBody(ground, Mobility::AxialRotation).setName("body1");
    Body& body1a = myPend.addRigidBody(ground, Mobility::AxialTranslation);
    Body& body1b = myPend.addRigidBody(ground, Mobility::Screw);
    Body& body2  = myPend.addRigidBody(body1,  Mobility::Orientation);
    Body& body3  = myPend.addRigidBody(body2,  Mobility::Translation);
    Body& ball   = myPend.addRigidBody(ground, Mobility::Free);
    Body& bead = myPend.addRigidBody();

    AnalyticCurve& curveOnBody3 = geometry.addCurve(body3, spline9);
    AnalyticFrame& beadFrame = geometry.addFrame(bead, Transform());
    bead.setMobility( BeadOnWire(curveOnBody3, beadFrame) );

    body1.setMassProperties(MassProperties());  // these can use any source of the right type
    body1.setJointFrame(Transform());           // these all return a reference to the body for convenience
    body1.setParentJointFrame(Transform());

    AnalyticGeometry& sphere = geometry.addSphere(body1).setRadius(1);
    masses.addUniformMassElement(sphere).setDensity(1);

    AnalyticGeometry& groundSphere = geometry.addSphere(ground, Vec3(1,2,3));

    RealMeasure& dist = measures.addDistanceMeasure(ground.getOrigin(), sphere.getCenter());

    forces.addContactElement(sphere, groundSphere);
    forces.addSpring(Station(body1, Vec3(0)), Station(ground, Vec3(1,1,1)), k, len0);
    forces.setGravity(Vec3(0,-9.8,0));

    viz.addDefaultMechanicalVisualization(myPend);

    State s;
    sys.realize(s, Stage::Modeled);
    density.set(state, 3.);
    sys.realize(s, Stage::Parametrized);
    body1.setQ(state, 30*DegreesToRadians);
    sys.realize(s, Stage::Configured);
    cout << "COM=" << body1.getStationLocation(state, body1.getCOMStation(state)) << endl;

    VTKReporter vtk(sys);
    vtk.report(state);


    }
    catch (const std::exception& e) {
        printf("EXCEPTION THROWN: %s\n", e.what());
    }
    return 0;
}
