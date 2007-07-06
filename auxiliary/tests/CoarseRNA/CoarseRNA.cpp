/* Portions copyright (c) 2006 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 * A big, chunky fake RNA built of cylinders and ball joints.
 */

#include "SimTKsimbody.h"

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

using namespace std;
using namespace SimTK;

static const Real Pi = (Real)SimTK_PI, Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN;
static const int  GroundBodyNum = 0; // ground is always body 0

static const Real g = 9.8; // meters/s^2; apply in –y direction

static const Real DuplexRadius = 3; // A
static const Real HalfHeight = 10;  // A
static const Real CylinderSlop = 1; // A

static const int  NAtoms = 20;
static const Real AtomMass = 12;    // Daltons
static const Real AtomRadius = 1;   // A

static const Real ConnectorRadius     = 1;  // A
static const Real ConnectorHalfHeight = 3;  // A
static const Real ConnectorEndSlop    = 0.2;// A
static const Real ConnectorDensity    = 10;  // Dalton/A^3

static int NSegments = 3;

class MyRNAExample : public SimbodyMatterSubsystem {
    struct PerBodyInfo {
        PerBodyInfo(MobilizedBodyId b, bool d) : bnum(b), isDuplex(d) { }
        MobilizedBodyId bnum;
        bool isDuplex;
    };
    std::vector<PerBodyInfo> bodyInfo;
    MobilizedBodyId end1, end2;
public:
    MyRNAExample(MultibodySystem& mbs, int nsegs, bool shouldFlop) : SimbodyMatterSubsystem(mbs)
    {
        bodyInfo.push_back(PerBodyInfo(GroundId, false)); // placeholder for ground
        end1 = makeChain(GroundId, Vec3(0), nsegs, shouldFlop);
        end2 = makeChain(GroundId, Vec3(20,0,0), nsegs, shouldFlop);

        if (true) {
            Constraint::Rod theConstraint2(updMobilizedBody(end1), Vec3(0, -HalfHeight,0),
                                           updMobilizedBody(end2), Vec3(0, -HalfHeight,0), 10);
        }

    }

    void decorateBody(MobilizedBodyId bodyNum, VTKReporter& display) const {
        assert(bodyInfo[bodyNum].bnum == bodyNum);
        if (bodyInfo[bodyNum].isDuplex)
            addDuplexDecorations(bodyNum, DuplexRadius, HalfHeight, CylinderSlop, 
                                 NAtoms, AtomRadius, display);
        else 
            addConnectorDecorations(bodyNum, ConnectorRadius, ConnectorHalfHeight, 
                                    ConnectorEndSlop, display);
    }

    void decorateGlobal(VTKReporter& display) const {
        DecorativeLine rbProto; rbProto.setColor(Black).setLineThickness(2);
        display.addRubberBandLine(end1, Vec3(0, -HalfHeight,0), end2, Vec3(0, -HalfHeight,0), rbProto);
    }

private:

    MobilizedBodyId makeChain(MobilizedBodyId startBodyId, const Vec3& startOrigin, int nSegs, bool shouldFlop) {
        //MobilizedBodyId baseBodyId = startBody;
        MobilizedBody baseBody = updMobilizedBody(startBodyId);
        Vec3 origin = startOrigin;
        MobilizedBody lastDup;
        for (int seg=0; seg < nSegs; ++seg) {

            MobilizedBody::Ball left1(
                baseBody, Transform(origin + Vec3(-DuplexRadius,-HalfHeight,0)),
                Body::Rigid(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity)),
                Transform(Vec3(0, ConnectorHalfHeight, 0)));
            bodyInfo.push_back(PerBodyInfo(left1, false));

            MobilizedBody::Ball left2(
                left1, Transform(Vec3(0, -ConnectorHalfHeight, 0)),
                Body::Rigid(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity)),
                Transform(Vec3(0, ConnectorHalfHeight, 0)));
            bodyInfo.push_back(PerBodyInfo(left2, false));

            MobilizedBody::Ball rt1(
                baseBody, Transform(origin + Vec3(DuplexRadius,-HalfHeight,0)),
                Body::Rigid(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity)),
                Transform(Vec3(0, ConnectorHalfHeight, 0)));
            bodyInfo.push_back(PerBodyInfo(rt1, false));

            MobilizedBody::Ball rt2(                             
                rt1, Transform(Vec3(0, -ConnectorHalfHeight, 0)),
                Body::Rigid(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity)),
                Transform(Vec3(0, ConnectorHalfHeight, 0)));
            bodyInfo.push_back(PerBodyInfo(rt2, false));

            MobilizedBody::Ball dup(
                rt2, Transform(Vec3(0, -ConnectorHalfHeight, 0)),
                Body::Rigid(calcDuplexMassProps(DuplexRadius, HalfHeight, NAtoms, AtomMass)),
                                Transform(Vec3(-DuplexRadius, HalfHeight, 0)));
            bodyInfo.push_back(PerBodyInfo(dup, true));

            if (!shouldFlop) {
                Constraint::Ball theConstraint(left2, Vec3(0, -ConnectorHalfHeight, 0),
                                               dup, Vec3(DuplexRadius, HalfHeight, 0));
            }

            baseBody = lastDup = dup;
            origin = Vec3(0);
        }
        return lastDup;
    }

    MassProperties calcDuplexMassProps(
        Real halfHeight, Real r, int nAtoms, Real atomMass)
    {
        const Real pitch = 2*Pi/halfHeight;
        const Real trans = (2*halfHeight)/(nAtoms-1);
        const Real rot = pitch*trans;
        Inertia iner(0);
        Vec3 com(0);
        Real mass = 0;
        for (int i=0; i<nAtoms; ++i) {
            const Real h = halfHeight - i*trans;
            const Real th = i*rot;
            const Vec3 p1(-r*cos(th),h,r*sin(th)), p2(r*cos(th),h,-r*sin(th));
            mass += 2*atomMass;
            iner += Inertia(p1, atomMass) + Inertia(p2, atomMass);
            com += atomMass*p1 + atomMass*p2;
        }
        return MassProperties(mass,com/mass,iner);
    }

    MassProperties calcConnectorMassProps(Real r, Real halfHeight, Real density)
    {
        const Real volume = Pi*r*r*halfHeight;
        const Real mass = volume*density;
        const Vec3 com = Vec3(0);
        const Inertia iner = mass*Inertia::cylinderAlongY(r, halfHeight);

        return MassProperties(mass,com,iner);
    }

    void addDuplexDecorations(MobilizedBodyId bodyNum, Real r, Real halfHeight, Real slop, int nAtoms,
                              Real atomRadius, VTKReporter& display) const
    {
        display.addDecoration(bodyNum, Transform(), 
            DecorativeCylinder(r+atomRadius+slop, halfHeight).setColor(Cyan).setOpacity(0.4));

        const Real pitch = 2*Pi/halfHeight;
        const Real trans = (2*halfHeight)/(nAtoms-1);
        const Real rot = pitch*trans;
        for (int i=0; i<nAtoms; ++i) {
            const Real h = halfHeight - i*trans;
            const Real th = i*rot;
            const Vec3 p1(-r*cos(th),h,r*sin(th)), p2(r*cos(th),h,-r*sin(th));
            display.addDecoration(bodyNum, Transform(Vec3(p1)), 
                DecorativeSphere(atomRadius).setColor(Red).setResolution(0.5));
            display.addDecoration(bodyNum, Transform(Vec3(p2)), 
                DecorativeSphere(atomRadius).setColor(Green).setResolution(0.5));
        }
    }

    void addConnectorDecorations(MobilizedBodyId bodyNum, Real r, Real halfHeight, Real endSlop,  
                                 VTKReporter& display) const
    {
        display.addDecoration(bodyNum, Transform(), 
            DecorativeCylinder(r, halfHeight-endSlop).setColor(Blue));
    }
};


int main(int argc, char** argv) {
    std::vector<State> saveEm;

try // If anything goes wrong, an exception will be thrown.
  { int nseg = NSegments;
    int shouldFlop = 0;
    if (argc > 1) sscanf(argv[1], "%d", &nseg);
    if (argc > 2) sscanf(argv[2], "%d", &shouldFlop);

    // Create a multibody system using Simbody.
    MultibodySystem mbs;
    MyRNAExample myRNA(mbs, nseg, shouldFlop != 0);
    GeneralForceElements forces(mbs);
    UniformGravitySubsystem ugs(mbs, Vec3(0, -g, 0));

    const Vec3 attachPt(150, -40, -50);

    forces.addTwoPointLinearSpring(GroundId, attachPt,
                                   MobilizedBodyId(myRNA.getNBodies()-1), Vec3(0),
                                   1000.,  // stiffness
                                   1.);    // natural length

    /* forces.addTwoPointLinearSpring(0, -attachPt,
                                   myRNA.getNBodies()-1, Vec3(0),
                                   1000.,  // stiffness
                                   1.);    // natural length
    */

    forces.addGlobalEnergyDrain(1000);


    State s = mbs.realizeTopology();
    //myRNA.setUseEulerAngles(s,true);
    mbs.realizeModel(s);

    printf("# quaternions in use = %d\n", myRNA.getNQuaternionsInUse(s));
    for (MobilizedBodyId i(0); i<myRNA.getNBodies(); ++i) {
        printf("body %2d: using quat? %s; quat index=%d\n",
            (int)i, myRNA.isUsingQuaternion(s,i) ? "true":"false", 
            myRNA.getQuaternionIndex(s,i));
    }

    ugs.updGravity(s) *= 10;
    ugs.disableGravity(s);
    ugs.enableGravity(s);
    ugs.updZeroHeight(s) = -0.8;

    // And a study using the Runge Kutta Merson integrator
    bool suppressProject = false;
    RungeKuttaMerson myStudy(mbs, s, suppressProject);
    //CPodesIntegrator myStudy(mbs, s);
    myStudy.setAccuracy(1e-2);
    myStudy.setConstraintTolerance(1e-3); 
    myStudy.setProjectEveryStep(false);

    VTKReporter display(mbs);
    for (MobilizedBodyId i(1); i<myRNA.getNBodies(); ++i)
        myRNA.decorateBody(i, display);
    myRNA.decorateGlobal(display);

    DecorativeLine rbProto; rbProto.setColor(Orange).setLineThickness(3);
    display.addRubberBandLine(GroundId, attachPt,MobilizedBodyId(myRNA.getNBodies()-1),Vec3(0), rbProto);
    //display.addRubberBandLine(GroundId, -attachPt,myRNA.getNBodies()-1,Vec3(0), rbProto);

    const Real dt = 0.05; // output intervals

    printf("time  nextStepSize\n");

    s.updTime() = 0;
    for (int i=0; i<50; ++i)
        saveEm.push_back(s);    // delay
    display.report(s);

    myStudy.initialize();
    saveEm.push_back(s);
    for (int i=0; i<50; ++i)
        saveEm.push_back(s);    // delay
    display.report(s);
    for (;;) {
        mbs.realize(s);
        printf("%5g qerr=%10.4g uerr=%10.4g hNext=%g\n", s.getTime(), 
            myRNA.calcQConstraintNorm(s), myRNA.calcUConstraintNorm(s),
            myStudy.getPredictedNextStep());
        printf("      E=%14.8g (pe=%10.4g ke=%10.4g)\n",
            mbs.getEnergy(s), mbs.getPotentialEnergy(s), mbs.getKineticEnergy(s));

        cout << "QERR=" << s.getQErr() << endl;
        cout << "UERR=" << s.getUErr() << endl;

        //if (s.getTime() - std::floor(s.getTime()) < 0.2)
        //    display.addEphemeralDecoration( DecorativeSphere(10).setColor(Green) );

        display.report(s);
        saveEm.push_back(s);

       // if (myStudy.getT() >= 10*expectedPeriod)
         //   break;

        if (s.getTime() >= 10)
            break;

        // TODO: should check for errors or have or teach RKM to throw. 
        myStudy.step(s.getTime() + dt);
    }
    while(true) {
        for (int i=0; i < (int)saveEm.size(); ++i) {
            display.report(saveEm[i]);
            //display.report(saveEm[i]); // half speed
        }
        getchar();
    }
  } 
catch (const exception& e)
  {
    printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);
  }
}

