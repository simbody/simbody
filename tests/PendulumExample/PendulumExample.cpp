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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 * The simple 2d pendulum example from the user's manual.
 */

#include "Simbody.h"
#include "simbody/internal/NumericalMethods.h"
#include "simbody/internal/DecorativeGeometry.h"
#include "simbody/internal/VTKReporter.h"

#include "windows.h"

#include <cmath>
#include <cstdio>
#include <exception>
#include <vector>

using namespace std;
using namespace SimTK;

static const Real Pi = std::acos(-1.), RadiansPerDegree = Pi/180;
static const int  GroundBodyNum = 0; // ground is always body 0

static const Real m = 5;   // kg
static const Real g = 9.8; // meters/s^2; apply in –y direction
static const Real d = 1.5; // meters
static const Real initialTheta   = 30;             // degrees
static const Real expectedPeriod = 2*Pi*sqrt(d/g); // s

class MySimbodyPendulum : public SimbodyMatterSubsystem {
public:
    MySimbodyPendulum() 
    {
        pendBodyNum =
            addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,-d/2,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                GroundBodyNum,           // parent body
                Transform(Vec3(1,1,1)),             // jt frame on parent             
                Mobilizer(Mobilizer::Ball, false)); // joint type; pin always aligns z axes

        int pendBodyNum2 =
            addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,-d/2,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                pendBodyNum,           // parent body
                Transform(Vec3(0,-d/2,0)),             // jt frame on parent (bottom)             
                Mobilizer(Mobilizer::Ball, false)); // joint type; pin always aligns z axes
        int pendBodyNum2a =
            addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,-d/2,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                pendBodyNum2,           // parent body
                Transform(Vec3(0,-d/2,0)),             // jt frame on parent (bottom)             
                Mobilizer(Mobilizer::Ball, false)); // joint type; pin always aligns z axes
        int pendBodyNum2b =
            addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,-d/2,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                pendBodyNum2a,           // parent body
                Transform(Vec3(0,-d/2,0)),             // jt frame on parent (bottom)             
                Mobilizer(Mobilizer::Ball, false)); // joint type; pin always aligns z axes
        int pendBodyNum2c =
            addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,-d/2,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                pendBodyNum2b,           // parent body
                Transform(Vec3(0,-d/2,0)),             // jt frame on parent (bottom)             
                Mobilizer(Mobilizer::Ball, false)); // joint type; pin always aligns z axes

        int pendBodyNum3 =
            addRigidBody(
                MassProperties(m,        // body mass, center of mass, inertia
                               Vec3(0,-d/2,0), 
                               InertiaMat(Vec3(0,-d/2,0), m)+InertiaMat(1e-3,1e-3,1e-3)),
                Transform(Vec3(0,d/2,0)),// jt frame on body (aligned w/body frame)
                pendBodyNum2c,           // parent body
                Transform(Vec3(0,-d/2,0)),             // jt frame on parent (bottom)             
                Mobilizer(Mobilizer::Ball, false)); // joint type; pin always aligns z axes
       
    int theConstraint =
      addCoincidentStationsConstraint(1, Vec3(0.3,-d/2,0),
                                           pendBodyNum3, Vec3(0,-d/2,0));

    int theConstraint2 =
       addConstantDistanceConstraint(0, Vec3(2,-3,0),
                                          3, Vec3(0,-d/2,0),1.5);

    }

    Real getPendulumAngle(const State& s) const {
        const Vec4 aa = getMobilizerConfiguration(s,pendBodyNum).R().convertToAngleAxis();
        return aa[0]/RadiansPerDegree;
    }

    // Assume rotation around z
    void setPendulumAngle(State& s, Real angleInDegrees) {
        const Vec4 aa(angleInDegrees*RadiansPerDegree,0, 0, 1);
        Quaternion q; q.setToAngleAxis(aa);
        setMobilizerConfiguration(s,pendBodyNum,Transform(RotationMat(q)));
    }
private:
    int pendBodyNum;
};

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
        PerBodyInfo(int b, bool d) : bnum(b), isDuplex(d) { }
        int  bnum;
        bool isDuplex;
    };
    std::vector<PerBodyInfo> bodyInfo;
    int end1, end2;
public:
    MyRNAExample(int nsegs, bool shouldFlop) 
    {
        bodyInfo.push_back(PerBodyInfo(0, false)); // placeholder for ground
        end1 = makeChain(GroundBodyNum, Vec3(0), nsegs, shouldFlop);
        end2 = makeChain(GroundBodyNum, Vec3(20,0,0), nsegs, shouldFlop);

        if (true) {
            int theConstraint2 =
               addConstantDistanceConstraint(end1, Vec3(0, -HalfHeight,0),
                                             end2, Vec3(0, -HalfHeight,0), 10);
        }

    }

    void decorateBody(int bodyNum, VTKReporter& display) const {
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

    int makeChain(int startBody, const Vec3& startOrigin, int nSegs, bool shouldFlop) {
        int baseBody = startBody;
        Vec3 origin = startOrigin;
        int lastDup = -1;
        for (int seg=0; seg < nSegs; ++seg) {
            int left1 = addRigidBody(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity),
                             Transform(Vec3(0, ConnectorHalfHeight, 0)),
                             baseBody,
                             Transform(origin + Vec3(-DuplexRadius,-HalfHeight,0)),
                             Mobilizer(Mobilizer::Ball, false));
            bodyInfo.push_back(PerBodyInfo(left1, false));

            int left2 = addRigidBody(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity),
                             Transform(Vec3(0, ConnectorHalfHeight, 0)),
                             left1,
                             Transform(Vec3(0, -ConnectorHalfHeight, 0)),
                             Mobilizer(Mobilizer::Ball, false));
            bodyInfo.push_back(PerBodyInfo(left2, false));

            int rt1 = addRigidBody(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity),
                             Transform(Vec3(0, ConnectorHalfHeight, 0)),
                             baseBody,
                             Transform(origin + Vec3(DuplexRadius,-HalfHeight,0)),
                             Mobilizer(Mobilizer::Ball, false));
            bodyInfo.push_back(PerBodyInfo(rt1, false));

            int rt2 = addRigidBody(calcConnectorMassProps(ConnectorRadius, ConnectorHalfHeight, ConnectorDensity),
                             Transform(Vec3(0, ConnectorHalfHeight, 0)),
                             rt1,
                             Transform(Vec3(0, -ConnectorHalfHeight, 0)),
                             Mobilizer(Mobilizer::Ball, false));
            bodyInfo.push_back(PerBodyInfo(rt2, false));

            int dup = addRigidBody(calcDuplexMassProps(DuplexRadius, HalfHeight, NAtoms, AtomMass),
                                Transform(Vec3(-DuplexRadius, HalfHeight, 0)),
                                rt2,
                                Transform(Vec3(0, -ConnectorHalfHeight, 0)),
                                Mobilizer(Mobilizer::Ball, false));
            bodyInfo.push_back(PerBodyInfo(dup, true));

            if (!shouldFlop) {
                int theConstraint =
                    addCoincidentStationsConstraint(left2, Vec3(0, -ConnectorHalfHeight, 0),
                                                    dup, Vec3(DuplexRadius, HalfHeight, 0));
                //int theConstraint =
                  //  addConstantDistanceConstraint(rt2, Vec3(0, -ConnectorHalfHeight, 0),
                  //                                dup, Vec3(DuplexRadius, HalfHeight, 0), .01);
            }

            baseBody = dup;
            origin = Vec3(0);
            lastDup = dup;
        }
        return lastDup;
    }

    MassProperties calcDuplexMassProps(
        Real halfHeight, Real r, int nAtoms, Real atomMass)
    {
        const Real pitch = 2*Pi/halfHeight;
        const Real trans = (2*halfHeight)/(nAtoms-1);
        const Real rot = pitch*trans;
        InertiaMat iner(0);
        Vec3 com(0);
        Real mass = 0;
        for (int i=0; i<nAtoms; ++i) {
            const Real h = halfHeight - i*trans;
            const Real th = i*rot;
            const Vec3 p1(-r*cos(th),h,r*sin(th)), p2(r*cos(th),h,-r*sin(th));
            mass += 2*atomMass;
            iner += InertiaMat(p1, atomMass) + InertiaMat(p2, atomMass);
            com += atomMass*p1 + atomMass*p2;
        }
        return MassProperties(mass,com/mass,iner);
    }

    MassProperties calcConnectorMassProps(Real r, Real halfHeight, Real density)
    {
        const Real volume = Pi*r*r*halfHeight;
        const Real mass = volume*density;
        const Vec3 com = Vec3(0);
        const InertiaMat iner = mass*InertiaMat::cylinderAlongY(r, halfHeight);

        return MassProperties(mass,com,iner);
    }

    void addDuplexDecorations(int bodyNum, Real r, Real halfHeight, Real slop, int nAtoms,
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

    void addConnectorDecorations(int bodyNum, Real r, Real halfHeight, Real endSlop,  
                                 VTKReporter& display) const
    {
        display.addDecoration(bodyNum, Transform(), 
            DecorativeCylinder(r, halfHeight-endSlop).setColor(Blue));
    }
};


int main(int argc, char** argv) {
    std::vector<State> saveEm;

    try { // If anything goes wrong, an exception will be thrown.
        int nseg = NSegments;
        int shouldFlop = 0;
        if (argc > 1) sscanf(argv[1], "%d", &nseg);
        if (argc > 2) sscanf(argv[2], "%d", &shouldFlop);
        //printf("Pendulum starting at angle +%g degrees from vertical.\n", start);

        // Create a multibody system using Simbody.
        MyRNAExample myRNA(nseg, shouldFlop != 0);
        const Vec3 attachPt(150, -40, -50);
        GeneralForceElements forces;

        MultibodySystem mbs;
        mbs.setMatterSubsystem(myRNA);

        forces.addLinearTwoPointSpring(0, attachPt,
                                       myRNA.getNBodies()-1, Vec3(0),
                                       1000.,  // stiffness
                                       1.);    // natural length

       /* forces.addLinearTwoPointSpring(0, -attachPt,
                                       myRNA.getNBodies()-1, Vec3(0),
                                       1000.,  // stiffness
                                       1.);    // natural length
        */

        forces.addGlobalMobilityDamping(1000);

        mbs.addForceSubsystem(forces);
        UniformGravitySubsystem ugs(Vec3(0, -g, 0));
        mbs.addForceSubsystem(ugs);

        State s;
        mbs.realize(s, Stage::Built);
        //myRNA.setUseEulerAngles(s,true);
        mbs.realize(s, Stage::Modeled);

        ugs.updGravity(s) *= 10;
        ugs.disableGravity(s);
        ugs.enableGravity(s);
        ugs.updZeroHeight(s) = -0.8;
        //cout << "STATE AS MODELED: " << s;
       
        //myPend.setPendulumAngle(s, start);

        // And a study using the Runge Kutta Merson integrator
        bool suppressProject = false;
        RungeKuttaMerson myStudy(mbs, s, suppressProject);
        myStudy.setAccuracy(1e-2);
        myStudy.setConstraintTolerance(1e-3);
        myStudy.setProjectEveryStep(false);

        VTKReporter display(mbs);
        for (int i=1; i<myRNA.getNBodies(); ++i)
            myRNA.decorateBody(i, display);
        myRNA.decorateGlobal(display);

        DecorativeLine rbProto; rbProto.setColor(Orange).setLineThickness(3);
        display.addRubberBandLine(0, attachPt,myRNA.getNBodies()-1,Vec3(0), rbProto);
        //display.addRubberBandLine(0, -attachPt,myRNA.getNBodies()-1,Vec3(0), rbProto);

        const Real dt = 0.05; // output intervals

        printf("time  nextStepSize\n");

        s.updTime() = 0;
        for (int i=0; i<100; ++i)
            saveEm.push_back(s);    // delay
        display.report(s);

        myStudy.initialize();
        saveEm.push_back(s);
        for (int i=0; i<100; ++i)
            saveEm.push_back(s);    // delay
        display.report(s);
        for (;;) {
            printf("%5g qerr=%10.4g uerr=%10.4g hNext=%g\n", s.getTime(), 
                myRNA.calcQConstraintNorm(s), myRNA.calcUConstraintNorm(s),
                myStudy.getPredictedNextStep());
            printf("      E=%14.8g (pe=%10.4g ke=%10.4g)\n",
                mbs.getEnergy(s), mbs.getPotentialEnergy(s), mbs.getKineticEnergy(s));

            display.report(s);
            saveEm.push_back(s);

           // if (myStudy.getT() >= 10*expectedPeriod)
             //   break;
    
            if (s.getTime() >= 20)
                break;

            // TODO: should check for errors or have or teach RKM to throw. 
            myStudy.step(s.getTime() + dt);
        }

        for (int i=0; i < (int)saveEm.size(); ++i)
            display.report(saveEm[i]);
    } 
    catch (const exception& e) {
        printf("EXCEPTION THROWN: %s\n", e.what());
        exit(1);
    }
}

