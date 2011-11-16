/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009-11 Stanford University and the Authors.        *
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
 * Adhoc main program for playing with prescribed motion.
 */

#include "SimTKsimbody.h"
#include "SimTKcommon/Testing.h"

#include <cstdio>
#include <exception>
#include <iostream>
using std::cout; using std::endl;
using std::clog; using std::cerr;

using namespace SimTK;

class MyReporter : public PeriodicEventReporter {
public:
    MyReporter(const MultibodySystem& system, 
        const Measure& power,
        const Measure& work,
        Real reportInterval)
    :   PeriodicEventReporter(reportInterval), 
        m_system(system), m_power(power), m_work(work)
    {}

    ~MyReporter() {}

    void handleEvent(const State& state) const {
        cout << state.getTime() << " " << m_system.calcEnergy(state);
        cout << " " << m_power.getValue(state) << " " << m_work.getValue(state) 
            << " " << 
            m_system.calcEnergy(state) - m_work.getValue(state);
        //for (int i=0; i < state.getNQ(); ++i)
        //    cout << " " << state.getQ()[i];
        cout << "\n";
        const SimbodyMatterSubsystem& matter = m_system.getMatterSubsystem();
        cout << "qperr=" << matter.calcMotionErrors(state,Stage::Position) << "\n";
        cout << "uperr=" << matter.calcMotionErrors(state,Stage::Velocity) << "\n";
        cout << "qerr=" << state.getQErr() << "\n";
        cout << "uerr=" << state.getUErr() << "\n";
    }
private:
    const MultibodySystem& m_system;
    const Measure          m_power, m_work;
};

class MySinusoid : public Function {
public:
    // Prescribe q=amplitude*sin(frequency*t + phase)
    MySinusoid(Real amplitude, Real frequency, Real phase=0) 
    :   a(amplitude), w(frequency), p(phase) {}

    void setFrequency(Real frequency) {w=frequency;}

    // These are the virtual methods you must write.

    virtual Real calcValue(const Vector& x) const {
        const Real t = x[0]; // we expect just one argument
        return a*std::sin(w*t + p);
    }

    virtual Real calcDerivative(const Array_<int>& derivComponents,
                                const Vector&      x) const {
        const Real t = x[0]; // time is the only argument
        const int  order = derivComponents.size();
        assert(1 <= order && order <= 2); // only 1st & 2nd derivs implemented
        if (order == 1) 
            return a*w*std::cos(w*t + p);
        else // order == 2
            return -a*w*w*std::sin(w*t + p);
    }

    virtual int getArgumentSize() const {return 1;} // just time
    virtual int getMaxDerivativeOrder() const {return 2;} // need two time derivs

private:
    Real a, w, p;
};

/* This Measure returns the instantaneous power being generated by the
constraints and motions. */
template <class T>
class PowerMeasure : public Measure_<T> {
public:
    SimTK_MEASURE_HANDLE_PREAMBLE(PowerMeasure, Measure_<T>);

    PowerMeasure(Subsystem& sub, const SimbodyMatterSubsystem& matter)
    :   Measure_<T>(sub, new Implementation(matter), AbstractMeasure::SetHandle()) {}
    SimTK_MEASURE_HANDLE_POSTSCRIPT(PowerMeasure, Measure_<T>);
};


template <class T>
class PowerMeasure<T>::Implementation : public Measure_<T>::Implementation {
public:
    Implementation(const SimbodyMatterSubsystem& matter) 
    :   Measure_<T>::Implementation(1), m_matter(matter) {}

    // Default copy constructor, destructor, copy assignment are fine.

    // Implementations of virtual methods.
    Implementation* cloneVirtual() const {return new Implementation(*this);}
    int getNumTimeDerivativesVirtual() const {return 0;}
    Stage getDependsOnStageVirtual(int order) const 
    {   return Stage::Acceleration; }

    void calcCachedValueVirtual(const State& s, int derivOrder, T& value) const
    {
        SimTK_ASSERT1_ALWAYS(derivOrder==0,
            "PowerMeasure::Implementation::calcCachedValueVirtual():"
            " derivOrder %d seen but only 0 allowed.", derivOrder);

        value = m_matter.calcMotionPower(s) + m_matter.calcConstraintPower(s);
    }
private:
    const SimbodyMatterSubsystem& m_matter;
};

const Real RodLength = 1.1;
int main() {
  try
  { // Create the system.
    
    MultibodySystem         system;
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::UniformGravity   gravity(forces, matter, Vec3(0, -9.8, 0));

    PowerMeasure<Real> powMeas(matter, matter);
    Measure::Zero zeroMeas(matter);
    Measure::Integrate workMeas(matter, powMeas, zeroMeas); 

    Body::Rigid pendulumBody(MassProperties(1.0, Vec3(0), Inertia(1)));
    pendulumBody.addDecoration(Transform(), DecorativeSphere(0.1));


    // Prescribed system.
    MobilizedBody::Pin pendulum(matter.Ground(), Transform(Vec3(0)), 
                                pendulumBody,    Transform(Vec3(0, 1, 0)));
    Motion::Sinusoid(pendulum, Motion::Position, Pi/8, 2*Pi, Pi/4); // amp, rate, phase
    MobilizedBody::Pin pendulum2(pendulum, Transform(Vec3(0)), 
                                 pendulumBody,    Transform(Vec3(0, 1, 0)));
    MobilizedBody::Pin pendulum3(pendulum2, Vec3(0),
                                 pendulumBody, Vec3(0,.5,0));
    Motion::Steady(pendulum3, 4*Pi); // rate
    Force::MobilityLinearSpring(forces, pendulum2, MobilizerUIndex(0),
        100, 0*(Pi/180));

    // Identical constrained system.
    MobilizedBody::Pin cpendulum(matter.Ground(), Transform(Vec3(2,0,0)), 
                                 pendulumBody,    Transform(Vec3(0, 1, 0)));
    Constraint::PrescribedMotion(matter, 
                                 new MySinusoid(Pi/8,2*Pi,Pi/4), //amp,rate,phase
                                 cpendulum, MobilizerQIndex(0));
    MobilizedBody::Pin cpendulum2(cpendulum, Transform(Vec3(0)), 
                                  pendulumBody, Transform(Vec3(0, 1, 0)));
    MobilizedBody::Pin cpendulum3(cpendulum2, Vec3(0),
                                 pendulumBody, Vec3(0,.5,0));
    Constraint::ConstantSpeed(cpendulum3, 4*Pi);
    Force::MobilityLinearSpring(forces, cpendulum2, MobilizerUIndex(0),
        100, 0*(Pi/180));

    Constraint::Rod(pendulum3, cpendulum3, RodLength);

    // Identical mixed system 1.
    MobilizedBody::Pin m1pendulum(matter.Ground(), Transform(Vec3(4,0,0)), 
                                  pendulumBody,    Transform(Vec3(0, 1, 0)));
    Motion::Sinusoid(m1pendulum, Motion::Position, Pi/8, 2*Pi, Pi/4); // amp, rate, phase
    MobilizedBody::Pin m1pendulum2(m1pendulum, Transform(Vec3(0)), 
                                   pendulumBody, Transform(Vec3(0, 1, 0)));
    MobilizedBody::Pin m1pendulum3(m1pendulum2, Vec3(0),
                                   pendulumBody, Vec3(0,.5,0));
    Constraint::ConstantSpeed(m1pendulum3, 4*Pi);
    Force::MobilityLinearSpring(forces, m1pendulum2, MobilizerUIndex(0),
        100, 0*(Pi/180));


    // Identical mixed system 2.
    MobilizedBody::Pin m2pendulum(matter.Ground(), Transform(Vec3(6,0,0)), 
                                 pendulumBody,    Transform(Vec3(0, 1, 0)));
    Constraint::PrescribedMotion(matter, 
                                 new MySinusoid(Pi/8,2*Pi,Pi/4), //amp,rate,phase
                                 m2pendulum, MobilizerQIndex(0));
    MobilizedBody::Pin m2pendulum2(m2pendulum, Transform(Vec3(0)), 
                                   pendulumBody, Transform(Vec3(0, 1, 0)));
    MobilizedBody::Pin m2pendulum3(m2pendulum2, Vec3(0),
                                   pendulumBody, Vec3(0,.5,0));
    Motion::Steady(m2pendulum3, 4*Pi); // rat
    Force::MobilityLinearSpring(forces, m2pendulum2, MobilizerUIndex(0),
        100, 0*(Pi/180));

    Constraint::Rod(m1pendulum3, m2pendulum3, RodLength);


    Visualizer viz(system);
    system.addEventReporter(new Visualizer::Reporter(viz, 0.01));
    system.addEventReporter(new MyReporter(system, powMeas, workMeas, 0.01));

   
    // Initialize the system and state.
    
    system.realizeTopology();
    State state = system.getDefaultState();
    system.realize(state, Stage::Instance);
    const int nq = state.getNQ();
    const int nu = state.getNU();
    const int m  = state.getNMultipliers();

    viz.report(state);
    system.realize(state, Stage::Velocity);
    clog << "Default state -- hit ENTER\n";
    clog << "t=" << state.getTime() 
         << "\nq=" << state.getQ() 
         << "\nu=" << state.getU() 
         << "\nqerr=" << state.getQErr()
         << "\nuerr=" << state.getUErr()
         << endl;
    char c=getchar();

    state.setTime(0);
    system.realize(state, Stage::Time);


    clog << "After realize(Time) motion qerr=" 
         << matter.calcMotionErrors(state, Stage::Position) << "\n";
    system.prescribeQ(state);

    clog << "After prescribe() motion qerr=" 
         << matter.calcMotionErrors(state, Stage::Position) << "\n";

    //if (matter.prescribe(state, Stage::Position)) clog << "Some PresQ\n";
    //else clog << "NO PresQ\n";
    //clog << "after prescribe q=" << state.getQ() << "\n";
    system.realize(state, Stage::Position);


    ProjectResults projResults;

    system.projectQ(state, Vector(), ProjectOptions(1e-10), projResults);
    //Assembler asmb(system);
    //asmb.setAccuracy(1e-10).assemble(state);

    viz.report(state);
    clog << "After prescribe(Position) & project -- hit ENTER\n";
    clog << "t=" << state.getTime() 
         << "\nq=" << state.getQ() 
         << "\nu=" << state.getU() 
         << "\nqerr=" << state.getQErr()
         << endl;

    const int nfq = matter.getFreeQIndex(state).size();
    clog << "freeQ:     " << matter.getFreeQIndex(state) << "\n";

    Vector q = state.getQ(), packedQ(nfq), unpackedQ(nq);
    clog << "allQ         =" << q << "\n";
    matter.packFreeQ(state, q, packedQ);
    clog << "packedFreeQ  =" << packedQ << "\n";
    matter.unpackFreeQ(state, packedQ, unpackedQ);
    clog << "unpackedFreeQ=" << unpackedQ << "\n";

    c=getchar();

    clog << "After realize(Position) motion uerr=" 
         << matter.calcMotionErrors(state, Stage::Velocity) << "\n";
    system.prescribeU(state);
    clog << "After prescribe() motion uerr=" 
         << matter.calcMotionErrors(state, Stage::Velocity) << "\n";

    system.realize(state, Stage::Velocity);
    system.projectU(state, Vector(), ProjectOptions(1e-10), projResults);
    viz.report(state);
    clog << "After prescribe(Velocity) & project -- hit ENTER\n";
    clog << "t=" << state.getTime() 
         << "\nq=" << state.getQ() 
         << "\nu=" << state.getU() 
         << "\nuerr=" << state.getUErr()
         << endl;

    const int nfu = matter.getFreeUIndex(state).size();
    clog << "freeU:     " << matter.getFreeUIndex(state) << "\n";

    Vector u = state.getU(), packedU(nfu), unpackedU(nu);
    clog << "allU         =" << u << "\n";
    matter.packFreeU(state, u, packedU);
    clog << "packedFreeU  =" << packedU << "\n";
    matter.unpackFreeU(state, packedU, unpackedU);
    clog << "unpackedFreeU=" << unpackedU << "\n";


    c=getchar();

    system.realize(state, Stage::Acceleration);
    clog << "After realize(Acceleration) motion udoterr=" 
         << matter.calcMotionErrors(state, Stage::Acceleration) << "\n";

    clog << "After realize(Acceleration) -- hit ENTER\n";
    clog << "t=" << state.getTime() 
         << "\nq=" << state.getQ()
         << "\nu=" << state.getU() 
         << "\nudot=" << state.getUDot() 
         << "\ntau=" << pendulum.getTauAsVector(state) << pendulum2.getTauAsVector(state) << pendulum3.getTauAsVector(state)
         << endl;

    Vector_<SpatialVec> reactionForces, reactionForcesFreebody;
    matter.calcMobilizerReactionForces(state, reactionForces);
    matter.calcMobilizerReactionForcesUsingFreebodyMethod(state, reactionForcesFreebody);

    clog << "reactions PA+z: " << reactionForces << endl;
    clog << "react freebody: " << reactionForcesFreebody << endl;
    clog << "diff=" << reactionForces-reactionForcesFreebody << "\n";
    SimTK_TEST_EQ_SIZE(reactionForces, reactionForcesFreebody, nu);

    clog << "tau=" << matter.getMotionMultipliers(state) << "\n";
    Vector motFrcs;
    matter.findMotionForces(state, motFrcs);
    clog << "motion frc=" << motFrcs << "\n";
    clog << "gen speeds=" << matter.getU(state) << "\n";
    clog << "motion pwr=" << matter.calcMotionPower(state) << "\n";

    clog << "lambda=" << matter.getConstraintMultipliers(state) << "\n";
    Vector_<SpatialVec> consBodyFrc;
    Vector consMobFrc;
    matter.findConstraintForces(state, consBodyFrc, consMobFrc);
    clog << "cons bfrc=" << consBodyFrc << "\n";
    clog << "cons mfrc=" << consMobFrc << "\n";
    clog << "cons pwr=" << matter.calcConstraintPower(state) << "\n";


    clog << "freeUDot:  " << matter.getFreeUDotIndex(state) << "\n";
    clog << "knownUDot: " << matter.getKnownUDotIndex(state) << "\n";

    Matrix GMInvGt_r;
    matter.calcProjectedMInv(state, GMInvGt_r);

    Matrix G, M, Minv_r;
    matter.calcG(state, G);
    matter.calcM(state, M);
    matter.calcMInv(state, Minv_r);
    Matrix Minv = M.invert();

    const Array_<UIndex>& freeUDot = matter.getFreeUDotIndex(state);
    const int nfudot = freeUDot.size();
    Matrix Gr(m,nfudot), Mr(nfudot,nfudot);

    for (int i=0; i < nfu; ++i) {
        Gr(i) = G(freeUDot[i]);
        matter.packFreeU(state, M(freeUDot[i]), Mr(i));
    }

    cout << std::fixed;
    cout << "G=" << G;
    cout << "Gr=" << Gr;
    cout << "M=" <<  M;
    cout << "Mr=" <<  Mr;

    cout << "inv(M) =" << Minv;
    cout << "Minv_r =" << Minv_r;
    cout << "inv(Mr)=" << Mr.invert();

    clog << "GMInvGt_r     =" << GMInvGt_r;
    clog << "GMInv_rGt     =" << G*Minv_r*~G;
    clog << "Gr inv(Mr) ~Gr=" << Gr*Mr.invert()*~Gr;
    clog << "GMInvGt       =" << G*Minv*~G;


    c=getchar();

    //pendulum.setOneU(state, 0, 1.0);

    
    // Simulate it.


    //ExplicitEulerIntegrator integ(system);
    //RungeKutta3Integrator integ(system);
    //VerletIntegrator integ(system);
    RungeKuttaMersonIntegrator integ(system);
    //RungeKuttaFeldbergIntegrator integ(system);
    //integ.setAllowInterpolation(false);
    //integ.setProjectInterpolatedStates(false);
    //integ.setMinimumStepSize(1e-1);
    integ.setAccuracy(1e-2);
    //integ.setConstraintTolerance(1e-3);
    TimeStepper ts(system, integ);
    ts.initialize(state);
    ts.stepTo(10.0);

  } catch (const std::exception& e) {
    cerr << "EXCEPTION THROWN: " << e.what() << "\n";
    exit(1);

  } catch (...) {
    cerr << "UNKNOWN EXCEPTION THROWN\n";
    exit(1);
  }

    return 0;
}
