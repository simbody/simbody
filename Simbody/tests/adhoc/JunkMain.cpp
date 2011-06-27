#include "SimTKsimbody.h"
#include "SimTKcommon/Testing.h"

#include <cstdio>
#include <iostream>
#include <exception>

using namespace SimTK;

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

bool firstTime;
State finalAnswer[20];

class EnergyReporter : public PeriodicEventReporter {
public:
    EnergyReporter(const MultibodySystem& system, Real reportInterval)
        : PeriodicEventReporter(reportInterval), system(system) {}
    void handleEvent(const State& state) const {
        system.realize(state, Stage::Velocity);
        printf("%8g\t%10g", state.getTime(),system.calcEnergy(state));
        const int tt = (int)state.getTime();
        if ((Real)tt != state.getTime()) {
            std::cout << std::endl;
            return;
        }
        if (firstTime) {
            finalAnswer[tt] = state;
            std::cout << std::endl;
        } else {
            Vector err = state.getQ()-finalAnswer[tt].getQ();
            printf("\t%10g", std::sqrt((~err*err)/err.size()));
            err = state.getU()-finalAnswer[tt].getU();
            printf("\t%10g\n", std::sqrt((~err*err)/err.size()));
        }
    }

private:
    const MultibodySystem& system;
};

int main() {
  try
  { // Create the system.
    
    MultibodySystem         system; system.setUseUniformBackground(true);
    SimbodyMatterSubsystem  matter(system);
    GeneralForceSubsystem   forces(system);
    Force::Gravity          gravity(forces, matter, -YAxis, 9.8);
    Force::GlobalDamper     damper(forces, matter, 10.0);

    const Real length = 3, mass = 10, rad = 0.5;
    Body::Rigid pendulumBody(MassProperties(mass,Vec3(0), mass*UnitInertia::sphere(rad)));
    pendulumBody.addDecoration(Transform(), DecorativeSphere(rad).setColor(Cyan));
    // Reverse so we slide along parent's x axis.
    MobilizedBody::BendStretch pendulum(matter.Ground(), Transform(Vec3(0)), 
                                        pendulumBody,    Transform(Vec3(0, length, 0))
                                        ,MobilizedBody::Reverse
                                        );

    MobilizedBody parent = pendulum;
    Vec3 offs(0,length,0), pos(0);
    int sign = -1;
    for (int i=1; i < 220; ++i) {
        Rotation R_BM = Test::randRotation();
        MobilizedBody::Pin mobod(parent, pos, 
                                 pendulumBody, Transform(R_BM,Vec3(0,length,0)));
        if ((i+1)%20 == 0) {
            parent = pendulum;
            if (sign == -1) {
                offs += Vec3(.5,0,0);
                sign = 1;
            } else sign = -1;
            pos = Vec3(sign*offs[0],length,0);
            pendulumBody.setDefaultRigidBodyMassProperties
                (MassProperties(mass+i/10, Vec3(0), (mass+i/10)*UnitInertia::sphere(rad)));
        } else {
            parent = mobod;
            pos = 0;
        }
    }
    
    Visualizer viz(system);

    // Add a menu for user interaction in the visualizer gui.
    Visualizer::InputSilo* silo = new Visualizer::InputSilo();
    viz.addInputListener(silo);
    Array_<std::pair<String,int> > runMenuItems;
    runMenuItems.push_back(std::make_pair("Go", 1));
    runMenuItems.push_back(std::make_pair("Quit", 2));
    viz.addMenu("Run", 1, runMenuItems);

    system.addEventReporter(new Visualizer::Reporter(viz, 1./10));
    system.addEventReporter(new EnergyReporter(system, 0.5));

    MySinusoid* func = new MySinusoid(4, .5, 0);// amp, freq, phase
    Constraint::PrescribedMotion motion(matter, 
        func, 
        pendulum, MobilizerQIndex(1));
    //Force::MobilityLinearSpring(forces,pendulum,1,100,0);

    RungeKuttaMersonIntegrator integ(system);
    //VerletIntegrator integ(system);
    //RungeKuttaFeldbergIntegrator integ(system);
    //RungeKutta3Integrator integ(system);
    //CPodesIntegrator integ(system);
        integ.setAccuracy(0.0001);
    TimeStepper ts(system, integ);

    system.realizeTopology();
    State startState = system.getDefaultState();
    pendulum.setOneQ(startState, 1, -0*0.5);   // translate along x
    pendulum.setOneQ(startState, 0, -0*Pi/4);  // rotate about z
    viz.report(startState);
    system.realize(startState, Stage::Position);
    LocalEnergyMinimizer::minimizeEnergy(system,startState,10);
    system.realize(startState, Stage::Position);
    viz.report(startState);
    for (int i=0; i < startState.getNU(); ++i)
        startState.updU()[i] = 0.2*Test::randReal();

    Real firstTimeAcc = 1e-9;
    firstTime=true;
    for (Real acc=1; acc >= .9e-7; acc /= 10) {
        Real freq = .5;
        //func->setFrequency(freq);
        std::cout << "Frequency is now " << freq << std::endl;

        //system.realizeTopology();
        State state = startState;
        system.realize(state, Stage::Velocity);

        //int menuId, item;
        //silo->waitForMenuPick(menuId, item);
        //if (item==2) break; // quit

        integ.resetAllStatistics();
        integ.setAccuracy(firstTime ? firstTimeAcc : acc);
        ts.initialize(state);
        Real start = threadCpuTime(), rt = realTime();
        ts.stepTo(20.0);
        printf("Took %gs CPU %gs REAL @acc=%g with \n%s\n",
            threadCpuTime()-start, realTime()-rt, integ.getAccuracyInUse(),
            integ.getMethodName());
        printf("nsteps=%d/%d avg=%gms\n", integ.getNumStepsTaken(),
            integ.getNumStepsAttempted(), 20000./integ.getNumStepsTaken());
        viz.report(ts.getState());
        /*
        if (firstTime) {
            finalAnswer=ts.getState();
            firstTime=false;
        } else {
            Vector err = ts.getState().getQ()-finalAnswer.getQ();
            printf("qerr=%g\n", std::sqrt((~err*err)/err.size()));
            err = ts.getState().getU()-finalAnswer.getU();
            printf("uerr=%g\n", std::sqrt((~err*err)/err.size()));
        }*/
        firstTime = false;
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