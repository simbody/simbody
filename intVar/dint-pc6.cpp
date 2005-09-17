
#include "dint-pc6.h"
#include "AtomTree.h"

#include "cdsMath.h"

#ifdef USE_CDS_NAMESPACE 
using namespace CDS;
#endif /* USE_CDS_NAMESPACE */

class TimeStepAdj {
    const IVM* ivm;
    double deltaE_ref;
    double tau;         // set tau <=0 to disable
    double energyPrev;
    double maxFactor;
    bool   active_;
public:
    TimeStepAdj(const IVM* ivm)
      : ivm(ivm), deltaE_ref(ivm->Etolerance()), tau(ivm->responseTime()), 
        energyPrev(ivm->Etotal()), maxFactor(ivm->maxTSFactor()), 
        active_(ivm->adjustTS())
    { }

    double getFactor(const double& energy) {
        double ret = 1.0;
        if ( !active() ) return ret;
        if (ivm->verbose()&InternalDynamics::printStepDebug)
            cout << "TimeStepAdj: energy: " << energy
                 << " pe: " << ivm->Epotential()
                 << " ke: " << ivm->Ekinetic()
                 << " energyPrev: " << energyPrev << '\n';
        double deltaE = fabs( energy - energyPrev );
        if ( deltaE > ivm->maxDeltaE() ) 
            throw InternalDynamics::Exception("large time-step");
        if ( deltaE > 0.0 ) 
            ret = sqrt(1 + (deltaE_ref-deltaE) / (tau*deltaE));
        ret = CDSMath::min(ret , maxFactor);
        if (ivm->verbose()&InternalDynamics::printStepDebug)
            cout << "TimeStepAdj: deltaE: " << deltaE 
                 << "   factor: " << ret << '\n';
        return ret;
    }
    void setEnergyPrev(const double& energy) { energyPrev = energy; }
    bool active() { return active_ && tau>0.0; }
};

class VelScale {
    // set tau <=0 to disable
    double tau;
    bool   active_;
public:
    VelScale(const double& tau, bool active)
      : tau(tau), active_(active)
    { }

    double getScale(const double& temp, const double& bathTemp) {
        double ret = 1.0;
        if (temp > 0.0)
            ret = sqrt( 1 + (bathTemp-temp) / (tau*temp) );
        return ret;
    }
    bool active() { return active_ && tau>0.0; }
};

void PC6::init(const RVec& pos,
               const RVec& vel,
               const RVec& acc)
{
    this->pos = pos;
    this->vel = vel;
    this->acc = acc;

    pos0 = pos; vel0 = vel;

    dq3.resize(pos.size());
    dq4.resize(pos.size());
    dq5.resize(pos.size());    //lazy initialization

    timeStepAdj = new TimeStepAdj(ivm);
    velScale    = new VelScale(ivm->responseTime(),
                ivm->scaleVel());
    dq3.set(0.0);
    dq4.set(0.0);
    dq5.set(0.0);
}

PC6::~PC6()
{ 
    delete timeStepAdj;
    delete velScale;
}

//
// reset all coordinates/velocities to values at beginning of step
//
void
PC6::stepUndo() {
    pos = pos0; vel = vel0; acc = acc0; dq3 = dq30; dq4 = dq40; dq5 = dq50;
}

//
// 6th order predictor corrector - internal coordinate version
//
void
PC6::step(double& timeStep)
{
    using namespace InternalDynamics;

    // save
    acc0 = acc; dq30 = dq3; dq40 = dq4; dq50 = dq5;

    //prediction
    vel *= ( timeStep );          //rescale vel/acc
    acc *= ( 0.5*CDSMath::sq(timeStep) );
    dq3 *= ( CDSMath::ipow(timeStep,3)/(2*3) );
    dq4 *= ( CDSMath::ipow(timeStep,4)/(2*3*4) );
    dq5 *= ( CDSMath::ipow(timeStep,5)/(2*3*4*5) );

    pos += vel + acc +     dq3 +     dq4 +      dq5;
    vel +=   2.0*acc + 3.0*dq3 + 4.0*dq4 +  5.0*dq5;
    acc +=             3.0*dq3 + 6.0*dq4 + 10.0*dq5;
    dq3 +=                       4.0*dq4 + 10.0*dq5;
    dq4 +=                                  5.0*dq5;

    prevAcc = acc;

    vel *= ( 1.0/timeStep );   //get unscaled velocity
    tree()->enforceConstraints(pos,vel);

    tree()->setPosVel(pos,vel);

    ivm->calcEnergy();  //calc energies, derivatives
    cout << "PE=" << std::setprecision(16) << ivm->Epotential() << endl;
    cout << "KE=" << ivm->Ekinetic() << endl;
    if ( ivm->verbose()&printStepDebug )
    cout << "PC6::step: before velocity scale: energy: "
         << ivm->Etotal() << " temp: " << ivm->currentTemp() << '\n';
    double newTimeStep = timeStep;
    try {
        newTimeStep = timeStep * timeStepAdj->getFactor(ivm->Etotal());
    }
    catch ( InternalDynamics::Exception e ) {
        if ( String(e.mess).contains("large time-step") )
            stepUndo();
        throw;
    }
    pos0 = pos;
    acc = tree()->getAccel();

    vel *= ( timeStep );          //rescale vel/acc
    acc *= ( 0.5*CDSMath::sq(timeStep) );

    //correction
    RVec dR = acc - prevAcc;

    pos +=   3.0/16.0  * dR;   //the value 3/20 should be used if dim=3*natom
    vel += 251.0/360.0 * dR;
    dq3 +=  11.0/18.0  * dR;
    dq4 +=   1.0/6.0   * dR;
    dq5 +=   1.0/60.0  * dR;

    vel *= ( 1.0/timeStep );      //undo vel/acc scaling
    acc *= ( 2.0/CDSMath::sq(timeStep) );
    dq3 *= ( (2*3) / CDSMath::ipow(timeStep,3) );
    dq4 *= ( (2*3*4) / CDSMath::ipow(timeStep,4) );
    dq5 *= ( (2*3*4*5) / CDSMath::ipow(timeStep,5) );

    if ( velScale->active() ) {  //scale velocity - for constant temp. run
        tree()->setVel(vel);
        ivm->calcTemperature();
        float_type scaleFac = velScale->getScale(ivm->currentTemp(),
                                                 ivm->bathTemp());
        vel *= scaleFac; acc *= scaleFac;
        dq3 *= scaleFac; dq4 *= scaleFac; dq5 *= scaleFac;
        tree()->setVel(vel);
    }

    //resetCM needs to done here to avoid interfering with the auto timestepper
    // note that resetCM leaves acc, dq3, dq4, dq5 slightly inconsistent.
    ivm->resetCM(); 
    //recalc energy
    vel0 = vel; 
    ivm->calcTemperature(); 
    ivm->setEtotal( ivm->Epotential() + ivm->Ekinetic() );  //FIX: need this?
    timeStepAdj->setEnergyPrev(ivm->Etotal());
    if ( ivm->verbose()&printStepDebug )
        cout << "PC6::step: after velocity scale: energy: "
             << ivm->Etotal() << " temp: " << ivm->currentTemp() << '\n';

    timeStep = newTimeStep;
}
