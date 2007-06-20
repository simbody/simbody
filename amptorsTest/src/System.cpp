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
 *
 * Implementation of System and SystemRep.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/System.h"
#include "simbody/internal/Subsystem.h"
#include "SystemRep.h"

#include <cassert>

namespace SimTK {

    ////////////
    // System //
    ////////////

bool System::isEmptyHandle() const {return rep==0;}
bool System::isOwnerHandle() const {return rep==0 || rep->myHandle==this;}

System::~System() {
    if (isOwnerHandle()) delete rep; 
    rep=0;
}

System::System(const System& src) : rep(0) {
    if (src.rep) {
        rep = src.rep->clone();
        rep->setMyHandle(*this);
    }
}

System& System::operator=(const System& src) {
    if (&src != this) {
        if (isOwnerHandle()) delete rep; 
        rep=0;
        if (src.rep) {
            rep = src.rep->clone();
            rep->setMyHandle(*this);
        }
    }
    return *this;
}


const String& System::getName()    const {return getRep().getName();}
const String& System::getVersion() const {return getRep().getVersion();}

int System::getNSubsystems() const {return getRep().getNSubsystems();}
const Subsystem& System::getSubsystem(int i) const {return getRep().getSubsystem(i);}
Subsystem& System::updSubsystem(int i) {return updRep().updSubsystem(i);}

void System::realize(const State& s, Stage g) const {
    getRep().realize(s,g);
}

void System::calcDecorativeGeometryAndAppend(const State& s, Stage g, Array<DecorativeGeometry>& dg) const {
    getRep().calcDecorativeGeometryAndAppend(s,g,dg);
}

Real System::calcTimescale(const State& s) const {
    return getRep().calcTimescale(s);
}

void System::calcYUnitWeights(const State& s, Vector& weights) const {
    return getRep().calcYUnitWeights(s,weights);
}

void System::calcYErrUnitTolerances(const State& s, Vector& tolerances) const {
    return getRep().calcYErrUnitTolerances(s,tolerances);
}

Real System::calcYErrorNorm(const State& s, const Vector& y_err) const {
    return getRep().calcYErrorNorm(s,y_err);
}

int System::takeOverSubsystem(Subsystem& src) {
    return updRep().takeOverSubsystem(src);
}

    ///////////////
    // SystemRep //
    ///////////////

void SystemRep::realize(const State& s, Stage g) const {
    Stage stageNow;
    while ((stageNow=s.getSystemStage()) < g) {
        switch (stageNow) {
        case Stage::Empty: {
            // Teach the State about the system & its subsystems.
            State& mutableState = const_cast<State&>(s);
            mutableState.setNSubsystems(getNSubsystems());
            realizeTopology(mutableState); 
            break;
        }
        case Stage::Topology:     realizeModel (const_cast<State&>(s)); break;
        case Stage::Model:        realizeInstance(s);     break;
        case Stage::Instance:     realizeTime(s);         break;
        case Stage::Time:         realizePosition(s);     break;
        case Stage::Position:     realizeVelocity(s);     break;
        case Stage::Velocity:     realizeDynamics(s);     break;
        case Stage::Dynamics:     realizeAcceleration(s); break;
        case Stage::Acceleration: realizeReport(s);       break;
        default: assert(!"System::realize(): bad stage");
        }
        // In case the concrete system didn't do anything with the
        // System Subsystem (subsystem 0), we'll just bump its stage
        // here.
        if (s.getSubsystemStage(0) == stageNow)
            s.advanceSubsystemToStage(0, stageNow.next());
        s.advanceSystemToStage(stageNow.next());
    }
}

void SystemRep::calcDecorativeGeometryAndAppend(const State& s, Stage stage, Array<DecorativeGeometry>& geom) const {
    assert(stage==Stage::Topology || s.getSystemStage() >= stage);
    for (int i=0; i<getNSubsystems(); ++i)
        subsystems[i].calcDecorativeGeometryAndAppend(s, stage, geom);
}


} // namespace SimTK

