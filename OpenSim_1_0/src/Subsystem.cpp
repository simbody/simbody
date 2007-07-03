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
 * Implementation of Subsystem, SubsystemRep and DefaultSystemSubsystem.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/System.h"
#include "simbody/internal/Subsystem.h"
#include "SubsystemRep.h"

#include <cassert>

namespace SimTK {

    ///////////////
    // Subsystem //
    ///////////////

bool Subsystem::isEmptyHandle() const {return rep==0;}
bool Subsystem::isOwnerHandle() const {return rep==0 || rep->myHandle==this;}


Subsystem::~Subsystem() {
    if (isOwnerHandle()) delete rep; 
    rep=0;
}

Subsystem::Subsystem(const Subsystem& src) : rep(0) {
    if (src.rep) {
		assert(!src.rep->isInSystem()); // TODO
        rep = src.rep->clone();	// create a new object
        rep->setMyHandle(*this);
    }
}

Subsystem& Subsystem::operator=(const Subsystem& src) {
    if (&src != this) {
        if (isOwnerHandle()) delete rep; 
        rep=0;
        if (src.rep) {
			assert(!src.rep->isInSystem()); // TODO
			rep = src.rep->clone();	// create a new object
			rep->setMyHandle(*this);
        }
    }
    return *this;
}



const String& Subsystem::getName()    const {return getRep().getName();}
const String& Subsystem::getVersion() const {return getRep().getVersion();}

void Subsystem::realize(const State& s, Stage g) const {
    getRep().realize(s,g);
}

void Subsystem::endConstruction() {updRep().endConstruction();}

bool Subsystem::isInSystem() const {return getRep().isInSystem();}
bool Subsystem::isInSameSystem(const System& sys) const {
	return getRep().isInSameSystem(sys);
}
const System& Subsystem::getSystem() const {return getRep().getSystem();}
System&       Subsystem::updSystem()	   {return updRep().updSystem();}
int Subsystem::getMySubsystemIndex() const {return getRep().getMySubsystemIndex();}

const Vector& Subsystem::getQ(const State& s) const {return getRep().getQ(s);}
const Vector& Subsystem::getU(const State& s) const {return getRep().getU(s);}
const Vector& Subsystem::getZ(const State& s) const {return getRep().getZ(s);}

const Vector& Subsystem::getQErr(const State& s)    const {return getRep().getQErr(s);}
const Vector& Subsystem::getUErr(const State& s)    const {return getRep().getUErr(s);}
const Vector& Subsystem::getUDotErr(const State& s) const {return getRep().getUDotErr(s);}

void Subsystem::calcQUnitWeights(const State& s, Vector& weights) const {
    getRep().calcQUnitWeights(s,weights);
}
void Subsystem::calcUUnitWeights(const State& s, Vector& weights) const {
    getRep().calcUUnitWeights(s,weights);
}
void Subsystem::calcZUnitWeights(const State& s, Vector& weights) const {
    getRep().calcZUnitWeights(s,weights);
}
void Subsystem::calcQErrUnitTolerances(const State& s, Vector& tolerances) const {
    getRep().calcQErrUnitTolerances(s,tolerances);
}
void Subsystem::calcUErrUnitTolerances(const State& s, Vector& tolerances) const {
    getRep().calcUErrUnitTolerances(s,tolerances);
}

    //////////////////
    // SubsystemRep //
    //////////////////

void SubsystemRep::realize(const State& s, Stage g) const {
    while (getStage(s) < g) {
        switch (getStage(s)) {
        case Stage::Empty: {
            State& mutableState = const_cast<State&>(s);
            mutableState.initializeSubsystem(
                getMySubsystemIndex(), getName(), getVersion());
            realizeTopology(mutableState); 
            break;
        }
        case Stage::Topology:     realizeModel(const_cast<State&>(s)); break;
        case Stage::Model:        realizeInstance(s);     break;
        case Stage::Instance:     realizeTime(s);         break;
        case Stage::Time:         realizePosition(s);     break;
        case Stage::Position:     realizeVelocity(s);     break;
        case Stage::Velocity:     realizeDynamics(s);     break;
        case Stage::Dynamics:     realizeAcceleration(s); break;
        case Stage::Acceleration: realizeReport(s);       break;
        default: assert(!"Subsystem::realize(): bad stage");
        }
        advanceToStage(s, getStage(s).next());
    }
}
    ////////////////////////////
    // DefaultSystemSubsystem //
    ////////////////////////////


DefaultSystemSubsystem::DefaultSystemSubsystem() {
    rep = new DefaultSystemSubsystemRep();
    rep->setMyHandle(*this);
}

DefaultSystemSubsystem::DefaultSystemSubsystem(const String& sysName, const String& sysVersion) {
    rep = new DefaultSystemSubsystemRep();
    rep->setMyHandle(*this);
}

/*static*/ bool
DefaultSystemSubsystem::isInstanceOf(const Subsystem& s) {
    return DefaultSystemSubsystemRep::isA(s.getRep());
}

/*static*/ const DefaultSystemSubsystem&
DefaultSystemSubsystem::downcast(const Subsystem& s) {
    assert(DefaultSystemSubsystemRep::isA(s.getRep()));
    return reinterpret_cast<const DefaultSystemSubsystem&>(s);
}

/*static*/ DefaultSystemSubsystem&
DefaultSystemSubsystem::updDowncast(Subsystem& s) {
    assert(DefaultSystemSubsystemRep::isA(s.getRep()));
    return reinterpret_cast<DefaultSystemSubsystem&>(s);
}

DefaultSystemSubsystemRep& DefaultSystemSubsystem::updRep() {
    return DefaultSystemSubsystemRep::downcast(*rep);
}

const DefaultSystemSubsystemRep& DefaultSystemSubsystem::getRep() const {
    return DefaultSystemSubsystemRep::downcast(*rep);
}

} // namespace SimTK

