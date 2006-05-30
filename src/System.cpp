/* Copyright (c) 2006 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


/**@file
 *
 * Implementation of System and SystemRep.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/System.h"
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

void System::realize(const State& s, Stage g) const {
    getRep().realize(s,g);
}

    ///////////////
    // SystemRep //
    ///////////////

void SystemRep::realize(const State& s, Stage g) const {
    while (s.getStage() < g) {
        switch (s.getStage()) {

        case Stage::Allocated: {
            // The State has nothing in it. We expect the
            // first few discrete variables to be numbered from
            // index=0 so they match the subsystem index.
            State& mutableState = const_cast<State&>(s);
            for (int i=0; i < getNSubsystems(); ++i) {
                int index = mutableState.allocateDiscreteVariable(Stage::Built,
                    new Value<SubsystemDescriptor>(SubsystemDescriptor(i, 
                                                    getSubsystem(i).getName(), 
                                                    getSubsystem(i).getVersion())));

                SimTK_ASSERT2_ALWAYS(index==i, 
                    "SystemRep::realize(): expected discrete variable index %d but got %d",
                    i, index);
            }
            realizeConstruction(mutableState); 
            break;
        }

        case Stage::Built:        realizeModeling    (const_cast<State&>(s)); break;
        case Stage::Modeled:      realizeParameters(s);    break;
        case Stage::Parametrized: realizeTime(s);          break;
        case Stage::Timed:        realizeConfiguration(s); break;
        case Stage::Configured:   realizeMotion(s);        break;
        case Stage::Moving:       realizeDynamics(s);      break;
        case Stage::Dynamics:     realizeReaction(s);      break;
        default: assert(!"System::realize(): bad stage");
        }
        s.advanceToStage(s.getStage().next());
    }
}

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

void Subsystem::endConstruction() {updRep().endConstruction();}
void Subsystem::realizeConstruction(State& s) const {
    getRep().realizeConstruction(s);
}
void Subsystem::realizeModeling(State& s) const {
    getRep().realizeModeling(s);
}
bool Subsystem::isInSystem() const {return getRep().isInSystem();}
bool Subsystem::isInSameSystem(const System& sys) const {
	return getRep().isInSameSystem(sys);
}
const System& Subsystem::getSystem() const {return getRep().getSystem();}
System&       Subsystem::updSystem()	   {return updRep().updSystem();}
int Subsystem::getMySubsystemIndex() const {return getRep().getMySubsystemIndex();}


    //////////////////
    // SubsystemRep //
    //////////////////

// nothing here yet

    /////////////////////////
    // SubsystemDescriptor //
    /////////////////////////

std::ostream& 
operator<<(std::ostream& o, const SubsystemDescriptor& sd) {
    o << "SubsystemDescriptor:\n"
      << "  index="    << sd.subsystemIndex << "\n" 
      << "  name='"    << sd.name << "'\n"
      << "  version='" << sd.version << "'\n"
      << "  constructionVarsIndex=" << sd.constructionVarsIndex << "\n"
      << "  modelingVarsIndex="     << sd.modelingVarsIndex << "\n";
    return o;
}



} // namespace SimTK

