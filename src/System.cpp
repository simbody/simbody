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

void System::realize(const State& s, Stage g) const {
    getRep().realize(s,g);
}

    ///////////////
    // SystemRep //
    ///////////////

void SystemRep::realize(const State& s, Stage g) const {
    while (s.getStage() < g) {
        switch (s.getStage()) {
        case Stage::Allocated:    realizeConstruction(const_cast<State&>(s)); break;
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
        rep = src.rep->clone();
        rep->setMyHandle(*this);
    }
}

Subsystem& Subsystem::operator=(const Subsystem& src) {
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

void Subsystem::endConstruction() {
    updRep().endConstruction();
}

void Subsystem::realizeConstruction(State& s) const {
    getRep().realizeConstruction(s);
}

void Subsystem::realizeModeling(State& s) const {
    getRep().realizeModeling(s);
}





} // namespace SimTK

