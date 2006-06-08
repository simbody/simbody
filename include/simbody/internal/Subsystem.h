#ifndef SimTK_SUBSYSTEM_H_
#define SimTK_SUBSYSTEM_H_

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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"

namespace SimTK {

class State;
class System;

/**
 * The abstract parent of all Subsystems.
 * A Subsystem is expected to be part of a larger System and to have
 * interdependencies with other subsystems of that same system. It
 * must NOT have dependencies on objects which are outside the system.
 * Consequently construction of any concrete subsystem requires
 * specification of a system at that time.
 * Subsystems go through an extended construction phase in which
 * their contents and interdependencies are created. Thus all
 * of a system's subsystems generally need to be available simultaneously 
 * during construction, so that they can reference each other.
 */
class SimTK_SIMBODY_API Subsystem {
public:
    Subsystem() : rep(0) { }
    ~Subsystem();
    Subsystem(const Subsystem&);
    Subsystem& operator=(const Subsystem&);

    const String& getName()    const;
    const String& getVersion() const;

    // Realize the Subsystem to the indicated Stage.
    void realize(const State& s, Stage g) const;

	bool isInSystem() const;
	bool isInSameSystem(const System&) const;
	const System& getSystem() const;
	System&       updSystem();

	int getMySubsystemIndex() const;

    void endConstruction();

    /// Is this handle the owner of this rep? This is true if the
    /// handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // Internal use only
    explicit Subsystem(class SubsystemRep* r) : rep(r) { }
    bool                hasRep() const {return rep!=0;}
    const SubsystemRep& getRep() const {assert(rep); return *rep;}
    SubsystemRep&       updRep() const {assert(rep); return *rep;}
	void setRep(SubsystemRep& r) {assert(!rep); rep = &r;}
protected:
    class SubsystemRep* rep;
};

/**
 * This is a concrete Subsystem used by default as the 0th Subsystem of
 * every System. Feel free to replace it with something useful!
 */
class SimTK_SIMBODY_API DefaultSystemSubsystem : public Subsystem {
public:
    DefaultSystemSubsystem();
    DefaultSystemSubsystem(const String& sysName, const String& sysVersion);

    SimTK_PIMPL_DOWNCAST(DefaultSystemSubsystem, Subsystem);
private:
    class DefaultSystemSubsystemRep& updRep();
    const DefaultSystemSubsystemRep& getRep() const;
};


} // namespace SimTK

#endif // SimTK_SUBSYSTEM_H_
