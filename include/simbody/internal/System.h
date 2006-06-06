#ifndef SimTK_SIMBODY_SYSTEM_H_
#define SimTK_SIMBODY_SYSTEM_H_

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
#include "simbody/internal/State.h"

namespace SimTK {

class AnalyticGeometry;
class DecorativeGeometry;
class System;
class Subsystem;
class Study;

/**
 * The abstract parent of all Systems.
 * A System serves as a mediator for a group of interacting Subsystems.
 * All will share a single system State, and typically subsystems will
 * need access to content in the state which is produced by other
 * subsystems. 
 *
 * A System provides a unique subsystem index (a small positive integer)
 * for each of its subsystems, and the subsystems are constructed
 * knowing their indices. The indices are used subsequently by the subsystems
 * to find their own entries in the system state, and by each subsystem
 * to refer to others within the same system. Index 0 is reserved for 
 * use by the System itself, e.g. for system-global state variables.
 *
 * Concrete Systems understand the kinds of subsystems they contain. 
 * For example, a MultibodySystem might contain a mechanical subsystem,
 * a force subsystem, and a geometry subsystem. At each computation
 * stage, a subsystem is realized in a single operation. That operation
 * can refer to computations from already-realized subsystems, but
 * cannot initiate computation in other subsystems. The System must
 * know the proper order with which to realize the subsystems at each
 * stage, and that ordering is likely to vary with stage. For example,
 * at configuration stage the mechanical positions must be realized
 * before the configuration-dependent force elements. However, at
 * reaction stage, the force elements must be realized before the
 * mechanical accelerations can be calculated.
 */
class SimTK_SIMBODY_API System {
public:
    System() : rep(0) { }
    ~System();
    System(const System&);
    System& operator=(const System&);

    const String& getName()    const;
    const String& getVersion() const;

    /// Realize the entire System to the indicated Stage.
    void realize(const State& s, Stage g) const;

    /// Take over ownership of the supplied subsystem and install it into 
    /// the indicated subsystem slot, which must already exist and not
    /// have anything in it. A reference to the new handle is returned,
    /// exactly as though updSubsystem(subsys) had been called.
    Subsystem& takeOverSubsystem(int subsys, Subsystem& src);

    /// How may Subsystems are in here?
    int getNSubsystems() const;
    /// Obtain read-only access to a particular subsystem by its index.
    const Subsystem& getSubsystem(int)   const;
    /// Obtain writable access to a particular subsystem by its index.
    Subsystem&       updSubsystem(int);

    // Internal use only
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;
    explicit System(class SystemRep* r) : rep(r) { }
    bool          hasRep() const {return rep!=0;}
    const SystemRep& getRep() const {assert(rep); return *rep;}
    SystemRep&       updRep() const {assert(rep); return *rep;}
protected:
    class SystemRep* rep;
};

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


/// The abstract parent of all Studies.
class SimTK_SIMBODY_API Study {
public:
    Study() : rep(0) { }
    ~Study();
    Study(const Study&);
    Study& operator=(const Study&);

    Study(const System& sys);

    const System& getSystem() const;
    const State&  getState()  const;
    State&        updState();

    /// Is this handle the owner of this rep? This is true if the
    /// handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // Internal use only
    explicit Study(class StudyRep* r) : rep(r) { }
    bool            hasRep() const {return rep!=0;}
    const StudyRep& getRep() const {assert(rep); return *rep;}
protected:
    class StudyRep* rep;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_SYSTEM_H_
