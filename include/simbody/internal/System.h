#ifndef SimTK_SYSTEM_H_
#define SimTK_SYSTEM_H_

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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"

namespace SimTK {

class State;
class Subsystem;

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
class SimTK_SIMBODY_EXPORT System {
public:
    System() : rep(0) { }
    ~System();
    System(const System&);
    System& operator=(const System&);

    const String& getName()    const;
    const String& getVersion() const;

    /// Realize the entire System to the indicated Stage.
    void realize(const State& s, Stage g = Stage::HighestValid) const;

    /// This operator can be called at Stage::Instance or higher and 
    /// returns a rough estimate of a length of time we consider significant
    /// for this system. For example, this could be the period of the highest-frequency
    /// oscillation that we care about. This can be used as a hint by 
    /// numerical integrators in choosing their initial step size, and 
    /// suggests how velocity variables should be scaled relative to their
    /// corresponding position variables.
    Real calcTimescale(const State&) const;

    /// This operator can be called at Stage::Position to take a vector
    /// of absolute state variable error estimates and return a weighted
    /// norm. This method is intended for use by numerical integration methods
    /// for step size control. This is a weighted norm, calculated so
    /// that a return value of 1 would indicate a "unit" error, which would 
    /// be huge. If your accuracy requirement is 0.1%, you would test that
    /// the norm return here is <= .001.
    Real calcYErrorNorm(const State&, const Vector& y_err) const;

    /// Take over ownership of the supplied subsystem and install it into 
    /// the next free subsystem slot. The new slot index is returned.
    int takeOverSubsystem(Subsystem& src);

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


/// The abstract parent of all Studies.
class SimTK_SIMBODY_EXPORT Study {
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

#endif // SimTK_SYSTEM_H_
