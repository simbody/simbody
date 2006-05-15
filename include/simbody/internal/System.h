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

/// The abstract parent of all Systems.
class SimTK_SIMBODY_API System {
public:
    System() : rep(0) { }
    ~System();
    System(const System&);
    System& operator=(const System&);

    void realize(const State& s, Stage g) const;

    /// Is this handle the owner of this rep? This is true if the
    /// handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // Internal use only
    explicit System(class SystemRep* r) : rep(r) { }
    bool          hasRep() const {return rep!=0;}
    const SystemRep& getRep() const {assert(rep); return *rep;}
protected:
    class SystemRep* rep;
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

/// The abstract parent of all Subsystems.
class SimTK_SIMBODY_API Subsystem {
public:
    Subsystem() : rep(0) { }
    ~Subsystem();
    Subsystem(const Subsystem&);
    Subsystem& operator=(const Subsystem&);

    void endConstruction();
    void realizeConstruction(State&) const;
    void realizeModeling(State&) const;

    /// Is this handle the owner of this rep? This is true if the
    /// handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // Internal use only
    explicit Subsystem(class SubsystemRep* r) : rep(r) { }
    bool                hasRep() const {return rep!=0;}
    const SubsystemRep& getRep() const {assert(rep); return *rep;}
    SubsystemRep&       updRep() const {assert(rep); return *rep;}
protected:
    class SubsystemRep* rep;
};



} // namespace SimTK

#endif // SimTK_SIMBODY_SYSTEM_H_
