#ifndef SimTK_SIMBODY_SYSTEM_REP_H_
#define SimTK_SIMBODY_SYSTEM_REP_H_

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
#include "simbody/internal/System.h"

namespace SimTK {

class AnalyticGeometry;
class DecorativeGeometry;

// Stubs.

class SystemRep {
public:
    SystemRep() : myHandle(0) { }

    virtual ~SystemRep() {clearMyHandle();}

    SystemRep* clone() const {
        SystemRep* dup = cloneSystemRep();
        dup->myHandle = 0;
        return dup;
    }

    virtual SystemRep* cloneSystemRep() const = 0;

    virtual void realizeConstruction (State& s)       const = 0;
    virtual void realizeModeling     (State& s)       const = 0;
    virtual void realizeParameters   (const State& s) const { }
    virtual void realizeTime         (const State& s) const { }
    virtual void realizeConfiguration(const State& s) const { }
    virtual void realizeMotion       (const State& s) const { }
    virtual void realizeDynamics     (const State& s) const { }
    virtual void realizeReaction     (const State& s) const { }

    void realize(const State& s, Stage g) const;

    void setMyHandle(System& h) {myHandle = &h;}
    void clearMyHandle() {myHandle=0;}
private:
    friend class System;
    System* myHandle;     // the owner of this rep
};

class StudyRep {
public:
    StudyRep(const System& sys)
      : myHandle(0), system(new System(sys))
    {
        system->realize(state, Stage::Built);
    }

    virtual ~StudyRep() {
        delete system;
    }

    StudyRep* clone() const {
        StudyRep* dup = cloneStudyRep();
        dup->myHandle = 0;
        return dup;
    }
    virtual StudyRep* cloneStudyRep() const = 0;

    const System& getSystem() const {return *system;}
    const State&  getState()  const {return state;}
    State&        updState()        {return state;}

    void setMyHandle(Study& h) {myHandle = &h;}
    void clearMyHandle() {myHandle=0;}
private:
    friend class Study;
    Study* myHandle;     // the owner of this rep

    System* system;
    State   state;
};

class SubsystemRep {
public:
    virtual ~SubsystemRep() { }

    SubsystemRep* clone() const {
        SubsystemRep* dup = cloneSubsystemRep();
        dup->myHandle = 0;
        return dup;
    }

    virtual SubsystemRep* cloneSubsystemRep() const = 0;
    virtual void endConstruction() { }
    virtual void realizeConstruction(State&) const = 0;
    virtual void realizeModeling(State&) const = 0;

    void setMyHandle(Subsystem& h) {myHandle = &h;}
    void clearMyHandle() {myHandle=0;}
private:
    friend class Subsystem;
    Subsystem* myHandle;     // the owner of this rep
};


} // namespace SimTK

#endif // SimTK_SIMBODY_SYSTEM_REP_H_
