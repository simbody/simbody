#ifndef SimTK_SIMBODY_FORCES_H_
#define SimTK_SIMBODY_FORCES_H_

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

/** @file
 * Define the public interface to the SimbodyForcesSubsystem, a collection
 * of basic force elements suitable for use with a SimbodyMatterSubsystem.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/State.h"
#include "simbody/internal/System.h"
#include "simbody/internal/ForceSubsystem.h"

#include <cassert>

namespace SimTK {

class SimTK_SIMBODY_API SimbodyForcesSubsystem : public ForceSubsystem {
public:
    SimbodyForcesSubsystem();

    // These are just the defaults but are nice to have explicitly for debugging.
    ~SimbodyForcesSubsystem() {
    }
    SimbodyForcesSubsystem(SimbodyForcesSubsystem& e) 
      : ForceSubsystem(e) {
    }
    SimbodyForcesSubsystem& operator=(SimbodyForcesSubsystem& e) {
        ForceSubsystem::operator=(e);
        return *this;
    }

    // This response is available at Stage::Configured.
    const Real& getPotentialEnergy(const State& s) const;

    // These responses are available at Stage::Moving
    const Vector& getMobilityForces(const State& s) const;
    const Vector_<SpatialVec>& getSpatialForces(const State& s) const;

    // Construction: add force elements.
    void addGravity(const Vec3& g);
    void addMobilitySpring(int body, int dof,
                           const Real& stiffness,
                           const Real& naturalLength);
    void addTwoPointSpring(int body1, const Vec3& s1,
                           int body2, const Vec3& s2,
                           const Real&    stiffness,
                           const Real&    naturalLength);
    void addTwoPointDamper(int body1, const Vec3& s1,
                           int body2, const Vec3& s2,
                           const Real&    damping);


    SimTK_PIMPL_DOWNCAST(SimbodyForcesSubsystem, ForceSubsystem);
private:
    class SimbodyForcesSubsystemRep& updRep();
    const SimbodyForcesSubsystemRep& getRep() const;
};

/**
 * This is a concrete subsystem which applies a single linear, two-point spring.
 */
class SimTK_SIMBODY_API TwoPointSpringSubsystem : public ForceSubsystem {
public:
    TwoPointSpringSubsystem(int b1, const Vec3& s1,
                            int b2, const Vec3& s2,
                            const Real& k, const Real& x0);

    // TODO: this doesn't belong here!
    const Vec3& getGravity(const State& s) const;
    Vec3&       updGravity(State& s)       const;

    const Real& getStiffness(const State& s) const;
    Real&       updStiffness(State& s)       const;

    const Real& getNaturalLength(const State& s) const;
    Real&       updNaturalLength(State& s)       const;

    const Real& getPotentialEnergy(const State& s) const;
    const Vec3& getForceOnStation1(const State& s) const;


    // These are just the defaults but are nice to have explicitly for debugging.
    ~TwoPointSpringSubsystem() {
    }
    TwoPointSpringSubsystem(const TwoPointSpringSubsystem& e) 
      : ForceSubsystem(e) {
    }
    TwoPointSpringSubsystem& operator=(const TwoPointSpringSubsystem& e) {
        ForceSubsystem::operator=(e);
        return *this;
    }

    SimTK_PIMPL_DOWNCAST(TwoPointSpringSubsystem, ForceSubsystem);
private:
    class TwoPointSpringSubsystemRep& updRep();
    const TwoPointSpringSubsystemRep& getRep() const;
};

/**
 * This is a concrete subsystem which applies no forces.
 */
class SimTK_SIMBODY_API EmptyForcesSubsystem : public ForceSubsystem {
public:
    EmptyForcesSubsystem();

    // These are just the defaults but are nice to have explicitly for debugging.
    ~EmptyForcesSubsystem() {
    }
    EmptyForcesSubsystem(const EmptyForcesSubsystem& e) 
      : ForceSubsystem(e) {
    }
    EmptyForcesSubsystem& operator=(const EmptyForcesSubsystem& e) {
        ForceSubsystem::operator=(e);
        return *this;
    }

    SimTK_PIMPL_DOWNCAST(EmptyForcesSubsystem, ForceSubsystem);
private:
    class EmptyForcesSubsystemRep& updRep();
    const EmptyForcesSubsystemRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCES_H_
