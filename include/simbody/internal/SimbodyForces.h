#ifndef SimTK_SIMBODY_FORCES_H_
#define SimTK_SIMBODY_FORCES_H_

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

    const Real& getDamping(const State& s) const;
    Real&       updDamping(State& s)       const;

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
 * This is a concrete subsystem which applies a uniform gravity field
 * to all the matter in the system. The gravity vector is given in
 * the ground frame. State variables (all Parameter stage):
 *    enabled    boolean, default true
 *    g          Vec3, the gravity vector, default 0
 *    zeroHeight Real, affects potential energy
 */
class SimTK_SIMBODY_API UniformGravitySubsystem : public ForceSubsystem {
public:
    UniformGravitySubsystem();
    explicit UniformGravitySubsystem(const Vec3& g, const Real& zeroHeight=0);

    /// State variables can be accessed after Modeled stage.
    const Vec3& getGravity(const State&) const;
    const Real& getZeroHeight(const State&) const;
    bool        isEnabled(const State&) const;

    /// Solver for setting the gravity vector in the state. Callable
    /// at Modeled stage or higher; backs stage up to Modeled.
    Vec3& updGravity(State& s) const;
    const UniformGravitySubsystem& setGravity(State& s, const Vec3& g) const {
        updGravity(s)=g;
        return *this;
    }

    Real& updZeroHeight(State& s) const;
    const UniformGravitySubsystem& setZeroHeight(State&s, const Real& z) const {
        updZeroHeight(s)=z;
        return *this;
    }

    /// Solvers for enabling or disabling gravity in the state. Callable
    /// at Modeled stage or higher; backs stage up to Modeled.
    bool& updIsEnabled(State& s) const;
    const UniformGravitySubsystem& enableGravity(State& s) const {
        updIsEnabled(s) = true;
        return *this;
    }
    const UniformGravitySubsystem& disableGravity(State& s) const{
        updIsEnabled(s) = false;
        return *this;
    }

    // These are just the defaults but are nice to have explicitly for debugging.
    ~UniformGravitySubsystem() {
    }
    UniformGravitySubsystem(const UniformGravitySubsystem& e) 
      : ForceSubsystem(e) {
    }
    UniformGravitySubsystem& operator=(const UniformGravitySubsystem& e) {
        ForceSubsystem::operator=(e);
        return *this;
    }

    SimTK_PIMPL_DOWNCAST(UniformGravitySubsystem, ForceSubsystem);
private:
    class UniformGravitySubsystemRep& updRep();
    const UniformGravitySubsystemRep& getRep() const;
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
