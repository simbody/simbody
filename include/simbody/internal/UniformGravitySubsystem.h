#ifndef SimTK_UNIFORM_GRAVITY_SUBSYSTEM_H_
#define SimTK_UNIFORM_GRAVITY_SUBSYSTEM_H_

/* Portions copyright (c) 2006-7 Stanford University and Michael Sherman.
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
 * This file declares the public interface to the UniformGravitySubsystem, a
 * concrete ForceSubsystem that implements a uniform gravitational acceleration
 * field applied to all the matter in the containing MultibodySystem.
 * Add this to any MultibodySystem as a force subsystem, set the gravity
 * vector (expressed in ground) and all the bodies will fall in that direction.
 * Each body will experience a force through its center of mass equal to
 * mg where m is the body's mass and g is the gravity vector.
 *
 * This subsystem also contributes potential energy to the system for each
 * body, with the scalar value m|g|h where |g| is the magnitude (length) of
 * the gravity vector and h is the height of the body's center of mass
 * measured in the -g direction from some arbitrary reference height. By
 * default the reference height is zero, but you can set it here to some
 * other value, which has the sole effect of subtracting m|g|z from each
 * body's potential energy, where z is the (constant) zero height.
 *
 * Finally, you can independently enable or disable gravity without changing
 * the value of g or the zero height.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/ForceSubsystem.h"

#include <cassert>

namespace SimTK {

class State;
class MultibodySystem;

/**
 * This is a concrete subsystem which applies a uniform gravity field
 * to all the matter in the system. The gravity vector is given in
 * the ground frame. State variables (all Instance stage):
 *    enabled    boolean, default true
 *    g          Vec3, the gravity vector, default 0
 *    zeroHeight Real, affects potential energy only, default 0
 */
class SimTK_SIMBODY_EXPORT UniformGravitySubsystem : public ForceSubsystem {
public:
    UniformGravitySubsystem();
    explicit UniformGravitySubsystem(MultibodySystem&);
    UniformGravitySubsystem(MultibodySystem&, const Vec3& g, const Real& zeroHeight=0);

    /// State variables can be accessed after Model stage.
    const Vec3& getGravity(const State&) const;
    const Real& getZeroHeight(const State&) const;
    bool        isEnabled(const State&) const;

    /// Solver for setting the gravity vector in the state. Callable
    /// at Model stage or higher; backs stage up to Model.
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
    /// at Model stage or higher; backs stage up to Model.
    bool& updIsEnabled(State& s) const;
    const UniformGravitySubsystem& enableGravity(State& s) const {
        updIsEnabled(s) = true;
        return *this;
    }
    const UniformGravitySubsystem& disableGravity(State& s) const{
        updIsEnabled(s) = false;
        return *this;
    }

    // private implementation
    SimTK_PIMPL_DOWNCAST(UniformGravitySubsystem, ForceSubsystem);
    class UniformGravitySubsystemRep& updRep();
    const UniformGravitySubsystemRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_UNIFORM_GRAVITY_SUBSYSTEM_H_
