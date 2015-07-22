#ifndef SimTK_SIMBODY_FORCE_MOBILITY_CONSTANT_FORCE_H_
#define SimTK_SIMBODY_FORCE_MOBILITY_CONSTANT_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-13 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon.h"
#include "simbody/internal/Force.h"

/** @file
This contains the user-visible API ("handle" class) for the SimTK::Force
subclass Force::MobilityConstantForce and is logically part of Force.h. The
file assumes that Force.h will have included all necessary declarations. **/

namespace SimTK {

/** A constant generalized force f (a scalar) applied to a mobility. The
mobility here selects a generalized speed (u), not a generalized coordinate (q),
and the meaning depends on the definition of the generalized speed. If that
speed is a translation then this is a force; if a rotation then this is a
torque; if anything more general then this is in compatible force units. This
force does not contribute to the potential energy, so adding it to a system will
cause energy not to be conserved unless you account for the power injected or
dissipated here. **/
class SimTK_SIMBODY_EXPORT Force::MobilityConstantForce : public Force {
public:
    /** Add a %MobilityConstantForce element to the indicated Subsystem.
    @param forces       Subsystem to which this force should be added.
    @param mobod        Mobilizer to which the force should be applied.
    @param whichU       Index within \a mobod of the generalized speed u to
                            which this force should be applied (first is zero).
    @param defaultForce Default value for the generalized force. **/
    MobilityConstantForce(GeneralForceSubsystem&    forces,
                          const MobilizedBody&      mobod,
                          MobilizerUIndex           whichU,
                          Real                      defaultForce);

    /** Alternate constructor signature for when the mobilizer has only
    a single generalized speed, in which case we'll use MobilizerUIndex(0).
    See the other signature for documentation. **/
    MobilityConstantForce(GeneralForceSubsystem&    forces,
                          const MobilizedBody&      mobod,
                          Real                      defaultForce)
    {   // Invoke the other constructor.
        new(this) MobilityConstantForce(forces, mobod, MobilizerUIndex(0),
                                        defaultForce);
    }

    /** Default constructor creates an empty handle. **/
    MobilityConstantForce() {}

    /** Provide a new value for the default generalied force to be applied by
    this force element. This is a topological change because it affects the
    value that the containing System's default state will have when
    realizeTopology() is called. This is for use during construction, not for
    during a simulation where you should be using setForce() to set the force
    in a State rather than in the System.

    @param defaultForce Default value for the generalized force.
    @returns A writable reference to this modified force element.
    @see setForce(), getDefaultForce() **/
    MobilityConstantForce& setDefaultForce(Real defaultForce);

    /** Return the default value for the generalized force. This is normally
    set at construction but may have been changed with setDefaultForce().
    @see setDefaultForce(), getForce() **/
    Real getDefaultForce() const;

    /** Change the value of the generalized force that is stored in the given
    \a state; this may differ from the default value supplied at construction.
    @param  state   The State in which the bounds are changed.
    @param  force   The value of the generalized force to be applied when the
                    \a state is subsequently used.

    Changing this force invalidates Stage::Dynamics and above in the \a state
    since it affects force generation.
    @see getForce(), setDefaultForce() **/
    void setForce(State& state, Real force) const;

    /** Return the value for the generalized force that is stored in the
    given \a state. Note that this is not the same thing as the default force
    that was supplied on construction or with setDefaultForce().
    @see setForce(), getDefaultForce() **/
    Real getForce(const State& state) const;

    /** @name                      Deprecated
    Methods here are for backwards compatibility but have been replaced with
    better ones that you should use. **/
    /**@{**/
    /** Deprecated: Alternate signature for backwards compatibilty -- for
    safety you should prefer using the other constructor signature that
    takes a MobilizerUIndex rather than a plain int. **/
    MobilityConstantForce(GeneralForceSubsystem&    forces,
                          const MobilizedBody&      mobod,
                          int                       whichU,
                          Real                      defaultForce)
    {   // Invoke the other constructor.
        new(this) MobilityConstantForce(forces, mobod, MobilizerUIndex(whichU),
                                        defaultForce);
    }
    /**@}**/

    /** @cond **/ // Hide from Doxygen.
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(MobilityConstantForce,
                                             MobilityConstantForceImpl, Force);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_MOBILITY_CONSTANT_FORCE_H_
