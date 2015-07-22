#ifndef SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_SPRING_H_
#define SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_SPRING_H_

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
subclass Force::MobilityLinearSpring and is logically part of Force.h. The file
assumes that Force.h will have included all necessary declarations. **/

namespace SimTK {

/** A linear spring that acts along or around a mobility coordinate to apply
a generalized force there.

The stiffness k is provided, along with an arbitrary "zero" coordinate value
q0 at which the spring generates no force. The generated force is k*(q-q0), and
potential energy is pe = 1/2 k (q-q0)^2.

@bug This is not meaningful unless the mobilizer coordinates are defined such
that qdot=u. In particular, do not use this on coordinates for Ball or Free
mobilizers. You can often replace those with Gimbal or Bushing mobilizers for
which qdot=u holds, or split up the mobilizer into simpler components separated
by massless bodies. **/
class SimTK_SIMBODY_EXPORT Force::MobilityLinearSpring : public Force {
public:
    /** Create a %MobilityLinearSpring force element on a particular generalized
    coordinate.

    @param[in,out]  forces
        The subsystem to which this force should be added.
    @param[in]      mobod
        Mobilizer to which the force should be applied.
    @param[in]      whichQ
        To which of the mobilizer's generalized coordinates q should this
        force be applied (first is 0)?
    @param[in]      defaultStiffness
        The default value for the spring constant k.
    @param[in]      defaultQZero
        The default for the value of the coordinate q0 at which the force is 0.
    **/
    MobilityLinearSpring(GeneralForceSubsystem& forces,
                         const MobilizedBody&   mobod,
                         MobilizerQIndex        whichQ,
                         Real                   defaultStiffness,
                         Real                   defaultQZero);

    /** Default constructor creates an empty handle that can be assigned to
    refer to any %MobilityLinearSpring object. **/
    MobilityLinearSpring() {}

    /** Provide a new value for the default stiffness k of this spring.
    This is a topological change because it affects the value that the
    containing System's default state will have when realizeTopology() is
    called. This is for use during construction, not for during a simulation
    where you should be using setStiffness() to set the stiffness in a State
    rather than in the System.
    @param[in]      defaultStiffness
        The default value for the spring constant k.
    @return
        A writable reference to this modified force element for convenience in
        chaining set methods.
    @see getDefaultStiffness(), setStiffness() **/
    MobilityLinearSpring& setDefaultStiffness(Real defaultStiffness);


    /** Provide a new value for the zero position q0 of this spring, at which
    position the spring force will be zero.
    This is a topological change because it affects the value that the
    containing System's default state will have when realizeTopology() is
    called. This is for use during construction, not for during a simulation
    where you should be using setQZero() to set the zero position in a State
    rather than in the System.
    @param[in]      defaultQZero
        The default for the value of the coordinate at which the force is 0.
    @return
        A writable reference to this modified force element for convenience in
        chaining set methods.
    @see getDefaultQZero(), setQZero() **/
    MobilityLinearSpring& setDefaultQZero(Real defaultQZero);

    /** Return the default value for the spring's stiffness k. This is
    normally set at construction but can be modified with setDefaultStiffness().
    @see setDefaultStiffness(), getStiffness() **/
    Real getDefaultStiffness() const;
    /** Return the default value for the spring's zero position q0. This is
    normally set at construction but can be modified with setDefaultQZero().
    @see setDefaultQZero(), getQZero() **/
    Real getDefaultQZero() const;

    /** Change the value of the spring stiffness in the given \a state; this may
    differ from the default value supplied at construction.
    @param[in,out]  state
        The State in which the stiffness is to be changed.
    @param[in]      stiffness
        The new stiffness k (>= 0) that overrides the default.
    @return
        A const reference to this %MobilityLinearSpring element for convenience
        in chaining set methods together.

    Changing the spring stiffness invalidates Stage::Dynamics and above in the
    \a state since it can affect force generation.
    @see setDefaultStiffness(), getStiffness() **/
    const MobilityLinearSpring& setStiffness(State&     state,
                                             Real       stiffness) const;

    /** Change the value of the spring zero length in the given \a state; this
    may differ from the default value supplied at construction.
    @param[in,out]  state
        The State in which the zero length is to be changed.
    @param[in]      qZero
        The value of the controlled coordinate q at which the generated spring
        force should be zero. This overrides the default.
    @return
        A const reference to this %MobilityLinearSpring element for convenience
        in chaining set methods together.

    Changing the spring stiffness invalidates Stage::Dynamics and above in the
    \a state since it can affect force generation.
    @see setDefaultStiffness(), getStiffness() **/
    const MobilityLinearSpring& setQZero(State&     state,
                                         Real       qZero) const;

    /** Return the value for the spring's stiffness k that is stored in
    the given \a state. Note that this is not the same thing as the default
    stiffness that was supplied on construction or in setDefaultStiffness().
    @see setStiffness(), getDefaultStiffness() **/
    Real getStiffness(const State& state) const;

    /** Return the value for the spring zero position q0 that is stored in
    the given \a state. Note that this is not the same thing as the default
    q0 that was supplied on construction or in setDefaultQZero().
    @see setQZero(), getDefaultQZero() **/
    Real getQZero(const State& state) const;

    /** @name                      Deprecated
    Methods here are for backwards compatibility but have been replaced with
    better ones that you should use. **/
    /**@{**/
    /** Deprecated: Alternate signature for backwards compatibilty -- for
    safety you should prefer using the other constructor signature that
    takes a MobilizerQIndex rather than a plain int. **/
    MobilityLinearSpring(GeneralForceSubsystem& forces,
                         const MobilizedBody&   mobod,
                         int                    whichQ,
                         Real                   defaultStiffness,
                         Real                   defaultQZero)
    {   // Invoke the other constructor.
        new(this) MobilityLinearSpring(forces, mobod, MobilizerQIndex(whichQ),
                                       defaultStiffness, defaultQZero);
    }
    /**@}**/

    /** @cond **/
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(MobilityLinearSpring,
                                             MobilityLinearSpringImpl, Force);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_SPRING_H_
