#ifndef SimTK_SIMBODY_FORCE_MOBILITY_DISCRETE_FORCE_H_
#define SimTK_SIMBODY_FORCE_MOBILITY_DISCRETE_FORCE_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
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
subclass Force::MobilityDiscreteForce and is logically part of Force.h. The 
file assumes that Force.h will have included all necessary declarations. **/

namespace SimTK {

/** A discrete mobility (generalized) force f applied to a particular mobility
that is specified at construction.\ Useful for applying
external forces or forces that are updated at discrete times due to the
occurrence of events. Note that a mobility is a generalized speed (u), not 
a generalized coordinate (q). The meaning of a generalized force depends on the
definition of the generalized speed. If that speed is a translation then this is
a force; if a rotation then this is a torque; if something else then f has a 
comparable definition (the defining condition is that f*u should always have 
physically meaningful units of power). This force does not contribute to the 
potential energy, so adding it to a system will cause potential+kinetic energy 
not to be conserved.

If you want to be able to apply discrete forces to any body or mobilizer 
without specifying which one in advance, see Force::DiscreteForces.
@see Force::DiscreteForces **/
class SimTK_SIMBODY_EXPORT Force::MobilityDiscreteForce : public Force {
public:
    /** Create a %MobilityDiscreteForce.
    
    @param forces        subsystem to which this force element should be added
    @param mobod         mobilizer to which the force should be applied
    @param whichU        to which of the mobilizer's mobilities (degrees of
                            freedom) should this force be applied (first is 0)?
    @param defaultForce  initial value for the generalized force to be applied 
                             (default 0)

    Note that if you have an integer value for the generalized speed (u) index,
    you have to cast it to a MobilizerUIndex here. The generalized speeds are
    numbered starting with 0 for each mobilizer. Here is an example:
    @code
        GeneralForceSubsystem forces;
        MobilizedBody::Pin    pinJoint(...);
        MobilityDiscreteForce myForce(forces, pinJoint, MobilizerUIndex(0));
    @endcode
    **/
    MobilityDiscreteForce(GeneralForceSubsystem&    forces, 
                          const MobilizedBody&      mobod, 
                          MobilizerUIndex           whichU, 
                          Real                      defaultForce=0);


    /** Alternate constructor signature for when the mobilizer has only
    a single generalized speed, in which case we'll use MobilizerUIndex(0). 
    See the other signature for documentation. **/
    MobilityDiscreteForce(GeneralForceSubsystem&    forces, 
                          const MobilizedBody&      mobod, 
                          Real                      defaultForce=0)    
    {   // Invoke the other constructor.
        new(this) MobilityDiscreteForce(forces, mobod, MobilizerUIndex(0),
                                        defaultForce);
    }
    
    /** Default constructor creates an empty handle. **/
    MobilityDiscreteForce() {}

    /** Provide a new value for the \a defaultForce, overriding the one 
    provided in the constructor. This is a topological change because it
    affects the value that the containing System's default state will have
    when realizeTopology() is called. This is for use during construction, not
    for during a simulation where you should be using setGeneralizedForce().

    @param defaultForce     the value this generalized force should have by 
                                default
    @returns a writable reference to this modified force element 
    @see setMobilityForce(), getDefaultMobilityForce()  **/
    MobilityDiscreteForce& setDefaultMobilityForce(Real defaultForce);

    /** Return the value that this generalized force will have by default. This
    is normally set in the constructor, or left to its default value of 0. It 
    can also be set in setDefaultMobilityForce(). Note that this is \e not
    the same as the value that may be set in any particular State. 
    @see getMobilityForce(), setDefaultMobilityForce() **/
    Real getDefaultMobilityForce() const;

    /** Change the value of the generalized force to be applied in the given
    \a state. Set this to zero if you don't want it to do anything. 
    @see getMobilityForce() **/
    void setMobilityForce(State& state, Real f) const;
    /** Return the value for this generalized force that is stored in the 
    given \a state. If no calls to setMobilityForce() have been made on this
    \a state then it will have the \a defaultForce value that was supplied on
    construction or via setDefaultMobilityForce(). 
    @see setMobilityForce() **/
    Real getMobilityForce(const State& state) const;

    /** @cond **/
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(MobilityDiscreteForce, 
                                             MobilityDiscreteForceImpl, Force);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_MOBILITY_DISCRETE_FORCE_H_
