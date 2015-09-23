#ifndef SimTK_SIMBODY_FORCE_DISCRETE_FORCES_H_
#define SimTK_SIMBODY_FORCE_DISCRETE_FORCES_H_

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
subclass Force::DiscreteForces and is logically part of Force.h. The file assumes that
Force.h will have included all necessary declarations. **/

namespace SimTK {

/** Arbitrary discrete body forces and mobility (generalized) forces.\ Useful
for applying external forces or forces that are updated at discrete times due
to the occurrence of events.
@see Force::MobilityDiscreteForce **/
class SimTK_SIMBODY_EXPORT Force::DiscreteForces : public Force {
public:
    /** Create a %DiscreteForces force element.

    @param forces   Subsystem to which this force element should be added.
    @param matter   Subsystem containing the mobilized bodies to which these
                    forces will be applied.
    **/
    DiscreteForces(GeneralForceSubsystem&           forces,
                   const SimbodyMatterSubsystem&    matter);


    /** Default constructor creates an empty handle. **/
    DiscreteForces() {}

    /** Set to zero all forces stored by this force element in the given
    \a state, including both generalized forces and body spatial forces. **/
    void clearAllForces(State& state) const {
        clearAllMobilityForces(state);
        clearAllBodyForces(state);
    }

    /** Set to zero all generalized forces stored by this force element in the
    given \a state. **/
    void clearAllMobilityForces(State& state) const;
    /** Set to zero all body spatial forces (force and torques) stored by this
    force element in the given \a state. **/
    void clearAllBodyForces(State& state) const;

    /** Change the value of a generalized force to be applied in the given
    \a state. Set this to zero if you don't want it to do anything. You can
    call this method any time after realizeTopology() since the force just needs
    to be saved in the \a state.

    @param  state   The state in which the generalized force is to be saved.
    @param  mobod   The mobilizer to which the force should be applied.
    @param  whichU  To which of the mobilizer's mobilities (degrees of
                    freedom) should this force be applied (first is 0)?
    @param  force   The scalar value of the generalized force.

    @see getOneMobilityForce() **/
    void setOneMobilityForce(State& state, const MobilizedBody& mobod,
                             MobilizerUIndex whichU, Real force) const;

    /** Change the value of the discrete spatial force (force and torque) to
    be applied to a body in the given \a state. Set this to zero if you don't
    want any force applied to the body. You can call this method any time after
    realizeTopology() since the force just needs to be saved in the \a state.

    @param      state
        The state in which the force is to be saved.
    @param      mobod
        The mobilized body to which the force should be applied.
    @param      spatialForceInG
        The (torque,force) pair as a SpatialVec. These are expressed in the
        Ground frame and the force is applied to the body origin.

    @see getOneBodyForce(), addForceToBodyPoint() **/
    void setOneBodyForce(State& state, const MobilizedBody& mobod,
                         const SpatialVec& spatialForceInG) const;

    /** Convenience method to \e add in a force applied at a point (station) on
    a particular body into the forces currently set to be applied to that body
    in this \a state. The station location is given in the body frame and the
    force is given in the Ground frame. Note that the \a state must have been
    realized to Stage::Position because we need to know the body orientation
    in order to shift the applied force to the body origin.

    Be sure to call setOneBodyForce() to initialize this force to zero before
    you start accumulating point forces using this method.

    @param  state   The state in which the force is saved. Must already have
                    been realized to at least Stage::Position.
    @param  mobod   The mobilized body to which the force will be applied.
    @param  pointInB The body station at which the force will be applied,
                     measured and expressed in the body frame.
    @param  forceInG The force vector to be applied, expressed in Ground.

    @see setOneBodyForce() **/
    void addForceToBodyPoint(State& state, const MobilizedBody& mobod,
                             const Vec3& pointInB,
                             const Vec3& forceInG) const;

    /** Set all the generalized forces at once. These will be applied until
    they are changed.
    @param  state   The state in which the generalized forces are saved.
    @param  f       The set of n generalized forces. This Vector must be either
                    empty (in which case it is treated as though it were
                    all zero) or the same length as the number of generalized
                    speeds (state.getNU()).
    @see clearAllMobilityForces() **/
    void setAllMobilityForces(State& state, const Vector& f) const;

    /** Set all the body spatial forces at once. These will be applied until
    they are changed.
    @param  state   The state in which the forces are saved.
    @param  bodyForcesInG
        The set of nb spatial forces (torque,force) pairs, expressed in the
        Ground frame, with the forces applied at each body frame origin. This
        Vector must be either empty (in which case it is treated as though it
        were all zero) or the same length as the number of mobilized bodies in
        the system (and don't forget that Ground is mobilized body 0).
    @see clearAllBodyForces() **/
    void setAllBodyForces(State& state,
                          const Vector_<SpatialVec>& bodyForcesInG) const;

    /** Return the value for a generalized force that is stored in the
    given \a state. If no calls to setGeneralizedForce() have been made on this
    \a state then zero will be returned.
    @see setOneMobilityForce() **/
    Real getOneMobilityForce(const State& state, const MobilizedBody& mobod,
                                MobilizerUIndex whichU) const;

    /** Return the value for the discrete spatial force (force and torque) on a
    particular body that is stored in the given \a state. The result is a force
    applied at the body origin and expressed in the Ground frame.
    @see setOneBodyForce() **/
    SpatialVec getOneBodyForce(const State& state,
                               const MobilizedBody& mobod) const;

    /** Get a reference to the internal Vector of all the generalized forces
    currently set to be applied by this force element. This will be length zero
    if no forces are being applied; otherwise, its length will be the number
    of generalized speeds u in the given \a state.
    The returned reference is a reference into the given \a state object. **/
    const Vector& getAllMobilityForces(const State& state) const;
    /** Get a reference to the internal Vector of all the body spatial forces
    currently set to be applied by this force element. This will be length zero
    if no body forces are being applied; otherwise, its length will be the
    number of mobilized bodies in the containing MultibodySystem.
    The returned reference is a reference into the given \a state object. **/
    const Vector_<SpatialVec>& getAllBodyForces(const State& state) const;

    /** @cond **/
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(DiscreteForces,
                                             DiscreteForcesImpl, Force);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_DISCRETE_FORCES_H_
