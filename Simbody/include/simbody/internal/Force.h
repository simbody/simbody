#ifndef SimTK_SIMBODY_FORCE_H_
#define SimTK_SIMBODY_FORCE_H_

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
#include "simbody/internal/common.h"
#include "simbody/internal/GeneralForceSubsystem.h"

namespace SimTK {

class SimbodyMatterSubsystem;
class GeneralForceSubsystem;
class MobilizedBody;
class Force;
class ForceImpl;

// We only want the template instantiation to occur once. This symbol is defined
// in the SimTK core compilation unit that defines the Force class but should 
// not be defined any other time.
#ifndef SimTK_SIMBODY_DEFINING_FORCE
    extern template class PIMPLHandle<Force, ForceImpl, true>;
#endif

/** This is the base class from which all Force element handle classes derive.
A Force object applies forces to some or all of the bodies, particles, and 
mobilities in a System. There are subclasses for various standard types of 
forces, or you can create your own forces by deriving from Force::Custom. **/
class SimTK_SIMBODY_EXPORT Force : public PIMPLHandle<Force, ForceImpl, true> {
public:
    /**@name                   Enabling and disabling
    These methods determine whether this force element is active in a given
    State. When disabled, the Force element is completely ignored and will
    not be updated during realization. Normally force elements are enabled
    when defined unless explicitly disabled; you can reverse that using the
    setDisabledByDefault() method below. **/
    /*@{*/
    /** Disable this force element, effectively removing it from the System
    for computational purposes (it is still using its ForceIndex, however).
    This is an Instance-stage change. **/
    void disable(State&) const;
    /** Enable this force element if it was previously disabled. This is an 
    Instance-stage change. Nothing happens if the force element was already
    enabled. **/
    void enable(State&) const;
    /** Test whether this force element is currently disabled in the supplied 
    State. If it is disabled you cannot depend on any computations it 
    normally performs being available. **/
    bool isDisabled(const State&) const;
    /** Normally force elements are enabled when defined and can be disabled 
    later. If you want to define this force element but have it be off by 
    default, use this method. Note that this is a Topology-stage (construction)
    change; you will have to call realizeTopology() before using the containing
    System after a change to this setting has been made. **/
    void setDisabledByDefault(bool shouldBeDisabled);
    /** Test whether this force element is disabled by default in which case it
    must be explicitly enabled before it will take effect.
    @see enable() **/
    bool isDisabledByDefault() const;
    /*@}*/

    /**@name                   Advanced methods
    Don't use these unless you're sure you know what you're doing. They aren't
    normally necessary but can be handy sometimes, especially when debugging
    newly-developed force elements. **/
    /*@{*/
    /** Calculate the force that would be applied by this force element if
    the given \a state were realized to Dynamics stage. This sizes the given
    arrays if necessary, zeroes them, and then calls the force element's
    calcForce() method which adds its force contributions if any to the
    appropriate array elements for bodies, particles, and mobilities. Note that
    in general we have no idea what elements of the system are affected by a 
    force element, and in fact that can change based on state and time (consider
    contact forces, for example). A disabled force element will return all 
    zeroes without invoking calcForce(), since that method may depend on
    earlier computations which may not have been performed in that case.
    @param[in]      state
        The State containing information to be used by the force element to
        calculate the current force. This must have already been realized to
        a high enough stage for the force element to get what it needs; if you
        don't know then realize it to Stage::Velocity.
    @param[out]     bodyForces
        This is a Vector of spatial forces, one per mobilized body in the 
        matter subsystem associated with this force element. This Vector is
        indexed by MobilizedBodyIndex so it has a 0th entry corresponding
        to Ground. A spatial force contains two Vec3's; index with [0] to get
        the moment vector, with [1] to get the force vector. This argument is
        resized if necessary to match the number of mobilized bodies and any
        unused entry will be set to zero on return.
    @param[out]     particleForces
        This is a Vector of force vectors, one per particle in the 
        matter subsystem associated with this force element. This vector is
        indexed by ParticleIndex; the 0th entry is the 1st particle, not Ground.
        This argument is resized if necessary to match the number of particles
        and any unused entry will be set to zero on return. (As of March 2010 
        Simbody treats particles as mobilized bodies so this is unused.)
    @param[out]     mobilityForces
        This is a Vector of scalar generalized forces, one per mobility in 
        the matter subsystem associated with this force element. This is the
        same as the number of generalized speeds u that collectively represent
        all the mobilities of the mobilizers. To determine the per-mobilizer
        correspondence, you must call methods of MobilizedBody; there is no
        hint here. 
    @note This method must zero out the passed in arrays, and in most cases
    almost all returned entries will be zero, so this is \e not the most
    efficent way to calculate forces; use it sparingly. **/
    void calcForceContribution(const State&          state,
                               Vector_<SpatialVec>&  bodyForces,
                               Vector_<Vec3>&        particleForces,
                               Vector&               mobilityForces) const;
    /** Calculate the potential energy contribution that is made by this
    force element at the given \a state. This calls the force element's
    calcPotentialEnergy() method. A disabled force element will return zero 
    without invoking calcPotentialEnergy().
    @param[in]      state
        The State containing information to be used by the force element to
        calculate the current potential energy. This must have already been 
        realized to a high enough stage for the force element to get what it 
        needs; if you don't know then realize it to Stage::Position.
    @return The potential energy contribution of this force element at this
    \a state value. **/
    Real calcPotentialEnergyContribution(const State& state) const;
    /*@}*/

    /**@name                   Bookkeeping
    These methods are not normally needed. They provide bookkeeping 
    information such as access to the parent force subsystem and the force
    index assigned to this force element. **/
    /*@{*/
    /** Default constructor for Force handle base class does nothing. **/
    Force() {}
    /** Implicit conversion to ForceIndex when needed. This will throw an 
    exception if the force element has not yet been adopted by a force 
    subsystem. **/
    operator ForceIndex() const {return getForceIndex();}
    /** Get the GeneralForceSubsystem of which this Force is an element. 
    This will throw an exception if the force element has not yet been
    adopted by a force subsystem. **/
    const GeneralForceSubsystem& getForceSubsystem() const;
    /** Get the index of this force element within its parent force subsystem.
    The returned index will be invalid if the force element has not yet been
    adopted by any subsystem (test with the index.isValid() method). **/
    ForceIndex getForceIndex() const;
    /*@}*/
    
    class TwoPointLinearSpring;
    class TwoPointLinearDamper;
    class TwoPointConstantForce;
    class MobilityLinearSpring;
    class MobilityLinearDamper;
    class MobilityConstantForce;
    class MobilityLinearStop;
    class MobilityDiscreteForce;
    class DiscreteForces;
    class LinearBushing;
    class ConstantForce;
    class ConstantTorque;
    class GlobalDamper;
    class Thermostat;
    class UniformGravity;
    class Gravity;
    class Custom;
    
    class TwoPointLinearSpringImpl;
    class TwoPointLinearDamperImpl;
    class TwoPointConstantForceImpl;
    class MobilityLinearSpringImpl;
    class MobilityLinearDamperImpl;
    class MobilityConstantForceImpl;
    class MobilityLinearStopImpl;
    class MobilityDiscreteForceImpl;
    class DiscreteForcesImpl;
    class LinearBushingImpl;
    class ConstantForceImpl;
    class ConstantTorqueImpl;
    class GlobalDamperImpl;
    class ThermostatImpl;
    class UniformGravityImpl;
    class GravityImpl;
    class CustomImpl;

protected:
    /** Use this in a derived Force handle class constructor to supply the 
    concrete implementation object to be stored in the handle base. **/
    explicit Force(ForceImpl* r) : HandleBase(r) { }
};


/**
 * A linear spring between two points, specified as a station on
 * each of two bodies. The stiffness k and 
 * unstretched length x0 are provided. Then if d is the unit vector
 * from point1 to point2, and x the current separation, we have
 * f = k(x-x0) and we apply a force f*d to point1 and -f*d to point2.
 * This contributes to potential energy: pe = 1/2 k (x-x0)^2.
 * It is an error if the two points become coincident, since we
 * are unable to determine a direction for the force in that case.
 */

class SimTK_SIMBODY_EXPORT Force::TwoPointLinearSpring : public Force {
public:
    /**
     * Create a TwoPointLinearSpring.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body1      the first body to which the force should be applied
     * @param station1   the location on the first body at which the force should be applied
     * @param body2      the second body to which the force should be applied
     * @param station2   the location on the second body at which the force should be applied
     * @param k          the spring constant
     * @param x0         the distance at which the force is 0
     */
    TwoPointLinearSpring(GeneralForceSubsystem& forces, const MobilizedBody& body1, const Vec3& station1, const MobilizedBody& body2, const Vec3& station2, Real k, Real x0);
    
    /** Default constructor creates an empty handle. **/
    TwoPointLinearSpring() {}

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(TwoPointLinearSpring, TwoPointLinearSpringImpl, Force);
};

/**
 * A force which resists changes in the distance between
 * two points, acting along the line between those points. The points
 * are specified as a station on each of two bodies. A damping constant
 * c >= 0 is given. If the relative (scalar) velocity between the points
 * is v, then we apply a force of magnitude f=c*|v| to each point in
 * a direction which opposes their separation. This is not a potential
 * force and thus does not contribute to the potential energy calculation.
 * It is an error if the two points become coincident, since we
 * are unable to determine a direction for the force in that case.
 */

class SimTK_SIMBODY_EXPORT Force::TwoPointLinearDamper: public Force {
public:
    /**
     * Create a TwoPointLinearDamper.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body1      the first body to which the force should be applied
     * @param station1   the location on the first body at which the force should be applied
     * @param body2      the second body to which the force should be applied
     * @param station2   the location on the second body at which the force should be applied
     * @param damping    the damping constant
     */
    TwoPointLinearDamper(GeneralForceSubsystem& forces, const MobilizedBody& body1, const Vec3& station1, const MobilizedBody& body2, const Vec3& station2, Real damping);
    
    /** Default constructor creates an empty handle. **/
    TwoPointLinearDamper() {}

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(TwoPointLinearDamper, TwoPointLinearDamperImpl, Force);
};

/**
 * A constant force f (a signed scalar) which acts along the line between
 * two points, specified as a station on each of two bodies. A positive
 * force acts to separate the points; negative pulls them together. The
 * force magnitude is independent of the separation between the points.
 * This force does not contribute to the potential energy, so adding it to
 * a system will cause energy not to be conserved.
 * It is an error if the two points become coincident, since we
 * are unable to determine a direction for the force in that case.
 */

class SimTK_SIMBODY_EXPORT Force::TwoPointConstantForce: public Force {
public:
    /**
     * Create a TwoPointConstantForce.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body1      the first body to which the force should be applied
     * @param station1   the location on the first body at which the force should be applied
     * @param body2      the second body to which the force should be applied
     * @param station2   the location on the second body at which the force should be applied
     * @param force      the magnitude of the force to apply
     */
    TwoPointConstantForce(GeneralForceSubsystem& forces, const MobilizedBody& body1, const Vec3& station1, const MobilizedBody& body2, const Vec3& station2, Real force);
    
    /** Default constructor creates an empty handle. **/
    TwoPointConstantForce() {}

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(TwoPointConstantForce, TwoPointConstantForceImpl, Force);
};


/**
 * A constant force applied to a body station. The force is a
 * vector fixed forever in the Ground frame.
 * This force does not contribute to the potential energy, so adding it to
 * a system will cause energy not to be conserved.
 */

class SimTK_SIMBODY_EXPORT Force::ConstantForce: public Force {
public:
    ConstantForce(GeneralForceSubsystem& forces, const MobilizedBody& body, const Vec3& station, const Vec3& force);
    
    /** Default constructor creates an empty handle. **/
    ConstantForce() {}

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(ConstantForce, ConstantForceImpl, Force);
};

/**
 * A constant torque to a body. The torque is a
 * vector fixed forever in the Ground frame.
 * This force does not contribute to the potential energy, so adding it to
 * a system will cause energy not to be conserved.
 */

class SimTK_SIMBODY_EXPORT Force::ConstantTorque: public Force {
public:
    ConstantTorque(GeneralForceSubsystem& forces, const MobilizedBody& body, const Vec3& torque);
    
    /** Default constructor creates an empty handle. **/
    ConstantTorque() {}

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(ConstantTorque, ConstantTorqueImpl, Force);
};

/**
 * A general energy "drain" on the system. This is 
 * done by effectively adding a damper to every generalized
 * speed (mobility) in the system. Each generalized speed 
 * u_i feels a force -dampingFactor*u_i.
 * This usually is not physically meaningful, but it can be useful in some
 * circumstances just to drain energy out of the model when
 * the specific energy-draining mechanism is not important.
 * You can have more than one of these in which case the
 * dampingFactors are added. No individual dampingFactor is
 * allowed to be negative. This is not a potential force and
 * hence does not contribute to potential energy.
 */

class SimTK_SIMBODY_EXPORT Force::GlobalDamper : public Force {
public:
    GlobalDamper(GeneralForceSubsystem& forces, const SimbodyMatterSubsystem& matter, Real damping);
    
    /** Default constructor creates an empty handle. **/
    GlobalDamper() {}

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(GlobalDamper, GlobalDamperImpl, Force);
};

/**
 * A uniform gravitational force applied to every body in the system.  The force
 * is specified by a vector in the Ground frame.  You can optionally specify
 * a height at which the gravitational potential energy is zero.
 */

class SimTK_SIMBODY_EXPORT Force::UniformGravity : public Force {
public:
    UniformGravity(GeneralForceSubsystem& forces, const SimbodyMatterSubsystem& matter, const Vec3& g, Real zeroHeight=0);
    
    /** Default constructor creates an empty handle. **/
    UniformGravity() {}

    Vec3 getGravity() const;
    void setGravity(const Vec3& g);
    Real getZeroHeight() const;
    void setZeroHeight(Real height);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(UniformGravity, UniformGravityImpl, Force);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_H_
