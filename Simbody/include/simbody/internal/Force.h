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
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
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
    class MobilityDiscreteForce;
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
    class MobilityDiscreteForceImpl;
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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(TwoPointConstantForce, TwoPointConstantForceImpl, Force);
};

/**
 * A linear spring along or around a mobility coordinate. The
 * stiffness k is provided, along with an arbitrary "zero" 
 * coordinate value q0 at which the spring generates no force.
 * The generated force is k*(q-q0), and potential energy is 
 * pe = 1/2 k (q-q0)^2.
 * This is not meaningful unless the mobility coordinate is such that qdot=u 
 * for that coordinate.  In particular, do not use this on a coordinate
 * which is part of a quaternion.
 */

class SimTK_SIMBODY_EXPORT Force::MobilityLinearSpring : public Force {
public:
    /**
     * Create a %MobilityLinearSpring.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body       the body to which the force should be applied
     * @param coordinate the index of the coordinate in the body's u vector to which the force should be applied
     * @param k          the spring constant
     * @param q0         the value of the coordinate at which the force is 0
     */
    MobilityLinearSpring(GeneralForceSubsystem& forces, const MobilizedBody& body, int coordinate, Real k, Real q0);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(MobilityLinearSpring, MobilityLinearSpringImpl, Force);
};

/**
 * A linear damper on a mobility coordinate. The
 * damping constant c is provided, with the generated force
 * being -c*u where u is the mobility's generalize speed.
 * This is meaningful on any mobility, since all our
 * generalized speeds have physical meaning. This is not
 * a potential force and hence does not contribute to
 * potential energy.
 */

class SimTK_SIMBODY_EXPORT Force::MobilityLinearDamper : public Force {
public:
    /**
     * Create a %MobilityLinearDamper.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body       the body to which the force should be applied
     * @param coordinate the index of the coordinate in the body's u vector to which the force should be applied
     * @param damping    the damping constant
     */
    MobilityLinearDamper(GeneralForceSubsystem& forces, const MobilizedBody& body, int coordinate, Real damping);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(MobilityLinearDamper, MobilityLinearDamperImpl, Force);
};

/**
 * A constant (scalar) "force" f applied to a mobility. The mobility here selects a
 * generalized speed (u), not a generalized coordinate (q), and the
 * meaning depends on the definition of the generalized speed. If that
 * speed is a translation then this is a force; if a rotation then
 * this is a torque.
 * This force does not contribute to the potential energy, so adding it to
 * a system will cause energy not to be conserved.
 */

class SimTK_SIMBODY_EXPORT Force::MobilityConstantForce : public Force {
public:
    /**
     * Create a %MobilityConstantForce.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body       the body to which the force should be applied
     * @param coordinate the index of the coordinate in the body's u vector to which the force should be applied
     * @param force      the force to apply
     */
    MobilityConstantForce(GeneralForceSubsystem& forces, const MobilizedBody& body, int coordinate, Real force);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(MobilityConstantForce, MobilityConstantForceImpl, Force);
};


/** A discrete generalized force f applied to a mobility.\ Useful for applying
external forces or forces that are updated at discrete times due to the
occurrence of events. Note that a mobility is a generalized speed (u), not 
a generalized coordinate (q). The meaning of a generalized force depends on the
definition of the generalized speed. If that speed is a translation then this is
a force; if a rotation then this is a torque; if something else then f has a 
comparable definition (the defining condition is that f*u should always have 
physically meaningful units of power). This force does not contribute to the 
potential energy, so adding it to a system will cause potential+kinetic energy 
not to be conserved. **/
class SimTK_SIMBODY_EXPORT Force::MobilityDiscreteForce : public Force {
public:
    /** Create a %MobilityDiscreteForce.
    
    @param forces        subsystem to which this force element should be added
    @param mobod         mobilizer to which the force should be applied
    @param whichU        to which of the mobilizer's mobilities (degrees of
                            freedom) should this force be applied (first is 0)?
    @param defaultForce  initial value for the force to be applied 
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

    /** Provide a new value for the \a defaultForce, overriding the one 
    provided in the constructor. This is a topological change because it
    affects the value that the containing System's default state will have
    when realizeTopology() is called. This is for use during construction, not
    for during a simulation where you should be using setGeneralizedForce().

    @param defaultForce     the value this generalized force should have by 
                                default
    @returns a writable reference to this modified force element 
    @see setGeneralizedForce(), getDefaultGeneralizedForce()  **/
    MobilityDiscreteForce& setDefaultGeneralizedForce(Real defaultForce);

    /** Return the value that this generalized force will have by default. This
    is normally set in the constructor, or left to its default value of 0. It 
    can also be set in setDefaultGeneralizedForce(). Note that this is \e not
    the same as the value that may be set in any particular State. 
    @see getGeneralizedForce(), setDefaultGeneralizedForce() **/
    Real getDefaultGeneralizedForce() const;

    /** Change the value of the generalized force to be applied in the given
    \a state. Set this to zero if you don't want it to do anything. 
    @see getGeneralizedForce() **/
    void setGeneralizedForce(State& state, Real f) const;
    /** Return the value for this generalized force that is stored in the 
    given \a state. If no calls to setGeneralizedForce() have been made on this
    \a state then it will have the \a defaultForce value that was supplied on
    construction or via setDefaultGeneralizedForce(). 
    @see setGeneralizedForce() **/
    Real getGeneralizedForce(const State& state) const;

    /** @cond **/
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(MobilityDiscreteForce, 
                                             MobilityDiscreteForceImpl, Force);
    /** @endcond **/
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
    Vec3 getGravity() const;
    void setGravity(const Vec3& g);
    Real getZeroHeight() const;
    void setZeroHeight(Real height);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(UniformGravity, UniformGravityImpl, Force);
};

/**
 * This class is used to define new force elements. To use it, you will 
 * create a class that extends Force::Custom::Implementation. Optionally, you
 * may also create a "handle" class extending Force::Custom that hides your 
 * implementation and provides a nicer API for users of your new force element.
 * Here we'll use as an example a force "MySpring" that we'll presume you will 
 * declare in a header file "MySpring.h". The MySpring constructor will take 
 * a reference to a MobilizedBody \a mobod and a spring constant \a k, and 
 * connect the origin OB of body \a mobod to the Ground origin O by a spring of
 * stiffness \a k. MySpring will apply a force of magnitude \a kx to OB, 
 * directed towards O, where x=|OB-O| is the current distance from OB to O. 
 * First we'll look at how MySpring would be used in a main program (whether by
 * you or some future user of your force element), then how you would 
 * implement it. 
 *
 * There are two possibilities: 
 *  - you supply only an %Implementation object MySpringImpl (less work 
 *    for you), or 
 *  - you supply both the %Implementation object and a handle MySpring
 *    (nicer for the user). 
 *
 * If you are planning only to use this force yourself, and perhaps just once
 * for a single application, the %Implementation-only approach is probably 
 * adequate. However, if you plan to have others use your new force object, it
 * is well worth the (minimal) extra effort to provide a handle class also to 
 * make your new force behave identically to the built-in forces that come 
 * with Simbody. In either case the %Implementation object is the same so you 
 * can add a handle later if you want. To reiterate: the handle class is 
 * completely optional; you \e must write an %Implementation class
 * but a handle class is an aesthetic addition whose main purpose is to
 * make a cleaner API for users of your force element.
 *
 * In the case where you write only the %Implementation class, a user will 
 * create an instance of that class and pass it to the generic Force::Custom 
 * constructor (this would be in main.cpp or some other piece of application
 * code):
 * @code
 *  // Using your custom force in a simulation application, no handle.
 *  #include "SimTKsimbody.h"
 *  #include "MySpring.h"
 *  using namespace SimTK;
 *  GeneralForceSubsystem forces;
 *  MobilizedBody mobod1 = MobilizedBody::Pin(...); // or whatever
 *  // ...
 *  Force::Custom spring1(forces, new MySpringImpl(mobod1,100.));
 * @endcode
 * If you also write a handle class, then the API seen by the user is the
 * same as for built-in force objects:
 * @code
 *  // Using your custom force in a simulation application, with handle.
 *  // Same as above except last line is now:
 *  MySpring spring1(forces, mobod1, 100.);
 * @endcode
 *
 * <h2>Writing the %Implementation class</h2>
 *
 * You would put the following code in MySpring.h, although you can hide it
 * elsewhere (such as in a separate .cpp file) if you provide a handle class:
 * @code
 *  // Sample implementation of a custom force Implementation class.
 *  class MySpringImpl: public Force::Custom::Implementation {
 *  public:
 *      MySpringImpl(MobilizedBody mobod, Real k)
 *      :   m_mobod(mobod), m_k(k) {}
 *
 *      virtual void calcForce(const State&         state, 
 *                             Vector_<SpatialVec>& bodyForcesInG, 
 *                             Vector_<Vec3>&       particleForcesInG, 
 *                             Vector&              mobilityForces) const 
 *      {
 *          Vec3 bodyPointInB(0,0,0); // body origin (in the body frame)
 *          Vec3 pos = m_mobod.getBodyOriginLocation(state); // in Ground frame
 *          // apply force towards origin, proportional to distance x
 *          m_mobod.applyForceToBodyPoint(state, bodyPointInB, -m_k*pos,
 *                                        bodyForcesInG);
 *      }
 *
 *      virtual Real calcPotentialEnergy(const State& state) const {
 *          Vec3 pos = m_mobod.getBodyOriginLocation(state); // in Ground
 *          Real x   = pos.norm(); // distance from Ground origin
 *          return m_k*x*x/2; // potential energy in the spring
 *      }
 *  private:
 *      MobilizedBody m_mobod;    // the body to which to apply the force
 *      Real          m_k;        // spring stiffness
 *  };
 * @endcode
 *
 * To write the code exactly as above, the compiler has to be told to look in
 * the %SimTK namespace for some of the symbols. There are a variety of ways
 * to do that; see the discussion below for details.
 * 
 * <h2>Writing the handle class</h2>
 *
 * Here is how to implement the handle class, assuming you've already 
 * written the %Implementation class MySpringImpl. You would put this 
 * code (at least the declarations) in MySpring.h, following the declaration
 * of the MySpringImpl class:
 * @code
 *  // Sample implementation of a handle class. Note: a handle class
 *  // may *not* have any data members or virtual methods.
 *  class MySpring : public Force::Custom {
 *  public:
 *      MySpring(GeneralForceSubsystem& forces, // note the "&"
 *               MobilizedBody          mobod,
 *               Real                   k) 
 *      :   Force::Custom(forces, new MySpringImpl(mobod,k)) {}
 *  };
 * @endcode
 * As you can see, the handle class is very simple and just hides the creation
 * of the implementation object. Since this removes any reference to the 
 * implementation object from the user's program, it also means you can hide
 * the implementation details completely (perhaps in a separately compiled
 * library), which has many advantages. You can add additional methods to the 
 * handle class to provide a clean API to users of your custom force; these 
 * methods will forward to the implementation object as necessary but will not 
 * expose any aspects of the implementation that are not needed by the user.
 *
 * @warning A handle class <em>must not</em> have any data members or virtual
 * methods and will not work properly if it does. Instead, store any data you
 * need in the %Implementation object.
 *
 * <h2>Getting access to symbols in the %SimTK namespace</h2>
 *
 * The examples above glossed over the naming of symbols in the SimTK namespace.
 * To write the code exactly as above you would need to precede it with
 * \c using statements to tell the compiler to look there to resolve symbols, 
 * for example:
 * @code
 *  using namespace SimTK; // search the entire namespace for symbols, or
 *  using SimTK::Real; using SimTK::Vector_; // define just particular symbols
 * @endcode
 * Either of those approaches will work, but note that this will have the same 
 * effect on user programs that include MySpring.h as it does within the header 
 * file. That may be unwanted behavior for some users who might prefer not to have
 * the namespace searched automatically, perhaps to avoid conflicts with their own
 * symbols that have the same names. If you want to avoid introducing unwanted 
 * symbols into users' compilation units, then in the header file you should 
 * refer to each symbol explicitly by its full name; that is, prepend the 
 * namespace each time the symbol is used, for example:
 * @code
 *  class MySpringImpl: public SimTK::Force::Custom::Implementation {
 *  void calcForce(const SimTK::State&                state, 
 *                 SimTK::Vector_<SimTK::SpatialVec>& bodyForcesInG, 
 *                 SimTK::Vector_<SimTK::Vec3>&       particleForcesInG, 
 *                 SimTK::Vector&                     mobilityForces) const ;
 *  };
 * @endcode
 * That is less convenient for you (and uglier) but avoids the possibility of
 * unwanted side effects in user code so should be considered if you expect
 * wide distribution of your new force element.
 *
 * Thanks to Nur Adila Faruk Senan (a.k.a. Adila Papaya) for help in clarifying
 * this documentation.
 */
class SimTK_SIMBODY_EXPORT Force::Custom : public Force {
public:
    class Implementation;
    /**
     * Create a Custom force.
     * 
     * @param forces         the subsystem to which this force should be added
     * @param implementation the object which implements the custom force.  The Force::Custom takes over
     *                       ownership of the implementation object, and deletes it when the Force itself
     *                       is deleted.
     */
    Custom(GeneralForceSubsystem& forces, Implementation* implementation);

    /** @cond **/ // don't let Doxygen see this
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Custom, CustomImpl, Force);
    /** @endcond **/
protected:
    const Implementation& getImplementation() const;
    Implementation& updImplementation();
};

/**
 * Every custom force requires implementation of a class that is derived from
 * this abstract class.\ See Force::Custom for details.
 */
class SimTK_SIMBODY_EXPORT Force::Custom::Implementation {
public:
    virtual ~Implementation() { }
    /**
     * Calculate the force for a given state.
     * 
     * @param state          the State for which to calculate the force
     * @param bodyForces     spatial forces on MobilizedBodies are accumulated in this.  To apply a force to a body,
     *                       add it to the appropriate element of this vector.
     * @param particleForces forces on particles are accumulated in this.  Since particles are not yet implemented,
     *                       this is ignored.
     * @param mobilityForces forces on individual mobilities (elements of the state's u vector) are accumulated in this.
     *                       To apply a force to a mobility, add it to the appropriate element of this vector.
     */
    virtual void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces) const = 0;
    /**
     * Calculate this force's contribution to the potential energy of the System.
     * 
     * @param state          the State for which to calculate the potential energy
     */
    virtual Real calcPotentialEnergy(const State& state) const = 0;
    /**
     * Get whether this force depends only on the position variables (q), not on the velocies (u) or auxiliary variables (z).
     * The default implementation returns false.  If the force depends only on positions, you should override this to return
     * true.  This allows force calculations to be optimized in some cases.
     */
    virtual bool dependsOnlyOnPositions() const {
        return false;
    }
    /** The following methods may optionally be overridden to do specialized 
    realization for a Force. **/
    //@{
    virtual void realizeTopology(State& state) const {}
    virtual void realizeModel(State& state) const {}
    virtual void realizeInstance(const State& state) const {}
    virtual void realizeTime(const State& state) const {}
    virtual void realizePosition(const State& state) const {}
    virtual void realizeVelocity(const State& state) const {}
    virtual void realizeDynamics(const State& state) const {}
    virtual void realizeAcceleration(const State& state) const {}
    virtual void realizeReport(const State& state) const {}
    //@}

    /** Override this if you want to generate some geometry for the visualizer
    to display. Be sure to \e append the geometry to the list. **/
    virtual void calcDecorativeGeometryAndAppend
       (const State& state, Stage stage, 
        Array_<DecorativeGeometry>& geometry) const {}

};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_H_
