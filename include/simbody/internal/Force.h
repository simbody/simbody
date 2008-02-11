#ifndef SimTK_SIMBODY_FORCE_H_
#define SimTK_SIMBODY_FORCE_H_
/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2008 Stanford University and the Authors.           *
 * Authors: Peter Eastman                                                     *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon.h"
#include "simbody/internal/GeneralForceSubsystem.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

namespace SimTK {

class Force;
class ForceImpl;

// We only want the template instantiation to occur once. This symbol is defined in the SimTK core
// compilation unit that defines the Force class but should not be defined any other time.
#ifndef SimTK_DEFINING_FORCE
    extern template class PIMPLHandle<Force, ForceImpl>;
#endif

/**
 * A Force object applies forces to the bodies in a system.  There are subclasses for various standard types
 * of forces, or you can create your own forces by using Custom.
 */
    
class SimTK_SIMBODY_EXPORT Force : public PIMPLHandle<Force, ForceImpl> {
public:
    Force() { }
    explicit Force(ForceImpl* r) : HandleBase(r) { }

    /// Get the index of this Force in its GeneralForceSubsystem.
    ForceIndex getForceIndex() const;
    
    class TwoPointLinearSpring;
    class TwoPointLinearDamper;
    class TwoPointConstantForce;
    class MobilityLinearSpring;
    class MobilityLinearDamper;
    class MobilityConstantForce;
    class ConstantForce;
    class ConstantTorque;
    class GlobalDamper;
    class Custom;
    
    class TwoPointLinearSpringImpl;
    class TwoPointLinearDamperImpl;
    class TwoPointConstantForceImpl;
    class MobilityLinearSpringImpl;
    class MobilityLinearDamperImpl;
    class MobilityConstantForceImpl;
    class ConstantForceImpl;
    class ConstantTorqueImpl;
    class GlobalDamperImpl;
    class CustomImpl;
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

class SimTK_SIMBODY_EXPORT Force::TwoPointLinearSpring : public PIMPLDerivedHandle<TwoPointLinearSpring, TwoPointLinearSpringImpl, Force> {
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

class SimTK_SIMBODY_EXPORT Force::TwoPointLinearDamper: public PIMPLDerivedHandle<TwoPointLinearDamper, TwoPointLinearDamperImpl, Force> {
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

class SimTK_SIMBODY_EXPORT Force::TwoPointConstantForce: public PIMPLDerivedHandle<TwoPointConstantForce, TwoPointConstantForceImpl, Force> {
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

class SimTK_SIMBODY_EXPORT Force::MobilityLinearSpring : public PIMPLDerivedHandle<MobilityLinearSpring, MobilityLinearSpringImpl, Force> {
public:
    /**
     * Create a MobilityLinearSpring.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body       the body to which the force should be applied
     * @param coordinate the index of the coordinate in the body's u vector to which the force should be applied
     * @param k          the spring constant
     * @param q0         the value of the coordinate at which the force is 0
     */
    MobilityLinearSpring(GeneralForceSubsystem& forces, const MobilizedBody& body, int coordinate, Real k, Real q0);
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

class SimTK_SIMBODY_EXPORT Force::MobilityLinearDamper : public PIMPLDerivedHandle<MobilityLinearDamper, MobilityLinearDamperImpl, Force> {
public:
    /**
     * Create a MobilityLinearDamper.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body       the body to which the force should be applied
     * @param coordinate the index of the coordinate in the body's u vector to which the force should be applied
     * @param damping    the damping constant
     */
    MobilityLinearDamper(GeneralForceSubsystem& forces, const MobilizedBody& body, int coordinate, Real damping);
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

class SimTK_SIMBODY_EXPORT Force::MobilityConstantForce : public PIMPLDerivedHandle<MobilityConstantForce, MobilityConstantForceImpl, Force> {
public:
    /**
     * Create a MobilityConstantForce.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param body       the body to which the force should be applied
     * @param coordinate the index of the coordinate in the body's u vector to which the force should be applied
     * @param force      the force to apply
     */
    MobilityConstantForce(GeneralForceSubsystem& forces, const MobilizedBody& body, int coordinate, Real force);
};

/**
 * A constant force applied to a body station. The force is a
 * vector fixed forever in the Ground frame.
 * This force does not contribute to the potential energy, so adding it to
 * a system will cause energy not to be conserved.
 */

class SimTK_SIMBODY_EXPORT Force::ConstantForce: public PIMPLDerivedHandle<ConstantForce, ConstantForceImpl, Force> {
public:
    ConstantForce(GeneralForceSubsystem& forces, const MobilizedBody& body, const Vec3& station, const Vec3& force);
};

/**
 * A constant torque to a body. The torque is a
 * vector fixed forever in the Ground frame.
 * This force does not contribute to the potential energy, so adding it to
 * a system will cause energy not to be conserved.
 */

class SimTK_SIMBODY_EXPORT Force::ConstantTorque: public PIMPLDerivedHandle<ConstantTorque, ConstantTorqueImpl, Force> {
public:
    ConstantTorque(GeneralForceSubsystem& forces, const MobilizedBody& body, const Vec3& torque);
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

class SimTK_SIMBODY_EXPORT Force::GlobalDamper : public PIMPLDerivedHandle<GlobalDamper, GlobalDamperImpl, Force> {
public:
    GlobalDamper(GeneralForceSubsystem& forces, const SimbodyMatterSubsystem& matter, Real damping);
};

/**
 * This class can be used to define new forces.  To use it, create a class that extends Force::Custom::Implementation.
 * You can then create an instance of it and pass it to the Force::Custom constructor:
 * 
 * <pre>
 * Force::Custom myForce(forces, new MyForceImplementation());
 * </pre>
 * 
 * Alternatively, you can create a subclass of Force::Custom which creates the Implementation itself:
 * 
 * <pre>
 * class MyForce : public Force::Custom {
 * public:
 *   MyForce(GeneralForceSubsystem& forces) : Force::Custom(forces, new MyForceImplementation()) {
 *   }
 * }
 * </pre>
 * 
 * This allows a user to simply write
 * 
 * <pre>
 * MyForce(forces);
 * </pre>
 * 
 * and not worry about implementation classes or creating objects on the heap.  If you do this, your Force::Custom
 * subclass must not have any data members or virtual methods.  If it does, it will not work correctly.  Instead,
 * store all data in the Implementation subclass.
 */

class SimTK_SIMBODY_EXPORT Force::Custom : public PIMPLDerivedHandle<Custom, CustomImpl, Force> {
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
protected:
    const Implementation& getImplementation() const;
    Implementation& updImplementation();
};

class SimTK_SIMBODY_EXPORT Force::Custom::Implementation {
public:
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
     * @param pe             the system's potential energy is accumulated in this.  If this force affects the potential
     *                       energy, add the energy change to this variable.
     */
    virtual void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, Vector_<Vec3>& particleForces, Vector& mobilityForces, Real& pe) const = 0;
    /**
     * Get whether this force depends only on the position variables (q), not on the velocies (u) or auxiliary variables (z).
     * The default implementation returns false.  If the force depends only on positions, you should override this to return
     * true.  This allows force calculations to be optimized in some cases.
     */
    virtual bool dependsOnlyOnPositions() const {
        return false;
    }
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_H_
