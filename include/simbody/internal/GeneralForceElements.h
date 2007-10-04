#ifndef SimTK_SIMBODY_GENERAL_FORCE_ELEMENTS_H_
#define SimTK_SIMBODY_GENERAL_FORCE_ELEMENTS_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
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

/** @file
 * Define the public interface to the SimbodyForcesSubsystem, a collection
 * of basic force elements suitable for use with a SimbodyMatterSubsystem.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/ForceSubsystem.h"

#include <cassert>

namespace SimTK {

class MultibodySystem;
class SimbodyMatterSubsystem;

/**
 * This is a concrete subsystem which can apply a variety of
 * simple force elements to the MatterSubsystem within a MultibodySystem.
 */
class SimTK_SIMBODY_EXPORT GeneralForceElements : public ForceSubsystem {
public:

    // See below for declaration of abstract class GeneralForceElements::CustomForce
    // which can be used to write your own force elements.
    class CustomForce;


public:
    GeneralForceElements();
    explicit GeneralForceElements(MultibodySystem&);

    /// Add a linear spring between two points, specified as a station on
    /// each of two bodies. The stiffness k and 
    /// unstretched length x0 are provided. Then if d is the unit vector
    /// from point1 to point2, and x the current separation, we have
    /// f = k(x-x0) and we apply a force f*d to point1 and -f*d to point2.
    /// This contributes to potential energy: pe = 1/2 k (x-x0)^2.
    /// It is an error if the two points become coincident, since we
    /// are unable to determine a direction for the force in that case.
    int addTwoPointLinearSpring(MobilizedBodyId body1, const Vec3& s1, // point 1
                                MobilizedBodyId body2, const Vec3& s2, // point 2
                                const Real& k, const Real& x0);

    /// Add a constant force f (a signed scalar) which acts along the line between
    /// two points, specified as a station on each of two bodies. A positive
    /// force acts to separate the points; negative pulls them together. The
    /// force magnitude is independent of the separation between the points,
    /// however the potential energy does depend on the separation. We allow
    /// an arbitrary "zero" distance x0 to be supplied at which the potential
    /// energy is considered zero; the default is that potential energy is zero
    /// at zero separation. Potential energy is then calculated as pe = -f*(x-x0).
    /// It is an error if the two points become coincident, since we
    /// are unable to determine a direction for the force in that case.
    int addTwoPointConstantForce(MobilizedBodyId body1, const Vec3& s1,
                                 MobilizedBodyId body2, const Vec3& s2,
                                 const Real& f, const Real& x0=0);

    /// Add a force which resists changes in the distance between
    /// two points, acting along the line between those points. The points
    /// are specified as a station on each of two bodies. A damping constant
    /// c >= 0 is given. If the relative (scalar) velocity between the points
    /// is v, then we apply a force of magnitude f=c*|v| to each point in
    /// a direction which opposes their separation. This is not a potential
    /// force and thus does not contribute to the potential energy calculation.
    /// It is an error if the two points become coincident, since we
    /// are unable to determine a direction for the force in that case.
    int addTwoPointLinearDamper(MobilizedBodyId body1, const Vec3& s1,
                                MobilizedBodyId body2, const Vec3& s2,
                                const Real& damping);

    /// Apply a constant force to a body station. The force is a
    /// vector fixed forever in the Ground frame. This force does not
    /// contribute to potential energy.
    int addConstantForce(MobilizedBodyId body, const Vec3& s_B, const Vec3& f_G);

    /// Apply a constant torque to a body. The torque is a
    /// vector fixed forever in the Ground frame. This torque does not
    /// contribute to potential energy.
    int addConstantTorque(MobilizedBodyId body, const Vec3& t_G);

    /// Apply a constant (scalar) "force" f to a mobility. This does not
    /// contribute to the potential energy. The axis here selects a
    /// generalized speed (u), not a generalized coordinate (q), and the
    /// meaning depends on the definition of the generalized speed. If that
    /// speed is a translation then this is a force; if a rotation then
    /// this is a torque.
    int addMobilityConstantForce(MobilizedBodyId body, int axis, const Real& f);

    /// Add a linear spring along or around a mobility coordinate. The
    /// stiffness k is provided, along with an arbitrary "zero" 
    /// coordinate value q0 at which the spring generates no force.
    /// The generated force is k*(q-q0), and potential energy is 
    /// pe = 1/2 k (q-q0)^2.
    /// Not meaningful unless the mobility coordinate is such that qdot=u 
    /// for that coordinate, in particular don't use this on a coordinate
    /// which is part of a quaternion.
    int addMobilityLinearSpring(MobilizedBodyId body, int axis,
                                const Real& stiffness,
                                const Real& neutralValue);

    /// Add a linear damper to a mobility coordinate. The
    /// damping constant c is provided, with the generated force
    /// being -c*u where u is the mobility's generalize speed.
    /// This is meaningful on any mobility, since all our
    /// generalized speeds have physical meaning. This is not
    /// a potential force and hence does not contribute to
    /// potential energy.
    int addMobilityLinearDamper(MobilizedBodyId body, int axis,
                                const Real& dampingFactor);

    /// Add a general energy "drain" to the system. This is 
    /// done by effectively adding a damper to every generalized
    /// speed (mobility) in the system. Each generalized speed 
    /// u_i feels a force -dampingFactor*u_i.
    /// This is not physically meaningful but can be useful in some
    /// circumstances just to drain energy out of the model when
    /// the specific energy-draining mechanism is not important.
    /// You can have more than one of these in which case the
    /// dampingFactors are added. No individual dampingFactor is
    /// allowed to be negative. This is not a potential force and
    /// hence does not contribute to potential energy.
    int addGlobalEnergyDrain(const Real& dampingFactor);

    /// Add a UserForce to the system. Note that a private copy of the UserForce
    /// is kept internally; the original is only used as a prototype and subsequent
    /// changes to it will have no effect.

    // (This routine must execute on the client side; it's implementation
    // is provided later in this header file when we have seen all the
    // needed declarations.) 
    inline int addCustomForce(const CustomForce& u);

    // No user-serviceable parts below here.

    // See the UserForce class declaration below for more information.
    typedef void (*CustomForceCalcMethod)
       (const CustomForce&             u,
        const SimbodyMatterSubsystem&  matter, 
        const State&                   state,
        Vector_<SpatialVec>&           bodyForces,
        Vector_<Vec3>&                 particleForces,
        Vector&                        mobilityForces,
        Real&                          pe);
    typedef CustomForce* (*CustomForceCloneMethod)(const CustomForce&);
    typedef void (*CustomForceDestructor)(CustomForce*);

    SimTK_PIMPL_DOWNCAST(GeneralForceElements, ForceSubsystem);
private:
    // This library-side routine is how user forces are actually
    // conveyed, ensuring binary compatibility when the client and
    // library are not at the same version.
    int addCustomForceMethods(const CustomForce& u, 
        CustomForceCalcMethod, CustomForceCloneMethod, CustomForceDestructor);

    class GeneralForceElementsRep& updRep();
    const GeneralForceElementsRep& getRep() const;
};


/**
 * This is an abstract class from which you can derive a concrete force element.
 * You must provide two methods:
 *    calc()  -- evaluates forces at a particular state
 *    clone() -- makes a copy of your UserForce object
 * If your class uses any heap-allocated private data, be sure to implement a
 * destructor also. In that case for your own sanity you *should* implement a
 * copy constructor and default assignment operator as well, although those
 * are not required for a UserForce to be used correctly by the GeneralForceElements
 * subsystem.
 *
 * A critical requirement on your UserForce, which no one but you can enforce,
 * is that it be *stateless*. That is, past history MUST NOT affect the results
 * of a call to your calc() routine with a given State. If you violate this
 * rule you may still get plausible-looking numbers and movies out of your
 * simulation, but you're not doing physics any more.
 *
 * TODO: the current implementation does not allow a CustomForce to allocate
 * any space in the State. It should be given such a chance in which case
 * it would no longer need to be stateless.
 */
class SimTK_SIMBODY_EXPORT GeneralForceElements::CustomForce {
public:
    virtual ~CustomForce() { }

    // The force arrays and potential energy are accumulated as we go,
    // so this routine must *add in* its own contributions, not replace
    // what is already there. 'bodyForces' is indexed by MobilizedBodyId, so the
    // first entry corresponds to Ground -- you can apply forces there
    // if you want but they won't do anything. Body torques and forces
    // are supplied in the Ground frame. 'mobilityForces' are
    // indexed by generalized speed number and take on the meaning
    // appropriate for that mobility. Note that mobility forces act
    // like a motor placed on that mobility -- they apply equal
    // and opposite forces or torques on the body and its parent body. 
    // Body (and particle) forces on the other hand act like the Hand
    // of God and apply just the single force or torque without regard
    // to how it would have to be implemented physically.
    // TODO: particleForces will be indexed by particle number; forces
    // will be expressed in the Ground frame.
    // If this is a potential (conservative) force, that is, a force
    // which depends only on configuration, then you should calculate
    // its current contribution to the system potential energy and add
    // it in to the indicated argument. Otherwise, just leave 'pe' alone;
    // don't set it zero!
    virtual void calc(const SimbodyMatterSubsystem& matter, const State& state,
                      Vector_<SpatialVec>& bodyForces,
                      Vector_<Vec3>&       particleForces,
                      Vector&              mobilityForces,
                      Real&                pe) const = 0;

    virtual CustomForce* clone() const = 0;

private:
    // This is tricky because no library-side code can depend on the ordering
    // of methods in the virtual function table of this abstract class. So
    // we call the virtual functions from private static methods here which
    // are generated on the client side and passed to us as though they were
    // C function addresses.

    // TODO: I'm not sure this works with static members -- they are supposed
    // to be shared across all instances of the containing class which means
    // the linker might try to use one from the library rather than here.
    // These should be made global statics just below the class definition, because
    // global statics are compilation-unit scoped and not shared with anyone.
    inline static void staticCalc(const CustomForce& u,
        const SimbodyMatterSubsystem& matter, const State& state,
        Vector_<SpatialVec>& bodyForces,
        Vector_<Vec3>&       particleForces,
        Vector&              mobilityForces,
        Real&                pe)
    {
        u.calc(matter,state,bodyForces,particleForces,mobilityForces,pe);
    }
    inline static CustomForce* staticClone(const CustomForce& u) {return u.clone();}
    inline static void staticDestructor(CustomForce* u) {delete u;}
    friend class GeneralForceElements;
};


// This routine must execute on the client side, so that the library side
// doesn't depend on the ordering of virtual function table entries.
inline int GeneralForceElements::addCustomForce(const CustomForce& u) {
    return addCustomForceMethods(u, CustomForce::staticCalc, CustomForce::staticClone,
                                    CustomForce::staticDestructor);
}

} // namespace SimTK

#endif // SimTK_SIMBODY_GENERAL_FORCE_ELEMENTS_H_
