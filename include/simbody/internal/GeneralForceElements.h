#ifndef SimTK_GENERAL_FORCE_ELEMENTS_H_
#define SimTK_GENERAL_FORCE_ELEMENTS_H_

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
 * This is a concrete subsystem which can apply a variety of
 * simple force elements to the MatterSubsystem within a MultibodySystem.
 */
class SimTK_SIMBODY_EXPORT GeneralForceElements : public ForceSubsystem {
public:

    // This is tricky because no library-side code can depend on the ordering
    // of methods in the virtual function table of this abstract class. So
    // we call the virtual functions from private static methods here which
    // are generated on the client side and passed to us as though they were
    // C function addresses.
    class UserForce {
    public:
        virtual ~UserForce() { }
        virtual void calc(const MatterSubsystem& matter, const State& state,
                          Vector_<SpatialVec>& bodyForces,
                          Vector_<Vec3>&       particleForces,
                          Vector&              mobilityForces,
                          Real&                pe) const = 0;

        virtual UserForce* clone() const = 0;

    private:
        // These routines provide a C-compatible interface to the virtual methods.
        inline static void staticCalc(const UserForce* u,
            const MatterSubsystem& matter, const State& state,
            Vector_<SpatialVec>& bodyForces,
            Vector_<Vec3>&       particleForces,
            Vector&              mobilityForces,
            Real&                pe)
        {
            u->calc(matter,state,bodyForces,particleForces,mobilityForces,pe);
        }
        inline static UserForce* staticClone(const UserForce* u) {
            return u->clone();
        }
        inline static void staticDestructor(UserForce* u) {
            delete u;
        }
        friend class GeneralForceElements;
    };
    typedef void (*UserForceCalcMethod)
       (const UserForce*        u,
        const MatterSubsystem&  matter, 
        const State&            state,
        Vector_<SpatialVec>&    bodyForces,
        Vector_<Vec3>&          particleForces,
        Vector&                 mobilityForces,
        Real&                   pe);

    typedef UserForce* (*UserForceCloneMethod)(const UserForce*);
    typedef void (*UserForceDestructor)(UserForce*);

public:
    GeneralForceElements();

    /// Add a linear spring between two points, specified as a station on
    /// each of two bodies. The stiffness k and 
    /// unstretched length x0 are provided. Then if d is the unit vector
    /// from point1 to point2, and x the current separation, we have
    /// f = k(x-x0) and we apply a force f*d to point1 and -f*d to point2.
    /// This contributes to potential energy: pe = 1/2 k (x-x0)^2.
    /// It is an error if the two points become coincident, since we
    /// are unable to determine a direction for the force in that case.
    int addTwoPointLinearSpring(int body1, const Vec3& s1, // point 1
                                int body2, const Vec3& s2, // point 2
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
    int addTwoPointConstantForce(int body1, const Vec3& s1,
                                 int body2, const Vec3& s2,
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
    int addTwoPointLinearDamper(int body1, const Vec3& s1,
                                int body2, const Vec3& s2,
                                const Real& damping);

    /// Apply a constant force to a body station. The force is a
    /// vector fixed forever in the Ground frame. This force contributes
    /// to the potential energy, with pe = |f|*(height-zeroHeight) where
    /// height=(s_G * -f/|f|) and zeroHeight is arbitrary.
    int addConstantForce(int body, const Vec3& s_B, const Vec3& f_G,
                         const Real& zeroEnergyHeight=0);

    /// Apply a constant (scalar) force f to a mobility. This force contributes
    /// to the potential energy, with pe = -f*(q-zeroEnergyValue) where
    /// q is the mobility coordinate and zeroEnergyValue is arbitrary. 
    /// Not meaningful unless the mobility coordinate is such that qdot=u 
    /// for that coordinate, in particular don't use this on a coordinate
    /// which is part of a quaternion.
    int addMobilityConstantForce(int body, int axis, const Real& f,
                                 const Real& zeroEnergyValue=0);

    /// Add a linear spring along or around a mobility coordinate. The
    /// stiffness k is provided, along with an arbitrary "zero" 
    /// coordinate value q0 at which the spring generates no force.
    /// The generated force is k*(q-q0), and potential energy is 
    /// pe = 1/2 k (q-q0)^2.
    /// Not meaningful unless the mobility coordinate is such that qdot=u 
    /// for that coordinate, in particular don't use this on a coordinate
    /// which is part of a quaternion.
    int addMobilityLinearSpring(int body, int axis,
                                const Real& stiffness,
                                const Real& neutralValue);

    /// Add a linear damper to a mobility coordinate. The
    /// damping constant c is provided, with the generated force
    /// being -c*u where u is the mobility's generalize speed.
    /// This is meaningful on any mobility, since all our
    /// generalized speeds have physical meaning. This is not
    /// a potential force and hence does not contribute to
    /// potential energy.
    int addMobilityLinearDamper(int body, int axis,
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

    // This routine execute on the client side.
    int addUserForce(UserForce* u) {
        return addUserForceMethods(u, UserForce::staticCalc, UserForce::staticClone,
            UserForce::staticDestructor);
    }

    SimTK_PIMPL_DOWNCAST(GeneralForceElements, ForceSubsystem);
private:
    int addUserForceMethods(UserForce* u, 
        UserForceCalcMethod calc, UserForceCloneMethod clone, 
        UserForceDestructor nuke);

    class GeneralForceElementsRep& updRep();
    const GeneralForceElementsRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_GENERAL_FORCE_ELEMENTS_H_
