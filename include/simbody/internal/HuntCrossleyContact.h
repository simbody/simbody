#ifndef SimTK_HUNT_CROSSLEY_CONTACT_H_
#define SimTK_HUNT_CROSSLEY_CONTACT_H_

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
 * Define the public interface to HuntCrossleyContact, a subsystem which
 * provides some minimal contact behavior.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/State.h"
#include "simbody/internal/System.h"
#include "simbody/internal/ForceSubsystem.h"

#include <cassert>

namespace SimTK {

/**
 * This is a concrete subsystem that handles simple, frictionless contact situations
 * with a model due to Hunt & Crossley: K. H. Hunt and F. R. E. Crossley, 
 * "Coefficient of Restitution Interpreted as Damping in Vibroimpact,"
 * ASME Journal of Applied Mechanics, pp. 440-445, June 1975. This is
 * a continuous model which correctly reproduces the empirically observed
 * dependence on velocity of coefficient of restitution, where e=(1-cv) for 
 * impact velocity v and a material property c with units 1/v. Note that
 * c can be measured right off the coefficient of restitution-vs.-velocity
 * curves: it is the absolute value of the slope at low velocities.
 *
 * Given a collision between two spheres, or a sphere and a plane,
 * we can generate a contact force from this equation
 *      f = -kx^n(1 + 3/2 cv)
 * where k is the (linear) material stiffness, x is penetration depth and 
 * v = dx/dt is penetration rate (positive during penetration and
 * negative during rebound). Exponent n depends on the surface
 * geometry. For Hertz contact where the geometry can be approximated
 * by sphere (or sphere-plane) interactions, n=3/2.
 *
 * When dissimilar materials collide, we use the following combining rules:
 *    k = k1k2/(k1+k2)
 *    c = (c1k2+c2k1)/(k1+k2)
 *
 */
class SimTK_SIMBODY_API HuntCrossleyContact : public ForceSubsystem {
public:
    HuntCrossleyContact();

    int addSphere(int body, const Vec3& center,
                  const Real& radius,
                  const Real& stiffness,     // f per unit distance
                  const Real& dissipation);  // 1/v

    int addHalfSpace(int body, const UnitVec3& normal,
                     const Real& height,
                     const Real& stiffness,     // f per unit distance
                     const Real& dissipation);  // 1/v

    SimTK_PIMPL_DOWNCAST(HuntCrossleyContact, ForceSubsystem);
private:
    class HuntCrossleyContactRep& updRep();
    const HuntCrossleyContactRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_HUNT_CROSSLEY_CONTACT_H_
