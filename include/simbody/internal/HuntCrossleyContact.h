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
 * a continuous model based on Hertz elastic contact theory,
 * which correctly reproduces the empirically observed dependence
 * on velocity of coefficient of restitution, where e=(1-cv) for 
 * (small) impact velocity v and a material property c with units 1/v. Note that
 * c can be measured right off the coefficient of restitution-vs.-velocity
 * curves: it is the absolute value of the slope at low velocities.
 *
 * Given a collision between two spheres, or a sphere and a plane,
 * we can generate a contact force from this equation
 *      f = kx^n(1 + 3/2 cv)
 * where k is a stiffness constant incorporating material properties
 * and geometry (to be defined below), x is penetration depth and 
 * v = dx/dt is penetration rate (positive during penetration and
 * negative during rebound). Exponent n depends on the surface
 * geometry. For Hertz contact where the geometry can be approximated
 * by sphere (or sphere-plane) interactions, which is all we are
 * currently handling here, n=3/2.
 *
 * Stiffness k is defined in terms of the relative radius of curvature R and
 * effective plane-strain modulus E, each of which is a combination of
 * the description of the two individual contacting elements. TODO: 
 * derivation of the following results should be in the SimTK Engr J; you'll
 * have to take my word for it now:
 *
 *          R1*R2                                         E2^(2/3)
 *     R = -------,  E = (s1 * E1^(2/3))^(3/2),  s1= ------------------- 
 *         R1 + R2                                   E1^(2/3) + E2^(2/3)
 *     
 *     c = c1*s1 + c2*(1-s1)
 *     k = (4/3) sqrt(R) E
 *     f = k x^(3/2) (1 + 3/2 c xdot)
 *     pe = 2/5 k x^(5/2)
 * Also, we can calculate the contact patch radius a as
 *     a = sqrt(R*x)
 *
 * In the above, E1 and E2 are the *plane strain* moduli. If you have instead
 * Young's modulus Y1 and Poisson's ratio p1, then E1=Y1/(1-p1^2). The interface
 * to this subsystem asks for E1 (pressure/%strain) and c1 (1/velocity), and
 * E2,c2 only.
 *
 */
class SimTK_SIMBODY_EXPORT HuntCrossleyContact : public ForceSubsystem {
public:
    HuntCrossleyContact();

    int addSphere(int body, const Vec3& center,
                  const Real& radius,
                  const Real& stiffness,     // E (plane strain)
                  const Real& dissipation);  // c (1/v)

    int addHalfSpace(int body, const UnitVec3& normal,
                     const Real& height,
                     const Real& stiffness,     // E (plane strain)
                     const Real& dissipation);  // c (1/v)

    SimTK_PIMPL_DOWNCAST(HuntCrossleyContact, ForceSubsystem);
private:
    class HuntCrossleyContactRep& updRep();
    const HuntCrossleyContactRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_HUNT_CROSSLEY_CONTACT_H_
