#ifndef SimTK_SIMBODY_HUNT_CROSSLEY_CONTACT_H_
#define SimTK_SIMBODY_HUNT_CROSSLEY_CONTACT_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-12 Stanford University and the Authors.        *
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

/** @file
 * Define the public interface to HuntCrossleyContact, a subsystem which
 * provides some minimal contact behavior.
 */

#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/ForceSubsystem.h"

#include <cassert>

namespace SimTK {

class MultibodySystem;

/**
 * This is a concrete subsystem that handles simple, frictionless contact situations
 * with a model due to Hunt & Crossley. See K. H. Hunt and F. R. E. Crossley, 
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
    explicit HuntCrossleyContact(MultibodySystem&);

    int addSphere(MobilizedBodyIndex body, const Vec3& center,
                  const Real& radius,
                  const Real& stiffness,     // E (plane strain)
                  const Real& dissipation);  // c (1/v)

    int addHalfSpace(MobilizedBodyIndex body, const UnitVec3& normal,
                     const Real& height,
                     const Real& stiffness,     // E (plane strain)
                     const Real& dissipation);  // c (1/v)

    SimTK_PIMPL_DOWNCAST(HuntCrossleyContact, ForceSubsystem);
private:
    class HuntCrossleyContactRep& updRep();
    const HuntCrossleyContactRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_HUNT_CROSSLEY_CONTACT_H_
