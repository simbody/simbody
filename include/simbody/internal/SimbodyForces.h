#ifndef SimTK_SIMBODY_FORCES_H_
#define SimTK_SIMBODY_FORCES_H_

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
 * This is a concrete subsystem which applies a single linear, two-point spring.
 */
class SimTK_SIMBODY_API TwoPointSpringSubsystem : public ForceSubsystem {
public:
    TwoPointSpringSubsystem();

    int addLinearTwoPointSpring(int body1, const Vec3& s1,
                                int body2, const Vec3& s2,
                                const Real& stiffness,
                                const Real& naturalLength);

    int addGlobalMobilityDamping(const Real& dampingFactor);

    int addMobilitySpring(int body, int axis,
                          const Real& stiffness,
                          const Real& neutralValue);

    int addMobilityDamper(int body, int axis,
                          const Real& dampingFactor);

    SimTK_PIMPL_DOWNCAST(TwoPointSpringSubsystem, ForceSubsystem);
private:
    class TwoPointSpringSubsystemRep& updRep();
    const TwoPointSpringSubsystemRep& getRep() const;
};

/**
 * This is a concrete subsystem which applies no forces.
 */
class SimTK_SIMBODY_API EmptyForcesSubsystem : public ForceSubsystem {
public:
    EmptyForcesSubsystem();

    // These are just the defaults but are nice to have explicitly for debugging.
    ~EmptyForcesSubsystem() {
    }
    EmptyForcesSubsystem(const EmptyForcesSubsystem& e) 
      : ForceSubsystem(e) {
    }
    EmptyForcesSubsystem& operator=(const EmptyForcesSubsystem& e) {
        ForceSubsystem::operator=(e);
        return *this;
    }

    SimTK_PIMPL_DOWNCAST(EmptyForcesSubsystem, ForceSubsystem);
private:
    class EmptyForcesSubsystemRep& updRep();
    const EmptyForcesSubsystemRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCES_H_
