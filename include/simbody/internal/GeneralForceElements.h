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
class SimTK_SIMBODY_API GeneralForceElements : public ForceSubsystem {
public:
    GeneralForceElements();

    int addLinearTwoPointSpring(int body1, const Vec3& s1,
                                int body2, const Vec3& s2,
                                const Real& stiffness,
                                const Real& naturalLength);

    // Each generalized speed u_i feels a force -dampingFactor*u_i.
    // This is not physically meaningful but can be useful in some
    // circumstances just to drain energy out of the model when
    // the specific energy-draining mechanism is not important.
    // You can have more than one of these in which case the
    // dampingFactors are added. No individual dampingFactor is
    // allowed to be negative.
    int addGlobalMobilityDamping(const Real& dampingFactor);

    int addMobilitySpring(int body, int axis,
                          const Real& stiffness,
                          const Real& neutralValue);

    int addMobilityDamper(int body, int axis,
                          const Real& dampingFactor);

    SimTK_PIMPL_DOWNCAST(GeneralForceElements, ForceSubsystem);
private:
    class GeneralForceElementsRep& updRep();
    const GeneralForceElementsRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_GENERAL_FORCE_ELEMENTS_H_
