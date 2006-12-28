#ifndef SimTK_MOLECULAR_MECHANICS_SYSTEM_H_
#define SimTK_MOLECULAR_MECHANICS_SYSTEM_H_

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

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/System.h"
#include "simbody/internal/MultibodySystem.h"

#include <vector>

namespace SimTK {

class MatterSubsystem;
class DuMMForceFieldSubsystem;

/**
 * This is a particular kind of MultibodySystem, one intended for use in
 * moleular mechanics (MM). The defining feature is that in addition to
 * the mandatory MatterSubsystem common to all MultibodySystems, this one
 * will also have a single MolecularMechanicsForceSubsystem.
 */
class SimTK_SIMBODY_EXPORT MolecularMechanicsSystem : public MultibodySystem {
public:
    MolecularMechanicsSystem();
    MolecularMechanicsSystem(MatterSubsystem&, DuMMForceFieldSubsystem&);

    // Steals ownership of the source; returns subsystem ID number.
    int setMolecularMechanicsForceSubsystem(DuMMForceFieldSubsystem&);
    const DuMMForceFieldSubsystem& getMolecularMechanicsForceSubsystem() const;
    DuMMForceFieldSubsystem&       updMolecularMechanicsForceSubsystem();

    SimTK_PIMPL_DOWNCAST(MolecularMechanicsSystem, System);
private:
    class MolecularMechanicsSystemRep& updRep();
    const MolecularMechanicsSystemRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_MOLECULAR_MECHANICS_SYSTEM_H_
