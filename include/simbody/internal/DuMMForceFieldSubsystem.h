#ifndef SimTK_DUMM_FORCE_FIELD_SUBSYSTEM_H_
#define SimTK_DUMM_FORCE_FIELD_SUBSYSTEM_H_

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
 * Define the public interface to DuMMForceFieldSubsystem, a subsystem which
 * provides some minimal molecular mechanics-like capability.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/State.h"
#include "simbody/internal/System.h"
#include "simbody/internal/ForceSubsystem.h"

#include <cassert>

namespace SimTK {

/**
 * This is a concrete subsystem that provides basic molecular mechanics 
 * functionality FOR DEMO AND PROOF OF CONCEPT only!!! It is not likely
 * to perform well on anything.
 *
 * UNITS: must be consistent!
 *
 */
class SimTK_SIMBODY_API DuMMForceFieldSubsystem : public ForceSubsystem {
public:
    DuMMForceFieldSubsystem();

    // This fails if the atom type already exists.
    void defineAtomType(int type, Real mass,
                        Real vdwRadius, Real vdwWellDepth,
                        Real partialCharge);

    // Bond stretch parameters (between 2 atom types). This
    // fails if (type1,type2) or (type2,type1) has already been assigned.
    void defineBondStretch(int type1, int type2,
                           Real stiffness, Real nominalLength);

    // Bending angle parameters (among 3 atom types). This fails
    // if (type1,type2,type3) or (type3,type2,type1) has already been seen.
    void defineBondBending(int type1, int type2, int type3,
                           Real stiffness, Real nominalAngle);

    // At least one amplitude must be non-zero. 1-2-3-4, 4-3-2-1 
    // are the same.
    void defineTorsion(int type1, int type2, int type3, int type4,
        Real amplitude, Real phase, int periodicity,
        Real amp2, Real phase2, int period2,
        Real amp3, Real phase3, int period3);

    // The third atom is the central one to which the other
    // three are bonded; this is not the same in reverse order.
    void defineImproperTorsion(int type1, int type2, int type3, int type4,
        Real amplitude, Real phase, int periodicity,
        Real amp2, Real phase2, int period2,
        Real amp3, Real phase3, int period3);

    int addAtom(int body, int type, const Vec3& station);
    int addBond(int atom1, int atom2);


    SimTK_PIMPL_DOWNCAST(DuMMForceFieldSubsystem, ForceSubsystem);
private:
    class DuMMForceFieldSubsystemRep& updRep();
    const DuMMForceFieldSubsystemRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_DUMM_FORCE_FIELD_SUBSYSTEM_H_
