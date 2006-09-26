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
 * provides some minimal molecular mechanics-like capability in a multi-rigid body
 * framework.
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
 * UNITS: must be as specified. TODO: allow different units.
 * Note: these are not consistent units; they are converted
 * internally so that correct energy and force units are applied
 * to the multibody system. Basic units are mass in Da (g/mol),
 * length in Angstroms, angles in radians, time in ps. Consistent
 * energy is then Da-A^2/ps^2, force is Da-A/ps^2. But here we expect
 * energy in Kcal/mol, force in (Kcal/mol)/A.
 *
 */

class SimTK_SIMBODY_API DuMMForceFieldSubsystem : public ForceSubsystem {
public:
    DuMMForceFieldSubsystem();

    // Atom classes are used for sets of atoms which share some properties.
    // These are: the element (as atomic number), expected valence, 
    // van der Waals parameters, and behavior in bonded situations.
    // Charge is not included in atom class but in a second classification
    // level called ChargedAtomType.
    //
    // This fails if the atom class already exists.
    //   mass in Da (g/mol)
    //   vdwRadius as Rmin, *not* Sigma, in Angstroms
    //     (i.e. 2*vdwRadius is the center-center separation
    //      at which the minimum energy occurs)
    //     To convert for LJ: Rmin = 2^(1/6) * Sigma
    //   vdwWellDepth potential minimum, in Kcal/mol

    void defineAtomClass(int atomClassId, const char* atomClassName,
                         int element, int valence,
                         Real vdwRadius, Real vdwWellDepth);

    //   partialCharge in units of e (charge on a proton) 
    void defineChargedAtomType(int atomTypeId, const char* atomTypeName,
                               int atomClassId, Real partialCharge);

    // Create an empty cluster (rigid group of atoms). The cluster Id number is returned;
    // you don't get to pick your own. The name is just for display; you must use the Id
    // to reference the cluster. Every cluster has its own reference frame.
    int createCluster(const char* clusterName);

    // Add a new atom to the model. The atom Id number is returned; you don't get to
    // pick your own.
    int addAtom(int chargedAtomTypeId);

    // Place an existing atom at a particular station in the local frame of a cluster. It
    // is fine for an atom to be in more than one cluster as long as only one of them ends up
    // attached to a body.
    void placeAtomInCluster(int atomId, int clusterId, const Vec3& station);

    // Place a cluster (the child) in another cluster (the parent). The child's
    // local frame is placed at a given transform with respect to the parent's frame.
    void placeClusterInCluster(int childClusterId, int parentClusterId, 
                               const Transform& placement);

    void attachClusterToBody(int clusterId, int body, const Transform& = Transform());
    void attachAtomToBody   (int atomId,    int body, const Vec3& station = Vec3(0));

    // Calcuate the composite mass properties of a cluster, either in its own reference
    // frame or transformed to the indicated frame.
    MassProperties calcClusterMassProperties(int clusterId, const Transform& = Transform()) const;

    // Bond stretch parameters (between 2 atom classes). This
    // fails if (class1,class2) or (class2,class1) has already been assigned.
    //   stiffness (energy per length^2) in (Kcal/mol)/A^2
    //     (note that energy is kx^2 using this definition,
    //      while force is 2kx; note factor of 2 in force)
    //   nominalLength in Angstroms
    void defineBondStretch(int class1, int class2,
                           Real stiffness, Real nominalLength);

    // Bending angle parameters (among 3 atom types). This fails
    // if (type1,type2,type3) or (type3,type2,type1) has already been seen.
    //   stiffness k (energy per degree^2) in (Kcal/mol)/Degree^2 (NOT Radians)
    //     Let k'=k*(180/pi)^2 (i.e. k' is in energy per radian^2). 
    //     Then energy is k' a^2 for angle a in radians,
    //     while torque is 2k'a; note factor of 2 in torque.
    //   nominalAngle in Degrees
    void defineBondBend(int class1, int class2, int class3,
                        Real stiffness, Real nominalAngle);


    // Only one term may have a given periodicity.
    void defineBondTorsion
       (int class1, int class2, int class3, int class4, 
        int periodicity1, Real amp1InKcal, Real phase1InDegrees);
    void defineBondTorsion
       (int class1, int class2, int class3, int class4,  
        int periodicity1, Real amp1InKcal, Real phase1InDegrees,
        int periodicity2, Real amp2InKcal, Real phase2InDegrees);
    void defineBondTorsion
       (int class1, int class2, int class3, int class4, 
        int periodicity1, Real amp1InKcal, Real phase1InDegrees,
        int periodicity2, Real amp2InKcal, Real phase2InDegrees,
        int periodicity3, Real amp3InKcal, Real phase3InDegrees);

    // The third atom is the central one to which the other
    // three are bonded; this is not the same in reverse order.
    // TODO: not implemented
    void defineImproperTorsion(int class1, int class2, int class3, int class4,
        Real amplitude, Real phase, int periodicity,
        Real amp2, Real phase2, int period2,
        Real amp3, Real phase3, int period3);

    void setVdw12ScaleFactor(Real); // default 0
    void setVdw13ScaleFactor(Real); // default 0
    void setVdw14ScaleFactor(Real); // default 1
    void setVdw15ScaleFactor(Real); // default 1

    void setCoulomb12ScaleFactor(Real); // default 0
    void setCoulomb13ScaleFactor(Real); // default 0
    void setCoulomb14ScaleFactor(Real); // default 1
    void setCoulomb15ScaleFactor(Real); // default 1


    // Note that these are atom Id numbers, not atom classes or types.
    int  addBond(int atomId1, int atomId2);
    int  getNBonds() const;

    // 'which' must be 0 or 1. 0 will return the lower-numbered atomId.
    int  getBondAtom(int bond, int which) const;

    int  getNAtoms() const;
    Real getAtomMass(int atomId) const;
    Real getAtomRadius(int atomId) const;
    Vec3 getAtomStationOnBody(int atomId) const;
    Vec3 getAtomStationInCluster(int atomId, int clusterId) const;
    int  getAtomBody(int atomId) const;
    Vec3 getAtomDefaultColor(int atomId) const;

    void dump() const; // to stdout

    SimTK_PIMPL_DOWNCAST(DuMMForceFieldSubsystem, ForceSubsystem);
private:
    class DuMMForceFieldSubsystemRep& updRep();
    const DuMMForceFieldSubsystemRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_DUMM_FORCE_FIELD_SUBSYSTEM_H_
