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
#include "simbody/internal/System.h"
#include "simbody/internal/ForceSubsystem.h"

#include <cassert>

namespace SimTK {

/**
 * This is a concrete subsystem that provides basic molecular mechanics 
 * functionality FOR DEMO AND PROOF OF CONCEPT only!!! It is not likely
 * to perform well on anything.
 *
 * UNITS: This subsystem requires that the system be modeled in "MD units"
 * of nanometers, daltons (g/mol), and picoseconds, yielding consistent
 * energy units of kJ/mol==(Da-nm^2/ps^2), and force in kJ/mol-nm. Charge
 * is in proton charge units e, and angles are in radians.
 * For convenience, we allow the force field to be defined in "KA" units,
 * that is, angstroms instead of nanometers, and energy in kcal rather
 * than kJ, and we also allow angles to be supplied in degrees. However,
 * these are immediately converted to the MD units described above.
 *
 */

class SimTK_SIMBODY_EXPORT DuMMForceFieldSubsystem : public ForceSubsystem {
public:
    enum VdwMixingRule {
        WaldmanHagler       = 1,    // Our default
        HalgrenHHG          = 2,    // MMFF, AMOEBA
        Jorgensen           = 3,    // OPLS
        LorentzBerthelot    = 4,    // AMBER, CHARMM
        Kong                = 5
    };
    const char* getVdwMixingRuleName(VdwMixingRule) const;

    static const Real Ang2Nm;       // multiply Angstroms by this to get nanometers
    static const Real Nm2Ang;       //   or multiply nanometers by this to get Angstroms
    static const Real Kcal2KJ;      // multiply kilocalories by this to get kilojoules
    static const Real KJ2Kcal;      //   or multiply kilojoules by this to get kilocalories
    static const Real Deg2Rad;      // multiply degrees by this to convert to radians
    static const Real Rad2Deg;      //   or multiply radians by this to get degrees
    static const Real Sigma2Radius; // multiply vdw sigma by this to get vdw radius
    static const Real Radius2Sigma; //   or multiply vdw radius by this to get vdw sigma

    DuMMForceFieldSubsystem();

        // MOLECULE

    // Add a new atom to the model. The atom Id number is returned; you don't get to
    // pick your own.
    int addAtom(int chargedAtomTypeId);

    // Note that these are atom Id numbers, not atom classes or types.
    int  addBond(int atom1Id, int atom2Id);

    int  getNAtoms() const;
    Real getAtomMass(int atomId) const;
    int  getAtomElement(int atomId) const;
    Real getAtomRadius(int atomId) const;
    Vec3 getAtomStationOnBody(int atomId) const;
    Vec3 getAtomStationInCluster(int atomId, int clusterId) const;
    int  getAtomBody(int atomId) const;
    Vec3 getAtomDefaultColor(int atomId) const;

    int  getNBonds() const;

    // 'which' must be 0 or 1. 0 will return the lower-numbered atomId.
    int  getBondAtom(int bond, int which) const;

        // CLUSTERS

    // Create an empty cluster (rigid group of atoms). The cluster Id number is returned;
    // you don't get to pick your own. The name is just for display; you must use the Id
    // to reference the cluster. Every cluster has its own reference frame.
    int createCluster(const char* clusterName);


    // Place an existing atom at a particular station in the local frame of a cluster. It
    // is fine for an atom to be in more than one cluster as long as only one of them ends up
    // attached to a body.
    void placeAtomInCluster(int atomId, int clusterId, const Vec3& station);

    // Place a cluster (the child) in another cluster (the parent). The child's
    // local frame is placed at a given transform with respect to the parent's frame.
    void placeClusterInCluster(int childClusterId, int parentClusterId, 
                               const Transform& placement);

    // Calcuate the composite mass properties of a cluster, either in its own reference
    // frame or transformed to the indicated frame.
    MassProperties calcClusterMassProperties(int clusterId, const Transform& = Transform()) const;

    int       getClusterBody(int clusterId) const;
    Transform getClusterPlacementOnBody(int clusterId) const;
    Transform getClusterPlacementInCluster(int childClusterId, int parentClusterId) const;

        // BODIES

    void attachClusterToBody(int clusterId, int body, const Transform& = Transform());
    void attachAtomToBody   (int atomId,    int body, const Vec3& station = Vec3(0));

        // DEFINE FORCE FIELD PARAMETERS
    
    // Atom classes are used for sets of atoms which share some properties.
    // These are: the element (as atomic number), expected valence, 
    // van der Waals parameters, and behavior in bonded situations.
    // Charge is not included in atom class but in a second more detailed
    // classification level called ChargedAtomType.
    //
    // This fails if the atom class already exists.
    // VdwRadius is given as Rmin, the radius at which the energy
    // well minimum is seen (actually it is 1/2 the distance between
    // atom centers for a pair of atoms of this class). This is *not*
    // Sigma, which is the radius (half distance) at which the 
    // energy crosses zero, that is, a little closer together than
    // when the energy well is at maximum depth.
    // To convert for LJ: Rmin = 2^(1/6) * Sigma.
    // The radius is in nm, the well depth in kJ/mol.

    void defineAtomClass(int atomClassId, const char* atomClassName,
                         int elementNumber, int expectedValence,
                         Real vdwRadiusInNm, Real vdwWellDepthInKJ);

    // Same routine in Kcal/Angstrom (KA) unit system, i.e., radius
    // (still not sigma) is in nm, and well depth in kcal/mol.
    void defineAtomClass_KA(int atomClassId, const char* atomClassName,
                            int element, int valence,
                            Real vdwRadiusInAng, Real vdwWellDepthInKcal)
    {
        defineAtomClass(atomClassId, atomClassName, element, valence,
                        vdwRadiusInAng*Ang2Nm, vdwWellDepthInKcal*Kcal2KJ);
    }

    // PartialCharge in units of e (charge on a proton); same in MD & KA
    void defineChargedAtomType(int atomTypeId, const char* atomTypeName,
                               int atomClassId, Real partialChargeInE);

    void defineChargedAtomType_KA(int atomTypeId, const char* atomTypeName,
                                  int atomClassId, Real partialChargeInE)
    {
        defineChargedAtomType(atomTypeId, atomTypeName, atomClassId, partialChargeInE); // easy!
    }

    // Bond stretch parameters (between 2 atom classes). This
    // fails if (class1,class2) or (class2,class1) has already been assigned.
    // Stiffness (energy per length^2) in (kJ/mol)/nm^2
    // (note that energy is kx^2 using this definition,
    // while force is 2kx; note factor of 2 in force)
    void defineBondStretch(int class1, int class2,
                           Real stiffnessInKJperNmSq, Real nominalLengthInNm);

    // Here stiffness is in (kcal/mol)/A^2, and nominal length is in A (angstroms).
    void defineBondStretch_KA(int class1, int class2,
                              Real stiffnessInKcalPerAngSq, Real nominalLengthInAng)
    {
        defineBondStretch(class1, class2, 
                          stiffnessInKcalPerAngSq * Kcal2KJ/(Ang2Nm*Ang2Nm),
                          nominalLengthInAng      * Ang2Nm);
    }

    // Bending angle parameters (among 3 atom types). This fails
    // if (type1,type2,type3) or (type3,type2,type1) has already been seen.
    // Stiffness k (energy per rad^2) in (kJ/mol)/rad^2
    // Then energy is k*a^2 for angle a in radians,
    // while torque is 2ka; note factor of 2 in torque.
    // Note that the nominal angle is in degrees while the stiffness
    // is in radians. Odd, I know, but that seems to be how it's done!
    void defineBondBend(int class1, int class2, int class3,
                        Real stiffnessInKJPerRadSq, Real nominalAngleInDeg);

    // Here the stiffness is given in (kcal/mol)/rad^2.
    void defineBondBend_KA(int class1, int class2, int class3,
        Real stiffnessInKcalPerRadSq, Real nominalAngleInDeg) 
    {
        defineBondBend(class1,class2,class3,
                       stiffnessInKcalPerRadSq * Kcal2KJ,
                       nominalAngleInDeg);
    }

    // Only one term may have a given periodicity. The amplitudes are
    // in kJ/mol, with no factor of 1/2 expected (as is sometimes
    // the convention). 
    void defineBondTorsion
       (int class1, int class2, int class3, int class4, 
        int periodicity1, Real amp1InKJ, Real phase1InDegrees);
    void defineBondTorsion
       (int class1, int class2, int class3, int class4,  
        int periodicity1, Real amp1InKJ, Real phase1InDegrees,
        int periodicity2, Real amp2InKJ, Real phase2InDegrees);
    void defineBondTorsion
       (int class1, int class2, int class3, int class4, 
        int periodicity1, Real amp1InKJ, Real phase1InDegrees,
        int periodicity2, Real amp2InKJ, Real phase2InDegrees,
        int periodicity3, Real amp3InKJ, Real phase3InDegrees);

    // Here the amplitudes are given in kcal/mol.
    void defineBondTorsion_KA
       (int class1, int class2, int class3, int class4, 
        int periodicity1, Real amp1InKcal, Real phase1InDegrees)
    { 
        defineBondTorsion(class1,class2,class3,class4,
                          periodicity1, amp1InKcal * Kcal2KJ, phase1InDegrees);
    }
    void defineBondTorsion_KA
       (int class1, int class2, int class3, int class4,  
        int periodicity1, Real amp1InKcal, Real phase1InDegrees,
        int periodicity2, Real amp2InKcal, Real phase2InDegrees)
    {
        defineBondTorsion(class1,class2,class3,class4,
                          periodicity1, amp1InKcal * Kcal2KJ, phase1InDegrees,
                          periodicity2, amp2InKcal * Kcal2KJ, phase2InDegrees);
    }
    void defineBondTorsion_KA
       (int class1, int class2, int class3, int class4, 
        int periodicity1, Real amp1InKcal, Real phase1InDegrees,
        int periodicity2, Real amp2InKcal, Real phase2InDegrees,
        int periodicity3, Real amp3InKcal, Real phase3InDegrees)
    {
        defineBondTorsion(class1,class2,class3,class4,
                          periodicity1, amp1InKcal * Kcal2KJ, phase1InDegrees,
                          periodicity2, amp2InKcal * Kcal2KJ, phase2InDegrees,
                          periodicity3, amp3InKcal * Kcal2KJ, phase3InDegrees);
    }

    // The third atom is the central one to which the other
    // three are bonded; this is not the same in reverse order.
    // TODO: not implemented
    void defineImproperTorsion(int class1, int class2, int class3, int class4,
        Real amplitude, Real phase, int periodicity,
        Real amp2, Real phase2, int period2,
        Real amp3, Real phase3, int period3);
    void defineImproperTorsion_KA(int class1, int class2, int class3, int class4,
        Real amplitude, Real phase, int periodicity,
        Real amp2, Real phase2, int period2,
        Real amp3, Real phase3, int period3);

    void setVdwMixingRule(VdwMixingRule); // default WaldmanHagler
    VdwMixingRule getVdwMixingRule() const;

    void setVdw12ScaleFactor(Real); // default 0
    void setVdw13ScaleFactor(Real); // default 0
    void setVdw14ScaleFactor(Real); // default 1
    void setVdw15ScaleFactor(Real); // default 1

    void setCoulomb12ScaleFactor(Real); // default 0
    void setCoulomb13ScaleFactor(Real); // default 0
    void setCoulomb14ScaleFactor(Real); // default 1
    void setCoulomb15ScaleFactor(Real); // default 1

    // These can be used to weaken or disable (or magnify)
    // individual force field terms.
    // These are always 1 for correct implementation of any force field.
    void setVdwGlobalScaleFactor(Real);
    void setCoulombGlobalScaleFactor(Real);
    void setBondStretchGlobalScaleFactor(Real);
    void setBondBendGlobalScaleFactor(Real);
    void setBondTorsionGlobalScaleFactor(Real);

    void dump() const; // to stdout
    SimTK_PIMPL_DOWNCAST(DuMMForceFieldSubsystem, Subsystem);
private:
    class DuMMForceFieldSubsystemRep& updRep();
    const DuMMForceFieldSubsystemRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_DUMM_FORCE_FIELD_SUBSYSTEM_H_
