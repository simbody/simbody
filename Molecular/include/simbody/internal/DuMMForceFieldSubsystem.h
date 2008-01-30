#ifndef SimTK_SIMBODY_DUMM_FORCE_FIELD_SUBSYSTEM_H_
#define SimTK_SIMBODY_DUMM_FORCE_FIELD_SUBSYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

/** @file
 * Define the public interface to DuMMForceFieldSubsystem, a subsystem which
 * provides some minimal molecular mechanics-like capability in a multi-rigid body
 * framework.
 */

#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/ForceSubsystem.h"

#include <cassert>

namespace SimTK {

class MolecularMechanicsSystem;

// Handy conversions. Note that these are compilation-unit statics, not members.
// That way we can be sure they are initialized before being used.
// To use these, multiply something in units on left of the "2" to get equivalent
// in units on right. E.g., 180*Deg2Rad gives Pi radians.
namespace DuMM {
static const Real Ang2Nm  = (Real)0.1L;
static const Real Nm2Ang  = (Real)10.L;
static const Real Deg2Rad = (Real)SimTK_DEGREE_TO_RADIAN;
static const Real Rad2Deg = (Real)SimTK_RADIAN_TO_DEGREE;
static const Real KJ2Kcal = (Real)SimTK_KJOULE_TO_KCAL;
static const Real Kcal2KJ = (Real)SimTK_KCAL_TO_KJOULE;

// There are several conventions for giving van der Waals parameters.
// Rmin is the radius at which the energy well minimum is seen
// (actually it is 1/2 the distance between atom centers for a pair
// of atoms of this class interacting with that minimum energy).
// This is *not* Sigma, which is the radius (half distance) at which the 
// energy crosses zero, that is, a little closer together than
// when the energy well is at maximum depth. To convert for LJ:
// Rmin = 2^(1/6) * Sigma.
static const Real Sigma2Radius = (Real)std::pow(2.L,  1.L/6.L); // sigma < radius
static const Real Radius2Sigma = (Real)std::pow(2.L, -1.L/6.L);
}


namespace DuMM {
    SimTK_DEFINE_UNIQUE_INDEX_TYPE(AtomIndex);
    SimTK_DEFINE_UNIQUE_INDEX_TYPE(AtomClassIndex);
    SimTK_DEFINE_UNIQUE_INDEX_TYPE(ChargedAtomTypeIndex);
    SimTK_DEFINE_UNIQUE_INDEX_TYPE(BondIndex);
    SimTK_DEFINE_UNIQUE_INDEX_TYPE(ClusterIndex);
} // namespace DuMM

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

    DuMMForceFieldSubsystem();
    explicit DuMMForceFieldSubsystem(MolecularMechanicsSystem&);

        // MOLECULE

    // Add a new atom to the model. The atom index number is returned; you don't get to
    // pick your own.
    DuMM::AtomIndex addAtom(DuMM::ChargedAtomTypeIndex chargedAtomTypeIx);

    // Note that these are atom index numbers, not atom classes or types.
    DuMM::BondIndex addBond(DuMM::AtomIndex atom1Ix, DuMM::AtomIndex atom2Ix);

    int    getNAtoms() const;
    Real   getAtomMass(DuMM::AtomIndex atomIx) const;
    int    getAtomElement(DuMM::AtomIndex atomIx) const;
    Real   getAtomRadius(DuMM::AtomIndex atomIx) const;
    Vec3   getAtomStationOnBody(DuMM::AtomIndex atomIx) const;
    Vec3   getAtomStationInCluster(DuMM::AtomIndex atomIx, DuMM::ClusterIndex clusterIx) const;
    MobilizedBodyIndex getAtomBody(DuMM::AtomIndex atomIx) const;
    Vec3   getAtomDefaultColor(DuMM::AtomIndex atomIx) const;

    int  getNBonds() const;

    // 'which' must be 0 or 1. 0 will return the lower-numbered atomIx.
    DuMM::AtomIndex  getBondAtom(DuMM::BondIndex bond, int which) const;

        // CLUSTERS

    // Create an empty cluster (rigid group of atoms). The cluster index number is returned;
    // you don't get to pick your own. The name is just for display; you must use the index
    // to reference the cluster. Every cluster has its own reference frame.
    DuMM::ClusterIndex createCluster(const char* clusterName);


    // Place an existing atom at a particular station in the local frame of a cluster. It
    // is fine for an atom to be in more than one cluster as long as only one of them ends up
    // attached to a body.
    void placeAtomInCluster(DuMM::AtomIndex atomIx, DuMM::ClusterIndex clusterIx, const Vec3& station);

    // Place a cluster (the child) in another cluster (the parent). The child's
    // local frame is placed at a given transform with respect to the parent's frame.
    void placeClusterInCluster(DuMM::ClusterIndex childClusterIndex, DuMM::ClusterIndex parentClusterIndex, 
                               const Transform& placement);

    // Calcuate the composite mass properties of a cluster, either in its own reference
    // frame or transformed to the indicated frame.
    MassProperties calcClusterMassProperties(DuMM::ClusterIndex clusterIx, const Transform& = Transform()) const;

    MobilizedBodyIndex    getClusterBody(DuMM::ClusterIndex clusterIx) const;
    Transform getClusterPlacementOnBody(DuMM::ClusterIndex clusterIx) const;
    Transform getClusterPlacementInCluster(DuMM::ClusterIndex childClusterIndex, DuMM::ClusterIndex parentClusterIndex) const;

        // BODIES

    void attachClusterToBody(DuMM::ClusterIndex clusterIx, MobilizedBodyIndex body, const Transform& = Transform());
    void attachAtomToBody   (DuMM::AtomIndex atomIx,    MobilizedBodyIndex body, const Vec3& station = Vec3(0));

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

    // Generate c++ code to reproduce forcefield parameters presently in memory
    void dumpCForcefieldParameters(std::ostream& os, const String& methodName = "loadParameters") const;

    void defineAtomClass(DuMM::AtomClassIndex atomClassIx, const char* atomClassName,
                         int elementNumber, int expectedValence,
                         Real vdwRadiusInNm, Real vdwWellDepthInKJ) {
         defineIncompleteAtomClass(atomClassIx, atomClassName, elementNumber, expectedValence);
         setAtomClassVdwParameters(atomClassIx, vdwRadiusInNm, vdwWellDepthInKJ);
    }

    // Same routine in Kcal/Angstrom (KA) unit system, i.e., radius
    // (still not sigma) is in nm, and well depth in kcal/mol.
    void defineAtomClass_KA(DuMM::AtomClassIndex atomClassIx, const char* atomClassName,
                            int element, int valence,
                            Real vdwRadiusInAng, Real vdwWellDepthInKcal)
    {
        defineAtomClass(atomClassIx, atomClassName, element, valence,
            vdwRadiusInAng*DuMM::Ang2Nm, vdwWellDepthInKcal*DuMM::Kcal2KJ);
    }
    // For backwards compatibility and compactness of expression, permit integer indices
    void defineAtomClass_KA(int atomClassIx, const char* atomClassName,
                            int element, int valence,
                            Real vdwRadiusInAng, Real vdwWellDepthInKcal)
    {
        defineAtomClass_KA((DuMM::AtomClassIndex)atomClassIx, atomClassName, element, valence, vdwRadiusInAng, vdwWellDepthInKcal);
    }


    // PartialCharge in units of e (charge on a proton); same in MD & KA
    void defineChargedAtomType(DuMM::ChargedAtomTypeIndex atomTypeIx, const char* atomTypeName,
        DuMM::AtomClassIndex atomClassIx, Real partialChargeInE) {
            defineIncompleteChargedAtomType(atomTypeIx, atomTypeName, atomClassIx);
            setChargedAtomTypeCharge(atomTypeIx, partialChargeInE);
    }

    void defineChargedAtomType_KA(DuMM::ChargedAtomTypeIndex atomTypeIx, const char* atomTypeName,
                                  DuMM::AtomClassIndex atomClassIx, Real partialChargeInE)
    {
        defineChargedAtomType(atomTypeIx, atomTypeName, atomClassIx, partialChargeInE); // easy!
    }
    // For backwards compatibility and compactness of expression, permit integer indices
    void defineChargedAtomType_KA(int atomTypeIx, const char* atomTypeName,
                                  int atomClassIx, Real partialChargeInE)
    {
        defineChargedAtomType_KA((DuMM::ChargedAtomTypeIndex) atomTypeIx, atomTypeName, (DuMM::AtomClassIndex)atomClassIx, partialChargeInE);
    }

    // Bond stretch parameters (between 2 atom classes). This
    // fails if (class1,class2) or (class2,class1) has already been assigned.
    // Stiffness (energy per length^2) in (kJ/mol)/nm^2
    // (note that energy is kx^2 using this definition,
    // while force is 2kx; note factor of 2 in force)
    void defineBondStretch(DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2,
                           Real stiffnessInKJperNmSq, Real nominalLengthInNm);

    // Here stiffness is in (kcal/mol)/A^2, and nominal length is in A (angstroms).
    void defineBondStretch_KA(DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2,
                              Real stiffnessInKcalPerAngSq, Real nominalLengthInAng)
    {
        defineBondStretch(class1, class2, 
                          stiffnessInKcalPerAngSq * DuMM::Kcal2KJ/square(DuMM::Ang2Nm),
                          nominalLengthInAng      * DuMM::Ang2Nm);
    }
    // For backwards compatibility and compactness of expression, permit integer indices
    void defineBondStretch_KA(int class1, int class2,
                              Real stiffnessInKcalPerAngSq, Real nominalLengthInAng)
    {
        defineBondStretch_KA((DuMM::AtomClassIndex)class1, (DuMM::AtomClassIndex)class2, stiffnessInKcalPerAngSq, nominalLengthInAng);
    }

    // Bending angle parameters (among 3 atom types). This fails
    // if (type1,type2,type3) or (type3,type2,type1) has already been seen.
    // Stiffness k (energy per rad^2) in (kJ/mol)/rad^2
    // Then energy is k*a^2 for angle a in radians,
    // while torque is 2ka; note factor of 2 in torque.
    // Note that the nominal angle is in degrees while the stiffness
    // is in radians. Odd, I know, but that seems to be how it's done!
    void defineBondBend(DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3,
                        Real stiffnessInKJPerRadSq, Real nominalAngleInDeg);

    // Here the stiffness is given in (kcal/mol)/rad^2.
    void defineBondBend_KA(DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3,
        Real stiffnessInKcalPerRadSq, Real nominalAngleInDeg) 
    {
        defineBondBend(class1,class2,class3,
                       stiffnessInKcalPerRadSq * DuMM::Kcal2KJ,
                       nominalAngleInDeg);
    }
    // For backwards compatibility and compactness of expression, permit integer indices
    void defineBondBend_KA(int class1, int class2, int class3,
        Real stiffnessInKcalPerRadSq, Real nominalAngleInDeg) 
    {
        defineBondBend_KA((DuMM::AtomClassIndex)class1, (DuMM::AtomClassIndex)class2, (DuMM::AtomClassIndex)class3,
                            stiffnessInKcalPerRadSq, nominalAngleInDeg);
    }

    // Only one term may have a given periodicity. The amplitudes are
    // in kJ/mol, with no factor of 1/2 expected (as is sometimes
    // the convention). 
    void defineBondTorsion
       (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
        int periodicity1, Real amp1InKJ, Real phase1InDegrees);
    void defineBondTorsion
       (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,  
        int periodicity1, Real amp1InKJ, Real phase1InDegrees,
        int periodicity2, Real amp2InKJ, Real phase2InDegrees);
    void defineBondTorsion
       (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
        int periodicity1, Real amp1InKJ, Real phase1InDegrees,
        int periodicity2, Real amp2InKJ, Real phase2InDegrees,
        int periodicity3, Real amp3InKJ, Real phase3InDegrees);

    // Here the amplitudes are given in kcal/mol.
    void defineBondTorsion_KA
       (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
        int periodicity1, Real amp1InKcal, Real phase1InDegrees)
    { 
        defineBondTorsion(class1,class2,class3,class4,
                          periodicity1, amp1InKcal * DuMM::Kcal2KJ, phase1InDegrees);
    }
    void defineBondTorsion_KA
       (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,  
        int periodicity1, Real amp1InKcal, Real phase1InDegrees,
        int periodicity2, Real amp2InKcal, Real phase2InDegrees)
    {
        defineBondTorsion(class1,class2,class3,class4,
                          periodicity1, amp1InKcal * DuMM::Kcal2KJ, phase1InDegrees,
                          periodicity2, amp2InKcal * DuMM::Kcal2KJ, phase2InDegrees);
    }
    void defineBondTorsion_KA
       (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4, 
        int periodicity1, Real amp1InKcal, Real phase1InDegrees,
        int periodicity2, Real amp2InKcal, Real phase2InDegrees,
        int periodicity3, Real amp3InKcal, Real phase3InDegrees)
    {
        defineBondTorsion(class1,class2,class3,class4,
                          periodicity1, amp1InKcal * DuMM::Kcal2KJ, phase1InDegrees,
                          periodicity2, amp2InKcal * DuMM::Kcal2KJ, phase2InDegrees,
                          periodicity3, amp3InKcal * DuMM::Kcal2KJ, phase3InDegrees);
    }
    // For backwards compatibility and compactness of expression, permit integer indices
    void defineBondTorsion_KA
       (int class1, int class2, int class3, int class4, 
        int periodicity1, Real amp1InKcal, Real phase1InDegrees)
    {
        defineBondTorsion_KA
               ((DuMM::AtomClassIndex)class1, (DuMM::AtomClassIndex)class2, (DuMM::AtomClassIndex)class3, (DuMM::AtomClassIndex)class4, 
                periodicity1, amp1InKcal, phase1InDegrees);
    }
    // For backwards compatibility and compactness of expression, permit integer indices
    void defineBondTorsion_KA
       (int class1, int class2, int class3, int class4, 
        int periodicity1, Real amp1InKcal, Real phase1InDegrees,
        int periodicity2, Real amp2InKcal, Real phase2InDegrees)
    {
        defineBondTorsion_KA
               ((DuMM::AtomClassIndex)class1, (DuMM::AtomClassIndex)class2, (DuMM::AtomClassIndex)class3, (DuMM::AtomClassIndex)class4, 
                periodicity1, amp1InKcal, phase1InDegrees,
                periodicity2, amp2InKcal, phase2InDegrees);
    }
    // For backwards compatibility and compactness of expression, permit integer indices
    void defineBondTorsion_KA
       (int class1, int class2, int class3, int class4, 
        int periodicity1, Real amp1InKcal, Real phase1InDegrees,
        int periodicity2, Real amp2InKcal, Real phase2InDegrees,
        int periodicity3, Real amp3InKcal, Real phase3InDegrees)
    {
        defineBondTorsion_KA
               ((DuMM::AtomClassIndex)class1, (DuMM::AtomClassIndex)class2, (DuMM::AtomClassIndex)class3, (DuMM::AtomClassIndex)class4, 
                periodicity1, amp1InKcal, phase1InDegrees,
                periodicity2, amp2InKcal, phase2InDegrees,
                periodicity3, amp3InKcal, phase3InDegrees);
    }

    // As with normal torsions, (see defineAmberImproperTorsion), only one term may have
    // a given periodicity. The amplitudes are in kJ/mol.
    void defineAmberImproperTorsion
       (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
        int periodicity1, Real amp1InKJ, Real phase1InDegrees);
    void defineAmberImproperTorsion
       (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
        int periodicity1, Real amp1InKJ, Real phase1InDegrees,
        int periodicity2, Real amp2InKJ, Real phase2InDegrees);
    void defineAmberImproperTorsion
       (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
        int periodicity1, Real amp1InKJ, Real phase1InDegrees,
        int periodicity2, Real amp2InKJ, Real phase2InDegrees,
        int periodicity3, Real amp3InKJ, Real phase3InDegrees);

    // Here the amplitudes are given in kcal/mol.
    void defineAmberImproperTorsion_KA
       (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
        int periodicity1, Real amp1InKcal, Real phase1InDegrees)
    {
        defineAmberImproperTorsion(class1,class2,class3,class4,
                          periodicity1, amp1InKcal * DuMM::Kcal2KJ, phase1InDegrees);
    }
    void defineAmberImproperTorsion_KA
       (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
        int periodicity1, Real amp1InKcal, Real phase1InDegrees,
        int periodicity2, Real amp2InKcal, Real phase2InDegrees)
    {
        defineAmberImproperTorsion(class1,class2,class3,class4,
                          periodicity1, amp1InKcal * DuMM::Kcal2KJ, phase1InDegrees,
                          periodicity2, amp2InKcal * DuMM::Kcal2KJ, phase2InDegrees);
    }
    void defineAmberImproperTorsion_KA
       (DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
        int periodicity1, Real amp1InKcal, Real phase1InDegrees,
        int periodicity2, Real amp2InKcal, Real phase2InDegrees,
        int periodicity3, Real amp3InKcal, Real phase3InDegrees)
    {
        defineAmberImproperTorsion(class1,class2,class3,class4,
                          periodicity1, amp1InKcal * DuMM::Kcal2KJ, phase1InDegrees,
                          periodicity2, amp2InKcal * DuMM::Kcal2KJ, phase2InDegrees,
                          periodicity3, amp3InKcal * DuMM::Kcal2KJ, phase3InDegrees);
    }

    // The third atom is the central one to which the other
    // three are bonded; this is not the same in reverse order.
    // TODO: not implemented
    // void defineImproperTorsion(DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
    //     Real amplitude, Real phase, int periodicity,
    //     Real amp2, Real phase2, int period2,
    //     Real amp3, Real phase3, int period3);
    // void defineImproperTorsion_KA(DuMM::AtomClassIndex class1, DuMM::AtomClassIndex class2, DuMM::AtomClassIndex class3, DuMM::AtomClassIndex class4,
    //     Real amplitude, Real phase, int periodicity,
    //     Real amp2, Real phase2, int period2,
    //     Real amp3, Real phase3, int period3);

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
    void setAmberImproperTorsionGlobalScaleFactor(Real);
    void setGbsaGlobalScaleFactor(Real);
    void setGbsaIncludeAceApproximation(bool);
    void setGbsaIncludeAceApproximationOn()  {setGbsaIncludeAceApproximation(true );}
    void setGbsaIncludeAceApproximationOff() {setGbsaIncludeAceApproximation(false);}

    void dump() const; // to stdout
    SimTK_PIMPL_DOWNCAST(DuMMForceFieldSubsystem, Subsystem);

protected:

    void defineIncompleteAtomClass(
        DuMM::AtomClassIndex classIx, 
        const char* name, 
        int elementNumber, 
        int valence);

    void defineIncompleteAtomClass_KA(
        DuMM::AtomClassIndex classIx, 
        const char* name, 
        int elementNumber, 
        int valence) 
    {
        defineIncompleteAtomClass(
            classIx, 
            name, 
            elementNumber, 
            valence
            );
    }

    void setAtomClassVdwParameters(DuMM::AtomClassIndex atomClassIx, Real vdwRadiusInNm, Real vdwWellDepthInKJPerMol);
    void setAtomClassVdwParameters_KA(DuMM::AtomClassIndex atomClassIx, Real radiusInAng, Real wellDepthInKcal) {
        setAtomClassVdwParameters(atomClassIx, radiusInAng*DuMM::Ang2Nm, wellDepthInKcal*DuMM::Kcal2KJ);
    }

    bool isValidAtomClass(DuMM::AtomClassIndex) const;
    
    void defineIncompleteChargedAtomType(
        DuMM::ChargedAtomTypeIndex typeIx, 
        const char* name,
        DuMM::AtomClassIndex classIx);

    void defineIncompleteChargedAtomType_KA(
        DuMM::ChargedAtomTypeIndex typeIx, 
        const char* name,
        DuMM::AtomClassIndex classIx) 
    {
        defineIncompleteChargedAtomType(typeIx, name, classIx);
    }

    void setChargedAtomTypeCharge(DuMM::ChargedAtomTypeIndex, Real charge);
    void setChargedAtomTypeCharge_KA(DuMM::ChargedAtomTypeIndex chargedAtomTypeIx, Real charge) {
        setChargedAtomTypeCharge(chargedAtomTypeIx, charge);
    }

private:
    class DuMMForceFieldSubsystemRep& updRep();
    const DuMMForceFieldSubsystemRep& getRep() const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_DUMM_FORCE_FIELD_SUBSYSTEM_H_
