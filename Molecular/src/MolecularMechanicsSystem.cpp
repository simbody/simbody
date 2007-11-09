/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2006-7 Stanford University and the Authors.         *
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

/**@file
 *
 * Implementation of MolecularMechanicsSystem, a kind of MultibodySystem.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/MolecularMechanicsSystem.h"

#include "MolecularMechanicsSystemRep.h"

namespace SimTK {


    ////////////////////////////////
    // MOLECULAR MECHANICS SYSTEM //
    ////////////////////////////////

class DuMMForceFieldSubsystem;

/*static*/ bool 
MolecularMechanicsSystem::isInstanceOf(const System& s) {
    return MolecularMechanicsSystemRep::isA(s.getSystemGuts());
}
/*static*/ const MolecularMechanicsSystem&
MolecularMechanicsSystem::downcast(const System& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const MolecularMechanicsSystem&>(s);
}
/*static*/ MolecularMechanicsSystem&
MolecularMechanicsSystem::updDowncast(System& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<MolecularMechanicsSystem&>(s);
}

const MolecularMechanicsSystemRep& 
MolecularMechanicsSystem::getRep() const {
    return dynamic_cast<const MolecularMechanicsSystemRep&>(getSystemGuts());
}
MolecularMechanicsSystemRep&       
MolecularMechanicsSystem::updRep() {
    return dynamic_cast<MolecularMechanicsSystemRep&>(updSystemGuts());
}

MolecularMechanicsSystem::MolecularMechanicsSystem() 
  : MultibodySystem(new MolecularMechanicsSystemRep())
{
}

MolecularMechanicsSystem::MolecularMechanicsSystem
   (SimbodyMatterSubsystem& matter, DuMMForceFieldSubsystem& mm)
  : MultibodySystem(new MolecularMechanicsSystemRep())
{
    setMatterSubsystem(matter);
    setMolecularMechanicsForceSubsystem(mm);
}

int MolecularMechanicsSystem::setMolecularMechanicsForceSubsystem(DuMMForceFieldSubsystem& mm) {
    return updRep().setMolecularMechanicsForceSubsystem(mm);
}

const DuMMForceFieldSubsystem&       
MolecularMechanicsSystem::getMolecularMechanicsForceSubsystem() const {
    return getRep().getMolecularMechanicsForceSubsystem();
}

DuMMForceFieldSubsystem&       
MolecularMechanicsSystem::updMolecularMechanicsForceSubsystem() {
    return updRep().updMolecularMechanicsForceSubsystem();
}

} // namespace SimTK

