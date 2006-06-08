/* Copyright (c) 2006 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


/**@file
 *
 * Implementation of MultibodySystem, a concrete System.
 */

#include "simbody/internal/common.h"
#include "simbody/internal/MultibodySystem.h"

#include "MultibodySystemRep.h"

namespace SimTK {

    /////////////////////
    // MultibodySystem //
    /////////////////////

// Default constructor is inline and creates an empty handle.
// Default copy & assignment just copy the parent class.
// Default destructor destructs the parent class.

MultibodySystem::MultibodySystem() {
    rep = new MultibodySystemRep();
    rep->setMyHandle(*this);
}

MultibodySystem::MultibodySystem(MatterSubsystem& m, 
                                 ForceSubsystem& f)
{
    rep = new MultibodySystemRep();
    rep->setMyHandle(*this);

    setMatterSubsystem(m);
    setForceSubsystem(f);
}

bool MultibodySystem::project(State& s, Vector& y_err, 
             const Real& tol,
             const Real& dontProjectFac,
             const Real& targetTol
             ) const
{
    return MultibodySystemRep::downcast(*rep).project(
                s,y_err,tol,dontProjectFac,targetTol);
}

Real MultibodySystem::calcYErrorNorm(const State& s, const Vector& y_err) const {
    return MultibodySystemRep::downcast(*rep).calcYErrorNorm(s,y_err);
}

MatterSubsystem&       
MultibodySystem::setMatterSubsystem(MatterSubsystem& m) {
    return MultibodySystemRep::downcast(*rep).setMatterSubsystem(m);
}
ForceSubsystem& 
MultibodySystem::setForceSubsystem(ForceSubsystem& f) {
    return MultibodySystemRep::downcast(*rep).setForceSubsystem(f);
}

const MatterSubsystem&       
MultibodySystem::getMatterSubsystem() const {
    return MultibodySystemRep::downcast(*rep).getMatterSubsystem();
}

const ForceSubsystem& 
MultibodySystem::getForceSubsystem() const {
    return MultibodySystemRep::downcast(*rep).getForceSubsystem();
}

MatterSubsystem&       
MultibodySystem::updMatterSubsystem() {
    return MultibodySystemRep::downcast(*rep).updMatterSubsystem();
}

ForceSubsystem& 
MultibodySystem::updForceSubsystem() {
    return MultibodySystemRep::downcast(*rep).updForceSubsystem();
}

// TODO: camera facing, screen fixed, calculated geometry (e.g. line between stations
// on two different bodies, marker at system COM)
void MultibodySystem::addAnalyticGeometry(int body, const Transform& X_BG, const AnalyticGeometry& g) {
    MultibodySystemRep::downcast(*rep).addAnalyticGeometry(body,X_BG,g);
}

void MultibodySystem::addDecorativeGeometry(int body, const Transform& X_BG, const DecorativeGeometry& g) {
    MultibodySystemRep::downcast(*rep).addDecorativeGeometry(body,X_BG,g);
}

const Array<AnalyticGeometry>&   
MultibodySystem::getBodyAnalyticGeometry(int body) {
    return MultibodySystemRep::downcast(*rep).getBodyAnalyticGeometry(body);
}

const Array<DecorativeGeometry>& 
MultibodySystem::getBodyDecorativeGeometry(int body) {
    return MultibodySystemRep::downcast(*rep).getBodyDecorativeGeometry(body);
}

} // namespace SimTK

