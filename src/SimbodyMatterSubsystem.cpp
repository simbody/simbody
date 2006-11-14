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


/**@file
 *
 * Implementation of SimbodyMatterSubsystem, a concrete 
 * MatterSubsystem.
 */

#include "Simbody.h"
#include "RigidBodyTree.h"
class RigidBodyNode;

#include <string>
#include <iostream>
using std::cout;
using std::endl;

namespace SimTK {


/*static*/ bool 
SimbodyMatterSubsystem::isInstanceOf(const Subsystem& s) {
    return RigidBodyTree::isA(s.getRep());
}
/*static*/ const SimbodyMatterSubsystem&
SimbodyMatterSubsystem::downcast(const Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const SimbodyMatterSubsystem&>(s);
}
/*static*/ SimbodyMatterSubsystem&
SimbodyMatterSubsystem::updDowncast(Subsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<SimbodyMatterSubsystem&>(s);
}

const RigidBodyTree& 
SimbodyMatterSubsystem::getRep() const {
    return dynamic_cast<const RigidBodyTree&>(*rep);
}
RigidBodyTree&       
SimbodyMatterSubsystem::updRep() {
    return dynamic_cast<RigidBodyTree&>(*rep);
}

SimbodyMatterSubsystem::SimbodyMatterSubsystem() : MatterSubsystem() {
    rep = new RigidBodyTree();
    rep->setMyHandle(*this);
}


int SimbodyMatterSubsystem::addRigidBody(
    const MassProperties&     mp,
    const Transform&          bodyJointFrameInB,    // X_BJ
    int                       parent,
    const Transform&          parentJointFrameInP,  // X_PJb
    const Mobilizer&          mobilizer)
{
    const int save = getRep().nextUSlot;

    RigidBodyNode& pn = updRep().updRigidBodyNode(parent);
    const int rbIndex = updRep().addRigidBodyNode(pn,
        mp, parentJointFrameInP, bodyJointFrameInB, 
        mobilizer.getMobilizerType(), mobilizer.isReversed(),
        updRep().nextUSlot, updRep().nextUSqSlot, updRep().nextQSlot);

    //cout << "CREATED BODY " << rbIndex << ": U states " << save << "-" << getRep().nextUSlot-1 << endl;
    return rbIndex;
}

int SimbodyMatterSubsystem::addConstantDistanceConstraint
    (int parent, const Vec3& stationInP,
     int child,  const Vec3& stationInC,
     const Real& distance)
{
    return updRep().addConstantDistanceConstraint(
       getRep().getRigidBodyNode(parent), stationInP,
       getRep().getRigidBodyNode(child),  stationInC,
       distance);
}

int SimbodyMatterSubsystem::addCoincidentStationsConstraint
    (int parent, const Vec3& stationInP,
     int child,  const Vec3& stationInC)
{
    return updRep().addCoincidentStationsConstraint(
       getRep().getRigidBodyNode(parent), stationInP,
       getRep().getRigidBodyNode(child),  stationInC);
}

int SimbodyMatterSubsystem::addWeldConstraint
    (int parent, const Transform& frameInP,
     int child,  const Transform& frameInC)
{
    return updRep().addWeldConstraint(
       getRep().getRigidBodyNode(parent), frameInP,
       getRep().getRigidBodyNode(child),  frameInC);
}

// Note the lack of a State argument when completing construction.
void SimbodyMatterSubsystem::endConstruction() {updRep().endConstruction();}

// Convert spatial forces to internal equivalent, ignoring velocity and
// constraints.
void SimbodyMatterSubsystem::calcInternalGradientFromSpatial(const State& s,
    const Vector_<SpatialVec>& dEdR,
    Vector&                    dEdQ) const
{
    getRep().calcInternalGradientFromSpatial(s,dEdR,dEdQ);
}

// Convert spatial forces and centrifugal forces to an equivalent set
// of mobilizer forces, ignoring constraints.
void SimbodyMatterSubsystem::calcTreeEquivalentMobilityForces(const State& s, 
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    mobForces) const
{
    getRep().calcTreeEquivalentMobilityForces(s,bodyForces,mobForces);
}

Real SimbodyMatterSubsystem::calcKineticEnergy(const State& s) const {
    return getRep().calcKineticEnergy(s);
}

void SimbodyMatterSubsystem::calcTreeUDot(const State& s,
    const Vector&              jointForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    udot,
    Vector_<SpatialVec>&       A_GB) const
{
    Vector              netHingeForces; // unwanted side effect

    getRep().calcTreeAccelerations(s,jointForces,bodyForces,
        netHingeForces, A_GB, udot);
}


void SimbodyMatterSubsystem::calcMInverseF(const State& s,
    const Vector&        f,
    Vector&              udot,
    Vector_<SpatialVec>& A_GB) const
{
    getRep().calcMInverseF(s,f, A_GB, udot);
}

void SimbodyMatterSubsystem::calcQDot(const State& s,
    const Vector& u,
    Vector&       qdot) const
{
    getRep().calcQDot(s, u, qdot);
}

void SimbodyMatterSubsystem::calcQDotDot(const State& s,
    const Vector& udot,
    Vector&       qdotdot) const
{
    getRep().calcQDotDot(s, udot, qdotdot);
}

// Topological info. Note the lack of a State argument.
int SimbodyMatterSubsystem::getNBodies()        const {return getRep().getNBodies();}
int SimbodyMatterSubsystem::getTotalDOF()       const {return getRep().getTotalDOF();}
int SimbodyMatterSubsystem::getTotalQAlloc()    const {return getRep().getTotalQAlloc();}
int SimbodyMatterSubsystem::getNConstraints()   const {return getRep().getNConstraints();}
int SimbodyMatterSubsystem::getTotalMultAlloc() const {return getRep().getTotalMultAlloc();}

int SimbodyMatterSubsystem::getQIndex(int body) const {return getRep().getQIndex(body);}
int SimbodyMatterSubsystem::getQAlloc(int body) const {return getRep().getQAlloc(body);}
int SimbodyMatterSubsystem::getUIndex(int body) const {return getRep().getUIndex(body);}
int SimbodyMatterSubsystem::getDOF   (int body) const {return getRep().getDOF(body);}

int SimbodyMatterSubsystem::getMultIndex(int constraint) const {return getRep().getMultIndex(constraint);}
int SimbodyMatterSubsystem::getMaxNMult (int constraint) const {return getRep().getMaxNMult(constraint);}

// Modeling info.
void SimbodyMatterSubsystem::setUseEulerAngles(State& s, bool useAngles) const
  { getRep().setUseEulerAngles(s,useAngles); }
void SimbodyMatterSubsystem::setMobilizerIsPrescribed(State& s, int body, bool prescribed) const
  { getRep().setMobilizerIsPrescribed(s,body,prescribed); }
void SimbodyMatterSubsystem::setConstraintIsEnabled(State& s, int constraint, bool enabled) const
  { getRep().setConstraintIsEnabled(s,constraint,enabled); }
bool SimbodyMatterSubsystem::getUseEulerAngles(const State& s) const
  { return getRep().getUseEulerAngles(s); }
bool SimbodyMatterSubsystem::isMobilizerPrescribed(const State& s, int body) const
  { return getRep().isMobilizerPrescribed(s,body); }
bool SimbodyMatterSubsystem::isConstraintEnabled(const State& s, int constraint) const
  { return getRep().isConstraintEnabled(s,constraint); }


const Vector& SimbodyMatterSubsystem::getQ(const State& s) const {return getRep().getQ(s);}
const Vector& SimbodyMatterSubsystem::getU(const State& s) const {return getRep().getU(s);}

const Vector&
SimbodyMatterSubsystem::getAppliedMobilityForces(const State& s) const {
    return getRep().getAppliedMobilityForces(s);
}
const Vector_<SpatialVec>&
SimbodyMatterSubsystem::getAppliedBodyForces(const State& s) const {
    return getRep().getAppliedBodyForces(s);
}

void SimbodyMatterSubsystem::setQ(State& s, const Vector& q) const {getRep().setQ(s,q);}
void SimbodyMatterSubsystem::setU(State& s, const Vector& u) const {getRep().setU(s,u);}
Vector& SimbodyMatterSubsystem::updQ(State& s) const {return getRep().updQ(s);}
Vector& SimbodyMatterSubsystem::updU(State& s) const {return getRep().updU(s);}

void SimbodyMatterSubsystem::setMobilizerQ(State& s, int body, int axis, const Real& r) const
  { return getRep().setMobilizerQ(s,body,axis,r); }
void SimbodyMatterSubsystem::setMobilizerU(State& s, int body, int axis, const Real& r) const
  { return getRep().setMobilizerU(s,body,axis,r); }

const Real& SimbodyMatterSubsystem::getMobilizerQ(const State& s, int body, int axis) const
  { return getRep().getMobilizerQ(s,body,axis); }
const Real& SimbodyMatterSubsystem::getMobilizerU(const State& s, int body, int axis) const
  { return getRep().getMobilizerU(s,body,axis); }


void SimbodyMatterSubsystem::enforcePositionConstraints(State& s, const Real& requiredTol, const Real& desiredTol) const
  { getRep().enforcePositionConstraints(s, requiredTol, desiredTol); }
void SimbodyMatterSubsystem::enforceVelocityConstraints(State& s, const Real& requiredTol, const Real& desiredTol) const
  { getRep().enforceVelocityConstraints(s, requiredTol, desiredTol); }

const Transform&
SimbodyMatterSubsystem::getBodyPosition(const State& s, int body) const
  { return getRep().getBodyPosition(s,body); }

const SpatialVec&
SimbodyMatterSubsystem::getBodyVelocity(const State& s, int body) const {
    return getRep().getBodyVelocity(s,body);
}

const SpatialVec&
SimbodyMatterSubsystem::getCoriolisAcceleration(const State& s, int body) const {
    return getRep().getCoriolisAcceleration(s,body);
}
const SpatialVec&
SimbodyMatterSubsystem::getTotalCoriolisAcceleration(const State& s, int body) const {
    return getRep().getTotalCoriolisAcceleration(s,body);
}
const SpatialVec&
SimbodyMatterSubsystem::getGyroscopicForce(const State& s, int body) const {
    return getRep().getGyroscopicForce(s,body);
}
const SpatialVec&
SimbodyMatterSubsystem::getCentrifugalForces(const State& s, int body) const {
    return getRep().getCentrifugalForces(s,body);
}

const SpatialMat& 
SimbodyMatterSubsystem::getArticulatedBodyInertia(const State& s, int body) const {
    return getRep().getArticulatedBodyInertia(s,body);
}

const SpatialVec&
SimbodyMatterSubsystem::getBodyAcceleration(const State& s, int body) const
  { return getRep().getBodyAcceleration(s,body); }

const Vector& SimbodyMatterSubsystem::getQDot   (const State& s) const {return getRep().getQDot(s);}
const Vector& SimbodyMatterSubsystem::getUDot   (const State& s) const {return getRep().getUDot(s);}
const Vector& SimbodyMatterSubsystem::getQDotDot(const State& s) const {return getRep().getQDotDot(s);}

} // namespace SimTK

