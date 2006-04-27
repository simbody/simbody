/**@file
 *
 * Implementation of SimbodySubsystem.
 */

#include "Simbody.h"
#include "RigidBodyTree.h"
class RigidBodyNode;

#include <string>
#include <iostream>
using std::cout;
using std::endl;

namespace SimTK {

SimbodySubsystem::SimbodySubsystem() : rep(0) {
    rep = new RigidBodyTree();
}

SimbodySubsystem::SimbodySubsystem(const SimbodySubsystem& src) : rep(0) {
    assert(src.rep);
    rep = new RigidBodyTree(*src.rep);
}

SimbodySubsystem& SimbodySubsystem::operator=(const SimbodySubsystem& src) {
    if (&src == this) return *this;
    assert(src.rep);
    delete rep; rep=0;
    rep = new RigidBodyTree(*src.rep);
    return *this;
}

SimbodySubsystem::~SimbodySubsystem() {
    delete rep; rep=0;
}

int SimbodySubsystem::addRigidBody(
    const MassProperties&     mp,
    const Transform&          bodyJointFrameInB,    // X_BJ
    int                       parent,
    const Transform&          parentJointFrameInP,  // X_PJb
    const JointSpecification& joint)
{
    const int save = rep->nextUSlot;

    RigidBodyNode& pn = rep->updRigidBodyNode(parent);
    const int rbIndex = rep->addRigidBodyNode(pn,
        mp, parentJointFrameInP, bodyJointFrameInB, joint.getJointType(), joint.isReversed(),
        rep->nextUSlot, rep->nextUSqSlot, rep->nextQSlot);

    //cout << "CREATED BODY " << rbIndex << ": U states " << save << "-" << rep->nextUSlot-1 << endl;
    return rbIndex;
}

int SimbodySubsystem::addConstantDistanceConstraint
    (int parent, const Vec3& stationInP,
     int child,  const Vec3& stationInC,
     const Real& distance)
{
    return rep->addConstantDistanceConstraint(
       rep->getRigidBodyNode(parent), stationInP,
       rep->getRigidBodyNode(child),  stationInC,
       distance);
}

int SimbodySubsystem::addCoincidentStationsConstraint
    (int parent, const Vec3& stationInP,
     int child,  const Vec3& stationInC)
{
    return rep->addCoincidentStationsConstraint(
       rep->getRigidBodyNode(parent), stationInP,
       rep->getRigidBodyNode(child),  stationInC);
}

int SimbodySubsystem::addWeldConstraint
    (int parent, const Transform& frameInP,
     int child,  const Transform& frameInC)
{
    return rep->addWeldConstraint(
       rep->getRigidBodyNode(parent), frameInP,
       rep->getRigidBodyNode(child),  frameInC);
}

// Note the lack of a State argument when completing construction.
void SimbodySubsystem::endConstruction()                          {rep->endConstruction();}

void SimbodySubsystem::realizeConstruction(State& s)        const {rep->realizeConstruction(s);}
void SimbodySubsystem::realizeModeling    (State& s)        const {rep->realizeModeling(s);}

void SimbodySubsystem::realizeParameters   (const State& s) const {rep->realizeParameters(s);}
void SimbodySubsystem::realizeTime         (const State& s) const {rep->realizeTime(s);}
void SimbodySubsystem::realizeConfiguration(const State& s) const {rep->realizeConfiguration(s);}
void SimbodySubsystem::realizeMotion       (const State& s) const {rep->realizeMotion(s);}
void SimbodySubsystem::realizeDynamics     (const State& s) const {rep->realizeDynamics(s);}
void SimbodySubsystem::realizeReaction     (const State& s) const {rep->realizeReaction(s);}

void SimbodySubsystem::realize(const State& s, Stage g) const {
    while (s.getStage() < g) {
        switch (s.getStage()) {
        case Stage::Allocated:    realizeConstruction(const_cast<State&>(s)); break;
        case Stage::Built:        realizeModeling    (const_cast<State&>(s)); break;
        case Stage::Modeled:      realizeParameters(s);    break;
        case Stage::Parametrized: realizeTime(s);          break;
        case Stage::Timed:        realizeConfiguration(s); break;
        case Stage::Configured:   realizeMotion(s);        break;
        case Stage::Moving:       realizeDynamics(s);      break;
        case Stage::Dynamics:     realizeReaction(s);      break;
        default: assert(!"SimbodySubsystem::realize(): bad stage");
        }
        s.advanceToStage(s.getStage().next());
    }
}

void SimbodySubsystem::calcInternalGradientFromSpatial(const State& s,
    const Vector_<SpatialVec>& dEdR,
    Vector&                    dEdQ) const
{
    rep->calcInternalGradientFromSpatial(s,dEdR,dEdQ);
}

void SimbodySubsystem::calcTreeEquivalentJointForces(const State& s, 
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    jointForces) const
{
    rep->calcTreeEquivalentJointForces(s,bodyForces,jointForces);
}

Real SimbodySubsystem::calcKineticEnergy(const State& s) const {
    return rep->calcKineticEnergy(s);
}

void SimbodySubsystem::calcTreeUDot(const State& s,
    const Vector&              jointForces,
    const Vector_<SpatialVec>& bodyForces,
    Vector&                    udot) const
{
    Vector              netHingeForces; // unwanted side effects
    Vector_<SpatialVec> A_GB;

    rep->calcTreeAccelerations(s,jointForces,bodyForces,
        netHingeForces, A_GB, udot);
}

void SimbodySubsystem::calcQDot(const State& s,
    const Vector& u,
    Vector&       qdot) const
{
    rep->calcQDot(s, u, qdot);
}

void SimbodySubsystem::calcQDotDot(const State& s,
    const Vector& udot,
    Vector&       qdotdot) const
{
    rep->calcQDotDot(s, udot, qdotdot);
}

// Topological info. Note the lack of a State argument.
int SimbodySubsystem::getNBodies()        const {return rep->getNBodies();}
int SimbodySubsystem::getTotalDOF()       const {return rep->getTotalDOF();}
int SimbodySubsystem::getTotalQAlloc()    const {return rep->getTotalQAlloc();}
int SimbodySubsystem::getNConstraints()   const {return rep->getNConstraints();}
int SimbodySubsystem::getTotalMultAlloc() const {return rep->getTotalMultAlloc();}

int SimbodySubsystem::getQIndex(int body) const {return rep->getQIndex(body);}
int SimbodySubsystem::getQAlloc(int body) const {return rep->getQAlloc(body);}
int SimbodySubsystem::getUIndex(int body) const {return rep->getUIndex(body);}
int SimbodySubsystem::getDOF   (int body) const {return rep->getDOF(body);}

int SimbodySubsystem::getMultIndex(int constraint) const {return rep->getMultIndex(constraint);}
int SimbodySubsystem::getMaxNMult (int constraint) const {return rep->getMaxNMult(constraint);}

// Modeling info.
void SimbodySubsystem::setUseEulerAngles(State& s, bool useAngles) const
  { rep->setUseEulerAngles(s,useAngles); }
void SimbodySubsystem::setJointIsPrescribed(State& s, int joint, bool prescribed) const
  { rep->setJointIsPrescribed(s,joint,prescribed); }
void SimbodySubsystem::setConstraintIsEnabled(State& s, int constraint, bool enabled) const
  { rep->setConstraintIsEnabled(s,constraint,enabled); }
bool SimbodySubsystem::getUseEulerAngles(const State& s) const
  { return rep->getUseEulerAngles(s); }
bool SimbodySubsystem::isJointPrescribed(const State& s, int joint) const
  { return rep->isJointPrescribed(s,joint); }
bool SimbodySubsystem::isConstraintEnabled(const State& s, int constraint) const
  { return rep->isConstraintEnabled(s,constraint); }


const VectorView SimbodySubsystem::getQ(const State& s) const {return rep->getQ(s);}
const VectorView SimbodySubsystem::getU(const State& s) const {return rep->getU(s);}

const Vector&
SimbodySubsystem::getAppliedJointForces(const State& s) const {
    return rep->getAppliedJointForces(s);
}
const Vector_<SpatialVec>&
SimbodySubsystem::getAppliedBodyForces(const State& s) const {
    return rep->getAppliedBodyForces(s);
}

void SimbodySubsystem::setQ(State& s, const Vector& q) const {rep->setQ(s,q);}
void SimbodySubsystem::setU(State& s, const Vector& u) const {rep->setU(s,u);}
VectorView SimbodySubsystem::updQ(State& s) const {return rep->updQ(s);}
VectorView SimbodySubsystem::updU(State& s) const {return rep->updU(s);}

void SimbodySubsystem::setJointQ(State& s, int body, int axis, const Real& r) const
  { return rep->setJointQ(s,body,axis,r); }
void SimbodySubsystem::setJointU(State& s, int body, int axis, const Real& r) const
  { return rep->setJointU(s,body,axis,r); }

const Real& SimbodySubsystem::getJointQ(const State& s, int body, int axis) const
  { return rep->getJointQ(s,body,axis); }
const Real& SimbodySubsystem::getJointU(const State& s, int body, int axis) const
  { return rep->getJointU(s,body,axis); }

void SimbodySubsystem::setPrescribedUdot(State& s, int body, int axis, const Real& r) const
  { return rep->setPrescribedUdot(s,body,axis,r); }

void SimbodySubsystem::clearAppliedForces(State& s) const {rep->clearAppliedForces(s);}
void SimbodySubsystem::applyGravity(State& s, const Vec3& g) const {rep->applyGravity(s,g);}
void SimbodySubsystem::applyPointForce(State& s, int body, const Vec3& stationInB, 
                                                        const Vec3& forceInG) const 
  { rep->applyPointForce(s,body,stationInB,forceInG); }

void SimbodySubsystem::applyBodyTorque(State& s, int body, const Vec3& torqueInG) const 
  { rep->applyBodyTorque(s,body,torqueInG); }
void SimbodySubsystem::applyJointForce(State& s, int body, int axis, const Real& d) const
  { rep->applyJointForce(s,body,axis,d); }

void SimbodySubsystem::enforceConfigurationConstraints(State& s) const
  { rep->enforceConfigurationConstraints(s); }
void SimbodySubsystem::enforceMotionConstraints(State& s) const
  { rep->enforceMotionConstraints(s); }

const Transform&
SimbodySubsystem::getBodyConfiguration(const State& s, int body) const
  { return rep->getBodyConfiguration(s,body); }

const SpatialVec&
SimbodySubsystem::getBodyVelocity(const State& s, int body) const
  { return rep->getBodyVelocity(s,body); }

const SpatialVec&
SimbodySubsystem::getBodyAcceleration(const State& s, int body) const
  { return rep->getBodyAcceleration(s,body); }

const VectorView SimbodySubsystem::getQDot   (const State& s) const {return rep->getQDot(s);}
const VectorView SimbodySubsystem::getUDot   (const State& s) const {return rep->getUDot(s);}
const VectorView SimbodySubsystem::getQDotDot(const State& s) const {return rep->getQDotDot(s);}

} // namespace SimTK

