#ifndef IVM_SIMBODY_INTERFACE_REP_H_
#define IVM_SIMBODY_INTERFACE_REP_H_

/** @file
 *
 * This is a Simbody-compatible interface to the IVM-derived rigid body
 * simulation toolset.
 */

#include "simmatrix/BigMatrix.h"
#include "simbody/Simbody.h"
#include "simbody/IVMSimbodyInterface.h"

#include "RigidBodyTree.h"

using namespace simtk;

#include <iostream>
#include <vector>

struct TreeMap {
    const Body* body;
    const Body* parent;
    const Joint* joint;
    int level;

    TreeMap() : body(0), parent(0), joint(0), level(-1) { }
    TreeMap(const Body* b, const Body* p, const Joint* j, int l)
        : body(b), parent(p), joint(j), level(l) { }

    const Body&  getBody()   const {assert(body);   return *body;}
    const Body&  getParent() const {assert(parent); return *parent;}
    const Joint& getJoint()  const {assert(joint);  return *joint;}
    int          getLevel()  const {return level;}
};


class IVMSimbodyInterfaceRep /* : public RigidBodyMechanicsResource(?) */ {
public:
    IVMSimbodyInterfaceRep(const Multibody&);

    int getNBodies() const;
    int getNParameters() const;
    int getNQ() const;
    int getNU() const;

    void realizeParameters(const State&) const;
    void realizeConfiguration(const State&) const;
    void realizeMotion(const State&) const;
    void realizeReaction(const State&) const;

    void enforceConfigurationConstraints(State&) const;
    void enforceMotionConstraints(State&) const;

    const Vector& getQDot(const State&) const;

    void setMyHandle(IVMSimbodyInterface& h) {handle=&h;}
    const IVMSimbodyInterface& getMyHandle() const {assert(handle); return *handle;}

    const Frame&         getBodyConfiguration(const State&, int body) const;
    const SpatialVector& getBodyVelocity     (const State&, int body) const;
    const SpatialVector& getBodyAcceleration (const State&, int body) const;
private:
    IVMSimbodyInterface* handle;

    Multibody               mbs;  // private copy
    std::vector<TreeMap>    mbs2tree;
    RigidBodyTree           tree;

};



#endif // IVM_SIMBODY_INTERFACE_REP_H_
