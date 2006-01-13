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

#include "IVMRigidBodyTree.h"
#include "RigidBodyTree.h"

using namespace simtk;

#include <iostream>
#include <vector>

class RBTreeMap {
public:
    RBTreeMap() : body(0), frame_BR(), frame_RJ(), mass(0), com_R(0), iner_OR_R(),
                  parentIndex(badSizeTValue()), joint(0), level(-1), rbIndex(-1) { }

    RBTreeMap(const Body* b, const Frame& ref, const Frame& jInRef,
              const Real& m, const Vec3& cm_R, const MatInertia& iner_R,
              size_t pix, const Joint* j, int l)
      : body(b), frame_BR(ref), frame_RJ(jInRef), 
        mass(m), com_R(cm_R), iner_OR_R(iner_R),
        parentIndex(pix), joint(j), level(l), rbIndex(-1)
    {
    }

    const Body&  getBody()           const          {assert(body); return *body;}
    const Frame& getRefFrameInBody() const          {assert(body); return frame_BR;}
    void         setRefFrameInBody(const Frame& f)  {frame_BR=f;}
    const Frame& getJointFrameInRef() const         {assert(body); return frame_RJ;}
    void         setJointFrameInRef(const Frame& f) {frame_RJ=f;}

    const Real&       getMass()            const {assert(body); return mass;}
    const Vec3&       getCOMInRef()        const {assert(body); return com_R;}
    const MatInertia& getInertiaAboutRef() const {assert(body); return iner_OR_R;}

    size_t       getParentIndex() const {assert(parentIndex != size_t(-1)); return parentIndex;}
    const Joint& getJoint()  const {assert(joint);  return *joint;}
    int          getLevel()  const {return level;}
    void         setRBIndex(int ix) {rbIndex=ix;}
    int          getRBIndex() const {assert(rbIndex != -1); return rbIndex;}

private:
    const Body*     body;
    Frame           frame_BR;   // reference frame R meas & expr in B; default=I
    Frame           frame_RJ;   // inboard joint frame meas & expr in ref frame R
    Real            mass;
    Vec3            com_R;      // COM, meas & expr in ref frame R
    MatInertia      iner_OR_R;  // Inertia about R origin, expr. in R
    size_t          parentIndex;
    const Joint*    joint;
    int             level;
    int             rbIndex;

    // This is to avoid compiler warnings
    size_t badSizeTValue() const {
        size_t x(0);
        return (x-1); // i.e., 0xfffffff etc.
    }
};


class IVMSimbodyInterfaceRep /* : public RigidBodyMechanicsResource(?) */ {
public:
    virtual int getNBodies() const = 0;
    virtual int getNParameters() const = 0;
    virtual int getNQ() const = 0;
    virtual int getNU() const = 0;

    virtual void realizeParameters(const State&) const = 0;
    virtual void realizeConfiguration(const State&) const = 0;
    virtual void realizeMotion(const State&) const = 0;
    virtual void realizeReaction(const State&) const = 0;

    virtual void enforceConfigurationConstraints(State&) const = 0;
    virtual void enforceMotionConstraints(State&) const = 0;

    virtual const Vector& getQDot(const State&) const = 0;
    virtual Vector        calcUDot(const State& s, 
                                   const Array<SpatialVec>& bodyForces,
                                   const Vector& hingeForces) const = 0;

    virtual Frame      getBodyConfiguration(const State&, const Body& body) const = 0;
    virtual SpatialVec getBodyVelocity     (const State&, const Body& body) const = 0;
    virtual SpatialVec getBodyAcceleration (const State&, const Body& body) const = 0;



    IVMSimbodyInterfaceRep(const Multibody&);

    const Multibody& getMultibody() const {return mbs;}

    const RBTreeMap& getBodyInfo(const Body& b) const {
        // TODO: info must be stored with (or indexed by) body directly for speed
        for (size_t i=0; i<mbs2tree.size(); ++i) 
            if (b.isSameSubsystem(mbs2tree[i].getBody()))
                return mbs2tree[i];
        assert(false); // where is it?
        return *reinterpret_cast<const RBTreeMap*>(0);
    }

    const RBTreeMap& getBodyInfoByIndex(int ix) const {
        return mbs2tree[ix];
    }

    
    void setMyHandle(IVMSimbodyInterface& h) {handle=&h;}
    const IVMSimbodyInterface& getMyHandle() const {assert(handle); return *handle;}

    static RBJointType mapToRBJointType(Joint::JointType jt);
    static IVMJointType mapToIVMJointType(Joint::JointType jt);

protected:
    Multibody               mbs;  // private copy
    std::vector<RBTreeMap>  mbs2tree;

private:
    IVMSimbodyInterface*    handle;

};

class OldIVMSimbodyInterfaceRep : public IVMSimbodyInterfaceRep {
public:
    OldIVMSimbodyInterfaceRep(const Multibody& m) : IVMSimbodyInterfaceRep(m) {
        buildTree();
    }

    // virtual methods
    int getNBodies()     const {return getRigidBodyTree().getNBodies();}
    int getNParameters() const {return 0;}
    int getNQ()          const {return getRigidBodyTree().getDim();}
    int getNU()          const {return getRigidBodyTree().getDim();}

    void realizeParameters(const State& s) const { }
    void realizeConfiguration(const State& s) const;
    void realizeMotion(const State& s) const;
    void realizeReaction(const State&) const { }

    void enforceConfigurationConstraints(State&) const{
        assert(false);
    }
    void enforceMotionConstraints(State&) const{
        assert(false);
    }

    const Vector& getQDot(const State& s) const {
        return s.getU();
    }

    Vector        calcUDot(const State& s, 
                           const Array<SpatialVec>& bodyForces,
                           const Vector& hingeForces) const;

    Frame getBodyConfiguration(const State& s, const Body& body) const;

    SpatialVec getBodyVelocity     (const State&, const Body& body) const {
        assert(false);
        return SpatialVec();
    }
    SpatialVec getBodyAcceleration (const State&, const Body& body) const {
        assert(false);
        return SpatialVec();
    }

    const IVMRigidBodyTree& getRigidBodyTree() const {return tree;}
    IVMRigidBodyTree&       updRigidBodyTree()       {return tree;}
private:
    IVMRigidBodyTree tree;

    void buildTree();
};


class NewIVMSimbodyInterfaceRep : public IVMSimbodyInterfaceRep {
public:
    NewIVMSimbodyInterfaceRep(const Multibody& m) : IVMSimbodyInterfaceRep(m) {
        buildTree();
    }


    // virtual methods

    int getNBodies()     const {return getRigidBodyTree().getNBodies();}
    int getNParameters() const {return 0;}
    int getNQ()          const {return getRigidBodyTree().getDim();}
    int getNU()          const {return getRigidBodyTree().getDim();}

    void realizeParameters(const State& s) const { }
    void realizeConfiguration(const State& s) const {
        RigidBodyTree& t = const_cast<RigidBodyTree&>(tree);
        t.setPos(s.getQ());
    }
    void realizeMotion(const State& s) const {
        RigidBodyTree& t = const_cast<RigidBodyTree&>(tree);
        t.setVel(s.getU());
    }
    void realizeReaction(const State&) const { }

    void enforceConfigurationConstraints(State&) const {
        assert(false);
    }
    void enforceMotionConstraints(State&) const {
        assert(false);
    }

    const Vector& getQDot(const State& s) const {
        return s.getU();
    }

    Vector        calcUDot(const State& s, 
                           const Array<SpatialVec>& bodyForces,
                           const Vector& hingeForces) const;

    Frame getBodyConfiguration(const State& s, const Body& body) const;

    SpatialVec getBodyVelocity     (const State&, const Body& body) const {
        assert(false);
        return SpatialVec();
    }
    SpatialVec getBodyAcceleration (const State&, const Body& body) const {
        assert(false);
        return SpatialVec();
    }

    const RigidBodyTree& getRigidBodyTree() const {return tree;}
    RigidBodyTree&       updRigidBodyTree()       {return tree;}
private:
    RigidBodyTree tree;

    void buildTree();
};


#endif // IVM_SIMBODY_INTERFACE_REP_H_
