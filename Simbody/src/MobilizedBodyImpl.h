#ifndef SimTK_SIMBODY_MOBILIZED_BODY_IMPL_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

/**@file
 *
 * Private implementation of Body and MobilizedBody, and their included 
 * subclasses which represent the built-in body and mobilizer types.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Body.h"
#include "simbody/internal/MobilizedBody.h"
#include "simbody/internal/MobilizedBody_BuiltIns.h"

#include "SimbodyMatterSubsystemRep.h"
#include "BodyRep.h"
#include "MotionImpl.h"

class RigidBodyNode;

namespace SimTK {

    /////////////////////////
    // MOBILIZED BODY REPS //
    /////////////////////////

// This is what a MobilizedBody handle points to.
class MobilizedBodyImpl : public PIMPLImplementation<MobilizedBody,MobilizedBodyImpl> {
public:
    explicit MobilizedBodyImpl(MobilizedBody::Direction d) 
    :   defaultLockLevel(Motion::NoLevel), 
        reversed(d==MobilizedBody::Reverse), 
        myMatterSubsystemRep(0), myLevel(-1), myRBnode(0), hasChildren(false) 
    {
    }

    void setDirection(MobilizedBody::Direction d) {
        const bool wantReversed = (d==MobilizedBody::Reverse);
        if (wantReversed != reversed) {
            invalidateTopologyCache();
            reversed = wantReversed;
        }
    }

    MobilizedBodyImpl(const MobilizedBodyImpl& src) {
        *this = src;
        myRBnode = 0;
    }


    void lock(State& state, Motion::Level level) const;
    void lockAt(State& state, int n, const Real* value, Motion::Level level) const;
    void unlock(State& state) const;
    Motion::Level getLockLevel(const State& state) const;
    Vector getLockValueAsVector(const State& state) const;
    void lockByDefault(Motion::Level level)
    {   invalidateTopologyCache(); defaultLockLevel=level; }
    Motion::Level getLockByDefaultLevel() const
    {   return defaultLockLevel; }

    // The matter subsystem must issue these MobilizedBody realize() calls in base-to-tip
    // order, because these methods are allowed to assume that their parent (and 
    // ancestors) have already been realized.

    // eventually calls realizeTopologyVirtual()
    const RigidBodyNode& realizeTopology(State& s, UIndex& nxtU, USquaredIndex& nxtUSq, QIndex& nxtQ) const;

    void realizeModel   (State&)       const; // eventually calls realizeModelVirtual()       
    void realizeInstance(const SBStateDigest&) const; // eventually calls realizeInstanceVirtual() 
    void realizeTime    (const SBStateDigest&) const; // eventually calls realizeTimeVirtual() 
    void realizePosition(const SBStateDigest&) const; // eventually calls realizePositionVirtual() 
    void realizeVelocity(const SBStateDigest&) const; // eventually calls realizeVelocityVirtual() 
    void realizeDynamics(const SBStateDigest&) const; // eventually calls realizeDynamicsVirtual() 
    void realizeAcceleration
                        (const SBStateDigest&) const; // eventually calls realizeAccelerationVirtual() 
    void realizeReport  (const SBStateDigest&) const; // eventually calls realizeReportVirtual() 

    virtual ~MobilizedBodyImpl() {
        delete myRBnode;
    }
    virtual MobilizedBodyImpl* clone() const = 0;
    
    // This creates a rigid body node using the appropriate mobilizer type.
    // Called during realizeTopology().
    virtual RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const = 0;

    virtual void realizeTopologyVirtual     (State&)        const {}
    virtual void realizeModelVirtual        (State&)        const {}
    virtual void realizeInstanceVirtual     (const State&)  const {}
    virtual void realizeTimeVirtual         (const State&)  const {}
    virtual void realizePositionVirtual     (const State&)  const {}
    virtual void realizeVelocityVirtual     (const State&)  const {}
    virtual void realizeDynamicsVirtual     (const State&)  const {}
    virtual void realizeAccelerationVirtual (const State&)  const {}
    virtual void realizeReportVirtual       (const State&)  const {}

    // Copy out nq default values for q, beginning at the indicated address.
    // The concrete class should assert if nq is not a reasonable
    // number for the kind of mobilizer; there is a bug somewhere in that case. 
    // This routine shouldn't be called directly -- call copyOutDefaultQ() below 
    // instead which has a nicer interface and does some error checking.
    virtual void copyOutDefaultQImpl(int nq, Real* q) const = 0;

    virtual void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
    {
        if (stage != Stage::Instance || !getMyMatterSubsystemRep().getShowDefaultGeometry())
            return;
        const Real scale = 1;
        DecorativeFrame axes(scale/2);
        axes.setLineThickness(2);
        axes.setBodyId(myMobilizedBodyIndex);
        geom.push_back(axes); // the body frame

        // Display the inboard joint frame (at half size), unless it is the
        // same as the body frame. Then find the corresponding frame on the
        // parent and display that in this body's color.
        if (myMobilizedBodyIndex != 0) {
            const Real pscale = 1;
            const Transform& X_BM = getDefaultOutboardFrame(); // TODO: get from state
            if (X_BM.p() != Vec3(0) || X_BM.R() != Mat33(1)) {
                DecorativeFrame frameOnChild(scale/4);
                frameOnChild.setBodyId(myMobilizedBodyIndex);
                frameOnChild.setColor(Red);
                frameOnChild.setTransform(X_BM);
                geom.push_back(frameOnChild);
                if (X_BM.p() != Vec3(0))
                    geom.push_back(DecorativeLine(Vec3(0), X_BM.p()).setBodyId(myMobilizedBodyIndex));
            }
            const Transform& X_PF = getDefaultInboardFrame(); // TODO: from state
            DecorativeFrame frameOnParent(pscale*Real(0.3)); // slightly larger than child
            frameOnParent.setBodyId(myParentIndex);
            frameOnParent.setColor(Blue);
            frameOnParent.setTransform(X_PF);
            geom.push_back(frameOnParent);
            if (X_PF.p() != Vec3(0))
                geom.push_back(DecorativeLine(Vec3(0),X_PF.p()).setBodyId(myParentIndex));
        }

        // Put a little purple wireframe sphere at the COM, and add a line from 
        // body origin to the com.

        DecorativeSphere com(scale*Real(.05));
        com.setBodyId(myMobilizedBodyIndex);
        com.setColor(Purple).setRepresentation(DecorativeGeometry::DrawPoints);
        const Vec3& comPos_B = theBody.getDefaultRigidBodyMassProperties().getMassCenter(); // TODO: from state
        com.setTransform(comPos_B);
        geom.push_back(com);
        if (comPos_B != Vec3(0))
            geom.push_back(DecorativeLine(Vec3(0), comPos_B).setBodyId(myMobilizedBodyIndex));
    }

    void calcDecorativeGeometryAndAppend
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
    {
        // We know how to deal with the topological (construction) geometry
        // here. For bodies we can just draw it at topology stage. For mobilizers,
        // we might not know the placement of the mobilizer frames on the parent
        // and child bodies until Instance stage, at which point we can transform
        // and then draw the topological geometry.
        if (stage == Stage::Topology) {
            appendTopologicalBodyGeometry(geom);
        } else if (stage == Stage::Instance) {
            const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
            const Transform& X_PF = getInboardFrame(s);
            const Transform& X_BM  = getOutboardFrame(s);
            appendTopologicalMobilizerGeometry(X_BM, X_PF, geom);
        }

        // Let the individual mobilizer deal with any complicated stuff.
        calcDecorativeGeometryAndAppendImpl(s,stage,geom);
    }

    int addOutboardDecoration(const Transform& X_MD, const DecorativeGeometry& g) {
        const int nxt = (int)outboardGeometry.size();
        outboardGeometry.push_back(g); // make a new copy
        // Combine the placement frame and the transform already in the geometry
        // so we end up with geometry expressed directly in the M frame.
        outboardGeometry.back().setTransform(X_MD*g.getTransform());
        // Record the assigned ordinal (not in the same ordinal space as
        // body decorations).
        outboardGeometry.back().setIndexOnBody(nxt);
        return nxt;
    }
    int addInboardDecoration(const Transform& X_FD, const DecorativeGeometry& g) {
        const int nxt = (int)inboardGeometry.size();
        inboardGeometry.push_back(g); // make a new copy
        // Combine the placement frame and the transform already in the geometry
        // so we end up with geometry expressed directly in the F frame.
        inboardGeometry.back().setTransform(X_FD*g.getTransform());
        // Record the assigned ordinal (not in the same ordinal space as
        // body decorations).
        inboardGeometry.back().setIndexOnBody(nxt);
        return nxt;
    }


    void findMobilizerQs(const State&, QIndex& qStart, int& nq) const;
    void findMobilizerUs(const State&, UIndex& uStart, int& nu) const;


    Motion::Method getQMotionMethod(const State&) const;
    Motion::Method getUMotionMethod(const State&) const;
    Motion::Method getUDotMotionMethod(const State&) const;

    // Given the Model stage state variables (if any are relevant), each
    // concrete mobilizer should return the number of q's and u's it expects
    // to be allotted.
    //virtual void getStateNeeds(const State& s, int& nq, int& nu) const = 0;

    // Given the Model stage variable values in the passed-in State (if
    // any of them are relevant) copy out the appropriate default values
    // to the appropriate slot in qDefault.
    void copyOutDefaultQ(const State& s, Vector& qDefault) const;

    // Generic State-access routines (that is, those which can be handled in
    // the MobilizedBody base class).


    // Invalidate Stage::Position.
    void setQToFitTransform(State& s, const Transform& X_FM) const;
    void setQToFitRotation(State& s, const Rotation& R_FM) const;
    void setQToFitTranslation(State& s, const Vec3& p_FM) const;

    // Invalidate Stage::Velocity.
    void setUToFitVelocity(State& s, const SpatialVec& V_FM) const;
    void setUToFitAngularVelocity(State& s, const Vec3& w_FM) const;
    void setUToFitLinearVelocity(State& s, const Vec3& v_FM)  const;

    const MassProperties& getBodyMassProperties(const State& s) const {
        // TODO: these should come from the state if the body has variable mass props
        const SBInstanceVars& iv = getMyMatterSubsystemRep().getInstanceVars(s);
        return getMyRigidBodyNode().getMassProperties_OB_B();
    }

    const SpatialInertia& getBodySpatialInertiaInGround(const State& s) const {
        const SBTreePositionCache& tpc = getMyMatterSubsystemRep().getTreePositionCache(s);
        return getMyRigidBodyNode().getMk_G(tpc);
    }

    const Transform& getInboardFrame (const State& s) const {
        // TODO: these should come from the state if the mobilizer has variable frames
        const SBInstanceVars& iv = getMyMatterSubsystemRep().getInstanceVars(s);
        return getMyRigidBodyNode().getX_PF();
    }
    const Transform& getOutboardFrame(const State& s) const {
        // TODO: these should come from the state if the mobilizer has variable frames
        const SBInstanceVars& iv = getMyMatterSubsystemRep().getInstanceVars(s);
        return getMyRigidBodyNode().getX_BM();
    }

    void setInboardFrame (State& s, const Transform& X_PF) const {
        assert(!"setInboardFrame(s) not implemented yet");
    }
    void setOutboardFrame(State& s, const Transform& X_BM) const {
        assert(!"setOutboardFrame(s) not implemented yet");
    }

    const Transform& getBodyTransform(const State& s) const {
        const SBTreePositionCache& pc = getMyMatterSubsystemRep().getTreePositionCache(s);
        return getMyRigidBodyNode().getX_GB(pc);
    }
    const SpatialVec& getBodyVelocity(const State& s) const {
        const SBTreeVelocityCache& vc = getMyMatterSubsystemRep().getTreeVelocityCache(s);
        return getMyRigidBodyNode().getV_GB(vc);
    }
    const SpatialVec& getBodyAcceleration(const State& s) const {
        const SBTreeAccelerationCache& ac = getMyMatterSubsystemRep().getTreeAccelerationCache(s);
        return getMyRigidBodyNode().getA_GB(ac);
    }

    const Transform& getMobilizerTransform(const State& s) const {
        const SBTreePositionCache& pc = getMyMatterSubsystemRep().getTreePositionCache(s);
        return getMyRigidBodyNode().getX_FM(pc);
    }
    const SpatialVec& getMobilizerVelocity(const State& s) const {
        const SBTreeVelocityCache& vc = getMyMatterSubsystemRep().getTreeVelocityCache(s);
        return getMyRigidBodyNode().getV_FM(vc);
    }

    SpatialVec getHCol(const State& s, MobilizerUIndex ux) const {
        const SBTreePositionCache& pc = getMyMatterSubsystemRep().getTreePositionCache(s);
        return getMyRigidBodyNode().getHCol(pc,ux);
    }

    SpatialVec getH_FMCol(const State& s, MobilizerUIndex ux) const {
        const SBTreePositionCache& pc = getMyMatterSubsystemRep().getTreePositionCache(s);
        return getMyRigidBodyNode().getH_FMCol(pc,ux);
    }

    void invalidateTopologyCache() const {
        delete myRBnode; myRBnode=0;
        if (myMatterSubsystemRep)
            myMatterSubsystemRep->invalidateSubsystemTopologyCache();
    }

    // This might get called *during* realizeTopology() so just make sure there is 
    // a node here without checking whether we're done with realizeTopology().
    const RigidBodyNode& getMyRigidBodyNode() const {
        SimTK_ASSERT(myRBnode && myMatterSubsystemRep,
          "An operation on a MobilizedBody was illegal because realizeTopology() has "
          "not been performed on the containing Subsystem since the last topological change."
        );
        return *myRBnode;
    }

    const Body& getBody() const {return theBody;}
    const Transform& getDefaultInboardFrame()  const {return defaultInboardFrame;}
    const Transform& getDefaultOutboardFrame() const {return defaultOutboardFrame;}
    const MassProperties& getDefaultRigidBodyMassProperties() const {
        return theBody.getDefaultRigidBodyMassProperties();
    }

    bool isReversed() const {return reversed;}

    void adoptMotion(Motion& ownerHandle) {
        SimTK_ERRCHK(!hasMotion(), "MobilizedBody::adoptMotion()",
            "This MobilizedBody already has a Motion object associated with it.\n"
            "Use clearMotion() first to replace an existing Motion object.");
        ownerHandle.disown(motion); // transfer ownership to handle "motion"
        invalidateTopologyCache();

        // Tell the Motion that it belongs to this MobilizedBody.
        motion.updImpl().setMobilizedBodyImpl(this);
    }

    void clearMotion() {
        motion.clearHandle();
        invalidateTopologyCache();
    }

    bool hasMotion() const {return !motion.isEmptyHandle();}

    const Motion& getMotion() const {
        SimTK_ERRCHK(!motion.isEmptyHandle(),
            "MobilizedBody::getMotion()",
            "There is no Motion object associated with this MobilizedBody.");
        return motion;
    }

    const SimbodyMatterSubsystem& getMySimbodyMatterSubsystem() const {
        return getMyMatterSubsystemRep().getMySimbodyMatterSubsystemHandle();
    }

    const SimbodyMatterSubsystemRep& getMyMatterSubsystemRep() const {
        SimTK_ASSERT(myMatterSubsystemRep,
            "An operation was illegal because a MobilizedBody was not in a Subsystem.");
        return *myMatterSubsystemRep;
    }
    SimbodyMatterSubsystemRep& updMyMatterSubsystemRep() {
        SimTK_ASSERT(myMatterSubsystemRep,
            "An operation was illegal because a MobilizedBody was not in a Subsystem.");
        return *myMatterSubsystemRep;
    }

    const SBModelPerMobodInfo& getMyModelInfo(const State& s) const {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        const SBModelCache& mc = matterRep.getModelCache(s);
        const SBModelPerMobodInfo& info = 
            mc.getMobodModelInfo(myMobilizedBodyIndex);
        return info;
    }

    const SBInstancePerMobodInfo& getMyInstanceInfo(const State& s) const {
        const SimbodyMatterSubsystemRep& matterRep = getMyMatterSubsystemRep();
        const SBInstanceCache& ic = matterRep.getInstanceCache(s);
        const SBInstancePerMobodInfo& info = 
            ic.getMobodInstanceInfo(myMobilizedBodyIndex);
        return info;
    }

    const MobilizedBody& getMyHandle() const {
        const MobilizedBody& mobod = 
            getMyMatterSubsystemRep().getMobilizedBody(getMyMobilizedBodyIndex());
        SimTK_ASSERT(&mobod.getImpl() == this,
            "A MobilizedBodyImpl's handle didn't refer back to the same Impl!");
        return mobod;
    }

    MobilizedBody& updMyHandle() {
        MobilizedBody& mobod = 
            updMyMatterSubsystemRep().updMobilizedBody(getMyMobilizedBodyIndex());
        SimTK_ASSERT(&mobod.getImpl() == this,
            "A MobilizedBodyImpl's handle didn't refer back to the same Impl!");
        return mobod;
    }

    MobilizedBodyIndex getMyMobilizedBodyIndex() const {
        assert(myMobilizedBodyIndex.isValid());
        return myMobilizedBodyIndex;
    }

    MobilizedBodyIndex getMyParentMobilizedBodyIndex() const {
        assert(myParentIndex.isValid());
        return myParentIndex;
    }

    MobilizedBodyIndex getMyBaseBodyMobilizedBodyIndex() const {
        assert(myBaseBodyIndex.isValid());
        return myBaseBodyIndex;
    }

    int getMyLevel() const {
        assert(myLevel >= 0);
        return myLevel;
    }

    bool isInSubsystem() const {return myMatterSubsystemRep != 0;}

    void setMyMatterSubsystem(SimbodyMatterSubsystem& matter,
                              MobilizedBodyIndex  parentIndex,
                              MobilizedBodyIndex  index)
    {
        // If the subsystem is already set it must be the same one.
        assert(!myMatterSubsystemRep || myMatterSubsystemRep==&matter.getRep());
        myMatterSubsystemRep = &matter.updRep();

        assert(index.isValid());
        assert(parentIndex.isValid() || index==GroundIndex);

        myParentIndex           = parentIndex; // invalid if this is Ground
        myMobilizedBodyIndex    = index;

        if (index != GroundIndex) {
            MobilizedBody& parent = matter.updMobilizedBody(parentIndex);
            myLevel = parent.getLevelInMultibodyTree() + 1;
            myBaseBodyIndex = (myLevel == 1 ? myMobilizedBodyIndex 
                                         : parent.getBaseMobilizedBody().getMobilizedBodyIndex());
            parent.updImpl().hasChildren = true;
        } else {
            myLevel = 0;
            myBaseBodyIndex = GroundIndex;
        }
    }

private:
    // Body topological geometry is defined with respect to the body frame so we
    // can draw it right away.
    void appendTopologicalBodyGeometry(Array_<DecorativeGeometry>& geom) const {
        getBody().getRep().appendDecorativeGeometry(getMyMobilizedBodyIndex(), geom);
    }

    // Mobilizer topological geometry is defined with respect to the M (outboard, child)
    // frame and the F (inboard, parent) frame. The placement of those frames with
    // respect to the body frame can be Instance variables, so we can't draw this geometry
    // until Instance stage. At that point we can find M and F, so they are passed in
    // here.
    void appendTopologicalMobilizerGeometry(const Transform& X_BM, const Transform& X_PF,
                                            Array_<DecorativeGeometry>& geom) const
    {
        for (int i=0; i<(int)outboardGeometry.size(); ++i) {
            geom.push_back(outboardGeometry[i]);
            geom.back().setBodyId(getMyMobilizedBodyIndex())
                       .setTransform(X_BM*outboardGeometry[i].getTransform());
        }
        for (int i=0; i<(int)inboardGeometry.size(); ++i) {
            geom.push_back(inboardGeometry[i]);
            geom.back().setBodyId(getMyParentMobilizedBodyIndex())
                       .setTransform(X_PF*inboardGeometry[i].getTransform());
        }
    }

private:
    friend class MobilizedBody;

        // TOPOLOGY "STATE"

    // Base class topological properties. Derived MobilizedBodies may
    // have further topological properties. Whenever these change, be
    // sure to set topologyRealized=false!

    // Body represents the mass structure, mass properties, and possibly some
    // decorative geometry for the body associated with this (body,mobilizer)
    // pair.
    Body            theBody;
    Transform       defaultInboardFrame;  // default for F (in Parent frame)
    Transform       defaultOutboardFrame; // default for M (in Body frame)
    Motion::Level   defaultLockLevel;

    // A Motion object, if present, defines how this mobilizer's motion is
    // to be calculated. Otherwise, the motion is determined dynamically
    // as a result of forces and constraints. A Motion always prescribes
    // \e all of the mobilizer's generalized accelerations udot (index 1); 
    // may also some of the prescribed generalized speeds u (index 2); and 
    // for speeds that are prescribed it may also prescribe the corresponding
    // generalized coordinates q (index 3).
    Motion          motion;

    bool        reversed; // is the mobilizer defined from M to F?

    Array_<DecorativeGeometry> outboardGeometry;
    Array_<DecorativeGeometry> inboardGeometry;

    // These data members are filled in once the MobilizedBody is added to
    // a MatterSubsystem. Note that this pointer is just a reference to
    // the Subsystem which owns this MobilizedBody -- we don't own the
    // heap space here so don't delete it!
    SimbodyMatterSubsystemRep* myMatterSubsystemRep;
    MobilizedBodyIndex  myMobilizedBodyIndex; // index within the subsystem
    MobilizedBodyIndex  myParentIndex;   // Invalid if this body is Ground, otherwise the parent's index
    MobilizedBodyIndex  myBaseBodyIndex; // GroundIndex if this is ground, otherwise a level-1 BodyIndex
    int                 myLevel;      // Distance from ground in multibody graph (0 if this is ground,
                                      //   1 if a base body, 2 if parent is a base body, etc.)

        // TOPOLOGY "CACHE"

    // A RigidBodyNode is created for each MobilizedBody during realizeTopology().
    // Think of it as the *computational* form of the MobilizedBody; the MobilizedBody
    // itself exists to make a nice API for the programmer. In fact is is possible
    // for several different MobilizedBody classes to use the same underlying
    // computational form while providing a different API.
    //
    // This is a pointer to an abstract object whose heap space is *owned* by
    // this MobilizedBody. Be sure to delete it upon destruction and whenever
    // topology is re-realized.
    mutable RigidBodyNode* myRBnode;

protected:
    // Keep track of whether this body has any children.
    bool hasChildren;
};

class MobilizedBody::PinImpl : public MobilizedBodyImpl {
public:
    explicit PinImpl(Direction d) : MobilizedBodyImpl(d), defaultQ(0) { }
    PinImpl* clone() const override { return new PinImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==1, 
            "MobilizedBody::PinImpl::copyOutDefaultQImpl(): wrong number of q's");
        *q = defaultQ;
    }

    SimTK_DOWNCAST(PinImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Pin;
    Real defaultQ;  // the angle in radians
};


class MobilizedBody::SliderImpl : public MobilizedBodyImpl {
public:
    explicit SliderImpl(Direction d) : MobilizedBodyImpl(d), defaultQ(0) { }
    SliderImpl* clone() const override { return new SliderImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==1, 
            "MobilizedBody::SliderImpl::copyOutDefaultQImpl(): wrong number of q's");
        *q = defaultQ;
    }

    SimTK_DOWNCAST(SliderImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Slider;
    Real defaultQ;  // the displacement
};

class MobilizedBody::UniversalImpl : public MobilizedBodyImpl {
public:
    explicit UniversalImpl(Direction d) : MobilizedBodyImpl(d), defaultQ(0) { }
    UniversalImpl* clone() const override { return new UniversalImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==2, 
            "MobilizedBody::UniversalImpl::copyOutDefaultQImpl(): wrong number of q's");
        Vec2::updAs(q) = defaultQ;
    }

    SimTK_DOWNCAST(UniversalImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Universal;
    Vec2 defaultQ;  // the two angles
};

class MobilizedBody::CylinderImpl : public MobilizedBodyImpl {
public:
    explicit CylinderImpl(Direction d) : MobilizedBodyImpl(d), defaultQ(0) { }
    CylinderImpl* clone() const override { return new CylinderImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==2, 
            "MobilizedBody::CylinderImpl::copyOutDefaultQImpl(): wrong number of q's");
        Vec2::updAs(q) = defaultQ;
    }

    SimTK_DOWNCAST(CylinderImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Cylinder;
    Vec2 defaultQ;  // the angle in radians followed by the displacement
};

class MobilizedBody::BendStretchImpl : public MobilizedBodyImpl {
public:
    explicit BendStretchImpl(Direction d) : MobilizedBodyImpl(d), defaultQ(0) { }
    BendStretchImpl* clone() const override { return new BendStretchImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==2, 
            "MobilizedBody::BendStretchImpl::copyOutDefaultQImpl(): wrong number of q's");
        Vec2::updAs(q) = defaultQ;
    }

    SimTK_DOWNCAST(BendStretchImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::BendStretch;
    Vec2 defaultQ;  // the angle in radians followed by the displacement
};

class MobilizedBody::PlanarImpl : public MobilizedBodyImpl {
public:
    explicit PlanarImpl(Direction d) : MobilizedBodyImpl(d), defaultQ(0) { }
    PlanarImpl* clone() const override { return new PlanarImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==3, 
            "MobilizedBody::PlanarImpl::copyOutDefaultQImpl(): wrong number of q's");
        Vec3::updAs(q) = defaultQ;
    }

    SimTK_DOWNCAST(PlanarImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Planar;
    Vec3 defaultQ;  // the angle in radians followed by the two displacements
};

class MobilizedBody::SphericalCoordsImpl : public MobilizedBodyImpl {
public:
    explicit SphericalCoordsImpl(Direction d) 
    :   MobilizedBodyImpl(d), az0(0), ze0(0), axisT(ZAxis), 
        negAz(false), negZe(false), negT(false), 
        defaultQ(0) {}

    SphericalCoordsImpl(Real az0, bool negAz, 
                        Real ze0, bool negZe,
                        CoordinateAxis axisT, bool negT,
                        Direction d) 
    :   MobilizedBodyImpl(d), az0(az0), ze0(ze0), axisT(axisT), 
        negAz(negAz), negZe(negZe), negT(negT), 
        defaultQ(0) {}

    SphericalCoordsImpl* clone() const override { return new SphericalCoordsImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==3, 
            "MobilizedBody::SphericalCoordsImpl::copyOutDefaultQImpl(): wrong number of q's");
        Vec3::updAs(q) = defaultQ;
    }

    SimTK_DOWNCAST(SphericalCoordsImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::SphericalCoords;

    Real              az0, ze0;               // angle offsets
    CoordinateAxis    axisT;                  // translation axis (X or Z)
    bool              negAz, negZe, negT;     // true if we should negate the corresponding coordinate

    Vec3 defaultQ;  // two angles in radians followed by a displacement
};

class MobilizedBody::GimbalImpl : public MobilizedBodyImpl {
public:
    explicit GimbalImpl(Direction d) 
    :   MobilizedBodyImpl(d), defaultRadius(Real(0.1)), defaultQ(0) { }
    GimbalImpl* clone() const override { return new GimbalImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==3, 
            "MobilizedBody::GimbalImpl::copyOutDefaultQImpl(): wrong number of q's");
        Vec3::updAs(q) = defaultQ;
    }

    void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const override;

    void setDefaultRadius(Real r) {
        assert(r>0);
        invalidateTopologyCache();
        defaultRadius=r;
    }
    Real getDefaultRadius() const {return defaultRadius;}

    SimTK_DOWNCAST(GimbalImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Gimbal;
    Real defaultRadius;   // used for visualization only
    Vec3 defaultQ;  // the three angles in radians
};

class MobilizedBody::BushingImpl : public MobilizedBodyImpl {
public:
    explicit BushingImpl(Direction d) 
    :   MobilizedBodyImpl(d), defaultQ(0) { }
    BushingImpl* clone() const override { return new BushingImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==6, 
            "MobilizedBody::BushingImpl::copyOutDefaultQImpl(): wrong number of q's");
        Vec6::updAs(q) = defaultQ;
    }

    SimTK_DOWNCAST(BushingImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Bushing;
    Vec6 defaultQ;  // 3 angles in radians, then p_FM
};

class MobilizedBody::BallImpl : public MobilizedBodyImpl {
public:
    explicit BallImpl(Direction d) 
    :   MobilizedBodyImpl(d), defaultRadius(Real(0.1)), defaultQ() {} // (1,0,0,0), the identity rotation
    BallImpl* clone() const override { return new BallImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==4||nq==3, 
            "MobilizedBody::BallImpl::copyOutDefaultQImpl(): wrong number of q's");
        if (nq==4)
            Vec4::updAs(q) = defaultQ.asVec4();
        else
            Vec3::updAs(q) = Rotation(defaultQ).convertRotationToBodyFixedXYZ();
    }

    void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const override;

    void setDefaultRadius(Real r) {
        assert(r>0);
        invalidateTopologyCache();
        defaultRadius=r;
    }
    Real getDefaultRadius() const {return defaultRadius;}

    SimTK_DOWNCAST(BallImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Ball;
    Real defaultRadius;   // used for visualization only
    Quaternion defaultQ;  // the default orientation
};

class MobilizedBody::EllipsoidImpl : public MobilizedBodyImpl {
public:
    explicit EllipsoidImpl(Direction d) 
    :   MobilizedBodyImpl(d), defaultRadii(Real(0.5),Real(1./3.),Real(0.25)), 
        defaultQ() { } // default is (1,0,0,0), the identity rotation
    EllipsoidImpl* clone() const override { return new EllipsoidImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==4||nq==3, 
            "MobilizedBody::EllipsoidImpl::copyOutDefaultQImpl(): wrong number of q's");
        if (nq==4)
            Vec4::updAs(q) = defaultQ.asVec4();
        else
            Vec3::updAs(q) = Rotation(defaultQ).convertRotationToBodyFixedXYZ();
    }

    void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const override;

    void setDefaultRadii(const Vec3& r) {
        assert(r[0]>0 && r[1]>0 && r[2]>0);
        invalidateTopologyCache();
        defaultRadii=r;
    }
    const Vec3& getDefaultRadii() const {return defaultRadii;}

    SimTK_DOWNCAST(EllipsoidImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Ellipsoid;
    Vec3 defaultRadii;    // used for visualization only
    Quaternion defaultQ;  // the default orientation
};

class MobilizedBody::TranslationImpl : public MobilizedBodyImpl {
public:
    explicit TranslationImpl(Direction d) : MobilizedBodyImpl(d), defaultQ(0) { }
    TranslationImpl* clone() const override { return new TranslationImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==3, 
            "MobilizedBody::TranslationImpl::copyOutDefaultQImpl(): wrong number of q's");
        Vec3::updAs(q) = defaultQ;
    }

    SimTK_DOWNCAST(TranslationImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Translation;
    Vec3 defaultQ;  // the three displacements
};

class MobilizedBody::FreeImpl : public MobilizedBodyImpl {
public:
    explicit FreeImpl(Direction d) : MobilizedBodyImpl(d), defaultQOrientation(), defaultQTranslation(0) { }
    FreeImpl* clone() const override { return new FreeImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==7||nq==6, 
            "MobilizedBody::FreeImpl::copyOutDefaultQImpl(): wrong number of q's");
        if (nq==7) {
            Vec4::updAs(q)   = defaultQOrientation.asVec4();
            Vec3::updAs(q+4) = defaultQTranslation;
        } else {
            Vec3::updAs(q)   = Rotation(defaultQOrientation).convertRotationToBodyFixedXYZ();
            Vec3::updAs(q+3) = defaultQTranslation;
        }
    }

    SimTK_DOWNCAST(FreeImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Free;
    Quaternion defaultQOrientation;
    Vec3       defaultQTranslation;
};

class MobilizedBody::LineOrientationImpl : public MobilizedBodyImpl {
public:
    explicit LineOrientationImpl(Direction d) : MobilizedBodyImpl(d), defaultQ() { } // 1,0,0,0
    LineOrientationImpl* clone() const override { return new LineOrientationImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==4||nq==3, 
            "MobilizedBody::LineOrientationImpl::copyOutDefaultQImpl(): wrong number of q's");
        if (nq==4)
            Vec4::updAs(q) = defaultQ.asVec4();
        else
            Vec3::updAs(q) = Rotation(defaultQ).convertRotationToBodyFixedXYZ();
    }

    SimTK_DOWNCAST(LineOrientationImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::LineOrientation;
    Quaternion defaultQ;
};

class MobilizedBody::FreeLineImpl : public MobilizedBodyImpl {
public:
    explicit FreeLineImpl(Direction d) : MobilizedBodyImpl(d), defaultQOrientation(), defaultQTranslation(0) { }
    FreeLineImpl* clone() const override { return new FreeLineImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==7||nq==6, 
            "MobilizedBody::FreeLineImpl::copyOutDefaultQImpl(): wrong number of q's");
        if (nq==7) {
            Vec4::updAs(q)   = defaultQOrientation.asVec4();
            Vec3::updAs(q+4) = defaultQTranslation;
        } else {
            Vec3::updAs(q)   = Rotation(defaultQOrientation).convertRotationToBodyFixedXYZ();
            Vec3::updAs(q+3) = defaultQTranslation;
        }
    }

    SimTK_DOWNCAST(FreeLineImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::FreeLine;
    Quaternion defaultQOrientation;
    Vec3       defaultQTranslation;
};

class MobilizedBody::WeldImpl : public MobilizedBodyImpl {
public:
    explicit WeldImpl(Direction d) : MobilizedBodyImpl(d) { }
    WeldImpl* clone() const override { return new WeldImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==0, 
            "MobilizedBody::WeldImpl::copyOutDefaultQImpl(): wrong number of q's");
    }

    SimTK_DOWNCAST(WeldImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Weld;
    // no q's
};

class MobilizedBody::GroundImpl : public MobilizedBodyImpl {
public:
    GroundImpl() : MobilizedBodyImpl(MobilizedBody::Forward) { }
    GroundImpl* clone() const override { return new GroundImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==0, 
            "MobilizedBody::GroundImpl::copyOutDefaultQImpl(): wrong number of q's");
    }

    SimTK_DOWNCAST(GroundImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Ground;
    // no q's
};

class MobilizedBody::ScrewImpl : public MobilizedBodyImpl {
public:
    ScrewImpl(Real p, Direction d) : MobilizedBodyImpl(d), defaultPitch(p), defaultQ(0) { }

    Real getDefaultPitch() const {return defaultPitch;}
    void setDefaultPitch(Real p) {
        invalidateTopologyCache();
        defaultPitch=p;
    }

    ScrewImpl* clone() const override { return new ScrewImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;

    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==1, 
            "MobilizedBody::ScrewImpl::copyOutDefaultQImpl(): wrong number of q's");
        *q = defaultQ;
    }

    SimTK_DOWNCAST(ScrewImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Screw;
    Real defaultPitch;
    Real defaultQ; // the angle in radians
};

/////////////////////////////////////////////////
// MOBILIZED BODY::CUSTOM::IMPLEMENTATION IMPL //
/////////////////////////////////////////////////

// This class exists primarily to allow the Custom::Implementation class to keep
// a pointer to its handle class's CustomImpl class which is derived from MobilizedBodyImpl
// which has all the goodies that are needed for defining a MobilizedBody.
//
// At first this class is the owner of the CustomImpl. Then when this is put in a
// Custom handle, that handle takes over ownership of the CustomImpl and the 
// CustomImpl takes over ownership of this ImplementationImpl object.
class MobilizedBody::Custom::ImplementationImpl : public PIMPLImplementation<Implementation, ImplementationImpl> {
public:
    // no default constructor
    explicit ImplementationImpl(CustomImpl* customImpl, int nu, int nq, int nAngles) : isOwner(true), builtInImpl(customImpl), nu(nu), nq(nq), nAngles(nAngles) { }
    inline ~ImplementationImpl(); // see below -- have to wait for CustomImpl's definition
    
    // Copying one of these just gives us a new one with a NULL CustomImpl pointer.
    ImplementationImpl(const ImplementationImpl& src) : isOwner(false), builtInImpl(0), nu(src.nu), nq(src.nq), nAngles(src.nAngles) { }
    
    ImplementationImpl* clone() const {return new ImplementationImpl(*this);}
    
    bool isOwnerOfCustomImpl() const {return builtInImpl && isOwner;}
    CustomImpl* removeOwnershipOfCustomImpl() {
        assert(isOwnerOfCustomImpl()); 
        isOwner=false; 
        return builtInImpl;
    }
    
    void setReferenceToCustomImpl(CustomImpl* cimpl) {
        assert(!builtInImpl); // you can only do this once
        isOwner=false;
        builtInImpl = cimpl;
    }
    
    bool hasCustomImpl() const {return builtInImpl != 0;}
    
    const CustomImpl& getCustomImpl() const {
        assert(builtInImpl);
        return *builtInImpl;
    }
    CustomImpl& updCustomImpl() {
        assert(builtInImpl);
        return *builtInImpl;
    }
    
    int getNU() const {
        return nu;
    }
    
    int getNQ() const {
        return nq;
    }
    
    int getNumAngles() const {
        return nAngles;
    }

private:
    bool isOwner;
    CustomImpl* builtInImpl; // just a reference; not owned
    int nu, nq, nAngles;

    // suppress assignment
    ImplementationImpl& operator=(const ImplementationImpl&);
};

/////////////////////////////////
// MOBILIZED BODY::CUSTOM IMPL //
/////////////////////////////////

class MobilizedBody::CustomImpl : public MobilizedBodyImpl {
public:
    explicit CustomImpl(Direction d) : MobilizedBodyImpl(d), implementation(0) { }
    
    void takeOwnershipOfImplementation(Custom::Implementation* userImpl);
    
    explicit CustomImpl(Custom::Implementation* userImpl, Direction d)
    :   MobilizedBodyImpl(d), implementation(0) { 
        assert(userImpl);
        implementation = userImpl;
        implementation->updImpl().setReferenceToCustomImpl(this);
    }    
    
    // Copy constructor
    CustomImpl(const CustomImpl& src)
    :   MobilizedBodyImpl(src), implementation(0) {
        if (src.implementation) {
            implementation = src.implementation->clone();
            implementation->updImpl().setReferenceToCustomImpl(this);
        }
    }
    
    ~CustomImpl() {
        if (implementation)
            delete implementation;
    }
    
    CustomImpl* clone() const override { return new CustomImpl(*this); }
    
    const Custom::Implementation& getImplementation() const {
        assert(implementation);
        return *implementation;
    }
    
    Custom::Implementation& updImplementation() {
        assert(implementation);
        return *implementation;
    }

    RigidBodyNode* createRigidBodyNode(
        UIndex&        nextUSlot,
        USquaredIndex& nextUSqSlot,
        QIndex&        nextQSlot) const override;
    
    void copyOutDefaultQImpl(int nq, Real* q) const override {
        SimTK_ASSERT(nq==getImplementation().getImpl().getNQ() || nq==getImplementation().getImpl().getNQ()-1, 
            "MobilizedBody::CustomImpl::copyOutDefaultQImpl(): wrong number of q's");
        for (int i = 0; i < nq; ++i)
            q[i] = 0.0;
        if (implementation->getImpl().getNumAngles() == 4)
            q[0] = 1.0;
    }

    // Forward all the virtuals to the Custom::Implementation virtuals.
    void realizeTopologyVirtual(State& s)       const override {getImplementation().realizeTopology(s);}
    void realizeModelVirtual   (State& s)       const override {getImplementation().realizeModel(s);}
    void realizeInstanceVirtual(const State& s) const override {getImplementation().realizeInstance(s);}
    void realizeTimeVirtual    (const State& s) const override {getImplementation().realizeTime(s);}
    void realizePositionVirtual(const State& s) const override {getImplementation().realizePosition(s);}
    void realizeVelocityVirtual(const State& s) const override {getImplementation().realizeVelocity(s);}
    void realizeDynamicsVirtual(const State& s) const override {getImplementation().realizeDynamics(s);}
    void realizeAccelerationVirtual
                               (const State& s) const override {getImplementation().realizeAcceleration(s);}
    void realizeReportVirtual  (const State& s) const override {getImplementation().realizeReport(s);}
        
    void calcDecorativeGeometryAndAppend(const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
       {getImplementation().calcDecorativeGeometryAndAppend(s,stage,geom);}
    
    SimTK_DOWNCAST(CustomImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Custom;
    
    Custom::Implementation* implementation;
    
    CustomImpl& operator=(const CustomImpl&); // suppress assignment
};

// Need definition for CustomImpl here in case we have to delete it.
inline MobilizedBody::Custom::ImplementationImpl::~ImplementationImpl() {
    if (isOwner) 
        delete builtInImpl; 
    builtInImpl=0;
}

////////////////////////////////////////
// MOBILIZED BODY::FUNCTIONBASED IMPL //
////////////////////////////////////////


class MobilizedBody::FunctionBasedImpl : public MobilizedBody::Custom::Implementation {
public:
    //Constructor that uses default axes
    FunctionBasedImpl(SimbodyMatterSubsystem& matter, int nmobilities, const Array_<const Function*>& functions, const Array_<Array_<int> >& coordIndices)
            : Implementation(matter, nmobilities, nmobilities, 0), subsystem(matter.getMySubsystemIndex()), nu(nmobilities), cacheIndex(0), functions(functions), coordIndices(coordIndices), referenceCount(new int[1]) {
        assert(functions.size() == 6);
        assert(coordIndices.size() == 6);
        referenceCount[0] = 1;
        for (int i = 0; i < (int)functions.size(); ++i) {
            assert(functions[i]->getArgumentSize() == coordIndices[i].size());
            assert(functions[i]->getMaxDerivativeOrder() >= 2);
        }
        Arot = Mat33(1);
        Atrans = Mat33(1);
    }
    FunctionBasedImpl(SimbodyMatterSubsystem& matter, int nmobilities, const Array_<const Function*>& functions, const Array_<Array_<int> >& coordIndices, const Array_<Vec3>& axes)
            : Implementation(matter, nmobilities, nmobilities, 0), subsystem(matter.getMySubsystemIndex()), nu(nmobilities), cacheIndex(0), functions(functions), coordIndices(coordIndices), referenceCount(new int[1]) {
        assert(functions.size() == 6);
        assert(coordIndices.size() == 6);
        assert(axes.size() == 6);
        referenceCount[0] = 1;
        for (int i = 0; i < (int)functions.size(); ++i) {
            assert(functions[i]->getArgumentSize() == coordIndices[i].size());
            assert(functions[i]->getMaxDerivativeOrder() >= 2);
        }
        double tol = 1e-5;
        // Verify that none of the rotation axes are colinear
        assert((axes[0]%axes[1]).norm() > tol);
        assert((axes[0]%axes[2]).norm() > tol);
        assert((axes[1]%axes[2]).norm() > tol);
        // Verify that none of the translational axes are colinear
        assert((axes[3]%axes[4]).norm() > tol);
        assert((axes[3]%axes[5]).norm() > tol);
        assert((axes[4]%axes[5]).norm() > tol);

        Arot = Mat33(axes[0].normalize(), axes[1].normalize(), axes[2].normalize());
        Atrans = Mat33(axes[3].normalize(), axes[4].normalize(), axes[5].normalize());
    }
    
    ~FunctionBasedImpl() {
        if (--referenceCount[0] == 0) {
            for (int i = 0; i < (int) functions.size(); i++)
                delete functions[i];
            delete[] referenceCount;
        }
    }

    MobilizedBody::Custom::Implementation* clone() const override {
        referenceCount[0]++;
        return new FunctionBasedImpl(*this);
    }

    Transform calcMobilizerTransformFromQ(const State& s, int nq, const Real* q) const override {
        // Initialize the tranformation to be returned
        Transform X(Vec3(0));
        Vec6 spatialCoords;
        
        // Get the spatial cooridinates as a function of the q's
        for(int i=0; i < 6; i++){
            //Coordinates for this function
            int nc = coordIndices[i].size();
            Vector fcoords(nc);
    
            for(int j=0; j < nc; j++)
                fcoords(j) = q[coordIndices[i][j]];            
            
            //default behavior of constant function should take a Vector of length 0
            spatialCoords(i) = functions[i]->calcValue(fcoords);
        }

/*
        UnitVec3 axis1;
        axis1.getAs(&Arot(0,0));
        UnitVec3 axis2;
        axis2.getAs(&Arot(0,1));
        UnitVec3 axis3;
        axis3.getAs(&Arot(0,2));
*/

        //X.updR().setRotationToBodyFixedXYZ(spatialCoords.getSubVec<3>(0));
        X.updR().setRotationFromMat33TrustMe(Rotation(spatialCoords(0), UnitVec3::getAs(&Arot(0,0)))*
            Rotation(spatialCoords(1), UnitVec3::getAs(&Arot(0,1)))*Rotation(spatialCoords(2), UnitVec3::getAs(&Arot(0,2))));
        X.updP() = spatialCoords(3)*UnitVec3::getAs(&Atrans(0,0))+spatialCoords(4)*UnitVec3::getAs(&Atrans(0,1))
                    +spatialCoords(5)*UnitVec3::getAs(&Atrans(0,2));

        return X;
    }

    SpatialVec multiplyByHMatrix(const State& s, int nu, const Real* u) const override {

        switch (nu) {
            case 1: {
                // Check that the H  matrices in the cache are valid
                if (!(Value<CacheInfo<1> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH))
                    updateH(s);

                Mat<2,1,Vec3> h = Value<CacheInfo<1> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().h;
                return h*Vec1::getAs(u);
            }
            case 2: {
                // Check that the H  matrices in the cache are valid
                if (!(Value<CacheInfo<2> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH))
                    updateH(s);
                
                Mat<2,2,Vec3> h = Value<CacheInfo<2> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().h;
                return h*Vec2::getAs(u);
            }
            case 3: {
                // Check that the H  matrices in the cache are valid
                if (!(Value<CacheInfo<3> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH))
                    updateH(s);

                Mat<2,3,Vec3> h = Value<CacheInfo<3> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().h;
                return h*Vec3::getAs(u);
            }
            case 4: {
                // Check that the H  matrices in the cache are valid
                if (!(Value<CacheInfo<4> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH))
                    updateH(s);

                Mat<2,4,Vec3> h = Value<CacheInfo<4> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().h;
                return h*Vec4::getAs(u);
            }
            case 5: {
                // Check that the H matrices in the cache are valid
                if (!(Value<CacheInfo<5> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH))
                    updateH(s);

                Mat<2,5,Vec3> h = Value<CacheInfo<5> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().h;
                return h*Vec5::getAs(u);
            }
            case 6: {
                // Check that the H  matrices in the cache are valid
                if (!(Value<CacheInfo<6> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH))
                    updateH(s);

                Mat<2,6,Vec3> h = Value<CacheInfo<6> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().h;
                return h*Vec6::getAs(u);
            }
        }
        SimTK_THROW5(SimTK::Exception::ValueOutOfRange, "nu", 1, nu, 6, "MobilizedBody::FunctionBasedImpl::multiplyByHMatrix");
    }

    void multiplyByHTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const override {
        
        switch (nu) {
            case 1: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<1> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH))
                    updateH(s);

                Mat<2,1,Vec3> h = Value<CacheInfo<1> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().h;
                Vec1::updAs(f) = ~h*F;
                return;
            }
            case 2: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<2> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH))
                    updateH(s);

                Mat<2,2,Vec3> h = Value<CacheInfo<2> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().h;
                Vec2::updAs(f) = ~h*F;
                return;
            }
            case 3: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<3> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH))
                    updateH(s);

                Mat<2,3,Vec3> h = Value<CacheInfo<3> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().h;
                Vec3::updAs(f) = ~h*F;
                return;
            }
            case 4: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<4> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH))
                    updateH(s);

                Mat<2,4,Vec3> h = Value<CacheInfo<4> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().h;
                Vec4::updAs(f) = ~h*F;
                return;
            }
            case 5: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<5> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH))
                    updateH(s);

                Mat<2,5,Vec3> h = Value<CacheInfo<5> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().h;
                Vec5::updAs(f) = ~h*F;
                return;
            }
            case 6: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<6> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH))
                    updateH(s);

                Mat<2,6,Vec3> h = Value<CacheInfo<6> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().h;
                Vec6::updAs(f) = ~h*F;
                return;
            }
        }
        SimTK_THROW5(SimTK::Exception::ValueOutOfRange, "nu", 1, nu, 6, "MobilizedBody::FunctionBasedImpl::multiplyByHTranspose");
    }

    SpatialVec multiplyByHDotMatrix(const State& s, int nu, const Real* u) const override {

        switch (nu) {
            case 1: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<1> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot))
                    updateHdot(s);

                Mat<2,1,Vec3> hdot = Value<CacheInfo<1> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().hdot;
                return hdot*Vec1::getAs(u);
            }
            case 2: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<2> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot))
                    updateHdot(s);

                Mat<2,2,Vec3> hdot = Value<CacheInfo<2> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().hdot;
                return hdot*Vec2::getAs(u);
            }
            case 3: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<3> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot))
                    updateHdot(s);

                Mat<2,3,Vec3> hdot = Value<CacheInfo<3> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().hdot;
                return hdot*Vec3::getAs(u);
            }
            case 4: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<4> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot))
                    updateHdot(s);

                Mat<2,4,Vec3> hdot = Value<CacheInfo<4> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().hdot;
                return hdot*Vec4::getAs(u);
            }
            case 5: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<5> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot))
                    updateHdot(s);

                Mat<2,5,Vec3> hdot = Value<CacheInfo<5> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().hdot;
                return hdot*Vec5::getAs(u);
            }
            case 6: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<6> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot))
                    updateHdot(s);

                Mat<2,6,Vec3> hdot = Value<CacheInfo<6> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().hdot;
                return hdot*Vec6::getAs(u);
            }
        }
        SimTK_THROW5(SimTK::Exception::ValueOutOfRange, "nu", 1, nu, 6, "MobilizedBody::FunctionBasedImpl::multiplyByHDotMatrix");
    }

    void multiplyByHDotTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const override {

        switch (nu) {
            case 1: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<1> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot))
                    updateHdot(s);

                Mat<2,1,Vec3> hdot = Value<CacheInfo<1> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().hdot;
                Vec1::updAs(f) = ~hdot*F;
                return;
            }
            case 2: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<2> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot))
                    updateHdot(s);

                Mat<2,2,Vec3> hdot = Value<CacheInfo<2> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().hdot;
                Vec2::updAs(f) = ~hdot*F;
                return;
            }
            case 3: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<3> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot))
                    updateHdot(s);

                Mat<2,3,Vec3> hdot = Value<CacheInfo<3> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().hdot;
                Vec3::updAs(f) = ~hdot*F;
                return;
            }
            case 4: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<4> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot))
                    updateHdot(s);

                Mat<2,4,Vec3> hdot = Value<CacheInfo<4> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().hdot;
                Vec4::updAs(f) = ~hdot*F;
                return;
            }
            case 5: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<5> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot))
                    updateHdot(s);

                Mat<2,5,Vec3> hdot = Value<CacheInfo<5> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().hdot;
                Vec5::updAs(f) = ~hdot*F;
                return;
            }
            case 6: {
                // Check that the H and Hdot matrices in the cache are valid
                if (!(Value<CacheInfo<6> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot))
                    updateHdot(s);

                Mat<2,6,Vec3> hdot = Value<CacheInfo<6> >::downcast(s.getCacheEntry(subsystem, cacheIndex)).get().hdot;
                Vec6::updAs(f) = ~hdot*F;
                return;
            }
        }
        SimTK_THROW5(SimTK::Exception::ValueOutOfRange, "nu", 1, nu, 6, "MobilizedBody::FunctionBasedImpl::multiplyByHDotTranspose");
    }

    void realizeTopology(State& s) const override {
        switch (nu) {
        case 1:
            cacheIndex = s.allocateCacheEntry(subsystem, Stage::Topology, new Value<CacheInfo<1> >());
            break;
        case 2:
            cacheIndex = s.allocateCacheEntry(subsystem, Stage::Topology, new Value<CacheInfo<2> >());
            break;
        case 3:
            cacheIndex = s.allocateCacheEntry(subsystem, Stage::Topology, new Value<CacheInfo<3> >());
            break;
        case 4:
            cacheIndex = s.allocateCacheEntry(subsystem, Stage::Topology, new Value<CacheInfo<4> >());
            break;
        case 5:
            cacheIndex = s.allocateCacheEntry(subsystem, Stage::Topology, new Value<CacheInfo<5> >());
            break;
        case 6:
            cacheIndex = s.allocateCacheEntry(subsystem, Stage::Topology, new Value<CacheInfo<6> >());
            break;
        }
    }

    void realizePosition(const State& s) const override {
        switch (nu) {
            case 1: {
                // invalidate H and Hdot matrices
                Value<CacheInfo<1> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH = false;
                break;
            }
            case 2: {
                // invalidate H and Hdot matrices
                Value<CacheInfo<2> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH = false;
                break;
            }
            case 3: {
                // invalidate H and Hdot matrices
                Value<CacheInfo<3> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH = false;
                break;
            }
            case 4: {
                // invalidate H and Hdot matrices
                Value<CacheInfo<4> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH = false;
                break;
            }
            case 5: {
                // invalidate H and Hdot matrices
                Value<CacheInfo<5> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH = false;
                break;
            }
            case 6: {
                // invalidate H and Hdot matrices
                Value<CacheInfo<6> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidH = false;
                break;
            }
        }
    }

    void realizeVelocity(const State& s) const override {
        switch (nu) {
            case 1: {
                // invalidate H and Hdot matrices
                Value<CacheInfo<1> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot = false;
                break;
            }
            case 2: {
                // invalidate H and Hdot matrices
                Value<CacheInfo<2> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot = false;
                break;
            }
            case 3: {
                // invalidate H and Hdot matrices
                Value<CacheInfo<3> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot = false;
                break;
            }
            case 4: {
                // invalidate H and Hdot matrices
                Value<CacheInfo<4> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot = false;
                break;
            }
            case 5: {
                // invalidate H and Hdot matrices
                Value<CacheInfo<5> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot = false;
                break;
            }
            case 6: {
                // invalidate H and Hdot matrices
                Value<CacheInfo<6> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd().isValidHdot = false;
                break;
            }
        }
    }

    void updateH(const State& s) const {
        // Get mobilizer kinematics 
        Vector q = getQ(s);
        Vector u = getU(s);
        switch (nu) {
            case 1: {
                CacheInfo<1>& cache = Value<CacheInfo<1> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd();
                cache.buildH(q, u, getMobilizerTransform(s), functions, coordIndices, Arot, Atrans);
                // H matrix is now valid 
                cache.isValidH = true;
                break;
            }
            case 2: {
                CacheInfo<2>& cache = Value<CacheInfo<2> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd();
                cache.buildH(q, u, getMobilizerTransform(s), functions, coordIndices, Arot, Atrans);
                // H matrix is now valid 
                cache.isValidH = true;
                break;
            }
            case 3: {
                CacheInfo<3>& cache = Value<CacheInfo<3> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd();
                cache.buildH(q, u, getMobilizerTransform(s), functions, coordIndices, Arot, Atrans);
                // H matrix is now valid  
                cache.isValidH = true;
                break;
            }
            case 4: {
                CacheInfo<4>& cache = Value<CacheInfo<4> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd();
                cache.buildH(q, u, getMobilizerTransform(s), functions, coordIndices, Arot, Atrans);
                // H matrix is now valid 
                cache.isValidH = true;
                break;
            }
            case 5: {
                CacheInfo<5>& cache = Value<CacheInfo<5> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd();
                cache.buildH(q, u, getMobilizerTransform(s), functions, coordIndices, Arot, Atrans);
                // H matrix is now valid  
                cache.isValidH = true;
                break;
            }
            case 6: {
                CacheInfo<6>& cache = Value<CacheInfo<6> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd();
                cache.buildH(q, u, getMobilizerTransform(s), functions, coordIndices, Arot, Atrans);
                // H matrix is now valid 
                cache.isValidH = true;
                break;
            }
        }
    }
 
    void updateHdot(const State& s) const {
        // Get mobilizer kinematics 
        Vector q = getQ(s);
        Vector u = getU(s);
        switch (nu) {
            case 1: {
                CacheInfo<1>& cache = Value<CacheInfo<1> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd();
                cache.buildHdot(q, u, getMobilizerTransform(s), functions, coordIndices, Arot, Atrans);
                // Hdot matrix is now valid 
                cache.isValidHdot = true;
                break;
            }
            case 2: {
                CacheInfo<2>& cache = Value<CacheInfo<2> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd();
                cache.buildHdot(q, u, getMobilizerTransform(s), functions, coordIndices, Arot, Atrans);
                // Hdot matrix is now valid 
                cache.isValidHdot = true;
                break;
            }
            case 3: {
                CacheInfo<3>& cache = Value<CacheInfo<3> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd();
                cache.buildHdot(q, u, getMobilizerTransform(s), functions, coordIndices, Arot, Atrans);
                // Hdot matrix is now valid  
                cache.isValidHdot = true;
                break;
            }
            case 4: {
                CacheInfo<4>& cache = Value<CacheInfo<4> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd();
                cache.buildHdot(q, u, getMobilizerTransform(s), functions, coordIndices, Arot, Atrans);
                // Hdot matrix is now valid 
                cache.isValidHdot = true;
                break;
            }
            case 5: {
                CacheInfo<5>& cache = Value<CacheInfo<5> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd();
                cache.buildHdot(q, u, getMobilizerTransform(s), functions, coordIndices, Arot, Atrans);
                // Hdot matrix is now valid 
                cache.isValidHdot = true;
                break;
            }
            case 6: {
                CacheInfo<6>& cache = Value<CacheInfo<6> >::downcast(s.updCacheEntry(subsystem, cacheIndex)).upd();
                cache.buildHdot(q, u, getMobilizerTransform(s), functions, coordIndices, Arot, Atrans);
                // Hdot matrix is now valid 
                cache.isValidHdot = true;
                break;
            }
        }
    }

private:
    const SubsystemIndex subsystem;
    const int nu;
    mutable CacheEntryIndex cacheIndex;
    const Array_<const Function*> functions;
    const Array_<Array_<int> > coordIndices;
    int* referenceCount;
    //const Array_<Vec3> axes;
    Mat33 Arot, Atrans;
    template <int N> class CacheInfo {
    public:
        CacheInfo() : isValidH(false), isValidHdot(false) { }

        void buildH(Vector& q, Vector& u, const Transform& X_FM, const Array_<const Function*>& functions, const Array_<Array_<int> >& coordIndices, const Mat33 Arot, const Mat33 Atrans)
        {
            // Build the Fq and Fqq matrices of partials of the spatial functions with respect to the gen coordinates, q    
            // Cycle through each row (function describing spatial coordinate)
            Fq = Mat<6,N>(0);
            Vec6 spatialCoords(0);
            Array_<int> deriv(1);
            Vector fcoords(coordIndices[0].size()); 

            for(int i=0; i < 6; i++){
                // Determine the number of coordinates for this function
                int nc = coordIndices[i].size();
                if (fcoords.size() != nc)
                    fcoords.resize(nc);

                if (nc > 0) {
                    // Get coordinate values to evaluate the function
                    for(int k = 0; k < nc; k++)
                        fcoords(k) = q(coordIndices[i][k]);

                    for (int j = 0; j < nc; j++) {
                        deriv[0] = j;
                        Fq(i, coordIndices[i][j]) = functions[i]->calcDerivative(deriv, fcoords);
                    }
                    //default behavior of constant function should take a Vector of length 0
                    spatialCoords(i) = functions[i]->calcValue(fcoords);
                }

            }
            // Note rotations are body fixed sequences, and so taking the derivative yields a body-fixed angular velocity
            // and therefore, must be transformed to the parent frame.

            // omega = [a1, R_F1*a2, R_F2*a3]*{Theta_dot(q)}, where Theta_dot(q) is described by the first three functions
            //       = [W*[Fq_theta]*qdot
            // vel = [A]*{X_dot(q)}, where X_dot(q) is described by the last three functions
            //     = [A]*[Fq_x]*qdot
            
            Rotation R_F1 = Rotation(spatialCoords(0), UnitVec3::getAs(&Arot(0,0)));
            Rotation R_F2 = R_F1*Rotation(spatialCoords(1), UnitVec3::getAs(&Arot(0,1)));

            Mat31 temp;
            Mat33 W(UnitVec3::getAs(&Arot(0,0)), R_F1*(UnitVec3::getAs(&Arot(0,1))), R_F2*(UnitVec3::getAs(&Arot(0,2))));

            for(int i=0; i < N; i++){
                temp = W*(Fq.template getSubMat<3,1>(0,i));
                h(0,i) = Vec3::getAs(&temp(0,0));
                temp = Atrans*(Fq.template getSubMat<3,1>(3,i));
                h(1,i) = Vec3::getAs(&temp(0,0));
            }
        }

        void buildHdot(Vector& q, Vector& u, const Transform& X_FM, const Array_<const Function*>& functions, const Array_<Array_<int> >& coordIndices, const Mat33 Arot, const Mat33 Atrans)
        {
            Mat<6,N> Fqdot(0);
            Vec6 spatialCoords;
            Array_<int> derivs(2);
            Vector fcoords(coordIndices[0].size()); 

            for(int i=0; i < 6; i++){
                // Determine the number of coordinates for this function
                int nc = coordIndices[i].size();
                if (fcoords.size() != nc)
                    fcoords.resize(nc);
                
                if (nc > 0) {
                    // Get coordinate values to evaluate the function
                    for(int k = 0; k < nc; k++)
                        fcoords(k) = q(coordIndices[i][k]);

                    // function is dependent on a mobility if its index is in the list of function coordIndices
                    // cycle through the mobilities
                    for (int j = 0; j < nc; j++) {
                        derivs[0] = j;
                        for (int k = 0; k < nc; k++) {
                            derivs[1] = k;
                            Fqdot(i, coordIndices[i][j]) += functions[i]->calcDerivative(derivs, fcoords)*u[coordIndices[i][k]];
                        }
                    }
                }
                //default behavior of constant function should take a Vector of length 0
                spatialCoords(i) = functions[i]->calcValue(fcoords);
            }

            Rotation R_F1 = Rotation(spatialCoords(0), UnitVec3::getAs(&Arot(0,0)));
            Rotation R_F2 = R_F1*Rotation(spatialCoords(1), UnitVec3::getAs(&Arot(0,1)));

            Mat33 W(UnitVec3::getAs(&Arot(0,0)), R_F1*(UnitVec3::getAs(&Arot(0,1))), R_F2*(UnitVec3::getAs(&Arot(0,2))));

            // Hdot_theta = [Wdot]*[Fq] + [W][Fqdot]
            // Hdot_x = [A]*[Fqdot]
            Real *up = &u[0];
            Vec<N> uv = Vec<N>::getAs(up);
            Vec<N> uv1 = uv;
            Vec<N> uv2 = uv;
            
            for(int i=1; i < N; i++){
                uv1(i) = 0;
            }
            if(N > 2)
                uv2(2) = 0;

            // Spatial velocity
            SpatialVec V = h*uv;
            SpatialVec V1 = h*uv1;
            SpatialVec V2 = h*uv2;
            // First V[0] is angular velocity, omega
            //Mat33 Wdot(Vec3(0), Vec3(V[0](0),0,0)%(Vec3::getAs(&W(0,1))), Vec3(V[0](0),V[0](1),0)%(Vec3::getAs(&W(0,2))));
            Mat33 Wdot(Vec3(0), V1[0]%(Vec3::getAs(&W(0,1))), V2[0]%(Vec3::getAs(&W(0,2))));
            
            //Sanity check Omega == V[0]
            Mat31 Omega = W*(Fq.template getSubMat<3,N>(0,0))*(Mat<N,1>::getAs(up));

            Mat31 temp;
            for(int i=0; i < N; i++){
                temp = Wdot*(Fq.template getSubMat<3,1>(0,i))+W*(Fqdot.template getSubMat<3,1>(0,i));
                hdot(0,i) = Vec3::getAs(&temp(0,0));
                temp = Atrans*(Fqdot.template getSubMat<3,1>(3,i));
                hdot(1,i) = Vec3::getAs(&temp(0,0));
            }
        }

        Mat<2,N,Vec3> h, hdot;
        Mat<6,N> Fq;
        bool isValidH;
        bool isValidHdot;
    };

};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_IMPL_H_
