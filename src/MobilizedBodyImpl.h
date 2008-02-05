#ifndef SimTK_SIMBODY_MOBILIZED_BODY_REP_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_REP_H_

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

/**@file
 *
 * Private implementation of Body and MobilizedBody, and their included subclasses which
 * represent the built-in body and mobilizer types.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Body.h"
#include "simbody/internal/MobilizedBody.h"

#include "SimbodyMatterSubsystemRep.h"
#include "BodyRep.h"

class RigidBodyNode;

namespace SimTK {

    /////////////////////////
    // MOBILIZED BODY REPS //
    /////////////////////////

class MobilizedBodyImpl : public PIMPLImplementation<MobilizedBody,MobilizedBodyImpl> {
public:
    MobilizedBodyImpl() : myHandle(0), myMatterSubsystemRep(0), myLevel(-1), myRBnode(0) {
    }
    MobilizedBodyImpl(const MobilizedBodyImpl& clone) {
        *this = clone;
        myRBnode = 0;
    }
    virtual ~MobilizedBodyImpl() {
        delete myRBnode;
    }
    virtual MobilizedBodyImpl* clone() const = 0;
    
    // This creates a rigid body node using the appropriate mobilizer type.
    // Called during realizeTopology().
    virtual RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const = 0; 

    virtual void realizeTopologyImpl(State&) const { }

    void realizeTopology(State& s) const {
        realizeTopologyImpl(s);
        //TODO: more? Move RigidBodyNode creation here?
    }

    // Copy out nq default values for q, beginning at the indicated address.
    // The concrete class should assert if nq is not a reasonable
    // number for the kind of mobilizer; there is a bug somewhere in that case. 
    // This routine shouldn't be called directly -- call copyOutDefaultQ() below 
    // instead which has a nicer interface and does some error checking.
    virtual void copyOutDefaultQImpl(int nq, Real* q) const = 0;

    virtual void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
    {
        if (stage != Stage::Instance || !getMyMatterSubsystemRep().getShowDefaultGeometry())
            return;
        const Real scale = 1.0;
        DecorativeFrame axes(scale*0.5);
        axes.setLineThickness(2);
        axes.setBodyId(myMobilizedBodyIndex);
        geom.push_back(axes); // the body frame

        // Display the inboard joint frame (at half size), unless it is the
        // same as the body frame. Then find the corresponding frame on the
        // parent and display that in this body's color.
        if (myMobilizedBodyIndex != 0) {
            const Real pscale = 1.0;
            const Transform& M = getDefaultOutboardFrame(); // TODO: get from state
            if (M.T() != Vec3(0) || M.R() != Mat33(1)) {
                geom.push_back(DecorativeFrame(scale*0.25).setBodyId(myMobilizedBodyIndex));
                if (M.T() != Vec3(0))
                    geom.push_back(DecorativeLine(Vec3(0), M.T()).setBodyId(myMobilizedBodyIndex));
            }
            const Transform& Mb = getDefaultInboardFrame(); // TODO: from state
            DecorativeFrame frameOnParent(pscale*0.25);
            frameOnParent.setBodyId(myParentIndex);
            frameOnParent.setColor(Black);
            frameOnParent.setTransform(Mb);
            geom.push_back(frameOnParent);
            if (Mb.T() != Vec3(0))
                geom.push_back(DecorativeLine(Vec3(0),Mb.T()).setBodyId(myParentIndex));
        }

        // Put a little purple wireframe sphere at the COM, and add a line from 
        // body origin to the com.

        DecorativeSphere com(scale*.05);
        com.setBodyId(myMobilizedBodyIndex);
        com.setColor(Purple).setRepresentation(DecorativeGeometry::DrawPoints);
        const Vec3& comPos_B = theBody.getDefaultRigidBodyMassProperties().getMassCenter(); // TODO: from state
        com.setTransform(comPos_B);
        geom.push_back(com);
        if (comPos_B != Vec3(0))
            geom.push_back(DecorativeLine(Vec3(0), comPos_B).setBodyId(myMobilizedBodyIndex));
    }

    void calcDecorativeGeometryAndAppend
       (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
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

    void addOutboardDecoration(const Transform& X_MD, const DecorativeGeometry& g) {
        outboardGeometry.push_back(g); // make a new copy
        // Combine the placement frame and the transform already in the geometry
        // so we end up with geometry expressed directly in the M frame.
        outboardGeometry.back().setTransform(X_MD*g.getTransform());
    }
    void addInboardDecoration(const Transform& X_FD, const DecorativeGeometry& g) {
        inboardGeometry.push_back(g); // make a new copy
        // Combine the placement frame and the transform already in the geometry
        // so we end up with geometry expressed directly in the F frame.
        inboardGeometry.back().setTransform(X_FD*g.getTransform());
    }


    void findMobilizerQs(const State&, int& qStart, int& nq) const;
    void findMobilizerUs(const State&, int& uStart, int& nu) const;

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
    void setQToFitTranslation(State& s, const Vec3& T_FM, 
                                 bool dontChangeOrientation) const;

    // Invalidate Stage::Velocity.
    void setUToFitVelocity(State& s, const SpatialVec& V_FM) const;
    void setUToFitAngularVelocity(State& s, const Vec3& w_FM) const;
    void setUToFitLinearVelocity(State& s, const Vec3& v_FM,
                                 bool dontChangeAngularVelocity)  const;

    const MassProperties& getBodyMassProperties(const State& s) const {
        // TODO: these should come from the state if the body has variable mass props
        const SBInstanceVars& iv = getMyMatterSubsystemRep().getInstanceVars(s);
        return getMyRigidBodyNode().getMassProperties_OB_B();
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
        const SBPositionCache& pc = getMyMatterSubsystemRep().getPositionCache(s);
        return getMyRigidBodyNode().getX_GB(pc);
    }
    const SpatialVec& getBodyVelocity(const State& s) const {
        const SBVelocityCache& vc = getMyMatterSubsystemRep().getVelocityCache(s);
        return getMyRigidBodyNode().getV_GB(vc);
    }
    const SpatialVec& getBodyAcceleration(const State& s) const {
        const SBAccelerationCache& ac = getMyMatterSubsystemRep().getAccelerationCache(s);
        return getMyRigidBodyNode().getA_GB(ac);
    }

    const Transform& getMobilizerTransform(const State& s) const {
        const SBPositionCache& pc = getMyMatterSubsystemRep().getPositionCache(s);
        return getMyRigidBodyNode().getX_FM(pc);
    }
    const SpatialVec& getMobilizerVelocity(const State& s) const {
        const SBVelocityCache& vc = getMyMatterSubsystemRep().getVelocityCache(s);
        return getMyRigidBodyNode().getV_FM(vc);
    }

    const RigidBodyNode& realizeTopology(int& nxtU, int& nxtUSq, int& nxtQ) const;

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
        assert(!myMatterSubsystemRep);
        myMatterSubsystemRep = &matter.updRep();

        assert(index.isValid());
        assert(parentIndex.isValid() || index==GroundIndex);

        myParentIndex           = parentIndex; // invalid if this is Ground
        myMobilizedBodyIndex    = index;

        if (index != GroundIndex) {
            const MobilizedBody& parent = matter.getMobilizedBody(parentIndex);
            myLevel = parent.getLevelInMultibodyTree() + 1;
            myBaseBodyIndex = (myLevel == 1 ? myMobilizedBodyIndex 
                                         : parent.getBaseMobilizedBody().getMobilizedBodyIndex());
        } else {
            myLevel = 0;
            myBaseBodyIndex = GroundIndex;
        }
    }

    void setMyHandle(MobilizedBody& h) {myHandle = &h;}
    const MobilizedBody& getMyHandle() const {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}

private:
    // Body topological geometry is defined with respect to the body frame so we
    // can draw it right away.
    void appendTopologicalBodyGeometry(std::vector<DecorativeGeometry>& geom) const {
        getBody().getRep().appendDecorativeGeometry(getMyMobilizedBodyIndex(), geom);
    }

    // Mobilizer topological geometry is defined with respect to the M (outboard, child)
    // frame and the F (inboard, parent) frame. The placement of those frames with
    // respect to the body frame can be Instance variables, so we can't draw this geometry
    // until Instance stage. At that point we can find M and F, so they are passed in
    // here.
    void appendTopologicalMobilizerGeometry(const Transform& X_BM, const Transform& X_PF,
                                            std::vector<DecorativeGeometry>& geom) const
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
    MobilizedBody* myHandle;	// the owner handle of this rep

        // TOPOLOGY "STATE"

    // Base class topological properties. Derived MobilizedBodies may
    // have further topological properties. Whenever these change, be
    // sure to set topologyRealized=false!
    Body theBody;
    Transform defaultInboardFrame;  // default for F (in Parent frame)
    Transform defaultOutboardFrame; // default for M (in Body frame)

    std::vector<DecorativeGeometry> outboardGeometry;
    std::vector<DecorativeGeometry> inboardGeometry;

    // These data members are filled in once the MobilizedBody is added to
    // a MatterSubsystem. Note that this pointer is just a reference to
    // the Subsystem which owns this MobilizedBody -- we don't own the
    // heap space here so don't delete it!
    SimbodyMatterSubsystemRep* myMatterSubsystemRep;
    MobilizedBodyIndex  myMobilizedBodyIndex; // index within the subsystem
    MobilizedBodyIndex  myParentIndex;   // Invalid if this body is Ground, otherwise the parent's index
    MobilizedBodyIndex  myBaseBodyIndex; // GroundIndex if this is ground, otherwise a level-1 BodyIndex
    int              myLevel;      // Distance from ground in multibody graph (0 if this is ground,
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
};

class MobilizedBody::PinImpl : public MobilizedBodyImpl {
public:
    PinImpl() : defaultQ(0) { }
    PinImpl* clone() const { return new PinImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==1, 
            "MobilizedBody::Pin::PinImpl::copyOutDefaultQImpl(): wrong number of q's expected");
        *q = defaultQ;
    }

    SimTK_DOWNCAST(PinImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Pin;
    Real defaultQ;  // the angle in radians
};


class MobilizedBody::SliderImpl : public MobilizedBodyImpl {
public:
    SliderImpl() : defaultQ(0) { }
    SliderImpl* clone() const { return new SliderImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==1, 
            "MobilizedBody::Slider::SliderImpl::copyOutDefaultQImpl(): wrong number of q's expected");
        *q = defaultQ;
    }

    SimTK_DOWNCAST(SliderImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Slider;
    Real defaultQ;  // the displacement
};

class MobilizedBody::UniversalImpl : public MobilizedBodyImpl {
public:
    UniversalImpl() : defaultQ(0) { }
    UniversalImpl* clone() const { return new UniversalImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==2, 
            "MobilizedBody::Universal::UniversalImpl::copyOutDefaultQImpl(): wrong number of q's expected");
        Vec2::updAs(q) = defaultQ;
    }

    SimTK_DOWNCAST(UniversalImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Universal;
    Vec2 defaultQ;  // the two angles
};

class MobilizedBody::CylinderImpl : public MobilizedBodyImpl {
public:
    CylinderImpl() : defaultQ(0) { }
    CylinderImpl* clone() const { return new CylinderImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==2, 
            "MobilizedBody::Cylinder::CylinderImpl::copyOutDefaultQImpl(): wrong number of q's expected");
        Vec2::updAs(q) = defaultQ;
    }

    SimTK_DOWNCAST(CylinderImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Cylinder;
    Vec2 defaultQ;  // the angle in radians followed by the displacement
};

class MobilizedBody::BendStretchImpl : public MobilizedBodyImpl {
public:
    BendStretchImpl() : defaultQ(0) { }
    BendStretchImpl* clone() const { return new BendStretchImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==2, 
            "MobilizedBody::BendStretch::BendStretchImpl::copyOutDefaultQImpl(): wrong number of q's expected");
        Vec2::updAs(q) = defaultQ;
    }

    SimTK_DOWNCAST(BendStretchImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::BendStretch;
    Vec2 defaultQ;  // the angle in radians followed by the displacement
};

class MobilizedBody::PlanarImpl : public MobilizedBodyImpl {
public:
    PlanarImpl() : defaultQ(0) { }
    PlanarImpl* clone() const { return new PlanarImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==3, 
            "MobilizedBody::Planar::PlanarImpl::copyOutDefaultQImpl(): wrong number of q's expected");
        Vec3::updAs(q) = defaultQ;
    }

    SimTK_DOWNCAST(PlanarImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Planar;
    Vec3 defaultQ;  // the angle in radians followed by the two displacements
};

class MobilizedBody::GimbalImpl : public MobilizedBodyImpl {
public:
    GimbalImpl() : defaultQ(0) { }
    GimbalImpl* clone() const { return new GimbalImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==3, 
            "MobilizedBody::Gimbal::GimbalImpl::copyOutDefaultQImpl(): wrong number of q's expected");
        Vec3::updAs(q) = defaultQ;
    }

    void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const;

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

class MobilizedBody::BallImpl : public MobilizedBodyImpl {
public:
    BallImpl() : defaultRadius(0.1), defaultQ() { } // default is (1,0,0,0), the identity rotation
    BallImpl* clone() const { return new BallImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==4||nq==3, 
            "MobilizedBody::Ball::BallImpl::copyOutDefaultQImpl(): wrong number of q's expected");
        if (nq==4)
            Vec4::updAs(q) = defaultQ.asVec4();
        else
            Vec3::updAs(q) = Rotation(defaultQ).convertRotationToBodyFixedXYZ();
    }

    void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const;

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
    EllipsoidImpl() : defaultRadii(0.5,1/3.,0.25), defaultQ() { } // default is (1,0,0,0), the identity rotation
    EllipsoidImpl* clone() const { return new EllipsoidImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==4||nq==3, 
            "MobilizedBody::Ellipsoid::EllipsoidImpl::copyOutDefaultQImpl(): wrong number of q's expected");
        if (nq==4)
            Vec4::updAs(q) = defaultQ.asVec4();
        else
            Vec3::updAs(q) = Rotation(defaultQ).convertRotationToBodyFixedXYZ();
    }

    void calcDecorativeGeometryAndAppendImpl
       (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const;

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
    TranslationImpl() : defaultQ(0) { }
    TranslationImpl* clone() const { return new TranslationImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==3, 
            "MobilizedBody::Translation::TranslationImpl::copyOutDefaultQImpl(): wrong number of q's expected");
        Vec3::updAs(q) = defaultQ;
    }

    SimTK_DOWNCAST(TranslationImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Translation;
    Vec3 defaultQ;  // the three displacements
};

class MobilizedBody::FreeImpl : public MobilizedBodyImpl {
public:
    FreeImpl() : defaultQOrientation(), defaultQTranslation(0) { }
    FreeImpl* clone() const { return new FreeImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==7||nq==6, 
            "MobilizedBody::Free::FreeImpl::copyOutDefaultQImpl(): wrong number of q's expected");
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
    LineOrientationImpl() : defaultQ() { } // 1,0,0,0
    LineOrientationImpl* clone() const { return new LineOrientationImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==4||nq==3, 
            "MobilizedBody::LineOrientation::LineOrientationImpl::copyOutDefaultQImpl(): wrong number of q's expected");
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
    FreeLineImpl() : defaultQOrientation(), defaultQTranslation(0) { }
    FreeLineImpl* clone() const { return new FreeLineImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==7||nq==6, 
            "MobilizedBody::FreeLine::FreeLineImpl::copyOutDefaultQImpl(): wrong number of q's expected");
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
    WeldImpl() { }
    WeldImpl* clone() const { return new WeldImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==0, 
            "MobilizedBody::Weld::WeldImpl::copyOutDefaultQImpl(): wrong number of q's expected");
    }

    SimTK_DOWNCAST(WeldImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Weld;
    // no q's
};

class MobilizedBody::GroundImpl : public MobilizedBodyImpl {
public:
    GroundImpl() { }
    GroundImpl* clone() const { return new GroundImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==0, 
            "MobilizedBody::Ground::GroundImpl::copyOutDefaultQImpl(): wrong number of q's expected");
    }

    SimTK_DOWNCAST(GroundImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Ground;
    // no q's
};

class MobilizedBody::ScrewImpl : public MobilizedBodyImpl {
public:
    explicit ScrewImpl(Real p) : defaultPitch(p), defaultQ(0) { }

    Real getDefaultPitch() const {return defaultPitch;}
    void setDefaultPitch(Real p) {
        invalidateTopologyCache();
        defaultPitch=p;
    }

    ScrewImpl* clone() const { return new ScrewImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int nq, Real* q) const {
        SimTK_ASSERT(nq==1, 
            "MobilizedBody::Screw::ScrewImpl::copyOutDefaultQImpl(): wrong number of q's expected");
        *q = defaultQ;
    }

    SimTK_DOWNCAST(ScrewImpl, MobilizedBodyImpl);
private:
    friend class MobilizedBody::Screw;
    Real defaultPitch;
    Real defaultQ; // the angle in radians
};

class MobilizedBody::CustomImpl : public MobilizedBodyImpl {
public:
    CustomImpl(int nMobilities, int nCoordinates) 
      : nu(nMobilities), nq(nCoordinates), defaultQ(nq) {
        defaultQ.setToZero();
    }
    CustomImpl* clone() const { return new CustomImpl(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    void copyOutDefaultQImpl(int expectedNq, Real* q) const {
        SimTK_ASSERT(expectedNq==nq, 
            "MobilizedBody::Custom::CustomImpl::copyOutDefaultQImpl(): wrong number of q's expected");
        for (int i=0; i<nq; ++i)
            q[i] = defaultQ[i];
    }

    SimTK_DOWNCAST(CustomImpl, MobilizedBodyImpl);
private:
    int nu;
    int nq;
    Vector defaultQ;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_REP_H_




