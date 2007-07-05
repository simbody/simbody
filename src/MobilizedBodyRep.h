#ifndef SimTK_SIMBODY_MOBILIZED_BODY_REP_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_REP_H_

/* Portions copyright (c) 2007 Stanford University and Michael Sherman.
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
 * Private implementation of Body and MobilizedBody, and their included subclasses which
 * represent the built-in body and mobilizer types.
 */

#include "SimTKsimbody.h"
#include "simbody/internal/MobilizedBody.h"

#include "SimbodyMatterSubsystemRep.h"
#include "ConstraintNode.h"

class RigidBodyNode;

namespace SimTK {

    ///////////////
    // BODY REPS //
    ///////////////

class Body::BodyRep {
public:
    BodyRep() : myHandle(0) {
    }
    virtual ~BodyRep() { }
    virtual BodyRep* clone() const = 0;

    virtual const MassProperties& getDefaultRigidBodyMassProperties() const = 0;

    void setMyHandle(Body& h) {myHandle = &h;}
    const Body& getMyBodyHandle() const {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}
private:
    friend class Body;
    Body* myHandle;
};

class Body::Ground::GroundRep : public Body::BodyRep {
public:
    GroundRep() : BodyRep() {
    }
    GroundRep* clone() const {
        return new GroundRep(*this);
    }
    const MassProperties& getDefaultRigidBodyMassProperties() const {
        static const MassProperties groundMass(NTraits<Real>::getInfinity(), 
                                               Vec3(0),
                                               NTraits<Real>::getInfinity()*Inertia(1));
        return groundMass;
    }

    const Body::Ground& getMyGroundBodyHandle() const {
        return reinterpret_cast<const Body::Ground&>(getMyBodyHandle());
    }

    SimTK_DOWNCAST(GroundRep, BodyRep);
private:
};

class Body::Rigid::RigidRep : public Body::BodyRep {
public:
    RigidRep() : BodyRep() {
    }
    RigidRep(const MassProperties& m) : BodyRep(), defaultMassProperties(m) {
    }
    const MassProperties& getDefaultRigidBodyMassProperties() const {
        return defaultMassProperties;
    }
    void setDefaultMassProperties(const MassProperties& m) {
        defaultMassProperties = m;
    }
    RigidRep* clone() const {
        return new RigidRep(*this);
    }
    const Body::Rigid& getMyRigidBodyHandle() const {
        return reinterpret_cast<const Body::Rigid&>(getMyBodyHandle());
    }

    SimTK_DOWNCAST(RigidRep, BodyRep);
private:
    MassProperties defaultMassProperties;
};

    /////////////////////////
    // MOBILIZED BODY REPS //
    /////////////////////////

class MobilizedBody::MobilizedBodyRep {
public:
    MobilizedBodyRep() : myHandle(0), myMatterSubsystemRep(0), myRBnode(0) {
    }
    virtual ~MobilizedBodyRep() {
        delete myRBnode;
    }
    virtual MobilizedBodyRep* clone() const = 0;
    
    // This creates a rigid body node using the appropriate mobilizer type.
    // Called during realizeTopology().
    virtual RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const = 0; 

    int getQIndex(const State&) const;
    int getUIndex(const State&) const;

    // Given the Model stage state variables (if any are relevant), each
    // concrete mobilizer should return the number of q's and u's it expects
    // to be allotted.
    //virtual void getStateNeeds(const State& s, int& nq, int& nu) const = 0;

    // Given the Model stage variable values in the passed-in State (if
    // any of them are relevant) copy out the appropriate default values
    // to the indicated memory location. The expected number of q's is passed
    // in for good hygiene -- it must match what the mobilizer thinks it should
    // have for these modeling options or something is terribly wrong and an 
    // exception will be thrown.
    //virtual void copyOutDefaultQs(const State&, int expectedNq, Real* qDefault) const = 0;

    // Generic State-access routines (that is, those which can be handled in
    // the MobilizedBody base class).
    const Transform& getBodyTransform(const State& s) const {
        const SBPositionCache& pc = getMyMatterSubsystemRep().getPositionCache(s);
        return getMyRigidBodyNode().getX_GB(pc);
    }
    const SpatialVec& getBodyVelocity(const State& s) const {
        const SBVelocityCache& vc = getMyMatterSubsystemRep().getVelocityCache(s);
        return getMyRigidBodyNode().getV_GB(vc);
    }
    const SpatialVec& getBodyAppliedForces(const State& s) const {
        const SBAccelerationVars& av = getMyMatterSubsystemRep().getAccelerationVars(s);
        return getMyRigidBodyNode().getAppliedBodyForce(av);
    }
    SpatialVec& updBodyAppliedForces(State& s) const {
        SBAccelerationVars& av = getMyMatterSubsystemRep().updAccelerationVars(s);
        return getMyRigidBodyNode().updAppliedBodyForce(av);
    }
    const SpatialVec& getBodyAcceleration(const State& s) const {
        const SBAccelerationCache& ac = getMyMatterSubsystemRep().getAccelerationCache(s);
        return getMyRigidBodyNode().getA_GB(ac);
    }
    const Transform& getMobilizerTransform(const State& s) const {
        const SBPositionCache& pc = getMyMatterSubsystemRep().getPositionCache(s);
        return getMyRigidBodyNode().getX_MbM(pc);
    }
    const SpatialVec& getMobilizerVelocity(const State& s) const {
        const SBVelocityCache& vc = getMyMatterSubsystemRep().getVelocityCache(s);
        return getMyRigidBodyNode().getV_MbM(vc);
    }

    const RigidBodyNode& realizeTopology(int& nxtU, int& nxtUSq, int& nxtQ) const {
        delete myRBnode;
        myRBnode = createRigidBodyNode(nxtU,nxtUSq,nxtQ);

        int level;
        if (!myParentId.isValid()) {
            // this is ground
            assert(myMobilizedBodyId == 0);
            level = 0;
        } else {
            // not ground
            const MobilizedBodyRep& parent = 
                myMatterSubsystemRep->getMobilizedBody(myParentId).getRep();
            level = parent.myRBnode->getLevel() + 1;
            parent.myRBnode->addChild(myRBnode);
            myRBnode->setParent(parent.myRBnode);
        }

        myRBnode->setLevel(level);
        myRBnode->setNodeNum(myMobilizedBodyId);
        return *myRBnode;
    }

    void invalidateTopologyCache() const {
        delete myRBnode; myRBnode=0;
        if (myMatterSubsystemRep)
            myMatterSubsystemRep->invalidateSubsystemTopologyCache();
    }

    const RigidBodyNode& getMyRigidBodyNode() const {
        SimTK_ASSERT(myRBnode && myMatterSubsystemRep && 
                     myMatterSubsystemRep->subsystemTopologyHasBeenRealized(),
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

    bool isInSubsystem() const {return myMatterSubsystemRep != 0;}

    void setMyMatterSubsystem(SimbodyMatterSubsystem& matter,
                              MobilizedBodyId  parentId,
                              MobilizedBodyId  id)
    {
        assert(!myMatterSubsystemRep);
        myMatterSubsystemRep = &matter.updRep();
        myParentId           = parentId;
        myMobilizedBodyId    = id;
    }

    void setMyHandle(MobilizedBody& h) {myHandle = &h;}
    const MobilizedBody& getMyHandle() const {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}

private:
    friend class MobilizedBody;
    MobilizedBody* myHandle;	// the owner handle of this rep

        // TOPOLOGY "STATE"

    // Base class topological properties. Derived MobilizedBodies may
    // have further topological properties. Whenever these change, be
    // sure to set topologyRealized=false!
    Body theBody;
    Transform defaultInboardFrame;  // default for Mb (in Parent frame)
    Transform defaultOutboardFrame; // default for M (in Body frame)

    // These data members are filled in once the MobilizedBody is added to
    // a MatterSubsystem. Note that this pointer is just a reference to
    // the Subsystem which owns this MobilizedBody -- we don't own the
    // heap space here so don't delete it!
    SimbodyMatterSubsystemRep* myMatterSubsystemRep;
    MobilizedBodyId  myMobilizedBodyId; // id within the subsystem
    MobilizedBodyId  myParentId;

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

class MobilizedBody::Pin::PinRep : public MobilizedBody::MobilizedBodyRep {
public:
    PinRep() : defaultQ(0) { }
    PinRep* clone() const { return new PinRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(PinRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::Pin;
    Real defaultQ;  // the angle in radians
};


class MobilizedBody::Slider::SliderRep : public MobilizedBody::MobilizedBodyRep {
public:
    SliderRep() : defaultQ(0) { }
    SliderRep* clone() const { return new SliderRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(SliderRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::Slider;
    Real defaultQ;  // the displacement
};

class MobilizedBody::Universal::UniversalRep : public MobilizedBody::MobilizedBodyRep {
public:
    UniversalRep() : defaultQ(0) { }
    UniversalRep* clone() const { return new UniversalRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(UniversalRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::Universal;
    Vec2 defaultQ;  // the two angles
};

class MobilizedBody::Cylinder::CylinderRep : public MobilizedBody::MobilizedBodyRep {
public:
    CylinderRep() : defaultQ(0) { }
    CylinderRep* clone() const { return new CylinderRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(CylinderRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::Cylinder;
    Vec2 defaultQ;  // the angle in radians followed by the displacement
};

class MobilizedBody::BendStretch::BendStretchRep : public MobilizedBody::MobilizedBodyRep {
public:
    BendStretchRep() : defaultQ(0) { }
    BendStretchRep* clone() const { return new BendStretchRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(BendStretchRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::BendStretch;
    Vec2 defaultQ;  // the angle in radians followed by the displacement
};

class MobilizedBody::Planar::PlanarRep : public MobilizedBody::MobilizedBodyRep {
public:
    PlanarRep() : defaultQ(0) { }
    PlanarRep* clone() const { return new PlanarRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(PlanarRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::Planar;
    Vec3 defaultQ;  // the angle in radians followed by the two displacements
};

class MobilizedBody::Gimbal::GimbalRep : public MobilizedBody::MobilizedBodyRep {
public:
    GimbalRep() : defaultQ(0) { }
    GimbalRep* clone() const { return new GimbalRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(GimbalRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::Gimbal;
    Vec3 defaultQ;  // the three angles in radians
};

class MobilizedBody::Ball::BallRep : public MobilizedBody::MobilizedBodyRep {
public:
    BallRep() : defaultQ() { } // default is (1,0,0,0), the identity rotation
    BallRep* clone() const { return new BallRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(BallRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::Ball;
    Quaternion defaultQ;  // the default orientation
};

class MobilizedBody::Translation::TranslationRep : public MobilizedBody::MobilizedBodyRep {
public:
    TranslationRep() : defaultQ(0) { }
    TranslationRep* clone() const { return new TranslationRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(TranslationRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::Translation;
    Vec3 defaultQ;  // the three displacements
};

class MobilizedBody::Free::FreeRep : public MobilizedBody::MobilizedBodyRep {
public:
    FreeRep() : defaultQOrientation(), defaultQTranslation(0) { }
    FreeRep* clone() const { return new FreeRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(FreeRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::Free;
    Quaternion defaultQOrientation;
    Vec3       defaultQTranslation;
};

class MobilizedBody::LineOrientation::LineOrientationRep : public MobilizedBody::MobilizedBodyRep {
public:
    LineOrientationRep() : defaultQ() { } // 1,0,0,0
    LineOrientationRep* clone() const { return new LineOrientationRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(LineOrientationRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::LineOrientation;
    Quaternion defaultQ;
};

class MobilizedBody::FreeLine::FreeLineRep : public MobilizedBody::MobilizedBodyRep {
public:
    FreeLineRep() : defaultQOrientation(), defaultQTranslation(0) { }
    FreeLineRep* clone() const { return new FreeLineRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(FreeLineRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::FreeLine;
    Quaternion defaultQOrientation;
    Vec3       defaultQTranslation;
};

class MobilizedBody::Weld::WeldRep : public MobilizedBody::MobilizedBodyRep {
public:
    WeldRep() { }
    WeldRep* clone() const { return new WeldRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(WeldRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::Weld;
    // no q's
};

class MobilizedBody::Ground::GroundRep : public MobilizedBody::MobilizedBodyRep {
public:
    GroundRep() { }
    GroundRep* clone() const { return new GroundRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(GroundRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::Ground;
    // no q's
};

class MobilizedBody::Screw::ScrewRep : public MobilizedBody::MobilizedBodyRep {
public:
    explicit ScrewRep(Real p) : defaultPitch(p), defaultQ(0) { }

    Real getDefaultPitch() const {return defaultPitch;}
    void setDefaultPitch(Real p) {defaultPitch=p;}

    ScrewRep* clone() const { return new ScrewRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(ScrewRep, MobilizedBodyRep);
private:
    friend class MobilizedBody::Screw;
    Real defaultPitch;
    Real defaultQ; // the angle in radians
};

class MobilizedBody::Custom::CustomRep : public MobilizedBody::MobilizedBodyRep {
public:
    CustomRep(int nMobilities, int nCoordinates) 
      : nu(nMobilities), nq(nCoordinates), defaultQ(nq) {
        defaultQ.setToZero();
    }
    CustomRep* clone() const { return new CustomRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(CustomRep, MobilizedBodyRep);
private:
    int nu;
    int nq;
    Vector defaultQ;
};

    /////////////////////
    // CONSTRAINT REPS //
    /////////////////////

class Constraint::ConstraintRep {
public:
    ConstraintRep() : myHandle(0), myMatterSubsystemRep(0), myConstraintNode(0) {
    }
    virtual ~ConstraintRep() {
        delete myConstraintNode;
    }
    virtual ConstraintRep* clone() const = 0;

    // This creates a constraint node using the appropriate constraint type.
    virtual ConstraintNode* createConstraintNode() const = 0; 

    void realizeTopology(int& nxtQErr, int& nxtUErr, int& nxtMult) const {
        delete myConstraintNode;
        myConstraintNode = createConstraintNode();
        myConstraintNode->setConstraintNum(myConstraintId);
        myConstraintNode->setQErrIndex(nxtQErr);
        myConstraintNode->setUErrIndex(nxtUErr);
        myConstraintNode->setMultIndex(nxtMult);
        const int nConsEqns = myConstraintNode->getNConstraintEquations();
        nxtQErr += nConsEqns;
        nxtUErr += nConsEqns;
        nxtMult += nConsEqns;
        myConstraintNode->finishConstruction(*myMatterSubsystemRep);
    }

    void invalidateTopologyCache() const {
        delete myConstraintNode; myConstraintNode=0;
        if (myMatterSubsystemRep)
            myMatterSubsystemRep->invalidateSubsystemTopologyCache();
    }

    const ConstraintNode& getMyConstraintNode() const {
        SimTK_ASSERT(myConstraintNode && myMatterSubsystemRep && 
                     myMatterSubsystemRep->subsystemTopologyHasBeenRealized(),
          "An operation on a Constraint was illegal because realizeTopology() has "
          "not been performed on the containing Subsystem since the last topological change."
        );
        return *myConstraintNode;
    }

    void setMyMatterSubsystem(SimbodyMatterSubsystem& matter,
                              ConstraintId  id)
    {
        assert(!isInSubsystem());
        myMatterSubsystemRep = &matter.updRep();
        myConstraintId = id;
    }

    bool isInSubsystem() const {
        return myMatterSubsystemRep != 0;
    }

    const SimbodyMatterSubsystemRep& getMyMatterSubsystemRep() const {
        SimTK_ASSERT(myMatterSubsystemRep,
            "Operation illegal on a Constraint that is not in a Subsystem.");
        return *myMatterSubsystemRep;
    }
    SimbodyMatterSubsystemRep& updMyMatterSubsystemRep() {
        SimTK_ASSERT(myMatterSubsystemRep,
            "Operation illegal on a Constraint that is not in a Subsystem.");
        return *myMatterSubsystemRep;
    }

    void setMyHandle(Constraint& h) {myHandle = &h;}
    const Constraint& getMyHandle() const {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}

private:
    friend class Constraint;
    Constraint* myHandle;	// the owner handle of this rep

        // TOPOLOGY "STATE"

    // These data members are filled in once the Constraint is added to
    // a MatterSubsystem.
    SimbodyMatterSubsystemRep* myMatterSubsystemRep;
    ConstraintId               myConstraintId; // id within the subsystem

        // TOPOLOGY "CACHE"
    // A ConstraintNode is created for each Constraint during realizeTopology().
    // Think of it as the *computational* form of the Constraint; the Constraint
    // itself exists to make a nice API for the programmer. In fact is is possible
    // for several different Constraint classes to use the same underling computational
    // form while providing a different API.
    //
    // This is a pointer to an abstract object whose heap space is *owned* by
    // this Constraint. Be sure to delete it upon destruction and whenever
    // topology is re-realized.
    mutable ConstraintNode* myConstraintNode;
};


class Constraint::Rod::RodRep : public Constraint::ConstraintRep {
public:
    RodRep() : defaultPoint1(0), defaultPoint2(0), defaultRodLength(1) { }
    RodRep* clone() const { return new RodRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    SimTK_DOWNCAST(RodRep, ConstraintRep);
private:
    friend class Constraint::Rod;

    MobilizedBodyId body1; // B1
    MobilizedBodyId body2; // B2
    Vec3            defaultPoint1; // on body 1, exp. in B1 frame
    Vec3            defaultPoint2; // on body 2, exp. in B2 frame
    Real            defaultRodLength;
};

class Constraint::Ball::BallRep : public Constraint::ConstraintRep {
public:
    BallRep() : defaultPoint1(0), defaultPoint2(0) { }
    BallRep* clone() const { return new BallRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    SimTK_DOWNCAST(BallRep, ConstraintRep);
private:
    friend class Constraint::Ball;

    MobilizedBodyId body1; // B1
    MobilizedBodyId body2; // B2
    Vec3            defaultPoint1; // on body 1, exp. in B1 frame
    Vec3            defaultPoint2; // on body 2, exp. in B2 frame
};



class Constraint::Weld::WeldRep : public Constraint::ConstraintRep {
public:
    WeldRep() { } // transforms are Identity
    WeldRep* clone() const { return new WeldRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    SimTK_DOWNCAST(WeldRep, ConstraintRep);
private:
    friend class Constraint::Weld;

    MobilizedBodyId body1; // B1
    MobilizedBodyId body2; // B2
    Transform       defaultFrame1; // on body 1, relative to B1 frame
    Transform       defaultFrame2; // on body 2, relative to B2 frame};
};


class Constraint::Custom::CustomRep : public Constraint::ConstraintRep {
public:
    CustomRep* clone() const { return new CustomRep(*this); }

    ConstraintNode* createConstraintNode() const; 

    SimTK_DOWNCAST(CustomRep, ConstraintRep);
private:
    friend class Constraint::Custom;

    // TODO: God only knows what goes here!
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_REP_H_