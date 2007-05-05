#ifndef SimTK_SIMBODY_MOBILIZER_REP_H_
#define SimTK_SIMBODY_MOBILIZER_REP_H_



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
 * Private implementation of Mobilizer, and its included subclasses which
 * represent the built-in mobilizer types.
 */

#include "SimTKsimbody.h"
#include "simbody/internal/Mobilizer.h"

class RigidBodyNode;

namespace SimTK {

class Mobilizer::MobilizerRep {
public:
    MobilizerRep() : myHandle(0) { }
    virtual ~MobilizerRep() { }
    virtual MobilizerRep* clone() const = 0;

    // This creates a rigid body node using the appropriate mobilizer type.
    virtual RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,            // mass properties in body frame
        const Transform&         X_PMb,        // parent's attachment frame for this mobilizer
        const Transform&         X_BM,         // inboard mobilizer frame M in body frame
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const = 0; 

    void setMyHandle(Mobilizer& h) {myHandle = &h;}
    const Mobilizer& getMyHandle() const {assert(myHandle); return *myHandle;}
    void clearMyHandle() {myHandle=0;}
private:
    friend class Mobilizer;
    Mobilizer* myHandle;	// the owner handle of this rep
};

class Mobilizer::Pin::PinRep : public Mobilizer::MobilizerRep {
public:
    PinRep* clone() const { return new PinRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(PinRep, MobilizerRep);
};


class Mobilizer::Slider::SliderRep : public Mobilizer::MobilizerRep {
public:
    SliderRep* clone() const { return new SliderRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(SliderRep, MobilizerRep);
};

class Mobilizer::Universal::UniversalRep : public Mobilizer::MobilizerRep {
public:
    UniversalRep* clone() const { return new UniversalRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(UniversalRep, MobilizerRep);
};

class Mobilizer::Cylinder::CylinderRep : public Mobilizer::MobilizerRep {
public:
    CylinderRep* clone() const { return new CylinderRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(CylinderRep, MobilizerRep);
};

class Mobilizer::BendStretch::BendStretchRep : public Mobilizer::MobilizerRep {
public:
    BendStretchRep* clone() const { return new BendStretchRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(BendStretchRep, MobilizerRep);
};

class Mobilizer::Planar::PlanarRep : public Mobilizer::MobilizerRep {
public:
    PlanarRep* clone() const { return new PlanarRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(PlanarRep, MobilizerRep);
};

class Mobilizer::Gimbal::GimbalRep : public Mobilizer::MobilizerRep {
public:
    GimbalRep* clone() const { return new GimbalRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(GimbalRep, MobilizerRep);
};

class Mobilizer::Orientation::OrientationRep : public Mobilizer::MobilizerRep {
public:
    OrientationRep* clone() const { return new OrientationRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(OrientationRep, MobilizerRep);
};

class Mobilizer::Translation::TranslationRep : public Mobilizer::MobilizerRep {
public:
    TranslationRep* clone() const { return new TranslationRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(TranslationRep, MobilizerRep);
};

class Mobilizer::Free::FreeRep : public Mobilizer::MobilizerRep {
public:
    FreeRep* clone() const { return new FreeRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(FreeRep, MobilizerRep);
};

class Mobilizer::LineOrientation::LineOrientationRep : public Mobilizer::MobilizerRep {
public:
    LineOrientationRep* clone() const { return new LineOrientationRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(LineOrientationRep, MobilizerRep);
};

class Mobilizer::FreeLine::FreeLineRep : public Mobilizer::MobilizerRep {
public:
    FreeLineRep* clone() const { return new FreeLineRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(FreeLineRep, MobilizerRep);
};

class Mobilizer::Weld::WeldRep : public Mobilizer::MobilizerRep {
public:
    WeldRep* clone() const { return new WeldRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(WeldRep, MobilizerRep);
};


class Mobilizer::Screw::ScrewRep : public Mobilizer::MobilizerRep {
public:
    ScrewRep(Real p) : pitch(p) { }

    ScrewRep* clone() const { return new ScrewRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(ScrewRep, MobilizerRep);
private:
    Real pitch;
};

class Mobilizer::User::UserRep : public Mobilizer::MobilizerRep {
public:
    UserRep* clone() const { return new UserRep(*this); }

    RigidBodyNode* createRigidBodyNode(
        const MassProperties&    m,
        const Transform&         X_PMb,
        const Transform&         X_BM,
        int&                     nxtU,
        int&                     nxtUSq,
        int&                     nxtQ) const;

    SimTK_DOWNCAST(UserRep, MobilizerRep);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZER_REP_H_