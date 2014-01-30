/* -------------------------------------------------------------------------- *
 *                       Simbody(tm) - Tim's Box PGS                          *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2012-14 Stanford University and the Authors.        *
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

/* Solve TimsBox contact & impact using the Projected Gauss Seidel iterative
solver rather than a direct solver. */

//#define NDEBUG 1

#include "Simbody.h"

#include <string>
#include <iostream>
#include <exception>

using std::cout;
using std::endl;

using namespace SimTK;

//#define USE_TIMS_PARAMS
#define ANIMATE // off to get more accurate CPU time (you can still playback)

//#define HERTZ
#define POISSON
#define NEWTON

//==============================================================================
//                           MY CONTACT ELEMENT
//==============================================================================
// This abstract class hides the details about which kind of contact constraint
// we're dealing with, while giving us enough to work with for deciding what's
// on and off and generating impulses.
//
// There is always a scalar associated with the constraint for making 
// decisions. There may be a friction element associated with this contact.
namespace {
class MyFrictionElement;
class MyContactElement {
public:
    MyContactElement(Constraint uni, Real multSign, Real coefRest) 
    :   m_uni(uni), m_multSign(multSign), m_coefRest(coefRest), 
        m_index(-1), m_friction(0),
        m_velocityDependentCOR(NaN), m_restitutionDone(false) 
    {   m_uni.setDisabledByDefault(true); }

    MultiplierIndex getMultIndex(const State& s) const {
        int mp, mv, ma;
        MultiplierIndex px0, vx0, ax0;
        m_uni.getNumConstraintEquationsInUse(s,mp,mv,ma);
        assert(mp==1 && mv==0 && ma==0); // don't call if not enabled
        m_uni.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
        assert(px0.isValid() && !vx0.isValid() && !ax0.isValid());
        return px0;
    }

    virtual ~MyContactElement() {}
    
    // (Re)initialize base & concrete class. If overridden, be sure to
    // invoke base class first.
    virtual void initialize() {
        setRestitutionDone(false); 
        m_velocityDependentCOR = NaN;
    }

    // Provide a human-readable string identifying the type of contact
    // constraint.
    virtual String getContactType() const = 0;

    // These must be constructed so that a negative value means the 
    // unilateral constraint condition is violated.
    virtual Real getPerr(const State& state) const = 0;
    virtual Real getVerr(const State& state) const = 0;
    virtual Real getAerr(const State& state) const = 0;

    // This returns a point in the ground frame at which you might want to
    // say the constraint is "located", for purposes of display. This should
    // return something useful even if the constraint is currently off.
    virtual Vec3 whereToDisplay(const State& state) const = 0;

    // This is used by some constraints to collect position information that
    // may be used later to set instance variables when enabling the underlying
    // Simbody constraint. All constraints zero impulses here.
    virtual void initializeForImpact(const State& state, Real captureVelocity) { 
        if (-captureVelocity <= getVerr(state) && getVerr(state) < 0) {
            m_velocityDependentCOR = 0;
            SimTK_DEBUG3("CAPTURING %d because %g <= v=%g < 0\n",
                m_index, -captureVelocity, getVerr(state));
        } else {
            m_velocityDependentCOR = m_coefRest;
        }
        
        setRestitutionDone(false);        
    }

    // Returns zero if the constraint is not currently enabled. Otherwise 
    // return the signed constraint force, with a negative value indicating
    // that the unilateral force condition is violated.
    Real getForce(const State& s) const {
        if (isDisabled(s)) return 0;
        return m_multSign*s.updMultipliers()[getMultIndex(s)];
    }

    // Append to geometry array.
    virtual void showContactForce(const State& s, 
                                  Array_<DecorativeGeometry>& geometry) const {}

    bool isProximal(const State& state, Real posTol) const
    {   return /*!isDisabled(state) || */getPerr(state) <= posTol; }
    bool isCandidate(const State& state, Real posTol, Real velTol) const
    {   return isProximal(state, posTol) && getVerr(state) <= velTol; }


    void enable(State& state) const {m_uni.enable(state);}
    void disable(State& state) const {m_uni.disable(state);}
    bool isDisabled(const State& state) const {return m_uni.isDisabled(state);}

    void setMyDesiredDeltaV(const State&    s,
                            Vector&         desiredDeltaV) const
    {   Vector myDesiredDV(1); myDesiredDV[0] = m_multSign*getVerr(s);
        m_uni.setMyPartInConstraintSpaceVector(s, myDesiredDV, 
                                                   desiredDeltaV); }

    Real getMyValueFromConstraintSpaceVector(const State& state,
                                             const Vector& lambda) const
    {   Vector myValue(1);
        m_uni.getMyPartFromConstraintSpaceVector(state, lambda, myValue);
        return -m_multSign*myValue[0]; }

    Real getMaxCoefRest() const {return m_coefRest;}
    Real getEffectiveCoefRest() const {return m_velocityDependentCOR;}
    void setRestitutionDone(bool isDone) {m_restitutionDone=isDone;}
    bool isRestitutionDone() const {return m_restitutionDone;}

    // Record position within the set of unilateral contact constraints.
    void setContactIndex(int index) {m_index=index;}
    int getContactIndex() const {return m_index;}
    // If there is a friction element for which this is the master contact,
    // record it here.
    void setFrictionElement(MyFrictionElement& friction)
    {   m_friction = &friction; }
    // Return true if there is a friction element associated with this contact
    // element.
    bool hasFrictionElement() const {return m_friction != 0;}
    // Get the associated friction element.
    const MyFrictionElement& getFrictionElement() const
    {   assert(hasFrictionElement()); return *m_friction; }
    MyFrictionElement& updFrictionElement() const
    {   assert(hasFrictionElement()); return *m_friction; }

protected:
    Constraint          m_uni;
    const Real          m_multSign; // 1 or -1
    const Real          m_coefRest;

    int                 m_index; // contact index in unilateral constraint set
    MyFrictionElement*  m_friction; // if any (just a reference, not owned)

    // Runtime -- initialized at start of impact handler.
    Real m_velocityDependentCOR; // Calculated at start of impact 
    bool m_restitutionDone;
};



//==============================================================================
//                           MY FRICTION ELEMENT
//==============================================================================
// Generated forces during sliding, and the force limit during rolling, depend 
// on a scalar normal force N that comes from a 
// separate "normal force master", which might be one of the following:
//  - a unilateral constraint
//  - a bilateral constraint 
//  - a mobilizer
//  - a compliant force element 
// If the master is an inactive unilateral constraint, or if N=0, then no 
// friction forces are generated. In this example, we're only going to use
// a unilateral contact constraint as the "normal force master".
//
// For all but the compliant normal force master, the normal force N is 
// acceleration-dependent and thus may be coupled to the force produced by a
// sliding friction element. This may require iteration to ensure consistency
// between the sliding friction force and its master contact's normal force.
//
// A Coulomb friction element depends on a scalar slip speed defined by the
// normal force master (this might be the magnitude of a generalized speed or
// slip velocity vector). When the slip velocity goes to zero, the stiction 
// constraint is enabled if its constraint force magnitude can be kept to
// mu_s*|N| or less. Otherwise, or if the slip velocity is nonzero, the sliding
// force is enabled instead and generates a force of constant magnitude mu_d*|N| 
// that opposes the slip direction, or impending slip direction, as defined by 
// the master.
class MyFrictionElement {
public:
    MyFrictionElement(Real mu_d, Real mu_s, Real mu_v)
    :   mu_d(mu_d), mu_s(mu_s), mu_v(mu_v), m_index(-1) {}

    virtual ~MyFrictionElement() {}

    // (Re)initialize base & concrete class. If overridden, be sure to
    // invoke base class first.
    virtual void initialize() {
    }

    Real getDynamicFrictionCoef() const {return mu_d;}
    Real getStaticFrictionCoef()  const {return mu_s;}
    Real getViscousFrictionCoef() const {return mu_v;}

    // This returns a point in the ground frame at which you might want to
    // say the friction is "located", for purposes of display.
    virtual Vec3 whereToDisplay(const State& state) const = 0;

    // Return true if the stiction constraint is enabled.
    virtual bool isEnabled(const State&) const = 0;

    virtual void setInstanceParameters(State&) const {}
    virtual void enable(State&) const = 0;
    virtual void disable(State&) const = 0;

    // Return true if the normal force master *could* be involved in an 
    // impact event (because it is touching).
    virtual bool isMasterProximal(const State&, Real posTol) const = 0;

    // Return true if the normal force master is currently generating a
    // normal force (or impulse) so that this friction element might be 
    // generating a force also.
    virtual bool isMasterActive(const State&) const = 0;


    // This is used by some stiction constraints to collect position information
    // that may be used later to set instance variables when enabling the 
    // underlying Simbody constraint. Recorded impulses should be zeroed.
    virtual void initializeFriction(const State& state) = 0; 

    // If this friction element's stiction constraint is enabled, set its
    // constraint-space velocity entry(s) in desiredDeltaV to the current
    // slip velocity (which might be a scalar or 2-vector).
    virtual void setMyDesiredDeltaV(const State& s,
                                    Vector&      desiredDeltaV) const = 0;

    // Output the status, friction force, slip velocity, prev slip direction
    // (scalar or vector) to the given ostream, indented as indicated and 
    // followed by a newline. May generate multiple lines.
    virtual std::ostream& writeFrictionInfo(const State& state,
                                            const String& indent,
                                            std::ostream& o) const = 0;

    // Optional: give some kind of visual representation for the friction force.
    virtual void showFrictionForce(const State& state, 
        Array_<DecorativeGeometry>& geometry, const Vec3& color) const {}


    void setFrictionIndex(int index) {m_index=index;}
    int getFrictionIndex() const {return m_index;}

private:
    Real mu_d, mu_s, mu_v;
    int  m_index; // friction index within unilateral constraint set
};



//==============================================================================
//                       MY UNILATERAL CONSTRAINT SET
//==============================================================================

// These are indices into the unilateral constraint set arrays.
struct MyElementSubset {
    void clear() {m_contact.clear();m_friction.clear();}
    Array_<int> m_contact;
    Array_<int> m_friction; // friction elements that might stick
};

class MyUnilateralConstraintSet {
public:
    // Capture velocity: if the normal approach velocity
    // is smaller, the coefficient of restitution is set to zero for the 
    // upcoming impact. Transition velocity: if a slip velocity is smaller than 
    // this use the static coefficient of friction, otherwise use dynamic
    // plus viscous.
    MyUnilateralConstraintSet(const MultibodySystem& mbs, 
                              Real captureVelocity, Real transitionVelocity)
    :   m_mbs(mbs), m_captureVelocity(captureVelocity),
        m_transitionVelocity(transitionVelocity)  {}

    // This class takes over ownership of the heap-allocated contact element.
    int addContactElement(MyContactElement* contact) {
        const int index = (int)m_contact.size();
        m_contact.push_back(contact);
        contact->setContactIndex(index);
        return index;
    }
    // This class takes over ownership of the heap-allocated friction element.
    int addFrictionElement(MyFrictionElement* friction) {
        const int index = (int)m_friction.size();
        m_friction.push_back(friction);
        friction->setFrictionIndex(index);
        return index;
    }

    Real getCaptureVelocity() const {return m_captureVelocity;}
    void setCaptureVelocity(Real v) {m_captureVelocity=v;}
    Real getTransitionVelocity() const {return m_transitionVelocity;}
    void setTransitionVelocity(Real v) {m_transitionVelocity=v;}

    int getNumContactElements() const {return (int)m_contact.size();}
    int getNumFrictionElements() const {return (int)m_friction.size();}
    const MyContactElement& getContactElement(int ix) const 
    {   return *m_contact[ix]; }
    const MyFrictionElement& getFrictionElement(int ix) const 
    {   return *m_friction[ix]; }

    // Allow writable access to elements from const set so we can record
    // runtime results (e.g. impulses).
    MyContactElement&  updContactElement(int ix) const {return *m_contact[ix];}
    MyFrictionElement& updFrictionElement(int ix) const {return *m_friction[ix];}

    // Initialize all runtime fields in the contact & friction elements.
    void initialize()
    {
        for (unsigned i=0; i < m_contact.size(); ++i)
            m_contact[i]->initialize();
        for (unsigned i=0; i < m_friction.size(); ++i)
            m_friction[i]->initialize();
    }

    // Return the contact and friction elements that might be involved in an
    // impact occurring in this configuration. They are the contact elements 
    // for which perr <= posTol, and friction elements whose normal force 
    // masters can be involved in the impact. State must be realized through 
    // Position stage.
    void findProximalElements(const State&      state,
                              Real              posTol,
                              MyElementSubset&  proximals,
                              MyElementSubset&  distals) const
    {
        proximals.clear(); distals.clear();
        for (unsigned i=0; i < m_contact.size(); ++i)
            if (m_contact[i]->isProximal(state,posTol)) 
                proximals.m_contact.push_back(i);
            else distals.m_contact.push_back(i);
        for (unsigned i=0; i < m_friction.size(); ++i)
            if (m_friction[i]->isMasterProximal(state,posTol))
                proximals.m_friction.push_back(i);
            else distals.m_friction.push_back(i);
        // Any friction elements might stick if they are proximal since
        // we'll be changing velocities.
    }

    // In Debug mode, produce a useful summary of the current state of the
    // contact and friction elements.
    void showConstraintStatus(const State& state, const String& place) const;

    ~MyUnilateralConstraintSet() {
        for (unsigned i=0; i < m_contact.size(); ++i)
            delete m_contact[i];
        for (unsigned i=0; i < m_friction.size(); ++i)
            delete m_friction[i];
    }

    const MultibodySystem& getMultibodySystem() const {return m_mbs;}
private:
    const MultibodySystem&      m_mbs;
    Real                        m_captureVelocity;
    Real                        m_transitionVelocity;
    Array_<MyContactElement*>   m_contact;
    Array_<MyFrictionElement*>  m_friction;
};


//==============================================================================
//                             MY POINT CONTACT
//==============================================================================
// Define a unilateral constraint to represent contact of a follower point on 
// one body with a plane fixed to another body.
class MyPointContact : public MyContactElement {
    typedef MyContactElement Super;
public:
    MyPointContact(
        MobilizedBody& planeBodyB, const UnitVec3& normal_B, Real height,
        MobilizedBody& followerBodyF, const Vec3& point_F, 
        Real           coefRest)
    :   MyContactElement( 
             Constraint::PointInPlane(planeBodyB, normal_B, height,
                                      followerBodyF, point_F),
             Real(-1), // multiplier sign
             coefRest),
        m_planeBody(planeBodyB), m_frame(normal_B, ZAxis), m_height(height), 
        m_follower(followerBodyF), m_point(point_F)
    {
    }

    Real getPerr(const State& s) const OVERRIDE_11 {
        const Vec3 p = m_follower.findStationLocationInAnotherBody
                                    (s, m_point, m_planeBody);
        return ~p*m_frame.z() - m_height;
    }
    Real getVerr(const State& s) const OVERRIDE_11 {
        const Vec3 v = m_follower.findStationVelocityInAnotherBody
                                    (s, m_point, m_planeBody);
        return ~v*m_frame.z(); // normal is constant in P
    }
    Real getAerr(const State& s) const OVERRIDE_11 {
        const Vec3 a = m_follower.findStationAccelerationInAnotherBody
                                    (s, m_point, m_planeBody);
        return ~a*m_frame.z(); // normal is constant in P
    }

    String getContactType() const OVERRIDE_11 {return "Point";}
    Vec3 whereToDisplay(const State& state) const OVERRIDE_11 {
        return m_follower.findStationLocationInGround(state,m_point);
    }

    // Will be zero if the stiction constraints are on.
    Vec2 getSlipVelocity(const State& s) const {
        const Vec3 v = m_follower.findStationVelocityInAnotherBody
                                                    (s, m_point, m_planeBody);
        return Vec2(~v*m_frame.x(), ~v*m_frame.y());
    }
    // Will be zero if the stiction constraints are on.
    Vec2 getSlipAcceleration(const State& s) const {
        const Vec3 a = m_follower.findStationAccelerationInAnotherBody
                                                    (s, m_point, m_planeBody);
        return Vec2(~a*m_frame.x(), ~a*m_frame.y());
    }

    Vec3 getContactPointInPlaneBody(const State& s) const
    {   return m_follower.findStationLocationInAnotherBody
                                                    (s, m_point, m_planeBody); }

    void showContactForce(const State& s, 
                          Array_<DecorativeGeometry>& geometry)
            const OVERRIDE_11
    {
        const Real Scale = .1;
        const Real f = getForce(s);
        if (std::abs(f) < SignificantReal)
            return;
        const Vec3 stationG = whereToDisplay(s);
        const Vec3 endG = stationG + Scale*f*m_frame.z();
        geometry.push_back(DecorativeLine(endG     + Vec3(0,.05,0),
                                          stationG + Vec3(0,.05,0))
                            .setColor(Magenta));
    }

    const MobilizedBody& getBody() const {return m_follower;}
    MobilizedBody& updBody() {return m_follower;}
    const Vec3& getBodyStation() const {return m_point;}

    const MobilizedBody& getPlaneBody() const  {return m_planeBody;}
    MobilizedBody& updPlaneBody() const {return m_planeBody;}
    const Rotation& getPlaneFrame() const {return m_frame;}
    Real getPlaneHeight() const {return m_height;}

private:
    MobilizedBody&      m_planeBody;    // body P
    const Rotation      m_frame;        // z is normal; expressed in P
    const Real          m_height;

    MobilizedBody&      m_follower;     // body F
    const Vec3          m_point;        // measured & expressed in F
};


//==============================================================================
//                        MY POINT CONTACT FRICTION
//==============================================================================
// This friction element expects its master to be a unilateral point contact 
// constraint. It provides slipping forces or stiction constraint forces acting
// in the plane, based on the normal force being applied by the point contact 
// constraint.
class MyPointContactFriction : public MyFrictionElement {
    typedef MyFrictionElement Super;
public:
    // The constructor allocates two NoSlip1D constraints.
    MyPointContactFriction(MyPointContact& contact,
        Real mu_d, Real mu_s, Real mu_v, Real vtol, //TODO: shouldn't go here
        GeneralForceSubsystem& forces)
    :   MyFrictionElement(mu_d,mu_s,mu_v), m_contact(contact),
        m_noslipX(contact.updPlaneBody(), Vec3(0), contact.getPlaneFrame().x(), 
                  contact.updPlaneBody(), contact.updBody()),
        m_noslipY(contact.updPlaneBody(), Vec3(0), contact.getPlaneFrame().y(), 
                  contact.updPlaneBody(), contact.updBody())
    {
        assert((0 <= mu_d && mu_d <= mu_s) && (0 <= mu_v));
        contact.setFrictionElement(*this);
        m_noslipX.setDisabledByDefault(true);
        m_noslipY.setDisabledByDefault(true);
        initializeRuntimeFields();
    }

    ~MyPointContactFriction() {}

    void initialize() OVERRIDE_11 {
        Super::initialize();
        initializeRuntimeFields();
    }

    Vec3 whereToDisplay(const State& state) const OVERRIDE_11 {
        return m_contact.whereToDisplay(state);
    }


    MultiplierIndex getMultIndexX(const State& s) const {
        int mp, mv, ma;
        MultiplierIndex px0, vx0, ax0;
        m_noslipX.getNumConstraintEquationsInUse(s,mp,mv,ma);
        assert(mp==0 && mv==1 && ma==0); // don't call if not enabled
        m_noslipX.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
        assert(!px0.isValid() && vx0.isValid() && !ax0.isValid());
        return vx0;
    }

    MultiplierIndex getMultIndexY(const State& s) const {
        int mp, mv, ma;
        MultiplierIndex px0, vx0, ax0;
        m_noslipY.getNumConstraintEquationsInUse(s,mp,mv,ma);
        assert(mp==0 && mv==1 && ma==0); // don't call if not enabled
        m_noslipY.getIndexOfMultipliersInUse(s, px0, vx0, ax0);
        assert(!px0.isValid() && vx0.isValid() && !ax0.isValid());
        return vx0;
    }

    // The way we constructed the NoSlip1D constraints makes the multipliers be
    // the force on Ground; we negate here so we'll get the force on the sliding
    // body instead.
    Vec2 getFrictionForce(const State& s) const {
        if (m_noslipX.isDisabled(s)) return Vec2(0);
        Vec2 fOnG(s.updMultipliers()[getMultIndexX(s)],
                  s.updMultipliers()[getMultIndexY(s)]);
        return -fOnG;
    }

    // Implement pure virtuals from MyFrictionElement base class.

    bool isEnabled(const State& s) const OVERRIDE_11
    {   return !m_noslipX.isDisabled(s); } // X,Z always on or off together

    // Note that initializeForStiction() must have been called first.
    void setInstanceParameters(State& s) const OVERRIDE_11
    {   m_noslipX.setContactPoint(s, m_contactPointInPlane);
        m_noslipY.setContactPoint(s, m_contactPointInPlane); }

    void enable(State& s) const OVERRIDE_11
    {   m_noslipX.setContactPoint(s, m_contactPointInPlane);
        m_noslipY.setContactPoint(s, m_contactPointInPlane);
        m_noslipX.enable(s); m_noslipY.enable(s); }

    void disable(State& s) const OVERRIDE_11
    {   m_noslipX.disable(s); m_noslipY.disable(s); }

    bool isMasterProximal(const State& s, Real posTol) const OVERRIDE_11
    {   return m_contact.isProximal(s, posTol); }

    bool isMasterActive(const State& s) const OVERRIDE_11
    {   return !m_contact.isDisabled(s); }


    // Set the friction application point to be the projection of the contact 
    // point onto the contact plane. This will be used the next time friction
    // is enabled. Requires state realized to Position stage.
    void initializeFriction(const State& s) OVERRIDE_11 {
        const Vec3 p = m_contact.getContactPointInPlaneBody(s);
        m_contactPointInPlane = p;
    }

    // Fill in deltaV to eliminate slip velocity using the stiction 
    // constraints.
    void setMyDesiredDeltaV(const State& s,
                            Vector& desiredDeltaV) const OVERRIDE_11
    {
        if (!isEnabled(s)) return;

        const Vec2 dv = -m_contact.getSlipVelocity(s); // X,Z
        Vector myDesiredDV(1); // Nuke translational velocity also.
        myDesiredDV[0] = dv[0];
        m_noslipX.setMyPartInConstraintSpaceVector(s, myDesiredDV, desiredDeltaV);
        myDesiredDV[0] = dv[1];
        m_noslipY.setMyPartInConstraintSpaceVector(s, myDesiredDV, desiredDeltaV);
    }

    Real getMyImpulseMagnitudeFromConstraintSpaceVector(const State& state,
                                                        const Vector& lambda) const
    {   Vector myImpulseX(1), myImpulseY(1);
        m_noslipX.getMyPartFromConstraintSpaceVector(state, lambda, myImpulseX);
        m_noslipY.getMyPartFromConstraintSpaceVector(state, lambda, myImpulseY);
        const Vec2 tI(myImpulseX[0], myImpulseY[0]);
        return tI.norm();
    }


    std::ostream& writeFrictionInfo(const State& s, const String& indent, 
                                    std::ostream& o) const OVERRIDE_11 
    {
        o << indent;
        if (!isMasterActive(s)) o << "OFF";
        else if (isEnabled(s)) o << "STICK";
        else o << "SLIP";

        const Vec2 v = m_contact.getSlipVelocity(s);
        const Vec2 f = getFrictionForce(s);
        o << " V=" << ~v << " F=" << ~f << endl;
        return o;
    }


    void showFrictionForce(const State& s, Array_<DecorativeGeometry>& geometry,
                           const Vec3& color) 
            const OVERRIDE_11
    {
        const Real Scale = 0.1;
        const Vec2 f = getFrictionForce(s);
        if (f.normSqr() < square(SignificantReal))
            return;
        const MobilizedBody& bodyB = m_contact.getBody();
        const Vec3& stationB = m_contact.getBodyStation();
        const Vec3 stationG = bodyB.getBodyTransform(s)*stationB;
        Vec3 F = f[0]*m_contact.getPlaneFrame().x()
                 + f[1]*m_contact.getPlaneFrame().y();
        const Vec3 endG = stationG - Scale*F;
        geometry.push_back(DecorativeLine(endG     + Vec3(0,.05,0),
                                          stationG + Vec3(0,.05,0))
                            .setColor(color));
    }

    const MyPointContact& getMyPointContact() const {return m_contact;}

private:
    void initializeRuntimeFields() {
        m_contactPointInPlane = Vec3(0); 
    }
    const MyPointContact&   m_contact;

    Constraint::NoSlip1D    m_noslipX, m_noslipY; // stiction

    // Runtime
    Vec3 m_contactPointInPlane; // point on plane body where friction will act
};



//==============================================================================
//                            BODY WATCHER
//==============================================================================
// Prior to rendering each frame, point the camera at the given body's
// origin.
class BodyWatcher : public Visualizer::FrameController {
public:
    explicit BodyWatcher(const MobilizedBody& body) : m_body(body) {}

    void generateControls(const Visualizer&             viz, 
                          const State&                  state, 
                          Array_< DecorativeGeometry >& geometry) OVERRIDE_11
    {
        const Vec3 Bo = m_body.getBodyOriginLocation(state);
        const Vec3 p_GC = Bo + Vec3(0, 1.5, 7); // above and back
        const Rotation R_GC(UnitVec3(0,1,0), YAxis,
                            p_GC-Bo, ZAxis);
        viz.setCameraTransform(Transform(R_GC,p_GC));
        //viz.pointCameraAt(Bo, Vec3(0,1,0));
    }
private:
    const MobilizedBody m_body;
};

//==============================================================================
//                            MY PUSH FORCE
//==============================================================================
// This is a force element that generates a constant force on a body for a
// specified time interval. It is useful to perturb the system to force it
// to transition from sticking to sliding, for example.
class MyPushForceImpl : public Force::Custom::Implementation {
public:
    MyPushForceImpl(const MobilizedBody& bodyB, 
                    const Vec3&          stationB,
                    const Vec3&          forceG, // force direction in Ground!
                    Real                 onTime,
                    Real                 offTime)
    :   m_bodyB(bodyB), m_stationB(stationB), m_forceG(forceG),
        m_on(onTime), m_off(offTime)
    {    }

    //--------------------------------------------------------------------------
    //                       Custom force virtuals
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const
                   OVERRIDE_11
    {
        if (!(m_on <= state.getTime() && state.getTime() <= m_off))
            return;

        m_bodyB.applyForceToBodyPoint(state, m_stationB, m_forceG, bodyForces);
    }

    // No potential energy.
    Real calcPotentialEnergy(const State& state) const OVERRIDE_11 {return 0;}

    void calcDecorativeGeometryAndAppend
       (const State& state, Stage stage, 
        Array_<DecorativeGeometry>& geometry) const OVERRIDE_11
    {
        const Real ScaleFactor = 0.1;
        if (stage != Stage::Time) return;
        if (!(m_on <= state.getTime() && state.getTime() <= m_off))
            return;
        const Vec3 stationG = m_bodyB.findStationLocationInGround(state, m_stationB);
        geometry.push_back(DecorativeLine(stationG-ScaleFactor*m_forceG, stationG)
                            .setColor(Yellow)
                            .setLineThickness(3));
    }
private:
    const MobilizedBody& m_bodyB; 
    const Vec3           m_stationB;
    const Vec3           m_forceG;
    Real                 m_on;
    Real                 m_off;
};


//==============================================================================
//                       PGS AUGMENTED MULTIBODY SYSTEM
//==============================================================================
/* This is a Simbody MultibodySystem able to provide some additional information
about its constraints, in a form suitable for a PGS solver. The extra 
information allows us to emulate conditional constraints using Simbody's 
unconditional constraints plus its constraint enable/disable feature.

We first construct the system with all possible constraints included, but with
conditional ones initially disabled. Then for any given configuration q, we 
determine which conditional constraints are "proximal", meaning that they may 
participate in applying forces to the system in that configuration. Those 
constraints are enabled and all other conditional constraints are disabled. 
Simbody can then calculate the constraint Jacobian G(q), and the 
constraint-space compliance matrix A(q)=GM\~G.

Extra constraint info:
Unconditional: 
  - which multipliers / holonomic or not
Bounded scalar (stop/contact/clutch/motor):
  - lower and upper bounds
  - which multipliers / holonomic or not
  - if holo, restitution coefficients, capture velocity
Length-limited vector:
  - function L(t,q,u) giving max length
  - which multipliers / holonomic or not
Frictional constraint:
  - master constraint (unconditional or bounded, holonomic, dim 1,2, or 3)
  - friction coefficients, transition velocity
  - which multipliers
*/
const Real DefaultCaptureVelocity    = .01,
           DefaultTransitionVelocity = .01;
class PGSAugmentedMultibodySystem : public MultibodySystem {
public:
    PGSAugmentedMultibodySystem() : m_matter(0), m_forces(0), m_unis(0) {
        m_matter = new SimbodyMatterSubsystem(*this);
        m_forces = new GeneralForceSubsystem(*this);
        m_unis   = new MyUnilateralConstraintSet(*this, 
                        DefaultCaptureVelocity, DefaultTransitionVelocity);   
        m_matter->setShowDefaultGeometry(false);
    }

    virtual ~PGSAugmentedMultibodySystem() 
    {   delete m_unis; delete m_forces; delete m_matter; }

    virtual const MobilizedBody& getBodyToWatch() const
    {   return m_matter->getGround(); }

    virtual void calcInitialState(State& state) const = 0;

    const SimbodyMatterSubsystem& getMatterSubsystem() const {return *m_matter;}
    SimbodyMatterSubsystem& updMatterSubsystem() {return *m_matter;}
    
    const GeneralForceSubsystem& getForceSubsystem() const {return *m_forces;}
    GeneralForceSubsystem& updForceSubsystem() {return *m_forces;}

    const MyUnilateralConstraintSet& getUnis() const {return *m_unis;}
    MyUnilateralConstraintSet& updUnis() {return *m_unis;}

private:
    //TODO: this shouldn't require pointers.
    SimbodyMatterSubsystem*     m_matter;
    GeneralForceSubsystem*      m_forces;
    MyUnilateralConstraintSet*  m_unis;
};

//==============================================================================
//                           PGS TIME STEPPER
//==============================================================================

struct Bounded {
    Bounded(MultiplierIndex ix, Real lb, Real ub, Real effCOR) 
    :   m_ix(ix), m_lb(lb), m_ub(ub), m_effCOR(effCOR), m_hitBound(false) 
    {   assert(m_lb<=m_ub); assert(0<=m_effCOR && m_effCOR<=1); }
    MultiplierIndex  m_ix;        // which constraint multiplier
    Real             m_lb, m_ub;  // lower, upper bounds; lb <= ub
    Real             m_effCOR;       // velocity-dependent COR
    bool             m_hitBound;
};

struct LengthLimited {
    LengthLimited(const Array_<MultiplierIndex>& components, Real maxLength)
    :   m_components(components), m_maxLength(maxLength), m_hitLimit(false) 
    {   assert(m_components.size()<=3); assert(m_maxLength>=0); }
    Array_<MultiplierIndex> m_components;
    Real                    m_maxLength;
    bool                    m_hitLimit;
};

struct Frictional {
    Frictional(const Array_<MultiplierIndex>& frictionComponents, 
               const Array_<MultiplierIndex>& normalComponents,
               Real                           muEff)
    :   m_Fk(frictionComponents), m_Nk(normalComponents), 
        m_effMu(muEff), m_wasLimited(false) 
    {   assert(m_Fk.size()<=3 && m_Nk.size()<=3); assert(m_effMu>=0); }
    Array_<MultiplierIndex> m_Fk, m_Nk;
    Real                    m_effMu;
    bool                    m_wasLimited;
};

/**
**/
class PGSTimeStepper {
public:
    explicit PGSTimeStepper(const PGSAugmentedMultibodySystem& ambs)
    :   m_ambs(ambs), 
        m_PGSConvergenceTol(1e-6), m_PGSMaxIters(100), m_PGSSOR(1),
        m_accuracy(1e-2), m_consTol(1e-3), m_useNewton(false) 
    {   resetPGSStats(); }

    void setUseNewtonRestitution(bool useNewton) {m_useNewton=useNewton;}
    bool getUseNewtonRestitution() const {return m_useNewton;}

    /** Set integration accuracy; requires variable length steps. **/
    void setAccuracy(Real accuracy) {m_accuracy=accuracy;}
    /** Set the tolerance to which constraints must be satisfied. **/
    void setConstraintTol(Real consTol) {m_consTol=consTol;}

    Real getAccuracy() const {return m_accuracy;}
    Real getConstraintTol() const {return m_consTol;}

    void resetPGSStats() const {
        for (int i=0; i<3; ++i) {
            m_PGSNumCalls[i] = 0;  // mutable
            m_PGSNumIters[i] = 0;
            m_PGSNumFailures[i] = 0;
        }
    }
    long long getPGSNumCalls(int phase) const {return m_PGSNumCalls[phase];}
    long long getPGSNumIters(int phase) const {return m_PGSNumIters[phase];}
    long long getPGSNumFailures(int phase) const 
    {   return m_PGSNumFailures[phase]; }

    void setPGSConvergenceTol(Real tol) {m_PGSConvergenceTol=tol;}
    void setPGSMaxIters(int mx) {m_PGSMaxIters=mx;}
    void setPGSSOR(Real sor) {assert(0<=sor && sor<=2); m_PGSSOR = sor;}

    Real getPGSConvergenceTol() const {return m_PGSConvergenceTol;}
    int  getPGSMaxIters() const {return m_PGSMaxIters;}
    Real getPGSSOR() const {return m_PGSSOR;}

    /** Initialize the PGSTimeStepper's internally maintained state to a copy
    of the given state. **/
    void initialize(const State& initState);
    const State& getState() const {return m_state;}
    Real getTime() const {return m_state.getTime();}

    /** Advance to the indicated time in one or more steps. **/
    Integrator::SuccessfulStepStatus stepTo(Real time);

    /** Advance to the indicated time in one or more steps. **/
    Integrator::SuccessfulStepStatus stepTo2(Real time);

private:
    // Determine which constraints will be involved for this step.
    void findProximalConstraints(const State&);
    // Enable all proximal constraints, disable all distal constraints, 
    // reassigning multipliers if needed. Returns true if anything changed.
    bool enableProximalConstraints(State&);
    // Calculate velocity-dependent coefficients of restitution and friction
    // and apply combining rules for dissimilar materials.
    void calcContactCoefficients(const State&, const Vector& verr);

    // Easy if there are no constraints active.
    void takeUnconstrainedStep(State& s, Real h);

    // Adjust given verr to reflect Newton restitution. 
    bool applyNewtonRestitutionIfAny(const State&, Vector& verr) const;

    // Adjust given compression impulse to include Poisson restitution impulse.
    bool applyPoissonRestitutionIfAny(const State&, const Vector& verr,
                                      Vector& impulse) const;

    // This phase uses all the proximal constraints and should use a starting
    // guess for impulse saved from the last step if possible.
    bool doCompressionPhase(const State&, const Vector& eps,
                            Vector& compressionImpulse);
    // This phase uses all the proximal constraints, but we expect the result
    // to be zero unless expansion causes new violations.
    bool doExpansionPhase(const State&, const Vector& eps,
                          Vector& reactionImpulse);

    bool anyPositionErrorsViolated(const State&, const Vector& perr) const;

    // This phase uses only holonomic constraints, and zero is a good initial
    // guess for the (hopefully small) position correction.
    bool doPositionCorrectionPhase(const State&, const Vector& eps,
                                   Vector& positionImpulse);

private:
    // Returns true if it converges.
    bool projGaussSeidel(int phase, // for stats
                         const Matrix& A, const Vector& eps, Vector& pi, 
                         const Array_<int>&     all, 
                         const Array_<int>&     unconditional, 
                         Array_<Bounded>&       bounded,
                         Array_<LengthLimited>& lengthLimited,
                         Array_<Frictional>&    frictional) const;

    Real    m_PGSConvergenceTol;
    int     m_PGSMaxIters;
    Real    m_PGSSOR;

    mutable long long   m_PGSNumCalls[3]; // phases 0=comp, 1=exp, 2=position
    mutable long long   m_PGSNumIters[3];
    mutable long long   m_PGSNumFailures[3];

private:
    const PGSAugmentedMultibodySystem&  m_ambs;
    Real                                m_accuracy;
    Real                                m_consTol;
    bool                                m_useNewton;

    // Runtime data.
    State                               m_state;
    MyElementSubset                     m_proximals, m_distals;
    Array_<Bounded>                     m_bounded, m_boundedPos;
    Array_<LengthLimited>               m_lengthLimited;
    Array_<Frictional>                  m_frictional;
    Array_<MultiplierIndex>             m_all, m_allPos, m_uncond, m_uncondPos;
    Matrix                              m_GMInvGt; // G M\ ~G

friend class ShowContact;
};



//==============================================================================
//                            SHOW CONTACT
//==============================================================================
// For each visualization frame, generate ephemeral geometry to show just 
// during this frame, based on the current State.
class ShowContact : public DecorationGenerator {
public:
    explicit ShowContact(const PGSTimeStepper& ts) 
    :   m_ts(ts) {}

    void generateDecorations(const State&                state, 
                             Array_<DecorativeGeometry>& geometry) OVERRIDE_11
    {
        const PGSAugmentedMultibodySystem& ambs = m_ts.m_ambs;
        const MyUnilateralConstraintSet& unis = ambs.getUnis();
        ambs.realize(state, Stage::Dynamics);
        const Real E = ambs.calcEnergy(state);
        DecorativeText energy; energy.setIsScreenText(true);
        energy.setText("Energy: " + String(E, "%.3f"));
        geometry.push_back(energy);

        for (unsigned i=0; i < m_ts.m_proximals.m_contact.size(); ++i) {
            const int id = m_ts.m_proximals.m_contact[i];
            const MyContactElement& contact = unis.getContactElement(id);
            const Vec3 loc = contact.whereToDisplay(state);
            if (!contact.isDisabled(state)) {
                geometry.push_back(DecorativeSphere(.1)
                    .setTransform(loc)
                    .setColor(Red).setOpacity(.25));
                contact.showContactForce(state, geometry);
                String text;
                if (!contact.hasFrictionElement())
                    text = "-ENABLED";
                geometry.push_back(DecorativeText(String(i)+text)
                    .setColor(White).setScale(.1)
                    .setTransform(loc+Vec3(0,.04,0)));
            } else {
                geometry.push_back(DecorativeText(String(i))
                    .setColor(White).setScale(.1)
                    .setTransform(loc+Vec3(0,.02,0)));
            }
        }

        for (unsigned i=0; i < m_ts.m_proximals.m_friction.size(); ++i) {
            const int id = m_ts.m_proximals.m_friction[i];
            const MyFrictionElement& felt = unis.getFrictionElement(id);
            const Vec3 loc = felt.whereToDisplay(state);
            const Frictional& fric = m_ts.m_frictional[i];
            String text = fric.m_wasLimited ? "slip" : "STICK";
            Vec3 color = fric.m_wasLimited ? Green : Orange;
            felt.showFrictionForce(state, geometry, color);
            geometry.push_back(DecorativeText(text)
                .setColor(color).setScale(.1)
                .setTransform(loc+Vec3(0.1,.04,0)));
        }
    }
private:
    const PGSTimeStepper& m_ts;
};

//==============================================================================
//                               TIM'S BOX
//==============================================================================
class TimsBox : public PGSAugmentedMultibodySystem {
public:
    TimsBox();

    void calcInitialState(State& state) const OVERRIDE_11;

    const MobilizedBody& getBodyToWatch() const OVERRIDE_11
    {   return m_brick; }

private:
    Force::Gravity          m_gravity;
    Force::GlobalDamper     m_damper;
    MobilizedBody::Free     m_brick;
    MobilizedBody::Ball     m_brick2;
};


//==============================================================================
//                              BOUNCING BALLS
//==============================================================================
class BouncingBalls : public PGSAugmentedMultibodySystem {
public:
    BouncingBalls();
    ~BouncingBalls() {delete m_contactForces; delete m_tracker;}

    void calcInitialState(State& state) const OVERRIDE_11;

    const MobilizedBody& getBodyToWatch() const OVERRIDE_11
    {   static const MobilizedBody nobod; return nobod; }

    const MobilizedBody::Slider& getHBall(int i) const {return m_Hballs[i];}
    const MobilizedBody::Slider& getPBall(int i) const {return m_Pballs[i];}
    const MobilizedBody::Slider& getNBall(int i) const {return m_Nballs[i];}

private:
    // Add subsystems for compliant contact. TODO: shouldn't need pointers
    ContactTrackerSubsystem*     m_tracker;
    CompliantContactSubsystem*   m_contactForces;

    static const int NBalls = 4;

    Force::Gravity          m_gravity;
    Force::GlobalDamper     m_damper;
    MobilizedBody::Slider   m_Hballs[NBalls];    // Hertz
    MobilizedBody::Slider   m_Pballs[NBalls];    // Poisson
    MobilizedBody::Slider   m_Nballs[NBalls];    // Newton
};
}

//==============================================================================
//                                   MAIN
//==============================================================================
int main(int argc, char** argv) {
    #ifdef USE_TIMS_PARAMS
        const Real RunTime=16;  // Tim's time
        const Real Accuracy = 1e-4;
    #else
        const Real RunTime=20;
        const Real Accuracy = 1e-2;
    #endif

    const bool UseNewton = false; // default is Poisson restitution

  try { // If anything goes wrong, an exception will be thrown.

    // Create the augmented multibody model.
    TimsBox mbs;
    //BouncingBalls mbs;
    PGSTimeStepper pgs(mbs);

    const SimbodyMatterSubsystem&    matter = mbs.getMatterSubsystem();
    const MyUnilateralConstraintSet& unis   = mbs.getUnis();

    Visualizer viz(mbs);
    viz.setShowSimTime(true);
    viz.setShowFrameNumber(true);
    viz.setShowFrameRate(true);
    viz.addDecorationGenerator(new ShowContact(pgs));

    if (!mbs.getBodyToWatch().isEmptyHandle())
        viz.addFrameController(new BodyWatcher(mbs.getBodyToWatch()));

    // Simulate it.
    State s;
    mbs.calcInitialState(s);

    printf("Initial state shown. ENTER to continue.\n");
    viz.report(s);
    getchar();

    const Real ConsTol = .001;
    const Real PGSConvergenceTol = 1e-5;
    const int  PGSMaxIters = 100;
    const Real PGSSor = 1.05; // successive over relaxation, 0..2, 1 is neutral

    pgs.setUseNewtonRestitution(UseNewton);
    pgs.setAccuracy(Accuracy); // integration accuracy
    pgs.setConstraintTol(ConsTol);

    pgs.setPGSConvergenceTol(PGSConvergenceTol);
    pgs.setPGSMaxIters(PGSMaxIters);
    pgs.setPGSSOR(PGSSor);

    pgs.initialize(s);
    mbs.resetAllCountersToZero();
    mbs.updUnis().initialize(); // reinitialize
        
    Array_<State> states; states.reserve(10000);

    int nSteps=0, nStepsWithEvent = 0;

    const double startReal = realTime();
    const double startCPU = cpuTime();

    const Real h = .001;
    const int SaveEvery = 33; // save every nth step ~= 33ms

    do {
        const State& pgsState = pgs.getState();
        if ((nSteps%SaveEvery)==0) {
            #ifdef ANIMATE
            viz.report(pgsState);
            #endif
            states.push_back(pgsState);
        }

        pgs.stepTo(pgsState.getTime() + h);

        ++nSteps;
    } while (pgs.getTime() < RunTime);


    const double timeInSec = realTime()-startReal;
    const double cpuInSec = cpuTime()-startCPU;
    cout << "Done -- took " << nSteps << " steps in " <<
        timeInSec << "s for " << pgs.getTime() << "s sim (avg step=" 
        << (1000*pgs.getTime())/nSteps << "ms) ";
    cout << "CPUtime " << cpuInSec << endl;

    printf("Used PGS solver (%s) at acc=%g consTol=%g"
           " convergenceTol=%g maxIters=%d SOR=%g.\n", 
           pgs.getUseNewtonRestitution() ? "Newton" : "Poisson",
           pgs.getAccuracy(), pgs.getConstraintTol(),
           pgs.getPGSConvergenceTol(), pgs.getPGSMaxIters(),
           pgs.getPGSSOR());
    printf("Compression: ncalls=%d, niters=%d (%g/call), nfail=%d\n",
           pgs.getPGSNumCalls(0), pgs.getPGSNumIters(0),
           (double)pgs.getPGSNumIters(0)/std::max(pgs.getPGSNumCalls(0),1LL),
           pgs.getPGSNumFailures(0));
    printf("Expansion: ncalls=%d, niters=%d (%g/call), nfail=%d\n",
           pgs.getPGSNumCalls(1), pgs.getPGSNumIters(1),
           (double)pgs.getPGSNumIters(1)/std::max(pgs.getPGSNumCalls(1),1LL),
           pgs.getPGSNumFailures(1));
    printf("Position: ncalls=%d, niters=%d (%g/call), nfail=%d\n",
           pgs.getPGSNumCalls(2), pgs.getPGSNumIters(2),
           (double)pgs.getPGSNumIters(2)/std::max(pgs.getPGSNumCalls(2),1LL),
           pgs.getPGSNumFailures(2));

    cout << "nstates =" << states.size() << endl;

    // Instant replay.
    while(true) {
        printf("Hit ENTER for replay (%d states) ...", 
                states.size());
        getchar();
        for (unsigned i=0; i < states.size(); ++i) {
            mbs.realize(states[i], Stage::Velocity);
            viz.report(states[i]);
        }
    }

  } 
  catch (const std::exception& e) {
    printf("EXCEPTION THROWN: %s\n", e.what());
    exit(1);
  }
  catch (...) {
    printf("UNKNOWN EXCEPTION THROWN\n");
    exit(1);
  }
}

//==============================================================================
//                           PGS TIME STEPPER
//==============================================================================
Integrator::SuccessfulStepStatus PGSTimeStepper::
stepTo(Real time) {
    // Abbreviations.
    const PGSAugmentedMultibodySystem&  mbs    = m_ambs;
    const SimbodyMatterSubsystem&       matter = mbs.getMatterSubsystem();
    const MyUnilateralConstraintSet&    unis   = mbs.getUnis();
    State&                              s      = m_state;

    const Real h = time - m_state.getTime();    // max timestep

    // Kinematics should already be realized so this won't do anything.
    mbs.realize(s, Stage::Position); 
    // Determine which constraints will be involved for this step.
    findProximalConstraints(s);

    // Enable all proximal constraints, reassigning multipliers if needed.
    enableProximalConstraints(s);
    mbs.realize(s, Stage::Velocity);

    const Vector& verr0 = s.getUErr();
    if (!verr0.size()) {
        takeUnconstrainedStep(s, h);
        return Integrator::SuccessfulStepStatus::ReachedScheduledEvent;
    }

    // Calculate velocity-dependent coefficients of restitution and friction
    // and apply combining rules for dissimilar materials.
    calcContactCoefficients(s, verr0);

    // We're going to accumulate velocity errors here, starting with the 
    // presenting violation.
    Vector verr = verr0;

    // If we're in Newton mode, or if a contact specifies Newton restitution,
    // then modify the appropriate verr's here.
    const bool anyNewton = applyNewtonRestitutionIfAny(s, verr);

    // Evaluate applied forces and get reference to them. These include
    // gravity but not centrifugal forces.
    mbs.realize(s, Stage::Dynamics);
    const Vector&              f = mbs.getMobilityForces(s, Stage::Dynamics);
    const Vector_<SpatialVec>& F = mbs.getRigidBodyForces(s, Stage::Dynamics);

    // Calculate the constraint compliance matrix GM\~G.
    matter.calcProjectedMInv(s, m_GMInvGt);

    // Calculate udotExt = M\(f + ~J*(F-C)) where C are centrifugal forces.
    // This is the unconstrained acceleration.
    Vector udotExt; Vector_<SpatialVec> A_GB;
    matter.calcAccelerationIgnoringConstraints(s,f,F,udotExt,A_GB);

    // There are some constraints to deal with.
    // Calculate verrExt = G*deltaU; the end-of-step constraint error due to 
    // external and centrifugal forces.
    Vector verrExt;
    matter.multiplyByG(s, h*udotExt, verrExt);

    // Update velocity error to include verrExt.
    verr += verrExt; // verr0+Newton restitution+external forces

    Vector impulse;
    //--------------------------------------------------------------------------
    // Initialize the impulse to a best guess from the previous step,
    // compute the compression impulse that eliminates verr.
    doCompressionPhase(s, verr, impulse);
    //--------------------------------------------------------------------------

    // Modify the compression impulse to add in expansion impulses.
    const bool anyPoisson = applyPoissonRestitutionIfAny(s, verr0, impulse);

    if (anyPoisson) {
        // Must calculate a reaction impulse. Move now-known impulse to RHS 
        // (convert to verr).
        Vector genImpulse, deltaU, expVerr, reacImpulse; 
        // constraint impulse -> gen impulse
        matter.multiplyByGTranspose(s, impulse, genImpulse);
        // gen impulse to delta-u, then to verr
        matter.multiplyByMInv(s, genImpulse, deltaU);
        matter.multiplyByG(s, deltaU, expVerr);

        // Update verr to include compression+expansion impulse.
        verr -= expVerr; // watch sign -- multipliers are opposite forces

        //----------------------------------------------------------------------
        doExpansionPhase(s, verr, reacImpulse); // initial guess=zero
        //----------------------------------------------------------------------

        impulse += reacImpulse;
    }

    // We need forces, not impulses for the next calculation.
    s.updMultipliers() = impulse/h;
    // Calculate constraint forces ~G*lambda (body frcs Fc, mobility frcs fc).
    Vector_<SpatialVec> Fc; Vector fc; 
    matter.calcConstraintForcesFromMultipliers
        (s,s.updMultipliers(), Fc, fc);

    // Now calculate udot = M\(f-fc + ~J*(F-Fc-C)).
    Vector udot;
    matter.calcAccelerationIgnoringConstraints(s,f-fc,F-Fc,udot,A_GB);

    // Update auxiliary states z, invalidating Stage::Dynamics.
    s.updZ() += h*s.getZDot();

    // Update u from deltaU, invalidating Stage::Velocity. 
    s.updU() += h*udot;

    // Done with velocity update. Now calculate qdot, possiblity including
    // an additional position error correction term.
    Vector qdot;
    const Vector& perr0 = s.getQErr();
    if (!anyPositionErrorsViolated(s, perr0)) {
        matter.multiplyByN(s,false,s.getU(),qdot);
    } else {
        // Don't include quaternions for position correction. 
        Vector posImpulse, genImpulse, posVerr, deltaU;
        posVerr.resize(s.getNUErr()); posVerr.setToZero();
        const int nQuat = matter.getNumQuaternionsInUse(s);
        posVerr(0, perr0.size()-nQuat) = perr0(0, perr0.size()-nQuat)/h;

        //----------------------------------------------------------------------
        // Calculate impulse and then deltaU=M\~G*impulse such that -h*deltaU 
        // will eliminate position errors, respecting only position constraints.
        doPositionCorrectionPhase(s, posVerr, posImpulse);
        //----------------------------------------------------------------------
        // constraint impulse -> gen impulse
        matter.multiplyByGTranspose(s, posImpulse, genImpulse);
        // gen impulse to deltaU (watch sign)
        matter.multiplyByMInv(s, genImpulse, deltaU);

        // convert corrected u to qdot (note we're not changing u)
        matter.multiplyByN(s,false,s.getU()-deltaU, qdot);
    }

    // We have qdot, now update q, fix quaternions, update time.
    s.updQ() += h*qdot; // invalidates Stage::Position
    matter.normalizeQuaternions(s);
    s.updTime() += h;   // invalidates Stage::Time

    // Return from step with kinematics realized. Note that we may have
    // broken the velocity constraints by updating q, but we won't fix that
    // until the next step. Also position constraints are only imperfectly
    // satisfied by the correction above.
    mbs.realize(s, Stage::Velocity);

    return Integrator::SuccessfulStepStatus::ReachedScheduledEvent;
}

//==============================================================================
//                           PROJECTED GAUSS SEIDEL
//==============================================================================
/* 
We are given 
    - A, square matrix of dimension mA 
    - b, rhs vector (length mA)
    - w, solution vector with initial value w=w0 (length mA)
representing mA scalar constraint equations A[i]*w=b[i].

A smaller square subset may be selected via
    - I, selection index set, a subset of IA={1,...,mA}

The selected subset I is partitioned into four disjoint index sets
    - IU: Unconditional
    - IB: Bounded scalar
    - IV: Length-limited vector
    - IF: Friction

Each bounded scalar constraint k provides
    - a single constraint index iB_k from IB, and 
    - lower and upper bounds lb_k, ub_k.

Each length-limited vector constraint k specifies 
    - a unique index set of 1-3 distinct constraints IV_k from IV,
    - a nonnegative scalar L_k specifing the maximum length of the vector.

Each friction constraint k specifies 
    - a unique index set of 1-3 distinct friction constraints IF_k from IF, 
    - an index set of 1-3 distinct normal constraints IN_k from IA-IF,
    - the effective coefficient of friction mu.
Note that the normal constraints in IN_k do not have to be selected from the 
active subset I; if not they will be fixed at w0[IN_k] on entry. 

Given those inputs, we attempt to solve: 
    A[I,I] w[I] = b[I]
    subject to lb_k <= w[iB_k] <= ub_k       for bounded constraints k in IB
    and        ||w[IV_k]|| <= L_k            for vector constraints k in IV
    and        ||w[IF_k]|| <= mu*||w[IN_k]|| for friction constraints k in IF

Implicitly, complementarity conditions must hold:
    w_i in interior of constraint -> A[i]*w == b[i]
    w_i on boundary of constraint -> A[i]*w != b[i]

*/

namespace { // local functions
/** Given a scalar s, ensure that lb <= s <= ub by moving s to the nearest
bound if necessary. Return true if any change is made. **/
bool boundScalar(Real lb, Real& s, Real ub) {
    assert(lb <= ub);
    if      (s > ub) {s=ub; return true;}
    else if (s < lb) {s=lb; return true;}
    return false;
}

/** Given an index set IV, ensure that ||w[IV]|| <= maxLen by scaling the
vector to that length if necessary. Return true if any change is made. **/
bool boundVector(Real maxLen, const Array_<int>& IV, Vector& w) {
    assert(maxLen >= 0);
    const Real maxLen2 = square(maxLen);
    Real wNorm2 = 0;
    for (unsigned i=0; i<IV.size(); ++i) wNorm2 += square(w[IV[i]]);
    if (wNorm2 <= maxLen2) 
        return false;
    const Real scale = std::sqrt(maxLen2/wNorm2); // 0 <= scale < 1
    for (unsigned i=0; i<IV.size(); ++i) w[IV[i]] *= scale;
    return true;
}

/** Given index set IN identifying the components of the normal force vector,
and index set IF identifying the components of the friction vector, ensure
that ||w[IF]|| <= mu*||w[IN]|| by scaling the friction vector if necessary.
Return true if any change is made. **/
bool boundFriction(Real mu, const Array_<int>& IN, 
                   const Array_<int>& IF, Vector& w) {
    assert(mu >= 0);
    Real N2=0, F2=0; // squares of normal and friction force magnitudes
    for (unsigned i=0; i<IN.size(); ++i) N2 += square(w[IN[i]]);
    for (unsigned i=0; i<IF.size(); ++i) F2 += square(w[IF[i]]);
    const Real mu2N2 = mu*mu*N2;
    if (F2 <= mu2N2) 
        return false;
    const Real scale = std::sqrt(mu2N2/F2); // 0 <= scale < 1
    for (unsigned i=0; i<IF.size(); ++i) w[IF[i]] *= scale;
    return true;
}
}

bool PGSTimeStepper::projGaussSeidel
                    (int phase,
                     const Matrix& A, const Vector& b, Vector& w,
                     const Array_<int>&     all, 
                     const Array_<int>&     unconditional, 
                     Array_<Bounded>&       bounded,
                     Array_<LengthLimited>& lengthLimited,
                     Array_<Frictional>&    frictional) const
{
    ++m_PGSNumCalls[phase];
    const int mA=A.nrow(), nA=A.ncol();
    assert(mA==nA); assert(b.nrow()==mA); assert(w.nrow()==nA);

    const int m = (int)all.size();
    assert(m<=mA);

    // Partitions of selected subset.
    const int mUncond  = (int)unconditional.size();
    const int mBounded = (int)bounded.size();
    const int mLength  = (int)lengthLimited.size();
    const int mFric    = (int)frictional.size();

    // If debugging, check for consistent constraint equation count.
    #ifndef NDEBUG
    {int mCount = mUncond + mBounded; // 1 each
    for (int k=0; k<mLength; ++k)
        mCount += lengthLimited[k].m_components.size();
    for (int k=0; k<mFric; ++k)
        mCount += frictional[k].m_Fk.size();
    assert(mCount == m);}
    #endif

    if (m == 0) {
        printf("PGS %d: nothing to do; converged in 0 iters.\n", phase);
        return true;
    }

    // Track total error for all included equations, and the error for just
    // those equations that are being enforced.
    bool converged = false;
    Real normRMSall = Infinity, normRMSenf = Infinity;
    for (int its=0; its < m_PGSMaxIters; ++its) {
        ++m_PGSNumIters[phase];
        const Real prevNormRMSall = normRMSall;
        Real sum2all = 0, sum2enf = 0; // track solution errors

        // UNCONDITIONAL: these are always on.
        for (int fx=0; fx < mUncond; ++fx) {
            const int rx = unconditional[fx];
            Real rowSum = 0;
            for (int c=0; c < m; ++c) {
                const int cx = all[c];
                rowSum += A(rx,cx)*w[cx];
            }
            const Real er = b[rx]-rowSum, er2=er*er;
            w[rx] += m_PGSSOR * er/A(rx,rx);
            sum2all += er2; sum2enf += er2;
        }

        // BOUNDED: conditional scalar constraints with constant bounds
        // on resulting w.
        for (int k=0; k < mBounded; ++k) {
            Bounded& bnd = bounded[k];
            const int rx = bnd.m_ix;
            Real rowSum = 0;
            for (int c=0; c < m; ++c) {
                const int cx = all[c];
                rowSum += A(rx,cx)*w[cx];
            }
            const Real er = b[rx]-rowSum, er2=er*er;
            w[rx] += m_PGSSOR * er/A(rx,rx);
            sum2all += er2;
            if (!(bnd.m_hitBound=boundScalar(bnd.m_lb, w[rx], bnd.m_ub)))
                sum2enf += er2;
        }

        // LENGTH: a set of constraint equations forming a vector whose
        // maximum length is limited.
        for (int k=0; k < mLength; ++k) {
            LengthLimited& len = lengthLimited[k];
            const Array_<int>& rows = len.m_components;
            Vec3 rowSums(0);
            for (int c=0; c < m; ++c) {
                const int cx = all[c];
                for (unsigned i=0; i<rows.size(); ++i)
                    rowSums[i] += A(rows[i],cx)*w[cx];
            }
            Real localEr2 = 0;
            for (unsigned i=0; i<rows.size(); ++i) {
                const int rx = rows[i];
                const Real er = b[rx]-rowSums[i];
                w[rx] += m_PGSSOR * er/A(rx,rx);
                localEr2 += square(er);
            }
            sum2all += localEr2;
            if (!(len.m_hitLimit=boundVector(len.m_maxLength, rows, w)))
                sum2enf += localEr2;
        }

        // FRICTIONAL: a set of constraint equations forming a vector whose
        // maximum length is limited by the norm of other multipliers w.
        for (int k=0; k < mFric; ++k) {
            Frictional& fric = frictional[k];
            const Array_<int>& Fk = fric.m_Fk; // friction components
            Vec3 rowSums(0);
            for (int c=0; c < m; ++c) {
                const int cx = all[c];
                for (unsigned i=0; i<Fk.size(); ++i)
                    rowSums[i] += A(Fk[i],cx)*w[cx];
            }
            Real localEr2 = 0;
            for (unsigned i=0; i<Fk.size(); ++i) {
                const int rx = Fk[i];
                const Real er = b[rx]-rowSums[i];
                w[rx] += m_PGSSOR * er/A(rx,rx);
                localEr2 += square(er);
            }
            sum2all += localEr2;
            if (!(fric.m_wasLimited=boundFriction(fric.m_effMu,fric.m_Nk,Fk,w)))
                sum2enf += localEr2;
        }
        normRMSall = std::sqrt(sum2all/m);
        normRMSenf = std::sqrt(sum2enf/m);

        SimTK_DEBUG4("%d/%d: EST rmsAll=%g rmsEnf=%g\n", phase, its,
                     normRMSall, normRMSenf);
        #ifdef NDEBUG // i.e., NOT debugging (TODO)
        if (its > 50)
            printf("%d/%d: EST rmsAll=%g rmsEnf=%g\n", phase, its,
                     normRMSall, normRMSenf);
        #endif
        //if (   normRMSall < m_PGSConvergenceTol 
        //    || (  std::abs(normRMSall-prevNormRMSall) 
        //        < m_PGSConvergenceTol*prevNormRMSall)) 
        if (normRMSenf < m_PGSConvergenceTol) //TODO: add failure-to-improve check
        {
            SimTK_DEBUG3("PGS %d converged to %g in %d iters\n", 
                         phase, normRMSenf, its);
            converged = true;
            break;
        }
    }

    if (!converged) {
        printf("PGS %d CONVERGENCE FAILURE: %d iters -> norm=%g\n",
               phase, m_PGSMaxIters, normRMSenf);
        ++m_PGSNumFailures[phase];
    }

    return converged;
}

void PGSTimeStepper::initialize(const State& initState) {
    m_state = initState;
    m_ambs.realize(m_state, Stage::Acceleration);
}

// Determine which constraints will be involved for this step.
void PGSTimeStepper::
findProximalConstraints(const State& s) { //TODO: redo
    const MyUnilateralConstraintSet& unis = m_ambs.getUnis();
    unis.findProximalElements(s, m_consTol, m_proximals, m_distals);
}

// Enable all proximal constraints, disable all distal constraints, 
// reassigning multipliers if needed. Returns true if any change was made.
bool PGSTimeStepper::
enableProximalConstraints(State& s) {
    const MyUnilateralConstraintSet& unis = m_ambs.getUnis();

    // Record friction application points. This has to be done while Position 
    // stage is still valid.
    for (unsigned i=0; i < m_proximals.m_friction.size(); ++i) {
        const int id = m_proximals.m_friction[i];
        unis.updFrictionElement(id).initializeFriction(s);
    }

    bool changed = false;

    // Disable non-proximal constraints if they were previously disabled.
    for (unsigned i=0; i < m_distals.m_friction.size(); ++i) {
        const int id = m_distals.m_friction[i];
        const MyFrictionElement& fric = unis.getFrictionElement(id);
        if (fric.isEnabled(s)) fric.disable(s), changed=true;
    }
    for (unsigned i=0; i < m_distals.m_contact.size(); ++i) {
        const int id = m_distals.m_contact[i];
        const MyContactElement& cont = unis.getContactElement(id);
        if (!cont.isDisabled(s)) cont.disable(s), changed=true;
    }

    for (unsigned i=0; i < m_proximals.m_contact.size(); ++i) {
        const int id = m_proximals.m_contact[i];
        const MyContactElement& cont = unis.getContactElement(id);
        if (cont.isDisabled(s)) cont.enable(s), changed=true;
    }
    for (unsigned i=0; i < m_proximals.m_friction.size(); ++i) {
        const int id = m_proximals.m_friction[i];
        const MyFrictionElement& fric = unis.getFrictionElement(id);
        fric.setInstanceParameters(s);
        if (!fric.isEnabled(s)) fric.enable(s), changed=true;
    }

    // TODO: Note that we always have to move the friction application points
    // which is an Instance stage change; shouldn't be.
    m_ambs.realize(s, Stage::Instance); // assign multipliers

    return changed;
}

void PGSTimeStepper::takeUnconstrainedStep(State& s, Real h) {
    const PGSAugmentedMultibodySystem&  mbs    = m_ambs;
    const SimbodyMatterSubsystem&       matter = mbs.getMatterSubsystem();
    mbs.realize(s, Stage::Acceleration);
    const Vector& udot = s.getUDot(); // grab before invalidated
    s.updZ() += h*s.getZDot(); // invalidates Stage::Dynamics
    s.updU() += h*udot;        // invalidates Stage::Velocity
    Vector qdot;
    matter.multiplyByN(s,false,s.getU(),qdot);
    s.updQ() += h*qdot;         // invalidates Stage::Position
    matter.normalizeQuaternions(s);
    s.updTime() += h;           // invalidates Stage::Time
    mbs.realize(s, Stage::Velocity);
}


// Calculate velocity-dependent coefficients of restitution and friction
// and apply combining rules for dissimilar materials.
void PGSTimeStepper::
calcContactCoefficients(const State& s, const Vector& verr) {
    const MyUnilateralConstraintSet& unis = m_ambs.getUnis();
    m_bounded.clear();
    for (unsigned i=0; i < m_proximals.m_contact.size(); ++i) {
        const int id = m_proximals.m_contact[i];
        const MyContactElement& elt = unis.getContactElement(id);
        const MultiplierIndex mx = elt.getMultIndex(s);
        Real cor = elt.getMaxCoefRest();
        if (verr[mx] >= -unis.getCaptureVelocity()) cor=0;
        SimTK_DEBUG3("Contact %d speed=%g -> N cor=%g\n", 
                     i, verr[mx], cor);
        m_bounded.push_back(
            Bounded(mx,-Infinity, Zero, cor)); // TODO: fix sign convention
    }

    m_frictional.clear();
    for (unsigned i=0; i < m_proximals.m_friction.size(); ++i) {
        const int id = m_proximals.m_friction[i];
        const MyFrictionElement& felt = unis.getFrictionElement(id);
        const Real v_trans = unis.getTransitionVelocity();
        const Real mu_d = felt.getDynamicFrictionCoef();
        const Real mu_s = felt.getStaticFrictionCoef();
        const Real mu_v = felt.getViscousFrictionCoef();
        const MyPointContactFriction& pelt = // TODO: generalize
            dynamic_cast<const MyPointContactFriction&>(felt);
        const MultiplierIndex mx = pelt.getMultIndexX(s);
        const MultiplierIndex my = pelt.getMultIndexY(s);
        const MultiplierIndex mN = pelt.getMyPointContact().getMultIndex(s);
        const Vec2 vSlip(verr[mx], verr[my]);
        const Real speed = vSlip.norm();
        const Real mu = (speed<=v_trans? mu_s : mu_d + speed*mu_v);
        SimTK_DEBUG3("Fric %d speed=%g -> mu=%g\n", i, speed, mu);
        Array_<MultiplierIndex> Fk; Fk.push_back(mx); Fk.push_back(my);
        Array_<MultiplierIndex> Nk; Nk.push_back(mN);
        m_frictional.push_back(Frictional(Fk,Nk,mu));
    }

    const int m = s.getNUErr();
    m_all.clear(); m_uncond.clear(); 
    for (int i=0; i<m; ++i) 
        m_all.push_back(MultiplierIndex(i));
    // TODO: add in unconditionals in m_uncond

    m_allPos.clear(); m_boundedPos.clear(); m_uncondPos.clear();
    for (unsigned i=0; i<m_bounded.size(); ++i) {
        m_allPos.push_back(m_bounded[i].m_ix);
        m_boundedPos.push_back(m_bounded[i]); // should be holonomics only
    }
    // TODO: add in unconditionals in m_uncondPos
}

bool PGSTimeStepper::
applyNewtonRestitutionIfAny(const State& s, Vector& verr) const {
    if (!m_useNewton) 
        return false; //TODO: check individual contacts
    bool anyRestitution = false;
    const PGSAugmentedMultibodySystem&  mbs    = m_ambs;
    const MyUnilateralConstraintSet&    unis   = mbs.getUnis();
    for (unsigned i=0; i < m_proximals.m_contact.size(); ++i) {
        const int id = m_proximals.m_contact[i];
        const MyContactElement& elt = unis.getContactElement(id);
        const MultiplierIndex mx = elt.getMultIndex(s);
        const Bounded& bounded = m_bounded[i];
        if (bounded.m_effCOR != 0) {
            verr[mx] *= (1+bounded.m_effCOR);
            anyRestitution = true;
        }
    }
    return anyRestitution;
}

bool PGSTimeStepper::
applyPoissonRestitutionIfAny(const State& s, const Vector& verr,
                             Vector& impulse) const {
    if (m_useNewton) 
        return false; //TODO: check individual contacts
    bool anyRestitution = false;
    const PGSAugmentedMultibodySystem&  mbs    = m_ambs;
    const MyUnilateralConstraintSet&    unis   = mbs.getUnis();
    for (unsigned i=0; i < m_proximals.m_contact.size(); ++i) {
        const int id = m_proximals.m_contact[i];
        const MyContactElement& elt = unis.getContactElement(id);
        const MultiplierIndex mx = elt.getMultIndex(s);
        const Bounded& bounded = m_bounded[i];
        if (bounded.m_effCOR != 0) {
            impulse[mx] *= (1+bounded.m_effCOR);
            anyRestitution = true;
        }
    }
    return anyRestitution;
}

// This phase uses all the proximal constraints and should use a starting
// guess for impulse saved from the last step if possible.
bool PGSTimeStepper::
doCompressionPhase(const State&, const Vector& eps, Vector& compImpulse) {
    // TODO: improve initial guess
    compImpulse.resize(m_GMInvGt.ncol()); compImpulse.setToZero();
    bool converged = projGaussSeidel(0, m_GMInvGt, eps, compImpulse, 
                                     m_all, m_uncond, m_bounded, 
                                     Array_<LengthLimited>(), m_frictional);
    return converged;
}
// This phase uses all the proximal constraints, but we expect the result
// to be zero unless expansion causes new violations.
bool PGSTimeStepper::
doExpansionPhase(const State&, const Vector& eps, Vector& reactionImpulse) {
    // TODO: improve initial guess
    reactionImpulse.resize(m_GMInvGt.ncol()); reactionImpulse.setToZero();
    bool converged = projGaussSeidel(1, m_GMInvGt, eps, reactionImpulse, 
                                     m_all, m_uncond, m_bounded, 
                                     Array_<LengthLimited>(), m_frictional);
    return converged;
}

// This phase uses only holonomic constraints, and zero is a good initial
// guess for the (hopefully small) position correction.
bool PGSTimeStepper::
doPositionCorrectionPhase(const State&, const Vector& eps,
                          Vector& positionImpulse) {
    positionImpulse.resize(m_GMInvGt.ncol()); positionImpulse.setToZero();
    bool converged = projGaussSeidel(2, m_GMInvGt, eps, positionImpulse, 
                                m_allPos, m_uncondPos, m_boundedPos, 
                                Array_<LengthLimited>(), Array_<Frictional>());
    return converged;
}

bool PGSTimeStepper::
anyPositionErrorsViolated(const State&, const Vector& perr) const {
    // TODO: no need to fix if large perrs satisfy inequalities.
    bool anyViolated = perr.normInf() > m_consTol;
    SimTK_DEBUG2("maxAbs(perr)=%g -> %s\n", perr.normInf(),
                anyViolated ? "VIOLATED" : "OK");
    return anyViolated;
}

//==============================================================================
//                               TIM'S BOX
//==============================================================================
TimsBox::TimsBox() {
    // Abbreviations.
    SimbodyMatterSubsystem&     matter = updMatterSubsystem();
    GeneralForceSubsystem&      forces = updForceSubsystem();
    MyUnilateralConstraintSet&  unis   = updUnis();
    MobilizedBody&              Ground = matter.updGround();

    // Build the multibody system.
    m_gravity = Force::Gravity(forces, matter, -YAxis, 9.8066);
    m_damper  = Force::GlobalDamper(forces, matter, .1);

    // Predefine some handy rotations.
    const Rotation Z90(Pi/2, ZAxis); // rotate +90 deg about z

    const Vec3 BrickHalfDims(.1, .25, .5);
    const Real BrickMass = /*10*/5;
    #ifdef USE_TIMS_PARAMS
        const Real RunTime=16;  // Tim's time
        const Real Stiffness = 2e7;
        const Real Dissipation = 1;
        const Real CoefRest = 0; 
        // Painleve problem with these friction coefficients.
        //const Real mu_d = 1; /* compliant: .7*/
        //const Real mu_s = 1; /* compliant: .7*/
        const Real mu_d = .5;
        const Real mu_s = .8;
        const Real mu_v = /*0.05*/0;
        const Real CaptureVelocity = 0.01;
        const Real TransitionVelocity = 0.01;
        const Inertia brickInertia(.1,.1,.1);
        const Real Radius = .02;
    #else
        const Real RunTime=20;
        const Real Stiffness = 1e6;
        const Real CoefRest = 0.3; 
        const Real TargetVelocity = 3; // speed at which to match coef rest
//        const Real Dissipation = (1-CoefRest)/TargetVelocity;
        const Real Dissipation = 0.1;
        const Real mu_d = .5;
        const Real mu_s = 1.0;
        const Real mu_v = 0*0.1;
        const Real CaptureVelocity = 0.01;
        const Real TransitionVelocity = 0.05;
        const Inertia brickInertia(BrickMass*UnitInertia::brick(BrickHalfDims));
        const Real Radius = BrickHalfDims[0]/3;
    #endif

    unis.setCaptureVelocity(CaptureVelocity);
    unis.setTransitionVelocity(TransitionVelocity);

    printf("\n******************** Tim's Box ********************\n");
    printf("USING RIGID CONTACT\n");
    #ifdef USE_TIMS_PARAMS
    printf("Using Tim's parameters:\n");
    #else
    printf("Using Sherm's parameters:\n");
    #endif
    printf("  coef restitution=%g\n", CoefRest);
    printf("  mu_d=%g mu_s=%g mu_v=%g\n", mu_d, mu_s, mu_v);
    printf("  transition velocity=%g\n", TransitionVelocity);
    printf("  radius=%g\n", Radius);
    printf("  brick inertia=%g %g %g\n",
        brickInertia.getMoments()[0], brickInertia.getMoments()[1], 
        brickInertia.getMoments()[2]); 
    printf("******************** Tim's Box ********************\n\n");

        // ADD MOBILIZED BODIES AND CONTACT CONSTRAINTS

    Body::Rigid brickBody = 
        Body::Rigid(MassProperties(BrickMass, Vec3(0), brickInertia));
    brickBody.addDecoration(Transform(), DecorativeBrick(BrickHalfDims)
                                   .setColor(Red).setOpacity(.3));
    m_brick = MobilizedBody::Free(Ground, Vec3(0),
                                  brickBody, Vec3(0));
    m_brick2 = MobilizedBody::Ball(m_brick, BrickHalfDims,
                                   brickBody, Vec3(-BrickHalfDims));
    //m_brick3 = MobilizedBody::Ball(brick2, BrickHalfDims,
    //                          brickBody, Vec3(-BrickHalfDims));

/*
1) t= 0.5, dt = 2 sec, pt = (0.05, 0.2, 0.4), fdir = (1,0,0), mag = 50N
2) t= 4.0, dt = 0.5 sec, pt = (0.03, 0.06, 0.09), fdir = (0.2,0.8,0), mag = 300N
3) t= 0.9, dt = 2 sec, pt = (0,0,0), fdir = (0,1,0), mag = 49.333N (half the weight of the block)
4) t= 13.0, dt = 1 sec, pt = (0 0 0), fdir = (-1,0,0), mag = 200N
*/
    Force::Custom(forces, new MyPushForceImpl(m_brick, Vec3(0.05,0.2,0.4),
                                                    50 * Vec3(1,0,0),
                                                    0.5, 0.5+2));
    Force::Custom(forces, new MyPushForceImpl(m_brick, Vec3(0.03, 0.06, 0.09),
                                                    300 * UnitVec3(0.2,0.8,0),
                                                    //300 * Vec3(0.2,0.8,0),
                                                    4, 4+0.5));
    Force::Custom(forces, new MyPushForceImpl(m_brick, Vec3(0),
                                                    1.25*49.033 * Vec3(0,1,0),
                                                    9., 9.+2));
    Force::Custom(forces, new MyPushForceImpl(m_brick, Vec3(0),
                                                    200 * Vec3(-1,0,0),
                                                    13, 13+1));

    #ifndef USE_TIMS_PARAMS
    // Extra late force.
    Force::Custom(forces, new MyPushForceImpl(m_brick, Vec3(.1, 0, .45),
                                                    20 * Vec3(-1,-1,.5),
                                                    15, Infinity));
    #endif

    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        const Vec3 pt = Vec3(i,j,k).elementwiseMultiply(BrickHalfDims);
        MyPointContact* contact = new MyPointContact
           (Ground, YAxis, 0., m_brick, pt, CoefRest);
        unis.addContactElement(contact);
        unis.addFrictionElement(
            new MyPointContactFriction(*contact, mu_d, mu_s, mu_v, 
                                       CaptureVelocity, // TODO: vtol?
                                       forces));
        if (i==-1 && j==-1 && k==-1)
            continue;
        MyPointContact* contact2 = new MyPointContact
           (Ground, YAxis, 0., m_brick2, pt, CoefRest);
        unis.addContactElement(contact2);
        unis.addFrictionElement(
            new MyPointContactFriction(*contact2, mu_d, mu_s, mu_v, 
                                       CaptureVelocity, // TODO: vtol?
                                       forces));
        //MyPointContact* contact3 = new MyPointContact
        //  (Ground, YAxis, 0., m_brick3, pt, CoefRest);
        //unis.addContactElement(contact3);
        //unis.addFrictionElement(
        //    new MyPointContactFriction(*contact3, mu_d, mu_s, mu_v, 
        //                               CaptureVelocity, // TODO: vtol?
        //                               forces));
    }
}

//---------------------------- CALC INITIAL STATE ------------------------------
void TimsBox::calcInitialState(State& s) const {
    s = realizeTopology(); // returns a reference to the the default state
    
    //matter.setUseEulerAngles(s, true);
    
    realizeModel(s); // define appropriate states for this System
    realize(s, Stage::Instance); // instantiate constraints if any


    /*
    rX_q = 0.7 rad
    rX_u = 1.0 rad/sec

    rY_q = 0.6 rad
    rY_u = 0.0 rad/sec

    rZ_q = 0.5 rad
    rZ_u = 0.2 rad/sec

    tX_q = 0.0 m
    tX_u = 10 m/s

    tY_q = 1.4 m
    tY_u = 0.0 m/s

    tZ_q = 0.0 m
    tZ_u = 0.0 m/s
    */

    #ifdef USE_TIMS_PARAMS
        m_brick.setQToFitTranslation(s, Vec3(0,10,0));
        m_brick.setUToFitLinearVelocity(s, Vec3(0,0,0));
    #else
        m_brick.setQToFitTranslation(s, Vec3(0,1.4,0));
        m_brick.setUToFitLinearVelocity(s, Vec3(10,0,0));
        const Rotation R_BC(SimTK::BodyRotationSequence,
                                    0.7, XAxis, 0.6, YAxis, 0.5, ZAxis);
        m_brick.setQToFitRotation(s, R_BC);
        m_brick.setUToFitAngularVelocity(s, Vec3(1,0,.2));
    #endif

    realize(s, Stage::Position);
    Assembler(*this).setErrorTolerance(1e-6).assemble(s);
}

//==============================================================================
//                              BOUNCING BALLS
//==============================================================================

BouncingBalls::BouncingBalls() {
    m_tracker       = new ContactTrackerSubsystem(*this);
    m_contactForces = new CompliantContactSubsystem(*this, *m_tracker);

    // Abbreviations.
    SimbodyMatterSubsystem&     matter = updMatterSubsystem();
    GeneralForceSubsystem&      forces = updForceSubsystem();
    MyUnilateralConstraintSet&  unis   = updUnis();
    MobilizedBody&              Ground = matter.updGround();


    // Build the multibody system.
    m_gravity = Force::Gravity(forces, matter, -YAxis, 9.8066);
    //m_damper  = Force::GlobalDamper(forces, matter, .1);

    const Real BallMass = 1;
    const Real BallRadius = .25;
    const Real CoefRest = 1;
    const Real CaptureVelocity = .001;
    const Real TransitionVelocity = .001;

    // Rubber
    const Real rubber_density = 1100.;  // kg/m^3
    const Real rubber_young   = 0.01e9; // pascals (N/m)
    const Real rubber_poisson = 0.5;    // ratio
    const Real rubber_planestrain = 
        ContactMaterial::calcPlaneStrainStiffness(rubber_young,rubber_poisson);
    const Real rubber_dissipation = 0.1;
    const ContactMaterial rubber(rubber_planestrain,rubber_dissipation,0,0,0);
    // Nylon
    const Real nylon_density = 1100.;  // kg/m^3
    const Real nylon_young   = 2.5e9;  // pascals (N/m)
    const Real nylon_poisson = 0.4;    // ratio
    const Real nylon_planestrain =
        ContactMaterial::calcPlaneStrainStiffness(nylon_young,nylon_poisson);
    const Real nylon_dissipation = 0*0.1;
    const ContactMaterial nylon(nylon_planestrain,nylon_dissipation,0,0,0);

    const Rotation X2Y(Pi/2, ZAxis); // rotate +90 deg about z
    const Rotation NegX2Y(-Pi/2,ZAxis); // -90

    Ground.updBody().addContactSurface(Transform(NegX2Y,Vec3(0)),
                ContactSurface(ContactGeometry::HalfSpace(),nylon));

    unis.setCaptureVelocity(CaptureVelocity);
    unis.setTransitionVelocity(TransitionVelocity);



        // ADD MOBILIZED BODIES AND CONTACT CONSTRAINTS

    Body::Rigid ballBody(MassProperties(BallMass, Vec3(0), 
                                        UnitInertia::sphere(BallRadius)));
    ballBody.addDecoration(Transform(), DecorativeSphere(BallRadius));

    const Vec3 HColor(Gray), PColor(Red), NColor(Orange);

#ifdef HERTZ
    m_Hballs[0] = MobilizedBody::Slider
       (Ground, Transform(X2Y,Vec3(-1,BallRadius,0)),
        ballBody, X2Y);
    m_Hballs[0].updBody().addContactSurface(Vec3(0),
            ContactSurface(ContactGeometry::Sphere(BallRadius), nylon));
    m_Hballs[0].updBody().updDecoration(0).setColor(HColor);
    for (int i=1; i<NBalls; ++i) {
        m_Hballs[i] = MobilizedBody::Slider
           (m_Hballs[i-1],Transform(X2Y,Vec3(0,2*BallRadius,0)),ballBody, X2Y);
        m_Hballs[i].updBody().updDecoration(0).setColor(HColor);
        m_Hballs[i].updBody().addContactSurface(Vec3(0),
                ContactSurface(ContactGeometry::Sphere(BallRadius), nylon));
    }
#endif
#ifdef POISSON
    m_Pballs[0] = MobilizedBody::Slider
       (Ground, Transform(X2Y,Vec3(0,BallRadius,0)),
        ballBody, X2Y);
    m_Pballs[0].updBody().updDecoration(0).setColor(PColor);
    unis.addContactElement(new MyPointContact
           (Ground, YAxis, 0., m_Pballs[0], Vec3(0,-BallRadius,0), CoefRest));
    for (int i=1; i<NBalls; ++i) {
        m_Pballs[i] = MobilizedBody::Slider
           (m_Pballs[i-1],Transform(X2Y,Vec3(0,2*BallRadius,0)),ballBody, X2Y);
        m_Pballs[i].updBody().updDecoration(0).setColor(PColor);
        unis.addContactElement(new MyPointContact
               (m_Pballs[i-1], YAxis, BallRadius, 
                m_Pballs[i], Vec3(0,-BallRadius,0), CoefRest));
    }

#endif
#ifdef NEWTON
    m_Nballs[0] = MobilizedBody::Slider
       (Ground, Transform(X2Y,Vec3(1,BallRadius,0)),
        ballBody, X2Y);
    m_Nballs[0].updBody().updDecoration(0).setColor(NColor);
    unis.addContactElement(new MyPointContact
           (Ground, YAxis, 0., m_Nballs[0], Vec3(0,-BallRadius,0), CoefRest));
    for (int i=1; i<NBalls; ++i) {
        m_Nballs[i] = MobilizedBody::Slider
           (m_Nballs[i-1],Transform(X2Y,Vec3(0,2*BallRadius,0)),ballBody, X2Y);
        m_Nballs[i].updBody().updDecoration(0).setColor(NColor);
        unis.addContactElement(new MyPointContact
               (m_Nballs[i-1], YAxis, BallRadius, 
                m_Nballs[i], Vec3(0,-BallRadius,0), CoefRest));
    }
#endif

}

static const Real Separation = .1;
void BouncingBalls::calcInitialState(State& s) const {
    s = realizeTopology(); // returns a reference to the the default state   
    realizeModel(s); // define appropriate states for this System
    realize(s, Stage::Instance); // instantiate constraints if any
    realize(s, Stage::Position);
    Assembler(*this).setErrorTolerance(1e-6).assemble(s);
    #ifdef HERTZ
        getHBall(0).setOneQ(s, MobilizerQIndex(0), 1);
        getHBall(0).setOneU(s, MobilizerUIndex(0), -2);
        for (int i=1; i<NBalls; ++i) 
            getHBall(i).setOneQ(s, MobilizerQIndex(0), Separation);
    #endif
    #ifdef POISSON
        getPBall(0).setOneQ(s, MobilizerQIndex(0), 1);
        getPBall(0).setOneU(s, MobilizerUIndex(0), -2);
        for (int i=1; i<NBalls; ++i) 
            getPBall(i).setOneQ(s, MobilizerQIndex(0), Separation);
    #endif
    #ifdef NEWTON
        getNBall(0).setOneQ(s, MobilizerQIndex(0), 1);
        getNBall(0).setOneU(s, MobilizerUIndex(0), -2);
        for (int i=1; i<NBalls; ++i) 
            getNBall(i).setOneQ(s, MobilizerQIndex(0), Separation);
    #endif
}


//-------------------------- SHOW CONSTRAINT STATUS ----------------------------
void MyUnilateralConstraintSet::
showConstraintStatus(const State& s, const String& place) const
{
#ifndef NDEBUG
    printf("\n%s: uni status @t=%.15g\n", place.c_str(), s.getTime());
    m_mbs.realize(s, Stage::Dynamics);
    for (int i=0; i < getNumContactElements(); ++i) {
        const MyContactElement& contact = getContactElement(i);
        const bool isActive = !contact.isDisabled(s);
        printf("  %6s %2d %s, h=%g dh=%g f=%g\n", 
                isActive?"ACTIVE":"off", i, contact.getContactType().c_str(), 
                contact.getPerr(s),contact.getVerr(s),
                isActive?contact.getForce(s):Zero);
    }
    for (int i=0; i < getNumFrictionElements(); ++i) {
        const MyFrictionElement& friction = getFrictionElement(i);
        if (!friction.isMasterActive(s))
            continue;
        const bool isEnabled = friction.isEnabled(s);
        printf("  %8s friction %2d\n", 
                isEnabled?"STICKING":"sliding", i);
        friction.writeFrictionInfo(s, "    ", std::cout);
    }
    printf("\n");
#endif
}
