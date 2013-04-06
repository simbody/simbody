/* -------------------------------------------------------------------------- *
 *               Simbody(tm) - Tim's Box (hybrid contact model)               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

/* Investigate a hybrid compliant contact/rigid stiction model, using Tim's
box as a test case. */

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

// Set to revert to the no-constraint stiction model for performance comparison
// with the new constraint-based one.
//#define USE_CONTINUOUS_STICTION


// This is the continuous stiction model used in Simbody's compliant contact
// system for comparison with the new hybrid model. See end of this file for 
// implementation.
static Real stribeck(Real us, Real ud, Real uv, Real v);

//==============================================================================
//                     MY HYBRID VERTEX CONTACT ELEMENT
//==============================================================================
// Given a contact material (with compliant material properties), this 
// represents a compliant point-and-rigid-halfspace contact element. We track
// the height h of a vertex V on body B over a halfspace H on body P. The 
// halfspace frame H serves as the contact frame, with Hz in the contact 
// normal direction and Hxy spanning the tangent plane. The pose X_PH is 
// a given constant. 
//
// Normal force
// ------------
// If the vertex height h is above the halfspace surface (h>=0), no forces are 
// generated. If below the surface (h<0) we generate a normal force of magnitude 
//      N=max(0, -k*h(1-d*hdot)), with N >= 0
// applied to both bodies at the contact point C, which is the point at h=0 just
// up the halfspace normal direction from the vertex (because we're considering
// the halfspace to be rigid). Denote by Cb the station (material point) of B 
// that is coincident with C, and Cp the station of P that is coincident with C.
//
// Sliding force
// -------------
// Then if we are in sliding mode, we also generate a tangential force of
// magnitude T=(mu_d+mu_v*|v|)*N, where mu_d, mu_v are the dynamic and viscous
// coefficients of friction, and v is the velocity of Cb relative to Cp, 
// projected into the tangent plane (a 2d vector). If |v|>tol then the 
// tangential slip direction is s=v/|v|, and we record this as the previous slip
// direction s_prev. Otherwise the slip direction remains s=s_prev. In any case 
// the tangential force applied is -T*s.
// 
// Stiction force
// --------------
// If we are instead in stiction mode, then no sliding force is generated.
// Instead a pair of no-slip constraints is active, generating constraint
// multipliers (tx,ty), and we record s_prev=-[tx,ty]/|[tx,ty]| as the 
// impending slip direction. Note that this will not be available until
// Acceleration stage, while the normal force and tangential sliding force can
// be calculated at Velocity stage.
//
// Witness function
// ----------------
// (trigger on + -> -):
//     sliding_to_stiction = dot(v, s_prev)
//     stiction_to_sliding = mu_s*N - |[tx,ty]|

// This is a velocity-stage cache entry.
struct MyHybridContactInfo {
    MyHybridContactInfo()
    :   h(NaN), Cp(NaN), Cb(NaN), v_HCb(NaN), vSlipMag(NaN), f_HCb(NaN) 
    {
    }
    // Position info.
    Real        h;    // signed distance, h<0 means contact
    Vec3        Cp;   // station on P coincident with the contact point
    Vec3        Cb;   // station on B coincident with the contact point
    Rotation    R_GH;

    // Velocity info. H frame has z along contact normal, x,y in tangent plane.
    Vec3        v_HCb;      // velocity of Cb in the contact frame H
    Real        vSlipMag;   // |(v_HCb.x, v_HCb.y)|
    Vec3        f_HCb;      // contact force on Cb in contact frame H
};

class MyHybridVertexContactElementImpl : public Force::Custom::Implementation {
public:
    MyHybridVertexContactElementImpl(const GeneralForceSubsystem& forces,
        MobilizedBody& hsmobod, const UnitVec3& hsn, Real hsh,
        MobilizedBody& vmobod, const Vec3& vertex,
        const ContactMaterial& material)
    :   m_matter(hsmobod.getMatterSubsystem()), m_forces(forces), 
        m_hsmobodx(hsmobod), m_X_PH(Rotation(hsn, ZAxis), hsh*hsn), 
        m_vmobodx(vmobod), m_vertex_B(vertex),
        m_material(material),
        m_noslipX(hsmobod, Vec3(NaN), m_X_PH.x(), hsmobod, vmobod),
        m_noslipY(hsmobod, Vec3(NaN), m_X_PH.y(), hsmobod, vmobod), 
        m_index(-1), m_vtrans(NaN), // assign later
        m_contactPointInP(NaN), m_recordedSlipDir(NaN)
    {
        m_noslipX.setDisabledByDefault(true);
        m_noslipY.setDisabledByDefault(true);
    }

    void setContactIndex(int index) {m_index=index;}
    void setTransitionVelocity(Real vtrans) {m_vtrans=vtrans;}

    void initialize() { // TODO
        m_contactPointInP = NaN;
    }

    // Set the friction application point to be the projection of the contact 
    // point onto the contact plane. This will be used the next time stiction
    // is enabled. Requires state realized to Position stage.
    void initializeForStiction(const State& s) {
        const Real h = findContactPointInP(s, m_contactPointInP);
    }

    bool isInContact(const State& s) const
    {   return findH(s) < 0; }

    bool isSticking(const State& s) const
    {   return !m_noslipX.isDisabled(s); } // X,Y always on or off together
    
    // Return a point in Ground coincident with the vertex.
    Vec3 whereToDisplay(const State& s) const {
        return getBodyB().findStationLocationInGround(s, m_vertex_B);
    }

    Real getActualSlipSpeed(const State& s) const {
        const MyHybridContactInfo& info = getContactInfo(s);
        return info.vSlipMag;
    }

    // Return the normal force N >= 0 currently being generated by this 
    // contact. State must be realized to Stage::Velocity.
    Real getNormalForce(const State& s) const {
        const MyHybridContactInfo& info = getContactInfo(s);
        if (info.h >= 0) return 0;
        const Real N = info.f_HCb[2]; // z component is normal
        return N;
    }

    Real getHeight(const State& s) const {
        const MyHybridContactInfo& info = getContactInfo(s);
        return info.h;
    }

    Real getHeightDot(const State& s) const {
        const MyHybridContactInfo& info = getContactInfo(s);
        return info.v_HCb[2]; // velocity in plane normal direction
    }

    // Return the sliding force 2-vector in the halfspace plane that is being
    // applied to the contact point station on body B. If we are sticking
    // then this force is zero; call getStickingForce() to get the real value.
    const Vec2& getSlidingForce(const State& s) const {
        const MyHybridContactInfo& info = getContactInfo(s);
        const Vec2& f = info.f_HCb.getSubVec<2>(0);
        return f;
    }

    // Return the sliding velocity 2-vector in the halfspace plane of the
    // contact point station on body B relative to the contact point station
    // on body P.
    const Vec2& getSlipVelocity(const State& s) const {
        const MyHybridContactInfo& info = getContactInfo(s);
        const Vec2& v = info.v_HCb.getSubVec<2>(0);
        return v;
    }

    // The way we constructed the NoSlip1D constraints makes the multipliers be
    // the force on the half space on body P; we negate here so we'll get the 
    // force on the vertex body B instead.
    Vec2 getStictionForce(const State& s) const {
        assert(isSticking(s));
        return Vec2(-m_noslipX.getMultiplier(s), -m_noslipY.getMultiplier(s));
    }

    void recordImpendingSlipDir(const State& s) const {
        const Vec2 f = getStictionForce(s);
        const Real fMag = f.norm();
        const Vec2 dir = fMag==0 ? Vec2(1,0) : Vec2(-f/fMag);
        m_recordedSlipDir = dir;
    }

    void recordActualSlipDir(const State& s) const {
        const Vec2 v = getSlipVelocity(s);
        const Real vMag = v.norm();
        const Vec2 dir = vMag==0 ? Vec2(1,0) : Vec2(v/vMag);
        m_recordedSlipDir = dir;
    }

    void updatePrevSlipDirFromRecorded(State& s) const {
        setPrevSlipDir(s, m_recordedSlipDir);
    }

    Real calcSlipSpeedWitness(const State& s) const {
        if (getHeight(s) >= 0 || isSticking(s)) return 0;
        const Vec2& slipDirPrev = getPrevSlipDir(s);
        const Vec2& vNow = getSlipVelocity(s);
        if (slipDirPrev.isNaN()) return vNow.norm();
        return dot(vNow, slipDirPrev);
    }

    Real calcStictionForceWitness(const State& s) const {
        if (getHeight(s) >= 0 || !isSticking(s)) return 0;
        const Real mu_s = m_material.getStaticFriction();
        const Real N = getNormalForce(s);
        const Vec2 f = getStictionForce(s);
        const Real fmag = f.norm();
        return mu_s*N - fmag;
    }



    // Note that initializeForStiction() must have been called first.
    void enableStiction(State& s) const
    {   m_noslipX.setContactPoint(s, m_contactPointInP);
        m_noslipY.setContactPoint(s, m_contactPointInP);
        m_noslipX.enable(s); m_noslipY.enable(s); }

    void disableStiction(State& s) const
    {   m_noslipX.disable(s); m_noslipY.disable(s); }

    // Return height of vertex over plane; negative if penetrated. You can
    // call this at Stage::Position.
    Real findH(const State& state) const {
        const Vec3 V_P = findVInP(state); // location of V in P
        const Real h = dot(V_P-m_X_PH.p(), m_X_PH.z());
        return h;
    }

    // Returns height of vertex over plane same as findH(); also returns 
    // contact point (projection of V onto plane along Hz).
    Real findContactPointInP(const State& state, Vec3& contactPointInP) const {
        const Vec3 V_P = findVInP(state);       // location of V in P
        const Real h = dot(V_P-m_X_PH.p(), m_X_PH.z());
        contactPointInP = V_P - h*m_X_PH.z();
        return h;
    }

    // Calculate v_eff, the direction to be opposed by the sliding force.
    Vec2 getEffectiveSlipDir(const State& s, const Vec2& vSlip, Real vSlipMag) const {
        const Vec2 prevVslipDir = getPrevSlipDir(s);
        if (shouldUpdate(vSlip, vSlipMag, prevVslipDir, m_vtol))
            return vSlip/vSlipMag;
        else return prevVslipDir;
    }

    // Return the slip velocity as recorded at the end of the last time step. 
    const Vec2& getPrevSlipDir(const State& state) const {
        const Vec2& prevSlipDir = Value<Vec2>::downcast
           (m_forces.getDiscreteVariable(state, m_prevSlipDirIx));
        return prevSlipDir;
    }
    // Modify the discrete state variable directly.
    void setPrevSlipDir(State& state, const Vec2& slipDir) const {
        Vec2& prevSlipDir = Value<Vec2>::updDowncast
           (m_forces.updDiscreteVariable(state, m_prevSlipDirIx));
        prevSlipDir = slipDir;
        SimTK_DEBUG3("STATE CHG %d: prevDir to %g %g\n",
            m_index, slipDir[0], slipDir[1]);
    }

    // Get access to the auto-update cache variable that will update the
    // previous slip direction state variable if the step is accepted.
    Vec2& updPrevSlipDirAutoUpdateValue(const State& state) const {
        Vec2& prevSlipUpdate = Value<Vec2>::updDowncast
            (m_forces.updDiscreteVarUpdateValue(state, m_prevSlipDirIx));
        return prevSlipUpdate;
    }
    // Mark the auto-update cache entry valid.
    void markPrevSlipDirAutoUpdateValueValid(const State& state) const {
        m_forces.markDiscreteVarUpdateValueRealized(state, m_prevSlipDirIx);
    }

    const MyHybridContactInfo& getContactInfo(const State& state) const {
        const MyHybridContactInfo& info =
            Value<MyHybridContactInfo>::downcast
                (m_forces.getCacheEntry(state, m_contactInfoIx));
        return info;
    }

    MyHybridContactInfo& updContactInfo(const State& state) const {
        MyHybridContactInfo& info =
            Value<MyHybridContactInfo>::updDowncast
                (m_forces.updCacheEntry(state, m_contactInfoIx));
        return info;
    }
    //--------------------------------------------------------------------------
    //                       Custom force virtuals
    // Apply the normal force, and the sliding friction force if it is enabled.
    // This is called during realize(Dynamics).
    void calcForce(const State& state, Vector_<SpatialVec>& bodyForces, 
                   Vector_<Vec3>& particleForces, Vector& mobilityForces) const
                   OVERRIDE_11
    {
        const MyHybridContactInfo& info = getContactInfo(state);
        if (info.h >= 0) 
            return; // no contact

        const Vec3 f_GCb = info.R_GH * info.f_HCb; // re-express in Ground
        getBodyB().applyForceToBodyPoint(state, info.Cb,  f_GCb, bodyForces);
        getBodyP().applyForceToBodyPoint(state, info.Cp, -f_GCb, bodyForces);
    }

    // The normal force stores energy as 2/5 k h^(5/2) when h<0.
    Real calcPotentialEnergy(const State& state) const OVERRIDE_11 {
        const MyHybridContactInfo& info = getContactInfo(state);
        if (info.h >= 0) 
            return 0; // no contact
        const Real h52 = square(info.h)*std::sqrt(-info.h);
        const Real k = m_material.getStiffness();
        return 0.4*k*h52;
    }

    // Allocate state variable for storing the previous sliding direction.
    void realizeTopology(State& state) const OVERRIDE_11 {
        // The previous sliding direction is used in an event witness that 
        // is evaluated at Velocity stage.
        m_prevSlipDirIx = m_forces.allocateAutoUpdateDiscreteVariable
           (state, Stage::Velocity, new Value<Vec2>(Vec2(NaN)), 
            Stage::Velocity);

        m_contactInfoIx = m_forces.allocateCacheEntry
           (state, Stage::Velocity, new Value<MyHybridContactInfo>());
    }


    // Calculate everything here and save in contact info cache entry where
    // it can be retrieved for generating forces, reporting, etc.
    void realizeVelocity(const State& state) const OVERRIDE_11 {
        MyHybridContactInfo& info = updContactInfo(state);

        // Forces generated only if h<0. Cp always be the projection of the
        // vertex onto the halfspace surface.
        info.h = findContactPointInP(state, info.Cp);

        const Rotation& R_PH = m_X_PH.R();
        const Rotation& R_GP = getBodyP().getBodyRotation(state);
        info.R_GH = R_GP*R_PH;

        // Station of B coincident with the contact point.
        info.Cb = info.h < 0 
                    ? findPointOfBCoincidentWithPointOfP(state, info.Cp)
                    : getVertexOnB();

        // Velocity of B's contact station in P.
        const Vec3 v_PCb = getBodyB().findStationVelocityInAnotherBody
                                                (state, info.Cb, getBodyP());
        info.v_HCb = ~R_PH*v_PCb; // re-express in H
        const Real hdot = info.v_HCb[2]; // z component
        const Vec2 vSlip_HC(info.v_HCb[0], info.v_HCb[1]);
        info.vSlipMag = vSlip_HC.norm();

        const Real k    = m_material.getStiffness(), 
                   c    = m_material.getDissipation(),
                   mu_d = m_material.getDynamicFriction(),
                   mu_s = m_material.getStaticFriction(),
                   mu_v = m_material.getViscousFriction();

        info.f_HCb = Vec3(0);
        if (info.h >= 0)
            return; // no contact

        const Real h32 = -info.h*sqrt(-info.h); // |h| ^ (3/2)
        const Real N = std::max(Real(0), k*h32*(1-c*hdot));
        if (N==0) 
            return; // no contact force

        // N is the Hz component of the force on Cb.
        Vec2 fT_HCb(0); // This will be the (Hx,Hy) components of force on Cb.
        #ifdef USE_CONTINUOUS_STICTION
        {
            // Make v unitless velocity and scale viscous coefficient to match.
            const Real v = info.vSlipMag / m_vtrans;
            const Real mu=stribeck(mu_s,mu_d,mu_v*m_vtrans,v);
            const Real T = mu*N;
            fT_HCb = (-T/info.vSlipMag)*vSlip_HC;
        }
        #else
        if (!isSticking(state)) {
            // Apply sliding force 
            const Real T = (mu_d + mu_v*info.vSlipMag)*N;
            fT_HCb = -T*getEffectiveSlipDir(state, vSlip_HC, // in Hxy
                                            info.vSlipMag); 
        }
        #endif

        info.f_HCb = Vec3(fT_HCb[0], fT_HCb[1], N); // force on Cb, in H
    }

    // If we're sticking, set the update value for the previous slip direction
    // to the opposite of the stiction force direction.
    // If we're sliding, set the update value for the previous slip direction
    // if the current slip velocity is usable.
    #ifndef USE_CONTINUOUS_STICTION
    void realizeAcceleration(const State& state) const OVERRIDE_11 {
        const MyHybridContactInfo& info = getContactInfo(state);
        const Vec2& prevSlipDir = getPrevSlipDir(state);

        if (info.h >= 0) {
            // No contact. Forget previous slip direction.
            if (!prevSlipDir.isNaN()) {
                SimTK_DEBUG1("%d LIFTOFF UPDATE, forget prevSlipDir\n", m_index);
                updPrevSlipDirAutoUpdateValue(state).setToNaN();
                markPrevSlipDirAutoUpdateValueValid(state);
            }
            return;
        }

        // Sticking.
        if (isSticking(state)) {
            const Vec2 f_HCb = getStictionForce(state);
            const Real fMag = f_HCb.norm();
            if (fMag > 0) {
                Vec2& prevSlipUpdate = updPrevSlipDirAutoUpdateValue(state);
                prevSlipUpdate = -f_HCb / fMag;
                #ifndef NDEBUG
                printf("%d STICKING UPDATE: prevSlipDir=%g %g; now=%g %g\n",
                    m_index, getPrevSlipDir(state)[0],getPrevSlipDir(state)[1],
                    prevSlipUpdate[0], prevSlipUpdate[1]);
                #endif
                markPrevSlipDirAutoUpdateValueValid(state);
            }
            return;
        }

        // Sliding.
        const Vec2& vSlip_HCb = info.v_HCb.getSubVec<2>(0); // x,y

        if (shouldUpdate(vSlip_HCb, info.vSlipMag, prevSlipDir, m_vtol)) {
            Vec2& prevSlipUpdate = updPrevSlipDirAutoUpdateValue(state);
            prevSlipUpdate = vSlip_HCb / info.vSlipMag;
            markPrevSlipDirAutoUpdateValueValid(state);

            #ifndef NDEBUG
            printf("%d SLIDING UPDATE: prevSlipDir=%g %g; now=%g %g; |v|=%g dot=%g vdot=%g\n",
                m_index, prevSlipDir[0],prevSlipDir[1],
                prevSlipUpdate[0],prevSlipUpdate[1], info.vSlipMag, 
                ~prevSlipUpdate*prevSlipDir, ~vSlip_HCb*prevSlipDir);
            #endif
        } else {
            #ifndef NDEBUG
            printf("%d SLIDING; NO UPDATE: prevSlipDir=%g %g; Vnow=%g %g; |v|=%g vdot=%g\n",
                m_index,
                prevSlipDir[0],prevSlipDir[1],vSlip_HCb[0],vSlip_HCb[1],
                info.vSlipMag, ~vSlip_HCb*prevSlipDir);
            #endif
        }
    }
    #endif

    std::ostream& writeFrictionInfo(const State& s, const String& indent, 
                                    std::ostream& o) const 
    {
        o << indent;
        if (!isInContact(s)) o << "OFF";
        else if (isSticking(s)) o << "STICK";
        else o << "SLIP";

        const Vec2 v = getSlipVelocity(s);
        const Vec2 pd = getPrevSlipDir(s);
        const Vec2 f = isSticking(s) ? getStictionForce(s)
                                     : getSlidingForce(s);
        o << " prevDir=" << ~pd << " V=" << ~v << " Vdot=" << ~v*pd 
          << " F=" << ~f << endl;
        return o;
    }


    void showContactForces(const State& s, Array_<DecorativeGeometry>& geometry) 
        const
    {
        const Real Scale = 0.01;
        const Real NH = getNormalForce(s);

        #ifdef USE_CONTINUOUS_STICTION
        const bool isInStiction = getActualSlipSpeed(s) <= m_vtrans;
        const Vec2 fH = getSlidingForce(s);
        #else
        const bool isInStiction = isSticking(s);
        const Vec2 fH = isSticking(s) ? getStictionForce(s) : getSlidingForce(s);
        #endif

        if (fH.normSqr() < square(SignificantReal) && NH < SignificantReal)
            return;

        const MyHybridContactInfo& info = getContactInfo(s);
        const Vec3 fG = info.R_GH * Vec3(fH[0],fH[1],0); // friction
        const Vec3 NG = info.R_GH * Vec3(0, 0, NH); // normal

        const MobilizedBody& bodyB = getBodyB();
        const Vec3& stationB = getVertexOnB();
        const Vec3 stationG = bodyB.getBodyTransform(s)*stationB;
        const Vec3 endfG = stationG - Scale*fG;
        const Vec3 endNG = stationG + Scale*NG;
        geometry.push_back(DecorativeLine(endfG    + Vec3(0,.05,0),
                                          stationG + Vec3(0,.05,0))
                            .setColor(isInStiction?Green:Orange));
        geometry.push_back(DecorativeLine(endNG    + Vec3(0,.05,0),
                                          stationG + Vec3(0,.05,0))
                            .setColor(Red));
    }

    int getIndex() const {return m_index;}
    const MobilizedBody& getBodyB() const 
    {   return m_matter.getMobilizedBody(m_vmobodx); }
    const Vec3& getVertexOnB() const {return m_vertex_B;}
    const MobilizedBody& getBodyP() const 
    {   return m_matter.getMobilizedBody(m_hsmobodx); }
    const Transform& getX_PH() const {return m_X_PH;}
    const ContactMaterial& getMaterial() const {return m_material;}

    //--------------------------------------------------------------------------
private:
    // Determine whether the current slip velocity is reliable enough that
    // we should use it to replace the previous slip velocity.
    static bool shouldUpdate(const Vec2& vSlip, Real vSlipMag, 
                             const Vec2& prevSlipDir, Real velTol) {
        if (prevSlipDir.isNaN())
            return vSlipMag > 0; // we'll take anything

        // Check for reversal.
        bool reversed = (~vSlip*prevSlipDir < 0);
        return !reversed && (vSlipMag > velTol);
    }

    // Find the location of the vertex on B, measured and expressed in P.
    Vec3 findVInP(const State& s) const {
        return getBodyB().findStationLocationInAnotherBody
                                                (s, m_vertex_B, getBodyP());
    }

    // Find the location and velocity of the vertex on B, measured from and
    // expressed in P.
    Vec3 findVDotInP(const State& s) const {
        return getBodyB().findStationVelocityInAnotherBody
            (s, m_vertex_B, getBodyP());
    }


    Vec3 findPointOfBCoincidentWithPointOfP
       (const State& s, const Vec3& r_P) const
    {
        return getBodyP().findStationLocationInAnotherBody(s, r_P, getBodyB());
    }


private:
    const GeneralForceSubsystem&    m_forces;
    const SimbodyMatterSubsystem&   m_matter;
    const MobilizedBodyIndex        m_hsmobodx; // body P with halfspace
    Transform                       m_X_PH;     // halfspace frame in P
    const MobilizedBodyIndex        m_vmobodx;  // body B with vertex
    Vec3                            m_vertex_B; // vertex location in B
    ContactMaterial                 m_material; // composite material props

    Real                            m_vtol;     // velocity tolerance
    Real                            m_vtrans;   // transition velocity for
                                                //   Stribeck stiction

    Constraint::NoSlip1D            m_noslipX;
    Constraint::NoSlip1D            m_noslipY;

    int                             m_index;    // unique id for this contact

    // This is recorded an Position stage prior to turning on stiction
    // constraints; then we use it to set the contact point in noslipX,Y.
    Vec3                            m_contactPointInP;

    // This is recorded during event handling and then used to set the
    // remembered previous slip direction at the end of the handler to make
    // sure all witness functions are positive then.
    mutable Vec2                    m_recordedSlipDir;

    // Set during realizeTopology().
    mutable DiscreteVariableIndex   m_prevSlipDirIx; // previous slip direction
    mutable CacheEntryIndex         m_contactInfoIx; 
};


//==============================================================================
//                       MY UNILATERAL CONSTRAINT SET
//==============================================================================

// These are indices into the unilateral constraint set arrays.
struct MyElementSubset {
    void clear() 
    {   m_contact.clear();m_sliding.clear();m_stiction.clear();
        m_ineligible.clear(); }
    Array_<int> m_contact;
    Array_<int> m_sliding;  // subset of above that can only be sliding
    Array_<int> m_stiction; // subset of above that might be in stiction
    Array_<int> m_ineligible; // stiction on, but not in contact any more
};

class MyUnilateralConstraintSet {
public:
    // Transition velocity: if a slip velocity is smaller than this the
    // contact is a candidate for stiction.
    MyUnilateralConstraintSet(const MultibodySystem& mbs, Real transitionVelocity)
    :   m_mbs(mbs), m_transitionVelocity(transitionVelocity) {}

    // Ownership of this force element belongs to the System; we're just keeping
    // a reference to it here.
    int addHybridElement(MyHybridVertexContactElementImpl* vertex) {
        const int index = (int)m_hybrid.size();
        m_hybrid.push_back(vertex);
        vertex->setContactIndex(index);
        vertex->setTransitionVelocity(m_transitionVelocity);
        return index;
    }

    Real getTransitionVelocity() const {return m_transitionVelocity;}
    void setTransitionVelocity(Real v) {m_transitionVelocity=v;}

    int getNumContactElements() const {return (int)m_hybrid.size();}
    const MyHybridVertexContactElementImpl& getContactElement(int ix) const 
    {   return *m_hybrid[ix]; }
    MyHybridVertexContactElementImpl& updContactElement(int ix) 
    {   return *m_hybrid[ix]; }

    // Initialize all runtime fields in the contact & friction elements.
    void initialize()
    {
        for (unsigned i=0; i < m_hybrid.size(); ++i)
            m_hybrid[i]->initialize();
    }

    // Return the contact and friction elements that might be involved in 
    // generating contact forces at the current state. Candidate contact
    // elements are those that are (a) already enabled, or (b) for which 
    // perr <= posTol and verr <= velTol. Candidate friction elements are those
    // whose normal force master is unconditional or a candidate and (a) which 
    // are already sticking, or (b) for which vslip <= velTol, or (c) for which
    // vslip opposes the previous slip direction, meaning it has reversed and 
    // must have passed through zero during the last step. These are the elements 
    // that can be activated without making any changes to the configuration or 
    // velocity state variables, except slightly for constraint projection. 
    //
    // We also record the friction elements that, if their masters are active, 
    // can only slide because they have a significant slip velocity. State must 
    // be realized through Velocity stage.
    void findCandidateElements(const State&     s,
                               Real             velTol,
                               MyElementSubset& candidates) const
    {
        candidates.clear();
        for (unsigned i=0; i < m_hybrid.size(); ++i) {
            if (!m_hybrid[i]->isInContact(s)) {
                if (m_hybrid[i]->isSticking(s)) {
                    candidates.m_ineligible.push_back(i); // must disable
                    SimTK_DEBUG2("%d NOW INELIGIBLE because h=%g.\n",
                        i, m_hybrid[i]->findH(s));
                }
                continue;
            }

            candidates.m_contact.push_back(i);

            if (m_hybrid[i]->isSticking(s) 
                || m_hybrid[i]->getActualSlipSpeed(s) <= velTol
                || m_hybrid[i]->calcSlipSpeedWitness(s) <= 0) 
            {
                m_hybrid[i]->initializeForStiction(s);
                candidates.m_stiction.push_back(i); // could stick or slide
            } else
                candidates.m_sliding.push_back(i); // can only slide
        }
    }

    // Look through the given constraint subset and enable any constraints
    // that are currently disabled, and disable any constraints that are still
    // on after liftoff. Returns true if any change was made.
    bool enableConstraintSubset(const MyElementSubset& subset,
                                State&                 state) const
    {
        bool changedSomething = false;

        // Disable all ineligible constraints.
        for (unsigned i=0; i < subset.m_ineligible.size(); ++i) {
            const int which = subset.m_ineligible[i];
            const MyHybridVertexContactElementImpl& fric = 
                getContactElement(which);
            SimTK_DEBUG1("%d DISABLING INELIGIBLE STICTION.\n", i);
            fric.disableStiction(state);
            changedSomething = true;
       }

        // Enable all stiction constraints.
        for (unsigned i=0; i < subset.m_stiction.size(); ++i) {
            const int which = subset.m_stiction[i];
            const MyHybridVertexContactElementImpl& fric = 
                getContactElement(which);
            if (!fric.isSticking(state)) {
                SimTK_DEBUG1("%d ENABLING CANDIDATE STICTION.\n", i);
                fric.enableStiction(state);
                changedSomething = true;
            }
        }

        m_mbs.realize(state, Stage::Instance);
        return changedSomething;
    }

    // All event handlers call this method before returning. Given a state for
    // which no (further) impulse is required, here we decide which contact and
    // stiction constraints are active, and ensure that they satisfy the 
    // required constraint tolerances to the given accuracy. For sliding 
    // contacts, we will have recorded the slip or impending slip direction and 
    // converged the normal forces.
    // TODO: in future this may return indicating that an impulse is required
    // after all, as in Painleve's paradox.
    void selectActiveConstraints(State& state, Real accuracy) const;

    // This is the inner loop of selectActiveConstraints(). Given a set of
    // candidates to consider, it finds an active subset and enables those
    // constraints.
    void findActiveCandidates(State&                 state, 
                              const MyElementSubset& candidates) const;

    // In Debug mode, produce a useful summary of the current state of the
    // contact and friction elements.
    void showConstraintStatus(const State& state, const String& place) const;

    ~MyUnilateralConstraintSet() {
        // Nothing to delete since we are only holding references.
    }

    const MultibodySystem& getMultibodySystem() const {return m_mbs;}
private:
    const MultibodySystem&                      m_mbs;
    Real                                        m_transitionVelocity;
    Array_<MyHybridVertexContactElementImpl*>   m_hybrid; // unowned ref
};



//==============================================================================
//                               STATE SAVER
//==============================================================================
// This reporter is called now and again to save the current state so we can
// play back a movie at the end.
class StateSaver : public PeriodicEventReporter {
public:
    StateSaver(const MultibodySystem&                   mbs,
               const MyUnilateralConstraintSet&         unis,
               const Integrator&                        integ,
               Real                                     reportInterval)
    :   PeriodicEventReporter(reportInterval), 
        m_mbs(mbs), m_unis(unis), m_integ(integ) 
    {   m_states.reserve(2000); }

    ~StateSaver() {}

    void clear() {m_states.clear();}
    int getNumSavedStates() const {return (int)m_states.size();}
    const State& getState(int n) const {return m_states[n];}

    void handleEvent(const State& s) const {
        const SimbodyMatterSubsystem& matter=m_mbs.getMatterSubsystem();
        const SpatialVec PG = matter.calcSystemMomentumAboutGroundOrigin(s);
        m_mbs.realize(s, Stage::Acceleration);

#ifndef NDEBUG
        printf("%3d: %5g mom=%g,%g E=%g", m_integ.getNumStepsTaken(),
            s.getTime(),
            PG[0].norm(), PG[1].norm(), m_mbs.calcEnergy(s));
        cout << " Triggers=" << s.getEventTriggers() << endl;
        m_unis.showConstraintStatus(s, "STATE SAVER");
#endif

        m_states.push_back(s);
    }
private:
    const MultibodySystem&                  m_mbs;
    const MyUnilateralConstraintSet&        m_unis;
    const Integrator&                       m_integ;
    mutable Array_<State>                   m_states;
};

// A periodic event reporter that does nothing; useful for exploring the
// effects of interrupting the simulation.
class Nada : public PeriodicEventReporter {
public:
    explicit Nada(Real reportInterval)
    :   PeriodicEventReporter(reportInterval) {} 

    void handleEvent(const State& s) const {
#ifndef NDEBUG
        printf("%7g NADA\n", s.getTime());
#endif
    }
};


//==============================================================================
//                            SHOW CONTACT
//==============================================================================
// For each visualization frame, generate ephemeral geometry to show just 
// during this frame, based on the current State.
class ShowContact : public DecorationGenerator {
public:
    ShowContact(const MyUnilateralConstraintSet& unis) 
    :   m_unis(unis) {}

    void generateDecorations(const State&                state, 
                             Array_<DecorativeGeometry>& geometry) OVERRIDE_11
    {
        for (int i=0; i < m_unis.getNumContactElements(); ++i) {
            const MyHybridVertexContactElementImpl& contact = 
                m_unis.getContactElement(i);
            const Vec3 loc = contact.whereToDisplay(state);
            if (contact.isInContact(state)) {
                geometry.push_back(DecorativeSphere(.1)
                    .setTransform(loc)
                    .setColor(Red).setOpacity(.25));
                String text = "LOCKED";
                text = contact.isSticking(state) ? "STICKING"
                                                    : "CONTACT";
                m_unis.getMultibodySystem().realize(state, Stage::Acceleration);
                contact.showContactForces(state, geometry);
                geometry.push_back(DecorativeText(String(i)+"-"+text)
                    .setColor(White).setScale(.1)
                    .setTransform(loc+Vec3(0,.04,0)));
            } else {
                geometry.push_back(DecorativeText(String(i))
                    .setColor(White).setScale(.1)
                    .setTransform(loc+Vec3(0,.02,0)));
            }
        }
    }
private:
    const MyUnilateralConstraintSet& m_unis;
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
        const Vec3 p_GC = Bo + Vec3(0, 1, 5); // above and back
        const Rotation R_GC(UnitVec3(0,1,0), YAxis,
                            p_GC-Bo, ZAxis);
        viz.setCameraTransform(Transform(R_GC,p_GC));
        //viz.pointCameraAt(Bo, Vec3(0,1,0));
    }
private:
    const MobilizedBody m_body;
};

//==============================================================================
//                          STICTION ON HANDLER
//==============================================================================
// Allocate one of these for each contact constraint that has friction. This 
// handler takes care of turning on the stiction constraints when the sliding 
// velocity drops to zero.
class StictionOn: public TriggeredEventHandler {
public:
    StictionOn(const MultibodySystem&       system,
        const MyUnilateralConstraintSet&    unis,
        unsigned                            which) 
    :   TriggeredEventHandler(Stage::Velocity), 
        m_mbs(system), m_unis(unis), m_which(which)
    { 
        getTriggerInfo().setTriggerOnRisingSignTransition(false);
    }

    // This is the witness function. It is positive as long as we continue
    // to slide in the same direction; negative means reversal.
    Real getValue(const State& state) const {
        const MyHybridVertexContactElementImpl& contact = 
            m_unis.getContactElement(m_which);
        if (!contact.isInContact(state)) return 0;
        const Real signedSlipSpeed = contact.calcSlipSpeedWitness(state);
        return signedSlipSpeed;
    }

    void handleEvent
       (State& s, Real accuracy, bool& shouldTerminate) const 
    {
        SimTK_DEBUG2("\nhandle %d slide->stick@%.17g\n", m_which, s.getTime());
        SimTK_DEBUG("\n----------------------------------------------------\n");
        SimTK_DEBUG2("STICTION ON triggered by friction element %d @t=%.15g\n", 
            m_which, s.getTime());
        m_mbs.realize(s, Stage::Acceleration);

        #ifndef NDEBUG
        m_unis.showConstraintStatus(s, "ENTER STICTION ON");
        cout << " entry triggers=" << s.getEventTriggers() << "\n";
        #endif

        m_unis.selectActiveConstraints(s, accuracy);

        #ifndef NDEBUG
        m_mbs.realize(s, Stage::Acceleration);
        cout << " exit triggers=" << s.getEventTriggers() << "\n";
        #endif

        SimTK_DEBUG("STICTION ON done.\n");
        SimTK_DEBUG("----------------------------------------------------\n");
    }

private:
    const MultibodySystem&              m_mbs; 
    const MyUnilateralConstraintSet&    m_unis;
    const int                           m_which; // one of the contact elements
};



//==============================================================================
//                          STICTION OFF HANDLER
//==============================================================================
// Allocate one of these for each contact constraint. This handler takes
// care of turning off stiction constraints when the stiction force magnitude
// exceeds mu*N.
class StictionOff: public TriggeredEventHandler {
public:
    StictionOff(const MultibodySystem&      system,
        const MyUnilateralConstraintSet&    unis,
        unsigned                            which) 
    :   TriggeredEventHandler(Stage::Acceleration), 
        m_mbs(system), m_unis(unis), m_which(which)
    {
        getTriggerInfo().setTriggerOnRisingSignTransition(false);
    }

    // This is the witness function. It is positive as long as mu_s*N is greater
    // than the friction force magnitude.
    Real getValue(const State& state) const {
        const MyHybridVertexContactElementImpl& contact = 
            m_unis.getContactElement(m_which);
        if (!contact.isInContact(state)) return 0;
        const Real capacity = contact.calcStictionForceWitness(state);
        return capacity; // how much stiction capacity is left
    }

    void handleEvent
       (State& s, Real accuracy, bool& shouldTerminate) const 
    {
        SimTK_DEBUG2("\nhandle %d stick->slide@%.17g\n", m_which, s.getTime());
        SimTK_DEBUG("\n----------------------------------------------------\n");
        SimTK_DEBUG2("STICTION OFF triggered by friction element %d @t=%.15g\n", 
            m_which, s.getTime());
        m_mbs.realize(s, Stage::Acceleration);

        #ifndef NDEBUG
        m_unis.showConstraintStatus(s, "ENTER STICTION OFF");
        cout << " triggers=" << s.getEventTriggers() << "\n";
        #endif

        m_unis.selectActiveConstraints(s, accuracy);

        #ifndef NDEBUG
        m_mbs.realize(s, Stage::Acceleration);
        cout << " exit triggers=" << s.getEventTriggers() << "\n";
        #endif

        SimTK_DEBUG("STICTION OFF done.\n");
        SimTK_DEBUG("----------------------------------------------------\n");
    }

private:
    const MultibodySystem&              m_mbs; 
    const MyUnilateralConstraintSet&    m_unis;
    const int                           m_which; // one of the friction elements
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

        //SimTK_DEBUG4("PUSHING @t=%g (%g,%g,%g)\n", state.getTime(),
        //    m_forceG[0], m_forceG[1], m_forceG[2]);
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
//                                   MAIN
//==============================================================================
int main(int argc, char** argv) {
  try { // If anything goes wrong, an exception will be thrown.
    const Real ReportInterval=1./30;
    const Vec3 BrickHalfDims(.1, .25, .5);
    const Real Radius = BrickHalfDims[0]/3;
    const Real BrickMass = 10;
    #ifdef USE_TIMS_PARAMS
        const Real RunTime=16;  // Tim's time
        const Real Stiffness = 2e7;
        const Real Dissipation = 0.1;
        const Real CoefRest = 0; 
        const Real mu_d = 1; /* compliant: .7*/
        const Real mu_s = 1; /* compliant: .7*/
        const Real mu_v = /*0.05*/0;
        const Real TransitionVelocity = 0.01;
        const Inertia brickInertia(.1,.1,.1);
    #else
        const Real RunTime=20;
        const Real Stiffness = 1e6;
        const Real CoefRest = 0; 
        const Real TargetVelocity = 3; // speed at which to match coef rest
//        const Real Dissipation = (1-CoefRest)/TargetVelocity;
        const Real Dissipation = 1;
        const Real mu_d = .5;
        const Real mu_s = .8;
        const Real mu_v = 0*0.05;
        const Real TransitionVelocity = 0.01;
        const Inertia brickInertia(BrickMass*UnitInertia::brick(BrickHalfDims));
    #endif

    printf("\n******************** Tim's Box Hybrid ********************\n");
    #ifdef USE_CONTINUOUS_STICTION
    printf("USING OLD MODEL: Continuous Stiction (Stribeck)\n");
    #else
    printf("USING NEW MODEL: Hybrid Compliant material/rigid stiction\n");
    #endif
    #ifdef USE_TIMS_PARAMS
    printf("Using Tim's parameters:\n");
    #else
    printf("Using Sherm's parameters:\n");
    #endif
    printf("  stiffness=%g dissipation=%g\n", Stiffness, Dissipation);
    printf("  mu_d=%g mu_s=%g mu_v=%g\n", mu_d, mu_s, mu_v);
    printf("  transition velocity=%g\n", TransitionVelocity);
    printf("  brick inertia=%g %g %g\n",
        brickInertia.getMoments()[0], brickInertia.getMoments()[1], 
        brickInertia.getMoments()[2]); 
    printf("******************** Tim's Box Hybrid ********************\n\n");

        // CREATE MULTIBODY SYSTEM AND ITS SUBSYSTEMS
    MultibodySystem             mbs;
    SimbodyMatterSubsystem      matter(mbs);
    GeneralForceSubsystem       forces(mbs);
    Force::Gravity              gravity(forces, matter, -YAxis, 9.81);

    MobilizedBody& Ground = matter.updGround();
    // Define a material to use for contact. This is not very stiff.
    ContactMaterial material(std::sqrt(Radius)*Stiffness,
                             Dissipation,
                             mu_s,  // mu_static
                             mu_d,  // mu_dynamic
                             mu_v); // mu_viscous

        // ADD MOBILIZED BODIES AND CONTACT CONSTRAINTS

    MyUnilateralConstraintSet unis(mbs, TransitionVelocity);

    Body::Rigid brickBody = 
        Body::Rigid(MassProperties(BrickMass, Vec3(0), brickInertia));

    // First body: cube
    MobilizedBody::Cartesian loc(Ground, MassProperties(0,Vec3(0),Inertia(0)));
    MobilizedBody::Ball brick(loc, Vec3(0),
                             brickBody, Vec3(0));
    brick.addBodyDecoration(Transform(), DecorativeBrick(BrickHalfDims)
                                                .setColor(Red).setOpacity(.3));
/*
1) t= 0.5, dt = 2 sec, pt = (0.05, 0.2, 0.4), fdir = (1,0,0), mag = 50N
2) t= 4.0, dt = 0.5 sec, pt = (0.03, 0.06, 0.09), fdir = (0.2,0.8,0), mag = 300N
3) t= 0.9, dt = 2 sec, pt = (0,0,0), fdir = (0,1,0), mag = 49.333N (half the weight of the block)
4) t= 13.0, dt = 1 sec, pt = (0 0 0), fdir = (-1,0,0), mag = 200N
*/
    Force::Custom(forces, new MyPushForceImpl(brick, Vec3(0.05,0.2,0.4),
                                                    50 * Vec3(1,0,0),
                                                    0.5, 0.5+2));
    Force::Custom(forces, new MyPushForceImpl(brick, Vec3(0.03, 0.06, 0.09),
                                                    300 * UnitVec3(0.2,0.8,0),
                                                    //300 * Vec3(0.2,0.8,0),
                                                    4, 4+0.5));
    Force::Custom(forces, new MyPushForceImpl(brick, Vec3(0),
                                                    49.333 * Vec3(0,1,0),
                                                    0.9, 0.9+2));
    Force::Custom(forces, new MyPushForceImpl(brick, Vec3(0),
                                                    200 * Vec3(-1,0,0),
                                                    13, 13+1));

    #ifndef USE_TIMS_PARAMS
    // Extra late force.
    Force::Custom(forces, new MyPushForceImpl(brick, Vec3(.1, 0, .45),
                                                    20 * Vec3(-1,-1,.5),
                                                    15, Infinity));
    #endif

    for (int i=-1; i<=1; i+=2)
    for (int j=-1; j<=1; j+=2)
    for (int k=-1; k<=1; k+=2) {
        const Vec3 pt = Vec3(i,j,k).elementwiseMultiply(BrickHalfDims);
        MyHybridVertexContactElementImpl* vertex =
            new MyHybridVertexContactElementImpl(forces,
                Ground, YAxis, 0, // halfplane
                brick, pt, material);
        Force::Custom(forces, vertex); // add force element to system
        unis.addHybridElement(vertex); // assign index, transition velocity
        #ifndef USE_CONTINUOUS_STICTION
        mbs.addEventHandler(new StictionOn(mbs, unis, vertex->getIndex()));
        mbs.addEventHandler(new StictionOff(mbs, unis, vertex->getIndex()));
        #endif
    }

    matter.setShowDefaultGeometry(false);
    Visualizer viz(mbs);
    //viz.setDesiredFrameRate(1000);
    viz.setShowSimTime(true);
    viz.setShowFrameNumber(true);
    viz.setShowFrameRate(true);
    viz.addDecorationGenerator(new ShowContact(unis));

    #ifdef ANIMATE
    mbs.addEventReporter(new Visualizer::Reporter(viz, ReportInterval));
    #else
    // This does nothing but interrupt the simulation.
    mbs.addEventReporter(new Nada(ReportInterval));
    #endif

    viz.addFrameController(new BodyWatcher(brick));

    Vec3 cameraPos(0, 1, 2);
    UnitVec3 cameraZ(0,0,1);
    viz.setCameraTransform(Transform(Rotation(cameraZ, ZAxis, 
                                              UnitVec3(YAxis), YAxis), 
                                     cameraPos));
    viz.pointCameraAt(Vec3(0,0,0), Vec3(0,1,0));

    #ifdef USE_TIMS_PARAMS
    Real accuracy = 1e-4;
    RungeKuttaMersonIntegrator integ(mbs);
    #else
    Real accuracy = 1e-2;
    RungeKutta3Integrator integ(mbs);
    //RungeKuttaMersonIntegrator integ(mbs);
    #endif

    integ.setAccuracy(accuracy);
    //integ.setMaximumStepSize(0.25);
    integ.setMaximumStepSize(0.1);
    //integ.setMaximumStepSize(0.001);

    StateSaver* stateSaver = new StateSaver(mbs,unis,integ,ReportInterval);
    mbs.addEventReporter(stateSaver);

    State s = mbs.realizeTopology(); // returns a reference to the the default state
    
    //matter.setUseEulerAngles(s, true);
    
    mbs.realizeModel(s); // define appropriate states for this System
    mbs.realize(s, Stage::Instance); // instantiate constraints if any


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
    loc.setQToFitTranslation(s, Vec3(0,10,0));
    loc.setUToFitLinearVelocity(s, Vec3(0,0,0));
    #else
    loc.setQToFitTranslation(s, Vec3(0,1.4,0));
    loc.setUToFitLinearVelocity(s, Vec3(10,0,0));
    const Rotation R_BC(SimTK::BodyRotationSequence,
                                0.7, XAxis, 0.6, YAxis, 0.5, ZAxis);
    brick.setQToFitRotation(s, R_BC);
    brick.setUToFitAngularVelocity(s, Vec3(1,0,.2));
    #endif

    mbs.realize(s, Stage::Velocity);
    viz.report(s);

    cout << "Initial configuration shown. Next? ";
    getchar();

    Assembler(mbs).setErrorTolerance(1e-6).assemble(s);
    viz.report(s);
    cout << "Assembled configuration shown. Ready? ";
    getchar();

    // Now look for enabled contacts that aren't sliding; turn on stiction
    // for those.
    mbs.realize(s, Stage::Velocity);
    Array_<int> enableTheseStictions;
    for (int i=0; i < unis.getNumContactElements(); ++i) {
        MyHybridVertexContactElementImpl& fric = unis.updContactElement(i);
        if (!fric.isInContact(s)) continue;
        const Real vSlip = fric.getActualSlipSpeed(s);
        fric.initializeForStiction(s); // just in case
        printf("friction element %d has v_slip=%g%s\n", i, vSlip,
            vSlip==0?" (ENABLING STICTION)":"");
        if (vSlip == 0)
            enableTheseStictions.push_back(i);
    }
    for (unsigned i=0; i < enableTheseStictions.size(); ++i)
        unis.getContactElement(enableTheseStictions[i]).enableStiction(s);

    // Make sure Lapack gets initialized.
    Matrix M(1,1); M(0,0)=1.23;
    FactorLU Mlu(M);

    
    // Simulate it.

    integ.setReturnEveryInternalStep(true);
    //integ.setAllowInterpolation(false);
    TimeStepper ts(mbs, integ);
    ts.setReportAllSignificantStates(true);

    const int NTries=1;
    Array_< Array_<State> > states(NTries);
    Array_< Array_<Real> > timeDiff(NTries-1);
    Array_< Array_<Vector> > yDiff(NTries-1);
    cout.precision(18);
    for (int ntry=0; ntry < NTries; ++ntry) {
        mbs.resetAllCountersToZero();
        unis.initialize(); // reinitialize
        ts.updIntegrator().resetAllStatistics();
        ts.initialize(s);
        int nStepsWithEvent = 0;

        const double startReal = realTime();
        const double startCPU = cpuTime();

        Integrator::SuccessfulStepStatus status;
        do {
            status=ts.stepTo(RunTime);
            //states[ntry].push_back(ts.getState());
            const int j = states[ntry].size()-1;
            if (ntry>0) {
                int prev = ntry-1;
                //prev=0;
                Real dt = states[ntry][j].getTime() 
                          - states[prev][j].getTime();
                volatile double vd;
                const int ny = states[ntry][j].getNY();
                Vector d(ny);
                for (int i=0; i<ny; ++i) {
                    vd = states[ntry][j].getY()[i] 
                           - states[prev][j].getY()[i];
                    d[i] = vd;
                }
                timeDiff[prev].push_back(dt);
                yDiff[prev].push_back(d);
                cout << "\n" << states[prev][j].getTime() 
                     << "+" << dt << ":" << d << endl;
            }

            const bool handledEvent = status==Integrator::ReachedEventTrigger;
            if (handledEvent) {
                ++nStepsWithEvent;
                SimTK_DEBUG3("EVENT %3d @%.17g status=%s\n\n", 
                    nStepsWithEvent, ts.getTime(),
                    Integrator::getSuccessfulStepStatusString(status).c_str());
            } else {
                SimTK_DEBUG3("step  %3d @%.17g status=%s\n", integ.getNumStepsTaken(),
                    ts.getTime(),
                    Integrator::getSuccessfulStepStatusString(status).c_str());
            }
            #ifndef NDEBUG
                viz.report(ts.getState());
            #endif
            //stateSaver->handleEvent(ts.getState());
        } while (ts.getTime() < RunTime);


        const double timeInSec = realTime()-startReal;
        const double cpuInSec = cpuTime()-startCPU;
        const int evals = integ.getNumRealizations();
        cout << "Done -- took " << integ.getNumStepsTaken() << " steps in " <<
            timeInSec << "s for " << ts.getTime() << "s sim (avg step=" 
            << (1000*ts.getTime())/integ.getNumStepsTaken() << "ms) " 
            << (1000*ts.getTime())/evals << "ms/eval\n";
        cout << "CPUtime " << cpuInSec << endl;

        printf("Used Integrator %s at accuracy %g:\n", 
            integ.getMethodName(), integ.getAccuracyInUse());
        printf("# STEPS/ATTEMPTS = %d/%d\n", integ.getNumStepsTaken(), integ.getNumStepsAttempted());
        printf("# ERR TEST FAILS = %d\n", integ.getNumErrorTestFailures());
        printf("# REALIZE/PROJECT = %d/%d\n", integ.getNumRealizations(), integ.getNumProjections());
        printf("# EVENT STEPS/HANDLER CALLS = %d/%d\n", 
            nStepsWithEvent, mbs.getNumHandleEventCalls());
    }

    for (int i=0; i<NTries; ++i)
        cout << "nstates " << i << " " << states[i].size() << endl;

    // Instant replay.
    while(true) {
        printf("Hit ENTER for replay (%d states) ...", 
                stateSaver->getNumSavedStates());
        getchar();
        for (int i=0; i < stateSaver->getNumSavedStates(); ++i) {
            mbs.realize(stateSaver->getState(i), Stage::Velocity);
            viz.report(stateSaver->getState(i));
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
//                       MY UNILATERAL CONSTRAINT SET
//==============================================================================

//------------------------ SELECT ACTIVE CONSTRAINTS ---------------------------
void MyUnilateralConstraintSet::
selectActiveConstraints(State& state, Real vtol) const {

    // Find all the contacts and stiction elements that might be active based
    // on kinematics.
    MyElementSubset candidates;

    bool needRestart;
    do {
        //TODO: this (mis)use of accuracy needs to be revisited.
        findCandidateElements(state, vtol, candidates);

        // Evaluate accelerations and reaction forces and check if 
        // any of the active contacts are generating negative ("pulling") 
        // forces; if so, inactivate them.
        findActiveCandidates(state, candidates);

        // Project active constraints to the constraint manifolds.
        const Real tol = vtol/1000;
        SimTK_DEBUG1("projecting active constraints to tol=%g\n", tol);
        m_mbs.project(state, tol);

        // It is possible that the projection above took some marginally-sliding
        // friction and slowed it down enough to make it a stiction candidate.
        needRestart = false;
        for (unsigned i=0; i < candidates.m_sliding.size(); ++i) {
            const int which = candidates.m_sliding[i];
            const MyHybridVertexContactElementImpl& contact = 
                getContactElement(which);
            assert(!contact.isSticking(state));
            if (   contact.getActualSlipSpeed(state) <= vtol
                || contact.calcSlipSpeedWitness(state) <= 0) 
            {
                SimTK_DEBUG3("***RESTART** selectActiveConstraints() because "
                    "sliding velocity %d is now |v|=%g or witness=%g\n",
                    which, contact.getActualSlipSpeed(state),
                    contact.calcSlipSpeedWitness(state));
                needRestart = true;
                break;
            }
        }
    } while (needRestart);
}




//-------------------------- FIND ACTIVE CANDIDATES ---------------------------
// Given a list of candidate stiction constraints,
// determine which ones are active in the least squares solution for the
// constraint multipliers. Candidates are those constraints that meet all 
// kinematic conditions -- for stiction, sliding velocity less than tolerance,
// or sliding direction reversed in the last step. Also, any
// constraint that is currently active is a candidate, regardless of its
// kinematics (might have drifted but that can't disable it).
//
// This method should be called only from within an event handler. For sliding
// friction it will have reset the "previous slip direction" to the current
// slip or impending slip direction.
//
// Algorithm
// ---------
// We're given a set of stiction candidates. If any are inactive, activate them.
// -- at this point all aerr==0, some ferr might be < 0
//
// loop
// - Realize(Acceleration) with the current active set
// - Calculate ferr for active constraints, aerr for inactive
// - If all ferr>=0, aerr>=0 -> break loop
// - Check for aerr < 0 [tol?]. Shouldn't happen but if it does must turn on the
//     associated constraint for the worst violation, then -> continue loop
// - Find worst (most negative) offender:
//    stiction offense = mu_s*max(0, fc) - |fs|
// - Inactivate chosen constraint
//     record impending slip direction stick->slide
// end loop 
//
void MyUnilateralConstraintSet::
findActiveCandidates(State& s, const MyElementSubset& candidates) const
{
    const int ReviseNormalNIters = 6;
    showConstraintStatus(s, "ENTER findActiveCandidates()");

    SimTK_DEBUG3(
        "findActiveCandidates() %d/%d/%d contact/stick/ineligible ...\n",
        candidates.m_contact.size(), candidates.m_stiction.size(),
        candidates.m_ineligible.size());

    // Disable any sticking constraints that are now ineligible due to
    // liftoff, and enable all other candidate stiction constraints if any 
    // are currently disabled.
    enableConstraintSubset(candidates, s);

    if (candidates.m_contact.empty()) {
        // Can't be any friction either, if there are no contacts.
        SimTK_DEBUG("EXIT findActiveCandidates: no candidates.\n");
        m_mbs.realize(s, Stage::Acceleration);
        return;
    }

    int pass=0, nStictionDisabled=0;
    while (true) {
        ++pass; 
        SimTK_DEBUG1("\nfindActiveCandidates(): pass %d\n", pass);

        // Given an active set, evaluate multipliers and accelerations.
        m_mbs.realize(s, Stage::Acceleration);
        if (pass==1) {
            // First time through record all the slip directions.
            for (unsigned i=0; i < candidates.m_contact.size(); ++i) {
                const int which = candidates.m_contact[i];
                const MyHybridVertexContactElementImpl& fric = 
                    getContactElement(which);
                if (fric.isSticking(s))
                    fric.recordImpendingSlipDir(s);
                else fric.recordActualSlipDir(s);
            }
        } 

        // Scan all candidate stictions to find the active one that has the
        // most negative capacity.

        int worstActiveStiction=-1; Real mostNegativeStictionCapacity=0;     
        SimTK_DEBUG("analyzing stiction constraints ...\n");
        for (unsigned i=0; i < candidates.m_stiction.size(); ++i) {
            const int which = candidates.m_stiction[i];
            SimTK_DEBUG1("  %d: ", which);
            const MyHybridVertexContactElementImpl& fric = 
                getContactElement(which);

            if (!fric.isSticking(s)) {
                SimTK_DEBUG("off\n");
                continue;
            }

            const Real capacity = fric.calcStictionForceWitness(s);
            SimTK_DEBUG2("on capacity=%g (N=%g)\n", 
                capacity, fric.getNormalForce(s));

            if (capacity < mostNegativeStictionCapacity) {
                worstActiveStiction = which;
                mostNegativeStictionCapacity = capacity;
            }
        }

        #ifndef NDEBUG
        if (mostNegativeStictionCapacity == 0)
            printf("  All active stiction constraints are satisfied.\n");
        else 
            printf("  Active stiction %d was worst violator with capacity=%g\n",
                worstActiveStiction, mostNegativeStictionCapacity);
        #endif

        if (mostNegativeStictionCapacity==0) {
            SimTK_DEBUG("DONE. Current active set is a winner.\n");
            break;
        }

        SimTK_DEBUG1("  Disable stiction %d\n", worstActiveStiction);
        const MyHybridVertexContactElementImpl& fric = 
            getContactElement(worstActiveStiction);

        ++nStictionDisabled;
        fric.disableStiction(s);

        // Go back for another pass.
    }

    // Reset all the slip directions so that all slip->stick event witnesses 
    // will be positive when integration resumes.
    for (unsigned i=0; i < candidates.m_contact.size(); ++i) {
        const int which = candidates.m_contact[i];
        const MyHybridVertexContactElementImpl& fric = 
            getContactElement(which);
        fric.updatePrevSlipDirFromRecorded(s);
    }

    // Always leave at acceleration stage.
    m_mbs.realize(s, Stage::Acceleration);

    SimTK_DEBUG1("... Done; %d stictions broken.\n", nStictionDisabled);

    showConstraintStatus(s, "EXIT findActiveCandidates()");
}



//-------------------------- SHOW CONSTRAINT STATUS ----------------------------
void MyUnilateralConstraintSet::
showConstraintStatus(const State& s, const String& place) const
{
#ifndef NDEBUG
    printf("\n%s: uni status @t=%.15g\n", place.c_str(), s.getTime());
    m_mbs.realize(s, Stage::Acceleration);
    for (int i=0; i < getNumContactElements(); ++i) {
        const MyHybridVertexContactElementImpl& contact = getContactElement(i);
        const bool isActive = contact.isInContact(s);
        printf("  %6s %2d %s, h=%g dh=%g f=%g\n", 
                isActive?"ACTIVE":"off", i, "hybrid", 
                contact.getHeight(s),contact.getHeightDot(s),
                isActive?contact.getNormalForce(s):Zero);
        if (!isActive) continue;

        const bool isSticking = contact.isSticking(s);
        printf("  %8s friction %2d, |v|=%g witness=%g\n", 
                isSticking?"STICKING":"sliding", i,
                contact.getActualSlipSpeed(s),
                isSticking?contact.calcStictionForceWitness(s)
                          :contact.calcSlipSpeedWitness(s));
        contact.writeFrictionInfo(s, "    ", std::cout);
    }
    printf("\n");
#endif
}

//------------------------ STRIBECK FRICTION STATICS ---------------------------
// This is extracted from Simbody's continuous friction model so that we can
// compare it with the new implementation.

// Input x goes from 0 to 1; output goes 0 to 1 but smoothed with an S-shaped 
// quintic with two zero derivatives at either end. Cost is 7 flops.
inline static Real step5(Real x) {
    assert(0 <= x && x <= 1);
    const Real x3=x*x*x;
    return x3*(10+x*(6*x-15)); //10x^3-15x^4+6x^5
}

// This is the sum of two curves:
// (1) a wet friction term mu_wet which is a linear function of velocity: 
//     mu_wet = uv*v
// (2) a dry friction term mu_dry which is a quintic spline with 4 segments:
//     mu_dry = 
//      (a) v=0..1: smooth interpolation from 0 to us
//      (b) v=1..3: smooth interp from us down to ud (Stribeck)
//      (c) v=3..Inf ud
// CAUTION: uv and v must be dimensionless in multiples of transition velocity.
// The function mu = mu_wet + mu_dry is zero at v=0 with 1st deriv (slope) uv
// and 2nd deriv (curvature) 0. At large velocities v>>0 the value is 
// ud+uv*v, again with slope uv and zero curvature. We want mu(v) to be c2
// continuous, so mu_wet(v) must have zero slope and curvature at v==0 and
// at v==3 where it takes on a constant value ud.
//
// Cost: stiction 12 flops
//       stribeck 14 flops
//       sliding 3 flops
// Curve looks like this:
//
//  us+uv     ***
//           *    *                     *
//           *     *               *____| slope = uv at Inf
//           *      *         *
// ud+3uv    *        *  *      
//          *          
//          *        
//         *
//  0  *____| slope = uv at 0
//
//     |    |           |
//   v=0    1           3 
//
// This calculates a composite coefficient of friction that you should use
// to scale the normal force to produce the friction force.
static Real stribeck(Real us, Real ud, Real uv, Real v) {
    const Real mu_wet = uv*v;
    Real mu_dry;
    if      (v >= 3) mu_dry = ud; // sliding
    else if (v >= 1) mu_dry = us - (us-ud)*step5((v-1)/2); // Stribeck
    else             mu_dry = us*step5(v); // 0 <= v < 1 (stiction)
    return mu_dry + mu_wet;
}