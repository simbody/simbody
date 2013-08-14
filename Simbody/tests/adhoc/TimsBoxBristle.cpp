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

/* Investigate a compliant contact/bristle friction model, using Tim's
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

// Control whether we use a dwell-time model. Otherwise we allow instantaneous
// sliding/stiction transition.
#define ENABLE_DWELL_TIME

// Define this to run the simulation NTries times, saving the states and
// comparing them bitwise to see if the simulations are perfectly repeatable
// as they should be. You should see nothing but exact zeroes print out for
// second and subsequent runs.
//#define TEST_REPEATABILITY
static const int NTries=3;

//==============================================================================
//                     MY HYBRID VERTEX CONTACT ELEMENT
//==============================================================================
/* Given a contact material (with compliant material properties), this 
represents a compliant point-and-rigid-halfspace contact element. We track
the height h of a vertex V on body B over a halfspace H on body P. The 
halfspace frame H serves as the contact frame, with Hz in the contact 
normal direction and Hxy spanning the tangent plane. The pose X_PH is 
a given constant. 

Normal force
------------
If the vertex height h is above the halfspace surface (h>=0), no forces are 
generated. If below the surface (h<0) we generate a normal force of magnitude 
     N=max(0, -k*h(1-d*hdot)), with N >= 0
applied to both bodies at the contact point C, which is the point at h=0 just
up the halfspace normal direction from the vertex (because we're considering
the halfspace to be rigid). Denote by Cb the station (material point) of B 
that is coincident with C, and Cp the station of P that is coincident with C.

Bristle friction force
----------------------
Ref: Gonthier, et al. 2004 Multibody Sys. Dyn. 11:209-33.

k_br     bristle stiffness     (= alpha0 from Gonthier)
t_br     bristle time constant (= alpha1/alpha0 from Gonthier)
t_dw     dwell time constant
v_trans  Stribeck transition velocity
v_tol    numerical velocity noise limit (< v_trans/10)

f        tangential friction force vector (2d)

z        bristle deformation vector (2d state)
v        tangential slip velocity vector (2d) [v=0 when not in contact]
r        rolling condition [0,1] (0=slipping, 1=rolling) (=s from Gonthier)
r_dw     r delayed by time constants t_dw,t_br (state) [0,1]

f = { sat(f_br,f_max) - N*mu_v*v,     if in contact
    { 0,                              otherwise

f_br = -k_br*(z + t_br*zdot)    [TODO: need area term here?]
zdot = r * zdot_s + (1-r) * zdot_d
zdot_s = v
zdot_d = -(z + f_d/k_br)/t_br
f_d = -mu_d*N*dir(v, v_tol)

f_max = mu_eff*N
mu_eff = mu_d + (mu_s-mu_d)*r_dw

r = exp(-|v|^2/v_trans^2)

r_dw_dot ={ (r - r_dw)/t_dw,   r >= r_dw
          { (r - r_dw)/t_br,   r <  r_dw

sat(f_br, f_max) = { f_br,                 |f_br| <= f_max
                   { f_max*f_br/|f_br|,    |f_br| > f_max

dir(v, v_tol) = { v/|v|,   |v| >= v_tol
                { v/v_tol, |v| < v_tol [Gonthier has *stepUp(|v|/v_tol)]

*/

// Return f, but with magnitude limited by fmax>=0.
static inline Vec2 sat(const Vec2& f, Real fmax) {
    assert(fmax >= 0);
    const Real fsq = f.normSqr();
    if (fsq <= fmax*fmax) 
        return f;

    return (fmax/std::sqrt(fsq))*f;
}

// Implement relaxed unit vector in direction of v, where the direction vector
// has magnitude less than 1 if |v|<vtol, dropping to 0 when |v|=0.
static inline Vec2 dir(const Vec2& v, Real vnorm, Real oovtol) {
    assert(oovtol > 0);
    const Real scale = vnorm*oovtol;
    if (scale >= 1)
        return v/vnorm;
    //const Real smooth = 1;
    //const Real smooth = stepUp(scale);  
    const Real smooth = 1.5*scale-0.5*cube(scale);// Gonthier
    return v*(smooth*oovtol);
}

// This is a velocity-stage cache entry.
struct MyBristleContactInfo {
    MyBristleContactInfo()
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

class MyBristleVertexContactElementImpl : public Force::Custom::Implementation {
public:
    MyBristleVertexContactElementImpl(const GeneralForceSubsystem& forces,
        MobilizedBody& hsmobod, const UnitVec3& hsn, Real hsh,
        MobilizedBody& vmobod, const Vec3& vertex,
        const ContactMaterial& material, Real vTrans)
    :   m_matter(hsmobod.getMatterSubsystem()), m_forces(forces), 
        m_hsmobodx(hsmobod), m_X_PH(Rotation(hsn, ZAxis), hsh*hsn), 
        m_vmobodx(vmobod), m_vertex_B(vertex),
        m_material(material),
        m_kBristle(NaN), m_tBristle(NaN), m_tDwell(NaN),
        m_vtrans(vTrans), m_vtol(NaN), // assign later
        m_index(-1)
    {
        m_kBristle = m_material.getStiffness()/10;
        m_tBristle = .01; //  10 ms
        m_tDwell   = .1;  // 100 ms
        m_vtol = m_vtrans/10;
    }

    void setContactIndex(int index) {m_index=index;}

    void initialize(State& s) {
        updZ(s) = Vec2(0);
        updRDwell(s) = 0;
        setPrevSlipDir(s, Vec2(NaN));
    }

    void initializeForStiction(State& s) {
        updZ(s) = Vec2(0);
        updRDwell(s) = 1; // rolling
        setPrevSlipDir(s, Vec2(NaN));
    }

    bool isInContact(const State& s) const
    {   return findH(s) < 0; }

    bool isSticking(const State& s) const
    { 
        return getActualSlipSpeed(s) <= m_vtrans;
    }
    
    // Return a point in Ground coincident with the vertex.
    Vec3 whereToDisplay(const State& s) const {
        return getBodyB().findStationLocationInGround(s, m_vertex_B);
    }

    Real getActualSlipSpeed(const State& s) const {
        const MyBristleContactInfo& info = getContactInfo(s);
        return info.vSlipMag;
    }

    // Return the normal force N >= 0 currently being generated by this 
    // contact. State must be realized to Stage::Velocity.
    Real getNormalForce(const State& s) const {
        const MyBristleContactInfo& info = getContactInfo(s);
        if (info.h >= 0) return 0;
        const Real N = info.f_HCb[2]; // z component is normal
        return N;
    }

    Real getHeight(const State& s) const {
        const MyBristleContactInfo& info = getContactInfo(s);
        return info.h;
    }

    Real getHeightDot(const State& s) const {
        const MyBristleContactInfo& info = getContactInfo(s);
        return info.v_HCb[2]; // velocity in plane normal direction
    }

    // Return the sliding force 2-vector in the halfspace plane that is being
    // applied to the contact point station on body B. If we are sticking
    // then this force is zero; call getStickingForce() to get the real value.
    const Vec2& getSlidingForce(const State& s) const {
        const MyBristleContactInfo& info = getContactInfo(s);
        const Vec2& f = info.f_HCb.getSubVec<2>(0);
        return f;
    }

    // Return the sliding velocity 2-vector in the halfspace plane of the
    // contact point station on body B relative to the contact point station
    // on body P.
    const Vec2& getSlipVelocity(const State& s) const {
        const MyBristleContactInfo& info = getContactInfo(s);
        const Vec2& v = info.v_HCb.getSubVec<2>(0);
        return v;
    }

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

    const Vec2& getZ(const State& s) const 
    {   return Vec2::getAs(&m_forces.getZ(s)[this->m_bristleDeformationIx]); }
    Vec2& updZ(State& s) const 
    {   return Vec2::updAs(&m_forces.updZ(s)[this->m_bristleDeformationIx]); }
    Vec2& updZDot(const State& s) const 
    {   return Vec2::updAs(&m_forces.updZDot(s)[m_bristleDeformationIx]); }

    const Real& getRDwell(const State& s) const 
    {   return m_forces.getZ(s)[this->m_dwellIx]; }
    Real& updRDwell(State& s) 
    {   return m_forces.updZ(s)[this->m_dwellIx]; }
    Real& updRDwellDot(const State& s) const 
    {   return m_forces.updZDot(s)[m_dwellIx]; }


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


    const MyBristleContactInfo& getContactInfo(const State& state) const {
        const MyBristleContactInfo& info =
            Value<MyBristleContactInfo>::downcast
                (m_forces.getCacheEntry(state, m_contactInfoIx));
        return info;
    }

    MyBristleContactInfo& updContactInfo(const State& state) const {
        MyBristleContactInfo& info =
            Value<MyBristleContactInfo>::updDowncast
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
        const MyBristleContactInfo& info = getContactInfo(state);
        if (info.h >= 0) 
            return; // no contact

        const Vec3 f_GCb = info.R_GH * info.f_HCb; // re-express in Ground
        getBodyB().applyForceToBodyPoint(state, info.Cb,  f_GCb, bodyForces);
        getBodyP().applyForceToBodyPoint(state, info.Cp, -f_GCb, bodyForces);
    }

    // The normal force stores energy as 2/5 k h^(5/2) when h<0.
    Real calcPotentialEnergy(const State& state) const OVERRIDE_11 {
        const MyBristleContactInfo& info = getContactInfo(state);
        if (info.h >= 0) 
            return 0; // no contact
        const Real h52 = square(info.h)*std::sqrt(-info.h);
        const Real k = m_material.getStiffness();
        return 0.4*k*h52;
    }

    // Allocate state variable for storing the previous sliding direction.
    void realizeTopology(State& state) const OVERRIDE_11 {
        m_bristleDeformationIx = m_forces.allocateZ(state, Vector(2, Real(0)));
        m_dwellIx = m_forces.allocateZ(state, Vector(1, Real(0))); 
        m_prevSlipDirIx = m_forces.allocateAutoUpdateDiscreteVariable
           (state, Stage::Velocity, new Value<Vec2>(Vec2(NaN)), 
            Stage::Velocity);

        m_contactInfoIx = m_forces.allocateCacheEntry
           (state, Stage::Velocity, new Value<MyBristleContactInfo>());
    }


    // Calculate everything here and save in contact info cache entry where
    // it can be retrieved for generating forces, reporting, etc.
    void realizeVelocity(const State& state) const OVERRIDE_11 {
        MyBristleContactInfo& info = updContactInfo(state);

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

        const Vec2& z = getZ(state);
        const Real& rdw = getRDwell(state);
        Vec2& zdot   = updZDot(state);      // to be computed
        Real& rdwdot = updRDwellDot(state); 

        info.f_HCb = Vec3(0); // Assume no contact force.
        Real N = 0;

        if (info.h < 0) { // contact
            const Real h32 = -info.h*sqrt(-info.h); // |h| ^ (3/2)
            N = std::max(Real(0), k*h32*(1-c*hdot));
        }

        if (N==0) {
            zdot   = -z / m_tBristle; // decay fast
            rdwdot = -rdw / m_tBristle;
            return; // no contact force
        }

        // N is the Hz component of the force on Cb.

        Vec2 fT_HCb(0); // This will be the (Hx,Hy) components of force on Cb.

        // Make v unitless slip ratio.
        const Real v = info.vSlipMag / m_vtrans;
        Real r = std::exp(-square(v)); // 1==rolling, 0==slipping

        const Vec2& dir_prev = getPrevSlipDir(state);
        Vec2 dir_eff = dir(vSlip_HC, info.vSlipMag, 1/m_vtol);
        if (~dir_eff*dir_prev < 0) {
            //cout << "REVERSED from " << dir_prev << " to " << dir_eff << endl;
            //dir_eff = dir_prev;
        }
        updPrevSlipDirAutoUpdateValue(state) = dir_eff;
        markPrevSlipDirAutoUpdateValueValid(state);
                                  
        // This is what the force would be if we were fully in sliding.
        //const Vec2 f_d = -mu_d*N*dir(vSlip_HC, info.vSlipMag, 1/m_vtol);
        const Vec2 f_d = -mu_d*N*dir_eff;
        // This is the bristle deformation that would generate f_d
        const Vec2 z_d = -f_d/m_kBristle;

        const Vec2 zdot_s = vSlip_HC;
        const Vec2 zdot_d = (z_d - z)/m_tBristle;

        #ifdef ENABLE_DWELL_TIME
        const Real r_dwell = clamp(0,rdw,1);
        // Update dwell state derivative.
        const Real tconst = r >= rdw ? m_tDwell : m_tBristle;
        rdwdot = (r - rdw) / tconst;
        #else
        const Real r_dwell = r;
        rdwdot = 0;
        #endif

        zdot = r*zdot_s + (1-r)*zdot_d; // Gonthier
        //zdot = r_dwell*zdot_s + (1-r_dwell)*zdot_d; // is this better?
        const Vec2 f_br = -m_kBristle*(z + m_tBristle*zdot);

        const Real mu_eff = mu_d + (mu_s-mu_d)*r_dwell;
        const Real f_max = mu_eff*N;
        fT_HCb = sat(f_br, f_max) - N * mu_v*vSlip_HC;

        info.f_HCb = Vec3(fT_HCb[0], fT_HCb[1], N); // force on Cb, in H
    }

    std::ostream& writeFrictionInfo(const State& s, const String& indent, 
                                    std::ostream& o) const 
    {
        o << indent;
        if (!isInContact(s)) o << "OFF";
        else if (isSticking(s)) o << "STICK";
        else o << "SLIP";

        const Vec2 v = getSlipVelocity(s);
        const Vec2 f = getSlidingForce(s);
        o << " V=" << ~v << " |v|=" << v.norm() << endl;
        o << indent << " F=" << ~f << " |f|=" << f.norm() << endl;
        o << indent << " z=" << ~getZ(s) << " |z|=" << getZ(s).norm() << endl;
        o << indent << " rdw=" << getRDwell(s) 
          << " prevDir=" << ~getPrevSlipDir(s) << endl;
        return o;
    }


    void showContactForces(const State& s, Array_<DecorativeGeometry>& geometry) 
        const
    {
        const Real Scale = 0.01;
        const Real NH = getNormalForce(s);

        const bool isInStiction = getActualSlipSpeed(s) <= m_vtrans;
        const Vec2 fH = getSlidingForce(s);

        if (fH.normSqr() < square(SignificantReal) && NH < SignificantReal)
            return;

        const MyBristleContactInfo& info = getContactInfo(s);
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
    Real                            m_kBristle; // bristle stiffness
    Real                            m_tBristle; // bristle time constant
    Real                            m_tDwell;   // dwell time constant

    Real                            m_vtrans;   // transition velocity for
                                                //   Stribeck stiction
    Real                            m_vtol;     // max meaningful slip velocity

    int                             m_index;    // unique id for this contact

    // Set during realizeTopology().
    mutable ZIndex                  m_bristleDeformationIx; // zx,zy
    mutable ZIndex                  m_dwellIx;              // r_dw
    mutable DiscreteVariableIndex   m_prevSlipDirIx; // previous slip direction

    mutable CacheEntryIndex         m_contactInfoIx; 
};


//==============================================================================
//                       MY UNILATERAL CONSTRAINT SET
//==============================================================================

class MyUnilateralConstraintSet {
public:
    // Transition velocity: if a slip velocity is smaller than this the
    // contact is a candidate for stiction.
    MyUnilateralConstraintSet(const MultibodySystem& mbs)
    :   m_mbs(mbs) {}

    // Ownership of this force element belongs to the System; we're just keeping
    // a reference to it here.
    int addBristleElement(MyBristleVertexContactElementImpl* vertex) {
        const int index = (int)m_bristle.size();
        m_bristle.push_back(vertex);
        vertex->setContactIndex(index);
        return index;
    }

    int getNumContactElements() const {return (int)m_bristle.size();}
    const MyBristleVertexContactElementImpl& getContactElement(int ix) const 
    {   return *m_bristle[ix]; }
    MyBristleVertexContactElementImpl& updContactElement(int ix) 
    {   return *m_bristle[ix]; }

    // Initialize all states and runtime fields in the contact elements.
    void initialize(State& state)
    {
        for (unsigned i=0; i < m_bristle.size(); ++i)
            m_bristle[i]->initialize(state);
    }

    // In Debug mode, produce a useful summary of the current state of the
    // contact and friction elements.
    void showContactStatus(const State& state, const String& place) const;

    ~MyUnilateralConstraintSet() {
        // Nothing to delete since we are only holding references.
    }

    const MultibodySystem& getMultibodySystem() const {return m_mbs;}
private:
    const MultibodySystem&                      m_mbs;
    Array_<MyBristleVertexContactElementImpl*>  m_bristle; // unowned ref
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
        m_unis.showContactStatus(s, "STATE SAVER");
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
            const MyBristleVertexContactElementImpl& contact = 
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
    const Real BrickMass = 10;
    #ifdef USE_TIMS_PARAMS
        const Real RunTime=16;  // Tim's time
        const Real Stiffness = 2e7;
        const Real Dissipation = 1;
        const Real CoefRest = 0; 
        const Real mu_d = .5; /* compliant: .7*/
        const Real mu_s = .8; /* compliant: .7*/
        const Real mu_v = /*0.05*/0; //TODO: fails with mu_v=1, vtrans=.01
        const Real TransitionVelocity = 0.01;
        const Inertia brickInertia(.1,.1,.1);
        //const Real Radius = .02;
        const Real Radius = 1;
    #else
        const Real RunTime=20;
        const Real Stiffness = 1e6;
        const Real CoefRest = 0; 
        const Real TargetVelocity = 3; // speed at which to match coef rest
//        const Real Dissipation = (1-CoefRest)/TargetVelocity;
        const Real Dissipation = .1;
        const Real mu_d = .5;
        const Real mu_s = .8;
        const Real mu_v = 0*0.05;
        const Real TransitionVelocity = 0.01;
        const Inertia brickInertia(BrickMass*UnitInertia::brick(BrickHalfDims));
        const Real Radius = BrickHalfDims[0]/3;
    #endif

    printf("\n******************** Tim's Box Bristle ********************\n");
    printf("USING BRISTLE MODEL: Compliant material/integrated stiction\n");
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
    printf("******************** Tim's Box Bristle ********************\n\n");

        // CREATE MULTIBODY SYSTEM AND ITS SUBSYSTEMS
    MultibodySystem             mbs;
    SimbodyMatterSubsystem      matter(mbs);
    GeneralForceSubsystem       forces(mbs);
    Force::Gravity              gravity(forces, matter, -YAxis, 9.81);
    //Force::Gravity              gravity(forces, matter, -UnitVec3(.3,1,0), 3*9.81);
    MobilizedBody& Ground = matter.updGround();

    // Define a material to use for contact. This is not very stiff.
    ContactMaterial material(std::sqrt(Radius)*Stiffness,
                             Dissipation,
                             mu_s,  // mu_static
                             mu_d,  // mu_dynamic
                             mu_v); // mu_viscous

        // ADD MOBILIZED BODIES AND CONTACT CONSTRAINTS

    MyUnilateralConstraintSet unis(mbs);

    Body::Rigid brickBody = 
        Body::Rigid(MassProperties(BrickMass, Vec3(0), brickInertia));

    MobilizedBody::Free brick(Ground, Vec3(0),
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
                                                    49.033 * Vec3(0,1,0),
                                                    9, 9+2));
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
        MyBristleVertexContactElementImpl* vertex =
            new MyBristleVertexContactElementImpl(forces,
                Ground, YAxis, 0, // halfplane
                brick, pt, material, TransitionVelocity);
        Force::Custom(forces, vertex); // add force element to system
        unis.addBristleElement(vertex); // assign index, transition velocity
    }

    matter.setShowDefaultGeometry(false);
    Visualizer viz(mbs);
    viz.setShowSimTime(true);
    viz.setShowFrameNumber(true);
    viz.setShowFrameRate(true);
    viz.addDecorationGenerator(new ShowContact(unis));

    #ifdef ANIMATE
    mbs.addEventReporter(new Visualizer::Reporter(viz, ReportInterval));
    #else
    // This does nothing but interrupt the simulation so that exact step
    // sequence will be maintained with animation off.
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
    //Real accuracy = 1e-4;
    Real accuracy = 1e-3;
    //Real accuracy = 1e-5;
    //ExplicitEulerIntegrator integ(mbs);
    //RungeKutta2Integrator integ(mbs);
    //RungeKutta3Integrator integ(mbs);
    //SemiExplicitEulerIntegrator integ(mbs, .005);
    SemiExplicitEuler2Integrator integ(mbs);
    //RungeKuttaFeldbergIntegrator integ(mbs);
    //RungeKuttaMersonIntegrator integ(mbs);
    //VerletIntegrator integ(mbs);
    //CPodesIntegrator integ(mbs);
    #endif

    integ.setAccuracy(accuracy);
    //integ.setMaximumStepSize(0.25);
    integ.setMaximumStepSize(0.05);
    //integ.setMaximumStepSize(0.002);
    //integ.setAllowInterpolation(false);

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
    brick.setQToFitTranslation(s, Vec3(0,10,0));
    brick.setUToFitLinearVelocity(s, Vec3(0,0,0));
    #else
    brick.setQToFitTranslation(s, Vec3(0,1.4,0));
    brick.setUToFitLinearVelocity(s, Vec3(10,0,0));
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
    for (int i=0; i < unis.getNumContactElements(); ++i) {
        MyBristleVertexContactElementImpl& fric = unis.updContactElement(i);
        if (!fric.isInContact(s)) continue;
        const Real vSlip = fric.getActualSlipSpeed(s);
        fric.initializeForStiction(s); // just in case
        printf("friction element %d has v_slip=%g%s\n", i, vSlip,
            vSlip==0?" (ENABLING STICTION)":"");
    }

    // Make sure Lapack gets initialized.
    Matrix M(1,1); M(0,0)=1.23;
    FactorLU Mlu(M);

    
    // Simulate it.

    integ.setReturnEveryInternalStep(true);
    TimeStepper ts(mbs, integ);
    ts.setReportAllSignificantStates(true);

    #ifdef TEST_REPEATABILITY
        const int tries = NTries;
    #else
        const int tries = 1;
    #endif

    Array_< Array_<State> > states(tries);
    Array_< Array_<Real> > timeDiff(tries-1);
    Array_< Array_<Vector> > yDiff(tries-1);
    cout.precision(18);
    for (int ntry=0; ntry < tries; ++ntry) {
        mbs.resetAllCountersToZero();
        unis.initialize(s); // reinitialize
        ts.updIntegrator().resetAllStatistics();
        ts.initialize(s);
        int nStepsWithEvent = 0;

        const double startReal = realTime();
        const double startCPU = cpuTime();

        Integrator::SuccessfulStepStatus status;
        do {
            status=ts.stepTo(RunTime);
            #ifdef TEST_REPEATABILITY
                states[ntry].push_back(ts.getState());
            #endif          
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
                SimTK_DEBUG3("step  %3d @%.17g status=%s\n",
                    integ.getNumStepsTaken(), ts.getTime(),
                    Integrator::getSuccessfulStepStatusString(status).c_str());
            }
            #ifndef NDEBUG
                viz.report(ts.getState());
            #endif
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

    for (int i=0; i<tries; ++i)
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


//--------------------------- SHOW CONTACT STATUS ------------------------------
void MyUnilateralConstraintSet::
showContactStatus(const State& s, const String& place) const
{
#ifndef NDEBUG
    printf("\n%s: uni status @t=%.15g\n", place.c_str(), s.getTime());
    m_mbs.realize(s, Stage::Acceleration);
    for (int i=0; i < getNumContactElements(); ++i) {
        const MyBristleVertexContactElementImpl& contact = getContactElement(i);
        const bool isActive = contact.isInContact(s);
        printf("  %6s %2d %s, h=%g dh=%g f=%g\n", 
                isActive?"ACTIVE":"off", i, "bristle", 
                contact.getHeight(s),contact.getHeightDot(s),
                isActive?contact.getNormalForce(s):Zero);
        if (!isActive) continue;

        const bool isSticking = contact.isSticking(s);
        printf("  %8s friction %2d, |v|=%g\n", 
                isSticking?"STICKING":"sliding", i,
                contact.getActualSlipSpeed(s));
        contact.writeFrictionInfo(s, "    ", std::cout);
    }
    printf("\n");
#endif
}
