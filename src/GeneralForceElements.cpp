/* Portions copyright (c) 2006-7 Stanford University and Michael Sherman.
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
 * Private implementation of GeneralForceElements.
 */

#include "SimTKcommon.h"

#include "simbody/internal/common.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/GeneralForceElements.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"
#include "simbody/internal/MultibodySystem.h"

#include "ForceSubsystemRep.h"


namespace SimTK {


// 
// Define linear springs between two stations s1 and s2 of a matter subsystem
// (a station is a point fixed on a particular body). A spring has a stiffness k, 
// and a natural length x0 at which it generates no force. Define the separation
// vector v=s2-s1, with x=|v| the spring's current length.
//
// We allocate parameters in the state to hold all the properties of the
// force elements so that they can be changed and thus the subject of studies.
//
// Then the potential energy stored in the spring is 
//    pe = k(x-x0)^2/2
// Forces are generated on both points, as the negative gradient of the
// potential energy at that point: 
//    f1 = d pe/d s1 =  k(x-x0)v/x
//    f2 = d pe/d s2 = -k(x-x0)v/x.
// Note that force is undefined when x=0; we'll return NaN vectors in that case.

class GeneralForceElementsRep : public ForceSubsystemRep {

    // state entries
    struct TwoPointLinearSpringParameters {
        TwoPointLinearSpringParameters() { 
            station1.setToNaN(); station2.setToNaN();
            stiffness = naturalLength = CNT<Real>::getNaN();
        }

        TwoPointLinearSpringParameters(
            MobilizedBodyId b1, const Vec3& s1, 
            MobilizedBodyId b2, const Vec3& s2,
            const Real& k, const Real& x0) 
          : body1(b1), body2(b2), station1(s1), station2(s2), 
            stiffness(k), naturalLength(x0) 
        { 
            assert(b1.isValid() && b2.isValid() && b1 != b2);
            assert(stiffness >= 0 && naturalLength >= 0);
        }

        MobilizedBodyId body1, body2;
        Vec3   station1, station2;    // in body frames
        Real   stiffness, naturalLength;
    };

    struct TwoPointConstantForceParameters {
        TwoPointConstantForceParameters() { 
            station1.setToNaN(); station2.setToNaN();
            force = zeroEnergyDistance = CNT<Real>::getNaN();
        }

        TwoPointConstantForceParameters(
            MobilizedBodyId b1, const Vec3& s1, 
            MobilizedBodyId b2, const Vec3& s2,
            const Real& f, const Real& x0) 
          : body1(b1), body2(b2), station1(s1), station2(s2), 
            force(f), zeroEnergyDistance(x0) 
        { 
            assert(b1.isValid() && b2.isValid() && b1 != b2);
            assert(zeroEnergyDistance >= 0);
        }

        MobilizedBodyId body1, body2;
        Vec3   station1, station2;    // in body frames
        Real   force, zeroEnergyDistance;
    };

    struct TwoPointLinearDamperParameters {
        TwoPointLinearDamperParameters() : body1(-1), body2(-1)
        { 
            station1.setToNaN(); station2.setToNaN();
            damping = CNT<Real>::getNaN();
        }

        TwoPointLinearDamperParameters(
            MobilizedBodyId b1, const Vec3& s1, 
            MobilizedBodyId b2, const Vec3& s2,
            const Real& c) 
          : body1(b1), body2(b2), station1(s1), station2(s2), 
            damping(c)
        { 
            assert(b1.isValid() && b2.isValid() && b1 != b2);
            assert(damping >= 0);
        }

        MobilizedBodyId body1, body2;
        Vec3   station1, station2;    // in body frames
        Real   damping;
    };

    struct ConstantForceParameters {
        ConstantForceParameters() : body(-1) { 
            station_B.setToNaN(); force_G.setToNaN();
            fmag = CNT<Real>::getNaN();
        }

        ConstantForceParameters(
            MobilizedBodyId b, const Vec3& s_B, const Vec3& f_G)
          : body(b), station_B(s_B), force_G(f_G),
            fmag(f_G.norm())
        { 
            assert(b.isValid());
        }

        MobilizedBodyId body;
        Vec3   station_B;   // in body frame
        Vec3   force_G;   // in ground

        // Pre-calculated
        Real   fmag;       // force magnitude
    };

    struct ConstantTorqueParameters {
        ConstantTorqueParameters() : body(-1) { 
            torque_G.setToNaN();
            tmag = CNT<Real>::getNaN();
        }

        ConstantTorqueParameters(
            MobilizedBodyId b, const Vec3& t_G)
          : body(b), torque_G(t_G), tmag(t_G.norm())
        { 
            assert(b.isValid());
        }

        MobilizedBodyId body;
        Vec3   torque_G;   // in ground

        // Pre-calculated
        Real   tmag;       // torque magnitude
    };

    struct MobilityLinearSpringParameters {
        MobilityLinearSpringParameters() : body(-1), axis(-1) { 
            stiffness = naturalLength = CNT<Real>::getNaN();
        }
        MobilityLinearSpringParameters(MobilizedBodyId b, int a, const Real& k, const Real& q0)
          : body(b), axis(a), stiffness(k), naturalLength(q0)
        { 
            assert(b.isValid() && a >= 0 && stiffness >= 0. && naturalLength >= 0.);
        }

        MobilizedBodyId body;
        int    axis;
        Real   stiffness, naturalLength;
    };

    struct MobilityLinearDamperParameters {
        MobilityLinearDamperParameters() : body(-1), axis(-1) { 
            damping = CNT<Real>::getNaN();
        }
        MobilityLinearDamperParameters(MobilizedBodyId b, int a, const Real& c)
          : body(b), axis(a), damping(c)
        { 
            assert(b.isValid() && a >= 0 && damping >= 0.);
        }

        MobilizedBodyId body;
        int    axis;
        Real   damping;
    };

    // means force or torque, depending on the meaning of the generalized speed
    struct MobilityConstantForceParameters {
        MobilityConstantForceParameters() : axis(-1) { 
            force = CNT<Real>::getNaN();
        }
        MobilityConstantForceParameters(MobilizedBodyId b, int a, const Real& f)
          : body(b), axis(a), force(f)
        { 
            assert(b.isValid() && a >= 0);
        }

        MobilizedBodyId body; 
        int    axis;
        Real   force;
    };

    struct GlobalEnergyDrainParameters {
        GlobalEnergyDrainParameters() : damping(0) { }
        GlobalEnergyDrainParameters(const Real& c)
          : damping(c) 
        { 
            assert(damping >= 0);
        }
        Real damping;   // 0 means "none"
    };

    struct CustomForceParameters {
        CustomForceParameters() : uforce(0), calc(0), clone(0), destruct(0) { }

        // Note that we will make a private copy of the user's force element.
        CustomForceParameters(const GeneralForceElements::CustomForce&        u, 
                              GeneralForceElements::CustomForceCalcMethod     ucalc,
                              GeneralForceElements::CustomForceCloneMethod    uclone,
                              GeneralForceElements::CustomForceDestructor     udestruct)
          : uforce(uclone(u)), calc(ucalc), clone(uclone), destruct(udestruct)
        {
        }
        CustomForceParameters(const CustomForceParameters& src) {
            uforce   = src.clone(*src.uforce);
            calc     = src.calc;
            clone    = src.clone;
            destruct = src.destruct;
        }
        CustomForceParameters& operator=(const CustomForceParameters& src) {
            if (&src == this) return *this;
            if (uforce) {destruct(uforce);} // out with the old
            if (src.uforce) {
                uforce = src.clone(*src.uforce);
                calc = src.calc; clone = src.clone; destruct = src.destruct;
            } else {
                uforce=0; calc=0; clone=0; destruct=0;
            }
            return *this;
        }


        // This is to separate construction from filling in the parameters.
        // It is only allowed if the parameters are empty currently.
        // Note that we make a private copy of the user force element.
        void setCustomForceParameters(const GeneralForceElements::CustomForce&        u, 
                                      GeneralForceElements::CustomForceCalcMethod     ucalc,
                                      GeneralForceElements::CustomForceCloneMethod    uclone,
                                      GeneralForceElements::CustomForceDestructor     udestruct) 
        {
            assert(uforce==0);
            uforce=uclone(u); calc=ucalc; clone=uclone; destruct=udestruct;
        }

        ~CustomForceParameters() {
            destruct(uforce); // toss out our copy
            uforce=0; calc=0; clone=0; destruct=0;
        }
        GeneralForceElements::CustomForce*           uforce;  // we own this object
        GeneralForceElements::CustomForceCalcMethod  calc;
        GeneralForceElements::CustomForceCloneMethod clone;
        GeneralForceElements::CustomForceDestructor  destruct;
    };

    struct Parameters {
        Parameters() : enabled(true) { }
        bool enabled;
        std::vector<TwoPointLinearSpringParameters>  twoPointLinearSprings;
        std::vector<TwoPointLinearDamperParameters>  twoPointLinearDampers;
        std::vector<TwoPointConstantForceParameters> twoPointConstantForces;
        std::vector<ConstantForceParameters>         constantForces;
        std::vector<ConstantTorqueParameters>        constantTorques;
        std::vector<MobilityLinearSpringParameters>  mobilityLinearSprings;
        std::vector<MobilityLinearDamperParameters>  mobilityLinearDampers;
        std::vector<MobilityConstantForceParameters> mobilityConstantForces;
        std::vector<GlobalEnergyDrainParameters>     globalEnergyDrains;
        std::vector<CustomForceParameters>           customForces;
    };

    // topological variables
    Parameters defaultParameters;

    // This must be filled in during realizeTopology and treated
    // as const thereafter.
    mutable int instanceVarsIndex;

    const Parameters& getParameters(const State& s) const {
        assert(subsystemTopologyHasBeenRealized());
        return Value<Parameters>::downcast(
            getDiscreteVariable(s,instanceVarsIndex)).get();
    }
    Parameters& updParameters(State& s) const {
        assert(subsystemTopologyHasBeenRealized());
        return Value<Parameters>::downcast(
            updDiscreteVariable(s,instanceVarsIndex)).upd();
    }

public:
    GeneralForceElementsRep()
     : ForceSubsystemRep("GeneralForceElements", "0.0.1"), 
       instanceVarsIndex(-1)
    {
    }

    int addTwoPointLinearSpring(MobilizedBodyId body1, const Vec3& s1,
                                MobilizedBodyId body2, const Vec3& s2,
                                const Real& stiffness,
                                const Real& naturalLength)
    {
        assert(stiffness >= 0);
        assert(naturalLength >= 0);
        assert(body1 != body2);

        invalidateSubsystemTopologyCache();
        defaultParameters.twoPointLinearSprings.push_back(
            TwoPointLinearSpringParameters(body1,s1,body2,s2,stiffness,naturalLength));
        return (int)defaultParameters.twoPointLinearSprings.size() - 1;
    }

    int addTwoPointLinearDamper(MobilizedBodyId body1, const Vec3& s1,
                                MobilizedBodyId body2, const Vec3& s2,
                                const Real& damping)
    {
        assert(damping >= 0);
        assert(body1 != body2);

        invalidateSubsystemTopologyCache();
        defaultParameters.twoPointLinearDampers.push_back(
            TwoPointLinearDamperParameters(body1,s1,body2,s2,damping));
        return (int)defaultParameters.twoPointLinearDampers.size() - 1;
    }

    int addTwoPointConstantForce(MobilizedBodyId body1, const Vec3& s1,
                                 MobilizedBodyId body2, const Vec3& s2,
                                 const Real& force, const Real& zeroEnergyDistance)
    {
        assert(zeroEnergyDistance >= 0);
        assert(body1 != body2);

        invalidateSubsystemTopologyCache();
        defaultParameters.twoPointConstantForces.push_back(
            TwoPointConstantForceParameters(body1,s1,body2,s2,force,zeroEnergyDistance));
        return (int)defaultParameters.twoPointConstantForces.size() - 1;
    }


    int addConstantForce(MobilizedBodyId body, const Vec3& station_B, const Vec3& force_G)
    {
        invalidateSubsystemTopologyCache();
        defaultParameters.constantForces.push_back(
            ConstantForceParameters(body,station_B,force_G));
        return (int)defaultParameters.constantForces.size() - 1;
    }

    int addConstantTorque(MobilizedBodyId body, const Vec3& torque_G)
    {
        invalidateSubsystemTopologyCache();
        defaultParameters.constantTorques.push_back(
            ConstantTorqueParameters(body,torque_G));
        return (int)defaultParameters.constantTorques.size() - 1;
    }

    int addMobilityLinearSpring(MobilizedBodyId body, int axis,
                                const Real& stiffness,
                                const Real& naturalLength)
    {
        assert(stiffness >= 0);
        assert(naturalLength >= 0);

        invalidateSubsystemTopologyCache();
        defaultParameters.mobilityLinearSprings.push_back(
            MobilityLinearSpringParameters(body,axis,stiffness,naturalLength));
        return (int)defaultParameters.mobilityLinearSprings.size() - 1;
    }

    int addMobilityLinearDamper(MobilizedBodyId body, int axis,
                                const Real& damping)
    {
        assert(damping >= 0);

        invalidateSubsystemTopologyCache();
        defaultParameters.mobilityLinearDampers.push_back(
            MobilityLinearDamperParameters(body,axis,damping));
        return (int)defaultParameters.mobilityLinearDampers.size() - 1;
    }


    int addMobilityConstantForce(MobilizedBodyId body, int axis,
                                 const Real& force)
    {
        invalidateSubsystemTopologyCache();
        defaultParameters.mobilityConstantForces.push_back(
            MobilityConstantForceParameters(body,axis,force));
        return (int)defaultParameters.mobilityConstantForces.size() - 1;
    }

    int addGlobalEnergyDrain(const Real& dampingFactor) {
        assert(dampingFactor >= 0);

        invalidateSubsystemTopologyCache();
        defaultParameters.globalEnergyDrains.push_back(
            GlobalEnergyDrainParameters(dampingFactor));
        return (int)defaultParameters.globalEnergyDrains.size() - 1;
    }

    int addCustomForce(const GeneralForceElements::CustomForce& u, 
        GeneralForceElements::CustomForceCalcMethod             calc, 
        GeneralForceElements::CustomForceCloneMethod            clone, 
        GeneralForceElements::CustomForceDestructor             destruct) 
    {
        assert(calc && clone && destruct);

        invalidateSubsystemTopologyCache();
        defaultParameters.customForces.push_back(
            CustomForceParameters(u,calc,clone,destruct));
        return (int)defaultParameters.customForces.size() - 1;
    }

    // These override default implementations of virtual methods in the Subsystem::Guts
    // class.

    GeneralForceElementsRep* cloneImpl() const {return new GeneralForceElementsRep(*this);}

    int realizeSubsystemTopologyImpl(State& s) const {
        instanceVarsIndex = s.allocateDiscreteVariable(getMySubsystemId(), Stage::Instance, 
            new Value<Parameters>(defaultParameters));
        return 0;
    }

    int realizeSubsystemModelImpl(State& s) const {
        // Sorry, no choices available at the moment.
        return 0;
    }

    int realizeSubsystemInstanceImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemTimeImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemPositionImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemVelocityImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemDynamicsImpl(const State& s) const {
        const Parameters& p = getParameters(s);
        if (!p.enabled) return 0;

        const MultibodySystem&        mbs    = getMultibodySystem(); // my owner
        const SimbodyMatterSubsystem& matter = mbs.getMatterSubsystem();

        // Get access to system-global cache entries.
        Real&                  pe              = mbs.updPotentialEnergy(s, Stage::Dynamics);
        Vector_<SpatialVec>&   rigidBodyForces = mbs.updRigidBodyForces(s, Stage::Dynamics);
        Vector_<Vec3>&         particleForces  = mbs.updParticleForces (s, Stage::Dynamics);
        Vector&                mobilityForces  = mbs.updMobilityForces (s, Stage::Dynamics);

        // Linear mobility springs
        for (int i=0; i < (int)p.mobilityLinearSprings.size(); ++i) {
            const MobilityLinearSpringParameters& f = p.mobilityLinearSprings[i];
            const MobilizedBody& body = matter.getMobilizedBody(f.body);
            const Real q = body.getOneQ(s,f.axis);
            const Real frc = -f.stiffness*(q-f.naturalLength);
            pe -= 0.5*frc*(q-f.naturalLength);
            body.applyOneMobilityForce(s,f.axis,frc,mobilityForces);
        }

        // Linear mobility dampers
        for (int i=0; i < (int)p.mobilityLinearDampers.size(); ++i) {
            const MobilityLinearDamperParameters& f = p.mobilityLinearDampers[i];
            const MobilizedBody& body = matter.getMobilizedBody(f.body);
            const Real u = body.getOneU(s,f.axis);
            const Real frc = -f.damping*u;
            // no PE contribution
            body.applyOneMobilityForce(s,f.axis,frc,mobilityForces);
        }

        // Constant mobility forces
        for (int i=0; i < (int)p.mobilityConstantForces.size(); ++i) {
            const MobilityConstantForceParameters& f = p.mobilityConstantForces[i];
            const MobilizedBody& body = matter.getMobilizedBody(f.body);
            // no PE contribution
            body.applyOneMobilityForce(s,f.axis,f.force,mobilityForces);
        }

        // Global energy drain (no PE contribution)
        for (int i=0; i < (int)p.globalEnergyDrains.size(); ++i) {
            const Real c = p.globalEnergyDrains[i].damping;
            mobilityForces -= c*matter.getU(s);
        }

        // Two-point linear springs
        for (int i=0; i < (int)p.twoPointLinearSprings.size(); ++i) {
            const TwoPointLinearSpringParameters& spring =
                p.twoPointLinearSprings[i];
            const Transform& X_GB1 = matter.getMobilizedBody(spring.body1).getBodyTransform(s);
            const Transform& X_GB2 = matter.getMobilizedBody(spring.body2).getBodyTransform(s);

            const Vec3 s1_G = X_GB1.R() * spring.station1;
            const Vec3 s2_G = X_GB2.R() * spring.station2;

            const Vec3 p1_G = X_GB1.T() + s1_G; // station measured from ground origin
            const Vec3 p2_G = X_GB2.T() + s2_G;

            const Vec3 r_G = p2_G - p1_G; // vector from point1 to point2
            const Real d   = r_G.norm();  // distance between the points
            const Real stretch   = d - spring.naturalLength; // + -> tension, - -> compression
            const Real frcScalar = spring.stiffness*stretch; // k(x-x0)

            pe += 0.5 * frcScalar * stretch; // 1/2 k (x-x0)^2

            const Vec3 f1_G = (frcScalar/d) * r_G;
            rigidBodyForces[spring.body1] +=  SpatialVec(s1_G % f1_G, f1_G);
            rigidBodyForces[spring.body2] -=  SpatialVec(s2_G % f1_G, f1_G);
        }

        // Two-point constant force (no PE contribution)
        for (int i=0; i < (int)p.twoPointConstantForces.size(); ++i) {
            const TwoPointConstantForceParameters& frc =
                p.twoPointConstantForces[i];
            const Transform& X_GB1 = matter.getMobilizedBody(frc.body1).getBodyTransform(s);
            const Transform& X_GB2 = matter.getMobilizedBody(frc.body2).getBodyTransform(s);

            const Vec3 s1_G = X_GB1.R() * frc.station1;
            const Vec3 s2_G = X_GB2.R() * frc.station2;

            const Vec3 p1_G = X_GB1.T() + s1_G; // station measured from ground origin
            const Vec3 p2_G = X_GB2.T() + s2_G;

            const Vec3 r_G = p2_G - p1_G; // vector from point1 to point2
            const Real x   = r_G.norm();  // distance between the points
            const UnitVec3 d(r_G/x, true);

            const Vec3 f2_G = frc.force * d;
            rigidBodyForces[frc.body1] -=  SpatialVec(s1_G % f2_G, f2_G);
            rigidBodyForces[frc.body2] +=  SpatialVec(s2_G % f2_G, f2_G);
        }

        // Two-point linear dampers (no PE contribution)
        for (int i=0; i < (int)p.twoPointLinearDampers.size(); ++i) {
            const TwoPointLinearDamperParameters& damper =
                p.twoPointLinearDampers[i];
            const Transform& X_GB1 = matter.getMobilizedBody(damper.body1).getBodyTransform(s);
            const Transform& X_GB2 = matter.getMobilizedBody(damper.body2).getBodyTransform(s);

            const Vec3 s1_G = X_GB1.R() * damper.station1;
            const Vec3 s2_G = X_GB2.R() * damper.station2;

            const Vec3 p1_G = X_GB1.T() + s1_G; // station measured from ground origin
            const Vec3 p2_G = X_GB2.T() + s2_G;

            const Vec3 v1_G = matter.getMobilizedBody(damper.body1).calcBodyFixedPointVelocityInGround(s, damper.station1);
            const Vec3 v2_G = matter.getMobilizedBody(damper.body2).calcBodyFixedPointVelocityInGround(s, damper.station2);
            const Vec3 vRel = v2_G - v1_G; // relative velocity

            const UnitVec3 d(p2_G - p1_G); // direction from point1 to point2
            const Real frc = damper.damping*dot(vRel,d); // c*v

            const Vec3 f1_G = frc*d;
            rigidBodyForces[damper.body1] +=  SpatialVec(s1_G % f1_G, f1_G);
            rigidBodyForces[damper.body2] -=  SpatialVec(s2_G % f1_G, f1_G);
        }

        // Constant forces (no PE contribution)
        for (int i=0; i < (int)p.constantForces.size(); ++i) {
            const ConstantForceParameters& f = p.constantForces[i];
            if (f.fmag == 0)
                continue;

            const Transform& X_GB = matter.getMobilizedBody(f.body).getBodyTransform(s);
            const Vec3 station_G = X_GB.R() * f.station_B;

            rigidBodyForces[f.body] += SpatialVec(station_G % f.force_G, f.force_G);
        }

        // Constant torques (no PE contribution)
        for (int i=0; i < (int)p.constantTorques.size(); ++i) {
            const ConstantTorqueParameters& t = p.constantTorques[i];
            if (t.tmag == 0)
                continue;

            // update only the angular component of the spatial force on this body
            rigidBodyForces[t.body][0] += t.torque_G;
        }

        // User forces
        for (int i=0; i < (int)p.customForces.size(); ++i) {
            const CustomForceParameters& u = p.customForces[i];
            u.calc(*u.uforce, matter, s, 
                   rigidBodyForces, particleForces, mobilityForces, pe);
        }

        return 0;
    }

    int realizeSubsystemAccelerationImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    int realizeSubsystemReportImpl(const State& s) const {
        // Nothing to compute here.
        return 0;
    }

    friend std::ostream& operator<<(std::ostream& o, 
                         const GeneralForceElementsRep::Parameters&); 
};
// Useless, but required by Value<T>.
std::ostream& operator<<(std::ostream& o, 
                         const GeneralForceElementsRep::Parameters&) 
{assert(false);return o;}

    //////////////////////////
    // GeneralForceElements //
    //////////////////////////


/*static*/ bool 
GeneralForceElements::isInstanceOf(const ForceSubsystem& s) {
    return GeneralForceElementsRep::isA(s.getRep());
}
/*static*/ const GeneralForceElements&
GeneralForceElements::downcast(const ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const GeneralForceElements&>(s);
}
/*static*/ GeneralForceElements&
GeneralForceElements::updDowncast(ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<GeneralForceElements&>(s);
}

const GeneralForceElementsRep& 
GeneralForceElements::getRep() const {
    return dynamic_cast<const GeneralForceElementsRep&>(ForceSubsystem::getRep());
}
GeneralForceElementsRep&       
GeneralForceElements::updRep() {
    return dynamic_cast<GeneralForceElementsRep&>(ForceSubsystem::updRep());
}

// Create Subsystem but don't associate it with any System. This isn't much use except
// for making std::vector's, which require a default constructor to be available.
GeneralForceElements::GeneralForceElements()
  : ForceSubsystem() 
{
    adoptSubsystemGuts(new GeneralForceElementsRep());
}

GeneralForceElements::GeneralForceElements(MultibodySystem& mbs)
  : ForceSubsystem() 
{
    adoptSubsystemGuts(new GeneralForceElementsRep());
    mbs.addForceSubsystem(*this); // steal ownership
}

int GeneralForceElements::addTwoPointLinearSpring
   (MobilizedBodyId body1, const Vec3& s1,
    MobilizedBodyId body2, const Vec3& s2,
    const Real& stiffness,
    const Real& naturalLength) 
{
    return updRep().addTwoPointLinearSpring(body1,s1,body2,s2,stiffness,naturalLength);
}

int GeneralForceElements::addTwoPointConstantForce
   (MobilizedBodyId body1, const Vec3& s1,
    MobilizedBodyId body2, const Vec3& s2,
    const Real& force,
    const Real& zeroEnergyDistance) 
{
    return updRep().addTwoPointConstantForce(body1,s1,body2,s2,force,zeroEnergyDistance);
}

int GeneralForceElements::addTwoPointLinearDamper
   (MobilizedBodyId body1, const Vec3& s1,
    MobilizedBodyId body2, const Vec3& s2,
    const Real& damping) 
{
    return updRep().addTwoPointLinearDamper(body1,s1,body2,s2,damping);
}

int GeneralForceElements::addConstantForce
   (MobilizedBodyId body, const Vec3& s_B, const Vec3& f_G) 
{
    return updRep().addConstantForce(body, s_B, f_G);
}

int GeneralForceElements::addConstantTorque
   (MobilizedBodyId body, const Vec3& t_G) 
{
    return updRep().addConstantTorque(body, t_G);
}

int GeneralForceElements::addMobilityConstantForce
   (MobilizedBodyId body, int axis, const Real& f) {
    return updRep().addMobilityConstantForce(body, axis, f);
}

int GeneralForceElements::addMobilityLinearSpring
    (MobilizedBodyId body, int axis, const Real& stiffness, const Real& neutralValue)
{
    return updRep().addMobilityLinearSpring(body, axis, stiffness, neutralValue);
}

int GeneralForceElements::addMobilityLinearDamper
    (MobilizedBodyId body, int axis, const Real& dampingFactor)
{
    return updRep().addMobilityLinearDamper(body, axis, dampingFactor);
}

int GeneralForceElements::addGlobalEnergyDrain(const Real& dampingFactor) {
    return updRep().addGlobalEnergyDrain(dampingFactor);
}

int GeneralForceElements::addCustomForceMethods(const CustomForce& u, 
    CustomForceCalcMethod calc, CustomForceCloneMethod clone, 
    CustomForceDestructor destruct)
{
    return updRep().addCustomForce(u,calc,clone,destruct);
}

} // namespace SimTK

