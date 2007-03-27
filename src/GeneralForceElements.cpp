/* Portions copyright (c) 2006 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */


/**@file
 *
 * Private implementation of GeneralForceElements.
 */

#include "SimTKsimbody.h"
#include "simbody/internal/ForceSubsystem.h"

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
        TwoPointLinearSpringParameters() : body1(-1), body2(-1)
        { 
            station1.setToNaN(); station2.setToNaN();
            stiffness = naturalLength = CNT<Real>::getNaN();
        }

        TwoPointLinearSpringParameters(
            int b1, const Vec3& s1, int b2, const Vec3& s2,
            const Real& k, const Real& x0) 
          : body1(b1), body2(b2), station1(s1), station2(s2), 
            stiffness(k), naturalLength(x0) 
        { 
            assert(b1 >= 0 && b2 >= 0 && b1 != b2);
            assert(stiffness >= 0 && naturalLength >= 0);
        }

        BodyId body1, body2;
        Vec3   station1, station2;    // in body frames
        Real   stiffness, naturalLength;
    };

    struct TwoPointConstantForceParameters {
        TwoPointConstantForceParameters() : body1(-1), body2(-1)
        { 
            station1.setToNaN(); station2.setToNaN();
            force = zeroEnergyDistance = CNT<Real>::getNaN();
        }

        TwoPointConstantForceParameters(
            int b1, const Vec3& s1, int b2, const Vec3& s2,
            const Real& f, const Real& x0) 
          : body1(b1), body2(b2), station1(s1), station2(s2), 
            force(f), zeroEnergyDistance(x0) 
        { 
            assert(b1 >= 0 && b2 >= 0 && b1 != b2);
            assert(zeroEnergyDistance >= 0);
        }

        BodyId body1, body2;
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
            int b1, const Vec3& s1, int b2, const Vec3& s2,
            const Real& c) 
          : body1(b1), body2(b2), station1(s1), station2(s2), 
            damping(c)
        { 
            assert(b1 >= 0 && b2 >= 0 && b1 != b2);
            assert(damping >= 0);
        }

        BodyId body1, body2;
        Vec3   station1, station2;    // in body frames
        Real   damping;
    };

    struct ConstantForceParameters {
        ConstantForceParameters() : body(-1) { 
            station_B.setToNaN(); force_G.setToNaN();
            zeroEnergyHeight = fmag = zeroEnergy = CNT<Real>::getNaN();
        }

        ConstantForceParameters(
            int b, const Vec3& s_B, const Vec3& f_G, const Real& z)
          : body(b), station_B(s_B), force_G(f_G), zeroEnergyHeight(z),
            fmag(f_G.norm()), zeroEnergy(fmag*zeroEnergyHeight)
        { 
            assert(b >= 0);
        }

        BodyId body;
        Vec3   station_B;   // in body frame
        Vec3   force_G;   // in ground
        Real   zeroEnergyHeight;

        // Pre-calculated
        Real   fmag;       // force magnitude
        Real   zeroEnergy; // fmag*zeroEnergyHeight (subtract this from PE)
    };

    struct MobilityLinearSpringParameters {
        MobilityLinearSpringParameters() : body(-1), axis(-1) { 
            stiffness = naturalLength = CNT<Real>::getNaN();
        }
        MobilityLinearSpringParameters(int b, int a, const Real& k, const Real& q0)
          : body(b), axis(a), stiffness(k), naturalLength(q0)
        { 
            assert(b >= 0 && a >= 0 && stiffness >= 0. && naturalLength >= 0.);
        }

        BodyId body;
        int    axis;
        Real   stiffness, naturalLength;
    };

    struct MobilityLinearDamperParameters {
        MobilityLinearDamperParameters() : body(-1), axis(-1) { 
            damping = CNT<Real>::getNaN();
        }
        MobilityLinearDamperParameters(int b, int a, const Real& c)
          : body(b), axis(a), damping(c)
        { 
            assert(b >= 0 && a >= 0 && damping >= 0.);
        }

        BodyId body;
        int    axis;
        Real   damping;
    };

    struct MobilityConstantForceParameters {
        MobilityConstantForceParameters() : body(-1), axis(-1) { 
            force = CNT<Real>::getNaN();
        }
        MobilityConstantForceParameters(int b, int a, const Real& f, const Real& z)
          : body(b), axis(a), force(f), zeroEnergyQ(z)
        { 
            assert(b >= 0 && a >= 0);
        }

        BodyId body; 
        int    axis;
        Real   force, zeroEnergyQ;
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

    struct UserForceParameters {
        UserForceParameters() : uforce(0), calc(0), clone(0), nuke(0) { }
        UserForceParameters(GeneralForceElements::UserForce*              u, 
                            GeneralForceElements::UserForceCalcMethod     ucalc,
                            GeneralForceElements::UserForceCloneMethod    uclone,
                            GeneralForceElements::UserForceDestructor     unuke)
          : uforce(u), calc(ucalc), clone(uclone), nuke(unuke)
        {
        }
        UserForceParameters(const UserForceParameters& src) {
            uforce = src.clone(src.uforce);
            calc   = src.calc;
            clone  = src.clone;
            nuke   = src.nuke;
        }
        UserForceParameters& operator=(const UserForceParameters& src) {
            if (&src == this) return *this;
            if (uforce) {nuke(uforce);} // out with the old
            if (src.uforce) {
                uforce = src.clone(uforce);
                calc = src.calc; clone = src.clone; nuke = src.nuke;
            } else {
                uforce=0; calc=0; clone=0; nuke=0;
            }
            return *this;
        }


        // This is to separate construction from filling in the parameters.
        // It is only allowed if the parameters are empty currently.
        void setUserForceParameters(GeneralForceElements::UserForce*              u, 
                                    GeneralForceElements::UserForceCalcMethod     ucalc,
                                    GeneralForceElements::UserForceCloneMethod    uclone,
                                    GeneralForceElements::UserForceDestructor     unuke) 
        {
            assert(uforce==0);
            uforce=u; calc=ucalc; clone=uclone; nuke=unuke;
        }

        ~UserForceParameters() {
            nuke(uforce);
            uforce=0; calc=0; clone=0; nuke=0;
        }
        GeneralForceElements::UserForce*           uforce;  // we own this object
        GeneralForceElements::UserForceCalcMethod  calc;
        GeneralForceElements::UserForceCloneMethod clone;
        GeneralForceElements::UserForceDestructor  nuke;
    };

    struct Parameters {
        Parameters() : enabled(true) { }
        bool enabled;
        std::vector<TwoPointLinearSpringParameters>  twoPointLinearSprings;
        std::vector<TwoPointLinearDamperParameters>  twoPointLinearDampers;
        std::vector<TwoPointConstantForceParameters> twoPointConstantForces;
        std::vector<ConstantForceParameters>         constantForces;
        std::vector<MobilityLinearSpringParameters>  mobilityLinearSprings;
        std::vector<MobilityLinearDamperParameters>  mobilityLinearDampers;
        std::vector<MobilityConstantForceParameters> mobilityConstantForces;
        std::vector<GlobalEnergyDrainParameters>     globalEnergyDrains;
        std::vector<UserForceParameters>             userForces;
    };

    // topological variables
    Parameters defaultParameters;

    // These must be filled in during realizeTopology and treated
    // as const thereafter. These are garbage unless built=true.
    mutable int instanceVarsIndex;
    mutable bool built;

    const Parameters& getParameters(const State& s) const {
        assert(built);
        return Value<Parameters>::downcast(
            getDiscreteVariable(s,instanceVarsIndex)).get();
    }
    Parameters& updParameters(State& s) const {
        assert(built);
        return Value<Parameters>::downcast(
            updDiscreteVariable(s,instanceVarsIndex)).upd();
    }

public:
    GeneralForceElementsRep()
     : ForceSubsystemRep("GeneralForceElements", "0.0.1"), 
       instanceVarsIndex(-1), built(false)
    {
    }

    int addTwoPointLinearSpring(int body1, const Vec3& s1,
                                int body2, const Vec3& s2,
                                const Real& stiffness,
                                const Real& naturalLength)
    {
        assert(stiffness >= 0);
        assert(naturalLength >= 0);
        assert(body1 != body2);
        defaultParameters.twoPointLinearSprings.push_back(
            TwoPointLinearSpringParameters(body1,s1,body2,s2,stiffness,naturalLength));
        return (int)defaultParameters.twoPointLinearSprings.size() - 1;
    }

    int addTwoPointLinearDamper(int body1, const Vec3& s1,
                                int body2, const Vec3& s2,
                                const Real& damping)
    {
        assert(damping >= 0);
        assert(body1 != body2);
        defaultParameters.twoPointLinearDampers.push_back(
            TwoPointLinearDamperParameters(body1,s1,body2,s2,damping));
        return (int)defaultParameters.twoPointLinearDampers.size() - 1;
    }

    int addTwoPointConstantForce(int body1, const Vec3& s1,
                                 int body2, const Vec3& s2,
                                 const Real& force, const Real& zeroEnergyDistance)
    {
        assert(zeroEnergyDistance >= 0);
        assert(body1 != body2);
        defaultParameters.twoPointConstantForces.push_back(
            TwoPointConstantForceParameters(body1,s1,body2,s2,force,zeroEnergyDistance));
        return (int)defaultParameters.twoPointConstantForces.size() - 1;
    }


    int addConstantForce(int body, const Vec3& station_B, const Vec3& force_G,
                         const Real& zeroEnergyHeight)
    {
        defaultParameters.constantForces.push_back(
            ConstantForceParameters(body,station_B,force_G,zeroEnergyHeight));
        return (int)defaultParameters.constantForces.size() - 1;
    }

    int addMobilityLinearSpring(int body, int axis,
                                const Real& stiffness,
                                const Real& naturalLength)
    {
        assert(stiffness >= 0);
        assert(naturalLength >= 0);
        defaultParameters.mobilityLinearSprings.push_back(
            MobilityLinearSpringParameters(body,axis,stiffness,naturalLength));
        return (int)defaultParameters.mobilityLinearSprings.size() - 1;
    }

    int addMobilityLinearDamper(int body, int axis,
                                const Real& damping)
    {
        assert(damping >= 0);
        defaultParameters.mobilityLinearDampers.push_back(
            MobilityLinearDamperParameters(body,axis,damping));
        return (int)defaultParameters.mobilityLinearDampers.size() - 1;
    }


    int addMobilityConstantForce(int body, int axis,
                                 const Real& force, const Real& zeroEnergyValue)
    {
        defaultParameters.mobilityConstantForces.push_back(
            MobilityConstantForceParameters(body,axis,force,zeroEnergyValue));
        return (int)defaultParameters.mobilityConstantForces.size() - 1;
    }

    int addGlobalEnergyDrain(const Real& dampingFactor) {
        assert(dampingFactor >= 0);
        defaultParameters.globalEnergyDrains.push_back(
            GlobalEnergyDrainParameters(dampingFactor));
        return (int)defaultParameters.globalEnergyDrains.size() - 1;
    }

    int addUserForce(GeneralForceElements::UserForce* u, 
        GeneralForceElements::UserForceCalcMethod calc, 
        GeneralForceElements::UserForceCloneMethod clone, 
        GeneralForceElements::UserForceDestructor nuke) 
    {
        assert(u && calc && clone && nuke);
        defaultParameters.userForces.push_back(
            UserForceParameters(u,calc,clone,nuke));
        return (int)defaultParameters.userForces.size() - 1;
    }

    void realizeTopology(State& s) const {
        instanceVarsIndex = s.allocateDiscreteVariable(getMySubsystemIndex(), Stage::Instance, 
            new Value<Parameters>(defaultParameters));
        built = true;
    }

    void realizeModel(State& s) const {
        // Sorry, no choices available at the moment.
    }

    void realizeInstance(const State& s) const {
        // Nothing to compute here.
    }

    void realizeTime(const State& s) const {
        // Nothing to compute here.
    }

    void realizePosition(const State& s) const {
        // Nothing to compute here.
    }

    void realizeVelocity(const State& s) const {
        // Nothing to compute here.
    }

    void realizeDynamics(const State& s) const {
        const Parameters& p = getParameters(s);
        if (!p.enabled) return;

        const MultibodySystem& mbs    = getMultibodySystem(); // my owner
        const MatterSubsystem& matter = mbs.getMatterSubsystem();

        // Get access to system-global cache entries.
        Real&                  pe              = mbs.updPotentialEnergy(s);
        Vector_<SpatialVec>&   rigidBodyForces = mbs.updRigidBodyForces(s);
        Vector_<Vec3>&         particleForces  = mbs.updParticleForces(s);
        Vector&                mobilityForces  = mbs.updMobilityForces(s);

        // Linear mobility springs
        for (int i=0; i < (int)p.mobilityLinearSprings.size(); ++i) {
            const MobilityLinearSpringParameters& f = p.mobilityLinearSprings[i];
            const Real q = matter.getMobilizerQ(s,f.body,f.axis);
            const Real frc = -f.stiffness*(q-f.naturalLength);
            pe -= 0.5*frc*(q-f.naturalLength);
            matter.addInMobilityForce(s,f.body,f.axis,frc,mobilityForces);
        }

        // Linear mobility dampers
        for (int i=0; i < (int)p.mobilityLinearDampers.size(); ++i) {
            const MobilityLinearDamperParameters& f = p.mobilityLinearDampers[i];
            const Real u = matter.getMobilizerU(s,f.body,f.axis);
            const Real frc = -f.damping*u;
            // no PE contribution
            matter.addInMobilityForce(s,f.body,f.axis,frc,mobilityForces);
        }

        // Constant mobility forces
        for (int i=0; i < (int)p.mobilityConstantForces.size(); ++i) {
            const MobilityConstantForceParameters& f = p.mobilityConstantForces[i];
            const Real q = matter.getMobilizerQ(s,f.body,f.axis);
            pe -= f.force*(q-f.zeroEnergyQ);
            matter.addInMobilityForce(s,f.body,f.axis,f.force,mobilityForces);
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
            const Transform& X_GB1 = matter.getBodyTransform(s, spring.body1);
            const Transform& X_GB2 = matter.getBodyTransform(s, spring.body2);

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

        // Two-point constant force
        for (int i=0; i < (int)p.twoPointConstantForces.size(); ++i) {
            const TwoPointConstantForceParameters& frc =
                p.twoPointConstantForces[i];
            const Transform& X_GB1 = matter.getBodyTransform(s, frc.body1);
            const Transform& X_GB2 = matter.getBodyTransform(s, frc.body2);

            const Vec3 s1_G = X_GB1.R() * frc.station1;
            const Vec3 s2_G = X_GB2.R() * frc.station2;

            const Vec3 p1_G = X_GB1.T() + s1_G; // station measured from ground origin
            const Vec3 p2_G = X_GB2.T() + s2_G;

            const Vec3 r_G = p2_G - p1_G; // vector from point1 to point2
            const Real x   = r_G.norm();  // distance between the points
            const UnitVec3 d(r_G/x, true);

            pe -= frc.force * (x - frc.zeroEnergyDistance);

            const Vec3 f2_G = frc.force * d;
            rigidBodyForces[frc.body1] -=  SpatialVec(s1_G % f2_G, f2_G);
            rigidBodyForces[frc.body2] +=  SpatialVec(s2_G % f2_G, f2_G);
        }

        // Two-point linear dampers (no PE contribution)
        for (int i=0; i < (int)p.twoPointLinearDampers.size(); ++i) {
            const TwoPointLinearDamperParameters& damper =
                p.twoPointLinearDampers[i];
            const Transform& X_GB1 = matter.getBodyTransform(s, damper.body1);
            const Transform& X_GB2 = matter.getBodyTransform(s, damper.body2);

            const Vec3 s1_G = X_GB1.R() * damper.station1;
            const Vec3 s2_G = X_GB2.R() * damper.station2;

            const Vec3 p1_G = X_GB1.T() + s1_G; // station measured from ground origin
            const Vec3 p2_G = X_GB2.T() + s2_G;

            const Vec3 v1_G = matter.calcStationVelocity(s, damper.body1, damper.station1);
            const Vec3 v2_G = matter.calcStationVelocity(s, damper.body2, damper.station2);
            const Vec3 vRel = v2_G - v1_G; // relative velocity

            const UnitVec3 d(p2_G - p1_G); // direction from point1 to point2
            const Real frc = damper.damping*dot(vRel,d); // c*v

            const Vec3 f1_G = frc*d;
            rigidBodyForces[damper.body1] +=  SpatialVec(s1_G % f1_G, f1_G);
            rigidBodyForces[damper.body2] -=  SpatialVec(s2_G % f1_G, f1_G);
        }

        // Constant forces
        for (int i=0; i < (int)p.constantForces.size(); ++i) {
            const ConstantForceParameters& f = p.constantForces[i];
            if (f.fmag == 0)
                continue;

            const Transform& X_GB = matter.getBodyTransform(s, f.body);
            const Vec3 station_G = X_GB.R() * f.station_B;
            const Vec3 point_G   = X_GB.T() + station_G;

            const Real rawPE = -dot(point_G, f.force_G);
            pe += rawPE-f.zeroEnergy;

            rigidBodyForces[f.body] += SpatialVec(station_G % f.force_G, f.force_G);
        }

        // User forces
        for (int i=0; i < (int)p.userForces.size(); ++i) {
            const UserForceParameters& u = p.userForces[i];
            u.calc(u.uforce, matter, s, 
                   rigidBodyForces, particleForces, mobilityForces, pe);
        }
    }

    void realizeAcceleration(const State& s) const {
        // Nothing to compute here.
    }

    GeneralForceElementsRep* cloneSubsystemRep() const {return new GeneralForceElementsRep(*this);}
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
    return dynamic_cast<const GeneralForceElementsRep&>(*rep);
}
GeneralForceElementsRep&       
GeneralForceElements::updRep() {
    return dynamic_cast<GeneralForceElementsRep&>(*rep);
}

GeneralForceElements::GeneralForceElements() {
    rep = new GeneralForceElementsRep();
    rep->setMyHandle(*this);
}

int GeneralForceElements::addTwoPointLinearSpring
   (int body1, const Vec3& s1,
    int body2, const Vec3& s2,
    const Real& stiffness,
    const Real& naturalLength) 
{
    return updRep().addTwoPointLinearSpring(body1,s1,body2,s2,stiffness,naturalLength);
}

int GeneralForceElements::addTwoPointConstantForce
   (int body1, const Vec3& s1,
    int body2, const Vec3& s2,
    const Real& force,
    const Real& zeroEnergyDistance) 
{
    return updRep().addTwoPointConstantForce(body1,s1,body2,s2,force,zeroEnergyDistance);
}

int GeneralForceElements::addTwoPointLinearDamper
   (int body1, const Vec3& s1,
    int body2, const Vec3& s2,
    const Real& damping) 
{
    return updRep().addTwoPointLinearDamper(body1,s1,body2,s2,damping);
}

int GeneralForceElements::addConstantForce
   (int body, const Vec3& s_B, const Vec3& f_G, const Real& zeroEnergyHeight) 
{
    return updRep().addConstantForce(body, s_B, f_G, zeroEnergyHeight);
}

int GeneralForceElements::addMobilityConstantForce
   (int body, int axis, const Real& f, const Real& zeroEnergyValue) {
    return updRep().addMobilityConstantForce(body, axis, f, zeroEnergyValue);
}

int GeneralForceElements::addMobilityLinearSpring
    (int body, int axis, const Real& stiffness, const Real& neutralValue)
{
    return updRep().addMobilityLinearSpring(body, axis, stiffness, neutralValue);
}

int GeneralForceElements::addMobilityLinearDamper
    (int body, int axis, const Real& dampingFactor)
{
    return updRep().addMobilityLinearDamper(body, axis, dampingFactor);
}

int GeneralForceElements::addGlobalEnergyDrain(const Real& dampingFactor) {
    return updRep().addGlobalEnergyDrain(dampingFactor);
}

int GeneralForceElements::addUserForceMethods(UserForce* u, 
    UserForceCalcMethod calc, UserForceCloneMethod clone, 
    UserForceDestructor nuke)
{
    return updRep().addUserForce(u,calc,clone,nuke);
}

} // namespace SimTK

