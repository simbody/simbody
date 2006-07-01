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
 * Private implementation of HuntCrossleyContact Subsystem.
 */

#include "Simbody.h"
#include "simbody/internal/ForceSubsystem.h"
#include "simbody/internal/HuntCrossleyContact.h"

#include "ForceSubsystemRep.h"


namespace SimTK {

class HuntCrossleyContactRep : public ForceSubsystemRep {

    // state entries
    struct SphereParameters {
        SphereParameters()
          : body(-1)
        { 
            center.setToNaN();
            radius = stiffness = dissipation = CNT<Real>::getNaN();
        }

        SphereParameters(
            int b, const Vec3& ctr,
            const Real& r, const Real& k, const Real& c) 
          : body(b), center(ctr), radius(r), stiffness(k), dissipation(c) 
        { 
            assert(b >= 0);
            assert(radius > 0 && stiffness >= 0 && dissipation >= 0);
        }

        int  body;
        Vec3 center;    // in body frame
        Real radius, stiffness, dissipation;    // r,k,c from H&C
    };

    struct HalfspaceParameters {
        HalfspaceParameters()
          : body(-1)
        { 
            height = stiffness = dissipation = CNT<Real>::getNaN();
        }

        HalfspaceParameters(
            int b, const UnitVec3& n,
            const Real& h, const Real& k, const Real& c) 
          : body(b), normal(n), height(h), stiffness(k), dissipation(c) 
        { 
            assert(b >= 0);
            assert(stiffness >= 0 && dissipation >= 0);
        }

        int  body;
        UnitVec3 normal;    // in body frame
        Real height, stiffness, dissipation;
    };


    struct Parameters {
        Parameters() : enabled(true) { }
        bool enabled;
        std::vector<SphereParameters>    spheres;
        std::vector<HalfspaceParameters> halfSpaces;
    };

    // topological variables
    Parameters defaultParameters;

    // These must be filled in during realizeConstruction and treated
    // as const thereafter. These are garbage unless built=true.
    mutable int parameterVarsIndex;
    mutable bool built;

    const Parameters& getParameters(const State& s) const {
        assert(built);
        return Value<Parameters>::downcast(
            getDiscreteVariable(s,parameterVarsIndex)).get();
    }
    Parameters& updParameters(State& s) const {
        assert(built);
        return Value<Parameters>::downcast(
            updDiscreteVariable(s,parameterVarsIndex)).upd();
    }

public:
    HuntCrossleyContactRep()
     : ForceSubsystemRep("HuntCrossleyContact", "0.0.1"), 
       parameterVarsIndex(-1), built(false)
    {
    }
    int addSphere(int body, const Vec3& center,
                  const Real& radius,
                  const Real& stiffness,
                  const Real& dissipation) 
    {
        assert(body >= 0 && radius > 0 && stiffness >= 0 && dissipation >= 0);
        defaultParameters.spheres.push_back(
            SphereParameters(body,center,radius,stiffness,dissipation));
        return (int)defaultParameters.spheres.size() - 1;    
    }

    int addHalfSpace(int body, const UnitVec3& normal,
                     const Real& height,
                     const Real& stiffness,
                     const Real& dissipation)
    {
        assert(body >= 0 && stiffness >= 0 && dissipation >= 0);
        defaultParameters.halfSpaces.push_back(
            HalfspaceParameters(body,normal,height,stiffness,dissipation));
        return (int)defaultParameters.halfSpaces.size() - 1;
    }

    void realizeConstruction(State& s) const {
        parameterVarsIndex = s.allocateDiscreteVariable(getMySubsystemIndex(), Stage::Parametrized, 
            new Value<Parameters>(defaultParameters));
        built = true;
    }

    void realizeModeling(State& s) const {
        // Sorry, no choices available at the moment.
    }

    void realizeParameters(const State& s) const {
        // Nothing to compute here.
    }

    void realizeTime(const State& s) const {
        // Nothing to compute here.
    }

    void realizeConfiguration(const State& s) const {
        // Nothing to compute here.
    }

    void realizeMotion(const State& s) const {
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

        // TODO: much of this should be precalculated

        for (int s1=0; s1 < (int)p.spheres.size(); ++s1) {
            const int body1 = p.spheres[s1].body;
            const Transform&  X_GB1 = matter.getBodyConfiguration(s,body1);
            const SpatialVec& V_GB1 = matter.getBodyVelocity(s,body1); // in G

            const Vec3& center1_B  = p.spheres[s1].center;
            const Vec3  center1_BG = X_GB1.R()*center1_B; // station, expr. in G
            const Vec3  center1_G  = X_GB1.T() +  center1_BG;

            const Real r1 = p.spheres[s1].radius;
            const Real k1 = p.spheres[s1].stiffness;
            const Real c1 = p.spheres[s1].dissipation;

            for (int s2=s1+1; s2 < (int)p.spheres.size(); ++s2) {
                const int body2 = p.spheres[s2].body;
                if (body2 == body1) continue;
                const Real           r2 = p.spheres[s2].radius;
                const Transform&  X_GB2 = matter.getBodyConfiguration(s,body2);
                const SpatialVec& V_GB2 = matter.getBodyVelocity(s,body2);

                const Vec3& center2_B  = p.spheres[s2].center;
                const Vec3  center2_BG = X_GB2.R()*center2_B; // station, expr. in G
                const Vec3  center2_G  = X_GB2.T() +  center2_BG;


                const Vec3 c2c1_G = center1_G - center2_G; // points at sphere1
                const Real dsq    = c2c1_G.normSqr();
                if (dsq >= (r1+r2)*(r1+r2)) continue;

                // There is a collision

                const Real k2 = p.spheres[s2].stiffness;
                const Real c2 = p.spheres[s2].dissipation;
                const Real k=(k1*k2)/(k1+k2);
                const Real c=(c1*k2 + c2*k1)/(k1+k2);

                const Real     d = std::sqrt(dsq);        // distance between the centers
                const UnitVec3 c2c1dir_G(c2c1_G/d, true); // direction from c2 center to c1 center

                const Real x = (r1+r2) - d; // penetration distance (total "squish"), > 0
                const Real squish1 = x*k2/(k1+k2);
                const Real squish2 = x - squish1; // i.e., x*k1/(k1+k2)
                const Vec3 contactPt_G = center2_G + (d-squish2)*c2c1dir_G;
                const Vec3 contactPt_B1_G = contactPt_G - X_GB1.T();   // meas from B1
                const Vec3 contactPt_B2_G = contactPt_G - X_GB2.T();   // meas from B2
                const Vec3 vContactB1_G = V_GB1[1] + V_GB1[0] % contactPt_B1_G ;
                const Vec3 vContactB2_G = V_GB2[1] + V_GB2[0] % contactPt_B2_G;
                const Real v = - ~(vContactB1_G-vContactB2_G)*c2c1dir_G; // dx/dt

                const Real x32 = x*std::sqrt(x); // x^(3/2)
                const Real f = k*x32*(1 + 1.5*c*v);
                const Vec3 fvec = f * c2c1dir_G; // points towards c1

                pe += 0.4*k*x32*x; // i.e., 2/5 k x^(5/2)
                rigidBodyForces[body1] += SpatialVec( contactPt_B1_G % fvec, fvec);
                rigidBodyForces[body2] -= SpatialVec( contactPt_B2_G % fvec, fvec);
            }

            // half spaces
            for (int h=0; h < (int)p.halfSpaces.size(); ++h) {
                const int body2 = p.halfSpaces[h].body;
                if (body2 == body1) continue;

                const Real        height = p.halfSpaces[h].height;
                const Transform&  X_GB2  = matter.getBodyConfiguration(s,body2);
                const SpatialVec& V_GB2  = matter.getBodyVelocity(s,body2);

                const Vec3& normal_B  = p.halfSpaces[h].normal;
                const Vec3  normal_G  = X_GB2.R()*normal_B;

                const Real h_G   = ~X_GB2.T()*normal_G + height; // height of half space over G
                const Real hc1_G = ~center1_G*normal_G; // height of s1's center over G
                const Real d = hc1_G - h_G; // height of center over half space
                if (d >= r1) continue;

                // There is a collision
   
                const Real k2 = p.halfSpaces[h].stiffness;
                const Real c2 = p.halfSpaces[h].dissipation;
                const Real k=(k1*k2)/(k1+k2);
                const Real c=(c1*k2 + c2*k1)/(k1+k2);

                const Real x = r1 - d; // penetration distance (total "squish"), > 0
                const Real squish1 = x*k2/(k1+k2);  // squish of sphere
                const Real squish2 = x - squish1;   // i.e., x*k1/(k1+k2) (squish of half space)
                const Vec3 contactPt_G = center1_G - (d-squish1)*normal_G;
                const Vec3 contactPt_B1_G = contactPt_G - X_GB1.T();   // meas from B1
                const Vec3 contactPt_B2_G = contactPt_G - X_GB2.T();   // meas from B2
                const Vec3 vContactB1_G = V_GB1[1] + V_GB1[0] % contactPt_B1_G ;
                const Vec3 vContactB2_G = V_GB2[1] + V_GB2[0] % contactPt_B2_G ;
                const Real v = - ~(vContactB1_G-vContactB2_G)*normal_G; // dx/dt

                const Real x32 = x*std::sqrt(x); // x^(3/2)
                const Real f = k*x32*(1 + 1.5*c*v);
                const Vec3 fvec = f * normal_G; // points towards c1

                pe += 0.4*k*x32*x; // i.e., 2/5 k x^(5/2)
                rigidBodyForces[body1] += SpatialVec( contactPt_B1_G % fvec, fvec);
                rigidBodyForces[body2] -= SpatialVec( contactPt_B2_G % fvec, fvec);
            }
        }
    }

    void realizeReaction(const State& s) const {
        // Nothing to compute here.
    }

    HuntCrossleyContactRep* cloneSubsystemRep() const {return new HuntCrossleyContactRep(*this);}
    friend std::ostream& operator<<(std::ostream& o, 
                         const HuntCrossleyContactRep::Parameters&); 
};
// Useless, but required by Value<T>.
std::ostream& operator<<(std::ostream& o, 
                         const HuntCrossleyContactRep::Parameters&) 
{assert(false);return o;}

    //////////////////////////
    // HuntCrossleyContact //
    //////////////////////////


/*static*/ bool 
HuntCrossleyContact::isInstanceOf(const ForceSubsystem& s) {
    return HuntCrossleyContactRep::isA(s.getRep());
}
/*static*/ const HuntCrossleyContact&
HuntCrossleyContact::downcast(const ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<const HuntCrossleyContact&>(s);
}
/*static*/ HuntCrossleyContact&
HuntCrossleyContact::updDowncast(ForceSubsystem& s) {
    assert(isInstanceOf(s));
    return reinterpret_cast<HuntCrossleyContact&>(s);
}

const HuntCrossleyContactRep& 
HuntCrossleyContact::getRep() const {
    return dynamic_cast<const HuntCrossleyContactRep&>(*rep);
}
HuntCrossleyContactRep&       
HuntCrossleyContact::updRep() {
    return dynamic_cast<HuntCrossleyContactRep&>(*rep);
}

HuntCrossleyContact::HuntCrossleyContact() {
    rep = new HuntCrossleyContactRep();
    rep->setMyHandle(*this);
}

int HuntCrossleyContact::addSphere(int body, const Vec3& center,
              const Real& radius,
              const Real& stiffness,
              const Real& dissipation)
{
    return updRep().addSphere(body,center,radius,stiffness,dissipation);
}

int HuntCrossleyContact::addHalfSpace(int body, const UnitVec3& normal,
                 const Real& height,
                 const Real& stiffness,
                 const Real& dissipation)
{
    return updRep().addHalfSpace(body,normal,height,stiffness,dissipation);
}

} // namespace SimTK

