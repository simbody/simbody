#ifndef SimTK_SIMBODY_MOBILIZED_BODY_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_H_

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
 * Contributors: Paul Mitiguy                                                 *
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

/** @file
 * This defines the MobilizedBody class, which associates some mass
 * (the "outboard" body) with a mobilizer and a reference frame (the
 * parent or "inboard" body), already present in a MatterSubsystem.
 *
 * MobilizedBody is an abstract base class handle, with concrete classes defined
 * for each kind of mobilizer. There are a set of built-in mobilizers
 * and a generic "Custom" mobilizer (an actual abstract base class) from
 * which advanced users may derive their own mobilizers.
 */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Body.h"

#include <cassert>

namespace SimTK {

class SimbodyMatterSubsystem;
class MobilizedBody;
class MobilizedBodyImpl;

/**
 * This is the base class for all MobilizedBody classes, which is just a handle for the underlying
 * hidden implementation. Each built-in MobilizedBody type is a local subclass within
 * MobilizedBody, so the built-ins have names like MobilizedBody::Pin. All concrete MobilizedBodies,
 * including the built-ins, are derived from MobilizedBody.
 *
 * SimTK Design Patterns used:
 *    - remote construction object
 *    - abstract private implementation
 *    - binary compatible interface
 *    - custom object interface
 */

// We only want the template instantiation to occur once. This symbol is defined in the SimTK core
// compilation unit that defines the mobilized body class but should not be defined any other time.
#ifndef SimTK_DEFINING_MOBILIZED_BODY
    extern template class PIMPLHandle<MobilizedBody, MobilizedBodyImpl>;
#endif

class SimTK_SIMBODY_EXPORT MobilizedBody : public PIMPLHandle<MobilizedBody, MobilizedBodyImpl> {
public:
    MobilizedBody() { }

    // These declarations are needed so we can have std::vectors of
    // these things, but they aren't implemented.
    //MobilizedBody(const MobilizedBody&) {assert(false);}
    //MobilizedBody& operator=(const MobilizedBody&){assert(false);return*this;}

    // These are the built-in MobilizedBody types. Types on the same line are
    // synonymous. Each of these has a known number of coordinates and speeds 
    // (at least a default number) so
    // can define routines which return and accept specific-size arguments, e.g.
    // Real (for 1-dof mobilizer) and Vec5 (for 5-dof mobilizer). Here is the
    // conventional interface that each built-in should provide. The base type
    // provides similar routines but using variable-sized or "one at a time"
    // arguments. (Vec<1> here will actually be a Real; assume the built-in
    // MobilizedBody class is "BuiltIn")
    //
    //    BuiltIn&       setDefaultQ(const Vec<nq>&);
    //    const Vec<nq>& getDefaultQ() const;
    //
    //    const Vec<nq>& getQ[Dot[Dot]](const State&) const;
    //    const Vec<nu>& getU[Dot](const State&) const;
    //
    //    void setQ(State&, const Vec<nq>&) const;
    //    void setU(State&, const Vec<nu>&) const;
    //
    //    const Vec<nq>& getMyPartQ(const State&, const Vector& qlike) const;
    //    const Vec<nu>& getMyPartU(const State&, const Vector& ulike) const;
    //   
    //    Vec<nq>& updMyPartQ(const State&, Vector& qlike) const;
    //    Vec<nu>& updMyPartU(const State&, Vector& ulike) const;
    //      


    class Pin;         typedef Pin    Torsion;
    class Slider;      typedef Slider Prismatic;
    class Universal;
    class Cylinder;
    class BendStretch;
    class Planar;
    class Gimbal;
    class Ball; typedef Ball Orientation, Spherical;
    class Translation; typedef Translation Cartesian;
    class Free;
    class LineOrientation;
    class FreeLine;
    class Weld;
    class Screw;
    class Ellipsoid;
    class Custom;
    class Ground;
    
    class BallImpl;
    class BendStretchImpl;
    class CustomImpl;
    class CylinderImpl;
    class EllipsoidImpl;
    class FreeImpl;
    class FreeLineImpl;
    class GimbalImpl;
    class GroundImpl;
    class LineOrientationImpl;
    class PinImpl;
    class PlanarImpl;
    class ScrewImpl;
    class SliderImpl;
    class TranslationImpl;
    class UniversalImpl;
    class WeldImpl;

    ///////////////////////////////
    // PAUL'S FRIENDLY INTERFACE //
    ///////////////////////////////

    /// @name High level interface
    /// In the API below, the current MobilizedBody B is the "object" or "main"
    /// body with which we are concerned. Often there will be an additional body mentioned
    /// in the argument list as a target for some conversion. That "auxiliary" body will
    /// be called "body A". The Ground body is abbreviated "G".
    ///
    /// We use OF to mean "the origin of frame F", CB is "the mass center of body B".
    /// R_AF is the rotation matrix giving frame F's
    /// orientation in frame A, such that a vector v expressed in F is reexpressed in
    /// A by v_A = R_AF * v_F. X_AF is the spatial transform giving frame F's origin
    /// location and orientation in frame A, such that a point P whose location is
    /// measured from F's origin OF and expressed in F by vector r_OF_P is remeasured from
    /// frame A's origin and reexpressed in A via r_OA_P = X_AF * r_OF_P. 

    //@{

        // MASS PROPERTIES //

    /// Return the mass properties of body B, measured in the B frame, but expressed
    /// in body A. That is, return the mass, mass center location r_OB_CB, and
    /// the inertia about the body origin OB, expressed in A. If body A is the
    /// same body as body B, then the mass properties can be obtained after realizing
    /// only to the Instance stage, otherwise the state must have been realized
    /// to Position stage.
    ///
    /// If inBodyA==Ground, the returned mass properties are equivalent to
    /// the Spatial Inertia matrix as used in the
    /// Spatial Operator Algebra formulation (that is, the local body mass
    /// properties but expressed in Ground). You can pull out the
    /// individual elements of MassProperties m with m.getMass(), 
    /// m.getMassCenter() and m.getInertia(). You can get them as
    /// a Spatial Inertia Matrix (2x2 x Mat33) with m.toSpatialMat()
    /// or as a 6x6 matrix with m.toMat66().
    ///
    /// @par Required stage
    ///   \c Stage::Instance, if \a inBodyA ==\a objectBodyB
    /// \n\c Stage::Position otherwise.
    MassProperties calcBodyMassPropertiesInBody(const State& s,
                                                const MobilizedBody& inBodyA) const
    {
        const MassProperties& mp = getBodyMassProperties(s);
        if (isSameMobilizedBody(inBodyA))
            return mp;

        // must be at Stage >= Position
        Rotation R_BA = ~getBodyRotation(s);  // R_BG (assume A==G)
        if (!inBodyA.isGround())
            R_BA *= inBodyA.getBodyRotation(s); // R_BA = R_BG*R_GA
        return mp.reexpress(R_BA); // i.e., reexpress from B to A
    }

    /**
     * Return the mass properties of body B, measured from and about
     * the B frame origin, but expressed in Ground and then returned
     * as a Spatial Inertia Matrix. The mass properties are arranged
     * in the SpatialMat like this:                  @verbatim
         M=[      I_OB      crossMat(m*CB) ]
           [ ~crossMat(m*CB)   diag(m)     ]         @endverbatim
     * where I_OB is the inertia taken about the B frame origin OB,
     * and CB is the vector r_OB_CB from B's origin to its mass center.
     *
     * The Spatial Inertia Matrix for Ground has infinite mass and
     * inertia, with the cross terms set to zero. That is, it looks
     * like a 6x6 diagonal matrix with Infinity on the diagonals.
     * 
     * @par Required stage
     *   \c Stage::Position, unless \a objectBodyB == \c GroundIndex
     */
    SpatialMat calcBodySpatialInertiaMatrixInGround(const State& s) const
    {
        if (isGround())
            return SpatialMat(Mat33(Infinity)); // sets diagonals to Inf

        const MassProperties& mp   = getBodyMassProperties(s);
        const Rotation&       R_GB = getBodyRotation(s);
         // re-express in Ground without shifting, convert to spatial mat.
        return mp.reexpress(~R_GB).toSpatialMat();
    }

    /// Calculate the location of body B's mass center, measured from the origin 
    /// of body A, and expressed in the A frame.
    ///
    /// @par Required stage
    ///   \c Stage::Instance, if \a inBodyA == \a thisBodyB
    /// \n\c Stage::Position otherwise.
    Vec3 calcBodyMassCenterLocationInBody(const State& s,
                                          const MobilizedBody& inBodyA) const
    {
        const Vec3& r_OB_CB = getBodyMassCenterStation(s);
        if (inBodyA.isGround()) return locateBodyPointOnGround(s, r_OB_CB);
        return locateBodyPointOnBody(s, r_OB_CB, inBodyA);
    }

    /// Return the central inertia for body B, that is, the inertia taken about
    /// body B's mass center CB, and expressed in B.
    ///
    /// @par Required stage
    ///   \c Stage::Instance
    Inertia calcBodyCentralInertia(const State& s, 
                                   MobilizedBodyIndex objectBodyB) const
    {
        return getBodyMassProperties(s).calcCentralInertia();
    }

    /// Return the inertia of this body B, taken about an arbitrary point PA of body A,
    /// and expressed in body A.
    /// TODO: this needs testing!
    Inertia calcBodyInertiaAboutBodyPoint(const State& s,
                                          const MobilizedBody& inBodyA, 
                                          const Vec3&          aboutLocationOnBodyA) const
    {
        // get B's mass props MB, measured about OB, exp. in B
        const MassProperties& MB_OB_B = getBodyMassProperties(s);

        // Calculate the vector from the body B origin (current "about" point) to the new "about" point PA,
        // expressed in B.
        const Vec3 r_OB_PA = inBodyA.calcBodyPointLocationInBody(s, aboutLocationOnBodyA, *this);

        // Now shift the "about" point for body B's inertia IB to PA, but still expressed in B.
        const Inertia IB_PA_B = MB_OB_B.calcShiftedInertia(r_OB_PA);
        
        // Finally reexpress the inertia in the A frame.
        const Rotation R_BA    = inBodyA.calcBodyRotationFromBody(s, *this);
        const Inertia  IB_PA_A = IB_PA_B.reexpress(R_BA);
        return IB_PA_A;
    }


        // POSITION //

    /// Return X_AB, the spatial transform to body B's frame from body A's frame.
    Transform calcBodyTransformFromBody(const State& s, 
                                        const MobilizedBody& fromBodyA) const
    {
        const Transform& X_GB = getBodyTransform(s);
        if (fromBodyA.isGround()) return X_GB;
        const Transform& X_GA = fromBodyA.getBodyTransform(s);
        return ~X_GA*X_GB;
    }

    /// Return R_AB, the rotation matrix to body B's x,y,z axes from body A's x,y,z axes.
    Rotation calcBodyRotationFromBody(const State& s, 
                                      const MobilizedBody& fromBodyA) const
    {
        if (isSameMobilizedBody(fromBodyA))
            return Rotation();  // Identity rotation; no access to State

        if      (fromBodyA.isGround())  return  getBodyRotation(s);             // R_GB
        else if (this->isGround())      return ~fromBodyA.getBodyRotation(s);   // R_AG (=~R_GA)
        else                            return ~fromBodyA.getBodyRotation(s)    // R_AB=R_AG*R_GB
                                              * getBodyRotation(s);
    }

    /// Return r_OA_OB, the location of body B's origin OB, measured from body A's
    /// origin OA, expressed in body A.
    Vec3 calcBodyOriginLocationInBody(const State& s, 
                                      const MobilizedBody& inBodyA) const
    {
        if (isSameMobilizedBody(inBodyA)) return Vec3(0);

        const Vec3& r_OG_OB = getBodyOriginLocation(s); // from G origin, exp. in G
        if (inBodyA.isGround()) return r_OG_OB;
        else                    return inBodyA.locateGroundPointOnBody(s, r_OG_OB);
    }

    /// Given a vector r_OB_P measured from body B's origin to a point P on body B, expressed in body B,
    /// return the vector r_OA_P measured from body A's origin to point P, expressed in body A.
    Vec3 calcBodyPointLocationInBody(const State& s,
                                     const Vec3&  locationOnBodyB, 
                                     const MobilizedBody& inBodyA) const
    {
        if      (isSameMobilizedBody(inBodyA))  return locationOnBodyB;
        else if (inBodyA.isGround()) return locateBodyPointOnGround(s,locationOnBodyB);
        else if (this->isGround())   return inBodyA.locateGroundPointOnBody(s,locationOnBodyB);
        else                         return locateBodyPointOnBody(s, locationOnBodyB, inBodyA);
    }

    /// Given a vector v_B expressed in body B, return v_A, that same vector reexpressed in body A.
    Vec3 calcBodyVectorInBody(const State& s, 
                              const Vec3&          vectorOnBodyB, 
                              const MobilizedBody& inBodyA) const
    {
        if      (isSameMobilizedBody(inBodyA))  return vectorOnBodyB;
        else if (inBodyA.isGround()) return expressBodyVectorInGround(s,vectorOnBodyB);
        else if (this->isGround())   return inBodyA.expressGroundVectorInBody(s,vectorOnBodyB);
        else                         return expressBodyVectorInBody(s, vectorOnBodyB, inBodyA);
    }

        // VELOCITY //

    /// Return the angular and linear velocity of body B's frame in body A's frame, expressed in body A,
    /// and arranged as a SpatialVec.
    SpatialVec calcBodySpatialVelocityInBody(const State& s,
                                             const MobilizedBody& inBodyA) const
    {
        const SpatialVec& V_GB = getBodyVelocity(s);
        if (inBodyA.isGround()) return V_GB;

        // Body A is not Ground so we'll have to compute relative velocity.

        const SpatialVec& V_GA   = inBodyA.getBodyVelocity(s);
        const Vec3        w_AB_G = V_GB[0]-V_GA[0]; // angular velocity of B in A, exp in G

        // Angular velocity was easy, but for linear velocity we have to add in an wXr term.
        const Transform&  X_GB      = getBodyTransform(s);
        const Transform&  X_GA      = inBodyA.getBodyTransform(s);
        const Vec3        r_OA_OB_G = X_GB.T() - X_GA.T(); // vector from OA to OB, exp in G

        const Vec3 v_AB_G = (V_GB[1]-V_GA[1]) + w_AB_G % r_OA_OB_G; // linear velocity of OB in A, exp in G

        // We're done, but the answer is expressed in Ground. Reexpress in A and return.
        const Rotation& R_GA = X_GA.R();
        return SpatialVec(~R_GA*w_AB_G, ~R_GA*v_AB_G);
    }

    /// Return the angular velocity w_AB of body B's frame in body A's frame, expressed in body A.
    Vec3 calcBodyAngularVelocityInBody(const State& s,
                                       const MobilizedBody& inBodyA) const 
    {
        const SpatialVec& V_GB = getBodyVelocity(s);
        if (inBodyA.isGround()) return V_GB[0];

        // Body A is not Ground so we'll have to compute relative angular velocity.
        const SpatialVec& V_GA   = inBodyA.getBodyVelocity(s);
        const Vec3        w_AB_G = V_GB[0]-V_GA[0]; // angular velocity of B in A, exp in G

        // Now reexpress in A.
        return inBodyA.expressGroundVectorInBody(s, w_AB_G);
    }

    /// Return the velocity of body B's origin point in body A's frame, expressed in body A.
    Vec3 calcBodyOriginVelocityInBody(const State& s,
                                      const MobilizedBody& inBodyA) const
    {
        // Doesn't save much to special case this one.
        return calcBodySpatialVelocityInBody(s,inBodyA)[1];
    }

    /// Return the velocity of a point P fixed on body B, in body A's frame, expressed in body A.
    Vec3 calcBodyFixedPointVelocityInBody(const State& s, 
                                          const Vec3&          locationOnBodyB, 
                                          const MobilizedBody& inBodyA) const;

    /// Return the velocity of a point P moving on body B, in body A's frame, expressed in body A.
    Vec3 calcBodyMovingPointVelocityInBody(const State& s,
                                           const Vec3& locationOnBodyB, 
                                           const Vec3& velocityOnBodyB,
                                           const MobilizedBody& inBodyA) const;

        // ACCELERATION //

    /// Return the angular and linear acceleration of body B's frame in body A's frame, expressed in body A,
    /// and arranged as a SpatialVec.
    SpatialVec calcBodySpatialAccelerationInBody(const State& s,
                                                 const MobilizedBody& inBodyA) const;

    /// Return the angular acceleration of body B's frame in body A's frame, expressed in body A.
    Vec3 calcBodyAngularAccelerationInBody(const State& s,
                                           const MobilizedBody& inBodyA) const;

    /// Return the acceleration of body B's origin point in body A's frame, expressed in body A.
    Vec3 calcBodyOriginAccelerationInBody(const State& s, 
                                          const MobilizedBody& inBodyA) const;

    /// Return the acceleration of a point P fixed on body B, in body A's frame, expressed in body A.
    Vec3 calcBodyFixedPointAccelerationInBody(const State& s,
                                              const Vec3&          locationOnBodyB, 
                                              const MobilizedBody& inBodyA) const;

    /// Return the velocity of a point P moving (and possibly accelerating) on body B,
    /// in body A's frame, expressed in body A.
    Vec3 calcBodyMovingPointAccelerationInBody(const State& s, 
                                               const Vec3&          locationOnBodyB, 
                                               const Vec3&          velocityOnBodyB, 
                                               const Vec3&          accelerationOnBodyB,
                                               const MobilizedBody& inBodyA) const;

        // SCALAR DISTANCE //

    /// Calculate the distance from a point PB on body B to a point PA on body A.
    /// We are given the location vectors (stations) r_OB_PB and r_OA_PA, expressed in
    /// their respective frames. We return |r_OB_OA|.
    Real calcPointToPointDistance(const State& s,
                                  const Vec3&          locationOnBodyB,
                                  const MobilizedBody& bodyA,
                                  const Vec3&          locationOnBodyA) const
    {
        if (isSameMobilizedBody(bodyA))
            return (locationOnBodyA-locationOnBodyB).norm();

        const Vec3 r_OG_PB = locateBodyPointOnGround(s,locationOnBodyB);
        const Vec3 r_OG_PA = bodyA.locateBodyPointOnGround(s,locationOnBodyA);
        return (r_OG_PA - r_OG_PB).norm();
    }

    /// Calculate the time rate of change of distance from a fixed point PB on body B to a fixed point
    /// PA on body A. We are given the location vectors r_OB_PB and r_OA_PA, expressed in their
    /// respective frames. We return d/dt |r_OB_OA|, under the assumption that the time derivatives
    /// of the two given vectors in their own frames is zero.
    Real calcFixedPointToPointDistanceTimeDerivative(const State& s,
                                                     const Vec3&          locationOnBodyB,
                                                     const MobilizedBody& bodyA,
                                                     const Vec3&          locationOnBodyA) const
    {
        if (isSameMobilizedBody(bodyA))
            return 0;

        Vec3 rB, rA, vB, vA;
        calcBodyFixedPointLocationAndVelocityInGround(s,locationOnBodyB,rB,vB);
        bodyA.calcBodyFixedPointLocationAndVelocityInGround(s,locationOnBodyA,rA,vA);
        const Vec3 r = rA-rB, v = vA-vB;
        const Real d = r.norm();

        // When the points are coincident, the rate of change of distance is just their relative speed.
        // Otherwise, it is the speed along the direction of separation. 
        if (d==0) return v.norm();
        else return dot(v, r/d);
    }


    /// Calculate the time rate of change of distance from a moving point PB on body B to a moving point
    /// PA on body A. We are given the location vectors r_OB_PB and r_OA_PA, and the velocities of
    /// PB in B and PA in A, all expressed in their respective frames. We return d/dt |r_OB_OA|,
    /// taking into account the (given) time derivatives of the locations in their local frames, as well
    /// as the relative velocities of the bodies.
    Real calcMovingPointToPointDistanceTimeDerivative(const State& s,
                                                      const Vec3&          locationOnBodyB,
                                                      const Vec3&          velocityOnBodyB,
                                                      const MobilizedBody& bodyA,
                                                      const Vec3&          locationOnBodyA,
                                                      const Vec3&          velocityOnBodyA) const;

    /// Calculate the second time derivative of distance from a fixed point PB on body B to a fixed point
    /// PA on body A. We are given the location vectors (stations) r_OB_PB and r_OA_PA, expressed in their
    /// respective frames. We return d^2/dt^2 |r_OB_OA|, under the assumption that the time derivatives
    /// of the two given vectors in their own frames is zero.
    Real calcFixedPointToPointDistance2ndTimeDerivative(const State& s,
                                                        const Vec3&          locationOnBodyB,
                                                        const MobilizedBody& bodyA,
                                                        const Vec3&          locationOnBodyA) const
    {
        if (isSameMobilizedBody(bodyA))
            return 0;

        Vec3 rB, rA, vB, vA, aB, aA;
        calcBodyFixedPointLocationVelocityAndAccelerationInGround(s,locationOnBodyB,rB,vB,aB);
        bodyA.calcBodyFixedPointLocationVelocityAndAccelerationInGround(s,locationOnBodyA,rA,vA,aA);

        const Vec3 r = rA-rB, v = vA-vB, a = aA-aB;
        const Real d = r.norm();
        
        // This method is the time derivative of calcFixedPointToPointDistanceTimeDerivative(), so it
        // must follow the same two cases. That is, when the points are coincident the change in 
        // separation rate is the time derivative of the speed |v|, otherwise it is the time
        // derivative of the speed along the separation vector.

        if (d==0) {
            // Return d/dt |v|. This has two cases: if |v| is zero, the rate of change of speed is
            // just the points' relative acceleration magnitude. Otherwise, it is the acceleration
            // in the direction of the current relative velocity vector.
            const Real s = v.norm(); // speed
            if (s==0) return a.norm();
            else return dot(a, v/s); // TODO: check with Paul
        }

        // Points are separated.
        const Vec3 u = r/d;             // u is the separation direction (a unit vector from B to A) 
        const Vec3 vp = v - dot(v,u)*u; // velocity perpendicular to separation direction
        return dot(a,u) + dot(vp,v)/d;
    }

    /// Calculate the second time derivative of distance from a moving point PB on body B to a moving point
    /// PA on body A. We are given the location vectors r_OB_PB and r_OA_PA, and the velocities and
    /// accelerations of PB in B and PA in A, all expressed in their respective frames. We return
    /// d^2/dt^2 |r_OA_OB|, taking into account the time derivatives of the locations in their
    /// local frames, as well as the relative velocities and accelerations of the bodies.
    Real calcMovingPointToPointDistance2ndTimeDerivative(const State& s,
                                                         const Vec3&          locationOnBodyB,
                                                         const Vec3&          velocityOnBodyB,
                                                         const Vec3&          accelerationOnBodyB,
                                                         const MobilizedBody& bodyA,
                                                         const Vec3&          locationOnBodyA,
                                                         const Vec3&          velocityOnBodyA,
                                                         const Vec3&          accelerationOnBodyA) const;

    // End of high-level interface.
    //@}



    // Add decorative geometry specified relative to the new (outboard) body's reference
    // frame B, or to the outboard mobilizer frame M attached to body B, or
    // to the inboard mobilizer frame F attached to the parent body P. Note that
    // the body itself may already have had some decorative geometry on it when
    // it was first put into this MobilizedBody; in that case this just adds more.
    MobilizedBody& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)updBody().addDecoration(X_BD,g);
        return *this;
    }
    MobilizedBody& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry&);
    MobilizedBody& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry&);

    // Topology stage (i.e., construction).
    // Calling these means you are (re)constructing the system and will have to do
    // realizeTopology() and extract a new State before doing any analysis.
    const Body&    getBody() const;
    Body&          updBody();
    MobilizedBody& setBody(const Body&);

    MobilizedBody& setDefaultMassProperties(const MassProperties& m) {
        updBody().setDefaultRigidBodyMassProperties(m); // might not be allowed
        return *this;
    }

    const MassProperties& getDefaultMassProperties() const {
        return getBody().getDefaultRigidBodyMassProperties(); // every body type can do this
    }

    MobilizedBody& setDefaultInboardFrame (const Transform& X_PF);
    MobilizedBody& setDefaultOutboardFrame(const Transform& X_BM);
    const Transform& getDefaultInboardFrame()  const; // X_PF
    const Transform& getDefaultOutboardFrame() const; // X_BM

        // Utilities //
    Real  getOneFromQPartition(const State&, int which, const Vector& qlike) const;
    Real& updOneFromQPartition(const State&, int which, Vector& qlike) const;

    Real  getOneFromUPartition(const State&, int which, const Vector& ulike) const;
    Real& updOneFromUPartition(const State&, int which, Vector& ulike) const;

    void applyOneMobilityForce(const State& s, int which, Real f, 
                               Vector& mobilityForces) const
    {
        updOneFromUPartition(s,which,mobilityForces) += f;
    }

    void applyBodyForce(const State& s, const SpatialVec& spatialForceInG, 
                        Vector_<SpatialVec>& bodyForces) const;

    void applyBodyTorque(const State& s, const Vec3& torqueInG, 
                         Vector_<SpatialVec>& bodyForces) const;

    void applyForceToBodyPoint(const State& s, const Vec3& pointInB, const Vec3& forceInG,
                               Vector_<SpatialVec>& bodyForces) const;

        // MODEL STAGE responses //
    int getNumQ(const State&) const;
    int getNumU(const State&) const;

    Real getOneQ(const State&, int which) const;
    Real getOneU(const State&, int which) const;

    Vector getQVector(const State&) const;
    Vector getUVector(const State&) const;

        // VELOCITY STAGE responses //
    Real   getOneQDot   (const State&, int which) const;
    Vector getQDotVector(const State&) const;

        // ACCELERATION STAGE responses //
    Real   getOneUDot   (const State&, int which) const;
    Real   getOneQDotDot(const State&, int which) const;
    Vector getUDotVector   (const State&) const;
    Vector getQDotDotVector(const State&) const;

        // MODEL STAGE solvers //
    void setOneQ(State&, int which, Real) const;
    void setOneU(State&, int which, Real) const;

    void setQVector(State& s, const Vector& v) const;
    void setUVector(State& s, const Vector& v) const;

    // These routines set the generalized coordinates, or speeds (state
    // variables) for just the mobilizer associated with this MobilizedBody
    // (ignoring all other mobilizers and constraints), without requiring knowledge
    // of the meanings of the individual state variables. The idea here
    // is to provide a physically-meaningful quantity relating the 
    // mobilizer's inboard and outboard frames, and then ask the mobilizer
    // to set its state variables to reproduce that quantity to the
    // extent it can.
    //
    // These routines can be called in Stage::Model, however the routines
    // may consult the current values of the state variables in some cases,
    // so you must make sure they have been set to reasonable, or at least
    // innocuous values (zero will work). In no circumstance will any of
    // these routines look at any state variables which belong to another
    // mobilizer; they are limited to working locally with one mobilizer.
    //
    // Routines which specify only translation (linear velocity) may use
    // rotational coordinates to help satisfy the translation requirement.
    // An alternate "Only" method is available to forbid modification of 
    // purely rotational coordinates in that case. When a mobilizer uses
    // state variables which have combined rotational and translational
    // character (e.g. a screw joint) consult the documentation for the
    // mobilizer to find out how it responds to these routines.
    //
    // There is no guarantee that the desired physical quantity will be
    // achieved by these routines; you can check on return if you're
    // worried. Individual mobilizers make specific promises about what
    // they will do; consult the documentation. These routines do not
    // throw exceptions even for absurd requests like specifying a
    // rotation for a sliding mobilizer. Nothing happens if
    // there are no mobilities here, i.e. Ground or a Weld mobilizer.

    void setQToFitTransform      (State&, const Transform& X_FM) const;
    void setQToFitRotation       (State&, const Rotation&  R_FM) const;
    void setQToFitTranslation    (State&, const Vec3&      r_FM) const;
    void setQToFitTranslationOnly(State&, const Vec3&      r_FM) const;

    // Routines which affect generalized speeds u depend on the generalized
    // coordinates q already having been set; they never change coordinates.
    void setUToFitVelocity          (State&, const SpatialVec& V_FM) const;
    void setUToFitAngularVelocity   (State&, const Vec3&       w_FM) const;
    void setUToFitLinearVelocity    (State&, const Vec3&       v_FM) const;
    void setUToFitLinearVelocityOnly(State&, const Vec3&       v_FM) const;


        // INSTANCE STAGE responses //

    const MassProperties& getBodyMassProperties(const State&) const;

    /// Return the mass of this body. Callable at Instance Stage or higher.
    Real getBodyMass(const State& s) const {
        return getBodyMassProperties(s).getMass();
    }

    /// Return this body's center of mass station (i.e., the vector fixed in the body,
    /// going from body origin to body mass center, expressed in the body frame.)
    /// Callable at Instance Stage or higher.
    const Vec3& getBodyMassCenterStation(const State& s) const {
        return getBodyMassProperties(s).getMassCenter();
    }

    /// Return this body's inertia matrix, taken about the body origin and 
    /// expressed in the body frame.
    /// Callable at Instance Stage or higher.
    const Inertia& getBodyInertiaAboutBodyOrigin(const State& s) const {
        return getBodyMassProperties(s).getInertia();
    }

    const Transform& getInboardFrame (const State&) const;  // X_PF
    const Transform& getOutboardFrame(const State&) const;  // X_BM

        // INSTANCE STAGE solvers //

    // Calling these reduces stage to Stage::Model.
    void setInboardFrame (State&, const Transform& X_PF) const;
    void setOutboardFrame(State&, const Transform& X_BM) const;

        // POSITION STAGE responses //

    /// Extract from the state cache the already-calculated spatial configuration of
    /// body B's body frame, measured with respect to the ground frame and expressed
    /// in the ground frame. That is, we return the location of the body frame's
    /// origin, and the orientation of its x, y, and z axes, as the transform X_GB.
    /// This response is available at Position stage.
    const Transform& getBodyTransform(const State&) const; // X_GB

    /// Extract from the state cache the already-calculated spatial orientation
    /// of body B's body frame x, y, and z axes expressed in the ground frame,
    /// as the rotation matrix R_GB. This response is available at Position stage.
    const Rotation& getBodyRotation(const State& s) const {
        return getBodyTransform(s).R();
    }
    /// Extract from the state cache the already-calculated spatial location
    /// of body B's body frame origin OB, measured from the ground origin and
    /// expressed in the ground frame, as the translation vector r_OG_OB.
    /// This response is available at Position stage.
    const Vec3& getBodyOriginLocation(const State& s) const {
        return getBodyTransform(s).T();
    }

    /// At stage Position or higher, return the cross-mobilizer transform.
    /// This is X_FM, the body's inboard mobilizer frame M measured and expressed in
    /// the parent body's corresponding outboard frame F.
    const Transform& getMobilizerTransform(const State&) const; // X_FM


        // VELOCITY STAGE responses //

    /// Extract from the state cache the already-calculated spatial velocity of this
    /// body's reference frame B, measured with respect to the ground frame and expressed
    /// in the ground frame. That is, we return the linear velocity v_GB of the body
    /// frame's origin in G, and the body's angular velocity w_GB as the spatial velocity
    /// vector V_GB = {w_GB, v_GB}. This response is available at Velocity stage.
    const SpatialVec& getBodyVelocity(const State&) const; // V_GB

    /// Extract from the state cache the already-calculated inertial angular
    /// velocity vector w_GB of this body B, measured with respect to the ground frame
    /// and expressed in the ground frame. This response is available at Velocity stage.
    const Vec3& getBodyAngularVelocity(const State& s) const {
        return getBodyVelocity(s)[0]; 
    }
    /// Extract from the state cache the already-calculated inertial linear
    /// velocity vector v_G_OB of this body B's origin point OB, measured with respect
    /// to the ground frame and expressed in the ground frame. This response
    /// is available at Velocity stage.
    const Vec3& getBodyOriginVelocity(const State& s) const {
        return getBodyVelocity(s)[1];
    }

    /// Extract from the state cache the already-calculated inertial angular
    /// acceleration vector aa_GB of this body B, measured with respect to the ground frame
    /// and expressed in the ground frame. This response is available at Acceleration stage.
    const Vec3& getBodyAngularAcceleration(const State& s) const {
        return getBodyAcceleration(s)[0]; 
    }
    /// Extract from the state cache the already-calculated inertial linear
    /// acceleration vector a_G_OB of this body B's origin point OB, measured with respect
    /// to the ground frame and expressed in the ground frame. This response
    /// is available at Acceleration stage.
    const Vec3& getBodyOriginAcceleration(const State& s) const {
        return getBodyAcceleration(s)[1];
    }

    /// Return the Cartesian (ground) location of a station fixed on body B. That is
    /// we return locationOnG = X_GB * locationOnB which means the result is measured from
    /// the ground origin and expressed in ground. Cost is 18 flops. This operator is
    /// available at Position stage.
    Vec3 locateBodyPointOnGround(const State& s, const Vec3& locationOnB) const {
        return getBodyTransform(s) * locationOnB;
    }

    /// Return the station fixed on this body B that is coincident with the given Ground location.
    /// That is we return locationOnB = X_BG * locationOnG, which means the result is measured
    /// from the body origin OB and expressed in the body frame. Cost is 18 flops. This operator
    /// is available at Position stage or higher.
    Vec3 locateGroundPointOnBody(const State& s, const Vec3& locationOnG) const {
        return ~getBodyTransform(s) * locationOnG;
    }

    /// Given a location on this body B, return the location on body A which is at the same location
    /// in space. That is, we return locationOnA = X_AB * locationOnB, which means the result
    /// is measured from the body A origin and expressed in body A. Cost is 36 flops.
    /// This operator is available at Position stage or higher.
    /// Note: if you know that one of the bodies is Ground, use one of the routines above
    /// which is specialized for Ground to avoid half the work.
    Vec3 locateBodyPointOnBody(const State& s, const Vec3& locationOnB, 
                               const MobilizedBody& toBodyA) const
    {
        return toBodyA.locateGroundPointOnBody(s, locateBodyPointOnGround(s,locationOnB));
    }

    /// Return the Cartesian (ground) location of this body B's mass center.
    Vec3 locateBodyMassCenterOnGround(const State& s) const {
        return locateBodyPointOnGround(s,getBodyMassCenterStation(s));
    }

    /// Re-express a vector expressed in this body B's frame into the same vector in G. That is,
    /// we return vectorInG = R_GB * vectorInB. Cost is 15 flops. 
    /// This operator is available at Position stage.
    Vec3 expressBodyVectorInGround(const State& s, const Vec3& vectorInB) const {
        return getBodyRotation(s)*vectorInB;
    }

    /// Re-express a vector expressed in Ground into the same vector expressed in this body B. That is,
    /// we return vectorInB = R_BG * vectorInG. Cost is 15 flops. 
    /// This operator is available at Position stage.
    Vec3 expressGroundVectorInBody(const State& s, const Vec3& vectorInG) const {
        return ~getBodyRotation(s)*vectorInG;
    }

    /// Re-express a vector expressed in this body B into the same vector expressed in body A.
    /// That is, we return vectorInA = R_AB * vectorInB. Cost is 30 flops.
    /// This operator is available at Position stage.
    /// Note: if you know one of the bodies is Ground, call one of the specialized methods
    /// above to save 15 flops.
    Vec3 expressBodyVectorInBody(const State& s, const Vec3& vectorInB,
                                 const MobilizedBody& inBodyA) const
    {
        return inBodyA.expressGroundVectorInBody(s, expressBodyVectorInGround(s,vectorInB));
    }

    /// Calculate this body B's mass properties, measured in the body B frame, taken about the body
    /// B origin OB, but reexpressed in Ground.
    MassProperties expressBodyMassPropertiesInGround(const State& s) {
            const MassProperties& M_OB_B = getBodyMassProperties(s);
            const Rotation&       R_GB   = getBodyRotation(s);
            return M_OB_B.reexpress(~R_GB);
    }

    /// Calculate body B's momentum (angular, linear) measured and expressed in ground, but taken about
    /// the body origin OB.
    SpatialVec calcBodyMomentumAboutBodyOriginInGround(const State& s) {
        const MassProperties M_OB_G = expressBodyMassPropertiesInGround(s);
        const SpatialVec&    V_GB   = getBodyVelocity(s);
        return M_OB_G.toSpatialMat() * V_GB;
    }

    /// Calculate body B's momentum (angular, linear) measured and expressed in ground, but taken about
    /// the body mass center CB.
    SpatialVec calcBodyMomentumAboutBodyMassCenterInGround(const State& s) const {
        const MassProperties& M_OB_B = getBodyMassProperties(s);
        const Rotation&       R_GB   = getBodyRotation(s);

        // Given a central inertia matrix I, angular velocity w, and mass center velocity v,
        // the central angular momentum is Iw and linear momentum is mv.
        const Inertia I_CB_B = M_OB_B.calcCentralInertia();
        const Inertia I_CB_G = I_CB_B.reexpress(~R_GB);
        const Real    mb     = M_OB_B.getMass();
        const Vec3&   w_GB   = getBodyAngularVelocity(s);
        Vec3          v_G_CB = calcBodyFixedPointVelocityInGround(s, M_OB_B.getMassCenter());

        return SpatialVec( I_CB_G*w_GB, mb*v_G_CB );
    }

    /// Given a station fixed on body B, return its inertial (Cartesian) velocity,
    /// that is, its velocity relative to the ground frame, expressed in the
    /// ground frame. Cost is 27 flops. This operator is available at Velocity stage.
    Vec3 calcBodyFixedPointVelocityInGround(const State& s, const Vec3& stationOnB) const {
        const Vec3& w = getBodyAngularVelocity(s); // in G
        const Vec3& v = getBodyOriginVelocity(s);  // in G
        const Vec3  r = expressBodyVectorInGround(s,stationOnB); // 15 flops
        return v + w % r;                                        // 12 flops
    }

    /// It is cheaper to calculate a station's ground location and velocity together
    /// than to do them separately. Here we can return them both in 30 flops, vs. 45 to
    /// do them in two calls.
    void calcBodyFixedPointLocationAndVelocityInGround(const State& s, const Vec3& locationOnB,
                                                       Vec3& locationOnGround, Vec3& velocityInGround) const
    {
        const Vec3& r_G_OB = getBodyOriginLocation(s);
        const Vec3  r      = expressBodyVectorInGround(s,locationOnB); // 15 flops
        locationOnGround = r_G_OB + r;   // 3 flops

        const Vec3& w = getBodyAngularVelocity(s); // in G
        const Vec3& v = getBodyOriginVelocity(s);  // in G
        velocityInGround = v + w % r; // 12 flops
    }


    /// Given a station fixed on body B, return its inertial (Cartesian) acceleration,
    /// that is, its acceleration relative to the ground frame, expressed in the
    /// ground frame. Cost is 48 flops. This operator is available at Acceleration stage.
    Vec3 calcBodyFixedPointAccelerationInGround(const State& s, const Vec3& stationOnB) const {
        const Vec3& w  = getBodyAngularVelocity(s);     // in G
        const Vec3& aa = getBodyAngularAcceleration(s); // in G
        const Vec3& a  = getBodyOriginAcceleration(s);  // in G

        const Vec3  r = expressBodyVectorInGround(s,stationOnB); // 15 flops
        return a + aa % r + w % (w % r);                         // 33 flops
    }

    /// It is cheaper to calculate a station's ground location, velocity, and acceleration together
    /// than to do them separately. Here we can return them all in 54 flops, vs. 93 to
    /// do them in three calls. This operator is available at Acceleration stage.
    void calcBodyFixedPointLocationVelocityAndAccelerationInGround
       (const State& s, const Vec3& locationOnB,
        Vec3& locationOnGround, Vec3& velocityInGround, Vec3& accelerationInGround) const
    {
        const Rotation&  R_GB   = getBodyRotation(s);
        const Vec3&      r_G_OB = getBodyOriginLocation(s);

        const Vec3 r = R_GB*locationOnB; // re-express station vector in G (15 flops)
        locationOnGround = r_G_OB + r;   // 3 flops

        const Vec3& w  = getBodyAngularVelocity(s); // in G
        const Vec3& v  = getBodyOriginVelocity(s);  // in G
        const Vec3& aa = getBodyAngularAcceleration(s); // in G
        const Vec3& a  = getBodyOriginAcceleration(s);  // in G

        const Vec3 wXr = w % r; // "whipping" velocity w X r due to angular velocity (9 flops)
        velocityInGround     = v + wXr;              // v + w X r (3 flops)
        accelerationInGround = a + aa % r + w % wXr; // 24 flops
    }

    /// Given a station fixed on body B, return its velocity relative to the body frame of
    /// body A, and expressed in body A's body frame. Cost is 54 flops.
    /// This operator is available at Velocity stage.
    /// TODO: UNTESTED!!
    /// TODO: maybe these between-body routines should return results in ground so that they
    /// can be easily combined. Easy to re-express vector afterwards.
    Vec3 calcStationVelocityInBody(const State& s, const Vec3& stationOnB, const MobilizedBody& bodyA) const {
        // If body B's origin were coincident with body A's, then Vdiff_AB would be the relative angular
        // and linear velocity of body B in body A, expressed in G. To get the point we're interested in,
        // we need the vector from body A's origin to stationB to account for the extra linear velocity
        // that will be created by moving away from the origin, due to the bodies' relative angular
        // velocity.
        const SpatialVec Vdiff_AB = getBodyVelocity(s) - bodyA.getBodyVelocity(s); // 6

        // This is a vector from body A's origin to the point of interest, expressed in G.
        const Vec3 stationA_G = locateBodyPointOnGround(s,stationOnB) - bodyA.getBodyOriginLocation(s); // 21
        const Vec3 v_AsB_G = Vdiff_AB[1] + Vdiff_AB[0] % stationA_G; // 12
        return ~bodyA.getBodyRotation(s) * v_AsB_G; // 15
    }

    /// At stage Velocity or higher, return the cross-mobilizer velocity.
    /// This is V_FM, the relative velocity of the body's inboard mobilizer
    /// frame M in the parent body's corresponding outboard frame F, 
    /// measured and expressed in F. Note that this isn't the usual 
    /// spatial velocity since it isn't expressed in G.
    const SpatialVec& getMobilizerVelocity(const State&) const; // V_FM

        // ACCELERATION STAGE responses //

    /// Extract from the state cache the already-calculated spatial acceleration of
    /// this body's reference frame B, measured with respect to the ground frame and expressed
    /// in the ground frame. That is, we return the linear acceleration a_GB of the body
    /// frame's origin in G, and the body's angular acceleration alpha_GB as the spatial acceleration
    /// vector A_GB = {alpha_GB, a_GB}. This response is available at Acceleration stage.
    const SpatialVec& getBodyAcceleration(const State& s) const; // A_GB

    // Implicit conversion to MobilizedBodyIndex when needed.
    operator MobilizedBodyIndex() const {return getMobilizedBodyIndex();}
    MobilizedBodyIndex        getMobilizedBodyIndex()     const; // id of this mobilized body
    const MobilizedBody&   getParentMobilizedBody() const; // the inboard body (not allowed if this is ground)
    const MobilizedBody&   getBaseMobilizedBody()   const; // the lowest numbered ancestor body on this branch
                                                           //   (returns Ground if if this is Ground)

    // These will fail unless this MobilizedBody is owned by a MatterSubsystem.
    const SimbodyMatterSubsystem& getMatterSubsystem()      const;
    SimbodyMatterSubsystem&       updMatterSubsystem();

    bool isInSubsystem() const;
    bool isInSameSubsystem(const MobilizedBody&) const;

    bool isSameMobilizedBody(const MobilizedBody&) const;
    bool isGround() const; // meaning mobilizedbody 0 of some subsystem, not just the Body type

    // Return this body's level in the tree of bodies, starting with ground at 0,
    // bodies directly connected to ground at 1, bodies directly connected to those at 2, 
    // etc. This is callable after realizeTopology(). This is the graph distance of
    // the body from Ground.
    int getLevelInMultibodyTree() const;
    
    /// Create a new MobilizedBody which is identical to this one, except that it has a
    /// different parent (and consequently might belong to a different MultibodySystem).
    MobilizedBody& cloneForNewParent(MobilizedBody& parent) const;

    // Internal use only

    explicit MobilizedBody(MobilizedBodyImpl* r) : HandleBase(r) { }
};

/// One mobility -- rotation about the common z axis of the inboard
/// and outboard mobilizer frames.
/// Synonym: Torsion
class SimTK_SIMBODY_EXPORT MobilizedBody::Pin : public PIMPLDerivedHandle<Pin, PinImpl, MobilizedBody> {
public:
        // SPECIALIZED INTERFACE FOR PIN MOBILIZER

    // "Angle" is just a nicer name for a pin joint's lone generalized coordinate q.
    Pin& setDefaultAngle(Real angleInRadians) {return setDefaultQ(angleInRadians);}
    Real getDefaultAngle() const              {return getDefaultQ();}

        // Friendly, mobilizer-specific access to generalized coordinates and speeds.

    void setAngle(State& s, Real angleInRadians) {setQ(s, angleInRadians);}
    Real getAngle(const State& s) const {return getQ(s);}

    void setRate(State& s, Real rateInRadiansPerTime) {setU(s, rateInRadiansPerTime);}
    Real getRate(const State& s) const {return getU(s);}

    // Mobility forces are "u-like", that is, one per dof.
    Real getAppliedPinTorque(const State& s, const Vector& mobilityForces) const {
        return getMyPartU(s,mobilityForces);
    }
    void applyPinTorque(const State& s, Real torque, Vector& mobilityForces) const {
        updMyPartU(s,mobilityForces) += torque;
    }

        // STANDARDIZED MOBILIZED BODY INTERFACE

        // required constructors
    Pin();
    Pin(MobilizedBody& parent, const Body&);
    Pin(MobilizedBody& parent, const Transform& inbFrame,
        const Body&,           const Transform& outbFrame);

        // access to generalized coordinates q and generalized speeds u
    Pin& setDefaultQ(Real);
    Real getDefaultQ() const;

    Real getQ(const State&) const;
    Real getQDot(const State&) const;
    Real getQDotDot(const State&) const;
    Real getU(const State&) const;
    Real getUDot(const State&) const;

    void setQ(State&, Real) const;
    void setU(State&, Real) const;

    Real getMyPartQ(const State&, const Vector& qlike) const;
    Real getMyPartU(const State&, const Vector& ulike) const;
   
    Real& updMyPartQ(const State&, Vector& qlike) const;
    Real& updMyPartU(const State&, Vector& ulike) const;

        // specialize return type for convenience
    Pin& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g)
      { (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this; }
    Pin& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g)
      { (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this; }
    Pin& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g)
      { (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this; }
    Pin& setDefaultInboardFrame(const Transform& X_PF)
      { (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this; }
    Pin& setDefaultOutboardFrame(const Transform& X_BM)
      { (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this; }
};

/// One mobility -- translation along the common x axis of the
/// inboard and outboard mobilizer frames.
/// Synonym: Prismatic
class SimTK_SIMBODY_EXPORT MobilizedBody::Slider : public PIMPLDerivedHandle<Slider, SliderImpl, MobilizedBody> {
public:
        // SPECIALIZED INTERFACE FOR SLIDER MOBILIZER

    // "Length" is just a nicer name for a sliding joint's lone generalized coordinate q.
    Slider& setDefaultLength(Real length) {return setDefaultQ(length);}
    Real getDefaultLength() const         {return getDefaultQ();}

        // Friendly, mobilizer-specific access to generalized coordinates and speeds.

    void setLength(State& s, Real length) {setQ(s, length);}
    Real getLength(const State& s) {return getQ(s);}

    void setRate(State& s, Real rateInLengthPerTime) {setU(s, rateInLengthPerTime);}
    Real getRate(const State& s) {return getU(s);}

    // Mobility forces are "u-like", that is, one per dof.
    Real getAppliedForce(const State& s, const Vector& mobilityForces) const {
        return getMyPartU(s,mobilityForces);
    }
    void applyForce(const State& s, Real force, Vector& mobilityForces) const {
        updMyPartU(s,mobilityForces) += force;
    }

        // STANDARDIZED MOBILIZED BODY INTERFACE

        // required constructors
    Slider();
    Slider(MobilizedBody& parent, const Body&);
    Slider(MobilizedBody& parent, const Transform& inbFrame,
        const Body&,           const Transform& outbFrame);

        // access to generalized coordinates q and generalized speeds u
    Slider& setDefaultQ(Real);
    Real getDefaultQ() const;

    Real getQ(const State&) const;
    Real getQDot(const State&) const;
    Real getQDotDot(const State&) const;
    Real getU(const State&) const;
    Real getUDot(const State&) const;

    void setQ(State&, Real) const;
    void setU(State&, Real) const;

    Real getMyPartQ(const State&, const Vector& qlike) const;
    Real getMyPartU(const State&, const Vector& ulike) const;
   
    Real& updMyPartQ(const State&, Vector& qlike) const;
    Real& updMyPartU(const State&, Vector& ulike) const;

        // specialize return type for convenience
    Slider& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g)
      { (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this; }
    Slider& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g)
      { (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this; }
    Slider& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g)
      { (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this; }
    Slider& setDefaultInboardFrame(const Transform& X_PF)
      { (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this; }
    Slider& setDefaultOutboardFrame(const Transform& X_BM)
      { (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this; }
};

/// One mobility -- coordinated rotation and translation along the
/// common z axis of the inboard and outboard mobilizer frames. A
/// "pitch" is specified relating the two. The generalized coordinate
/// q is the rotation angle in radians, the translation is always
/// pitch*q.
class SimTK_SIMBODY_EXPORT MobilizedBody::Screw : public PIMPLDerivedHandle<Screw, ScrewImpl, MobilizedBody> {
public:
    Screw(Real pitch);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Screw(MobilizedBody& parent, const Body&, Real pitch);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Screw(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame,
         Real pitch);

    Screw& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Screw& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Screw& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Screw& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Screw& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    Screw& setDefaultPitch(Real pitch);
    Real   getDefaultPitch() const;

    Screw& setDefaultQ(Real);
    Real   getDefaultQ() const;

    Real getQ(const State&) const;
    Real getQDot(const State&) const;
    Real getQDotDot(const State&) const;
    Real getU(const State&) const;
    Real getUDot(const State&) const;

    void setQ(State&, Real) const;
    void setU(State&, Real) const;

    Real getMyPartQ(const State&, const Vector& qlike) const;
    Real getMyPartU(const State&, const Vector& ulike) const;
   
    Real& updMyPartQ(const State&, Vector& qlike) const;
    Real& updMyPartU(const State&, Vector& ulike) const;
};

/// Two mobilities -- rotation about the x axis, followed by a rotation
/// about the new y axis. This mobilizer is badly behaved when the
/// second rotation is near 90 degrees.
class SimTK_SIMBODY_EXPORT MobilizedBody::Universal : public PIMPLDerivedHandle<Universal, UniversalImpl, MobilizedBody> {
public:
    Universal();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Universal(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Universal(MobilizedBody& parent, const Transform& inbFrame,
              const Body&,           const Transform& outbFrame);

    Universal& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Universal& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Universal& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Universal& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Universal& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }
};

/// Two mobilities -- rotation and translation along the common z axis
/// of the inboard and outboard mobilizer frames.
class SimTK_SIMBODY_EXPORT MobilizedBody::Cylinder : public PIMPLDerivedHandle<Cylinder, CylinderImpl, MobilizedBody> {
public:
    Cylinder();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Cylinder(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Cylinder(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame);

    Cylinder& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Cylinder& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Cylinder& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Cylinder& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Cylinder& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }
};

/// Two mobilities: The z axis of the parent's F frame is 
/// used for rotation (and that is always aligned with the M frame z axis).
/// The x axis of the *M* (outboard) frame is then used for translation;
/// that is, first we rotate around z, which moves M's x with respect to F's x. Then
/// we slide along the rotated x axis. The two generalized coordinates are the
/// rotation and the translation, in that order.
class SimTK_SIMBODY_EXPORT MobilizedBody::BendStretch : public PIMPLDerivedHandle<BendStretch, BendStretchImpl, MobilizedBody> {
public:
    BendStretch();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    BendStretch(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    BendStretch(MobilizedBody& parent, const Transform& inbFrame,
                const Body&,           const Transform& outbFrame);

    BendStretch& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    BendStretch& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    BendStretch& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    BendStretch& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    BendStretch& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }
};

/// Three mobilities -- z rotation and x,y translation. The generalized
/// coordinates are rotation about the shared z axis of the F and M
/// frame, translation along the F frame's x axis, and translation along
/// its y axis, in that order.
class SimTK_SIMBODY_EXPORT MobilizedBody::Planar : public PIMPLDerivedHandle<Planar, PlanarImpl, MobilizedBody> {
public:
    Planar();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Planar(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Planar(MobilizedBody& parent, const Transform& inbFrame,
           const Body&,           const Transform& outbFrame);

    Planar& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Planar& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Planar& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Planar& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Planar& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    // Friendly, mobilizer-specific access to coordinates and speeds.
    Planar& setDefaultAngle(Real a) {
        Vec3 q = getDefaultQ(); q[0] = a; setDefaultQ(q);
        return *this;
    }
    Planar& setDefaultTranslation(const Vec2& r) {
        Vec3 q = getDefaultQ(); q.updSubVec<2>(1) = r; setDefaultQ(q);
        return *this;
    }

    Real getDefaultAngle() const {return getDefaultQ()[0];}
    const Vec2& getDefaultTranslation() const {return getDefaultQ().getSubVec<2>(1);}

    void setAngle      (State& s, Real        a) {setOneQ(s,0,a);}
    void setTranslation(State& s, const Vec2& r) {setOneQ(s,1,r[0]); setOneQ(s,2,r[1]);}

    Real getAngle(const State& s) const {return getQ(s)[0];}
    const Vec2& getTranslation(const State& s) const {return getQ(s).getSubVec<2>(1);}

    // Generic default state Topology methods.
    const Vec3& getDefaultQ() const;
    Planar& setDefaultQ(const Vec3& q);

    const Vec3& getQ(const State&) const;
    const Vec3& getQDot(const State&) const;
    const Vec3& getQDotDot(const State&) const;
    const Vec3& getU(const State&) const;
    const Vec3& getUDot(const State&) const;

    void setQ(State&, const Vec3&) const;
    void setU(State&, const Vec3&) const;

    const Vec3& getMyPartQ(const State&, const Vector& qlike) const;
    const Vec3& getMyPartU(const State&, const Vector& ulike) const;
   
    Vec3& updMyPartQ(const State&, Vector& qlike) const;
    Vec3& updMyPartU(const State&, Vector& ulike) const;
};

/// Three mobilities -- unrestricted orientation modeled as a 1-2-3
/// body-fixed Euler angle sequence. This is singular when the middle
/// angle is 90 degrees.
class SimTK_SIMBODY_EXPORT MobilizedBody::Gimbal : public PIMPLDerivedHandle<Gimbal, GimbalImpl, MobilizedBody> {
public:
    Gimbal();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Gimbal(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Gimbal(MobilizedBody& parent, const Transform& inbFrame,
           const Body&,           const Transform& outbFrame);

    Gimbal& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Gimbal& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Gimbal& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Gimbal& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Gimbal& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    // This is just a nicer name for the generalized coordinate.
    Gimbal& setDefaultRotation(const Rotation& R_FM) {
        return setDefaultQ(R_FM.convertRotationToBodyFixedXYZ());
    }
    Rotation getDefaultRotation() const {
        const Vec3& q = getDefaultQ();
        return Rotation(BodyRotationSequence,
            q[0], XAxis, q[1], YAxis, q[2], ZAxis);
    }

    // This is used only for visualization.
    Gimbal& setDefaultRadius(Real r);
    Real getDefaultRadius() const;

    // Generic default state Topology methods.
    const Vec3& getDefaultQ() const; // X,Y,Z body-fixed Euler angles
    Gimbal& setDefaultQ(const Vec3& q);

    const Vec3& getQ(const State&) const;
    const Vec3& getQDot(const State&) const;
    const Vec3& getQDotDot(const State&) const;
    const Vec3& getU(const State&) const;
    const Vec3& getUDot(const State&) const;

    void setQ(State&, const Vec3&) const;
    void setU(State&, const Vec3&) const;

    const Vec3& getMyPartQ(const State&, const Vector& qlike) const;
    const Vec3& getMyPartU(const State&, const Vector& ulike) const;
   
    Vec3& updMyPartQ(const State&, Vector& qlike) const;
    Vec3& updMyPartU(const State&, Vector& ulike) const;
};

/// Three mobilities -- unrestricted orientation modeled with a
/// quaternion which is never singular. A modeling option allows the
/// joint to use a 1-2-3 Euler sequence (identical to a Gimbal) 
/// instead.
class SimTK_SIMBODY_EXPORT MobilizedBody::Ball : public PIMPLDerivedHandle<Ball, BallImpl, MobilizedBody> {
public:
    explicit Ball();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Ball(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Ball(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame);

    Ball& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Ball& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Ball& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Ball& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Ball& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    // This is just a nicer name for the generalized coordinate.
    Ball& setDefaultRotation(const Rotation& R_FM) {
        return setDefaultQ(R_FM.convertRotationToQuaternion());
    }
    Rotation getDefaultRotation() const {return Rotation(getDefaultQ());}

    // This is used only for visualization.
    Ball& setDefaultRadius(Real r);
    Real getDefaultRadius() const;

    // Generic default state Topology methods.
    const Quaternion& getDefaultQ() const;
    Ball& setDefaultQ(const Quaternion& q);

    const Vec4& getQ(const State&) const;
    const Vec4& getQDot(const State&) const;
    const Vec4& getQDotDot(const State&) const;
    const Vec3& getU(const State&) const;
    const Vec3& getUDot(const State&) const;

    void setQ(State&, const Vec4&) const;
    void setU(State&, const Vec3&) const;

    const Vec4& getMyPartQ(const State&, const Vector& qlike) const;
    const Vec3& getMyPartU(const State&, const Vector& ulike) const;
   
    Vec4& updMyPartQ(const State&, Vector& qlike) const;
    Vec3& updMyPartU(const State&, Vector& ulike) const;
};

/// Three mobilities -- coordinated rotation and translation along the
/// surface of an ellipsoid fixed to the parent (inboard) body.
/// The generalized coordinates are the same as for a Ball (Orientation)
/// joint, that is, a quaternion or 1-2-3 Euler sequence.
class SimTK_SIMBODY_EXPORT MobilizedBody::Ellipsoid : public PIMPLDerivedHandle<Ellipsoid, EllipsoidImpl, MobilizedBody> {
public:
    // The ellipsoid is placed on the mobilizer's inboard frame F, with
    // half-axis dimensions along F's x,y,z respectively.
    Ellipsoid(); // not very useful until radii are set, but has some defaults
    explicit Ellipsoid(const Vec3& radii);
    Ellipsoid(Real a, Real b, Real c);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Ellipsoid(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Ellipsoid(MobilizedBody& parent, const Transform& inbFrame,
              const Body&,           const Transform& outbFrame);

    Ellipsoid& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Ellipsoid& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Ellipsoid& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Ellipsoid& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Ellipsoid& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    // This is just a nicer name for the generalized coordinate.
    Ellipsoid& setDefaultRotation(const Rotation& R_FM) {
        return setDefaultQ(R_FM.convertRotationToQuaternion());
    }
    Rotation getDefaultRotation() const {return Rotation(getDefaultQ());}

    Ellipsoid& setDefaultRadii(const Vec3& r);
    const Vec3& getDefaultRadii() const;

    // Generic default state Topology methods.
    const Quaternion& getDefaultQ() const;
    Quaternion& updDefaultQ();
    Ellipsoid& setDefaultQ(const Quaternion& q) {updDefaultQ()=q; return *this;}
};

/// Three translational mobilities. The generalized coordinates are
/// x,y,z translations along the parent (inboard) F frame axes.
class SimTK_SIMBODY_EXPORT MobilizedBody::Translation : public PIMPLDerivedHandle<Translation, TranslationImpl, MobilizedBody> {
public:
    Translation();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Translation(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Translation(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame);

    Translation& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Translation& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Translation& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Translation& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Translation& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    // This is just a nicer name for setting this mobilizer's generalized coordinates,
    // which together constitute the vector from the F frame's origin to the M
    // frame's origin, expressed in F.

    // Set the topological default values for the initial q's.
    Translation& setDefaultTranslation(const Vec3& p_FM) {
        return setDefaultQ(p_FM);
    }

    // Get the topological default values for the initial q's.
    const Vec3& getDefaultTranslation() const {
        return getDefaultQ();
    }

    // Set the current value of q's in the given State. Note that this is
    // the *cross-mobilizer* translation, not location in the Ground frame.
    void setMobilizerTranslation(State& s, const Vec3& p_FM) const {
        setQ(s,p_FM);
    }

    // Get the current value of the q's for this mobilizer from the given State.
    const Vec3& getMobilizerTranslation(const State& s) const {
        return getQ(s);
    }


    // Set the current value of u's in the given State. Note that this is
    // the *cross-mobilizer* velocity v_FM, not velocity in the Ground frame.
    void setMobilizerVelocity(State& s, const Vec3& v_FM) const {
        setU(s,v_FM);
    }

    // Get the current value of the u's for this mobilizer from the given State.
    const Vec3& getMobilizerVelocity(const State& s) const {
        return getU(s);
    }

    // Get the value of the udot's for this mobilizer from the given State.
    const Vec3& getMobilizerAcceleration(const State& s) const {
        return getUDot(s);
    }

    // Generic default state Topology methods.
    const Vec3& getDefaultQ() const;
    Translation& setDefaultQ(const Vec3& q);

    const Vec3& getQ(const State&) const;
    const Vec3& getQDot(const State&) const;
    const Vec3& getQDotDot(const State&) const;
    const Vec3& getU(const State&) const;
    const Vec3& getUDot(const State&) const;

    void setQ(State&, const Vec3&) const;
    void setU(State&, const Vec3&) const;

    const Vec3& getMyPartQ(const State&, const Vector& qlike) const;
    const Vec3& getMyPartU(const State&, const Vector& ulike) const;
   
    Vec3& updMyPartQ(const State&, Vector& qlike) const;
    Vec3& updMyPartU(const State&, Vector& ulike) const;
};

/// Unrestricted motion for a rigid body (six mobilities). Orientation
/// is modeled the same as for the Orientation mobilizer, that is, using
/// quaternions to avoid singularities. A modeling option exists to 
/// have the joint modeled with a 1-2-3 body fixed Euler sequence like
/// a Gimbal mobilizer. Translational generalized coordinates are
/// x,y,z translations along the F (inboard) axes.
class SimTK_SIMBODY_EXPORT MobilizedBody::Free : public PIMPLDerivedHandle<Free, FreeImpl, MobilizedBody> {
public:
    Free();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Free(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Free(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame);

    Free& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Free& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Free& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Free& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Free& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    // Leaves rotation unchanged.
    Free& setDefaultTranslation(const Vec3&);

    // Leaves translation unchanged. The internal representation is a quaternion
    // so we guarantee that the stored value is numerically identical to the
    // supplied one.
    Free& setDefaultQuaternion(const Quaternion&);

    // Leaves translation unchanged. The Rotation matrix will be converted to
    // a quaternion for storage.
    Free& setDefaultRotation(const Rotation&);
    // Sets both translation and rotation. The Rotation part of the Transform 
    // will be converted to a quaternion for storage.
    Free& setDefaultTransform(const Transform&);

    // These return references to the stored default values.
    const Vec3& getDefaultTranslation() const;
    const Quaternion& getDefaultQuaternion() const;

    // These next two are derived from the stored values.
    Rotation getDefaultRotation() const {
        return Rotation(getDefaultQuaternion());
    }
    Transform getDefaultTransform() const {
        return Transform(Rotation(getDefaultQuaternion()), getDefaultTranslation());
    }

    // Generic default state Topology methods.

    // Returns (Vec4,Vec3) where the Vec4 is a normalized quaternion.
    const Vec7& getDefaultQ() const;

    // Interprets the supplied q as (Vec4,Vec3) where the Vec4 is a possibly
    // unnormalized quaternion. The quaternion will be normalized before it is
    // stored here, so you may not get back exactly the value supplied here if
    // you call getDefaultQ().
    Free& setDefaultQ(const Vec7& q);

    // Note that there is no guarantee that the quaternion part of the returned Q is normalized.
    const Vec7& getQ(const State&) const;
    const Vec7& getQDot(const State&) const;
    const Vec7& getQDotDot(const State&) const;

    const Vec6& getU(const State&) const;
    const Vec6& getUDot(const State&) const;

    // The Q's in the state are set exactly as supplied without normalization.
    void setQ(State&, const Vec7&) const;
    void setU(State&, const Vec6&) const;

    const Vec7& getMyPartQ(const State&, const Vector& qlike) const;
    const Vec6& getMyPartU(const State&, const Vector& ulike) const;
   
    Vec7& updMyPartQ(const State&, Vector& qlike) const;
    Vec6& updMyPartU(const State&, Vector& ulike) const;
};


// These are special "ball" and "free" joints designed to allow arbitrary orientations
// for "linear" bodies, such as a CO2 molecule consisting only of point masses arranged
// along a straight line. Such bodies have no inertia about the line and cause singularities
// in the equations of motion if attached to Orientation or Free mobilizers. Instead, use the
// LineOrientation and LineFree moblizers, making sure that the inertialess direction is
// along the outboard body's z axis (that is, Mz). These mobilizers introduce only two
// mobilities (generalized speeds u), being incapable of representing non-zero angular
// velocity of M in F about Mz. The generalized speeds are in fact the wx and wy 
// components of w_FM_M, that is, the x and y components of the angular velocity of M
// in F *expressed in M*. However, at least three generalized coordinates (q's)
// are required to represent the orientation. By default we use four quaternions for
// unconditional stability. Alternatively, you can request a 1-2-3 body fixed 
// Euler angle sequence (that is, about x, then new y, then new z) which will
// suffer a singularity when the y rotation is 90 degrees since that aligns the
// first rotation axis (x) with the last (z) which is the inertialess direction.

/// Two mobilities, representing unrestricted orientation for a body which is
/// inertialess along its own z axis. The generalized coordinates are the same
/// as for the general Orientation (Ball) mobilizer, but there are only
/// two generalized speeds. These are the x,y components of the angular velocity
/// of frame M in F, but expressed in the *M* (outboard frame).
class SimTK_SIMBODY_EXPORT MobilizedBody::LineOrientation : public PIMPLDerivedHandle<LineOrientation, LineOrientationImpl, MobilizedBody> {
public:
    LineOrientation();


    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    LineOrientation(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    LineOrientation(MobilizedBody& parent, const Transform& inbFrame,
                    const Body&,           const Transform& outbFrame);

    LineOrientation& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    LineOrientation& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    LineOrientation& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    LineOrientation& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    LineOrientation& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }
};

/// Five mobilities, representing unrestricted motion for a body which is
/// inertialess along its own z axis. The rotational generalized coordinates are the same
/// as for the LineOrientation mobilizer. The translational coordinates are
/// the same as in a Free mobilizer, or a Cartesian (Translation) mobilizer.
class SimTK_SIMBODY_EXPORT MobilizedBody::FreeLine : public PIMPLDerivedHandle<FreeLine, FreeLineImpl, MobilizedBody> {
public:
    FreeLine();

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    FreeLine(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    FreeLine(MobilizedBody& parent, const Transform& inbFrame,
             const Body&,           const Transform& outbFrame);

    FreeLine& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    FreeLine& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    FreeLine& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    FreeLine& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    FreeLine& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }
};

/// Zero mobilities. This degenerate "mobilizer" serves only to weld together
/// the M frame of a body to the F frame on its parent.
/// TODO: not implemented yet.
class SimTK_SIMBODY_EXPORT MobilizedBody::Weld : public PIMPLDerivedHandle<Weld, WeldImpl, MobilizedBody> {
public:
    Weld();


    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Weld(MobilizedBody& parent, const Body&);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Weld(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame);

    Weld& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    Weld& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    Weld& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    Weld& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    Weld& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }
};


/// This is a special type of "mobilized" body used as a placeholder for Ground
/// in the 0th slot for a MatterSubsystem's mobilized bodies.
/// The body type will also be Ground.
class SimTK_SIMBODY_EXPORT MobilizedBody::Ground : public PIMPLDerivedHandle<Ground, GroundImpl, MobilizedBody> {
public:
    Ground();
    Ground& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
};

/// 1-6 mobilities. TODO: this will be an abstract class with virtual methods
/// that a user's derived class can implement to define a custom mobilizer.
/// TODO: not implemented yet.
class SimTK_SIMBODY_EXPORT MobilizedBody::Custom : public PIMPLDerivedHandle<Custom, CustomImpl, MobilizedBody> {
public:
    Custom(int nMobilities, int nCoordinates);

    // Get calculations through Stage::Instance from State.
    virtual void calcTransform(const State&, const Vector& q, 
                               Transform& X_FM) const = 0;
    //TODO: should H be a nuX2 Matrix_<Vec3> instead? or Vector_<SpatialVec>?
    //      or nuX6 Matrix?
    virtual void calcTransitionMatrix(const State& s, 
                               Vector_<SpatialRow>& H_FM) const = 0;
    virtual void calcTransitionMatrixTimeDerivative(const State& s,  
                               Vector_<SpatialRow>& H_FM_Dot) const = 0;

    // get q and calculations through Stage::Position from State if needed
    virtual void calcQDot(const State&, const Vector& u, Vector& qdot) const {
        qdot = u; //TODO: only if sizes match
    }
    // get q,u and calculations through Stage::Dynamics from State if needed
    virtual void calcQDotDot(const State&, const Vector& udot, Vector& qdotdot) const {
        qdotdot = udot; //TODO: only if sizes match
    }
protected:
    // Utilities for use by Custom mobilized body implementation.

    // Be sure to call this whenever you make a change to any data contained
    // in a concrete Custom MobilizedBody class. This method ensures that the
    // containing matter subsystem will have its topology invalidated so that
    // a subsequent call to realizeTopology() will recalculate the topology 
    // cache. A good rule of thumb is that any method you provide which is
    // non-const should start by calling invalidateTopologyCache().
    void invalidateTopologyCache() const;
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_H_



