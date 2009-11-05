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
 * Portions copyright (c) 2007-9 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors: Paul Mitiguy, Peter Eastman                                  *
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
 * This defines the MobilizedBody class, which associates a new body (the 
 * "child", "outboard", or "successor" body) with a Mobilizer and a reference 
 * frame on an existing body (the "parent", "inboard", or "predecessor" body)
 * that is already part of a MatterSubsystem.
 *
 * MobilizedBody is an abstract base class handle, with concrete classes 
 * defined for each kind of mobilizer. There are a set of built-in mobilizers
 * and a generic "Custom" mobilizer (an actual abstract base class) from
 * which advanced users may derive their own mobilizers.
 *
 * A Mobilizer may be associated with a Motion object which defines how
 * it is to move; otherwise its motion is calculated as a result of the 
 * application of forces (either directly applied or resulting from constraint
 * forces generated to satisfy restrictions imposed by Constraint objects).
 */

#include "SimTKmath.h"
#include "simbody/internal/common.h"
#include "simbody/internal/Body.h"
#include "simbody/internal/Motion.h"

#include <cassert>
#include <vector>

namespace SimTK {

class SimbodyMatterSubsystem;
class Motion;
class MobilizedBody;
class MobilizedBodyImpl;

// We only want the template instantiation to occur once. This symbol is 
// defined in the SimTK core compilation unit that instantiates the mobilized 
// body class but should not be defined any other time.
#ifndef SimTK_SIMBODY_DEFINING_MOBILIZED_BODY
    extern template class PIMPLHandle<MobilizedBody, MobilizedBodyImpl, true>;
#endif

/**
 * This is the base class for all MobilizedBody classes, just a handle for the 
 * underlying hidden implementation. Each built-in MobilizedBody type is a local 
 * subclass within MobilizedBody, so the built-ins have names like 
 * MobilizedBody::Pin. All concrete MobilizedBodies, including the built-ins, 
 * are derived from MobilizedBody.
 *
 * There are three sets of methods used for obtaining MobilizedBody-specific 
 * data from the containing System's State. These are:
 *    - State Access
 *    - Basic Operators
 *    - High Level Operators
 *
 * <em>State Access</em> methods simply extract already-calculated data from the 
 * State or State Cache, or set State values. They involve no additional 
 * computation, have names beginning with "get" and "upd" (update) and return 
 * references to the requested quantities rather than calculated values. We 
 * divide these into routines which deal with bodies and routines which deal 
 * with mobilizers and mobilities.
 *
 * <em>Basic Operators</em> use State Access methods to compute basic quantities
 * which cannot be precomputed, such as the velocity of an arbitrary point, 
 * using an inline combination of basic floating point operations which can be 
 * reliably determined at compile time. These have names beginning with "find" 
 * or a more specific verb, as a reminder that they do not require a great deal 
 * of computation. 
 *
 * <em>High Level Operators</em> combine responses and basic operators with 
 * run-time tests to calculate more complex quantities, with more complicated 
 * implementations that can exploit special cases at run time. These begin with 
 * "calc" (calculate) as a reminder that they may involve substantial run time 
 * computation.
 *
 * There is also a set of methods used for construction, and miscellaneous 
 * utilities. These methods are primarly intended for use by concrete 
 * MobilizedBody classes and are not generally used by end users.
 *
 * In the API below, we'll refer to the current ("this") MobilizedBody as "body
 * B". It is the "object" or "main" body with which we are concerned. Often 
 * there will be another body mentioned in the argument list as a target for 
 * some conversion. That "another" body will be called "body A". The Ground 
 * body is abbreviated "G".
 *
 * We use OF to mean "the origin of frame F", CB is "the mass center of body 
 * B". R_AF is the rotation matrix giving frame F's orientation in frame A, 
 * such that a vector v expressed in F is reexpressed in A by v_A = R_AF * v_F.
 * X_AF is the spatial transform giving frame F's origin location and 
 * orientation in frame A, such that a point P whose location is measured 
 * from F's origin OF and expressed in F by position vector p_FP (or more 
 * explicitly p_OF_P) is remeasured from frame A's origin OA and reexpressed 
 * in A via p_AP = X_AF * p_FP, where p_AP==p_OA_P. 
 */

class SimTK_SIMBODY_EXPORT MobilizedBody : public PIMPLHandle<MobilizedBody, MobilizedBodyImpl, true> {
public:

    /// Constructors can take an argument of this type to indicate that the 
    /// mobilizer is being defined in the reverse direction, meaning from 
    /// child to parent. That means that the mobilizer coordinates and speeds
    /// will be defined as though the tree had been built in the opposite
    /// direction. This is a topological setting and can't be changed dynamically.
    enum Direction {
        Forward = 0,
        Reverse = 1
    };

    /// The default behavior of this mobilizer will normally be determined
    /// by whether you provide a Motion object for it. However, you can override
    /// that afterwards.
    MobilizedBody& setDefaultMotionType(Motion::Level, Motion::Method=Motion::Prescribed);

    /// This is an Instance stage setting.
    void setMotionType(State&, Motion::Level, Motion::Method=Motion::Prescribed) const;

    bool isAccelerationAlwaysZero(const State&) const;
    bool isVelocityAlwaysZero(const State&) const;


        //////////////////////////
        // STATE ACCESS METHODS //
        //////////////////////////

    /// @name State Access - Bodies
    /// These methods extract already-computed information from the State or 
    /// State cache, or set values in the State.
    //@{

    /// Extract from the state cache the already-calculated spatial 
    /// configuration X_GB of body B's body frame, measured with respect to the 
    /// Ground frame and expressed in the Ground frame. That is, we return the 
    /// location of the body frame's origin, and the orientation of its x, y, 
    /// and z axes, as the Transform X_GB. This notation is intended to convey 
    /// unambiguously the sense of this transform, which is as follows: if you 
    /// have a station (body fixed point) S on body B, represented by position 
    /// vector p_BS (a.k.a. p_OB_S) from the origin OB of B to the point S and 
    /// expressed in the B frame, then p_GS=X_GB*p_BS where p_GS (== p_OG_S) is 
    /// the position vector from the Ground origin OG to the point in space 
    /// currently coincident with S and expressed in the Ground frame. The 
    /// inverse transformation is obtained using the "~" operator where 
    /// ~X_GB=X_BG, so that p_BS = ~X_GB*p_GS. This response is available at 
    /// Position stage.
    const Transform& getBodyTransform(const State&) const; // X_GB

    /// Extract from the state cache the already-calculated spatial orientation
    /// R_GB of body B's body frame x, y, and z axes expressed in the Ground 
    /// frame, as the Rotation matrix R_GB. The sense of this rotation matrix 
    /// is such that if you have a vector v fixed on body B, represented by the
    /// vector v_B expressed in the B frame, then v_G=R_GB*v_B where v_G is the 
    /// same vector but re-expressed in the Ground frame. The inverse 
    /// transformation is obtained using the "~" operator where ~R_GB=R_BG, so 
    /// that v_B = ~R_GB*v_G. This response is available at Position stage.
    const Rotation& getBodyRotation(const State& s) const {
        return getBodyTransform(s).R();
    }
    /// Extract from the state cache the already-calculated spatial location
    /// of body B's body frame origin OB, measured from the Ground origin OG and
    /// expressed in the Ground frame, as the position vector p_GB (== p_OG_OB).
    /// This response is available at Position stage.
    const Vec3& getBodyOriginLocation(const State& s) const {
        return getBodyTransform(s).p();
    }

    /// At stage Position or higher, return the cross-mobilizer transform X_FM,
    /// the body's inboard mobilizer frame M measured and expressed in
    /// the parent body's corresponding outboard frame F.
    const Transform& getMobilizerTransform(const State&) const; // X_FM

    /// Extract from the state cache the already-calculated spatial velocity 
    /// V_GB of this body's reference frame B, measured with respect to the 
    /// Ground frame and expressed in the Ground frame. That is, we return the 
    /// linear velocity v_GB of the body frame's origin in G, and the body's 
    /// angular velocity w_GB as the spatial velocity vector V_GB = {w_GB, v_GB}.
    /// This response is available at Velocity stage.
    const SpatialVec& getBodyVelocity(const State&) const;          // V_GB

    /// Extract from the state cache the already-calculated inertial angular
    /// velocity vector w_GB of this body B, measured with respect to the Ground frame
    /// and expressed in the Ground frame. This response is available at Velocity stage.
    const Vec3& getBodyAngularVelocity(const State& s) const {      // w_GB
        return getBodyVelocity(s)[0]; 
    }
    /// Extract from the state cache the already-calculated inertial linear
    /// velocity vector v_GB (more explicitly, v_G_OB) of this body B's origin
    /// point OB, measured with respect to the Ground frame and expressed in
    /// the Ground frame. This response
    /// is available at Velocity stage.
    const Vec3& getBodyOriginVelocity(const State& s) const {       // v_GB
        return getBodyVelocity(s)[1];
    }

    /// At stage Velocity or higher, return the cross-mobilizer velocity
    /// V_FM, the relative velocity of this body's "moving" mobilizer
    /// frame M in the parent body's corresponding "fixed" frame F, 
    /// measured and expressed in F. Note that this isn't the usual 
    /// spatial velocity since it isn't expressed in G.
    const SpatialVec& getMobilizerVelocity(const State&) const;     // V_FM

    /// Extract from the state cache the already-calculated spatial acceleration A_GB of
    /// this body's reference frame B, measured with respect to the Ground frame and expressed
    /// in the Ground frame. That is, we return the linear acceleration a_GB of the body
    /// frame's origin in G, and the body's angular acceleration b_GB as the spatial acceleration
    /// vector A_GB = {b_GB, a_GB}. This response is available at Acceleration stage.
    const SpatialVec& getBodyAcceleration(const State& s) const;    // A_GB

    /// Extract from the state cache the already-calculated inertial angular
    /// acceleration vector b_GB of this body B, measured with respect to the Ground frame
    /// and expressed in the Ground frame. This response is available at Acceleration stage.
    const Vec3& getBodyAngularAcceleration(const State& s) const {  // b_GB
        return getBodyAcceleration(s)[0]; 
    }

    /// Extract from the state cache the already-calculated inertial linear
    /// acceleration vector a_GB (more explicitly, a_G_OB) of this body B's origin
    /// point OB, measured with respect to the Ground frame and expressed in the
    /// Ground frame. This response is available at Acceleration stage.
    const Vec3& getBodyOriginAcceleration(const State& s) const {   // a_GB
        return getBodyAcceleration(s)[1];
    }

    /// TODO: Not implemented yet -- any volunteers? At stage Acceleration, return the cross-mobilizer
    /// acceleration A_FM, the relative acceleration of body B's "moving" mobilizer
    /// frame M in the parent body's corresponding "fixed" frame F, 
    /// measured and expressed in F. Note that this isn't the usual 
    /// spatial acceleration since it isn't expressed in G.
    const SpatialVec& getMobilizerAcceleration(const State&) const { // A_FM
        SimTK_ASSERT_ALWAYS(!"unimplemented method", 
            "MobilizedBody::getMobilizerAcceleration() is not yet implemented -- any volunteers?");
        return *(new SpatialVec());
    }

    /// Return a reference to this body's mass properties in the State cache.
    /// The State must have been realized to Stage::Instance or higher.
    const MassProperties& getBodyMassProperties(const State&) const;

    /// Return the mass of this body. The State must have been realized to Stage::Instance.
    Real getBodyMass(const State& s) const {
        return getBodyMassProperties(s).getMass();
    }

    /// Return this body's center of mass station (i.e., the vector fixed in the body,
    /// going from body origin to body mass center, expressed in the body frame.)
    /// The State must have been realized to Stage::Instance or higher.
    const Vec3& getBodyMassCenterStation(const State& s) const {
        return getBodyMassProperties(s).getMassCenter();
    }

    /// Return a reference to this body's inertia matrix in the State cache, taken
    /// about the body origin and expressed in the body frame.
    /// The State must have been realized to Stage::Instance or higher.
    const Inertia& getBodyInertiaAboutBodyOrigin(const State& s) const {
        return getBodyMassProperties(s).getInertia();
    }

    /// Return a reference to this mobilizer's frame F fixed on the parent body P, as the fixed Transform
    /// from P's body frame to the frame F fixed to P. If this frame is changeable, the result comes
    /// from the State cache, otherwise it is from the MobilizedBody object itself.
    /// The State must have been realized to Stage::Instance or higher.
    const Transform& getInboardFrame (const State&) const;  // X_PF
    /// Return a reference to this MobilizedBody's mobilizer frame M, as the fixed Transform
    /// from this body B's frame to the frame M fixed on B. If this frame is changeable, the result comes
    /// from the State cache, otherwise it is from the MobilizedBody object itself.
    /// The State must have been realized to Stage::Instance or higher.
    const Transform& getOutboardFrame(const State&) const;  // X_BM

    /// TODO: not implemented yet. Set the location and orientation of the inboard (parent) mobilizer
    /// frame F, fixed to this mobilizer's parent body P.
    void setInboardFrame (State&, const Transform& X_PF) const;
    /// TODO: not implemented yet. Set the location and orientation of the outboard mobilizer
    /// frame M, fixed to this body B.
    void setOutboardFrame(State&, const Transform& X_BM) const;

    // End of State Access - Bodies
    //@}

    /// @name State Access - Mobilizer generalized coordinates q and speeds u
    /// These methods extract q- or u-related information from the State or State cache, or set
    /// q or u values in the State.
    //@{
    /// Return the number of generalized coordinates q currently in use by this mobilizer.
    /// State must have been realized to Stage::Model.
    int getNumQ(const State&) const;
    /// Return the number of generalized speeds u currently in use by this mobilizer.
    /// State must have been realized to Stage::Model.
    int getNumU(const State&) const;

    /// Return one of the generalized coordinates q from this mobilizer's partition of the matter
    /// subsystem's full q vector in the State. The particular coordinate is selected using the \p which
    /// parameter, numbering from zero to getNumQ()-1.
    Real getOneQ(const State&, int which) const;

    /// Return one of the generalized speeds u from this mobilizer's partition of the matter
    /// subsystem's full u vector in the State. The particular coordinate is selected using the \p which
    /// parameter, numbering from zero to getNumU()-1.
    Real getOneU(const State&, int which) const;

    /// Return as a Vector of length getNumQ() all the generalized coordinates q
    /// currently in use by this mobilizer, from this mobilizer's partion in 
    /// the matter subsystem's full q vector in the State.
    Vector getQAsVector(const State&) const;
    /// Return as a Vector of length getNumU() all the generalized speeds u
    /// currently in use by this mobilizer, from this mobilizer's partion in 
    /// the matter subsystem's full u vector in the State.
    Vector getUAsVector(const State&) const;

    /// Return one of the generalized coordinate derivatives qdot from this mobilizer's partition of the matter
    /// subsystem's full qdot vector in the State cache. The particular coordinate is selected using the \p which
    /// parameter, numbering from zero to getNumQ()-1.
    Real getOneQDot   (const State&, int which) const;
    /// Return as a Vector of length getNumQ() all the generalized coordinate derivatives qdot
    /// currently in use by this mobilizer, from this mobilizer's partion in 
    /// the matter subsystem's full qdot vector in the State cache.
    Vector getQDotAsVector(const State&) const;

    /// Return one of the generalized accelerations udot from this mobilizer's 
    /// partition of the matter subsystem's full udot vector in the State 
    /// cache. The particular coordinate is selected using the \p which
    /// parameter, numbering from zero to getNumU()-1.
    Real getOneUDot(const State&, int which) const;
    /// Return one of the generalized coordinate second derivatives qdotdot 
    /// from this mobilizer's partition of the matter subsystem's full 
    /// qdotdot vector in the State cache. The particular coordinate is 
    /// selected using the \p which parameter, numbering from zero to 
    /// getNumQ()-1.
    Real getOneQDotDot(const State&, int which) const;
    /// Return as a Vector of length getNumU() all the generalized 
    /// accelerations udot currently in use by this mobilizer, from this 
    /// mobilizer's partion in the matter subsystem's full udot vector in 
    /// the State cache.
    Vector getUDotAsVector(const State&) const;
    /// Return as a Vector of length getNumQ() all the generalized coordinate 
    /// second derivatives qdotdot currently in use by this mobilizer, from 
    /// this mobilizer's partion in the matter subsystem's full qdotdot vector 
    /// in the State cache.
    Vector getQDotDotAsVector(const State&) const;

    /// Return the tau forces resulting from known (prescribed) acceleration, 
    /// corresponding to each of this mobilizer's mobilities, as a Vector 
    /// of length getNumU().
    ///
    /// If this mobilizer has known accelerations (UDots) due to an active
    /// Motion object, the set of generalized forces tau that must be added
    /// in order to produce those accelerations is calculated at Acceleration
    /// stage. There is one scalar tau per mobility and they can be returned
    /// individually or as a Vector. The return value is zero if the
    /// accelerations are free.
    Vector getTauAsVector(const State&) const;
    /// Return one of the tau forces resulting from known (prescribed) 
    /// acceleration, corresponding to one of this mobilizer's mobilities 
    /// as selected here using the \p which parameter, numbered from 
    /// zero to getNumU()-1.
    /// @see getTauAsVector() for more information
    Real getOneTau(const State&, MobilizerUIndex which) const;


    /// Set one of the generalized coordinates q to value \p v, in this mobilizer's partition of the matter
    /// subsystem's full q vector in the State. The particular coordinate is selected using the \p which
    /// parameter, numbering from zero to getNumQ()-1.
    void setOneQ(State&, int which, Real v) const;
    /// Set one of the generalized speeds u to value \p v, in this mobilizer's partition of the matter
    /// subsystem's full u vector in the State. The particular coordinate is selected using the \p which
    /// parameter, numbering from zero to getNumU()-1.
    void setOneU(State&, int which, Real v) const;

    /// Set all of the generalized coordinates q to value \p v (a Vector of length getNumQ()),
    /// in this mobilizer's partition of the matter subsystem's full q vector in the State.
    void setQFromVector(State& s, const Vector& v) const;
    /// Set all of the generalized speeds u to value \p v (a Vector of length getNumU()),
    /// in this mobilizer's partition of the matter subsystem's full u vector in the State.
    void setUFromVector(State& s, const Vector& v) const;

    /// Adjust this mobilizer's q's to best approximate the supplied Transform
    /// which requests a particular relative orientation and translation between
    /// the "fixed" and "moving" frames connected by this mobilizer.
    ///
    /// This set of methods sets the generalized coordinates, or speeds (state
    /// variables) for just the mobilizer associated with this MobilizedBody
    /// (ignoring all other mobilizers and constraints), without requiring knowledge
    /// of the meanings of the individual state variables. The idea here
    /// is to provide a physically-meaningful quantity relating the 
    /// mobilizer's inboard and outboard frames, and then ask the mobilizer
    /// to set its state variables to reproduce that quantity to the
    /// extent it can.
    ///
    /// These routines can be called in Stage::Model, however the routines
    /// may consult the current values of the state variables in some cases,
    /// so you must make sure they have been set to reasonable, or at least
    /// innocuous values (zero will work). In no circumstance will any of
    /// these routines look at any state variables which belong to another
    /// mobilizer; they are limited to working locally with one mobilizer.
    ///
    /// Routines which specify only translation (linear velocity) may use
    /// rotational coordinates to help satisfy the translation requirement.
    /// An alternate "Only" method is available to forbid modification of 
    /// purely rotational coordinates in that case. When a mobilizer uses
    /// state variables which have combined rotational and translational
    /// character (e.g. a screw joint) consult the documentation for the
    /// mobilizer to find out how it responds to these routines.
    ///
    /// There is no guarantee that the desired physical quantity will be
    /// achieved by these routines; you can check on return if you're
    /// worried. Individual mobilizers make specific promises about what
    /// they will do; consult the documentation. These routines do not
    /// throw exceptions even for absurd requests like specifying a
    /// rotation for a sliding mobilizer. Nothing happens if
    /// there are no mobilities here, i.e. Ground or a Weld mobilizer.

    void setQToFitTransform      (State&, const Transform& X_FM) const;
    /// Adjust this mobilizer's q's to best approximate the supplied Rotation matrix
    /// which requests a particular relative orientation between the "fixed"
    /// and "moving" frames connected by this mobilizer.
    /// @see setQToFitTransform()
    void setQToFitRotation       (State&, const Rotation&  R_FM) const;
    /// Adjust this mobilizer's q's to best approximate the supplied position vector
    /// which requests a particular offset between the origins of the "fixed"
    /// and "moving" frames connected by this mobilizer, with <em>any</em> q's (rotational
    /// or translational) being modified if doing so helps satisfy the request.
    /// @see setQToFitTransform()
    void setQToFitTranslation    (State&, const Vec3&      p_FM) const;

    /// Adjust this mobilizer's u's (generalized speeds) to best approximate
    /// the supplied spatial velocity \p V_FM which requests the relative angular
    /// and linear velocity between the "fixed" and "moving" frames connected by
    /// this mobilizer. Routines which affect generalized speeds u depend on the generalized
    /// coordinates q already having been set; they never change these coordinates.
    /// @see setQToFitTransform()
    void setUToFitVelocity          (State&, const SpatialVec& V_FM) const;
    /// Adjust this mobilizer's u's (generalized speeds) to best approximate
    /// the supplied angular velocity \p w_FM which requests a particular relative angular
    /// between the "fixed" and "moving" frames connected by
    /// this mobilizer.
    /// @see setQToFitTransform()
    /// @see setUToFitVelocity()
    void setUToFitAngularVelocity   (State&, const Vec3&       w_FM) const;
    /// Adjust <em>any</em> of this mobilizer's u's (generalized speeds) to best approximate
    /// the supplied linear velocity \p v_FM which requests a particular velocity for
    /// the "moving" frame M origin in the "fixed" frame F on the parent where these
    /// are the frames connected by this mobilizer.
    /// @see setQToFitTransform()
    /// @see setUToFitVelocity()
    void setUToFitLinearVelocity    (State&, const Vec3&       v_FM) const;

    /// Expert use only: obtain a column of the hinge matrix H corresponding to
    /// one of this mobilizer's mobilities.
    SpatialVec getHCol(const State& s, UIndex ux) const;

    // End of State Access Methods.
    //@} 

        /////////////////////
        // BASIC OPERATORS //
        /////////////////////

    /// @name Basic Operators

    /// These methods use state variables and Response methods to compute basic quantities
    /// which cannot be precomputed, but which can be implemented with an inline combination
    /// of basic floating point operations which can be reliably determined at compile time.
    /// The method names and descriptions use the following terms:
    ///     - Body or ThisBody: the Body B associated with the current MobilizedBody. ThisBody
    ///       is implied when no other Body is mentioned.
    ///     - Ground: the "MobilizedBody" G representing the Ground reference frame which
    ///       never moves.
    ///     - AnotherBody: the Body A being referenced, which in general is neither ThisBody
    ///       nor Ground.
    ///     - Station: a point S fixed on ThisBody B, located by a position vector
    ///       p_BS (or more explicitly, p_OB_S) from the B-frame origin OB to the point S,
    ///       expressed in the B-frame coordinate system.
    ///     - Vector: a vector v fixed on ThisBody B, given by a vector v_B expressed in
    ///       the B-frame coordinate system.
    ///     - Origin: the Station located at (0,0,0) in ThisBody frame B, that is, body B's origin
    ///       point.
    ///     - MassCenter: the Station on ThisBody B which is the center of mass for B.
    ///     - GroundPoint, GroundVector: a Point P or Vector v on the Ground "Body" G. These
    ///       are measured and expressed in the Ground frame, as p_GP or v_G.
    ///     - AnotherBodyStation, AnotherBodyVector, etc.: a Station S or Vector v on AnotherBody A.
    ///       These are measured and expressed in the A frame, as p_AS or v_A. 

    //@{

    /// Return X_AB, the spatial transform giving this body B's frame in body A's frame.
    /// Cost is 63 flops. If you know that one of the bodies is Ground, use the 0-cost
    /// Response getBodyTransform() instead of this operators.
    /// This operator is available in Position stage.
    /// @see getBodyTransform()
    Transform findBodyTransformInAnotherBody(const State& s, 
                                             const MobilizedBody& inBodyA) const
    {
        const Transform& X_GA = inBodyA.getBodyTransform(s);
        const Transform& X_GB = this->getBodyTransform(s);

        return ~X_GA*X_GB; // X_AB=X_AG*X_GB
    }

    /// Return R_AB, the rotation matrix giving this body B's axes in body A's frame.
    /// Cost is 45 flops. If you know that one of the bodies is Ground, use the 0-cost
    /// response getBodyRotation() instead of this operators.
    /// This operator is available in Position stage.
    /// @see getBodyRotation()
    Rotation findBodyRotationInAnotherBody(const State& s, 
                                            const MobilizedBody& inBodyA) const
    {
        const Rotation& R_GA = inBodyA.getBodyRotation(s);
        const Rotation& R_GB = this->getBodyRotation(s);

        return ~R_GA*R_GB; // R_AB=R_AG*R_GB
    }

    /// Return the station on another body A (that is, a point measured and expressed in A) that is 
    /// currently coincident in space with the origin OB of this body B. Cost is 18 flops.
    /// This operator is available at Position stage. Note: "findBodyOriginLocationInGround" 
    /// doesn't exist because it would be the same as the Response getBodyOriginLocation().
    /// @see getBodyOriginLocation()
    Vec3 findBodyOriginLocationInAnotherBody(const State& s, const MobilizedBody& toBodyA) const {
        return toBodyA.findStationAtGroundPoint(s,getBodyOriginLocation(s));
    }

    /// Return the angular and linear velocity of body B's frame in body A's frame, expressed in body A,
    /// and arranged as a SpatialVec.
    /// Cost is 51 flops. If you know inBodyA is Ground, don't use this routine; use the response
    /// method getBodyVelocity() which is free.
    /// This operator is available in Velocity stage.
    SpatialVec findBodyVelocityInAnotherBody(const State& s,
                                             const MobilizedBody& inBodyA) const
    {
        const SpatialVec& V_GB   = this->getBodyVelocity(s);
        const SpatialVec& V_GA   = inBodyA.getBodyVelocity(s);
        const Vec3        w_AB_G = V_GB[0]-V_GA[0]; // angular velocity of B in A, exp in G   ( 3 flops)

        // Angular velocity was easy, but for linear velocity we have to add in an wXr term.
        const Transform&  X_GB       = getBodyTransform(s);
        const Transform&  X_GA       = inBodyA.getBodyTransform(s);
        const Vec3        p_AB_G     = X_GB.p() - X_GA.p(); // vector from OA to OB, exp in G ( 3 flops)
        const Vec3        p_AB_G_dot = V_GB[1]  - V_GA[1];  // d/dt p taken in G              ( 3 flops)

        const Vec3 v_AB_G = p_AB_G_dot - V_GA[0] % p_AB_G;  // d/dt p taken in A, exp in G    (12 flops)

        // We're done, but the answer is expressed in Ground. Reexpress in A and return.
        return ~X_GA.R()*SpatialVec(w_AB_G, v_AB_G);        //                                (30 flops)
    }

    /// Return the angular velocity w_AB of body B's frame in body A's frame, expressed in body A.
    /// Cost is 18 flops. If you know inBodyA is Ground, don't use this routine; use the response
    /// method getBodyAngularVelocity() which is free.
    /// This operator is available in Velocity stage.
    Vec3 findBodyAngularVelocityInAnotherBody(const State& s,
                                       const MobilizedBody& inBodyA) const 
    {
        const Vec3& w_GB   = getBodyAngularVelocity(s);
        const Vec3& w_GA   = inBodyA.getBodyAngularVelocity(s);
        const Vec3  w_AB_G = w_GB-w_GA; // angular velocity of B in A, exp in G ( 3 flops)

        // Now reexpress in A.
        return inBodyA.expressGroundVectorInBodyFrame(s, w_AB_G); //            (15 flops)
    }

    /// Return the velocity of body B's origin point in body A's frame, expressed in body A.
    /// Cost is 51 flops. If you know inBodyA is Ground, don't use this routine; use the response
    /// method getBodyOriginVelocity() which is free.
    /// This operator is available in Velocity stage.
    Vec3 findBodyOriginVelocityInAnotherBody(const State& s,
                                      const MobilizedBody& inBodyA) const
    {
        // Doesn't save much to special case this one.
        return findBodyVelocityInAnotherBody(s,inBodyA)[1];
    }

    /// Return the angular and linear acceleration of body B's frame in body A's frame, expressed in body A,
    /// and arranged as a SpatialVec. Cost is 105 flops. If you know that inBodyA is Ground, don't
    /// use this operator; instead, use the response method getBodyAcceleration() which is free.
    /// This operator is available in Acceleration stage.
    SpatialVec findBodyAccelerationInAnotherBody(const State& s,
                                                 const MobilizedBody& inBodyA) const
    {
        const Vec3&       p_GB = this->getBodyOriginLocation(s);
        const Transform&  X_GA = inBodyA.getBodyTransform(s);
        const SpatialVec& V_GB = this->getBodyVelocity(s);
        const SpatialVec& V_GA = inBodyA.getBodyVelocity(s);
        const SpatialVec& A_GB = this->getBodyAcceleration(s);
        const SpatialVec& A_GA = inBodyA.getBodyAcceleration(s);
        const Vec3&       p_GA = X_GA.p();
        const Vec3&       w_GA = V_GA[0];
        const Vec3&       w_GB = V_GB[0];
        const Vec3&       b_GA = A_GA[0];
        const Vec3&       b_GB = A_GB[0];

        const Vec3 p_AB_G        = p_GB     - p_GA;         // vector from OA to OB, in G   ( 3 flops)
        const Vec3 p_AB_G_dot    = V_GB[1]  - V_GA[1];      // d/dt p taken in G            ( 3 flops)
        const Vec3 p_AB_G_dotdot = A_GB[1]  - A_GA[1];      // d^2/dt^2 taken in G          ( 3 flops)

        const Vec3 w_AB_G     = w_GB - w_GA;                // relative ang. vel. of B in A, exp. in G (3 flops)
        const Vec3 v_AB_G     = p_AB_G_dot - w_GA % p_AB_G; // d/dt p taken in A, exp in G  (12 flops)

        const Vec3 w_AB_G_dot = b_GB - b_GA;                // d/dt of w_AB_G taken in G    ( 3 flops)
        const Vec3 v_AB_G_dot = p_AB_G_dotdot - (b_GA % p_AB_G + w_GA % p_AB_G_dot); // d/dt v_AB_G taken in G
                                                                                     //     (24 flops)

        // We have the derivative in G; change it to derivative in A by adding in contribution caused
        // by motion of G in A, that is w_AG X w_AB_G. (Note that w_AG=-w_GA.)
        const Vec3 b_AB_G = w_AB_G_dot - w_GA % w_AB_G; // ang. accel. of B in A            (12 flops)
        const Vec3 a_AB_G = v_AB_G_dot - w_GA % v_AB_G; // taken in A, exp. in G            (12 flops)

        return ~X_GA.R() * SpatialVec(b_AB_G, a_AB_G); // taken in A, expressed in A        (30 flops)
    }

    /// Return the angular acceleration of body B's frame in body A's frame, expressed in body A.
    /// Cost is 33 flops. If you know \p inBodyA is Ground, don't use this operator; instead use
    /// the response method getBodyAngularAccleration() which is free. This operator is available
    /// at AccelerationStage.
    Vec3 findBodyAngularAccelerationInAnotherBody(const State& s,
                                                  const MobilizedBody& inBodyA) const
    {
        const Rotation& R_GA = inBodyA.getBodyRotation(s);
        const Vec3&     w_GA = inBodyA.getBodyAngularVelocity(s);
        const Vec3&     w_GB = this->getBodyAngularVelocity(s);
        const Vec3&     b_GA = inBodyA.getBodyAngularAcceleration(s);
        const Vec3&     b_GB = this->getBodyAngularAcceleration(s);

        const Vec3 w_AB_G     = w_GB - w_GA;                // relative ang. vel. of B in A, exp. in G (3 flops)
        const Vec3 w_AB_G_dot = b_GB - b_GA;                // d/dt of w_AB_G taken in G    ( 3 flops)

        // We have the derivative in G; change it to derivative in A by adding in contribution caused
        // by motion of G in A, that is w_AG X w_AB_G. (Note that w_AG=-w_GA.)
        const Vec3 b_AB_G = w_AB_G_dot - w_GA % w_AB_G; // ang. accel. of B in A            (12 flops)

        return ~R_GA * b_AB_G; // taken in A, expressed in A                                (15 flops)
    }

    /// Return the acceleration of body B's origin point in body A's frame, expressed in body A.
    /// Cost is 105 flops. If you know that inBodyA is Ground, don't
    /// use this operator; instead, use the response method getBodyOriginAcceleration() which is free.
    /// This operator is available in Acceleration stage.
    Vec3 findBodyOriginAccelerationInAnotherBody(const State& s, 
                                          const MobilizedBody& inBodyA) const
    {
        // Not much to be saved by trying to optimize this since the linear part
        // is the most expensive to calculate.
        return findBodyAccelerationInAnotherBody(s,inBodyA)[1];
    }

    /// Return the Cartesian (ground) location that is currently coincident with
    /// a station (point) S fixed on body B. That is, we return locationInG = X_GB * stationOnB
    /// which means the result is measured from the Ground origin and expressed in Ground.
    /// In more precise notation, we're calculating p_GS = X_GB * p_BS for position vectors
    /// p_GS and p_BS. Cost is 18 flops. This operator is available at Position stage.
    Vec3 findStationLocationInGround(const State& s, const Vec3& stationOnB) const {
        return getBodyTransform(s) * stationOnB;
    }


    /// Given a station S on this body B, return the location on another body A which is at
    /// the same location in space. That is, we return locationOnA = X_AB * locationOnB,
    /// which means the result is measured from the body A origin and expressed in body A. In
    /// more precise notation, we're calculating p_AS = X_AB * p_BS, which we actually
    /// calculate as p_AS = X_AG*(X_GB*p_BS). Cost is 36 flops.
    /// Note: if you know that one of the bodies is Ground, use one of the routines above
    /// which is specialized for Ground to avoid half the work.
    /// This operator is available at Position stage or higher.
    Vec3 findStationLocationInAnotherBody(const State& s, const Vec3& stationOnB, 
                               const MobilizedBody& toBodyA) const
    {
        return toBodyA.findStationAtGroundPoint(s, findStationLocationInGround(s,stationOnB));
    }

    /// Given a station fixed on body B, return its inertial (Cartesian) velocity,
    /// that is, its velocity relative to the Ground frame, expressed in the
    /// Ground frame. Cost is 27 flops. This operator is available at Velocity stage.
    Vec3 findStationVelocityInGround(const State& s, const Vec3& stationOnB) const {
        const Vec3& w = getBodyAngularVelocity(s); // in G
        const Vec3& v = getBodyOriginVelocity(s);  // in G
        const Vec3  r = expressVectorInGroundFrame(s,stationOnB); // 15 flops
        return v + w % r;                                         // 12 flops
    }


    /// Return the velocity of a station S fixed on body B, in body A's frame, expressed in body A.
    /// Cost is 93 flops. If you know \p inBodyA is Ground, don't use this operator; instead use
    /// the operator findStationVelocityInGround() which is much cheaper.
    /// This operator is available in Velocity stage.
    Vec3 findStationVelocityInAnotherBody(const State& s, 
                                          const Vec3&          stationOnBodyB, // p_BS
                                          const MobilizedBody& inBodyA) const
    {
        const SpatialVec V_AB   = findBodyVelocityInAnotherBody(s,inBodyA); // (51 flops)
         // OB->S rexpressed in A but not shifted to OA
        const Vec3       p_BS_A = expressVectorInAnotherBodyFrame(s, stationOnBodyB, inBodyA);
                                                                            // (30 flops)
        return V_AB[1] + (V_AB[0] % p_BS_A);                                // (12 flops)
    }

      
    /// Given a station fixed on body B, return its inertial (Cartesian) acceleration,
    /// that is, its acceleration relative to the ground frame, expressed in the
    /// ground frame. Cost is 48 flops. This operator is available at Acceleration stage.
    Vec3 findStationAccelerationInGround(const State& s, const Vec3& stationOnB) const {
        const Vec3& w = getBodyAngularVelocity(s);     // in G
        const Vec3& b = getBodyAngularAcceleration(s); // in G
        const Vec3& a = getBodyOriginAcceleration(s);  // in G

        const Vec3  r = expressVectorInGroundFrame(s,stationOnB); // 15 flops
        return a + b % r + w % (w % r);                           // 33 flops
    }

    /// Return the acceleration of a station S fixed on body B, in another body A's frame, expressed in body A.
    /// Cost is 186 flops.  If you know that \p inBodyA is Ground, don't
    /// use this operator; instead, use the operator findStationAccelerationInGround() which is 
    /// much cheaper. This operator is available in Acceleration stage.
    Vec3 findStationAccelerationInAnotherBody(const State& s,
                                              const Vec3&          stationOnBodyB, 
                                              const MobilizedBody& inBodyA) const 
    {
        const Vec3       w_AB = findBodyAngularVelocityInAnotherBody(s,inBodyA);  // ( 18 flops)
        const SpatialVec A_AB = findBodyAccelerationInAnotherBody(s,inBodyA);     // (105 flops)
         // OB->S rexpressed in A but not shifted to OA
        const Vec3       p_BS_A = expressVectorInAnotherBodyFrame(s, stationOnBodyB, inBodyA);
                                                                                  // ( 30 flops)

        return A_AB[1] + (A_AB[0] % p_BS_A) + w_AB % (w_AB % p_BS_A);             // ( 33 flops)
    }

    /// It is cheaper to calculate a station's ground location and velocity together
    /// than to do them separately. Here we can return them both in 30 flops, vs. 45 to
    /// do them in two calls. This operator is available at Velocity stage.
    void findStationLocationAndVelocityInGround(const State& s, const Vec3& locationOnB,
                                                Vec3& locationOnGround, Vec3& velocityInGround) const
    {
        const Vec3& p_GB   = getBodyOriginLocation(s);
        const Vec3  p_BS_G = expressVectorInGroundFrame(s,locationOnB); // 15 flops
        locationOnGround = p_GB + p_BS_G;                               //  3 flops

        const Vec3& w_GB = getBodyAngularVelocity(s);
        const Vec3& v_GB = getBodyOriginVelocity(s);
        velocityInGround = v_GB + w_GB % p_BS_G;                        // 12 flops
    }


    /// It is cheaper to calculate a station's ground location, velocity, and acceleration together
    /// than to do them separately. Here we can return them all in 54 flops, vs. 93 to
    /// do them in three calls. This operator is available at Acceleration stage.
    void findStationLocationVelocityAndAccelerationInGround
       (const State& s, const Vec3& locationOnB,
        Vec3& locationOnGround, Vec3& velocityInGround, Vec3& accelerationInGround) const
    {
        const Rotation&  R_GB = getBodyRotation(s);
        const Vec3&      p_GB = getBodyOriginLocation(s);

        const Vec3 r = R_GB*locationOnB; // re-express station vector p_BS in G (15 flops)
        locationOnGround  = p_GB + r;    // 3 flops

        const Vec3& w = getBodyAngularVelocity(s);      // in G
        const Vec3& v = getBodyOriginVelocity(s);       // in G
        const Vec3& b = getBodyAngularAcceleration(s);  // in G
        const Vec3& a = getBodyOriginAcceleration(s);   // in G

        const Vec3 wXr = w % r; // "whipping" velocity w X r due to angular velocity (9 flops)
        velocityInGround     = v + wXr;                 // v + w X r (3 flops)
        accelerationInGround = a + b % r + w % wXr;     // 24 flops
    }

    /// Return the Cartesian (ground) location of this body B's mass center. Cost is 18 flops.
    /// This operator is available at Position stage.
    Vec3 findMassCenterLocationInGround(const State& s) const {
        return findStationLocationInGround(s,getBodyMassCenterStation(s));
    }

    /// Return the point of another body A that is currently coincident in space with the
    /// mass center CB of this body B. Cost is 36 flops. This operator is available at
    /// Position stage.
    Vec3 findMassCenterLocationInAnotherBody(const State& s, const MobilizedBody& toBodyA) const {
        return findStationLocationInAnotherBody(s,getBodyMassCenterStation(s),toBodyA);
    }

    /// Return the station (point) S of this body B that is coincident with the given Ground location.
    /// That is we return locationOnB = X_BG * locationInG, which means the result is measured
    /// from the body origin OB and expressed in the body frame. In more precise notation,
    /// we're calculating p_BS = X_BG * p_GS. Cost is 18 flops. This operator
    /// is available at Position stage.
    Vec3 findStationAtGroundPoint(const State& s, const Vec3& locationInG) const {
        return ~getBodyTransform(s) * locationInG;
    }

    /// Return the station (point) on this body B that is coincident with the given station
    /// on another body A. That is we return stationOnB = X_BA * stationOnA, which means
    /// the result is measured from the body origin OB and expressed in the body frame.
    /// Cost is 36 flops. This operator is available at Position stage.
    /// @see findStationLocationInAnotherBody()
    Vec3 findStationAtAnotherBodyStation(const State& s, const MobilizedBody& fromBodyA, 
                                const Vec3& stationOnA) const {
        return fromBodyA.findStationLocationInAnotherBody(s, stationOnA, *this);
    }

    /// Return the station S of this body that is currently coincident in space with the
    /// origin OA of another body A. Cost is 18 flops. This operator is available at
    /// Position stage.
    Vec3 findStationAtAnotherBodyOrigin(const State& s, const MobilizedBody& fromBodyA) const {
        return findStationAtGroundPoint(s,fromBodyA.getBodyOriginLocation(s));
    }

    /// Return the station S of this body that is currently coincident in space with the
    /// mass center CA of another body B. Cost is 36 flops. This operator is available at
    /// Position stage.
    Vec3 findStationAtAnotherBodyMassCenter(const State& s, const MobilizedBody& fromBodyA) const {
        return fromBodyA.findStationLocationInAnotherBody(s,getBodyMassCenterStation(s),*this);
    }

    /// Re-express a vector expressed in this body B's frame into the same vector in G, by applying
    /// only a rotation. That is, we return vectorInG = R_GB * vectorInB. Cost is 15 flops. 
    /// This operator is available at Position stage.
    Vec3 expressVectorInGroundFrame(const State& s, const Vec3& vectorInB) const {
        return getBodyRotation(s)*vectorInB;
    }

    /// Re-express a vector expressed in Ground into the same vector expressed in this body B, by
    /// applying only rotation. That is, we return vectorInB = R_BG * vectorInG. Cost is 15 flops. 
    /// This operator is available at Position stage.
    Vec3 expressGroundVectorInBodyFrame(const State& s, const Vec3& vectorInG) const {
        return ~getBodyRotation(s)*vectorInG;
    }

    /// Re-express a vector expressed in this body B into the same vector expressed in body A, by
    /// applying only a rotation. That is, we return vectorInA = R_AB * vectorInB. Cost is 30 flops.
    /// This operator is available at Position stage.
    /// Note: if you know one of the bodies is Ground, call one of the specialized methods
    /// above to save 15 flops.
    Vec3 expressVectorInAnotherBodyFrame(const State& s, const Vec3& vectorInB,
                                         const MobilizedBody& inBodyA) const
    {
        return inBodyA.expressGroundVectorInBodyFrame(s, expressVectorInGroundFrame(s,vectorInB));
    }

    /// Re-express this body B's mass properties in Ground by applying only a rotation, not a shift
    /// of reference point. The mass properties remain measured in the body B frame, taken about the body
    /// B origin OB, but are reexpressed in Ground.
    MassProperties expressMassPropertiesInGroundFrame(const State& s) {
            const MassProperties& M_OB_B = getBodyMassProperties(s);
            const Rotation&       R_GB   = getBodyRotation(s);
            return M_OB_B.reexpress(~R_GB);
    }

    /// Re-express this body B's mass properties in another body A's frame by applying only a rotation, not a shift
    /// of reference point. The mass properties remain measured in the body B frame, taken about the body
    /// B origin OB, but are reexpressed in A.
    MassProperties expressMassPropertiesInAnotherBodyFrame(const State& s, const MobilizedBody& inBodyA) {
            const MassProperties& M_OB_B = getBodyMassProperties(s);
            const Rotation        R_AB   = findBodyRotationInAnotherBody(s,inBodyA);
            return M_OB_B.reexpress(~R_AB);
    }

    // End of Basic Operators.
    //@}

    /// @name High-Level Operators
    /// High level operators combine State Access and Basic Operators with run-time tests 
    /// to calculate more complex MoblizedBody-specific quantities, with more complicated
    /// implementations that can exploit special cases at run time.


    //@{

    /**
     * Return the mass properties of body B, measured from and about
     * the B frame origin, but expressed in Ground and then returned
     * as a Spatial Inertia Matrix. The mass properties are arranged
     * in the SpatialMat like this:
     * <pre>
     *       M=[      I_OB      crossMat(m*CB) ]
     *         [ ~crossMat(m*CB)   diag(m)     ]
     * </pre>
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
    Inertia calcBodyInertiaAboutAnotherBodyStation(const State& s,
                                                   const MobilizedBody& inBodyA, 
                                                   const Vec3&          aboutLocationOnBodyA) const
    {
        // get B's mass props MB, measured about OB, exp. in B
        const MassProperties& MB_OB_B = getBodyMassProperties(s);

        // Calculate the vector from the body B origin (current "about" point) to the new "about" point PA,
        // expressed in B.
        const Vec3 r_OB_PA = findStationAtAnotherBodyStation(s, inBodyA, aboutLocationOnBodyA);

        // Now shift the "about" point for body B's inertia IB to PA, but still expressed in B.
        const Inertia IB_PA_B = MB_OB_B.calcShiftedInertia(r_OB_PA);
        
        // Finally reexpress the inertia in the A frame.
        const Rotation R_BA    = inBodyA.findBodyRotationInAnotherBody(s, *this);
        const Inertia  IB_PA_A = IB_PA_B.reexpress(R_BA);
        return IB_PA_A;
    }


    /// Calculate body B's momentum (angular, linear) measured and expressed in ground, but taken about
    /// the body origin OB.
    SpatialVec calcBodyMomentumAboutBodyOriginInGround(const State& s) {
        const MassProperties M_OB_G = expressMassPropertiesInGroundFrame(s);
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
        Vec3          v_G_CB = findStationVelocityInGround(s, M_OB_B.getMassCenter());

        return SpatialVec( I_CB_G*w_GB, mb*v_G_CB );
    }

    /// Calculate the distance from a station PB on body B to a station PA on body A.
    /// We are given the location vectors (stations) p_OB_PB and p_OA_PA, expressed in
    /// their respective frames. We return |p_PB_PA|.
    Real calcStationToStationDistance(const State& s,
                                      const Vec3&          locationOnBodyB,
                                      const MobilizedBody& bodyA,
                                      const Vec3&          locationOnBodyA) const
    {
        if (isSameMobilizedBody(bodyA))
            return (locationOnBodyA-locationOnBodyB).norm();

        const Vec3 r_OG_PB = this->findStationLocationInGround(s,locationOnBodyB);
        const Vec3 r_OG_PA = bodyA.findStationLocationInGround(s,locationOnBodyA);
        return (r_OG_PA - r_OG_PB).norm();
    }

    /// Calculate the time rate of change of distance from a fixed point PB on body B to a fixed point
    /// PA on body A. We are given the location vectors r_OB_PB and r_OA_PA, expressed in their
    /// respective frames. We return d/dt |r_OB_OA|, under the assumption that the time derivatives
    /// of the two given vectors in their own frames is zero.
    Real calcStationToStationDistanceTimeDerivative(const State& s,
                                                    const Vec3&          locationOnBodyB,
                                                    const MobilizedBody& bodyA,
                                                    const Vec3&          locationOnBodyA) const
    {
        if (isSameMobilizedBody(bodyA))
            return 0;

        Vec3 rB, rA, vB, vA;
        this->findStationLocationAndVelocityInGround(s,locationOnBodyB,rB,vB);
        bodyA.findStationLocationAndVelocityInGround(s,locationOnBodyA,rA,vA);
        const Vec3 r = rA-rB, v = vA-vB;
        const Real d = r.norm();

        // When the points are coincident, the rate of change of distance is just their relative speed.
        // Otherwise, it is the speed along the direction of separation. 
        if (d==0) return v.norm();
        else return dot(v, r/d);
    }


    /// Calculate the second time derivative of distance from a fixed point PB on body B to a fixed point
    /// PA on body A. We are given the location vectors (stations) r_OB_PB and r_OA_PA, expressed in their
    /// respective frames. We return d^2/dt^2 |p_PB_PA|, under the assumption that the time derivatives
    /// of the two given vectors in their own frames is zero.
    Real calcStationToStationDistance2ndTimeDerivative(const State& s,
                                                       const Vec3&          locationOnBodyB,
                                                       const MobilizedBody& bodyA,
                                                       const Vec3&          locationOnBodyA) const
    {
        if (isSameMobilizedBody(bodyA))
            return 0;

        Vec3 rB, rA, vB, vA, aB, aA;
        this->findStationLocationVelocityAndAccelerationInGround(s,locationOnBodyB,rB,vB,aB);
        bodyA.findStationLocationVelocityAndAccelerationInGround(s,locationOnBodyA,rA,vA,aA);

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
            else return dot(a, v/s);
        }

        // Points are separated.
        const Vec3 u = r/d;             // u is the separation direction (a unit vector from B to A) 
        const Vec3 vp = v - dot(v,u)*u; // velocity perpendicular to separation direction
        return dot(a,u) + dot(vp,v)/d;
    }


    /// TODO: not implemented yet -- any volunteers? Return the velocity of a point P moving on
    /// body B, in body A's frame, expressed in body A.
    Vec3 calcBodyMovingPointVelocityInBody(const State& s,
                                           const Vec3& locationOnBodyB, 
                                           const Vec3& velocityOnBodyB,
                                           const MobilizedBody& inBodyA) const
    {
        SimTK_ASSERT_ALWAYS(!"unimplemented method", 
            "MobilizedBody::calcBodyMovingPointVelocityInBody() is not yet implemented -- any volunteers?");
        return Vec3::getNaN();
    }


    /// TODO: not implemented yet -- any volunteers? Return the velocity of a point
    /// P moving (and possibly accelerating) on body B, in body A's frame, expressed in body A.
    Vec3 calcBodyMovingPointAccelerationInBody(const State& s, 
                                               const Vec3&          locationOnBodyB, 
                                               const Vec3&          velocityOnBodyB, 
                                               const Vec3&          accelerationOnBodyB,
                                               const MobilizedBody& inBodyA) const
    {
        SimTK_ASSERT_ALWAYS(!"unimplemented method", 
            "MobilizedBody::calcBodyMovingPointAccelerationInBody() is not yet implemented -- any volunteers?");
        return Vec3::getNaN();
    }

    /// TODO: not implemented yet -- any volunteers? Calculate the time rate of change of distance from a moving
    /// point PB on body B to a moving point
    /// PA on body A. We are given the location vectors r_OB_PB and r_OA_PA, and the velocities of
    /// PB in B and PA in A, all expressed in their respective frames. We return d/dt |r_OB_OA|,
    /// taking into account the (given) time derivatives of the locations in their local frames, as well
    /// as the relative velocities of the bodies.
    Real calcMovingPointToPointDistanceTimeDerivative(const State& s,
                                                      const Vec3&          locationOnBodyB,
                                                      const Vec3&          velocityOnBodyB,
                                                      const MobilizedBody& bodyA,
                                                      const Vec3&          locationOnBodyA,
                                                      const Vec3&          velocityOnBodyA) const
    {
        SimTK_ASSERT_ALWAYS(!"unimplemented method", 
            "MobilizedBody::calcMovingPointToPointDistanceTimeDerivative() is not yet implemented -- any volunteers?");
        return NaN;
    }

    /// TODO: not implemented yet -- any volunteers? Calculate the second time derivative of distance
    /// from a moving point PB on body B to a moving point
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
                                                         const Vec3&          accelerationOnBodyA) const
    {
        SimTK_ASSERT_ALWAYS(!"unimplemented method", 
            "MobilizedBody::calcMovingPointToPointDistance2ndTimeDerivative() is not yet implemented -- any volunteers?");
        return NaN;
    }


    // End of High Level Operators.
    //@}


            //////////////////////////
            // CONSTRUCTION METHODS //
            //////////////////////////

    /// The default constructor initializes the base class so that it contains a null
    /// implementation. This should be called only from concrete MobilizedBody 
    /// constructors.
    MobilizedBody();

    /// Internal use only
    explicit MobilizedBody(MobilizedBodyImpl* r);

    /// @name Construction and Misc Methods
    /// These methods are the base class services which are used while building a concrete
    /// MobilizedBody, or to query a MobilizedBody to find out how it was built. These are
    /// unlikely to be used by end users of MobilizedBodies.
    //@{

    /// Add decorative geometry specified relative to the new (outboard) body's reference
    /// frame B, or to the outboard mobilizer frame M attached to body B, or
    /// to the inboard mobilizer frame F attached to the parent body P. Note that
    /// the body itself may already have had some decorative geometry on it when
    /// it was first put into this MobilizedBody; in that case this just adds more.
    MobilizedBody& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)updBody().addDecoration(X_BD,g);
        return *this;
    }

    /// Add decorative geometry specified relative to the outboard mobilizer frame M
    /// attached to body B. If body B already has decorative geometry on it,
    /// this just adds some more.
    MobilizedBody& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry&);

    /// Add decorative geometry specified relative to the inboard mobilizer frame F
    /// attached to the parent body P. If body P already has decorative geometry on it,
    /// this just adds some more.
    MobilizedBody& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry&);

    /// Return a reference to the Body contained within this MobilizedBody.
    const Body& getBody() const;
    /// Return a writable reference to the Body contained within this MobilizedBody.
    /// Calling this method invalidates the MobilizedBody's topology, so the containing
    /// matter subsystem's realizeTopology() method must be called again.
    Body& updBody();

    /// Replace the Body contained within this MobilizedBody with a new one. 
    /// Calling this method invalidates the MobilizedBody's topology, so the containing
    /// matter subsystem's realizeTopology() method must be called again. A reference
    /// to this MobilizedBody is returned so that this can be chained like an assignment
    /// operator.
    MobilizedBody& setBody(const Body&);

    /// If the contained Body can have its mass properties set to the supplied value \p m
    /// its mass properties are changed, otherwise the method fails. Calling this method
    /// invalidates the MobilizedBody's topology, so the containing
    /// matter subsystem's realizeTopology() method must be called again. A reference
    /// to this MobilizedBody is returned so that this can be chained like an assignment
    /// operator.
    MobilizedBody& setDefaultMassProperties(const MassProperties& m) {
        updBody().setDefaultRigidBodyMassProperties(m); // might not be allowed
        return *this;
    }

    /// Return the mass properties of the Body stored within this MobilizedBody.
    const MassProperties& getDefaultMassProperties() const {
        return getBody().getDefaultRigidBodyMassProperties(); // every body type can do this
    }

    /// Provide a unique Motion object for this MobilizedBody. The MobilizedBody takes
    /// over ownership of the Motion object and is responsible for cleaning up 
    /// its heap space when the time comes. This is a Topology-changing operation and
    /// consequently requires write access to the MobilizedBody which will propagate
    /// to invalidate the containing Subsystem and System's topology. There can only
    /// be one Motion object per mobilizer; this method will throw an exception if
    /// there is already one here.
    void adoptMotion(Motion& ownerHandle);

    /// If there is a Motion object associated with this MobilizedBody it is removed;
    /// otherwise, nothing happens. If a Motion is deleted, the containing System's
    /// topology is invalidated.
    void clearMotion();

    /// Check whether this MobilizedBody has an associated Motion object. This does
    /// not tell you whether the Motion object is currently enabled or in use; just
    /// whether it is available.
    bool hasMotion() const;

    /// If there is a Motion object assocated with this MobilizedBody, this returns
    /// a const reference to it. Otherwise it will throw an exception. You can check first
    /// using hasMotion(). Note that there is no provision to obtain a writable
    /// reference to the contained Motion object; if you want to change it clear the
    /// existing object instead and replace it with a new one.
    /// @see hasMotion()
    const Motion& getMotion() const;

    /// Change this mobilizer's frame F on the parent body P. Calling this method
    /// invalidates the MobilizedBody's topology, so the containing
    /// matter subsystem's realizeTopology() method must be called again. A reference
    /// to this MobilizedBody is returned so that this can be chained like an assignment
    /// operator.
    MobilizedBody& setDefaultInboardFrame (const Transform& X_PF);
    /// Change this mobilizer's frame M fixed on this (the outboard) body B. Calling this method
    /// invalidates the MobilizedBody's topology, so the containing
    /// matter subsystem's realizeTopology() method must be called again. A reference
    /// to this MobilizedBody is returned so that this can be chained like an assignment
    /// operator.
    MobilizedBody& setDefaultOutboardFrame(const Transform& X_BM);

    /// Return a reference to this mobilizer's default for the frame F fixed on
    /// the parent body P, as the fixed Transform
    /// from P's body frame to the frame F fixed to P. This default Transform is stored
    /// with the MobilizedBody object, not the State.
    const Transform& getDefaultInboardFrame()  const; // X_PF
    /// Return a reference to this MobilizedBody's default for mobilizer frame M, as the fixed Transform
    /// from this body B's frame to the frame M fixed on B. This default Transform is stored
    /// with the MobilizedBody object, not the State.
    const Transform& getDefaultOutboardFrame() const; // X_BM

    /// This is an implicit conversion from MobilizedBody to MobilizedBodyIndex 
    /// when needed. This will fail unless this MobilizedBody is owned by some 
    /// SimbodyMatterSubsystem. We guarantee that the MobilizedBodyIndex of a
    /// mobilized body is numerically larger than the MobilizedBodyIndex of its 
    /// parent.
    operator MobilizedBodyIndex() const {return getMobilizedBodyIndex();}

    /// Return the MobilizedBodyIndex of this MobilizedBody within the owning 
    /// SimbodyMatterSubsystem. This will fail unless this MobilizedBody is 
    /// owned by some SimbodyMatterSubsystem. We guarantee that the 
    /// MobilizedBodyIndex of a mobilized body is numerically larger than the 
    /// MobilizedBodyIndex of its parent.
    MobilizedBodyIndex     getMobilizedBodyIndex()  const;

    /// Return a reference to the MobilizedBody serving as the parent body of 
    /// the current MobilizedBody. This call will fail if the current 
    /// MobilizedBody is Ground, since Ground has no parent.
    const MobilizedBody&   getParentMobilizedBody() const;

    /// Return a reference to this MobilizedBody's oldest ancestor other than 
    /// Ground, or return Ground if this MobilizedBody is Ground. That is, we 
    /// return the "base" MobilizedBody for this MobilizedBody, meaning the one 
    /// which connects this branch of the multibody tree directly to Ground.
    const MobilizedBody&   getBaseMobilizedBody()   const;

    /// Obtain a reference to the SimbodyMatterSubsystem which contains this 
    /// MobilizedBody. This will fail unless this MobilizedBody is owned by 
    /// some SimbodyMatterSubsystem.
    const SimbodyMatterSubsystem& getMatterSubsystem() const;
    /// Obtain a writable reference to the SimbodyMatterSubsystem which 
    /// contains this MobilizedBody. This will fail unless this MobilizedBody 
    /// is owned by some SimbodyMatterSubsystem.
    SimbodyMatterSubsystem&       updMatterSubsystem();

    /// Determine whether the current MobilizedBody object is owned by a matter 
    /// subsystem.
    bool isInSubsystem() const;

    /// Determine whether a given MobilizedBody \p mBody is in the same matter 
    /// subsystem as the current body. If the bodies are not in a subsystem, 
    /// this routine will return \c false.
    bool isInSameSubsystem(const MobilizedBody&) const;

    /// Determine whether a given MobilizedBody \p mBody is the same  
    /// MobilizedBody as this one. For this to be true the handles must not be 
    /// empty, and the implementation objects must be <em>the same object</em> 
    /// not separate objects with identical contents.
    bool isSameMobilizedBody(const MobilizedBody& mBody) const;

    /// Determine whether this body is Ground, meaning that it is actually 
    /// body 0 of some matter subsytem, not just that its body type is Ground.
    bool isGround() const;

    /// Return this body's level in the tree of bodies, starting with ground 
    /// at 0, bodies directly connected to ground at 1, bodies directly 
    /// connected to those at 2, etc. This is callable after realizeTopology(). 
    /// This is the graph distance of the body from Ground.
    int getLevelInMultibodyTree() const;
    
    /// Create a new MobilizedBody which is identical to this one, except that 
    /// it has a different parent (and consequently might belong to a different 
    /// MultibodySystem).
    MobilizedBody& cloneForNewParent(MobilizedBody& parent) const;


        // Utility operators //

    /// This utility selects one of the q's (generalized coordinates) associated 
    /// with this mobilizer from a supplied "q-like" Vector, meaning a Vector 
    /// which is the same length as the Vector of q's for the containing matter 
    /// subsystem.
    Real  getOneFromQPartition(const State&, int which, const Vector& qlike) const;

    /// This utility returns a writable reference to one of the q's (generalized 
    /// coordinates) associated with this mobilizer from a supplied "q-like" 
    /// Vector, meaning a Vector which is the same length as the Vector of q's 
    /// for the containing matter subsystem.
    Real& updOneFromQPartition(const State&, int which, Vector& qlike) const;

    /// This utility selects one of the u's (generalized speeds) associated with 
    /// this mobilizer from a supplied "u-like" Vector, meaning a Vector which is 
    /// the same length as the Vector of u's for the containing matter subsystem.
    Real  getOneFromUPartition(const State&, int which, const Vector& ulike) const;

    /// This utility returns a writable reference to one of the u's (generalized 
    /// speeds) associated with this mobilizer from a supplied "u-like" Vector, 
    /// meaning a Vector which is the same length as the Vector of u's for the 
    /// containing matter subsystem.
    Real& updOneFromUPartition(const State&, int which, Vector& ulike) const;

    /// This utility adds in the supplied generalized force \p force (a scalar) to the
    /// appropriate slot of the supplied \p mobilityForces Vector, which is a "u-like"
    /// Vector. Note that we are <em>adding</em> this not <em>setting</em> it so it
    /// important that \p mobilityForces be initialized to zero before making a set
    /// of calls to applyOneMobilityForce().
    void applyOneMobilityForce(const State& s, int which, Real f, 
                               Vector& mobilityForces) const
    {
        updOneFromUPartition(s,which,mobilityForces) += f;
    }

    /// This utility adds in the supplied spatial force \p spatialForceInG (consisting
    /// of a torque vector, and a force vector to be applied at the current body's
    /// origin) to the appropriate slot of the supplied \p bodyForcesInG Vector. 
    /// Note that we are <em>adding</em> this not <em>setting</em> it so it
    /// important that \p mobilityForces be initialized to zero before making a set
    /// of calls to applyBodyForce().
    void applyBodyForce(const State& s, const SpatialVec& spatialForceInG, 
                        Vector_<SpatialVec>& bodyForcesInG) const;

    /// This utility adds in the supplied pure torque \p torqueInG to the appropriate
    /// slot of the supplied \p bodyForcesInG Vector. Note that we are <em>adding</em>
    /// this not <em>setting</em> it so it
    /// important that \p bodyForcesInG be initialized to zero before making a set
    /// of calls to applyBodyTorque().
    void applyBodyTorque(const State& s, const Vec3& torqueInG, 
                         Vector_<SpatialVec>& bodyForcesInG) const;

    /// This utility adds in the supplied force \p forceInG applied at a point \p pointInB
    /// to the appropriate slot of the supplied \p bodyForcesInG Vector. Notes: 
    ///    - we are <em>adding</em> this not <em>setting</em> it so it
    ///      important that \p bodyForcesInG be initialized to zero before making a set
    ///      of calls to applyForceToBodyPoint().
    ///    - \p pointInB represents a fixed station of B and is provided by giving the
    ///      vector from body B's origin to the point, expressed in the B frame, while
    ///      the applied force (and resulting body forces and torques) are expressed
    ///      in the ground frame.
    void applyForceToBodyPoint(const State& s, const Vec3& pointInB, const Vec3& forceInG,
                               Vector_<SpatialVec>& bodyForcesInG) const;

    // End of Construction and Misc Methods.
    //@}

        /////////////////////////////////////
        // BUILT IN MOBILIZER DECLARATIONS //
        /////////////////////////////////////

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


    class Pin;              typedef Pin             Torsion;
    class Universal;
    class Cylinder;

    class Weld;
    class Slider;           typedef Slider          Prismatic;
    class Translation2D;    typedef Translation2D   Cartesian2D, CartesianCoords2D;
    class Translation;      typedef Translation     Cartesian, CartesianCoords;
    class BendStretch;      typedef BendStretch     PolarCoords;
    class TorsionStretch;   typedef TorsionStretch  ConicalCoords2D;
    class SphericalCoords;
    class CylindricalCoords;
    class LineOrientation;


    class Planar;
    class Gimbal;
    class Ball;             typedef Ball            Orientation, Spherical;
    class Free;
    class FreeLine;
    class Screw;
    class Ellipsoid;
    class Custom;
    class Ground;
    class FunctionBased;
    
    class PinImpl;
    class SliderImpl;
    class UniversalImpl;
    class CylinderImpl;
    class BendStretchImpl;
    class TorsionStretchImpl;
    class PlanarImpl;
    class GimbalImpl;
    class BallImpl;
    class TranslationImpl;
    class SphericalCoordsImpl;
    class FreeImpl;
    class LineOrientationImpl;
    class FreeLineImpl;
    class WeldImpl;
    class ScrewImpl;
    class EllipsoidImpl;
    class CustomImpl;
    class GroundImpl;
    class FunctionBasedImpl;

};

/// One mobility -- rotation about the common z axis of the inboard
/// and outboard mobilizer frames.
/// Synonym: Torsion
class SimTK_SIMBODY_EXPORT MobilizedBody::Pin : public MobilizedBody {
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
    explicit Pin(Direction=Forward);
    Pin(MobilizedBody& parent, const Body&, Direction=Forward);
    Pin(MobilizedBody& parent, const Transform& inbFrame,
        const Body&,           const Transform& outbFrame,
        Direction=Forward);

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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Pin, PinImpl, MobilizedBody);
};

/// One mobility -- translation along the common x axis of the
/// inboard and outboard mobilizer frames.
/// Synonym: Prismatic
class SimTK_SIMBODY_EXPORT MobilizedBody::Slider : public MobilizedBody {
public:
        // SPECIALIZED INTERFACE FOR SLIDER MOBILIZER

    // "Length" is just a nicer name for a sliding joint's lone generalized coordinate q.
    Slider& setDefaultLength(Real length) {return setDefaultQ(length);}
    Real getDefaultLength() const         {return getDefaultQ();}

        // Friendly, mobilizer-specific access to generalized coordinates and speeds.

    void setLength(State& s, Real length) {setQ(s, length);}
    Real getLength(const State& s) const {return getQ(s);}

    void setRate(State& s, Real rateInLengthPerTime) {setU(s, rateInLengthPerTime);}
    Real getRate(const State& s) const {return getU(s);}

    // Mobility forces are "u-like", that is, one per dof.
    Real getAppliedForce(const State& s, const Vector& mobilityForces) const {
        return getMyPartU(s,mobilityForces);
    }
    void applyForce(const State& s, Real force, Vector& mobilityForces) const {
        updMyPartU(s,mobilityForces) += force;
    }

        // STANDARDIZED MOBILIZED BODY INTERFACE

        // required constructors
    explicit Slider(Direction=Forward);
    Slider(MobilizedBody& parent, const Body&, Direction=Forward);
    Slider(MobilizedBody& parent, const Transform& inbFrame,
           const Body&,           const Transform& outbFrame, Direction=Forward);

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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Slider, SliderImpl, MobilizedBody);
};

/// One mobility -- coordinated rotation and translation along the
/// common z axis of the inboard and outboard mobilizer frames. A
/// "pitch" is specified relating the two. The generalized coordinate
/// q is the rotation angle in radians, the translation is always
/// pitch*q.
class SimTK_SIMBODY_EXPORT MobilizedBody::Screw : public MobilizedBody {
public:
    explicit Screw(Real pitch, Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Screw(MobilizedBody& parent, const Body&, Real pitch, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Screw(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame,
         Real pitch, Direction=Forward);

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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Screw, ScrewImpl, MobilizedBody);
};

/// Two mobilities -- rotation about the x axis, followed by a rotation
/// about the new y axis. This mobilizer is badly behaved when the
/// second rotation is near 90 degrees.
class SimTK_SIMBODY_EXPORT MobilizedBody::Universal : public MobilizedBody {
public:
    explicit Universal(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Universal(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Universal(MobilizedBody& parent, const Transform& inbFrame,
              const Body&,           const Transform& outbFrame, Direction=Forward);

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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Universal, UniversalImpl, MobilizedBody);
};

/// Two mobilities -- rotation and translation along the common z axis
/// of the inboard and outboard mobilizer frames.
class SimTK_SIMBODY_EXPORT MobilizedBody::Cylinder : public MobilizedBody {
public:
    explicit Cylinder(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Cylinder(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Cylinder(MobilizedBody& parent, const Transform& inbFrame,
             const Body&,           const Transform& outbFrame, Direction=Forward);

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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Cylinder, CylinderImpl, MobilizedBody);
};

/// Two mobilities: The z axis of the parent's F frame is 
/// used for rotation (and that is always aligned with the M frame z axis).
/// The x axis of the *M* (outboard) frame is then used for translation;
/// that is, first we rotate around z, which moves M's x with respect to F's x. Then
/// we slide along the rotated x axis. The two generalized coordinates are the
/// rotation and the translation, in that order.
/// This can also be viewed a a 2D polar coordinate mobilizer since the coordinates
/// are (theta, r) about perpendicular axes.
class SimTK_SIMBODY_EXPORT MobilizedBody::BendStretch : public MobilizedBody {
public:
    explicit BendStretch(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    BendStretch(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    BendStretch(MobilizedBody& parent, const Transform& inbFrame,
                const Body&,           const Transform& outbFrame, Direction=Forward);

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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(BendStretch, BendStretchImpl, MobilizedBody);
};

/// Three mobilities -- z rotation and x,y translation. The generalized
/// coordinates are rotation about the shared z axis of the F and M
/// frame, translation along the F frame's x axis, and translation along
/// its y axis, in that order.
class SimTK_SIMBODY_EXPORT MobilizedBody::Planar : public MobilizedBody {
public:
    explicit Planar(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Planar(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Planar(MobilizedBody& parent, const Transform& inbFrame,
           const Body&,           const Transform& outbFrame,
           Direction=Forward);

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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Planar, PlanarImpl, MobilizedBody);
};

/**
 * Three mobilities -- body fixed 3-2 (z-y) rotation followed by translation
 * along body z or body x. Interpreted as spherical coordinates the first rotation
 * is the azimuth angle, the second is the zenith, and the translation is
 * the radius. We permit a simple mapping from generalized coordinates to
 * (azimuth, zenith, radius):
 * <pre>
 *     azimuth = s0*q0 + az0   (about Fz==Mz)
 *     zenith  = s1*q1 + ze0   (about My)
 *     radius  = s2*q2         (along Mz or Mx; Mz is default)
 * </pre>
 * where s0,s1,s2 are signs (1 or -1) and az0 and ze0 are offset angles. The
 * F and M frames are coincident when azimuth==zenith==radius==0. But note
 * that with non-zero offsets the F and M frames will not be aligned in
 * the reference configuration where q0==q1==q2==0. The F and M origins
 * will always be coincident when q2==0, however.
 *
 * This mobilizer can be used to give unrestricted 3-d motion to inertialess 
 * particles (as with a Cartesian mobilizer but parameterized torsion,bend,stretch
 * instead of x,y,z) but in this case you must watch for two possible 
 * singularities: (1) radius==0, and (2) zenith==n*Pi (or equivalently 
 * q1==n*Pi-s1*ze0). If your operating range steers clear of those singularities, 
 * you're fine.
 */
class SimTK_SIMBODY_EXPORT MobilizedBody::SphericalCoords : public MobilizedBody {
public:
    explicit SphericalCoords(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    SphericalCoords(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames. This gives you a pure spherical
    /// coordinate system in which q0=azimuth about Fz(==Mz), q1=zenith about My, 
    /// and q2=radius along Mz.
    SphericalCoords(MobilizedBody& parent, const Transform& inbFrame,
                    const Body&,           const Transform& outbFrame,
                    Direction=Forward);

    /// Use this constructor to specify the general case described above.
    SphericalCoords(MobilizedBody& parent,      const Transform& inbFrame,
                    const Body&,                const Transform& outbFrame,
                    Real azimuthOffset,         bool azimuthNegated,
                    Real zenithOffset,          bool zenithNegated,
                    CoordinateAxis radialAxis,  bool radialNegated,
                    Direction=Forward);

    SphericalCoords& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    SphericalCoords& addOutboardDecoration(const Transform& X_MD,  const DecorativeGeometry& g) {
        (void)MobilizedBody::addOutboardDecoration(X_MD,g); return *this;
    }
    SphericalCoords& addInboardDecoration (const Transform& X_FD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addInboardDecoration(X_FD,g); return *this;
    }

    SphericalCoords& setDefaultInboardFrame(const Transform& X_PF) {
        (void)MobilizedBody::setDefaultInboardFrame(X_PF); return *this;
    }

    SphericalCoords& setDefaultOutboardFrame(const Transform& X_BM) {
        (void)MobilizedBody::setDefaultOutboardFrame(X_BM); return *this;
    }

    // Friendly, mobilizer-specific access to coordinates and speeds.
    SphericalCoords& setDefaultAngles(const Vec2& a) {
        Vec3 q = getDefaultQ(); q.updSubVec<2>(0) = a; setDefaultQ(q);
        return *this;
    }
    SphericalCoords& setDefaultRadius(Real r) {
        Vec3 q = getDefaultQ(); q[2] = r; setDefaultQ(q);
        return *this;
    }
    SphericalCoords& setRadialAxis(CoordinateAxis);
    SphericalCoords& setNegateAzimuth(bool);
    SphericalCoords& setNegateZenith(bool);
    SphericalCoords& setNegateRadial(bool);

    const Vec2&    getDefaultAngles()      const {return getDefaultQ().getSubVec<2>(0);}
    Real           getDefaultTranslation() const {return getDefaultQ()[2];}
    
    CoordinateAxis getRadialAxis()    const;
    bool           isAzimuthNegated() const;
    bool           isZenithNegated()  const;
    bool           isRadialNegated()  const;

    void setAngles(State& s, const Vec2& a) {setOneQ(s,0,a[0]); setOneQ(s,1,a[1]);}
    void setRadius(State& s, Real        r) {setOneQ(s,2,r);}

    const Vec2& getAngles(const State& s) const {return getQ(s).getSubVec<2>(0);}
    Real        getRadius(const State& s) const {return getQ(s)[2];}

    // Generic default state Topology methods.
    const Vec3& getDefaultQ() const;
    SphericalCoords& setDefaultQ(const Vec3& q);

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

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(SphericalCoords, SphericalCoordsImpl, MobilizedBody);
};

/// Three mobilities -- unrestricted orientation modeled as a 1-2-3
/// body-fixed Euler angle sequence. This is singular when the middle
/// angle is 90 degrees.
class SimTK_SIMBODY_EXPORT MobilizedBody::Gimbal : public MobilizedBody {
public:
    explicit Gimbal(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Gimbal(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Gimbal(MobilizedBody& parent, const Transform& inbFrame,
           const Body&,           const Transform& outbFrame, Direction=Forward);

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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Gimbal, GimbalImpl, MobilizedBody);
};

/// Three mobilities -- unrestricted orientation modeled with a
/// quaternion which is never singular. A modeling option allows the
/// joint to use a 1-2-3 Euler sequence (identical to a Gimbal) 
/// instead.
class SimTK_SIMBODY_EXPORT MobilizedBody::Ball : public MobilizedBody {
public:
    explicit Ball(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Ball(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Ball(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame, Direction=Forward);

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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Ball, BallImpl, MobilizedBody);
};

/// Three mobilities -- coordinated rotation and translation along the
/// surface of an ellipsoid fixed to the parent (inboard) body.
/// The generalized coordinates are the same as for a Ball (Orientation)
/// joint, that is, a quaternion or 1-2-3 Euler sequence.
class SimTK_SIMBODY_EXPORT MobilizedBody::Ellipsoid : public MobilizedBody {
public:
    /// The ellipsoid is placed on the mobilizer's inboard frame F, with
    /// half-axis dimensions along F's x,y,z respectively.
    explicit Ellipsoid(Direction=Forward); // not very useful until radii are set, but has some defaults

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Ellipsoid(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Ellipsoid(MobilizedBody& parent, const Transform& inbFrame,
              const Body&,           const Transform& outbFrame,
              Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames, and give the radii at
    /// the same time.
    Ellipsoid(MobilizedBody& parent, const Transform& inbFrame,
              const Body&,           const Transform& outbFrame,
              const Vec3& radii,
              Direction=Forward);

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

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Ellipsoid, EllipsoidImpl, MobilizedBody);
};

/// Three translational mobilities. The generalized coordinates are
/// x,y,z translations along the parent (inboard) F frame axes.
class SimTK_SIMBODY_EXPORT MobilizedBody::Translation : public MobilizedBody {
public:
    explicit Translation(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Translation(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Translation(MobilizedBody& parent, const Transform& inbFrame,
                const Body&,           const Transform& outbFrame, Direction=Forward);

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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Translation, TranslationImpl, MobilizedBody);
};

/// Unrestricted motion for a rigid body (six mobilities). Orientation
/// is modeled the same as for the Orientation mobilizer, that is, using
/// quaternions to avoid singularities. A modeling option exists to 
/// have the joint modeled with a 1-2-3 body fixed Euler sequence like
/// a Gimbal mobilizer. Translational generalized coordinates are
/// x,y,z translations along the F (inboard) axes.
class SimTK_SIMBODY_EXPORT MobilizedBody::Free : public MobilizedBody {
public:
    explicit Free(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    Free(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    Free(MobilizedBody& parent, const Transform& inbFrame,
         const Body&,           const Transform& outbFrame,
         Direction=Forward);

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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Free, FreeImpl, MobilizedBody);
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
class SimTK_SIMBODY_EXPORT MobilizedBody::LineOrientation : public MobilizedBody {
public:
    explicit LineOrientation(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    LineOrientation(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    LineOrientation(MobilizedBody& parent, const Transform& inbFrame,
                    const Body&,           const Transform& outbFrame, Direction=Forward);

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

    // This is just a nicer name for the generalized coordinate.
    LineOrientation& setDefaultRotation(const Rotation& R_FM) {
        return setDefaultQ(R_FM.convertRotationToQuaternion());
    }
    Rotation getDefaultRotation() const {return Rotation(getDefaultQ());}

    // Generic default state Topology methods.
    const Quaternion& getDefaultQ() const;
    LineOrientation& setDefaultQ(const Quaternion& q);

    const Vec4& getQ(const State&) const;
    const Vec4& getQDot(const State&) const;
    const Vec4& getQDotDot(const State&) const;
    const Vec2& getU(const State&) const;
    const Vec2& getUDot(const State&) const;

    void setQ(State&, const Vec4&) const;
    void setU(State&, const Vec2&) const;

    const Vec4& getMyPartQ(const State&, const Vector& qlike) const;
    const Vec2& getMyPartU(const State&, const Vector& ulike) const;
   
    Vec4& updMyPartQ(const State&, Vector& qlike) const;
    Vec2& updMyPartU(const State&, Vector& ulike) const;
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(LineOrientation, LineOrientationImpl, MobilizedBody);
};

/// Five mobilities, representing unrestricted motion for a body which is
/// inertialess along its own z axis. The rotational generalized coordinates are the same
/// as for the LineOrientation mobilizer. The translational coordinates are
/// the same as in a Free mobilizer, or a Cartesian (Translation) mobilizer.
class SimTK_SIMBODY_EXPORT MobilizedBody::FreeLine : public MobilizedBody {
public:
    explicit FreeLine(Direction=Forward);

    /// By default the parent body frame and the body's own frame are
    /// used as the inboard and outboard mobilizer frames, resp.
    FreeLine(MobilizedBody& parent, const Body&, Direction=Forward);

    /// Use this constructor to specify mobilizer frames which are
    /// not coincident with the body frames.
    FreeLine(MobilizedBody& parent, const Transform& inbFrame,
             const Body&,           const Transform& outbFrame, Direction=Forward);

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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(FreeLine, FreeLineImpl, MobilizedBody);
};

/// Zero mobilities. This degenerate "mobilizer" serves only to weld together
/// the M frame of a body to the F frame on its parent.
class SimTK_SIMBODY_EXPORT MobilizedBody::Weld : public MobilizedBody {
public:
    /// Note: there is no "reverse" weld, because "reverse" refers to
    /// how the q's and u's are defined and there are none.
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
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Weld, WeldImpl, MobilizedBody);
};


/// This is a special type of "mobilized" body used as a placeholder for Ground
/// in the 0th slot for a MatterSubsystem's mobilized bodies.
/// The body type will also be Ground.
class SimTK_SIMBODY_EXPORT MobilizedBody::Ground : public MobilizedBody {
public:
    /// There is no "reverse" Ground.
    Ground();
    Ground& addBodyDecoration(const Transform& X_BD, const DecorativeGeometry& g) {
        (void)MobilizedBody::addBodyDecoration(X_BD,g); return *this;
    }
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Ground, GroundImpl, MobilizedBody);
};


/**
 * The handle class MobilizedBody::Custom (dataless) and its companion class MobilizedBody::Custom::Implementation
 * can be used together to define new MobilizedBody types with arbitrary properties. To use it, create a class
 * that extends MobilizedBody::Custom::Implementation. You can then create an instance of it and pass it to the
 * MobilizedBody::Custom constructor:
 * 
 * <pre>
 * MobilizedBody::Custom myMobilizedBody(new MyMobilizedBodyImplementation(args));
 * </pre>
 * 
 * Alternatively, you can also create a new Handle class which is a subclass of MobilizedBody::Custom
 * and which creates the Implementation itself in its constructors.
 * 
 * <pre>
 * class MyMobilizedBody : public MobilizedBody::Custom {
 * public:
 *   MyMobilizedBody(args) : MobilizedBody::Custom(new MyForceImplementation(args)) {
 *   }
 * }
 * </pre>
 * 
 * This allows an end user to simply write
 * 
 * <pre>
 * MyMobilizedBody(args);
 * </pre>
 * 
 * and not worry about implementation classes or creating objects on the heap.  If you do this, your MobilizedBody::Custom
 * subclass must not have any data members or virtual methods.  If it does, it will not work correctly.  Instead,
 * store all data in the Implementation subclass.
 */
class SimTK_SIMBODY_EXPORT MobilizedBody::Custom : public MobilizedBody {
public:
    class Implementation;
    class ImplementationImpl;

    /* Create a Custom MobilizedBody.
     * 
     * @param parent         the MobilizedBody's parent body
     * @param implementation the object which implements the custom mobilized body.  The MobilizedBody::Custom takes over
     *                       ownership of the implementation object, and deletes it when the MobilizedBody itself
     *                       is deleted.
     * @param body           describes this MobilizedBody's physical properties
     * @param direction      whether you want the coordinates defined as though parent & child were swapped
     */
    Custom(MobilizedBody& parent, Implementation* implementation, const Body& body, Direction direction=Forward);
    /* Create a Custom MobilizedBody.
     * 
     * @param parent         the MobilizedBody's parent body
     * @param implementation the object which implements the custom mobilized body.  The MobilizedBody::Custom takes over
     *                       ownership of the implementation object, and deletes it when the MobilizedBody itself
     *                       is deleted.
     * @param inbFrame       the MobilizedBody's inboard reference frame
     * @param body           describes this MobilizedBody's physical properties
     * @param outbFrame      the MobilizedBody's outboard reference frame
     * @param direction      whether you want the coordinates defined as though parent & child were swapped
     */
    Custom(MobilizedBody& parent, Implementation* implementation, 
           const Transform& inbFrame, const Body& body, const Transform& outbFrame,
           Direction direction=Forward);
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Custom, CustomImpl, MobilizedBody);
protected:
    const Implementation& getImplementation() const;
    Implementation&       updImplementation();
};

// We only want the template instantiation to occur once. This symbol is defined in the SimTK core
// compilation unit that defines the MobilizedBody class but should not be defined any other time.
// BE SURE TO DEAL WITH THIS IN MobilizedBody_Instantiation.cpp.
#ifndef SimTK_SIMBODY_DEFINING_MOBILIZED_BODY
    extern template class PIMPLHandle<MobilizedBody::Custom::Implementation, MobilizedBody::Custom::ImplementationImpl>;
#endif


class SimTK_SIMBODY_EXPORT MobilizedBody::Custom::Implementation 
  : public PIMPLHandle<Implementation,ImplementationImpl> 
{
public:
    // No default constructor because you have to supply at least the SimbodyMatterSubsystem
    // to which this MobilizedBody belongs.

    /// Destructor is virtual so derived classes get a chance to clean up if necessary.
    virtual ~Implementation();

    /// This method should produce a deep copy identical to the concrete derived Implementation
    /// object underlying this Implementation base class object.
    /// Note that the result is new heap space; the caller must be sure to take ownership
    /// of the returned pointer and call delete on it when done.
    virtual Implementation* clone() const = 0;

    /// This Implementation base class constructor sets the topological defaults for
    /// the number of mobilities (generalized speeds) u, the number of generalized
    /// coordinates q, and the number of those q's that are angles. There can be up
    /// to 3 angular coordinates (which must be measured in radians). You also can
    /// specify 4 as the number of angles, which is interpreted to mean the the
    /// mobilizer uses a quaternion to represent orientation. Because quaternions
    /// are not appropriate for some calculations, however, the user may globally
    /// disable them by calling setUseEulerAngles() on the SimbodyMatterSubsystem.
    /// Therefore, if you specify nAngles=4, the actual number of angular state variables
    /// may be either 3 (a set of Euler angles) or 4 (quaternion components), and the
    /// total number of state variables could be either nq-1 or nq. Before
    /// interpreting the state variables, you must first call getUseEulerAngles() to
    /// determine which representation is in use.
    ///
    /// In any case, if there are any angular coordinates they must be the <i>first</i>
    /// coordinates in the array of q's associated with this mobilizer. Translational
    /// or other q's will immediately follow the angular ones. This permits Simbody
    /// to handle quaternion normalization and conversion automatically, and to find angles which
    /// need to have their sines and cosines calculated.
    ///
    /// NOTE: if you don't say there are any angles, you can mange things yourself.
    /// However, there is no way to get quaternions normalized and converted if you don't tell
    /// Simbody about them.
    Implementation(SimbodyMatterSubsystem&, int nu, int nq, int nAngles=0);

    /// Return a Vector containing all the generalized coordinates q currently in use by this mobilizer.
    /// Note that if this mobilizer uses quaternions, the number of q's will depend on whether
    /// quaternions are currently enabled.  Call getUseEulerAngles() to check this.
    Vector getQ(const State& s) const;
    
    /// Return a Vector containing all the generalized speeds u currently in use by this mobilizer.
    Vector getU(const State& s) const;

    /// Return a Vector containing all the generalized coordinate derivatives qdot currently in use by this mobilizer.
    /// Note that if this mobilizer uses quaternions, the number of q's will depend on whether
    /// quaternions are currently enabled.  Call getUseEulerAngles() to check this.
    Vector getQDot(const State& s) const;

    /// Return a Vector containing all the generalized accelerations udot currently in use by this mobilizer.
    Vector getUDot(const State& s) const;
    
    /// Return a Vector containing all the generalized coordinate second derivatives qdotdot currently in use by this mobilizer.
    /// Note that if this mobilizer uses quaternions, the number of q's will depend on whether
    /// quaternions are currently enabled.  Call getUseEulerAngles() to check this.
    Vector getQDotDot(const State& s) const;

    /// Get the cross-mobilizer transform X_FM, the body's "moving" mobilizer frame M measured and expressed in
    /// the parent body's corresponding "fixed" frame F.  The state must have been realized to at least
    /// Position stage. Note: this refers to F and M <em>as defined</em>, not as they are if the 
    /// mobilizer has been reversed (that is, we're really returning X_F0M0 here).
    Transform getMobilizerTransform(const State& s) const;

    /// Get the cross-mobilizer velocity V_FM, the relative velocity of this body's "moving" mobilizer
    /// frame M in the parent body's corresponding "fixed" frame F, measured and expressed in F.
    /// Note that this isn't the usual spatial velocity since it isn't expressed in G.
    /// The state must have been realized to at least Velocity stage. Note: this refers to 
    /// F and M <em>as defined</em>, not as they are if the 
    /// mobilizer has been reversed (that is, we're really returning V_F0M0 here).
    SpatialVec getMobilizerVelocity(const State& s) const;

    /// Get whether rotations are being represented as quaternions or Euler angles.
    /// This method is only relevant if the constructor was invoked with nAngles==4.
    /// If this returns false, the first four q's should be interpreted as the
    /// components of a (possibly not normalized) quaternion.  If it returns true, the
    /// first three q's should be interpreted as Euler angles.
    ///
    /// Note that the total number of state variables is one less when using Euler
    /// angles than when using quaternions.
    bool getUseEulerAngles(const State& s) const;

    /// Call this if you want to make sure that the next realizeTopology() call does
    /// something. This is done automatically when you modify the MobilizedBody in ways
    /// understood by Simbody. But if you are just
    /// changing some of your own topology and want to make sure you get a chance to
    /// recompute something in realizeTopology(), make this call at the time of 
    /// modification.
    void invalidateTopologyCache() const;

    /// @name MobilizedBody Virtuals
    /// These must be defined for any Custom MobilizedBody.
    /// Note that the numbers nu, nq, and nAngles are passed in to these routines for
    /// redundancy -- you should make sure they have the values you are expecting!
    //@{

    /// Given values for this mobilizer's nq generalized coordinates q, compute X_FM(q), that is,
    /// the cross-mobilizer spatial Transform giving the configuration of the "moving" frame M
    /// fixed to the outboard (child) body B in the "fixed" frame F attached to the inboard (parent)
    /// body P. The state is guaranteed to have been realized to at least Instance stage.
    virtual Transform calcMobilizerTransformFromQ(const State& s, int nq, const Real* q) const = 0;


    /// Calculate V_FM(u) = H*u where H=H(q) is the joint transition matrix mapping the mobilities to
    /// the relative spatial velocity between the F frame on the parent to the M frame on the child.
    /// The state is guaranteed to have been realized to at least Position stage.
    ///
    /// IMPORTANT -- H should depend only on X_FM(q), not directly on q, since different sets
    /// of q's can generate the same Transform (e.g. quaternions and Euler angles). You can
    /// call getMobilizerTransform(s) to get the already calculated Transform.
    ///
    /// EVEN MORE IMPORTANT -- H here must be the same as the H^T used in multiplyByHTranspose(),
    /// and the HDot methods must use the time derivative of H.
    ///
    /// Note: the "H" we're using here is the transpose of what is used in Schwieter's IVM
    /// paper and in all of Abhi Jain's papers. That's because Jain used H^T as the joint 
    /// kinematics Jacobian, with H being the force transmission matrix which no 
    /// mobilizer-writing user is going to be thinking about.
    /// @see multiplyByHTranspose()
    virtual SpatialVec multiplyByHMatrix(const State& s, int nu, const Real* u) const = 0;

    /// Calculate f = ~H*F where F is a spatial force (torque+force) and f is its mapping onto
    /// the mobilities.
    ///
    /// IMPORTANT -- H should depend only on X_FM(q), not directly on q, since different sets
    /// of q's can generate the same Transform (e.g. quaternions and Euler angles). You can
    /// call getMobilizerTransform(s) to get the already calculated Transform.
    /// H here must match H and HDot in the other methods for this mobilizer.
    /// @see multiplyByHMatrix()
    virtual void multiplyByHTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const = 0;

    /// Calculate A0_FM = HDot*u where HDot=HDot(q,u) is the time derivative of H. This calculates
    /// the "bias acceleration" due to coriolis effects, such that the full cross-mobilizer
    /// acceleration is A_FM=A0_FM + H*udot.
    /// The state is guaranteed to have been realized to at least Velocity stage.
    ///
    /// IMPORTANT -- HDot should depend only on X_FM(q) and V_FM(q,u), not directly on q or u,
    /// since different choices of coordinates can generate the same X and V, but all such choices
    /// must produce the same H and HDot. You can call getMobilizerTransform(s) to get the already
    /// calculated Transform, and getMobilizerVelocity(s) to get the already calculated velocity.
    virtual SpatialVec multiplyByHDotMatrix(const State& s, int nu, const Real* u) const = 0;

    /// Calculate f = ~HDot*F where F is a spatial vector and f is its mapping onto
    /// the mobilities.
    /// The state is guaranteed to have been realized to at least Velocity stage.
    ///
    /// IMPORTANT -- HDot should depend only on X_FM(q) and V_FM(q,u), not directly on q or u,
    /// since different choices of coordinates can generate the same X and V, but all such choices
    /// must produce the same H and HDot. You can call getMobilizerTransform(s) to get the already
    /// calculated Transform, and getMobilizerVelocity(s) to get the already calculated velocity.
    virtual void multiplyByHDotTranspose(const State& s, const SpatialVec& F, int nu, Real* f) const = 0;

    /// Calculate out_q = N(q)*in_u (e.g., qdot=N*u)
    /// or out_u = ~N*in_q. Note that one of "in" and "out" is always "q-like" while
    /// the other is "u-like", but which is which changes if the matrix is transposed.
    /// Note that the transposed operation here is the same as multiplying by N on
    /// the right, with the Vectors viewed as RowVectors instead.
    /// The default implementation assumes that N is an identity matrix, and will
    /// only work if nq=nu=nIn=nOut and nAngles < 4 (i.e., no quaternions). If this
    /// is true for your mobilizer, you do not need to implement this method.
    ///
    /// The state is guaranteed to have been realized to at least Position stage.
    virtual void multiplyByN(const State& s, bool transposeMatrix, 
                             int nIn, const Real* in, int nOut, Real* out) const;

    /// Calculate out_u = NInv(q)*in_q (e.g., u=NInv*qdot)
    /// or out_q = ~NInv*in_u. Note that one of "in" and "out" is always "q-like" while
    /// the other is "u-like", but which is which changes if the matrix is transposed.
    /// Note that the transposed operation here is the same as multiplying by NInv on
    /// the right, with the Vectors viewed as RowVectors instead.
    /// The default implementation assumes that NInv is an identity matrix, and will
    /// only work if nq=nu=nIn=nOut and nAngles < 4 (i.e., no quaternions). If this
    /// is true for your mobilizer, you do not need to implement this method.
    ///
    /// The state is guaranteed to have been realized to at least Position stage.
    virtual void multiplyByNInv(const State& s, bool transposeMatrix, 
                                int nIn, const Real* in, int nOut, Real* out) const;

    /// Calculate out_q = NDot(q)*in_u
    /// or out_u = ~NDot*in_q. Note that one of "in" and "out" is always "q-like" while
    /// the other is "u-like", but which is which changes if the matrix is transposed.
    /// Note that the transposed operation here is the same as multiplying by NDot on
    /// the right, with the Vectors viewed as RowVectors instead.
    /// The default implementation assumes that NDot is zero, and will
    /// only work if nq=nu=nIn=nOut and nAngles < 4 (i.e., no quaternions) If this
    /// is true for your mobilizer, you do not need to implement this method.
    ///
    /// The state is guaranteed to have been realized to at least Position stage.
    virtual void multiplyByNDot(const State& s, bool transposeMatrix, 
                                int nIn, const Real* in, int nOut, Real* out) const;

        // Methods for setting Mobilizer initial conditions. Note -- I've stripped this
        // down to the two basic routines but the built-ins have 8 so that you can 
        // specify only rotations or translations. I'm not sure that's needed here and
        // I suppose you could add more routines later if needed.
        // Eventually it might be nice to provide default implementation here that would
        // use a root finder to attempt to solve these initial condition problems. For
        // most joints there is a much more direct way to do it, and sometimes there
        // are behavioral choices to make, which is why it is nice to have mobilizer-specific
        // overrides for these.

    /// Find a set of q's for this mobilizer that best approximate the supplied Transform
    /// which requests a particular relative orientation and translation between
    /// the "fixed" and "moving" frames connected by this mobilizer.
    /// The state is guaranteed to have been realized to at least Instance stage.
    ///
    /// The default implementation uses a nonlinear optimizer to search for the best
    /// fit.  Whenever possible, subclasses should override this to provide a faster
    /// and more robust implementation.
    virtual void setQToFitTransform(const State&, const Transform& X_FM, int nq, Real* q) const;

    /// Find a set of u's (generalized speeds) for this mobilizer that best approximate
    /// the supplied spatial velocity \p V_FM which requests the relative angular
    /// and linear velocity between the "fixed" and "moving" frames connected by
    /// this mobilizer. Routines which affect generalized speeds u depend on the generalized
    /// coordinates q already having been set; they never change these coordinates.
    /// The state is guaranteed to have been realized to at least Position stage.
    /// @see setQToFitTransform()
    ///
    /// The default implementation uses a nonlinear optimizer to search for the best
    /// fit.  Whenever possible, subclasses should override this to provide a faster
    /// and more robust implementation.
    virtual void setUToFitVelocity(const State&, const SpatialVec& V_FM, int nu, Real* u) const;

    /// Implement this optional method if you would like your MobilizedBody to generate any suggestions
    /// for geometry that could be used as default visualization as an aid to understanding a system
    /// containing this MobilizedBody. For example, if your mobilizer connects two points, you might
    /// want to draw a line between those points. You can also generate text labels, and you can
    /// provide methods for controlling the presence or appearance of your generated geometry.
    /// If you don't implement this routine no extra geometry will be generated here.
    virtual void calcDecorativeGeometryAndAppend
       (const State& s, Stage stage, std::vector<DecorativeGeometry>& geom) const
    {
    }
    //@}


    /// @name Optional realize() Virtual Methods
    /// Provide implementations of these methods if you want to allocate State variables (such
    /// as modeling options or parameters) or want to pre-calculate some expensive quantities and
    /// store them in the State cache for your future use. Note that the Position and Velocity
    /// realize methods will be called <em>before</em> calling the matrix operator methods
    /// for this MobilizedBody. That way if you want to precalculate the H or HDot matrix,
    /// for example, you can do so in realizePosition() or realizeVelocity() and then use it
    /// in multiplyByHMatrix(), etc.

    //@{
    /// The Matter Subsystem's realizeTopology() method will call this method along with the built-in
    /// MobilizedBodies' realizeTopology() methods. This gives the MobilizedBody a chance to 
    ///   - pre-calculate Topology stage "cache" values (mutable values which are stored
    ///     in the derived Implementation class directly), and
    ///   - allocate Model-stage state variables for later use, and
    ///   - allocate Model-stage cache entries in the State.
    /// The indices to the Model-stage state & cache entries are stored locally as part of 
    /// the Topology-stage cache.
    virtual void realizeTopology(State&) const { }

    /// The Matter Subsystem's realizeModel() method will call this method along with the built-in
    /// MobilizedBodies' realizeModel() methods. This gives the MobilizedBody a chance to 
    ///   - pre-calculate Model stage cache values according to the settings of the Model variables,
    ///   - allocate any later-Stage variables that may be needed (typically these will be 
    ///     Instance stage variables containing geometric information or parameters
    ///     like lengths or pitch for a Screw.
    /// The indices to any of the State entries allocated here are stored in the State as part
    /// of the Model-stage cache.
    virtual void realizeModel(State&) const { }

    /// The Matter Subsystem's realizeInstance() method will call this method along with the built-in
    /// MobilizedBodies' realizeInstance() methods. This gives the MobilizedBody a chance to 
    ///   - pre-calculate Instance stage cache values according to the settings of the Instance variables.
    virtual void realizeInstance(const State&) const { }

    /// The Matter Subsystem's realizeTime() method will call this method along with the built-in
    /// MobilizedBodies' realizeTime() methods. This gives the MobilizedBody a chance to 
    ///   - pre-calculate Time stage cache values according to the current value of time found
    ///     in the State.
    virtual void realizeTime(const State&) const { }

    /// The Matter Subsystem's realizePosition() method will call this method along with the built-in
    /// MobilizedBodies' realizePosition() methods. This gives the MobilizedBody a chance to 
    ///   - pre-calculate Position stage cache values according to the current values of positions found
    ///     in the State.
    /// Note that this is called <em>before</em> methods which implement operators involving position-dependent
    /// matrices N and H.
    virtual void realizePosition(const State&) const { }

    /// The Matter Subsystem's realizeVelocity() method will call this method along with the built-in
    /// MobilizedBodies' realizeVelocity() methods. This gives the MobilizedBody a chance to 
    ///   - pre-calculate Velocity stage cache values according to the current values of velocities found
    ///     in the State.
    /// Note that this is called <em>before</em> methods which implement operators involving velocity-dependent
    /// matrices NDot and HDot.
    virtual void realizeVelocity(const State&) const { }

    /// The Matter Subsystem's realizeDynamics() method will call this method along with the built-in
    /// MobilizedBodies' realizeDynamics() methods. This gives the MobilizedBody a chance to 
    ///   - pre-calculate Dynamics stage cache values according to the current values found
    ///     in the State.
    /// Computations at Dynamics stage cannot affect the behavior of the MobilizedBody since that
    /// is completely determined by the Position and Velocity stage operators.
    virtual void realizeDynamics(const State&) const { }

    /// The Matter Subsystem's realizeAcceleration() method will call this method along with the built-in
    /// MobilizedBodies' realizeAcceleration() methods. This gives the MobilizedBody a chance to 
    ///   - pre-calculate Acceleration stage cache values according to the current values of body
    ///     and mobility accelerations found in the State.
    /// Computations at Acceleration stage cannot affect the behavior of the MobilizedBody since that
    /// is completely determined by the Position and Velocity stage operators.
    virtual void realizeAcceleration(const State&) const { }

    /// The Matter Subsystem's realizeReport() method will call this method along with the built-in
    /// MobilizedBodies' realizeReport() methods. This gives the MobilizedBody a chance to 
    ///   - calculate Report stage cache values according to the current values found
    ///     in the State.
    /// Computations at Report stage cannot affect the progress of a simulation in any way.
    virtual void realizeReport(const State&) const { }
    //@}

    friend class MobilizedBody::CustomImpl;
};

/**
 * This is a subclass of MobilizedBody::Custom which uses a set of Function objects to define the behavior of the
 * MobilizedBody.  When you create it, you specify the number of generalized coordinates, and six Functions which
 * calculate the spatial rotations and translations based on those coordinates.  It assumes there is a one to one
 * correspondence between generalized coordinates and generalized speeds, so qdot == u.
 * 
 * Each of the Function objects must take some subset of the generalized coordinates as inputs, and produce a single
 * number as its output.  It also must support derivatives up to second order.  Taken together, the six Functions
 * define a SpatialVec giving the body's mobilizer transform.
 */

class SimTK_SIMBODY_EXPORT MobilizedBody::FunctionBased : public MobilizedBody::Custom {
public:
    /* Create a FunctionBased MobilizedBody.
     * 
     * @param parent         the MobilizedBody's parent body
     * @param body           describes this MobilizedBody's physical properties
     * @param nmobilities    the number of generalized coordinates belonging to this MobilizedBody
     * @param functions      the Functions describing how the body moves based on its generalized coordinates.
     *                       This must be of length 6.  The elements correspond to, in order, x rotation, y rotation, z rotation,
     *                       x translation, y translation, and z translation.  The MobilizedBody takes over ownership of the functions,
     *                       and automatically deletes them when the MobilizedBody is deleted.
     * @param coordIndices   the indices of the generalized coordinates that are inputs to each function.  For example, if coordIndices[2] = {0, 1},
     *                       that means that functions[2] takes two input arguments, and q[0] and q[1] respectively should be passed as those arguments.
     * @param direction      whether you want the coordinates defined as though parent & child were swapped
     */
    FunctionBased(MobilizedBody& parent, const Body& body, 
                  int nmobilities, const std::vector<const Function*>& functions,
                  const std::vector<std::vector<int> >& coordIndices,
                  Direction direction=Forward);
    /* Create a FunctionBased MobilizedBody.
     * 
     * @param parent         the MobilizedBody's parent body
     * @param inbFrame       the default inboard frame
     * @param body           describes this MobilizedBody's physical properties
     * @param outbFrame      the default outboard frame
     * @param nmobilities    the number of generalized coordinates belonging to this MobilizedBody
     * @param functions      the Functions describing how the body moves based on its generalized coordinates.
     *                       This must be of length 6.  The elements correspond to, in order, x rotation, y rotation, z rotation,
     *                       x translation, y translation, and z translation.  The MobilizedBody takes over ownership of the functions,
     *                       and automatically deletes them when the MobilizedBody is deleted.
     * @param coordIndices   the indices of the generalized coordinates that are inputs to each function.  For example, if coordIndices[2] = {0, 1},
     *                       that means that functions[2] takes two input arguments, and q[0] and q[1] respectively should be passed as those arguments.
     * @param direction      whether you want the coordinates defined as though parent & child were swapped
     */
    FunctionBased(MobilizedBody& parent, const Transform& inbFrame, 
                  const Body& body, const Transform& outbFrame, 
                  int nmobilities, const std::vector<const Function*>& functions,
                  const std::vector<std::vector<int> >& coordIndices,
                  Direction direction=Forward);
    /* Create a FunctionBased MobilizedBody.
     * 
     * @param parent         the MobilizedBody's parent body
     * @param body           describes this MobilizedBody's physical properties
     * @param nmobilities    the number of generalized coordinates belonging to this MobilizedBody
     * @param functions      the Functions describing how the body moves based on its generalized coordinates.
     *                       This must be of length 6.  The elements correspond to, in order, x rotation, y rotation, z rotation,
     *                       x translation, y translation, and z translation.  The MobilizedBody takes over ownership of the functions,
     *                       and automatically deletes them when the MobilizedBody is deleted.
     * @param coordIndices   the indices of the generalized coordinates that are inputs to each function.  For example, if coordIndices[2] = {0, 1},
     *                       that means that functions[2] takes two input arguments, and q[0] and q[1] respectively should be passed as those arguments.
	 * @param axes			 the axes directions (as Vec3's) for each spatial coordinate, which each function describes, and is therefore length 6.
	 *						 First 3 and last 3 axes must be linearly independent, otherwise there will be redundant speeds for the same motion.
     * @param direction      whether you want the coordinates defined as though parent & child were swapped
     */
    FunctionBased(MobilizedBody& parent, const Body& body, 
                  int nmobilities, const std::vector<const Function*>& functions,
                  const std::vector<std::vector<int> >& coordIndices, const std::vector<Vec3>& axes,
                  Direction direction=Forward);
    /* Create a FunctionBased MobilizedBody.
     * 
     * @param parent         the MobilizedBody's parent body
     * @param inbFrame       the default inboard frame
     * @param body           describes this MobilizedBody's physical properties
     * @param outbFrame      the default outboard frame
     * @param nmobilities    the number of generalized coordinates belonging to this MobilizedBody
     * @param functions      the Functions describing how the body moves based on its generalized coordinates.
     *                       This must be of length 6.  The elements correspond to, in order, x rotation, y rotation, z rotation,
     *                       x translation, y translation, and z translation.  The MobilizedBody takes over ownership of the functions,
     *                       and automatically deletes them when the MobilizedBody is deleted.
     * @param coordIndices   the indices of the generalized coordinates that are inputs to each function.  For example, if coordIndices[2] = {0, 1},
     *                       that means that functions[2] takes two input arguments, and q[0] and q[1] respectively should be passed as those arguments.
     * @param axes			 the axes directions (as Vec3's) for each spatial coordinate, which each function describes, and is therefore length 6.
	 *						 First 3 and last 3 axes must be linearly independent, otherwise there will be redundant speeds for the same motion.
     * @param direction      whether you want the coordinates defined as though parent & child were swapped
	 */
    FunctionBased(MobilizedBody& parent, const Transform& inbFrame, 
                  const Body& body, const Transform& outbFrame, 
                  int nmobilities, const std::vector<const Function*>& functions,
                  const std::vector<std::vector<int> >& coordIndices, const std::vector<Vec3>& axes,
                  Direction direction=Forward);
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_H_



