#ifndef SimTK_SIMBODY_FORCE_LINEAR_BUSHING_H_
#define SimTK_SIMBODY_FORCE_LINEAR_BUSHING_H_
/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2009 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
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

#include "SimTKcommon.h"
#include "simbody/internal/Force.h"

/** @file
 * This contains the user-visible API ("handle" class) for the SimTK::Force 
 * subclass Force::LinearBushing and is logically part of Force.h. The file 
 * assumes that Force.h will have included all necessary declarations.
 */

namespace SimTK {
/**
 * This force element represents a bushing acting to connect a frame F fixed 
 * on one body (B1) to a frame M fixed on a second body (B2), by a massless,
 * compliant element with linear stiffness and damping properties.
 * 
 * The relative orientation R_FM of frame M given in frame F is parametrized by
 * an x-y-z (1-2-3) B2-fixed Euler angle sequence and their time derivatives, 
 * as well as a position vector p_FM from F's origin OF to M's origin OM 
 * expressed in F, and its time derivative (velocity) v_FM (taken in F). For 
 * small orientation displacements, the Euler angles can be considered 
 * independent rotations about x, y, and z. Stiffness and damping parameters 
 * (6 of each) are provided for each direction's rotation and translation. The
 * 6 coordinates q are defined as the three Euler angles [qx,qy,qz] followed 
 * by the three translations [px,py,pz]=p_FM.
 * 
 * This force element is intended for small-displacement use, but is defined
 * nonlinearly so is physically correct for large displacements also. However,
 * be aware that the inferred q's cannot "wrap" so you must keep the motion 
 * small enough that the set of q's inferred from the transform X_FM (M's 
 * configuration in F) is continuous during a simulation. Also, the Bushing is 
 * singular in a configuration in which the middle rotation angle is near 90 
 * degrees, because the time derivative of that angle is unbounded; you should 
 * stay far away from that configuration.
 *
 * If you would like a force element like this one but suited for very large 
 * rotations (e.g., multiple revolutions along one of the axes) you should 
 * arrange to have a Bushing Mobilizer connected between frames F and M so that
 * you can apply mobility forces rather than body forces as we're doing here. 
 * In that case the q's are the defining coordinates rather than the frames, 
 * and mobilizer Euler angle q's can wrap continuously without limit. However,
 * with Euler angles there is necessarily a singular configuration so one of
 * the axes must be kept within a modest range even if you use mobilizer
 * coordinates.
 * 
 * \par Theory:
 *
 * The LinearBushing calculates the relative orientation X_FM between the 
 * frames, and the relative spatial velocity V_FM. From these it infers a set
 * of q's and qdot's as though there were a gimbal-like mechanism connecting
 * frame F to frame M. This is done using standard formulas to decompose
 * a rotation matrix into an Euler sequence (1-2-3 body fixed in this case)
 * and to convert angular velocity to the time derivatives of the Euler
 * angles. (There are also translational q's but they are trivial.)
 *
 * The generated force is then calculated for each of the 6 inferred 
 * coordinates and their 6 time derivatives as
 * <pre>    f_i = -(k_i*q_i + c_i*qdot_i).  </pre>
 * Each contribution to potential energy is
 * <pre>    e_i = k_i*q_i^2/2.              </pre>
 * The damping terms contribute no potential energy, but dissipate power at 
 * a rate 
 * <pre>    p_i = c_i*qdot_i^2.             </pre>
 * Numerically integrating the dissipated power over time gives the energy 
 * dissipated by the Bushing since some arbitrary starting time and can be used
 * to check conservation of energy in the presence of damping. The 
 * LinearBushing force allocates a state variable for this purpose so tracks
 * the energy it dissipates; see the getDissipatedEnergy() method.
 *
 * The scalar rotational moments f_0, f_1, and f_2 act about rotated axes so
 * do not constitute a vector; they are transformed internally here to produce 
 * the appropriate moments on the bodies. The scalar translational forces
 * f_3, f_4, f_5 on the other hand are aligned with frame F's axes so 
 * constitute a vector in F.
 */
class SimTK_SIMBODY_EXPORT Force::LinearBushing : public Force {
public:
    /// Create a LinearBushing between a pair of arbitrary frames, one fixed 
    /// to each of two bodies.
    /// 
    /// @param[in,out]      forces       
    ///     The subsystem to which this force should be added.
    /// @param[in]          body1        
    ///     The first body to which the force should be applied.
    /// @param[in]          frameFOnB1   
    ///     A frame F fixed to body 1 given by its constant transform X_B1F 
    ///     from the B1 frame.
    /// @param[in]          body2        
    ///     The other body to which the force should be applied.
    /// @param[in]          frameMOnB2   
    ///     A frame M fixed to body 2 given by its constant transform X_B2M 
    ///     from the B2 frame.
    /// @param[in]          stiffness    
    ///     The six nonnegative spring constants, torsional followed by 
    ///     translational.
    /// @param[in]          damping      
    ///     The six nonnegative damping coefficients, torsional followed by
    ///     translational.
    LinearBushing(GeneralForceSubsystem& forces, 
                  const MobilizedBody& body1, const Transform& frameFOnB1, 
                  const MobilizedBody& body2, const Transform& frameMOnB2, 
                  const Vec6& stiffness, const Vec6& damping);

    /// Create a LinearBushing connecting the body frames of two bodies.
    /// This is the same as the more general constructor except it assumes
    /// identity transforms for the two frames, meaning they are coincident
    /// with the body frames.
    /// 
    /// @param[in,out]      forces
    ///     The subsystem to which this force should be added.
    /// @param[in]          body1
    ///     The first body to which the force should be applied.
    /// @param[in]          body2
    ///     The other body to which the force should be applied.
    /// @param[in]          stiffness    
    ///     The six nonnegative spring constants, torsional followed by 
    ///     translational.
    /// @param[in]          damping      
    ///     The six nonnegative damping coefficients, torsional followed by
    ///     translational.
    LinearBushing(GeneralForceSubsystem& forces, 
                  const MobilizedBody& body1, 
                  const MobilizedBody& body2, 
                  const Vec6& stiffness, const Vec6& damping);



    //--------------------------------------------------------------------------
    /// @name Default Parameters
    ///
    /// These refer to Topology-stage parameters normally set in the
    /// constructor; these determine the initial values assigned to the 
    /// corresponding Instance-stage state variables.
    ///
    /// \par Notes
    /// - Changing one of these parameters invalidates the containing System's 
    ///   topology, meaning that realizeTopology() will have to be called 
    ///   before subsequent use. 
    /// - The set...() methods return a reference to "this" LinearBushing (in
    ///   the manner of an assignment operator) so they can be chained in a 
    ///   single expression.
    //@{
    /// Set the frame F on body B1 that will be used by default as one of the
    /// frames connected by the bushing, by giving the transform X_B1F locating 
    /// frame F with respect to the body frame B1.
    /// @param[in]          X_B1F    
    ///     The Transform giving the default placement of the first Bushing 
    ///     frame F in the body 1 frame B1.
    /// @return
    ///     A writable reference to "this" Linear Bushing which will now have
    ///     the new default F frame
    LinearBushing& setDefaultFrameOnBody1(const Transform& X_B1F);
    /// Set the frame M on body B2 that will be used by default as one of the
    /// frames connected by the bushing, by giving the transform X_B2M locating 
    /// frame M with respect to the body frame B2.
    /// @param[in]          X_B2M    
    ///     The Transform giving the default placement of the second Bushing 
    ///     frame M in the body 2 frame B2.
    /// @return
    ///     A writable reference to "this" Linear Bushing which will now have 
    ///     the new default M frame.
    LinearBushing& setDefaultFrameOnBody2(const Transform& X_B2M);
    /// Set the stiffnesses (spring constants) k that will be used by default 
    /// for this Bushing; these must be nonnegative.
    /// @param[in]          stiffness    
    ///     The six nonnegative spring constants, one for each axis in the
    ///     order qx,qy,qz,px,py,pz (i.e., rotations first).
    /// @return 
    ///     A writable reference to "this" Linear Bushing which will now have 
    ///     the new default stiffnesses.
    LinearBushing& setDefaultStiffness(const Vec6& stiffness);
    /// Set the damping coefficients c that will be used by default for this 
    /// Bushing; these must be nonnegative.
    /// @param[in]          damping
    ///     The six nonnegative damping coefficients, one for each axis in the
    ///     order qx,qy,qz,px,py,pz (i.e., rotations first).
    /// @return
    ///     A writable reference to "this" Linear Bushing which will now have 
    ///     the new default damping coefficients
    LinearBushing& setDefaultDamping(const Vec6& damping);

    /// Return the frame F on body B1 that will be used by default as one of the
    /// frames connected by the bushing, as the transform X_B1F locating frame F
    /// with respect to the body frame B1.
    /// @return
    ///     The Transform X_B1F giving the default placement of the first 
    ///     Bushing frame F in the body 1 frame B1.
    const Transform& getDefaultFrameOnBody1() const;
    /// Return the frame M on body B2 that will be used by default as one of the
    /// frames connected by the bushing, as the transform X_B2M locating frame M
    /// with respect to the body frame B2.
    /// @return
    ///     The Transform X_B2M giving the default placement of the second 
    ///     Bushing frame M in the body 2 frame B2.
    const Transform& getDefaultFrameOnBody2() const;
    /// Return the stiffnesses k that will be used by default for this Bushing.
    /// @return
    ///     The six default spring constants, one for each axis in the
    ///     order qx,qy,qz,px,py,pz (i.e., rotations first).
    const Vec6& getDefaultStiffness() const;
    /// Return the damping coefficients c that will be used by default for this 
    /// Bushing.
    /// @return
    ///     The six default damping coefficients, one for each axis in the
    ///     order qx,qy,qz,px,py,pz (i.e., rotations first).
    const Vec6& getDefaultDamping()   const;
    //@}......................... Default Parameters ...........................



    //--------------------------------------------------------------------------
    /// @name Instance Parameters
    ///
    /// These refer to the Instance-stage state variables that determine the
    /// geometry and material properties that will be used in computations
    /// involving this Bushing when performed with the given State. If these
    /// are not set explicitly, the default values are set to those provided in
    /// the constructor or via the correponding setDefault...() methods.
    ///
    /// \par Notes
    /// - Changing one of these parameters invalidates the given State's 
    ///   Instance stage, meaning that realize(Instance) will have to be called 
    ///   (explicitly or implicitly by realizing a higher stage) before 
    ///   subsequent use.
    /// - The set...() methods here return a const reference to "this" 
    ///   LinearBushing (in the manner of an assignment operator, except
    ///   read-only) so they can be chained in a single expression. 
    //@{
    /// Set the frame F on body B1 that will be used as one of the frames 
    /// connected by the bushing when evaluating with this State, by giving 
    /// the transform X_B1F locating frame F with respect to the body frame B1.
    /// @pre \a state realized to Stage::Topology
    /// @param[in,out]      state
    ///     The State object that is modified by this method.
    /// @param[in]          X_B1F
    ///     The Transform giving the placement of the first Bushing frame F 
    ///     in the body 1 frame B1.
    /// @return
    ///     A const reference to "this" Linear Bushing for convenient
    ///     chaining of set...() methods in a single expression.
    const LinearBushing& setFrameOnBody1(State&             state, 
                                         const Transform&   X_B1F) const;
    /// Set the frame M on body B2 that will be used as one of the frames 
    /// connected by the bushing when evaluating with this State, by giving 
    /// the transform X_B2M locating frame M with respect to the body frame B2.
    /// @pre \a state realized to Stage::Topology
    /// @param[in,out]      state
    ///     The State object that is modified by this method.
    /// @param[in]          X_B2M
    ///     The Transform giving the placement of the second Bushing frame M 
    ///     in the body 2 frame B2.
    /// @return
    ///     A const reference to "this" Linear Bushing for convenient
    ///     chaining of set...() methods in a single expression.
    const LinearBushing& setFrameOnBody2(State&             state, 
                                         const Transform&   X_B2M) const;
    /// Set the stiffnesses (spring constants) k that will be used for this 
    /// Bushing when evaluated using this State.
    /// @pre \a state realized to Stage::Topology
    /// @param[in,out]      state
    ///     The State object that is modified by this method.
    /// @param[in]          stiffness    
    ///     Six nonnegative spring constants, one for each axis in the
    ///     order qx,qy,qz,px,py,pz (i.e., rotations first).
    /// @return
    ///     A const reference to "this" Linear Bushing for convenient
    ///     chaining of set...() methods in a single expression.
    const LinearBushing& setStiffness(State&        state, 
                                      const Vec6&   stiffness) const;
    /// Set the damping coefficients c that will be used for this Bushing when 
    /// evaluated using this State.
    /// @pre \a state realized to Stage::Topology
    /// @param[in,out]      state
    ///     The State object that is modified by this method.
    /// @param[in]          damping
    ///     Six nonnegative damping coefficients, one for each axis in the
    ///     order qx,qy,qz,px,py,pz (i.e., rotations first).
    /// @return
    ///     A const reference to "this" Linear Bushing for convenient
    ///     chaining of set...() methods in a single expression.
    const LinearBushing& setDamping(State&      state, 
                                    const Vec6& damping) const;

    /// Return the frame F on body B1 currently being used by this State as one 
    /// of the frames connected by the bushing, as the transform X_B1F locating 
    /// frame F with respect to the body frame B1.
    /// @pre \a state realized to Stage::Topology
    /// @param[in]          state
    ///     The State object from which we obtain the frame.
    /// @return
    ///     The Transform X_B1F giving the given state's current placement 
    ///     of the first Bushing frame F in the body 1 frame B1.
    const Transform& getFrameOnBody1(const State& state) const;
    /// Return the frame M on body B2 currently being used by this State as one 
    /// of the frames connected by the bushing, as the transform X_B2M locating 
    /// frame M with respect to the body frame B2.
    /// @pre \a state realized to Stage::Topology
    /// @param[in]          state
    ///     The State object from which we obtain the frame.
    /// @return
    ///     The Transform X_B2M giving the given state's current placement 
    ///     of the second Bushing frame M in the body 2 frame B2.
    const Transform& getFrameOnBody2(const State& state) const;
    /// Return the stiffnesses (spring constants) k currently being used for 
    /// this Bushing by this State. 
    /// @pre \a state realized to Stage::Topology
    /// @param[in]          state
    ///     The State object from which we obtain the stiffnesses.
    /// @return
    ///     The current values in the given state of the six spring constants
    ///     k, one for each axis in the order qx,qy,qz,px,py,pz (i.e., 
    ///     rotations first).
    const Vec6& getStiffness(const State& state) const;
    /// Return the damping coefficients c currently being used for this Bushing 
    /// by this State. 
    /// @pre \a state realized to Stage::Topology
    /// @param[in]          state
    ///     The State object from which we obtain the damping coefficients.
    /// @return
    ///     The current values in the given state of the six damping
    ///     coefficients c, one for each axis in the order qx,qy,qz,px,py,pz 
    ///     (i.e., rotations first).
    const Vec6& getDamping(const State& state)   const;
    //@}...................... Instance Parameters .............................



    //--------------------------------------------------------------------------
    /// @name Position-related Quantities
    ///
    /// These methods return position-dependent quantities that are calculated 
    /// by the Bushing as part of its force calculations and stored in the 
    /// State cache. These can be obtained at no cost, although the first call
    /// after a position change may initiate computation if forces haven't 
    /// already been computed.
    ///
    /// @pre These methods may be called only after the supplied state has 
    ///      been realized to Position stage.
    //@{
    /// Obtain the generalized coordinates last calculated by this force
    /// element. These are the body2-fixed x-y-z Euler angles qx,qy,qz and
    /// p_FM=[px,py,pz], the vector from the frame F origin OF to
    /// the frame M origin OM, expressed in the F frame. The full generalized
    /// coordinate vector is q=[qx,qy,qz,px,py,pz].
    /// @pre \a state realized to Stage::Position
    /// @param[in]          state    
    ///     The State from whose cache the coordinates are retrieved.
    /// @return
    ///     The six generalized coordinates relating the frames in order
    ///     qx,qy,qz,px,py,pz (i.e., rotational coordinates first).
    /// @see getQDot()
    const Vec6& getQ(const State& state) const;

    /// Obtain the spatial transform X_GF giving the location and orientation of
    /// body 1's frame F in Ground.
    /// @pre \a state realized to Stage::Position
    /// @param[in]          state
    ///     The State from whose cache the transform is retrieved.
    /// @return
    ///     The transform X_GF.
    const Transform& getX_GF(const State& state) const;

    /// Obtain the spatial transform X_GM giving the location and orientation of
    /// body 2's frame M in Ground.
    /// @pre \a state realized to Stage::Position
    /// @param[in]          state    
    ///     The State from whose cache the transform is retrieved.
    /// @return
    ///     The transform X_GM.
    const Transform& getX_GM(const State& state) const;

    /// Obtain the spatial transform X_FM giving the location and orientation of
    /// body 2's frame M in body 1's frame F.
    /// @pre \a state realized to Stage::Position
    /// @param[in]          state    
    ///     The State from whose cache the transform is retrieved.
    /// @return
    ///     The transform X_FM.
    const Transform& getX_FM(const State& state) const;
    //@}.................. Position-related quantities .........................



    //--------------------------------------------------------------------------
    /// @name Velocity-related Quantities
    ///
    /// These methods return velocity-dependent quantities that are calculated 
    /// by the Bushing as part of its force calculations and stored in the 
    /// State cache. These can be obtained at no cost, although the first call
    /// after a velocity change may initiate computation if forces haven't 
    /// already been computed.
    ///
    /// @pre These methods may be called only after the supplied state has 
    ///      been realized to Velocity stage.
    //@{
    /// Obtain the generalized coordinate derivatives last calculated by this 
    /// force element. These are the bodyB-fixed x-y-z Euler angle derivatives
    /// qdotx,qdoty,qdotz and v_FM=[vx,vy,vz], the velocity of point OM in 
    /// frame F, expressed in F. That is, v_FM = d/dt p_FM with the derivative
    /// taken in the F frame. The full generalized coordinate derivative vector 
    /// qdot=[qdotx,qdoty,qdotz,vx,vy,vz].
    /// @pre \a state realized to Stage::Velocity
    /// @param[in]          state
    ///     The State from whose cache the coordinate derivatives are retrieved.
    /// @return
    ///     The six generalized coordinate derivatives.
    /// @see getQ()
    const Vec6& getQDot(const State& state) const;

    /// Obtain the spatial velocity V_GF giving the velocity of body 1's frame 
    /// F in the Ground frame. Note that this is the time derivative of X_GF 
    /// <em>taken in G</em>, that is, the ordinary spatial velocity.
    /// @pre \a state realized to Stage::Velocity
    /// @param[in]          state    
    ///     The State from whose cache the velocity is retrieved.
    /// @return
    ///     The spatial velocity V_GF as a spatial vector, that is, V_GF[0]
    ///     is the angular velocity of F in G; V_GF[1] is the linear velocity
    ///     of F's origin point OF in G. Both vectors are expressed in G.
    /// @see getV_FM() to get the local bushing deformation velocity
    const SpatialVec& getV_GF(const State& state) const;

    /// Obtain the spatial velocity V_GM giving the velocity of body 2's frame 
    /// M in the Ground frame. Note that this is the time derivative of X_GM 
    /// <em>taken in G</em>, that is, the ordinary spatial velocity.
    /// @pre \a state realized to Stage::Velocity
    /// @param[in]          state
    ///     The State from whose cache the velocity is retrieved.
    /// @return
    ///     The spatial velocity V_GM as a spatial vector, that is, V_GM[0]
    ///     is the angular velocity of M in G; V_GM[1] is the linear velocity
    ///     of M's origin point OM in G. Both vectors are expressed in G.
    /// @see getV_FM() to get the local bushing deformation velocity
    const SpatialVec& getV_GM(const State& state) const;

    /// Obtain the spatial velocity V_FM giving the velocity of body 2's frame 
    /// M in body 1's frame F, expressed in the F frame. Note that this is 
    /// the time derivative of X_FM <em>taken in F</em>, which means that V_FM
    /// is a local relationship between F and M and is not affected by F's 
    /// motion with respect to Ground.
    /// @pre \a state realized to Stage::Velocity
    /// @param[in]          state    
    ///     The State from whose cache the velocity is retrieved.
    /// @return 
    ///     The spatial velocity V_FM as a spatial vector, that is, V_FM[0]
    ///     is the angular velocity of F in M; V_FM[1] is the linear velocity
    ///     of M's origin point OM in F. Both vectors are expressed in F.
    const SpatialVec& getV_FM(const State& state) const;
    //@}.................. Velocity-related quantities .........................



    //--------------------------------------------------------------------------
    /// @name Forces
    ///
    /// These methods return the forces being applied by this Bushing in the
    /// configuration and velocities contained in the supplied State. These are
    /// evaluated during force calculation and available at no computatinal cost
    /// afterwards, although the first call after a velocity change may initiate 
    /// computation if forces haven't already been computed.
    ///
    /// @pre These methods may be called only after the supplied state has 
    ///      been realized to Velocity stage.
    //@{
    /// Obtain the generalized forces f being applied by this Bushing force
    /// element on each of its six axes. The sign is such that it would be 
    /// appropriate to apply to a Bushing mobilizer connecting the same two
    /// frames; that is, these are the generalized forces acting on the
    /// "outboard" body 2; the negative of these forces acts on the "inboard"
    /// body 1.
    /// @pre \a state realized to Stage::Velocity
    /// @param[in]          state
    ///     The State from whose cache the forces are retrieved.
    /// @return
    ///     The six generalized forces as a Vec6 in the order mx,my,mz,
    ///     fx,fy,fz where the m's are moments about the rotated Euler axes
    ///     and the f's are forces along the translational axes of the F frame. 
    const Vec6& getF(const State& state) const;

    /// Return the spatial force F_GM being applied by this Bushing to body 2
    /// at its M frame, expressed in the Ground frame.
    /// @pre \a state realized to Stage::Velocity
    /// @param[in]          state    
    ///     The State from whose cache the force is retrieved.
    /// @return
    ///     The spatial force F_GM as a spatial vector, expressed in G. That
    ///     is, F_GM[0] is the torque vector applied to body 2; F_GM[1] is 
    ///     the force vector applied to body 2 at frame M's origin OM.
    const SpatialVec& getF_GM(const State& state) const;

    /// Return the spatial force F_GF being applied by this Bushing to body 1
    /// at its F frame, expressed in the Ground frame. This is equal and 
    /// opposite to the spatial force applied on body 2.
    /// @pre \a state realized to Stage::Velocity
    /// @param[in]          state
    ///     The State from whose cache the force is retrieved.
    /// @return
    ///     The spatial force F_GF as a spatial vector, expressed in G. That
    ///     is, F_GF[0] is the torque vector applied to body 1; F_GF[1] is 
    ///     the force vector applied to body 1 at frame F's origin OF.
    const SpatialVec& getF_GF(const State& state) const;
    //@}............................ Forces ....................................



    //--------------------------------------------------------------------------
    /// @name Energy, Power, and Work
    ///
    /// These methods return the energy, power, and work-related quantities
    /// associated with this Bushing element for the values in the supplied
    /// State.
    //@{
    /// Obtain the potential energy stored in this Bushing in the current
    /// configuration.
    /// @pre \a state realized to Stage::Position
    /// @param[in]          state    
    ///     The State from whose cache the potential energy is retrieved.
    /// @return
    ///     The potential energy currently contained in the Bushing in the
    ///     configuration contained in \a state (a nonnegative scalar).
    Real getPotentialEnergy(const State& state) const;

    /// Obtain the rate at which energy is being dissipated by this Bushing,
    /// that is, the power being lost. This is in units of energy/time which
    /// is watts in MKS. 
    /// @pre \a state realized to Stage::Velocity
    /// @param[in]          state    
    ///     The State from which to obtain the current value of the power 
    ///     dissipation.
    /// @return
    ///     The dissipated power (a nonnegative scalar).
    /// @see getDissipatedEnergy() for the time-integrated power loss
    Real getPowerDissipation(const State& state) const;

    /// Obtain the total amount of energy dissipated by this Bushing since some
    /// arbitrary starting point. This is the time integral of the power
    /// dissipation. For a system whose only non-conservative forces are 
    /// Bushings, the sum of potential, kinetic, and dissipated energies should 
    /// be conserved. This is a State variable so you can obtain its value any 
    /// time after it is allocated.
    /// @pre \a state realized to Stage::Model
    /// @param[in]          state    
    ///     The State from which to obtain the current value of the dissipated
    ///     energy.
    /// @return
    ///     The total dissipated energy (a nonnegative scalar).
    /// @see getPowerDissipation() for the instantaneous power loss
    Real getDissipatedEnergy(const State& state) const;

    /// Set the accumulated dissipated energy to an arbitrary value. Typically
    /// this is used only to reset the dissipated energy to zero, but non-zero
    /// values can be useful if you are trying to match some existing data or
    /// continuing a simulation. This is a State variable so you can set its 
    /// value any time after it is allocated.
    /// @pre \a state realized to Stage::Model
    /// @param[in,out]      state    
    ///     The State whose dissipated energy variable for this Bushing is set
    ///     to zero.
    /// @param[in]          energy   
    ///     The new value for the accumulated dissipated energy (must be a 
    ///     nonnegative scalar).
    void setDissipatedEnergy(State& state, Real energy) const;
    //@}..................... Energy, Work, and Power ..........................

    // Don't show this in Doxygen.
    /// @cond
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS
       (LinearBushing, LinearBushingImpl, Force);
    /// @endcond
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_LINEAR_BUSHING_H_
