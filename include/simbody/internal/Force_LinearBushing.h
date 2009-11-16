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


/** @file
 * This contains the user-visible API ("handle" class) for the SimTK::Force 
 * subclass Force::LinearBushing and is logically part of Force.h. The file 
 * assumes that Force.h will have included all necessary declarations.
 */

/**
 * This force element represents a bushing with linear material properties, 
 * acting to connect two frames F ("fixed") and M ("moving"), with F fixed to 
 * body A and M fixed to body B. The relative orientation R_FM is parametrized 
 * by an x-y-z (1-2-3) B-fixed Euler angle sequence and its time derivatives, as 
 * well as a position vector p_FM from OF to OM expressed in F, and its time 
 * derivative v_FM (taken in F). For small orientation displacements, the Euler 
 * angles can be considered independent rotations about x, y, and z. Stiffness 
 * and damping parameters (6 of each) are provided for each direction's rotation 
 * and translation. The 6 coordinates q are defined as the three Euler angles 
 * [qx,qy,qz] followed by the three translations [px,py,pz]=p_FM.
 * 
 * The generated force is calculated for each of the 6 coordinates as
 * f_i = -(k_i*q_i + c_i*qdot_i). Each contribution to potential energy is
 * e_i = k_i*q_i^2/2 (damping term contributes no potential energy).
 * The scalar rotational moments f_0, f_1, and f_2 act about rotated axes so
 * do not constitute a vector; the are transformed internally here to produce 
 * the appropriate moments on the bodies. The scalar translational forces
 * f_3, f_4, f_5 on the other hand are aligned with frame F's axes so 
 * constitute a vector in F.
 * 
 * This force element is intended for small-displacement use, but is defined
 * nonlinearly so physically correct for large displacements also. However, be
 * aware that the q's cannot "wrap" so you must keep the motion small enough
 * that the set of q's inferred from the transform X_FM (M's configuration in
 * F) is continuous during a simulation. Also, the Bushing is singular in a 
 * configuration in which the middle rotation angle is near 90 degrees, because 
 * the time derivative of that angle is unbounded; you should stay far away
 * from that configuration.
 *
 * If you would like a force element like this one but suited for very large 
 * rotations (e.g., multiple revolutions) you should arrange to have a Bushing
 * Mobilizer connected between frames F and M so that you can apply mobility
 * forces rather than body forces as we're doing here. In that case the q's are
 * the defining coordinates rather than the frames, and mobilizer q's can wrap
 * continuously without limit.
 */

class SimTK_SIMBODY_EXPORT Force::LinearBushing : public Force {
public:
    /**
     * Create a LinearBushing between arbitrary frames fixed to two bodies.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param bodyA      the first body to which the force should be applied
     * @param frameFOnA  a frame F fixed to body A given by its constant
     *                      transform X_AF from the A frame
     * @param bodyB      the other body to which the force should be applied
     * @param frameMOnB  a frame M fixed to body B given by its constant
     *                      transform X_BM from the B frame
     * @param stiffness  the six spring constants, torsional followed by 
     *                      translational
     * @param damping    the six damping coefficients, torsional followed by
     *                      translational
     */
    LinearBushing(GeneralForceSubsystem& forces, 
                  const MobilizedBody& bodyA, const Transform& frameFOnA, 
                  const MobilizedBody& bodyB, const Transform& frameMOnB, 
                  const Vec6& stiffness, const Vec6& damping);

    /**
     * Create a LinearBushing connecting the body frames of two bodies.
     * This is the same as the more general constructor except it assumes
     * identity transforms for the two frames, meaning they are coincident
     * with the body frames.
     * 
     * @param forces     the subsystem to which this force should be added
     * @param bodyA      the first body to which the force should be applied
     * @param bodyB      the other body to which the force should be applied
     * @param stiffness  the six spring constants, torsional followed by 
     *                      translational
     * @param damping    the six damping coefficients, torsional followed by
     *                      translational
     */
    LinearBushing(GeneralForceSubsystem& forces, 
                  const MobilizedBody& bodyA, 
                  const MobilizedBody& bodyB, 
                  const Vec6& stiffness, const Vec6& damping);

    /**
     * Obtain the generalized coordinates last calculated by this force
     * element. These are the bodyB-fixed x-y-z Euler angles qx,qy,qz and
     * p_FM=[px,py,pz], the vector from the frame F origin OF to
     * the frame M origin OM, expressed in the F frame. The full generalized
     * coordinate vector q=[qx,qy,qz,px,py,pz]. You can call this 
     * after the supplied State has been realized to Position stage.
     *
     * @param state     the State from whose cache the coordinates are retrieved
     * @return          the six generalized coordinates relating the frames
     */
    const Vec6& getQ(const State& state) const;

    /**
     * Obtain the generalized coordinate derivatives last calculated by this 
     * force element. These are the bodyB-fixed x-y-z Euler angle derivatives
     * qdotx,qdoty,qdotz and v_FM=[vx,vy,vz], the velocity of point OM in 
     * frame F, expressed in F. That is, v_FM = d/dt p_FM with the derivative
     * taken in the F frame. The full generalized coordinate derivative vector 
     * qdot=[qdotx,qdoty,qdotz,vx,vy,vz]. You can call this after the supplied 
     * State has been realized to Velocity stage.
     *
     * @param state     the State from whose cache the coordinate derivatives
     *                      are retrieved
     * @return          the six generalized coordinate derivatives
     * @see getQ()
     */
    const Vec6& getQDot(const State& state) const;

    /**
     * Obtain the spatial transform X_FM giving the location and orientation of
     * body B's frame M in body A's frame F. You can call this after the 
     * supplied State has been realized to Position stage.
     *
     * @param state     the State from whose cache the transform is retrieved
     * @return          the transform X_FM
     */
    const Transform& getX_FM(const State& state) const;

    /**
     * Obtain the spatial velocity V_FM giving the velocity of body B's frame 
     * M in body A's frame F, expressed in the F frame. Note that this is 
     * the time derivative of X_FM <em>taken in F</em>, which means that V_FM
     * is a local relationship between F and M and is not affected by F's 
     * motion with respect to Ground. You can call this after the supplied 
     * State has been realized to Velocity stage.
     *
     * @param state     the State from whose cache the velocity is retrieved
     * @return          the spatial velocity V_FM as a spatial vector
     */
    const SpatialVec& getV_FM(const State& state) const;

    /**
     * Obtain the generalized forces f begin applied by this Bushing force
     * element on each of its six axes. The sign is such that it would be 
     * appropriate to apply to a Bushing mobilizer connecting the same two
     * frames; that is, these are the generalized forces acting on the
     * "outboard" body B; the negative of these forces acts on body A. You can 
     * call this after the supplied State has been realized to Velocity stage.
     *
     * @param state     the State from whose cache the forces are retrieved
     * @return          the six generalized forces as a Vec6 in the order
     *                      mx,my,mz,fx,fy,fz where the m's are moments
     *                      along the rotated Euler axes
     */
    const Vec6& getF(const State& state) const;

    /**
     * Return the spatial force F_GM being applied by this Bushing to body B
     * at its M frame, expressed in the Ground frame. You can call this after 
     * the supplied State has been realized to Velocity stage.
     *
     * @param state     the State from whose cache the force is retrieved
     * @return          the spatial force F_GM as a spatial vector, expressed
     *                      in G
     */
    const SpatialVec& getF_GM(const State& state) const;

    /**
     * Return the spatial force F_GF being applied by this Bushing to body A
     * at its F frame, expressed in the Ground frame. This is equal and 
     * opposite to the spatial force applied on body B.
     *
     * @param state     the State from whose cache the force is retrieved
     * @return          the spatial force F_GF as a spatial vector, expressed
     *                      in G
     */
    const SpatialVec& getF_GF(const State& state) const;

    /**
     * Obtain the potential energy stored in this Bushing in the current
     * configuration. You can call this after the supplied State has been 
     * realized to Position stage.
     *
     * @param state     the State from whose cache the potential energy is
     *                      retrieved
     * @return          the potential energy (a scalar)
     */
    const Real& getPotentialEnergy(const State& state) const;

    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(LinearBushing, LinearBushingImpl, Force);
};


#endif // SimTK_SIMBODY_FORCE_LINEAR_BUSHING_H_
