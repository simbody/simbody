#ifndef SimTK_SIMBODY_MOBILIZED_BODY_CUSTOM_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_CUSTOM_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman, Michael Sherman                                    *
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

/** @file
Declares the MobilizedBody::Custom and MobilizedBody::Custom::Implementation
subclasses. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

//==============================================================================
//                        MOBILIZED BODY :: CUSTOM  
//==============================================================================
/** The handle class MobilizedBody::Custom (dataless) and its companion class 
MobilizedBody::Custom::Implementation can be used together to define new 
MobilizedBody types with arbitrary properties. To use it, create a class that 
extends MobilizedBody::Custom::Implementation. You can then create an instance 
of it and pass it to the MobilizedBody::Custom constructor:

@code
MobilizedBody::Custom myMobilizedBody(new MyMobilizedBodyImplementation(args));
@endcode
("args" here and below stands for whatever arguments are needed for your
particular mobilizer; it isn't meant literally.)

Alternatively, you can also create a new Handle class which is a subclass of 
MobilizedBody::Custom and which creates the Implementation itself in its 
constructors.
@code
class MyMobilizedBody : public MobilizedBody::Custom {
public:
  MyMobilizedBody(args) 
  :   MobilizedBody::Custom(new MyMobilizedBodyImplementation(args)) {}
};
@endcode

This allows an end user to simply write
@code
MyMobilizedBody(args);
@endcode

and not worry about implementation classes or creating objects on the heap.  
If you do this, your MobilizedBody::Custom subclass must not have any data 
members or virtual methods. If it does, it will not work correctly. Instead,
store all data in the Implementation subclass.

@see MobilizedBody::Custom::Implementation **/
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

    /** @cond **/ // hide from Doxygen
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Custom, CustomImpl, MobilizedBody);
    /** @endcond **/
protected:
    const Implementation& getImplementation() const;
    Implementation&       updImplementation();

    Custom() {}
};

// We only want the template instantiation to occur once. This symbol is 
// defined in the Simbody compilation unit that defines the MobilizedBody class
// but should not be defined any other time.
#ifndef SimTK_SIMBODY_DEFINING_MOBILIZED_BODY
    extern template class PIMPLHandle<MobilizedBody::Custom::Implementation, 
                                      MobilizedBody::Custom::ImplementationImpl>;
#endif


//==============================================================================
//                  MOBILIZED BODY :: CUSTOM :: IMPLEMENTATION
//==============================================================================
/** This is the implementation class for Custom mobilizers.
@see MobilizedBody::Custom **/
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
    /// NOTE: if you don't say there are any angles, you can manage things yourself.
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
    /// Caution: if your mobilizer has a quaternion, the four q's will not necessarily be
    /// normalized here but you \e must normalize them before converting them to a Rotation.
    /// Casting them to a Quaternion will do that automatically.
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
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const
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

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_CUSTOM_H_



