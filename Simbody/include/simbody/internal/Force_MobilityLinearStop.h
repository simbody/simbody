#ifndef SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_STOP_H_
#define SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_STOP_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
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

#include "SimTKcommon.h"
#include "simbody/internal/Force.h"

/** @file
This contains the user-visible API ("handle" class) for the SimTK::Force 
subclass Force::MobilityLinearStop and is logically part of Force.h. The file
assumes that Force.h will have included all necessary declarations. **/

namespace SimTK {

/** Model a compliant stop element that acts to keep a mobilizer coordinate q
within specified bounds. This force element generates no force when the 
coordinate is in bounds, but when either the lower or upper bound is exceeded
it generates a generalized force opposing further violation of the bound. The
generated force is composed of a stiffness force (or torque, or generalized
force) that is linear in the violation of the bound, and dissipation that is 
linear in the rate qdot.

@bug MobilityLinearStop currently works only for coordinates q whose time 
derivatives qdot are just the corresponding generalized speed u. That is the
case for translational mobilities, pin, universal, and gimbal mobilizers but
not free or ball mobilizers.

<h3>Theory:</h3>
Given a mobilizer coordinate q and limits q_low and q_high, the generalized
force generated here is:
<pre>
      {           0,             q_low <= q <= q_high
  f = { min(0, -k*x*(1+d*qdot)), q > q_high, x=q-q_high (x>0, f<=0)
      { max(0, -k*x*(1-d*qdot)), q < q_low,  x=q-q_low  (x<0, f>=0)
</pre>
where k is a stiffness parameter, and d is a dissipation coefficient.

Note that dissipation occurs both during compression and expansion, but we
will never generate a "sticking" force. This is a Hunt and Crossley-like
dissipation model. It has the nice property that the damping force is zero
when you first touch the stop.
**/
class SimTK_SIMBODY_EXPORT Force::MobilityLinearStop : public Force {
public:
    /** Create a %MobilityLinearStop force element on a particular generalized 
    coordinate.
    
    @param forces   subsystem to which this force element should be added
    @param mobod    mobilizer to which the force should be applied
    @param whichQ   to which of the mobilizer's generalized coordinates q
                      should this force be applied (first is 0)?
    @param defaultStiffness     default stop stiffness (>= 0)
    @param defaultDissipation   default stop dissipation coefficient (>= 0)
    @param defaultQLow          default lower bound (-Infinity)
    @param defaultQHigh         default upper bound (+Infinity)

    Note that if you have an integer value for the generalized coordinate (q)
    index, you have to cast it to a MobilizerQIndex here. The generalized 
    coordinates are numbered starting with 0 for each mobilizer. Here is an 
    example: @code
        GeneralForceSubsystem forces;
        MobilizedBody::Pin pinJoint(...);
        MobilityLinearStop myStop(forces, pinJoint, MobilizerQIndex(0),
                                  k, d, -Pi/4, Pi/4);
    @endcode
    **/
    MobilityLinearStop(GeneralForceSubsystem&    forces, 
                       const MobilizedBody&      mobod, 
                       MobilizerQIndex           whichQ, 
                       Real                      defaultStiffness,
                       Real                      defaultDissipation,
                       Real                      defaultQLow =-Infinity,
                       Real                      defaultQHigh= Infinity);
    
    /** Default constructor creates an empty handle that can be assigned to
    refer to any %MobilityLinearStop object. **/
    MobilityLinearStop() {}

    /** Provide new values for the default lower and upper bounds of this stop. 
    This is a topological change because it affects the value that the 
    containing System's default state will have when realizeTopology() is 
    called. This is for use during construction, not for during a simulation 
    where you should be using setBounds() to set the bounds in a State rather 
    than in the System.

    @param      defaultQLow      
        Default lower bound (generalized coordinate value); <= defaultQHigh.
        Set to <code>-Infinity</code> to disable the lower bound by default.          
    @param      defaultQHigh     
        Default upper bound (generalized coordinate value); >= defaultQLow.
        Set to <code>Infinity</code> to disable the upper bound by default.          

    @returns a writable reference to this modified force element
    @see setBounds(), getDefaultLowerBound(), getDefaultUpperBound() **/
    MobilityLinearStop& setDefaultBounds(Real defaultQLow, Real defaultQHigh);


    /** Provide new values for the default material properties of this stop, 
    which are assumed to be the same for the upper and lower stops. This is
    a topological change because it affects the value that the containing 
    System's default state will have when realizeTopology() is called. This is 
    for use during construction, not for during a simulation where you should 
    be using setMaterialProperties() to set the material properties in a State 
    rather than in the System.

    @param defaultStiffness     default stop stiffness (>= 0)
    @param defaultDissipation   default stop dissipation coefficient (>= 0)
         
    @returns a writable reference to this modified force element
    @see setMaterialProperties(), getDefaultStiffness(), 
         getDefaultDissipation() **/
    MobilityLinearStop& setDefaultMaterialProperties
                            (Real defaultStiffness, Real defaultDissipation);

    /** Return the default value for the stop's lower bound (a generalized
    coordinate value). This is normally set at construction. 
    @see getDefaultUpperBound(), getLowerBound() **/
    Real getDefaultLowerBound() const;

    /** Return the default value for the stop's upper bound (a generalized
    coordinate value). This is normally set at construction. 
    @see getDefaultLowerBound(), getUpperBound() **/
    Real getDefaultUpperBound() const;

    /** Return the default value for the stop material's stiffness k. This is 
    normally set at construction. 
    @see getDefaultDissipation(), getStiffness() **/
    Real getDefaultStiffness() const;

    /** Return the default value for the stop's material's dissipation 
    coefficient d. This is normally set at construction. 
    @see getDefaultStiffness(), getDissipation() **/
    Real getDefaultDissipation() const;

    /** Change the values of the lower and upper bounds in the given \a state;
    these may differ from the default values supplied at construction.

    @param      state    
        The State in which the bounds are changed.
    @param      qLow     
        Lower bound (generalized coordinate value); <= qHigh. Set to
        \c -Infinity to disable the lower bound.             
    @param      qHigh    
        Upper bound (generalized coordinate value); >= qLow. Set to
        \c Infinity to disable the upper bound.

    Changing these bounds invalidates Stage::Dynamics and above in the \a state
    since it can affect force generation. 
    @see setDefaultBounds(), getLowerBound(), getUpperBound() **/
    void setBounds(State& state, Real qLow, Real qHigh) const;

    /** Change the values of this stop's material properties in the given 
    \a state; these may differ from the default values supplied at construction.

    @param      state    
        The State in which the material properties are changed.
    @param      stiffness     
        The stop material stiffness k (>= 0).
    @param      dissipation   
        The stop material dissipation coefficient d (>= 0).

    Changing these properties invalidates Stage::Dynamics and above in the 
    \a state since it can affect force generation. 
    @see setDefaultMaterialProperties(), getStiffness(), getDissipation() **/
    void setMaterialProperties
                        (State& state, Real stiffness, Real dissipation) const;

    /** Return the value for the lower bound that is stored in the 
    given \a state. Note that this is not the same thing as the default lower
    bound that was supplied on construction. 
    @see getUpperBound(), setBounds(), getDefaultLowerBound() **/
    Real getLowerBound(const State& state) const;

    /** Return the value for the upper bound that is stored in the 
    given \a state. Note that this is not the same thing as the default upper
    bound that was supplied on construction. 
    @see getLowerBound(), setBounds(), getDefaultUpperBound() **/
    Real getUpperBound(const State& state) const;

    /** Return the value for the stop material's stiffness k that is stored in 
    the given \a state. Note that this is not the same thing as the default 
    stiffness that was supplied on construction. 
    @see getDissipation(), setMaterialProperties(), getDefaultStiffness() **/
    Real getStiffness(const State& state) const;

    /** Return the value for the stop material's dissipation coefficient d that
    is stored in the given \a state. Note that this is not the same thing as the
    default dissipation coefficient that was supplied on construction. 
    @see getStiffness(), setMaterialProperties(), getDefaultDissipation() **/
    Real getDissipation(const State& state) const;

    /** @cond **/
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(MobilityLinearStop, 
                                             MobilityLinearStopImpl, Force);
    /** @endcond **/
};



} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_STOP_H_
