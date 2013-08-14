#ifndef SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_DAMPER_H_
#define SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_DAMPER_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-13 Stanford University and the Authors.        *
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

#include "SimTKcommon.h"
#include "simbody/internal/Force.h"

/** @file
This contains the user-visible API ("handle" class) for the SimTK::Force 
subclass Force::MobilityLinearDamper and is logically part of Force.h. The file
assumes that Force.h will have included all necessary declarations. **/

namespace SimTK {

/** A linear damper that acts along or around a mobility coordinate to apply
a generalized force there. 

The damping constant c is provided, with the generated force being -c*u where 
u is the mobility's generalized speed. This is meaningful on any mobility, since 
all our generalized speeds have physical meaning. This is not a potential force 
and hence does not contribute to potential energy.
**/

class SimTK_SIMBODY_EXPORT Force::MobilityLinearDamper : public Force {
public:
    /** Create a %MobilityLinearDamper force element on a particular mobility
    (generalized speed).
    
    @param[in,out]  forces
        The subsystem to which this force should be added.
    @param[in]      mobod    
        Mobilizer to which the force should be applied.
    @param[in]      whichU   
        To which of the mobilizer's mobilities (generalized speeds) u should 
        this force be applied (first is 0)?
    @param[in]      defaultDamping     
        The default value for the damping constant c.
    **/
    MobilityLinearDamper(GeneralForceSubsystem&     forces, 
                         const MobilizedBody&       mobod, 
                         MobilizerUIndex            whichU, 
                         Real                       defaultDamping);
    
    /** Default constructor creates an empty handle that can be assigned to
    refer to any %MobilityLinearDamper object. **/
    MobilityLinearDamper() {}

    /** Provide a new value for the default damping constant c of this damper.     
    This is a topological change because it affects the value that the 
    containing System's default state will have when realizeTopology() is 
    called. This is for use during construction, not for during a simulation 
    where you should be using setDamping() to set the damping constant in a 
    State rather than in the System.
    @param[in]      defaultDamping     
        The default value for the damping constant c.        
    @return
        A writable reference to this modified force element for convenience in
        chaining set methods.
    @see getDefaultDamping(), setDamping() **/
    MobilityLinearDamper& setDefaultDamping(Real defaultDamping);

    /** Return the default value for the damper's damping constant c. This is 
    normally set at construction but can be modified with setDefaultDamping(). 
    @see setDefaultDamping(), getDamping() **/ 
    Real getDefaultDamping() const;

    /** Change the value of the damping constant c in the given \a state; this
    may differ from the default value supplied at construction.
    @param[in,out]  state    
        The State in which the damping constant is to be changed.
    @param[in]      damping     
        The new damping constant c (>= 0) that overrides the default.
    @return 
        A const reference to this %MobilityLinearDamper element for convenience
        in chaining set methods together.

    Changing the damping constant invalidates Stage::Dynamics and above in the 
    \a state since it can affect force generation. 
    @see setDefaultDamping(), getDamping() **/
    const MobilityLinearDamper& setDamping(State&     state, 
                                           Real       damping) const;

    /** Return the value for the damping constant c that is stored in 
    the given \a state. Note that this is not the same thing as the default 
    damping constant that was supplied on construction or in 
    setDefaultDamping(). 
    @see setDamping(), getDefaultDamping() **/
    Real getDamping(const State& state) const;

    /** @name                      Deprecated
    Methods here are for backwards compatibility but have been replaced with
    better ones that you should use. **/
    /**@{**/
    /** Deprecated: Alternate signature for backwards compatibilty -- for 
    safety you should prefer using the other constructor signature that
    takes a MobilizerUIndex rather than a plain int. **/
    MobilityLinearDamper(GeneralForceSubsystem&     forces, 
                         const MobilizedBody&       mobod, 
                         int                        whichU, 
                         Real                       defaultDamping)
    {   // Invoke the other constructor.
        new(this) MobilityLinearDamper(forces, mobod, MobilizerUIndex(whichU),
                                       defaultDamping);
    }
    /**@}**/

    /** @cond **/    
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(MobilityLinearDamper, 
                                             MobilityLinearDamperImpl, Force);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_FORCE_MOBILITY_LINEAR_DAMPER_H_
