#ifndef SimTK_SIMBODY_MOBILIZED_BODY_GROUND_H_
#define SimTK_SIMBODY_MOBILIZED_BODY_GROUND_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-13 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Paul Mitiguy, Peter Eastman                                  *
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
Declares the MobilizedBody::Ground class. **/

#include "simbody/internal/MobilizedBody.h"

namespace SimTK {

/** This is a special type of "mobilized" body generated automatically by
Simbody as a placeholder for Ground in the 0th slot for a 
SimbodyMatterSubsystem's mobilized bodies; don't create this yourself. The body 
type will also be Ground. You can think of this as a Weld-like mobilizer that 
connects the Ground body to the world, located at the Ground origin. The 
reaction force in this mobilizer represents the total torque and force 
applied by the System to Ground. This mobilizer is not available for
users -- if you want to weld something to Ground use MobilizedBody::Weld
instead. 
@see MobilizedBody::Weld, SimbodyMatterSubsystem::getGround(),
     SimbodyMatterSubsystem::updGround() **/
class SimTK_SIMBODY_EXPORT MobilizedBody::Ground : public MobilizedBody {
public:
    Ground();

    /** Add some artwork to Ground where the Visualizer can find it. **/
    Ground& addBodyDecoration(const Transform& X_GD, 
                              const DecorativeGeometry& artwork) 
    {
        (void)MobilizedBody::addBodyDecoration(X_GD,artwork); return *this;
    }

    /** @cond **/ // hide from Doxygen
    SimTK_INSERT_DERIVED_HANDLE_DECLARATIONS(Ground, GroundImpl, MobilizedBody);
    /** @endcond **/
};

} // namespace SimTK

#endif // SimTK_SIMBODY_MOBILIZED_BODY_GROUND_H_



