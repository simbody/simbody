#ifndef SimTK_SIMBODY_MULTIBODY_SYSTEM_H_
#define SimTK_SIMBODY_MULTIBODY_SYSTEM_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-15 Stanford University and the Authors.        *
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
#include "simbody/internal/common.h"

#include <vector>

namespace SimTK {

class SimbodyMatterSubsystem;
class ForceSubsystem;
class DecorationSubsystem;
class GeneralContactSubsystem;
class ImpactEvent;
class ContactEvent;


/** The job of the MultibodySystem class is to coordinate the activities of 
various subsystems which can be part of a multibody system. We insist on 
having exactly one SimbodyMatterSubsystem, and expect to have one or more
ForceSubsystems.

There will also be a generic System-level Subsystem for global variables that
span multiple subsystems.
**/
class SimTK_SIMBODY_EXPORT MultibodySystem : public System {
public:
    MultibodySystem();
    explicit MultibodySystem(SimbodyMatterSubsystem& m);

    //--------------------------------------------------------------------------
    /** @name                       System Energy **/
    /**@{**/
    /** Calculate the total potential energy of the system. The state must be at
    Dynamics stage or later. **/
    const Real calcPotentialEnergy(const State&) const;

    /** Calculate the total kinetic energy of the system. The state must be at 
    Velocity stage or later. **/
    const Real calcKineticEnergy(const State&) const;

    /** Calculate the total energy of the system. The state must be at Dynamics 
    stage or later. **/
    Real calcEnergy(const State& s) const {
        return calcPotentialEnergy(s)+calcKineticEnergy(s);
    }
    /**@}**/

    //--------------------------------------------------------------------------
    /** @name                       System Forces
    These cache entries belong to the global subsystem, which zeroes them at
    the start of the corresponding stage. They are filled in by the force 
    subsystems when they are realized to each stage. Forces are cumulative 
    from stage to stage, so the Dynamics stage includes everything. That may
    then be accessed by the matter subsystem in Acceleration stage to 
    generate the accelerations. **/
    /**@{**/
    const Vector_<SpatialVec>& getRigidBodyForces(const State&, Stage) const;
    const Vector_<Vec3>&       getParticleForces (const State&, Stage) const;
    const Vector&              getMobilityForces (const State&, Stage) const;

    // These routines are for use by force subsystems during Dynamics stage.
    Vector_<SpatialVec>& updRigidBodyForces(const State&, Stage) const;
    Vector_<Vec3>&       updParticleForces (const State&, Stage) const;
    Vector&              updMobilityForces (const State&, Stage) const;
    /**@}**/

    //--------------------------------------------------------------------------
    /** @name                       System Events
    **/
    /**@{**/
    const ImpactEvent& getImpactEvent() const;
    ImpactEvent& updImpactEvent();

    const ContactEvent& getContactEvent() const;
    ContactEvent& updContactEvent();
    /**@}**/

    //--------------------------------------------------------------------------
    /** @name                       Subsystems 
    These are normally invoked from withing Subsystem constructors rather
    than by end users. **/
    /**@{**/
    /** Steals ownership of the source; returns Subsystem ID number. **/
    int addForceSubsystem(ForceSubsystem&);

    int setMatterSubsystem(SimbodyMatterSubsystem&);
    const SimbodyMatterSubsystem& getMatterSubsystem() const;
    SimbodyMatterSubsystem&       updMatterSubsystem();
    bool hasMatterSubsystem() const;

    int setDecorationSubsystem(DecorationSubsystem&);
    const DecorationSubsystem& getDecorationSubsystem() const;
    DecorationSubsystem&       updDecorationSubsystem();
    bool hasDecorationSubsystem() const;

    int setContactSubsystem(GeneralContactSubsystem&);
    const GeneralContactSubsystem& getContactSubsystem() const;
    GeneralContactSubsystem&       updContactSubsystem();
    bool hasContactSubsystem() const;
    /**@}**/


    //--------------------------------------------------------------------------
    /**@name              Advanced/obscure/debugging/obsolete
    You probably don't want to use these methods. **/
    /**@{**/
    // Private implementation.
    SimTK_PIMPL_DOWNCAST(MultibodySystem, System);
    class MultibodySystemRep& updRep();
    const MultibodySystemRep& getRep() const;
    /**@}**/

protected:
    explicit MultibodySystem(MultibodySystemRep*);

private:
    void initRep(); // shared by constructors
};



//==============================================================================
//                                IMPACT EVENT
//==============================================================================
// This is the Event that is triggered whenever any currently-inactive, 
// position-level unilateral constraint is about to be violated. We expect
// an EventAction to correct that by making a discontinuous velocity change.
// This should trigger a follow-up ContactEvent.
class ImpactEvent : public Event {
public:
    ImpactEvent() 
    :   Event("Impact") {}

private:
    ImpactEvent* cloneVirtual() const override 
    {   return new ImpactEvent(*this); }
};



//==============================================================================
//                               CONTACT EVENT
//==============================================================================
// This is the Event that is triggered whenever a change of contact active set 
// is needed. Triggers include liftoff, contact without impact, persistent 
// contact after an impact, sliding-to-stiction, stiction-to-impending slip.
// TODO: Come here for Painleve problems also?
class ContactEvent : public Event {
public:
    ContactEvent() 
    :   Event("Contact change") {}

private:
    ContactEvent* cloneVirtual() const override 
    {   return new ContactEvent(*this); }
};


} // namespace SimTK

#endif // SimTK_SIMBODY_MULTIBODY_SYSTEM_H_
