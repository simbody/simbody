#ifndef SimTK_SIMBODY_CONTACT_EVENT_ACTIONS_H_
#define SimTK_SIMBODY_CONTACT_EVENT_ACTIONS_H_

/* -------------------------------------------------------------------------- *
 *                               Simbody(tm)                                  *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2015 Stanford University and the Authors.           *
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

// Collects declarations of event actions associated with event-based 
// conditional constraint handling. Currently:
//      ImpactEventAction
//      ContactEventAction
// 

#include "simbody/internal/common.h"
#include "simbody/internal/ImpulseSolver.h"

#include <memory>
#include <utility>

namespace SimTK {
class SimbodyMatterSubsystemRep;
class ImpulseSolver;

//==============================================================================
//                           IMPACT EVENT ACTION
//==============================================================================
/* This is the EventAction we will perform when an ImpactEvent triggers.
This action assumes that the given Study's internal state has velocities that
need to be corrected using an impact model. The state is corrected by applying
impulses that correct the velocities in a manner that satisfies the equations
of motion and the impact model equations.

The action here is:
(1) Considering only currently-active constraints, or constraints that are 
    inactive but proximal,
(2) determine an acceptable impulse that corrects hard-constraint velocity 
    errors while respecting restitution and friction laws,
(3) follow up with a ContactEventAction to choose the active set and project
    active constraints onto the constraint manifolds.

Note that on entry we will accept as proximal unilateral contacts that are
active, *even if* they are not being maintained to position tolerance (for
example, they may be Baumgarte stabilized). When we choose the final active set, 
those position errors should remain acceptable since impact doesn't change them. 
However, the constraints cannot be considered fully active since their 
velocities will have changed arbitrarily due to the applied impulse.
*/
class ImpactEventAction : public EventAction {
public:
    explicit ImpactEventAction(const SimbodyMatterSubsystemRep& matter);

    // A unilateral contact constraint is proximal if it is currently active
    // or satisfies a proximity condition.
    static void findProximalConstraints
       (const MultibodySystem&          mbs,
        const State&                    state,
        Real                            proximityTol,
        Array_<UnilateralContactIndex>& proximalUniContacts,
        Array_<UnilateralContactIndex>& needToActivate);

private:
    void changeVirtual(Study&                  study,
                       const Event&            event,
                       const EventTriggers&    triggers,
                       EventChangeResult&      result) const override;

    ImpactEventAction* cloneVirtual() const override 
    {   return new ImpactEventAction(*this); }

    void collectConstraintInfo
       (const MultibodySystem&                  mbs,
        const State&                            state,
        const Array_<UnilateralContactIndex>&   proximalUniContacts,
        Array_<ImpulseSolver::UncondRT>&        unconditional,
        Array_<ImpulseSolver::UniContactRT>&    uniContact,
        Array_<MultiplierIndex>&                allParticipating) const; 

    void classifyUnilateralContactsForSimultaneousImpact
       (Real                                    posTol,
        Real                                    velTol,
        const Vector&                           verr,
        const Vector&                           expansionImpulse,
        Array_<ImpulseSolver::UniContactRT>&    uniContacts, 
        Array_<int>&                            impacters,
        Array_<int>&                            expanders,
        Array_<int>&                            observers,
        Array_<MultiplierIndex>&                participaters,
        Array_<MultiplierIndex>&                expanding) const;

    void calcCoefficientsOfRestitution
       (const SimbodyMatterSubsystem& matter,
        const State& s, const Vector& verr,
        bool disableRestitution,
        Real consTol, Real defCaptureVel, Real defMinCORVel,
        Array_<ImpulseSolver::UniContactRT>& uniContact) const;

    bool calcExpansionImpulseIfAny
       (const Array_<ImpulseSolver::UniContactRT>&  uniContacts,
        const Array_<int>&                          impacters,
        const Vector&                               compressionImpulse,
        Vector&                                     expansionImpulse,
        Array_<int>&                                expanders) const; 

    const SimbodyMatterSubsystemRep&            m_matter;
    ResetOnCopy<std::unique_ptr<ImpulseSolver>> m_solver;
};


//==============================================================================
//                           CONTACT EVENT ACTION
//==============================================================================
/* This is the EventAction we will perform when a ContactEvent triggers.
This action assumes that the conditional constraint active set in the given
Study's internal state needs to be revised. The current set generates a 
constraint force solution that is inconsistent with one or more of the 
conditional constraint inequality conditions. The active set
is corrected and new, consistent constraint forces are computed.

The action here is
(1) Considering only currently-active constraints, or constraints that are 
    inactive but proximal and lingering,
(2) choose an acceptable active set for contact, sticking, and sliding 
    constraint equations,
(3) project active, hard position constraints onto the position manifold, and
(4) project active, hard velocity constraints onto the velocity manifold.

Soft constraints need only be satisfied at the acceleration level.
*/
class ContactEventAction : public EventAction {
public:
    explicit ContactEventAction(const SimbodyMatterSubsystemRep& matter) 
    :   EventAction(Change), m_matter(matter) {}

    /* Given a state and tolerances, update the state's conditional constraint
    active set. The only constraints considered are those that are
      - already active in the given state, or
      - inactive position constraints that are proximal and not separating,
      - inactive velocity constraints that have the right sign, 
      - all acceleration/force constraints.
    
    Constraint projection is done with the final active set. */
    static void updateActiveSet
       (const MultibodySystem&  mbs,
        State&                  state, // in/out
        Real                    proximalTol,
        Real                    velocityTol)
    {
        Array_<UnilateralContactIndex> lingering;
        findLingeringConstraintsFromScratch(mbs,state,proximalTol,velocityTol,
                                              lingering);
        updateActiveSetFromLingering(mbs,state,lingering);
    }

    static void updateActiveSetFromProximals
       (const MultibodySystem&                mbs,
        State&                                state, // in/out
        const Array_<UnilateralContactIndex>& proximals,
        Real                                  velocityTol)
    {
        Array_<UnilateralContactIndex> lingering;
        findLingeringConstraintsFromProximals(mbs,state,proximals,velocityTol,
                                              lingering);
        updateActiveSetFromLingering(mbs,state,lingering);
    }

    static Array_<UnilateralContactIndex> findActiveSet
       (const MultibodySystem&          mbs,
        const State&                    state); 

private:
    static void updateActiveSetFromLingering
       (const MultibodySystem&                mbs,
        State&                                state, // in/out
        const Array_<UnilateralContactIndex>& lingering); 

    // A unilateral contact constraint is lingering if it is (1) currently 
    // active, or (2) proximal and not separating to given tolerances.
    static void findLingeringConstraintsFromScratch
       (const MultibodySystem&          mbs,
        const State&                    state,
        Real                            proximityTol,
        Real                            velocityTol,
        Array_<UnilateralContactIndex>& lingeringUniContacts); 

    // Here we look only at contact constraints that have already been deemed
    // proximal (ignoring tolerance). A proximal constraints is lingering
    // if it is (1) currently active, or (2) not separating to given tolerance.
    static void findLingeringConstraintsFromProximals
       (const MultibodySystem&                mbs,
        const State&                          state,
        const Array_<UnilateralContactIndex>& proximals,
        Real                                  velocityTol,
        Array_<UnilateralContactIndex>&       lingeringUniContacts); 

    static void activateActiveSet
       (const MultibodySystem&                              mbs,
        State&                                              state,
        const Array_<UnilateralContactIndex>&               fullSet,
        const Array_<std::pair<UnilateralContactIndex,
                               CondConstraint::Condition>>& activeSubset);

    static bool scoreActiveSet
       (const MultibodySystem&                              mbs,
        State&                                              state,
        const Array_<UnilateralContactIndex>&               lingering,
        const Array_<std::pair<UnilateralContactIndex,
                               CondConstraint::Condition>>& activeSubset,
        Real&                                               norm2Lambda,
        std::pair<UnilateralContactIndex,Real>&             worstForce,
        std::pair<UnilateralContactIndex,Real>&             worstAcc); 

    static void solveActive
       (const Array_<int, MultiplierIndex>&  mult2participating, // m->mp
		const Matrix&                  Wp,            // Gp M\ ~Gp
        const Vector&                  Dp,
        const Vector&                  rhsp,          // e.g. aerrp
		const Array_<ConstraintIndex>       unconditional,
		const Array_<std::pair<UnilateralContactIndex,
                     CondConstraint::Condition>>&  activeSubset,
        Vector&                        lambdap
        );

    void changeVirtual(Study&                  study,
                       const Event&            event,
                       const EventTriggers&    triggers,
                       EventChangeResult&      result) const override;

    ContactEventAction* cloneVirtual() const override 
    {   return new ContactEventAction(*this); }

    const SimbodyMatterSubsystemRep&    m_matter;
};

} // namespace


#endif // SimTK_SIMBODY_CONTACT_EVENT_ACTIONS_H_


