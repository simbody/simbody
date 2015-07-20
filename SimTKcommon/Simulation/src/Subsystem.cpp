/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2006-15 Stanford University and the Authors.        *
 * Authors: Michael Sherman                                                   *
 * Contributors: Peter Eastman                                                *
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


/**@file
 *
 * Implementation of Subsystem and Subsystem::Guts.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/Measure.h"
#include "SimTKcommon/internal/State.h"
#include "SimTKcommon/internal/EventHandler.h"
#include "SimTKcommon/internal/EventReporter.h"
#include "SimTKcommon/internal/System.h"
#include "SimTKcommon/internal/Subsystem.h"
#include "SimTKcommon/internal/Study.h"

#include "SimTKcommon/internal/MeasureImplementation.h"

#include <cassert>

using namespace SimTK;

namespace {
//------------------------ INITIALIZE MEASURES ACTION --------------------------
// Each Subsystem adds one of these Actions to the System's Initialization
// Event that occurs at the start of a simulation so that its Measures will get
// initialized. 
class InitializeMeasuresAction: public EventAction {
public:
    InitializeMeasuresAction(const Subsystem::Guts& subsys) 
    :   EventAction(Change), m_subsys(subsys) {}

private:
    InitializeMeasuresAction* cloneVirtual() const override 
    {   return new InitializeMeasuresAction(*this); }

    void changeVirtual
       (Study&                  study,
        const Event&            /*ignored*/,
        const EventTriggers&    /*ignored*/,
        EventChangeResult&      result) const override 
    {
        m_subsys.initializeMeasures(study.updInternalState());
        result.reportExitStatus(EventChangeResult::Succeeded);
    }

private:
    const Subsystem::Guts& m_subsys;
};

} // anonymous file-scope namespace


//==============================================================================
//                                 SUBSYSTEM
//==============================================================================

Subsystem::Subsystem(const Subsystem& src) : guts(0) {
    if (src.guts) {
        guts = src.guts->clone();
        guts->setOwnerSubsystemHandle(*this);
    }

}

Subsystem& Subsystem::operator=(const Subsystem& src) {
    if (!isSameSubsystem(src)) {
        if (isOwnerHandle())
            delete guts;
        guts=0;
        if (src.guts) {
            guts = src.guts->clone();
            guts->setOwnerSubsystemHandle(*this);
        }
    }
    return *this;
}

Subsystem::~Subsystem() {
    if (guts && isOwnerHandle())
        delete guts;
    guts=0;
}

void Subsystem::adoptSubsystemGuts(Subsystem::Guts* g) {
    SimTK_ASSERT_ALWAYS(g, "Subsystem::adoptSubsystemGuts(): can't adopt null Guts");
    SimTK_ASSERT_ALWAYS(!guts,
        "Subsystem::adoptSubsystemGuts(): this Subsystem handle is already in use");
    guts = g;
    guts->setOwnerSubsystemHandle(*this);
}


//==============================================================================
//                           SUBSYSTEM :: GUTS
//==============================================================================

// This serves as default constructor.
Subsystem::Guts::Guts(const String& name, const String& version)
:   m_subsystemName(name), m_subsystemVersion(version),
    m_mySystem(nullptr), m_mySubsystemIndex(InvalidSubsystemIndex), 
    m_myHandle(nullptr), m_subsystemTopologyRealized(false)
{ 
}

// Copy constructor isn't very useful. Note that it doesn't copy Measures.
Subsystem::Guts::Guts(const Subsystem::Guts& src) 
:   m_subsystemName(src.m_subsystemName), 
    m_subsystemVersion(src.m_subsystemVersion),
    m_mySystem(nullptr), m_mySubsystemIndex(InvalidSubsystemIndex), 
    m_myHandle(nullptr), m_subsystemTopologyRealized(false)
{
}

// Destructor must unreference and possibly delete measures.
Subsystem::Guts::~Guts() {
    for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
        if (m_measures[mx]->decrRefCount()==0) delete m_measures[mx];
    m_myHandle = 0;
    invalidateSubsystemTopologyCache();
}

// Subsystem is being added to a System. Can now set up any System-level
// resources, such as Event actions. 
// TODO: concrete subsystem needs to get a chance to do this too.
void Subsystem::Guts::setSystem(System& sys, SubsystemIndex index) {
    SimTK_ASSERT(!isInSystem(), "Subsystem::setSystem()");
    SimTK_ASSERT(index.isValid(), "Subsystem::setSystem()");
    m_mySystem = &sys;
    m_mySubsystemIndex = index;

    // Now that we're part of a System we can help ourselves to some
    // System goodies.
    acquireSystemResources();

}

MeasureIndex Subsystem::Guts::adoptMeasure(AbstractMeasure& m) {
    SimTK_ASSERT(m.hasImpl(), "Subsystem::Guts::adoptMeasure()");

    // In Debug mode check that this measure hasn't already been adopted.
    // This is an expensive check if there are lots of measures.
    SimTK_ASSERT(std::find(m_measures.begin(), m_measures.end(), &m.getImpl())
                 == m_measures.end(), "Subsystem::Guts::adoptMeasure()");

    invalidateSubsystemTopologyCache();
    const MeasureIndex mx(m_measures.size());
    m_measures.push_back(&m.updImpl());
    m_measures.back()->incrRefCount();
    m_measures.back()->setSubsystem(updOwnerSubsystemHandle(), mx);
    return mx;
}

void Subsystem::Guts::initializeMeasures(State& state) const {
    for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
        m_measures[mx]->initialize(state);
}

bool Subsystem::Guts::isInSameSystem(const Subsystem& otherSubsystem) const {
    return isInSystem() && otherSubsystem.isInSystem()
        && getSystem().isSameSystem(otherSubsystem.getSystem());
}

// Invalidating a Subsystem's topology cache forces invalidation of the
// whole System's topology cache, which will in turn invalidate all the other
// Subsystem's topology caches.
void Subsystem::Guts::invalidateSubsystemTopologyCache() const {
    if (m_subsystemTopologyRealized) {
        auto mThis = const_cast<Subsystem::Guts*>(this);
        mThis->m_subsystemTopologyRealized = false;
        if (isInSystem()) 
            getSystem().invalidateSystemTopologyCache();
    }
}

    // wrappers for Subsystem::Guts virtuals

//------------------------------------------------------------------------------
//                                  CLONE
//------------------------------------------------------------------------------
Subsystem::Guts* Subsystem::Guts::clone() const {
    return cloneImpl();
}

//------------------------------------------------------------------------------
//                     REALIZE SUBSYSTEM TOPOLOGY
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemTopology(State& s) const {
    SimTK_STAGECHECK_EQ_ALWAYS(getStage(s), Stage::Empty, 
        "Subsystem::Guts::realizeSubsystemTopology()");
    realizeSubsystemTopologyImpl(s);

    // Realize this Subsystem's Measures.
    for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
        m_measures[mx]->realizeTopology(s);

    auto mThis = const_cast<Subsystem::Guts*>(this);
    mThis->m_subsystemTopologyRealized = true; // mark subsys itself (mutable)
    advanceToStage(s, Stage::Topology);  // mark the State as well
}

//------------------------------------------------------------------------------
//                        REALIZE SUBSYSTEM MODEL
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemModel(State& s) const {
    SimTK_STAGECHECK_TOPOLOGY_REALIZED_ALWAYS(subsystemTopologyHasBeenRealized(),
        "Subsystem", getName(), "Subsystem::Guts::realizeSubsystemModel()");

    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage::Topology, 
        "Subsystem::Guts::realizeSubsystemModel()");
    if (getStage(s) < Stage::Model) {
        realizeSubsystemModelImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizeModel(s);

        advanceToStage(s, Stage::Model);
    }
}

//------------------------------------------------------------------------------
//                     REALIZE SUBSYSTEM INSTANCE
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemInstance(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Instance).prev(), 
        "Subsystem::Guts::realizeSubsystemInstance()");
    if (getStage(s) < Stage::Instance) {
        realizeSubsystemInstanceImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizeInstance(s);

        advanceToStage(s, Stage::Instance);
    }
}

//------------------------------------------------------------------------------
//                         REALIZE SUBSYSTEM TIME
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemTime(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Time).prev(), 
        "Subsystem::Guts::realizeTime()");
    if (getStage(s) < Stage::Time) {
        realizeSubsystemTimeImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizeTime(s);

        advanceToStage(s, Stage::Time);
    }
}

//------------------------------------------------------------------------------
//                     REALIZE SUBSYSTEM POSITION
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemPosition(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Position).prev(), 
        "Subsystem::Guts::realizeSubsystemPosition()");
    if (getStage(s) < Stage::Position) {
        realizeSubsystemPositionImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizePosition(s);

        advanceToStage(s, Stage::Position);
    }
}

//------------------------------------------------------------------------------
//                       REALIZE SUBSYSTEM VELOCITY
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemVelocity(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Velocity).prev(), 
        "Subsystem::Guts::realizeSubsystemVelocity()");
    if (getStage(s) < Stage::Velocity) {
        realizeSubsystemVelocityImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizeVelocity(s);

        advanceToStage(s, Stage::Velocity);
    }
}

//------------------------------------------------------------------------------
//                       REALIZE SUBSYSTEM DYNAMICS
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemDynamics(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Dynamics).prev(), 
        "Subsystem::Guts::realizeSubsystemDynamics()");
    if (getStage(s) < Stage::Dynamics) {
        realizeSubsystemDynamicsImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizeDynamics(s);

        advanceToStage(s, Stage::Dynamics);
    }
}

//------------------------------------------------------------------------------
//                     REALIZE SUBSYSTEM ACCELERATION
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemAcceleration(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Acceleration).prev(), 
        "Subsystem::Guts::realizeSubsystemAcceleration()");
    if (getStage(s) < Stage::Acceleration) {
        realizeSubsystemAccelerationImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizeAcceleration(s);

        advanceToStage(s, Stage::Acceleration);
    }
}

//------------------------------------------------------------------------------
//                         REALIZE SUBSYSTEM REPORT
//------------------------------------------------------------------------------
void Subsystem::Guts::realizeSubsystemReport(const State& s) const { 
    SimTK_STAGECHECK_GE_ALWAYS(getStage(s), Stage(Stage::Report).prev(), 
        "Subsystem::Guts::realizeSubsystemReport()");
    if (getStage(s) < Stage::Report) {
        realizeSubsystemReportImpl(s);

        // Realize this Subsystem's Measures.
        for (MeasureIndex mx(0); mx < m_measures.size(); ++mx)
            m_measures[mx]->realizeReport(s);

        advanceToStage(s, Stage::Report);
    }
}

//------------------------------------------------------------------------------
//                  CALC DECORATIVE GEOMETRY AND APPEND
//------------------------------------------------------------------------------
void Subsystem::Guts::calcDecorativeGeometryAndAppend
   (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const 
{
    calcDecorativeGeometryAndAppendImpl(s,stage,geom);
}

//------------------------------------------------------------------------------
//                      ACQUIRE SYSTEM RESOURCES
//------------------------------------------------------------------------------
void Subsystem::Guts::acquireSystemResources() 
{
    SimTK_ASSERT_ALWAYS(isInSystem(), 
        "Subsystem::Guts::acquireSystemResources(): "
        "This Subsystem has not been put into a System.");

    // Every Subsystem has as its first Initialization action to call its
    // Measures' initialize() methods.
    m_mySystem->updInitializationEvent()
        .adoptEventAction(new InitializeMeasuresAction(*this));

    acquireSystemResourcesImpl();
}



