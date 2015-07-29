#ifndef SimTK_SimTKCOMMON_SIMBODY_ELEMENT_IMPL_H_
#define SimTK_SimTKCOMMON_SIMBODY_ELEMENT_IMPL_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
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

/** @file
Declares the abstract SimbodyElementImpl class, the base class for all 
element implementations. **/

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/ReferencePtr.h"

#include <ostream>
#include <cassert>
#include <algorithm>
#include <memory>
#include <cstdint>

namespace SimTK {
class Subsystem;
class State;
class DecorativeGeometry;

//==============================================================================
//                          SIMBODY ELEMENT IMPL
//==============================================================================
/** TODO: just a sketch for now.
Abstract base class for lowest-level components of a Simbody System. Elements
are always owned by a particular Subsystem.
**/
class SimbodyElementImpl {
public:
    virtual ~SimbodyElementImpl() {}

    /** Duplicate this element and return a pointer to the new copy. The new
    element does not belong to any Subsystem even if this one does. **/
    SimbodyElementImpl* clone() const 
    {   return cloneVirtual(); }

    // acquireSubsystemResources() is only virtual; it is invoked through
    // nonvirtual setSubsystem() below.

    void realizeTopology(State& state) const 
    {   realizeTopologyVirtual(state); }
    void realizeModel(State& state) const 
    {   realizeModelVirtual(state); }
    void realizeInstance(const State& state) const 
    {   realizeInstanceVirtual(state); }
    void realizeTime(const State& state) const 
    {   realizeTimeVirtual(state); }
    void realizePosition(const State& state) const 
    {   realizePositionVirtual(state); }
    void realizeVelocity(const State& state) const 
    {   realizeVelocityVirtual(state); }
    void realizeDynamics(const State& state) const 
    {   realizeDynamicsVirtual(state); }
    void realizeAcceleration(const State& state) const 
    {   realizeAccelerationVirtual(state); }
    void realizeReport(const State& state) const 
    {   realizeReportVirtual(state); }

    void calcDecorativeGeometryAndAppend
       (const State& state, Stage stage, Array_<DecorativeGeometry>& geom) const 
    {   calcDecorativeGeometryAndAppendVirtual(state, stage, geom); }

    /** Transfer the internal state from a corresponding element of the same
    type to the internal state of this element. Typically the "from" element
    is in a different System, which is fine as long as the "from" state is 
    compatible with that System. **/
    void transferState(const SimbodyElementImpl& fromElement,
                       const State&              fromState,
                       State&                    toState) const
    {   transferStateVirtual(fromElement, fromState, toState); }


    /** Set this element's Subsystem and index within the Subsystem, immediately
    after this element has been adopted by `subsystem`. It is an error if
    this element already belongs to a Subsystem. The meaning of the index may 
    vary for different types of elements, but it must be nonnegative and should 
    uniquely identify the element within the Subsystem. After these are set, the
    element will be given a chance to acquire Subsystem resources via a call to 
    the virtual method acquireSubsystemResourcesVirtual(). **/
    void setSubsystem(Subsystem& subsystem, int indexInSubsystem) {
        assert(!isInSubsystem()); assert(indexInSubsystem >= 0);
        m_subsystem         = &subsystem;
        m_indexInSubsystem  = indexInSubsystem;
        acquireSubsystemResourcesVirtual();
    }

    /** Returns `true` if this element belongs to a Subsystem. **/
    bool isInSubsystem() const {return m_subsystem;}

    /** Returns `true` if this element, and the supplied `other` element each 
    belong to a Subsystem, and they both belong to the *same* subsystem. **/
    bool isInSameSubsystem(const SimbodyElementImpl& other);

    /** Return a const reference to the containing Subsystem. Throws an 
    exception if the element doesn't belong to a Subsystem; use isInSubsystem()
    if you aren't sure. **/
    const Subsystem& getSubsystem() const
    {   assert(isInSubsystem()); return *m_subsystem;}

    /** Return a writable reference to the containing Subsystem. Throws an 
    exception if the element doesn't belong to a Subsystem; use isInSubsystem()
    if you aren't sure. **/
    Subsystem& updSubsystem()
    {   assert(isInSubsystem()); return *m_subsystem;}

    /** Return the element index within its containing Subsystem. Throws an 
    exception if the element doesn't belong to a Subsystem; use isInSubsystem()
    if you aren't sure. **/
    int getIndexInSubsystem() const 
    {   assert(isInSubsystem()); return m_indexInSubsystem;}
protected:
    /** Default constructor initializes subsystem to null and index to -1. **/
    SimbodyElementImpl() : m_indexInSubsystem(-1) {}

    /** Base class copy constructor does *not* preserve Subsystem and index;
    the source is ignored and the result is the same as default 
    construction. **/
    SimbodyElementImpl(const SimbodyElementImpl&) : SimbodyElementImpl() {}

    /** Base class move constructor does *not* preserve Subsystem and index;
    the source is ignored and the result is the same as default 
    construction. **/
    SimbodyElementImpl(SimbodyElementImpl&&) : SimbodyElementImpl() {}

    /** Copy assignment is suppressed. **/
    SimbodyElementImpl& operator=(const SimbodyElementImpl&) = delete;

    /** Move assignment is suppressed. **/
    SimbodyElementImpl& operator=(SimbodyElementImpl&&) = delete;

    /** Invoke the copy constructor for your concrete element implementation to
    create a new instance of this element implementation on the heap. Be careful
    not to copy any data members that refer to other objects because this is
    most likely being copied to make an entirely new System, meaning that all
    the other objects will have different addresses. If you use ReferencePtr
    smart pointers for your external pointers those will be set to `nullptr`
    automatically when copied. **/
    virtual SimbodyElementImpl* cloneVirtual() const = 0;

    /** This method is called just once when the element object is first
    adopted by its Subsystem. The owning Subsystem will have been set in the
    element so this method can talk to it. The default implementation assumes
    the element requires no resources from the Subsystem. **/
    virtual void acquireSubsystemResourcesVirtual() {}

    virtual void realizeTopologyVirtual    (State& state) const {}
    virtual void realizeModelVirtual       (State& state) const {}
    virtual void realizeInstanceVirtual    (const State& state) const {}
    virtual void realizeTimeVirtual        (const State& state) const {}
    virtual void realizePositionVirtual    (const State& state) const {}
    virtual void realizeVelocityVirtual    (const State& state) const {}
    virtual void realizeDynamicsVirtual    (const State& state) const {}
    virtual void realizeAccelerationVirtual(const State& state) const {}
    virtual void realizeReportVirtual      (const State& state) const {}

    virtual void calcDecorativeGeometryAndAppendVirtual
       (const State& s, Stage stage, Array_<DecorativeGeometry>& geom) const {}

    /** Default implementation assumes element uses no state resources. **/
    virtual void transferStateVirtual
       (const SimbodyElementImpl& fromElement,
        const State&              fromState,
        State&                    toState) const {}

private:
    ReferencePtr<Subsystem> m_subsystem;
    int                     m_indexInSubsystem;
};

} // namespace SimTK


#endif // SimTK_SimTKCOMMON_SIMBODY_ELEMENT_IMPL_H_
