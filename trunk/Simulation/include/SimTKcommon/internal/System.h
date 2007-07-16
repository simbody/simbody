#ifndef SimTK_SYSTEM_H_
#define SimTK_SYSTEM_H_

/* Portions copyright (c) 2006-7 Stanford University and Michael Sherman.
 * Contributors:
 * 
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including 
 * without limitation the rights to use, copy, modify, merge, publish, 
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 * 
 * The above copyright notice and this permission notice shall be included 
 * in all copies or substantial portions of the Software.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"
#include "SimTKcommon/internal/State.h"

namespace SimTK {

class State;
class Subsystem;
class DecorativeGeometry;

/**
 * The abstract parent of all Systems.
 * A System serves as a mediator for a group of interacting Subsystems.
 * All will share a single system State, and typically subsystems will
 * need access to content in the state which is produced by other
 * subsystems. 
 *
 * A System provides a unique subsystem index (a small positive integer)
 * for each of its subsystems, and the subsystems are constructed
 * knowing their indices. The indices are used subsequently by the subsystems
 * to find their own entries in the system state, and by each subsystem
 * to refer to others within the same system. Index 0 is reserved for 
 * use by the System itself, e.g. for system-global state variables.
 *
 * Concrete Systems understand the kinds of subsystems they contain. 
 * For example, a MultibodySystem might contain a mechanical subsystem,
 * a force subsystem, and a geometry subsystem. At each computation
 * stage, a subsystem is realized in a single operation. That operation
 * can refer to computations from already-realized subsystems, but
 * cannot initiate computation in other subsystems. The System must
 * know the proper order with which to realize the subsystems at each
 * stage, and that ordering is likely to vary with stage. For example,
 * at Position stage the mechanical positions must be realized
 * before the configuration-dependent force elements. However, at
 * Acceleration stage, the force elements must be realized before the
 * mechanical accelerations can be calculated.
 */
class SimTK_SimTKCOMMON_EXPORT System {
public:
    System() : rep(0) { }
    ~System();
    System(const System&);
    System& operator=(const System&);

    const String& getName()    const;
    const String& getVersion() const;

    /// The following call must be made after any topological change
    /// has been made to this System, before the System can be used
    /// to perform any computations. Perhaps surprisingly, the method
    /// is const. That's because the topology cannot be changed by this method.
    /// Various mutable "cache" entries will get calculated, including the
    /// default State, a reference to which is returned.
    /// The returned State has already been realized through the highest
    /// Stage, using the defaults for the Model-stage variables, and
    /// default values for all later stage variables as well. You can access
    /// this same default State again using getDefaultState().
    /// If the current topology has already been realized, this call does 
    /// nothing but return a refernece to the already-built default State.
    const State& realizeTopology() const;

    /// You can check whether realizeTopology() has been called since the last
    /// topological change to this Syatem. If you don't check and just plunge
    /// ahead you are likely to encounter an exception since very few things
    /// will work without topology having been realized.
    bool topologyHasBeenRealized() const;

    /// This is available after realizeTopology(), and will throw an
    /// exception if realizeTopology() has not been called since the
    /// most recent topological change to this System. This method returns the
    /// same reference returned by realizeTopology(). The State to which
    /// a reference is returned was created by the most recent
    /// realizeTopology() call. It has default values for all the 
    /// Model-stage variables, and has already been realized through
    /// the highest Stage, so you can use it directly to obtain information
    /// about the System in its default state or you can use this state
    /// to initialize other States to which you have write access. Those
    /// States are then suitable for further computation with this System.
    const State& getDefaultState() const;

    /// This call is required if Model-stage variables are changed from
    /// their default values. The System topology must already have been
    /// realized (that is, realizeTopology() must have been called since
    /// the last topological change made to the System). Also, the supplied
    /// State must already have been initialized to work with this System
    /// either by copying the default state or some other State of this System.
    /// If it has already been realized to Stage::Model or higher, nothing
    /// happens here. Otherwise, all the state variables at Stage::Instance or
    /// higher are allocated or reallocated (if necessary), and reinitialized
    /// to their default values. NOTE: any State information at Stage::Instance
    /// or higher in the passed-in State is *destroyed* here. The number, types
    /// and memory locations of those state variables will change, so any
    /// existing references or pointers to them are invalid after this call.
    /// Note that this routine modifies its argument, but makes no changes
    /// at all to the System itself.
    void realizeModel(State&) const;

    /// Realize the entire System to the indicated Stage. The passed-in
    /// State must have been initialized to work with this System, and
    /// it must already have been realized through Stage::Model, since
    /// the realize() method doesn't have write access to the State.
    /// If the state has already been realized to the requested stage
    /// or higher, nothing happens here. Otherwise, the state is realized
    /// one stage at a time until it reaches the requested stage. 
    void realize(const State& s, Stage g = Stage::HighestValid) const;

    // Generate all decorative geometry computable at a specific stage. This will
    // throw an exception if the state hasn't already been realized
    // to that stage. Note that the list is not inclusive -- you have to
    // request geometry from each stage to get all of it.
    // This routine asks each subsystem in succession to generate its decorative geometry
    // and append it to the end of the Array.
    // If the stage is Stage::Topology, realizeTopology() must already have
    // been called but the State is ignored.
    void calcDecorativeGeometryAndAppend(const State&, Stage, Array<DecorativeGeometry>&) const;

    /// This operator can be called at Stage::Instance or higher and 
    /// returns a rough estimate of a length of time we consider significant
    /// for this system. For example, this could be the period of the highest-frequency
    /// oscillation that we care about. This can be used as a hint by 
    /// numerical integrators in choosing their initial step size, and 
    /// suggests how velocity variables should be scaled relative to their
    /// corresponding position variables.
    Real calcTimescale(const State&) const;

    /// This operator can be called at Stage::Position to calculate a weighting
    /// vector w, with one entry for each state variable y={q,u,z}, ordered
    /// the same as Y in the State and calculated specifically for the current
    /// values of Y in the State. Weight wi is proportional to the "importance"
    /// of state variable yi with respect to some criteria determined by
    /// the System, such that wi*dyi=1 indicates that a change dyi in state
    /// yi produces approximately a unit change in the weighting criteria. This
    /// is intended for use by numerical integration methods for step size control.
    /// The idea is to allow creation of a weighted RMS norm which returns 1 just when
    /// all the state variable changes have a unit effect. The norm is RMS(W*dy) where
    /// W=diag(w). A value of 1 for this norm would typically be a huge error. 
    /// For example, if your accuracy requirement is 0.1%, you would test that
    /// the weighted RMS norm is <= .001.
    /// We expect this operation to be fairly expensive and thus the integrator
    /// is expected to invoke it only occasionally.
    void calcYUnitWeights(const State&, Vector& weights) const;

    /// This provides scaling information for each of the position and velocity
    /// constraints (YErr) in the State. The tolerance is the absolute error in the
    /// constraint which is considered a "unit violation" of that state. Then
    /// if T=diag(tol) and c the vector of constraint errors, we can use a
    /// weighted RMS norm condition like RMS(T*c) <= accuracy to define
    /// when constraints have been adequately met.
    /// This is expected to be a cheap operation but not to change during
    /// a study. State must be realized to Stage::Model.
    void calcYErrUnitTolerances(const State&, Vector& tolerances) const;

    static Real calcWeightedRMSNorm(const Vector& values, const Vector& weights) {
        assert(weights.size() == values.size());
        if (values.size()==0) return 0;
        Real sumsq = 0;
        for (int i=0; i<values.size(); ++i) {
            const Real wv = weights[i]*values[i];
            sumsq += wv*wv;
        }
        return std::sqrt(sumsq/weights.size());
    }

    static Real calcWeightedInfinityNorm(const Vector& values, const Vector& weights) {
        assert(weights.size() == values.size());
        if (values.size()==0) return 0;
        Real maxval = 0;
        for (int i=0; i<values.size(); ++i) {
            const Real wv = std::abs(weights[i]*values[i]);
            if (wv > maxval) maxval=wv;
        }
        return maxval;
    }

    /// This operator can be called at Stage::Position to take a vector
    /// of absolute state variable error estimates and return a weighted
    /// norm. This method is intended for use by numerical integration methods
    /// for step size control. This is a weighted norm, calculated so
    /// that a return value of 1 would indicate a "unit" error, which would 
    /// be huge. If your accuracy requirement is 0.1%, you would test that
    /// the norm return here is <= .001.
    /// XXX OBSOLETE? TODO
    Real calcYErrorNorm(const State&, const Vector& y_err) const;

    /// Take over ownership of the supplied subsystem and install it into 
    /// the next free subsystem slot. The new slot index is returned.
    SubsystemId adoptSubsystem(Subsystem& child);

    /// How may Subsystems are in here?
    int getNSubsystems() const;
    /// Obtain read-only access to a particular subsystem by its index.
    const Subsystem& getSubsystem(SubsystemId)   const;
    /// Obtain writable access to a particular subsystem by its index.
    Subsystem&       updSubsystem(SubsystemId);

    // Internal use only
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // There can be multiple handles on the same System.
    bool isSameSystem(const System& otherSystem) const;

    explicit System(class SystemRep* r) : rep(r) { }
    bool          hasRep() const {return rep!=0;}
    const SystemRep& getRep() const {assert(rep); return *rep;}
    SystemRep&       updRep() const {assert(rep); return *rep;}
protected:
    class SystemRep* rep;
};


/// The abstract parent of all Studies.
/// TODO
class SimTK_SimTKCOMMON_EXPORT Study {
public:
    Study() : rep(0) { }
    ~Study();
    Study(const Study&);
    Study& operator=(const Study&);

    Study(const System& sys);

    const System& getSystem() const;
    const State&  getState()  const;
    State&        updState();

    /// Is this handle the owner of this rep? This is true if the
    /// handle is empty or if its rep points back here.
    bool isOwnerHandle() const;
    bool isEmptyHandle() const;

    // Internal use only
    explicit Study(class StudyRep* r) : rep(r) { }
    bool            hasRep() const {return rep!=0;}
    const StudyRep& getRep() const {assert(rep); return *rep;}
protected:
    class StudyRep* rep;
};

} // namespace SimTK

#endif // SimTK_SYSTEM_H_
