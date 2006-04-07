#ifndef SimTK_SIMBODY_STATE_H_
#define SimTK_SIMBODY_STATE_H_

#include "simbody/internal/common.h"

#include <cassert>

namespace SimTK {

/** 
 * Coordinate allocation:
 * Body 0 is ground and has no coordinates. Welded bodies have a higher body number
 * but no coordinates. We still set qoff and uoff as we would if these had dofs,
 * but nq==nu==0 for these so the next body has the same qoff & uoff.
 *
 *
 *   Body 0     1       2    3      4      ...
 *         ---------- ----- --- ---------- ---
 *      q |          |     |   |          |
 *         ---------- ----- --- ---------- ---
 *                    ^
 *                    qoff[2], nq[2], nqmax[2]
 *
 *   Body 0     1     2    3    4      ...
 *         -------- ----- --- -------- ---
 *      u |        |     |   |        |
 *         -------- ----- --- -------- ---
 *                  ^
 *                  uoff[2], nu[2] (==ndof[2])
 *
 * qdot, qdotdot are allocated exactly as for q. Derivatives of unused coordinates
 * are set to 0.
 * udot, prescribedUdot are allocated exactly as for u. All udots are set.
 * Entries of prescribedUdot which do not correspond to prescribed joints will
 * be ignored and need not be set.
 */

class SBStation {
public:
    SBStation(int bodyNum, const Vec3& s) : body(bodyNum), station(s) { }

    int  body;
    Vec3 station;
};

class SBDirection {
public:
    SBDirection(int bodyNum, const UnitVec3& d) : body(bodyNum), direction(d) { }

    int  body;
    UnitVec3 direction;
};

// This is the handle class for the hidden SBState implementation.
class SimTK_SIMBODY_API SBState {
public:
    SBState() : rep(0) { }
    ~SBState();
    SBState(const SBState&);
    SBState& operator=(const SBState&);

    Stage getCurrentStage() const;
    void backUpToStage(Stage) const;  // cache is mutable
    void advanceToStage(Stage) const; // can only advance 1 stage

    int allocateQRange(int nq); // qdot, qdotdot also allocated in cache
    int allocateURange(int nu); // udot                    "
    int allocateZRange(int nz); // zdot                    "
    int allocateDiscreteVariable(Stage, AbstractValue* v);
    int allocateCacheEntry(Stage, AbstractValue* v);

    const Vector& getQ() const;
    const Vector& getU() const;
    const Vector& getZ() const;

    // Stage >= dv.stage
    const AbstractValue& getDiscreteVariable(int index) const;

    // Stage >= g, g=min(dv.stage-1, Modeled); then back up to g
    AbstractValue&       updDiscreteVariable(int index);

    // Stage >= ce.stage
    const AbstractValue& getCacheEntry(int index) const;

    // Stage >= ce.stage-1; does not change stage
    AbstractValue&       updCacheEntry(int index) const; // mutable

// ignore everything below here, please.
    class SBStateRep* rep;
    const SBStateRep& getRep() const {assert(rep); return *rep;}
    SBStateRep&       updRep()       {assert(rep); return *rep;}
};

} // namespace SimTK

#endif // SimTK_SIMBODY_STATE_H_
