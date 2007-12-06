#ifndef SimTK_SIMBODY_MATTER_SUBSYSTEM_SUBTREE_H_
#define SimTK_SIMBODY_MATTER_SUBSYSTEM_SUBTREE_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007 Stanford University and the Authors.           *
 * Authors: Michael Sherman                                                   *
 * Contributors:                                                              *
 *                                                                            *
 * Permission is hereby granted, free of charge, to any person obtaining a    *
 * copy of this software and associated documentation files (the "Software"), *
 * to deal in the Software without restriction, including without limitation  *
 * the rights to use, copy, modify, merge, publish, distribute, sublicense,   *
 * and/or sell copies of the Software, and to permit persons to whom the      *
 * Software is furnished to do so, subject to the following conditions:       *
 *                                                                            *
 * The above copyright notice and this permission notice shall be included in *
 * all copies or substantial portions of the Software.                        *
 *                                                                            *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR *
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,   *
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL    *
 * THE AUTHORS, CONTRIBUTORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,    *
 * DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR      *
 * OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE  *
 * USE OR OTHER DEALINGS IN THE SOFTWARE.                                     *
 * -------------------------------------------------------------------------- */

#include "SimTKcommon.h"
#include "simbody/internal/common.h"
#include "simbody/internal/SimbodyMatterSubsystem.h"

#include <cassert>
#include <vector>
#include <iosfwd>

namespace SimTK {

/**
 * A Subtree is a view of a connected subgraph of the tree of mobilized bodies
 * in a SimbodyMatterSubsystem. It is used to perform kinematic operations on
 * the subgraph to facilitate the handling of constraints, which typically
 * involve only small subgraphs.
 *
 * A Subtree is characterized by a single ancestor body A and a set of terminal
 * mobilized bodies T={Ti}, where A is the outmost body which is on the inboard
 * path of each Ti. Note that a Subtree's "terminal" bodies do not have to be terminal
 * in the full tree. The Subtree includes T and all "branch" mobilized bodies
 * B={Bij} found on any path from a Ti to A, and A itself which serves as Ground.
 * A will not be the same as any terminal body unless one of the terminal bodies
 * is Ground. A's mobilizer is *not* part of the Subtree. The path from Ti to A
 * is called the ith branch of the Subtree; branches can overlap.
 *                                                          @verbatim
 *                           . .
 *      .    .                .
 *       .  .                 .
 *        T0      T1          T2     }
 *         *       *          *      }
 *     B0   *     * B1       *       }
 *           *   *          *        }   A Subtree with
 *             *           *  B2     }   three branches.
 *              *        *           }
 *                *     *            }
 *           B0,B1  *  *             }
 *                    A              }
 *                    .
 *                    .
 *                   ...
 *                  Ground                                  @endverbatim
 *
 * Each body in the Subtree is assigned an index called a SubtreeBodyId,
 * with the Ancestor being SubtreeBodyId 0 and other ids assigned such
 * that ids increase going outwards along a branch. Maps are
 * kept in the Subtree object to track its relationship to the full tree.
 *
 * A Subtree can be constructed at Topology stage and needed ones can
 * thus be precalculated and stored in the SimbodyMatterSubsystem Topology
 * Cache (i.e., in the System not the State). Calculations done on the Subtree,
 * on the other hand, require further state information and cannot be stored
 * as part of the System. For those, we define a companion class below called
 * SubtreeResults.
 *
 * A SubtreeResults object is initialized at Model stage, at which point
 * we can determine the mobilities u and generalized coordinates q. These are
 * assigned SubtreeUId's and SubtreeQId's in the same order that the Subtree
 * bodies are numbered. Maps are kept in the SubtreeResults object to track
 * the relationship between the Subtree mobilities and those in the full tree.
 *
 * Note that Subtree operations are elaborate *operators*, not *responses*.
 * That means the results are not stored in the State, but rather in the
 * private SubtreeResults objects.
 *
 * Operators here perform kinematic operations based on perturbations of
 * the global System State values. The supported perturbations are: 
 *   General
 *   1a same as global System state (except answers are in A rather than G)
 *   1b all mobility variables set
 *   2 all mobility variables from 1a or 1b, except for one which is perturbed (q,u,udot)
 *
 *   Linear
 *   3 all mobility variables are zero (u,udot)
 *   4 all mobility variables are zero *again*, except for one which is 1 (u,udot)
 * Steps 1 and 2 are designed to work together, as are 3 and 4: first evaluate
 * nominal kinematics; then perturb.
 */

class SimTK_SIMBODY_EXPORT SimbodyMatterSubsystem::Subtree {
public:
    Subtree();
    Subtree(const Subtree&);
    Subtree& operator=(const Subtree&);
    ~Subtree();

    explicit Subtree(const SimbodyMatterSubsystem&);
    Subtree(const SimbodyMatterSubsystem&, 
            const std::vector<MobilizedBodyId>& terminalBodies);

    void setSimbodyMatterSubsystem(const SimbodyMatterSubsystem& matter);
    const SimbodyMatterSubsystem& getSimbodyMatterSubsystem() const;

    // This doesn't change the associated SimbodyMatterSubsystem if there
    // is one, but does remove all the bodies from the Subtree.
    void clear();

    Subtree& addTerminalBody(MobilizedBodyId);

    void realizeTopology();

    int getNumSubtreeBodies() const; // includes ancestor
    MobilizedBodyId getAncestorBody() const;

    // These are in the same order they were added; body[i] is the terminus
    // of branch i.
    const std::vector<MobilizedBodyId>& getTerminalBodies() const;

    // These are indexed by SubtreeBodyId starting with 0 for the ancestor body
    // and monotonically increasing outwards along a branch.
    const std::vector<MobilizedBodyId>& getAllBodies() const;

    SubtreeBodyId getParentSubtreeBodyId(SubtreeBodyId) const; // 0 returns an invalid Id
    const std::vector<SubtreeBodyId>& getChildSubtreeBodyIds(SubtreeBodyId) const;

        // MODEL STAGE

    // State must be realized to at least Stage::Model for this call to work. 
    // The supplied SubtreeResults object is allocated and properly initialized to
    // be able to hold computation results from this Subtree.
    void initializeSubtreeResults(const State&, SubtreeResults&) const;

    // This can be used as a sanity check that initializeSubtreeResults() was already called
    // in this Subtree to produce these SubtreeResults. It is by no means exhaustive but
    // will catch egregious errors.
    bool isCompatibleSubtreeResults(const SubtreeResults&) const;

        // POSITION STAGE

    // State must be realized to at least Stage::Position for this to work. SubtreeResults
    // must have already been initialized to work with this Subtree. SubtreeResults stage
    // will be Stage::Position after this call. All body transforms will be the same as
    // the corresponding ones in the state, except they will be measured from the ancestor
    // frame instead of ground. Subtree q's will be identical to corresponding State q's.
    void copyPositionsFromState(const State&, SubtreeResults&) const;

    // State must be realized to Stage::Instance. subQ must be the right length for this
    // Subtree, and SubtreeResults must have been properly initialized. SubtreeResults
    // stage will be Stage::Position after this call.
    void calcPositionsFromSubtreeQ(const State&, const Vector& subQ, SubtreeResults&) const;

    // Calculates a perturbed position result starting with the subQ's and position results
    // which must already be in SubtreeResults.
    void perturbPositions(const State&, SubtreeQId subQIndex, Real perturbation, SubtreeResults&) const;


        // VELOCITY STAGE

    // State must be realized to at least Stage::Velocity for this to work. SubtreeResults
    // must already be at Stage::Position. SubtreeResults stage
    // will be Stage::Velocity after this call. All subtree body spatial velocities will be
    // the same as in the State, except measured relative to A and expressed in A. Subtree u's
    // will be identical to corresponding State u's.
    void copyVelocitiesFromState(const State&, SubtreeResults&) const;

    // State must be realized to Stage::Instance. subU must be the right length for this
    // Subtree, and SubtreeResults must already be at Stage::Position. SubtreeResults
    // stage will be Stage::Velocity after this call.
    void calcVelocitiesFromSubtreeU(const State&, const Vector& subU, SubtreeResults&) const;

    // State must be realized to Stage::Instance and SubtreeResults must already be at
    // Stage::Position. SubtreeResults stage will be Stage::Velocity after this call, but
    // all Subtree u's and body velocities will be zero.
    void calcVelocitiesFromZeroU(const State&, SubtreeResults&) const;

    // Calculates a perturbed velocity result starting with the subU's and velocity results
    // which must already be in SubtreeResults.
    void perturbVelocities(const State&, SubtreeUId subUIndex, Real perturbation, SubtreeResults&) const;


        // ACCELERATION STAGE

    // State must be realized to at least Stage::Acceleration for this to work. SubtreeResults
    // must already be at Stage::Velocity. SubtreeResults stage
    // will be Stage::Acceleration after this call. All subtree body spatial accelerations will be
    // the same as in the State, except measured relative to A and expressed in A. Subtree udots
    // will be identical to corresponding State udots.
    void copyAccelerationsFromState(const State&, SubtreeResults&) const;

    // State must be realized to Stage::Instance. subUDot must be the right length for this
    // Subtree, and SubtreeResults must already be at Stage::Velocity. SubtreeResults
    // stage will be Stage::Acceleration after this call.
    void calcAccelerationsFromSubtreeUDot(const State&, const Vector& subUDot, SubtreeResults&) const;

    // State must be realized to Stage::Instance and SubtreeResults must already be at
    // Stage::Velocity. SubtreeResults stage will be Stage::Acceleration after this call.
    // All Subtree udots's will be zero, body accelerations will have only their bias values
    // (coriolis accelerations from nonzero u's).
    void calcAccelerationsFromZeroUDot(const State&, SubtreeResults&) const;

    // Calculates a perturbed velocity result starting with the subUDot's and acceleration results
    // which must already be in SubtreeResults.
    void perturbAccelerations(const State&, SubtreeUId subUDotIndex, Real perturbation, SubtreeResults&) const;

    class SubtreeRep;
private:
    SubtreeRep* rep;
    const SubtreeRep& getRep() const {assert(rep);return *rep;}
    SubtreeRep&       updRep()       {assert(rep);return *rep;}
};

SimTK_SIMBODY_EXPORT std::ostream& 
operator<<(std::ostream&, const SimbodyMatterSubsystem::Subtree&);

/*
 * This is the writable "cache" for a Subtree. Once the full State has
 * been realized to the Model stage, a Subtree can initialize one of these
 * objects and then use it to hold operator results.
 */
class SimTK_SIMBODY_EXPORT SimbodyMatterSubsystem::SubtreeResults {
public:
    SubtreeResults();
    SubtreeResults(const SubtreeResults&);
    SubtreeResults& operator=(const SubtreeResults&);
    ~SubtreeResults();

    void clear();

    void reallocateBodies(int nBodies);
    void addMobilities(SubtreeBodyId, QId qStart, int nq, UId uStart, int nu);
    void realizeModel(const Vector& stateQ, const Vector& stateU);

    Stage getStage() const;

    int getNumSubtreeBodies() const;
    int getNumSubtreeQs() const;
    int getNumSubtreeUs() const;

    const Vector&     getSubtreeQ() const;
    const Transform&  getSubtreeBodyTransform(SubtreeBodyId) const; // from ancestor frame

    const Vector&     getSubtreeU() const;
    const SpatialVec& getSubtreeBodyVelocity(SubtreeBodyId) const; // measured & expressed  in ancestor frame

    const Vector&     getSubtreeUDot() const;
    const SpatialVec& getSubtreeBodyAcceleration(SubtreeBodyId) const; // measured & expressed in ancestor frame

    // These are indexed by SubtreeQId and SubtreeUId.
    const std::vector<QId>& getQSubset() const; // subset of Subsystem Qs used by this Subtree
    const std::vector<UId>& getUSubset() const; // subset of Subsystem Us used by this Subtree

    void findSubtreeBodyQ(SubtreeBodyId, SubtreeQId& qStart, int& nq) const; // indices into QSubset
    void findSubtreeBodyU(SubtreeBodyId, SubtreeUId& uStart, int& nu) const; // indices into USubset

    class SubtreeResultsRep;
private:
    friend class Subtree;
    SubtreeResultsRep* rep;
    const SubtreeResultsRep& getRep() const {assert(rep);return *rep;}
    SubtreeResultsRep&       updRep()       {assert(rep);return *rep;}
};

SimTK_SIMBODY_EXPORT std::ostream& 
operator<<(std::ostream&, const SimbodyMatterSubsystem::SubtreeResults&);

} // namespace SimTK

#endif // SimTK_SIMBODY_MATTER_SUBSYSTEM_SUBTREE_H_
