#ifndef SimTK_SIMBODY_MATTER_SUBTREE_H_
#define SimTK_SIMBODY_MATTER_SUBTREE_H_

/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTK Simbody(tm)                         *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2007-9 Stanford University and the Authors.         *
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

#include <cassert>
#include <iosfwd>

namespace SimTK {

class SimbodyMatterSubsystem;
class MobilizedBody;
class SimbodyMatterSubtree;
class SimbodyMatterSubtreeResults;

/** A SimbodyMatterSubtree is a view of a connected subgraph of the tree of 
mobilized bodies in a SimbodyMatterSubsystem. It is used to perform kinematic 
operations on the subgraph to facilitate the handling of constraints, which 
typically involve only small subgraphs.

A SimbodyMatterSubtree is characterized by a single ancestor body A and a set
of terminal mobilized bodies T={Ti}, where A is the outmost body that is on 
the inboard path of each Ti. Note that a SimbodyMatterSubtree's "terminal" 
bodies do not have to be terminal in the full tree. The SimbodyMatterSubtree 
includes T and all "branch" mobilized bodies B={Bij} found on any path from a 
Ti to A, and A itself which serves as Ground. A may be one of the terminal 
bodies. A's mobilizer is *not* part of the SimbodyMatterSubtree. The path 
from Ti to A is called the ith branch of the SimbodyMatterSubtree; branches 
can overlap.
@verbatim
                           . .
      .    .                .
       .  .                 .
        T0      T1          T2     }
         *       *       .  *      }
     B0   *     * B1       *       }
           *   *          *        }   A SimbodyMatterSubtree with
        .    *           *  B2     }   three branches.
          . . *        *           }
                *     *            }
           B0,B1  *  *             }
                    A              }
                    .
                    .
                   ...
                 Ground
@endverbatim
Each body in the SimbodyMatterSubtree is assigned an index called a 
SubtreeBodyIndex, with the Ancestor being SubtreeBodyIndex 0 and other ids 
assigned such that ids increase going outwards along a branch. Maps are kept
in the SimbodyMatterSubtree object to track its relationship to the full tree. 

A SimbodyMatterSubtree can be constructed at Topology stage and needed ones 
can thus be precalculated and stored in the SimbodyMatterSubsystem Topology
Cache (i.e., in the System not the State). Calculations done on the 
SimbodyMatterSubtree, on the other hand, require further state information 
and cannot be stored as part of the System. For those, we define a companion 
class below called SimbodyMatterSubtreeResults.

A SimbodyMatterSubtreeResults object is initialized at Model stage, at which 
point we can determine the mobilities u and generalized coordinates q. These 
are assigned SubtreeUIndex's and SubtreeQIndex's in the same order that the 
SimbodyMatterSubtree bodies are numbered. Maps are kept in the 
SimbodyMatterSubtreeResults object to track the relationship between the 
SimbodyMatterSubtree mobilities and those in the full tree.

Note that SimbodyMatterSubtree operations are elaborate \e operators, not 
\e responses. That means the results are not stored in the State, but rather 
in the private SimbodyMatterSubtreeResults objects.

Operators here perform kinematic operations based on perturbations of the
global System State values. The supported perturbations are: 
  General
  1a same as global System state (except answers are in A rather than G)
  1b all mobility variables set
  2  all mobility variables from 1a or 1b, except for one which is perturbed 
     (q,u,udot)

  Linear
  3  all mobility variables are zero (u,udot)
  4  all mobility variables are zero *again*, except for one which is 1 
     (u,udot)
Steps 1 and 2 are designed to work together, as are 3 and 4: first evaluate
nominal kinematics; then perturb. **/
class SimTK_SIMBODY_EXPORT SimbodyMatterSubtree {
public:
    SimbodyMatterSubtree();
    SimbodyMatterSubtree(const SimbodyMatterSubtree&);
    SimbodyMatterSubtree& operator=(const SimbodyMatterSubtree&);
    ~SimbodyMatterSubtree();

    explicit SimbodyMatterSubtree(const SimbodyMatterSubsystem&);
    SimbodyMatterSubtree(const SimbodyMatterSubsystem&, 
            const Array_<MobilizedBodyIndex>& terminalBodies);

    void setSimbodyMatterSubsystem(const SimbodyMatterSubsystem& matter);
    const SimbodyMatterSubsystem& getSimbodyMatterSubsystem() const;

    // This doesn't change the associated SimbodyMatterSubsystem if there
    // is one, but does remove all the bodies from the SimbodyMatterSubtree.
    void clear();

    SimbodyMatterSubtree& addTerminalBody(MobilizedBodyIndex);

    void realizeTopology();

    int getNumSubtreeBodies() const; // includes ancestor
    MobilizedBodyIndex getAncestorMobilizedBodyIndex() const;

    // These are in the same order they were added; body[i] is the terminus
    // of branch i.
    const Array_<MobilizedBodyIndex>& getTerminalBodies() const;

    // These are indexed by SubtreeBodyIndex starting with 0 for the ancestor
    // body and monotonically increasing outwards along a branch.
    const Array_<MobilizedBodyIndex>& getAllBodies() const;

    // 0 returns an invalid Index
    SubtreeBodyIndex getParentSubtreeBodyIndex(SubtreeBodyIndex) const;
    const Array_<SubtreeBodyIndex>& 
        getChildSubtreeBodyIndices(SubtreeBodyIndex) const;

        // MODEL STAGE

    // State must be realized to at least Stage::Model for this call to work. 
    // The supplied SimbodyMatterSubtreeResults object is allocated and properly initialized to
    // be able to hold computation results from this SimbodyMatterSubtree.
    void initializeSubtreeResults(const State&, SimbodyMatterSubtreeResults&) const;

    // This can be used as a sanity check that initializeSubtreeResults() was 
    // already called in this SimbodyMatterSubtree to produce these 
    // SimbodyMatterSubtreeResults. It is by no means exhaustive but will catch
    // egregious errors.
    bool isCompatibleSubtreeResults(const SimbodyMatterSubtreeResults&) const;

        // POSITION STAGE

    // State must be realized to at least Stage::Position for this to work. SimbodyMatterSubtreeResults
    // must have already been initialized to work with this SimbodyMatterSubtree. SimbodyMatterSubtreeResults stage
    // will be Stage::Position after this call. All body transforms will be the same as
    // the corresponding ones in the state, except they will be measured from the ancestor
    // frame instead of ground. SimbodyMatterSubtree q's will be identical to corresponding State q's.
    void copyPositionsFromState(const State&, SimbodyMatterSubtreeResults&) const;

    // State must be realized to Stage::Instance. subQ must be the right length for this
    // SimbodyMatterSubtree, and SimbodyMatterSubtreeResults must have been properly initialized. SimbodyMatterSubtreeResults
    // stage will be Stage::Position after this call.
    void calcPositionsFromSubtreeQ(const State&, const Vector& subQ, SimbodyMatterSubtreeResults&) const;

    // Calculates a perturbed position result starting with the subQ's and position results
    // which must already be in SimbodyMatterSubtreeResults.
    void perturbPositions(const State&, SubtreeQIndex subQIndex, Real perturbation, SimbodyMatterSubtreeResults&) const;


        // VELOCITY STAGE

    // State must be realized to at least Stage::Velocity for this to work. SimbodyMatterSubtreeResults
    // must already be at Stage::Position. SimbodyMatterSubtreeResults stage
    // will be Stage::Velocity after this call. All subtree body spatial velocities will be
    // the same as in the State, except measured relative to A and expressed in A. SimbodyMatterSubtree u's
    // will be identical to corresponding State u's.
    void copyVelocitiesFromState(const State&, SimbodyMatterSubtreeResults&) const;

    // State must be realized to Stage::Instance. subU must be the right length for this
    // SimbodyMatterSubtree, and SimbodyMatterSubtreeResults must already be at Stage::Position. SimbodyMatterSubtreeResults
    // stage will be Stage::Velocity after this call.
    void calcVelocitiesFromSubtreeU(const State&, const Vector& subU, SimbodyMatterSubtreeResults&) const;

    // State must be realized to Stage::Instance and SimbodyMatterSubtreeResults must already be at
    // Stage::Position. SimbodyMatterSubtreeResults stage will be Stage::Velocity after this call, but
    // all SimbodyMatterSubtree u's and body velocities will be zero.
    void calcVelocitiesFromZeroU(const State&, SimbodyMatterSubtreeResults&) const;

    // Calculates a perturbed velocity result starting with the subU's and velocity results
    // which must already be in SimbodyMatterSubtreeResults.
    void perturbVelocities(const State&, SubtreeUIndex subUIndex, Real perturbation, SimbodyMatterSubtreeResults&) const;


        // ACCELERATION STAGE

    // State must be realized to at least Stage::Acceleration for this to work. SimbodyMatterSubtreeResults
    // must already be at Stage::Velocity. SimbodyMatterSubtreeResults stage
    // will be Stage::Acceleration after this call. All subtree body spatial accelerations will be
    // the same as in the State, except measured relative to A and expressed in A. SimbodyMatterSubtree udots
    // will be identical to corresponding State udots.
    void copyAccelerationsFromState(const State&, SimbodyMatterSubtreeResults&) const;

    // State must be realized to Stage::Instance. subUDot must be the right length for this
    // SimbodyMatterSubtree, and SimbodyMatterSubtreeResults must already be at Stage::Velocity. SimbodyMatterSubtreeResults
    // stage will be Stage::Acceleration after this call.
    void calcAccelerationsFromSubtreeUDot(const State&, const Vector& subUDot, SimbodyMatterSubtreeResults&) const;

    // State must be realized to Stage::Instance and SimbodyMatterSubtreeResults must already be at
    // Stage::Velocity. SimbodyMatterSubtreeResults stage will be Stage::Acceleration after this call.
    // All SimbodyMatterSubtree udots's will be zero, body accelerations will have only their bias values
    // (coriolis accelerations from nonzero u's).
    void calcAccelerationsFromZeroUDot(const State&, SimbodyMatterSubtreeResults&) const;

    // Calculates a perturbed velocity result starting with the subUDot's and acceleration results
    // which must already be in SimbodyMatterSubtreeResults.
    void perturbAccelerations(const State&, SubtreeUIndex subUDotIndex, Real perturbation, SimbodyMatterSubtreeResults&) const;

    class SubtreeRep;
private:
    SubtreeRep* rep;
    const SubtreeRep& getRep() const {assert(rep);return *rep;}
    SubtreeRep&       updRep()       {assert(rep);return *rep;}
};

SimTK_SIMBODY_EXPORT std::ostream& 
operator<<(std::ostream&, const SimbodyMatterSubtree&);

/*
 * This is the writable "cache" for a SimbodyMatterSubtree. Once the full State has
 * been realized to the Model stage, a SimbodyMatterSubtree can initialize one of these
 * objects and then use it to hold operator results.
 */
class SimTK_SIMBODY_EXPORT SimbodyMatterSubtreeResults {
public:
    SimbodyMatterSubtreeResults();
    SimbodyMatterSubtreeResults(const SimbodyMatterSubtreeResults&);
    SimbodyMatterSubtreeResults& operator=(const SimbodyMatterSubtreeResults&);
    ~SimbodyMatterSubtreeResults();

    void clear();

    void reallocateBodies(int nBodies);
    void addMobilities(SubtreeBodyIndex, QIndex qStart, int nq, UIndex uStart, int nu);
    void realizeModel(const Vector& stateQ, const Vector& stateU);

    Stage getStage() const;

    int getNumSubtreeBodies() const;
    int getNumSubtreeQs() const;
    int getNumSubtreeUs() const;

    const Vector&     getSubtreeQ() const;
    const Transform&  getSubtreeBodyTransform(SubtreeBodyIndex) const; // from ancestor frame

    const Vector&     getSubtreeU() const;
    const SpatialVec& getSubtreeBodyVelocity(SubtreeBodyIndex) const; // measured & expressed  in ancestor frame

    const Vector&     getSubtreeUDot() const;
    const SpatialVec& getSubtreeBodyAcceleration(SubtreeBodyIndex) const; // measured & expressed in ancestor frame

    // These are indexed by SubtreeQIndex and SubtreeUIndex.
    const Array_<QIndex>& getQSubset() const; // subset of Subsystem Qs used by this SimbodyMatterSubtree
    const Array_<UIndex>& getUSubset() const; // subset of Subsystem Us used by this SimbodyMatterSubtree

    void findSubtreeBodyQ(SubtreeBodyIndex, SubtreeQIndex& qStart, int& nq) const; // indices into QSubset
    void findSubtreeBodyU(SubtreeBodyIndex, SubtreeUIndex& uStart, int& nu) const; // indices into USubset

    class SubtreeResultsRep;
private:
    friend class SimbodyMatterSubtree;
    SubtreeResultsRep* rep;
    const SubtreeResultsRep& getRep() const {assert(rep);return *rep;}
    SubtreeResultsRep&       updRep()       {assert(rep);return *rep;}
};

SimTK_SIMBODY_EXPORT std::ostream& 
operator<<(std::ostream&, const SimbodyMatterSubtreeResults&);

} // namespace SimTK

#endif // SimTK_SIMBODY_MATTER_SUBTREE_H_
