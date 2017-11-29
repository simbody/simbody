#ifndef SimTK_SIMMATH_MULTIBODY_GRAPH_MAKER_H_
#define SimTK_SIMMATH_MULTIBODY_GRAPH_MAKER_H_

/* -------------------------------------------------------------------------- *
 *                     Simbody(tm): Multibody Graph Maker                     *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013-4 Stanford University and the Authors.         *
 * Authors: Michael Sherman                                                   *
 * Contributors: Kevin He                                                     *
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
Declares the SimTK::MultibodyGraphMaker class for use in constructing a
spanning-tree-plus-constraints representation of a multibody system from a
list of its bodies and joints. **/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"

#include <utility>
#include <string>
#include <vector>
#include <map>
#include <iosfwd>

namespace SimTK {

//==============================================================================
//                          MULTIBODY GRAPH MAKER
//==============================================================================
/** Construct a reasonably good spanning-tree-plus-constraints structure for
modeling a given set of bodies and joints with a generalized coordinate 
multibody system like Simbody. Each body has a unique name; each joint 
connects two distinct bodies with one designated as the "parent" body and
the other the "child". We will output a spanning tree with Ground as the root,
containing a mobilizer for every body. Each mobilizer corresponds to one of
the joints from the input and has an "inboard" (topologically closer to Ground)
and "outboard" (further from Ground) body, chosen from the parent and child
bodies of the joint but possibly reordered in which case the mobilizer is 
marked as "reversed". Additional "free" mobilizers are added as needed between 
bodies and Ground so that there is a path from every body in the inboard 
direction all the way to Ground.

<h3>Results</h3>
The output is 
  -# A sequence of mobilizers (tree joints) in ascending order from
Ground to a set of terminal bodies. Every input body will be mobilized by one
of these, and there may also be some additional mobilized "slave" bodies.
  -# A set of constraints, representing either welds to attach slave bodies
to their masters, or loop joint constraints that correspond to particular input
joints.
  -# A correspondence between input bodies and mobilized bodies, with some
mobilized bodies designated as slaves to particular master input bodies. That
is, more than one mobilized body may correspond to the same input body.
  -# A correspondence between mobilizers and input joints, with some extra
mobilizers having been added to connect base bodies to Ground.
  -# Statistics and diagnostics.

Then to build a multibody model
  -# Run through the mobilizers in order, adding
one mobilized body to the tree each time, setting the "reversed" flag if 
appropriate (that just affects the interpretation of the generalized 
coordinates). Mass properties and joint information are obtained from the 
original bodies and joints; they are not stored here.
  -# Run through the list of constraints adding master-slave welds or loop 
joint constraints as indicated.

<h3>Inputs</h3>
  - bodies: name, mass, mustBeBaseBody flag
  - joints: name, type, parent body, child body, mustBeLoopJoint flag
  - joint type: name, # dofs, haveGoodLoopJoint flag
  - names for the Ground body and important joint types weld and free

The first body you add is assumed to be Ground and its name is used for 
recognizing Ground connections later in joints. The default names for weld
and free joints are, not surprisingly, "weld" and "free" but you can change
them.

<h3>Loop joints</h3>
Normally every joint produces a corresponding mobilizer. A joint that would
form a loop is marked as such, and then its child body is split to form a new
"slave" body that can be mobilized by the joint. We expect that in the 
constructed multibody model, a master body and its slaves will be reconnected
via weld constraints. This provides for uniform treatment of all joints, such
as providing pin joint coordinates that can wrap. However, for some joint types
that doesn't matter because the multibody system has equally good "loop joint"
constraints. For example, Simbody's Ball mobilizer uses quaternions for 
orientation and hence has no wrapping; it can be replaced with a Ball constraint
which removes the 3-dof mobilizer and 6 weld constraint equations, leaving just
3 translational constraints and indistinguishable behavior. Similarly a loop-
closing Weld joint can just be replaced by a Weld constraint, and a loop-
closing Free joint can simply be ignored.

The algorithm normally decides which joints are the loop-breakers, however you
can specify in the input that certain joints must be loop joints.

<h3>Massless bodies</h3>
Massless bodies can be very useful in defining compound joints, by composing
pin and slider joints, for example. This is fine in an internal coordinate code
as long as no massless (or inertialess) body is a terminal (most outboard) body
in the spanning tree. (Terminal massless bodies are OK if their inboard joint
has no dofs, i.e. is a weld mobilizer.) Bodies can have a mass provided; if a
movable body has zero mass the algorithm will try to avoid ending any branch of 
the spanning tree there. If it fails, the resulting spanning tree is no good and
an exception will be thrown. If no mass is provided we assume it is 1.

<h3>Base bodies</h3>
A body in the spanning tree that is directly connected to Ground is called a
"base" body. If its mobilizer is a free joint then it can be given an arbitrary
pose with respect to Ground and everything outboard of it moves together. This
is a nice property and you can influence this choice by providing an explicit
joint to Ground or designating some bodies as base bodies. Note that although 
you will then have a set of generalized coordinates permitting you to place 
these bodies arbitrarily, there is no guarantee that such placement will satisfy
constraints. If you don't designate base bodies, the algorithm will pick
them heuristically, which often leads to a good choice. The heuristic is to
pick the body that has the most children but does not itself appear as a child.
Failing that, pick the body that has the most children. In case of tie, pick
the first body among the tied ones.

<h3>Ground body</h3>
The first body you supply in the input will be interpreted as the Ground body.
Its name will be used to recognize joints that connect other bodies to Ground.
Typical names are "ground" or "world".

<h3>Levels</h3>
Each body in the spanning tree will be assigned a level, an integer that 
measures how many edges must be traversed to get from the body to Ground.
Ground is at level 0, bodies that have a direct connection to Ground (base
bodies) are at level 1, bodies whose mobilizer connects them to base bodies
are at level 2, and so on. Multibody models need to be constructed in order
of increasing level, so that the inboard body is already in the tree when the
outboard body is added. We consider the level of a mobilizer to be
the same as the level of its outboard body. Loop joints are not mobilizers and
do not have a level.
**/
class SimTK_SIMMATH_EXPORT MultibodyGraphMaker {
public:
    // Local classes.
    class Body;
    class Joint;
    class JointType;
    class Mobilizer;
    class LoopConstraint;

    /** Construct an empty %MultibodyGraphMaker object and set the default
    names for weld and free joints to "weld" and "free". **/
    MultibodyGraphMaker();


    /** Specify relevant properties of a joint type. Type name must be unique.
    Weld and free types are predefined and their names are reserved (though you
    can change the names to be used).
    @returns Small integer joint type number used for referencing. **/
    int addJointType(const std::string& name,
                     int                numMobilities,
                     bool               haveGoodLoopJointAvailable = false,
                     void*              userRef                    = 0);


    /** Add a new body (link) to the set of input bodies. 
    @param[in]      name     
        A unique string identifying this body. There are no other restrictions
        on the contents of \a name. The first body you add is considered to be
        Ground, and its name is remembered and recognized when used in joints.
    @param[in]      mass     
        The mass here is used as a graph-building hint. If the body is massless,
        be sure to set mass=0 so that the algorithm can avoid making this a 
        terminal body. The algorithm might also use \a mass to preferentially
        choose heavier bodies as terminal to improve conditioning.
    @param[in]      mustBeBaseBody
        If you feel strongly that this body should be able to move freely with
        respect to Ground, set this flag so that the algorithm will connect it
        to Ground by a free joint before attempting to build the rest of the
        tree. Alternatively, provide a joint that connects this body directly
        to Ground, in which case you should \e not set this flag. 
    @param[in]      userRef
        This is a generic user reference pointer that is kept with the body
        and can be used by the caller to map back to his or her own data 
        structure containing body information.
    @see deleteBody() **/
    void addBody(const std::string&  name, 
                 double              mass, 
                 bool                mustBeBaseBody,
                 void*               userRef = 0);

    /** Delete a body (link) from the set of input bodies. All the joints that
    reference this body will be deleted too.
    @param[in]      name     
        A unique string identifying this body. There are no other restrictions
        on the contents of \a name. Don't delete the Ground body.
    @returns \c true if the body is successfully deleted, \c false if it
        didn't exist. **/
    bool deleteBody(const std::string&  name);

    /** Add a new joint to the set of input joints.
    @param[in]      name
        A string uniquely identifying this joint. There are no other 
        restrictions on the contents of \a name.
    @param[in]      type
        A string designating the type of this joint, such as "revolute" or
        "ball". This must be chosen from the set of joint types previously
        specified.
    @param[in]      parentBodyName
        This must be the name of a body that was already specified in an earlier
        addBody() call, or it must be the designated name for the Ground body. 
        If possible, this will be used as the inboard body for the corresponding
        mobilizer.
    @param[in]      childBodyName
        This must be the name of a body that was already specified in an earlier
        addBody() call, or it must be the designated name for the Ground body.
        It must be distinct from \a parentBodyName. If possible, this will be
        used as the outboard body for the corresponding mobilizer.
    @param[in]      mustBeLoopJoint
        If you feel strongly that this joint should be chosen as a loop joint,
        set this flag. In that case the joint will not appear in the list of
        joints that are candidates for mobilizers (tree joints). Only after the
        tree has been successfully built will this joint be added, either using
        a loop joint equivalent if one is available, or by splitting the child
        body into master and slave otherwise. In the latter case this joint 
        \e will be made into a mobilizer but a loop weld joint will be added
        to attach the child slave body to its master.  
    @param[in]      userRef
        This is a generic user reference pointer that is kept with the joint
        and can be used by the caller to map back to his or her own data 
        structure containing joint information.
    @see deleteJoint() **/
    void addJoint(const std::string& name,
                  const std::string& type,
                  const std::string& parentBodyName,
                  const std::string& childBodyName,
                  bool               mustBeLoopJoint,
                  void*              userRef = 0);

    /** Delete an existing joint from the set of input joints. The bodies(links)
    referenced by the joint are expected to exist and their references to this
    joint will be removed as well.
    @param[in]      name
        A string uniquely identifying this joint. There are no other 
        restrictions on the contents of \a name. 
    @returns \c true if the joint is successfully deleted, \c false if it
        didn't exist. **/
    bool deleteJoint(const std::string& name);

    /** Generate a new multibody graph from the input data. Throws an std::exception
    if it fails, with a message in the what() string. **/
    void generateGraph();
    /** Throw away the multibody graph produced by generateGraph(). **/
    void clearGraph();

    /** Output a text representation of the multibody graph for debugging. **/
    void dumpGraph(std::ostream& out) const;

    /** Returns the number of mobilizers (tree joints) in the spanning tree.
    These are numbered in order of level. The 0th mobilizer has level 0 and
    is just a placeholder for Ground's immobile connection to the universe.
    After that come the base mobilizers at level 1, then mobilizers that 
    connect children to base bodies at level 2, and so on. This is also the 
    number of mobilized bodies, including Ground and slave bodies. **/
    int getNumMobilizers() const {return (int)mobilizers.size();}
    /** Get a Mobilizer object by its mobilizer number, ordered outwards by
    topological distance from Ground. **/
    const Mobilizer& getMobilizer(int mobilizerNum) const
    {   return mobilizers[mobilizerNum]; }

    /** Return the number of loop joint constraints that were used to close 
    loops in the graph topology. These do not include loops that were broken
    by cutting a body to make a slave body, just those where the joint itself 
    was implemented using a constraint rather than a mobilizer plus a slave. 
    The latter occurs only if were told there is a perfectly good loop joint
    constraint available; typically that applies for ball joints and not much
    else. **/
    int getNumLoopConstraints() const {return (int)constraints.size();}
    /** Get a loop constraint by its assigned number. These are assigned in
    an arbitrary order. **/
    const LoopConstraint& getLoopConstraint(int loopConstraintNum) const
    {   return constraints[loopConstraintNum]; }

    /** Return the number of bodies, including all input bodies, a ground body,
    and any slave bodies. **/
    int getNumBodies() const {return (int)bodies.size();}
    /** Get a Body object by its assigned number. These are assigned first to
    input bodies, then we add one for Ground, then we add slave bodies created
    by body splitting after that. **/
    const Body& getBody(int bodyNum) const {return bodies[bodyNum];}
    /** Return the body number assigned to the input body with the given name.
    Returns -1 if the body name is not recognized. You can't look up by name
    slave bodies that were added by the graph-making algorithm. **/
    int getBodyNum(const std::string& bodyName) const {
        std::map<std::string,int>::const_iterator p = 
            bodyName2Num.find(bodyName);
        return p==bodyName2Num.end() ? -1 : p->second;
    }

    /** Return the number of joints, including all input joints, and all joints
    added to connect otherwise disconnected bodies to Ground. Don't confuse
    these with mobilizers, which are an ordered subset of the joints that are
    chosen to form a spanning tree connecting all the bodies. **/
    int getNumJoints() const {return (int)joints.size();}
    /** Get a Joint object by its assigned number. These are assigned first to
    input joints, then we add additional joints to Ground as needed. **/
    const Joint& getJoint(int jointNum) const {return joints[jointNum];}
    /** Return the joint number assigned to the input joint with the given name.
    Returns -1 if the joint name is not recognized. You can't look up by name
    extra joints that were added by the graph-making algorithm. **/
    int getJointNum(const std::string& jointName) const {
        std::map<std::string,int>::const_iterator p = 
            jointName2Num.find(jointName);
        return p==jointName2Num.end() ? -1 : p->second;
    }

    /** Return the number of registered joint types. **/
    int getNumJointTypes() const {return (int)jointTypes.size();}
    /** Get a JointType object by its assigned number. **/
    const JointType& getJointType(int jointTypeNum) const 
    {   return jointTypes[jointTypeNum]; }
    /** Get the assigned number for a joint type from the type name. **/
    int getJointTypeNum(const std::string& jointTypeName) const {
        std::map<std::string,int>::const_iterator p = 
            jointTypeName2Num.find(jointTypeName);
        return p==jointTypeName2Num.end() ? -1 : p->second;
    }

    /** Change the name to be used to identify the weld joint type (0 dof) and 
    weld loop constraint type (6 constraints). The default is "weld".  Changing 
    this name clears and reinitializes this %MultibodyGraphMaker object. **/
    void setWeldJointTypeName(const std::string& name) 
    {   weldTypeName=name; initialize(); }
    /** Return the name currently being used to identify the weld joint type
    and weld loop constraint type. **/
    const std::string& getWeldJointTypeName() const {return weldTypeName;}

    /** Change the name to be used to identify the free (6 dof) joint type and 
    free (0 constraints) loop constraint type. The default is "free". Changing 
    this name clears and reinitializes this %MultibodyGraphMaker object. **/
    void setFreeJointTypeName(const std::string& name) 
    {   freeTypeName=name; initialize(); }
    /** Return the name currently being used to identify the free joint type
    and free loop constraint type. **/
    const std::string& getFreeJointTypeName() const {return freeTypeName;}

    /** Return the name we recognize as the Ground (or World) body. This is
    the name that was provided in the first addBody() call. **/
    const std::string& getGroundBodyName() const;
private:
    // Get writable access to bodies and joints.
    Body& updBody(int bodyNum) {return bodies[bodyNum];}
    Joint& updJoint(int jointNum) {return joints[jointNum];}
    Joint& updJoint(const std::string& name) {return joints[jointName2Num[name]];}

    void initialize();
    int splitBody(int bodyNum);
    int chooseNewBaseBody() const;
    void connectBodyToGround(int bodyNum);
    int addMobilizerForJoint(int jointNum);
    int findHeaviestUnassignedForwardJoint(int inboardBody) const;
    int findHeaviestUnassignedReverseJoint(int inboardBody) const;
    void growTree();
    void breakLoops();
    bool bodiesAreConnected(int b1, int b2) const;

    // Clear everything except for default names.
    void clear() {
        bodies.clear(); joints.clear(); jointTypes.clear();
        bodyName2Num.clear(); jointName2Num.clear(); jointTypeName2Num.clear();
        mobilizers.clear(); constraints.clear();
    }

    std::string                 weldTypeName, freeTypeName;
    std::vector<Body>           bodies; // ground + input bodies + slaves
    std::vector<Joint>          joints; // input joints + added joints
    std::vector<JointType>      jointTypes;
    std::map<std::string,int>   bodyName2Num;
    std::map<std::string,int>   jointName2Num;
    std::map<std::string,int>   jointTypeName2Num;

    // Calculated by generateGraph()
    std::vector<Mobilizer>      mobilizers; // mobilized bodies
    std::vector<LoopConstraint> constraints;
};

//------------------------------------------------------------------------------
//                      MULTIBODY GRAPH MAKER :: BODY
//------------------------------------------------------------------------------
/** Local class that collects information about bodies. **/
class MultibodyGraphMaker::Body {
public:
    explicit Body(const std::string&    name, 
                    double              mass, 
                    bool                mustBeBaseBody,
                    void*               userRef) 
    :   name(name), mass(mass), mustBeBaseBody(mustBeBaseBody), 
        userRef(userRef), level(-1), mobilizer(-1), master(-1) {}

    void forgetGraph(MultibodyGraphMaker& graph);
    int getNumFragments() const {return 1 + getNumSlaves();}
    int getNumSlaves() const {return (int)slaves.size();}
    int getNumJoints() const 
    {   return int(jointsAsChild.size() + jointsAsParent.size()); }
    bool isSlave() const {return master >= 0;}
    bool isMaster() const {return getNumSlaves()>0;}
    bool isInTree() const {return level>=0;}

    // Inputs
    std::string name;
    double      mass;
    bool        mustBeBaseBody;
    void*       userRef;

    // How this body appears in joints (input and added).
    std::vector<int>    jointsAsChild;  // where this body is the child
    std::vector<int>    jointsAsParent; // where this body is the parent

    // Disposition of this body in the spanning tree.

    int level; // Ground=0, connected to Ground=1, contact to that=2, etc.
    int mobilizer; // the unique mobilizer where this is the outboard body

    int                 master; // >=0 if this is a slave
    std::vector<int>    slaves; // slave links, if this is a master
};

//------------------------------------------------------------------------------
//                      MULTIBODY GRAPH MAKER :: JOINT
//------------------------------------------------------------------------------
/** Local class that collects information about joints. **/
class MultibodyGraphMaker::Joint {
public:
    Joint(const std::string& name, int jointTypeNum, 
          int parentBodyNum, int childBodyNum,
          bool mustBeLoopJoint, void* userRef)
    :   name(name), 
        mustBeLoopJoint(mustBeLoopJoint), 
        userRef(userRef),
        parentBodyNum(parentBodyNum), 
        childBodyNum(childBodyNum),
        jointTypeNum(jointTypeNum), 
        isAddedBaseJoint(false),
        mobilizer(-1), 
        loopConstraint(-1) {}

    /** Return true if the joint is deleted as a result of restoring it
        to the state prior to generateGraph(). **/
    bool forgetGraph(MultibodyGraphMaker& graph);

    // Only one of these will be true -- we don't consider it a LoopConstraint
    // if we split a body and weld it back.
    bool hasMobilizer() const {return mobilizer>=0;}
    bool hasLoopConstraint() const {return loopConstraint>=0;}

    // Inputs
    std::string name;
    bool        mustBeLoopJoint;
    void*       userRef;

    // Mapping of strings to indices for fast lookup.
    int parentBodyNum, childBodyNum;
    int jointTypeNum;

    bool isAddedBaseJoint; // true if this wasn't one of the input joints

    // Disposition of this joint in the multibody graph.
    int mobilizer;      // if this joint is part of the spanning tree, else -1
    int loopConstraint; // if this joint used a loop constraint, else -1
};

//------------------------------------------------------------------------------
//                   MULTIBODY GRAPH MAKER :: JOINT TYPE
//------------------------------------------------------------------------------
/** Local class that defines the properties of a known joint type. **/
class MultibodyGraphMaker::JointType {
public:
    JointType(const std::string& name, int numMobilities, 
              bool haveGoodLoopJointAvailable, void* userRef)
    :   name(name), numMobilities(numMobilities), 
        haveGoodLoopJointAvailable(haveGoodLoopJointAvailable),
        userRef(userRef) {}
    std::string name;
    int         numMobilities;
    bool        haveGoodLoopJointAvailable;
    void*       userRef;
};

//------------------------------------------------------------------------------
//                   MULTIBODY GRAPH MAKER :: MOBILIZER
//------------------------------------------------------------------------------
/** Local class that represents one of the mobilizers (tree joints) in the 
generated spanning tree. There is always a corresponding joint, although that
joint might be a ground-to-body free joint that was added automatically. **/
class MultibodyGraphMaker::Mobilizer {
public:
    Mobilizer() 
    :   joint(-1), level(-1), inboardBody(-1), outboardBody(-1),
        isReversed(false), mgm(0) {}
    Mobilizer(int jointNum, int level, int inboardBodyNum, int outboardBodyNum, 
              bool isReversed, MultibodyGraphMaker* graphMaker)
    :   joint(jointNum), level(level), inboardBody(inboardBodyNum), 
        outboardBody(outboardBodyNum), isReversed(isReversed),
        mgm(graphMaker) {}

    /** Return true if this mobilizer does not represent one of the input
    joints, but is instead a joint we added connecting a base body to ground. 
    If this returns true then there will be no user reference pointer returned
    from getJointRef(). Also, the inboard body is always ground. When you 
    create this mobilizer, the joint frames should be identity, that is, the
    joint should connect the ground frame to the outboard body frame. **/
    bool isAddedBaseMobilizer() const
    {   return mgm->getJoint(joint).isAddedBaseJoint; }
    /** Get the user reference pointer for the joint associated with this
    mobilizer, if there is such a joint. If this mobilizer doesn't correspond
    to one of the input joints then a null pointer is returned. **/
    void* getJointRef() const
    {   return mgm->getJoint(joint).userRef; } 
    /** Get the user reference pointer for the inboard body of this mobilizer. 
    The inboard body is always one of the input bodies so this will not be
    returned null unless no reference pointer was supplied in the addBody()
    call that defined this body. **/
    void* getInboardBodyRef() const
    {   return mgm->getBody(inboardBody).userRef; }   
    /** Get the user reference pointer for the outboard body of this mobilizer. 
    The outboard body may be one of the input bodies, but could also be a 
    slave body, in which case a null pointer will be returned. You can use
    getOutboardMasterBodyRef() instead to ensure that you will get a reference
    to one of the input bodies. **/
    void* getOutboardBodyRef() const
    {   return mgm->getBody(outboardBody).userRef; }
    /** Get the user reference pointer for the outboard body of this mobilizer,
    if it is one of the input bodes, or to the master body for the outboard
    body if the outboard body is a slave body. This ensures that you will get a 
    reference to one of the input bodies. **/
    void* getOutboardMasterBodyRef() const
    {   return mgm->getBody(getOutboardMasterBodyNum()).userRef; }
    /** Get the joint type name of the joint that this mobilizer represents. **/
    const std::string& getJointTypeName() const
    {   return mgm->getJointType(mgm->getJoint(joint).jointTypeNum).name; }
    /** Get the reference pointer (if any) that was provided when this 
    mobilizer's joint type was defined in an addJointType() call. **/
    void* getJointTypeRef() const
    {   return mgm->getJointType(mgm->getJoint(joint).jointTypeNum).userRef; }
    /** Return true if the outboard body of this mobilizer is a slave we 
    created in order to cut a loop, rather than one of the input bodies. **/
    bool isSlaveMobilizer() const
    {   return mgm->getBody(outboardBody).isSlave(); }
    /** Return the number of fragments into which we chopped the outboard body
    of this mobilizer. There is one fragment for the master body plus however
    many slaves of that body were created. Thus you should divide the master
    body's mass by this number to obtain the mass to be assigned to each of
    the body fragments. **/
    int getNumFragments() const 
    {   return mgm->getBody(getOutboardMasterBodyNum()).getNumFragments(); }
    /** Return true if this mobilizer represents one of the input joints but
    the sense of inboard->outboard is reversed from the parent->child sense
    defined in the input joint. In that case you should use a reverse joint
    when you build the system. **/
    bool isReversedFromJoint() const {return isReversed;}
    /** Return the level of the outboard body (Ground is level 0) **/
    int getLevel() const {return level;}

private:
friend class MultibodyGraphMaker;

    int getOutboardMasterBodyNum() const
    {   const Body& outb = mgm->getBody(outboardBody);
        return outb.isSlave() ? outb.master : outboardBody; }

    int  joint;         ///< corresponding joint (not necessarily from input)
    int  level;         ///< level of the outboard body; distance from ground
    int  inboardBody;   ///< might be ground
    int  outboardBody;  ///< might be a slave body; can't be ground
    bool isReversed;    ///< if so, inboard=child, outboard=parent

    MultibodyGraphMaker*    mgm; // just a reference to container
};


//------------------------------------------------------------------------------
//                  MULTIBODY GRAPH MAKER :: LOOP CONSTRAINT
//------------------------------------------------------------------------------
/** Local class that represents one of the constraints that were added to close 
topological loops that were cut to form the spanning tree. **/
class MultibodyGraphMaker::LoopConstraint {
public:
    LoopConstraint() : joint(-1), parentBody(-1), childBody(-1), mgm(0) {}
    LoopConstraint(const std::string& type, int jointNum, 
                   int parentBodyNum, int childBodyNum,
                   MultibodyGraphMaker* graphMaker) 
    :   type(type), joint(jointNum), 
        parentBody(parentBodyNum), childBody(childBodyNum), 
        mgm(graphMaker) {}

    /** Get the user reference pointer for the joint associated with this
    loop constraint. **/
    void* getJointRef() const
    {   return mgm->getJoint(joint).userRef; } 
    /** Get the loop constraint type name of the constraint that should be
    used here. **/
    const std::string& getJointTypeName() const
    {   return type; }
    /** Get the user reference pointer for the parent body defined by the
    joint associated with this loop constraint. **/
    void* getParentBodyRef() const
    {   return mgm->getBody(parentBody).userRef; }   
    /** Get the user reference pointer for the child body defined by the
    joint associated with this loop constraint. **/
    void* getChildBodyRef() const
    {   return mgm->getBody(childBody).userRef; }

private:
friend class MultibodyGraphMaker;

    std::string type;        // e.g., ball
    int         joint;       // always one of the input joints
    int         parentBody;  // parent from the joint
    int         childBody;   // child from the joint

    MultibodyGraphMaker*    mgm; // just a reference to container
};


} // namespace SimTK

#endif // SimTK_SIMMATH_MULTIBODY_GRAPH_MAKER_H_

