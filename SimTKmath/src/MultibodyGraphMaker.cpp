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

/* This is the implementation of MultibodyGraphMaker and its auxiliary classes. 
*/

#include "SimTKcommon.h"
#include "simmath/internal/common.h"
#include "simmath/MultibodyGraphMaker.h"

#include <exception>
#include <set>
#include <stdexcept>
#include <string>
#include <iostream>
#include <cstdio>

using std::cout; using std::endl;

#ifdef _MSC_VER
#pragma warning(disable:4996) // don't warn about strcat, sprintf, etc.
#endif

namespace SimTK {

//------------------------------------------------------------------------------
//                              CONSTRUCTOR
//------------------------------------------------------------------------------
MultibodyGraphMaker::MultibodyGraphMaker()
:   weldTypeName("weld"), freeTypeName("free")
{   initialize(); }
 

//------------------------------------------------------------------------------
//                          GET GROUND BODY NAME
//------------------------------------------------------------------------------
const std::string& MultibodyGraphMaker::getGroundBodyName() const {
    if (bodies.empty())
        throw std::logic_error(
            "getGroundBodyName(): you can't call this until you have called "
            "addBody() at least once -- the first body is Ground.");
    return bodies[0].name;
}


//------------------------------------------------------------------------------
//                             ADD JOINT TYPE
//------------------------------------------------------------------------------
int MultibodyGraphMaker::
addJointType(const std::string&     name,
             int                    numMobilities,
             bool                   haveGoodLoopJointAvailable,
             void*                  userRef)
{
    if (!(0 <= numMobilities && numMobilities <= 6)) throw std::runtime_error
       ("addJointType(): Illegal number of mobilities for joint type '" 
        + name + "'");

    // Reject duplicate type name, and reserved type names.
    if (name==getWeldJointTypeName() || name==getFreeJointTypeName())
        throw std::runtime_error("addJointType(): Joint type '" + name 
            + "' is reserved (you can change the reserved names).");   
    std::map<std::string,int>::const_iterator p = jointTypeName2Num.find(name);
    if (p != jointTypeName2Num.end()) throw std::runtime_error
       ("addJointType(): Duplicate joint type '" + name + "'");

    const int jointTypeNum = (int)jointTypes.size(); // next available
    jointTypeName2Num[name] = jointTypeNum; // provide fast name lookup

    jointTypes.push_back(JointType(name, numMobilities, 
                                   haveGoodLoopJointAvailable, userRef));

    return jointTypeNum;
}


//------------------------------------------------------------------------------
//                                 ADD BODY
//------------------------------------------------------------------------------
void MultibodyGraphMaker::addBody(const std::string& name, 
                                  double             mass, 
                                  bool               mustBeBaseBody,
                                  void*              userRef)
{
    if (mass < 0) throw std::invalid_argument
        ("addBody(): Body '" + name + "' specified negative mass");

    // Reject duplicate body name.
    std::map<std::string,int>::const_iterator p = bodyName2Num.find(name);
    if (p != bodyName2Num.end()) throw std::runtime_error
       ("addBody(): Duplicate body name '" + name + "'");

    const int bodyNum = (int)bodies.size(); // next available slot
    bodyName2Num[name] = bodyNum; // provide fast name lookup

    if (bodyNum == 0) { // First body is Ground; use only the name and ref.
        Body ground(name, Infinity, false, userRef);
        ground.level = 0; // already in tree
        bodies.push_back(ground);
    } else { // This is a real body.
        bodies.push_back(Body(name, mass, mustBeBaseBody, userRef));
    }
}

//------------------------------------------------------------------------------
//                                 DELETE BODY
//------------------------------------------------------------------------------
bool MultibodyGraphMaker::deleteBody(const std::string& name)
{
    // Reject non-existing body name.
    std::map<std::string,int>::iterator p = bodyName2Num.find(name);
    if (p == bodyName2Num.end())
        return false;

    const int bodyNum = p->second;

    std::vector<int>& jointAsParent = updBody(bodyNum).jointsAsParent;
    while (jointAsParent.size() > 0)
        deleteJoint(joints[jointAsParent[0]].name);

    std::vector<int>& jointAsChild = updBody(bodyNum).jointsAsChild;
    while (jointAsChild.size() > 0)
        deleteJoint(joints[jointAsChild[0]].name);

    bodies.erase(bodies.begin() + bodyNum);
    bodyName2Num.erase(p);

    // Update body indices due to the deletion of this body
    for (unsigned int i = 0; i < joints.size(); ++i) {
        if (joints[i].parentBodyNum > bodyNum)
            --(updJoint(i).parentBodyNum);
        if (joints[i].childBodyNum > bodyNum)
            --(updJoint(i).childBodyNum);
    }

    for (unsigned int i = bodyNum; i < bodies.size(); ++i)
        bodyName2Num[bodies[i].name] = i;

    return true;
}

//------------------------------------------------------------------------------
//                                ADD JOINT
//------------------------------------------------------------------------------
void MultibodyGraphMaker::addJoint(const std::string&  name,
                                   const std::string&  type,
                                   const std::string&  parentBodyName,
                                   const std::string&  childBodyName,
                                   bool                mustBeLoopJoint,
                                   void*               userRef)
{
    // Reject duplicate joint name, unrecognized type or body names.
    std::map<std::string,int>::const_iterator p = jointName2Num.find(name);
    if (p != jointName2Num.end()) throw std::runtime_error
       ("addJoint(): Duplicate joint name '" + name + "'");

    const int typeNum = getJointTypeNum(type);
    if (typeNum < 0) throw std::runtime_error
       ("addJoint(): Joint " + name + " had unrecognized joint type '" 
        + type + "'");

    const int parentBodyNum = getBodyNum(parentBodyName);
    if (parentBodyNum < 0) throw std::runtime_error
       ("addJoint(): Joint " + name + " had unrecognized parent body '" 
        + parentBodyName + "'");

    const int childBodyNum = getBodyNum(childBodyName);
    if (childBodyNum < 0) throw std::runtime_error
       ("addJoint(): Joint " + name + " had unrecognized child body '" 
        + childBodyName + "'");

    const int jointNum = (int)joints.size(); // next available slot
    jointName2Num[name] = jointNum; // provide fast name lookup
  
    joints.push_back(Joint(name, typeNum, parentBodyNum, childBodyNum,
                           mustBeLoopJoint, userRef));

    updBody(parentBodyNum).jointsAsParent.push_back(jointNum);
    updBody(childBodyNum).jointsAsChild.push_back(jointNum);
}

//------------------------------------------------------------------------------
//                                DELETE JOINT
//------------------------------------------------------------------------------
bool MultibodyGraphMaker::deleteJoint(const std::string&  name)
{
    // Reject duplicate joint name, unrecognized type or body names.
    std::map<std::string,int>::iterator p = jointName2Num.find(name);
    if (p == jointName2Num.end())
        return false;

    const int jointNum = p->second;
    jointName2Num.erase(p);

    std::vector<int>& jointsAsParent = 
        updBody(joints[jointNum].parentBodyNum).jointsAsParent;
    std::vector<int>::iterator it = 
        std::find(jointsAsParent.begin(), jointsAsParent.end(), jointNum);
    if (it == jointsAsParent.end()) throw std::runtime_error
        ("deleteJoint(): Joint " + name + 
         " doesn't exist in jointsAsParent of parent body ");
    jointsAsParent.erase(it);

    std::vector<int>& jointsAsChild = 
        updBody(joints[jointNum].childBodyNum).jointsAsChild;
    it = std::find(jointsAsChild.begin(), jointsAsChild.end(), jointNum);
    if (it == jointsAsChild.end()) throw std::runtime_error
        ("deleteJoint(): Joint " + name + 
         " doesn't exist in jointsAsChild of child body ");
    jointsAsChild.erase(it);

    joints.erase(joints.begin() + jointNum);

    // Update indices due to the deletion of this joint
    for (unsigned int i = 0; i < bodies.size(); ++i) {

        for (it = updBody(i).jointsAsParent.begin(); it != updBody(i).jointsAsParent.end(); ++it)
            if (*it > jointNum) {
                --(*it);
                jointName2Num[joints[*it].name] = *it;
            }

        for (it = updBody(i).jointsAsChild.begin(); it != updBody(i).jointsAsChild.end(); ++it)
            if (*it > jointNum) {
                --(*it);
                // no need to adjust jointName2Num here since the first loops has done it already
            }

    }
    return true;
}

//------------------------------------------------------------------------------
//                              GENERATE GRAPH
//------------------------------------------------------------------------------
void MultibodyGraphMaker::generateGraph() {
    // Body and Joint inputs have been supplied. Body 0 is Ground.

    // 1. Add free joints to Ground for any body that didn't appear
    // in any joint, and for any bodies that said they must be base bodies and
    // don't already have a connection to Ground.
    //
    // While we're at it, look for dangling massless bodies -- they are only
    // allowed if they were originally connected to at least two other bodies
    // by given tree-eligible joints, or if they are connected by a weld 
    // joint and thus have no mobility.
    // TODO: "must be loop joint" joints shouldn't count but we're not
    // checking.
    for (int bn=1; bn < getNumBodies(); ++bn) { // skip Ground
        const Body& body = getBody(bn);
        const int nJoints = body.getNumJoints();
        if (body.mass == 0) {
            if (nJoints == 0) {
                throw std::runtime_error(
                    "generateGraph(): body " + body.name + 
                    " is massless but free (no joint). Must be internal or"
                    " welded to a massful body.");
            }
            if (nJoints == 1) {
                const int jnum = body.jointsAsChild.empty() 
                                    ? body.jointsAsParent[0] 
                                    : body.jointsAsChild[0];
                const Joint& joint = getJoint(jnum);
                const JointType& jtype = getJointType(joint.jointTypeNum);
                if (jtype.numMobilities > 0) {
                    throw std::runtime_error(
                        "generateGraph(): body " + body.name + 
                        " is massless but not internal and not welded to"
                        " a massful body.");
                }                               
            }
        }
        if (nJoints==0 || (body.mustBeBaseBody && !bodiesAreConnected(bn, 0))) 
            connectBodyToGround(bn);
    }

    // 2. Repeat until all bodies are in the tree:
    //   - try to build the tree from Ground outwards
    //   - if incomplete, add one missing connection to Ground
    // This terminates because we add at least one body to the tree each time.
    while (true) {
        growTree();
        const int newBaseBody = chooseNewBaseBody();
        if (newBaseBody < 0) 
            break; // all bodies are in the tree
        // Add joint to Ground.
        connectBodyToGround(newBaseBody);
    }

   // 3. Split the loops
   // This requires adding new "slave" bodies to replace the child bodies in
   // the loop-forming joints. Then each slave is welded back to its master
   // body (the original child).
   breakLoops();
}


//------------------------------------------------------------------------------
//                              CLEAR GRAPH
//------------------------------------------------------------------------------
void MultibodyGraphMaker::clearGraph() {
    for (int bn=1; bn < getNumBodies(); ++bn) {  // skip Ground
        if (bodies[bn].isSlave()) {
            // Assumption: all slave bodies are clustered at the end of the body array
            bodies.erase(bodies.begin() + bn, bodies.begin() + bodies.size());
            break;
        }
        updBody(bn).forgetGraph(*this);
    }

    for (int jn=0; jn < getNumJoints(); ++jn)
        if (updJoint(jn).forgetGraph(*this))
            --jn;

    mobilizers.clear();
    constraints.clear();
}

//------------------------------------------------------------------------------
//                               DUMP GRAPH
//------------------------------------------------------------------------------
void MultibodyGraphMaker::dumpGraph(std::ostream& o) const {
    char buf[1024];
    o << "\nMULTIBODY GRAPH\n";
    o <<   "---------------\n";
    o << "\n" << getNumBodies() << " BODIES:\n";
    for (int i=0; i < getNumBodies(); ++i) {
        const MultibodyGraphMaker::Body& body  = getBody(i);
        sprintf(buf, "%2d %2d: %s mass=%g mob=%d master=%d %s\n", 
            i, body.level,
            body.name.c_str(), body.mass, body.mobilizer, body.master,
            body.mustBeBaseBody?"MUST BE BASE BODY":"");
        o << buf;
        o << "  jointsAsParent=[";
        for (unsigned j=0; j<body.jointsAsParent.size(); ++j)
            o << " " << body.jointsAsParent[j];
        o << "]\t  jointsAsChild=[";
        for (unsigned j=0; j<body.jointsAsChild.size(); ++j)
            o << " " << body.jointsAsChild[j];
        o << "]\t  slaves=[";
        for (unsigned j=0; j<body.slaves.size(); ++j)
            o << " " << body.slaves[j];
        o << "]\n";
    }

    o << "\n" << getNumJoints() << " JOINTS:\n";
    for (int i=0; i < getNumJoints(); ++i) {
        const MultibodyGraphMaker::Joint& joint = getJoint(i);
        sprintf(buf, "%2d %2d: %20s %20s->%-20s %10s loopc=%2d %s %s\n",
            i, joint.mobilizer, joint.name.c_str(),
            getBody(joint.parentBodyNum).name.c_str(),
            getBody(joint.childBodyNum).name.c_str(),
            getJointType(joint.jointTypeNum).name.c_str(),
            joint.loopConstraint,
            joint.mustBeLoopJoint?"MUST BE LOOP":"",
            joint.isAddedBaseJoint?"ADDED BASE JOINT":"");
        o << buf;
    }

    o << "\n" << getNumMobilizers() << " MOBILIZERS:\n";
    for (int i=0; i < getNumMobilizers(); ++i) {
        const MultibodyGraphMaker::Mobilizer& mo    = getMobilizer(i);
        const MultibodyGraphMaker::Joint&     joint = getJoint(mo.joint);
        const MultibodyGraphMaker::Body&      inb   = getBody(mo.inboardBody);
        const MultibodyGraphMaker::Body&      outb  = getBody(mo.outboardBody);
        sprintf(buf, "%2d %2d: %20s %20s->%-20s %10s %2d %3s\n",
            i, mo.level, joint.name.c_str(), 
            inb.name.c_str(), outb.name.c_str(), 
            mo.getJointTypeName().c_str(), 
            mo.joint, mo.isReversed?"REV":"");
        o << buf;
    }

    o << "\n" << getNumLoopConstraints() << " LOOP CONSTRAINTS:\n";
    for (int i=0; i < getNumLoopConstraints(); ++i) {
        const MultibodyGraphMaker::LoopConstraint& lc = getLoopConstraint(i);
        const MultibodyGraphMaker::Body parent = getBody(lc.parentBody);
        const MultibodyGraphMaker::Body child  = getBody(lc.childBody);
        sprintf(buf, "%d: %s parent=%s child=%s jointNum=%d\n",
            i, lc.type.c_str(), parent.name.c_str(), child.name.c_str(),
            lc.joint);
        o << buf;
    }
    o << "\n----- END OF MULTIBODY GRAPH.\n\n";
}


//------------------------------------------------------------------------------
//                               INITIALIZE
//------------------------------------------------------------------------------
// Clear all data and prime with free and weld joint types.
void MultibodyGraphMaker::initialize() {
    clear();
    jointTypes.push_back(JointType(weldTypeName, 0, true, NULL));
    jointTypes.push_back(JointType(freeTypeName, 6, true, NULL));
    jointTypeName2Num[weldTypeName] = 0;
    jointTypeName2Num[freeTypeName] = 1;
}


//------------------------------------------------------------------------------
//                                SPLIT BODY
//------------------------------------------------------------------------------
// Create a new slave body for the given master, and add it to the list of
// bodies. Does not create the related loop constraint. The body
// number assigned to the slave is returned.
int MultibodyGraphMaker::splitBody(int masterBodyNum) {
    const int slaveBodyNum = (int)bodies.size(); // next available
    Body& master = updBody(masterBodyNum);
    // First slave is number 1, slave 0 is the master.
    std::stringstream ss;
    ss << "#" << master.name << "_slave_" << master.getNumSlaves()+1;
    Body slave(ss.str(), NaN, false, NULL);
    slave.master = masterBodyNum;
    master.slaves.push_back(slaveBodyNum);
    bodies.push_back(slave);
    // No name lookup for slave bodies.
    return slaveBodyNum;
}


//------------------------------------------------------------------------------
//                          CHOOSE NEW BASE BODY
//------------------------------------------------------------------------------
// We've tried to build the tree but might not have succeeded in using
// all the bodies. That means we'll have to connect one of them to Ground. 
// Select the best unattached body to use as a base body, or -1 if all bodies
// are already included in the spanning tree. We hope to find a body that is
// serving as a parent but has never been listed as a child; that is crying out
// to be connected directly to Ground. Failing that, we'll pick a body that
// has been used as a parent often.
int MultibodyGraphMaker::chooseNewBaseBody() const {
    bool parentOnlyBodySeen = false;
    int bestBody = -1; int nChildren=-1;
    // Skip the Ground body.
    for (int bx=1; bx < getNumBodies(); ++bx) {
        const Body& body = bodies[bx];
        if (body.isInTree()) 
            continue; // unavailable
        if (parentOnlyBodySeen && !body.jointsAsChild.empty()) 
            continue; // already seen parent-only bodies; no need for this one
        if (!parentOnlyBodySeen && body.jointsAsChild.empty()) {
            // This is our first parent-only body; it is automatically the
            // best candidate now.
            parentOnlyBodySeen = true;
            bestBody = bx; nChildren = (int)body.jointsAsParent.size();
        } else { // Keep the body that has the most children.
            if ((int)body.jointsAsParent.size() > nChildren) {
                bestBody  = bx; 
                nChildren = (int)body.jointsAsParent.size();
            }
        }
    }
    return bestBody;
}


//------------------------------------------------------------------------------
//                          CONNECT BODY TO GROUND
//------------------------------------------------------------------------------
// Connect the given body to Ground by a free joint with parent Ground and
// child the given body. This makes the body a base body (level 1 body)
// for some subtree of the multibody graph.
void MultibodyGraphMaker::connectBodyToGround(int bodyNum) {
    const Body& body = getBody(bodyNum);
    assert(!body.isInTree());
    const std::string jointName = "#" + getGroundBodyName() + "_" + body.name;
    addJoint(jointName, getFreeJointTypeName(), 
                    getGroundBodyName(), body.name, false, NULL);
    updJoint(jointName).isAddedBaseJoint = true;
}


//------------------------------------------------------------------------------
//                          ADD MOBILIZER FOR JOINT
//------------------------------------------------------------------------------
// We've already determined this joint is eligible to become a mobilizer; now
// make it so. Given a joint for which either the parent or child is in the 
// tree, but not both, create a mobilizer implementing this joint and attaching 
// the unattached body (which will be the mobilizer's outboard body) to the
// tree. The mobilizer number is returned and also recorded in the joint
// and outboard body.
int MultibodyGraphMaker::addMobilizerForJoint(int jointNum) {
    Joint& joint = updJoint(jointNum);
    assert(!joint.mustBeLoopJoint);   // can't be a tree joint
    assert(!joint.hasMobilizer());      // already done
    const int pNum = joint.parentBodyNum, cNum = joint.childBodyNum;
    Body& parent = updBody(pNum);
    Body& child  = updBody(cNum);
    // Exactly one of these must already be in the tree.
    assert(parent.isInTree() ^ child.isInTree());

    const int mobNum = (int)mobilizers.size(); // next available
    if (parent.isInTree()) { 
        // Child is outboard body (forward joint)
        child.level = parent.level + 1;
        mobilizers.push_back(Mobilizer(jointNum,child.level,
                                       pNum,cNum,false,this));
        child.mobilizer = mobNum;
    } else if (child.isInTree()) { 
        // Parent is outboard body (reverse joint)
        parent.level = child.level + 1;
        mobilizers.push_back(Mobilizer(jointNum,parent.level,
                                       cNum,pNum,true,this));
        parent.mobilizer = mobNum;
    }
    joint.mobilizer = mobNum;
    return mobNum;
}


//------------------------------------------------------------------------------
//                   FIND HEAVIEST UNASSIGNED FORWARD JOINT
//------------------------------------------------------------------------------
// Starting with a given body that is in the tree already, look at its 
// unassigned children (meaning bodies connected by joints that aren't 
// mobilizers yet) and return the joint connecting it to the child with the 
// largest mass. (The idea is that this would be a good direction to extend the 
// chain.) Return -1 if this body has no children.
int MultibodyGraphMaker::
findHeaviestUnassignedForwardJoint(int inboardBody) const {
    const Body& inb = getBody(inboardBody);
    int jointNum = -1;
    double maxMass=0;
    // Search joints for which this is the parent body.
    for (unsigned i=0; i < inb.jointsAsParent.size(); ++i) {
        const int jfwd = inb.jointsAsParent[i];
        const Joint& joint = joints[jfwd];
        if (joint.hasMobilizer()) continue; // already in the tree
        if (joint.mustBeLoopJoint) continue; // can't be a tree joint
        const Body& child = getBody(joint.childBodyNum);
        if (child.isInTree()) continue; // this is a loop joint
        if (child.mass > maxMass) {
            jointNum = jfwd;
            maxMass = child.mass;
        }
    }
    return jointNum;
}


//------------------------------------------------------------------------------
//                   FIND HEAVIEST UNASSIGNED REVERSE JOINT
//------------------------------------------------------------------------------
// Same as previous method, but now we are looking for unassigned joints where 
// the given body serves as the child so we would have to reverse the joint to 
// add it to the tree. (This is less desirable so is a fallback.)
int MultibodyGraphMaker::
findHeaviestUnassignedReverseJoint(int inboardBody) const {
    const Body& inb = getBody(inboardBody);
    int jointNum = -1;
    double maxMass=0;
    // Search joints for which this is the child body.
    for (unsigned i=0; i < inb.jointsAsChild.size(); ++i) {
        const int jrev = inb.jointsAsChild[i];
        const Joint& joint = joints[jrev];
        if (joint.hasMobilizer()) continue; // already in the tree
        if (joint.mustBeLoopJoint) continue; // can't be a tree joint
        const Body& parent = getBody(joint.parentBodyNum);
        if (parent.isInTree()) continue; // this is a loop joint
        if (parent.mass > maxMass) {
            jointNum = jrev;
            maxMass = parent.mass;
        }
    }
    return jointNum;
}


//------------------------------------------------------------------------------
//                                 GROW TREE
//------------------------------------------------------------------------------
// Process unused joints for which one body is in the tree (at level h) and the 
// other is not. Add the other body to the tree at level h+1, marking the joint 
// as forward (other body is child) or reverse (other body is parent). Repeat 
// until no changes are made. Does not assign loop joints or any bodies that 
// don't have a path to Ground. Extend the multibody tree starting with the 
// lowest-level eligible joints and moving outwards. We're doing this 
// breadth-first so that we get roughly even-length chains and we'll stop at 
// the first level at which we fail to find any joint. If we fail to consume 
// all the bodies, the caller will have to add a level 1 joint to attach a 
// previously-disconnected base body to Ground and then call this again.
//
// We violate the breadth-first heuristic to avoid ending a branch with a 
// massless body, unless it is welded to its parent. If we add a mobile massless
// body, we'll keep going out along its branch until we hit a massful body. It 
// is a disaster if we fail to find a massful body because a tree that 
// terminates in a mobile massless body will generate a singular mass matrix. 
// We'll throw an exception in that case, but note that this may just be a 
// failure of the heuristic; there may be some tree that could have avoided the 
// terminal massless body but we failed to discover it.
//
// TODO: keep a list of unprocessed joints so we don't have to go through
// all of them again at each level.
void MultibodyGraphMaker::growTree() {
    // Record the joints for which we added mobilizers during this subtree
    // sweep. That way if we jumped levels ahead due to massless bodies we
    // can take credit and keep going rather than quit.
    std::set<int> jointsAdded;
    for (int level=1; ;++level) { // level of outboard (mobilized) body
        bool anyMobilizerAdded = false;
        for (int jNum=0; jNum<getNumJoints(); ++jNum) {
            // See if this joint is at the right level (meaning its inboard
            // body is at level-1) and is available to become a mobilizer.
            const Joint& joint = getJoint(jNum);
            if (joint.hasMobilizer()) {
                // Already done -- one of ours?
                if (jointsAdded.count(jNum)) {
                    // We added it during this growTree() call -- is it at
                    // the current level?
                    const Mobilizer& mob = mobilizers[joint.mobilizer];
                    if (mob.getLevel() == level)
                        anyMobilizerAdded = true; // Added one at level.
                }
                continue; 
            }
            if (joint.mustBeLoopJoint) continue; // can't be a tree joint
            const Body& parent = getBody(joint.parentBodyNum);
            const Body& child  = getBody(joint.childBodyNum);
            // Exactly one of parent or child must already be in the tree.
            if (!(parent.isInTree() ^ child.isInTree())) continue;
            if (parent.isInTree()) {
                if (parent.level != level-1) continue; // not time yet
            } else { // child is in tree
                if (child.level != level-1) continue; // not time yet
            } 
            addMobilizerForJoint(jNum);
            jointsAdded.insert(jNum);
            anyMobilizerAdded = true;

            // We just made joint jNum a mobilizer. If its outboard body 
            // is massless and the mobilizer was not a weld, we need to keep 
            // extending this branch of the tree until we can end it with a 
            // massful body.
            const Body& outboard = getBody(mobilizers.back().outboardBody);
            const JointType& jtype = getJointType(joint.jointTypeNum);
            if (jtype.numMobilities == 0 || outboard.mass > 0)
                continue;

            // Pick a further-outboard body with mass and extend branch to it.
            // Prefer forward-direction joints if possible.
            // If only massless outboard bodies are available, add one and
            // keep going.
            while (true) {
                const int bNum = mobilizers.back().outboardBody;
                const int jfwd = findHeaviestUnassignedForwardJoint(bNum);
                if (jfwd>=0 && getBody(getJoint(jfwd).childBodyNum).mass > 0) {
                    addMobilizerForJoint(jfwd);
                    jointsAdded.insert(jfwd);
                    break;
                }
                const int jrev = findHeaviestUnassignedReverseJoint(bNum);
                if (jrev>=0 && getBody(getJoint(jrev).parentBodyNum).mass > 0) {
                    addMobilizerForJoint(jrev);
                    jointsAdded.insert(jrev);
                    break;
                }

                // Couldn't find a massful body to add. Add another massless 
                // body (if there is one) and keep trying.
                if (jfwd >= 0) {
                    addMobilizerForJoint(jfwd);
                    jointsAdded.insert(jfwd);
                    continue;
                }
                if (jrev >= 0) {
                    addMobilizerForJoint(jrev);
                    jointsAdded.insert(jrev);
                    continue;
                }

                // Uh oh. Nothing outboard of the massless body we just added.
                // We're ending this chain with a massless body; not good.
                throw std::runtime_error("growTree(): algorithm produced an"
                    " invalid tree containing a terminal massless body ("
                    + getBody(bNum).name + "). This may be due to an invalid" 
                    " model or failure of heuristics. In the latter case you"
                    " may be able to get a valid tree by forcing a different"
                    " loop break or changing parent->child ordering.");
            }
        }
        if (!anyMobilizerAdded) 
            break;
    }
}


//------------------------------------------------------------------------------
//                                 BREAK LOOPS
//------------------------------------------------------------------------------
// Find any remaining unused joints, which will have both parent and
// child bodies already in the tree. This includes joints that the user
// told us must be loop joints, and ones picked by the algorithm. 
// For each of those, implement the joint with a loop constraint if one
// is provided, otherwise implement it with a mobilizer and split off a 
// slave body from the child body, use that slave as the outboard body
// for the mobilizer, and add a weld constraint to reconnect the slave to
// its master (the original child body).
void MultibodyGraphMaker::breakLoops() {
    for (int jx=0; jx < getNumJoints(); ++jx) {
        Joint& jinfo = joints[jx];
        if (jinfo.hasMobilizer()) continue; // already done
        const int px = jinfo.parentBodyNum, cx = jinfo.childBodyNum;
        assert(getBody(px).isInTree() && getBody(cx).isInTree());

        const JointType& jtype = getJointType(jinfo.jointTypeNum);
        if (jtype.haveGoodLoopJointAvailable) {
            const int loopNum = (int)constraints.size();
            constraints.push_back(LoopConstraint(jtype.name,jx,px,cx,this));
            jinfo.loopConstraint = loopNum;
            continue;
        }

        // No usable loop constraint for this type of joint. Add a new slave 
        // body so we can use a mobilizer. (Body references are invalidated here
        // since we're adding a new body -- don't reuse them.)
        const int sx = splitBody(cx);
        Body& parent = updBody(px);
        Body& child  = updBody(cx);
        Body& slave  = updBody(sx);

        // Mobilize the slave body.
        const int mobNum = (int)mobilizers.size(); // next available
        const int level = parent.level+1;
        jinfo.mobilizer = slave.mobilizer = mobNum;
        slave.jointsAsChild.push_back(jx);
        mobilizers.push_back(Mobilizer(jx,level,px,sx,false,this));
        slave.level = level;
    }
}


//------------------------------------------------------------------------------
//                         BODIES ARE CONNECTED
//------------------------------------------------------------------------------
// Return true if there is a joint between two bodies given by body number,
// regardless of parent/child ordering.
bool MultibodyGraphMaker::bodiesAreConnected(int b1, int b2) const {
    const Body& body1 = getBody(b1);
    for (unsigned i=0; i < body1.jointsAsParent.size(); ++i)
        if (getJoint(body1.jointsAsParent[i]).childBodyNum == b2)
            return true;
    for (unsigned i=0; i < body1.jointsAsChild.size(); ++i)
        if (getJoint(body1.jointsAsChild[i]).parentBodyNum == b2)
            return true;
    return false;
}

//------------------------------------------------------------------------------
//                         FORGET GRAPH
//------------------------------------------------------------------------------
// Restore the Body to its state prior to generateGraph()
void MultibodyGraphMaker::Body::forgetGraph(MultibodyGraphMaker& graph)
{   
    level = -1;
    mobilizer = -1;
    master = -1;
    slaves.clear();
}

//------------------------------------------------------------------------------
//                         FORGET GRAPH
//------------------------------------------------------------------------------
// Restore the Joint to its state prior to generateGraph()
bool MultibodyGraphMaker::Joint::forgetGraph(MultibodyGraphMaker& graph)
{
    if (isAddedBaseJoint) {
        graph.deleteJoint(name);
        return true;
    }
    mobilizer = -1;
    loopConstraint = -1;
    return false;
}
} // namespace SimTK

