/* -------------------------------------------------------------------------- *
 *                   Simbody(tm) Example: Gazebo 2 Simbody                    *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2013 Stanford University and the Authors.           *
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

/* This example reads input in the form of a Gazebo ".sdf" file and generates
a roughly-equivalent Simbody model. Gazebo is a robot simulator from the Open 
Source Robotics Foundation (http://osrfoundation.org).

An sdf file describes a model in XML format using the following objects:
world   the Ground frame. May have associated geometry and lights. Contains one
        or more "model" objects. Normally the ground is the x-y plane, with
        +z the "up" direction and a ground plane at z=0. The world also acts
        as a link whose name is "world".
model   a named grouping of physical objects. Provides a model frame given
        relative to the world frame, however this frame is still fixed in the
        world. A model contains links and joints.
link    a named mass- and geometry-carrying object (a body). There is a link frame,
        and inertial, visual, and collision objects each with their own frame
        given relative to the link frame. An initial pose for the link frame
        is given, relative to the model frame.
joint   a connection between two links or between a link and the world
        There is a parent and child link, given by name, with "world"
        the name for the ground link. There is a joint type and then appropriate
        parameters depending on the type. There is a child frame and a
        parent frame. However, if only a single frame is given for the joint
        it is interpreted as the child frame and then the default link frame
        poses are used to find the coincident frame on the parent.

A <pose> element (containing six scalars) is equivalent to a Simbody Transform,
interpreted as x,y,z translation vector followed by an x-y-z body-fixed Euler 
angle sequence given in radians.

A <visual> element has <geometry> and <material> subelements; we'll use the
material only to pick a color.

A <collision> element has <geometry> and <surface> subelements, with the
latter giving the physical properties used in contact. We'll generate a gray
visual for these.

We have to take the links and joints and construct a tree whose root is the
world, with every link appearing exactly once. Any link that does not appear
in any joint will be given six degrees of freedom relative to the world (3 dofs 
if the link is a point mass). If there are disconnected branches we will pick
one of the bodies as a "base" body and connect it to world by a free (6dof)
joint.

We can't take the parent-child specification too seriously when building the
tree since these can be arbitrary in the Gazebo description. To avoid confusion, 
we'll speak of the "inboard" and "outboard" bodies of a mobilizer; usually 
inboard==parent and outboard==child but the relationship is opposite for a 
reverse mobilizer.For example, a three-body chain might be defined
     link1--j1--link2--j2--link3
     joint1 parent=link1 child=link2
     joint2 parent=link3 child=link2
There is no way to build a tree where inboard/outboard matches parent/child
for both joints in this example. So we will try to
build the tree respecting the specified parent-child ordering, but if we 
have to reverse we'll do that with the Simbody reverse mobilizer capability so
that the meaning of the joint is unchanged. 

After the tree is built, all the links will have been used but there may still
be some unused joints. Those joints will involve parent and child links both
of which are already in the tree. That means these joints form topological
loops in the multibody graph. For now we'll handle that uniformly by 
splitting the child link to create a new "shadow" body and then introducing a 
weld constraint to hold the master link and its shadows together. Note that
the visual and contact geometry stays with the master; the mass and inertia
are divided equally among the master and shadows.
TODO: if the loop-closing joint is a Ball, use a Ball constraint directly
(note: a loop-closing Free joint doesn't do anything)

We will construct the Simbody model here and maintain a mapping between the
Gazebo objects and the Simbody implementation of them. 
*/
#include "Simbody.h"

#include <utility>
#include <map>
#include <iostream>
using std::cout; using std::endl; using std::string;
using std::map; using std::vector;

using namespace SimTK;

//==============================================================================
//                                LINK INFO
//==============================================================================
// There is one of these for every mobilized body in the multibody tree that
// we build. That includes links from the Gazebo input file, and "shadow" links
// that we add to deal with closed topological loops. Shadow bodies are always
// terminal bodies, appearing as the child body in just one joint. Every shadow
// is connected to its master body by a weld constraint. The master and each
// of its shadows have identical mass properties.
class LinkInfo {
public:
    explicit LinkInfo(const string& name) : name(name) {clear();}


    bool isInTree() const {return level>=0;}
    int getNumShadows() const 
    {   assert(!isShadow()); return (int)shadowLinks.size(); }
    int isShadow() const {return masterIx >= 0;}

    void clear() {
        jointsAsParent.clear(); jointsAsChild.clear(); 
        level=mobilizer=masterIx=-1; 
        shadowLinks.clear();
        mass=NaN; com_L=Vec3(0); uinertia_L=UnitInertia(0);
    }

    Xml::Element    element;        // if this was in the Xml file
    const string    name;

    vector<int>     jointsAsParent; // joints where this link is the parent
    vector<int>     jointsAsChild;  // joints where this link is the child

    int             level; // world=0, connected to world=1, contact to that=2, etc.
    int             mobilizer; // this link's connection to tree

    int             masterIx;    // >=0 if this is a shadow
    vector<int>     shadowLinks; // shadow links, if this is a master

    // Data converted to Simbody terms. The Gazebo link frame L is the same 
    // as the Simbody body frame B. Gazebo's inertia has to be transformed to
    // be about Bo and expressed in B; it comes in with its own inertial frame
    // I with origin Io placed at the mass center. We are given pose X_LI.
    Real            mass; // total mass, master+shadows
    Vec3            com_L;  // com location in link frame
    UnitInertia     uinertia_L; // inertia about L origin, in L
                                // master,shadow have inertia=(mass/n)*uinertia

    Transform       defX_GL;    // default pose in world frame
};

//==============================================================================
//                               JOINT INFO
//==============================================================================
// There is one of these for every joint in the system. Most are read in from
// the Gazebo input file, but free joints will be added where there are missing 
// connections to the world. Joints that we decide are the loop closers will
// have their Gazebo child links replaced with shadow links of that child. Every 
// joint here will correspond to a Simbody mobilizer; the outboard body listed
// here is the Simbody mobilized body. (Note that "outboard" is not always the
// same as "child".)
struct JointInfo {
    explicit JointInfo(const string& name, const string& type) 
    :   name(name), type(type) {clear();}
    void clear() {
        parentIx=childIx=masterIx=level=-1; 
        isReversed=false; 
    }
    int getInboardLink()  const {return isReversed?childIx:parentIx;}
    int getOutboardLink() const {return isReversed?parentIx:childIx;}

    bool isInTree() const {return level>=0;}
    bool isLoopCloser() const {return masterIx >= 0;}

    Xml::Element    element; // if this was in the Xml file
    const string    name, type;
    int             parentIx, childIx;
    int             masterIx; // if childIx is shadow

    int             level; // level of outboard body
    bool            isReversed; // true if parent is outboard body

    // Normally A=F, B=M. But if reversed, then B=F, A=M.
    Transform       X_PA; // parent body frame to mobilizer frame
    Transform       X_CB; // child body frame to mobilizer frame
    Transform    defX_AB; // default mobilizer pose
};

//==============================================================================
//                                 LINKS
//==============================================================================
// This collects all the links and supports the steps to building the
// multibody tree.
struct Links {
    // Add a new link to the collection and index it by name.
    int addLink(const LinkInfo& info);

    // Split off a shadow link for the link given by its index. The index
    // of the shadow is returned.
    int splitLink(int lx);

    // We've tried to build the tree but might not have succeeded in using
    // all the bodies. That means we'll have to connect one of them to world. 
    // Here we try to pick the best one: 
    //  - if any unused links were only used as parents, pick the one
    //    that serves as parent for the most links
    //  - otherwise pick the one that serves as parent most often
    int chooseNewBaseBody() const;

    bool hasLink(const string& name) const 
    {   return name2index.find(name) != name2index.end(); }
    int getLinkIndex(const string& name) const
    {   assert(hasLink(name)); return name2index.find(name)->second; }
    LinkInfo& getLink(const string& name) {return getLink(getLinkIndex(name));}
    LinkInfo& getLink(int index) {return linkByIndex[index];}

    int size() const {return (int)linkByIndex.size();}

    std::vector<LinkInfo>   linkByIndex;
    std::map<string,int>    name2index;
};

//==============================================================================
//                                 JOINTS
//==============================================================================
class Joints {
public:
    // Add a new joint to the collection and index it by name.
    int addJoint(const JointInfo& info);

    // Connect the given link to world by a free joint with parent world and
    // child the new link.
    // TODO: check mass properties and use translation joint if inertialess.
    int addJointToWorld(Links& links, int linkIndex);

    // Process unused joints for which one link is in the tree (at level h)
    // and the other is not. Add the other body to the tree at level h+1,
    // marking the joint as forward (other body is child) or reverse (other
    // body is parent). Repeat until no changes are made. Does not assign
    // loop joints or any links that don't have a path to world.
    void makeTree(Links& links);

    // Find any remaining unused joints, which will have both parent and
    // child links already in the tree. For each of those, split off a shadow
    // for the child link and use that shadow as the outboard body.
    void breakLoops(Links& links);

    bool hasJoint(const string& name) const 
    {   return name2index.find(name) != name2index.end(); }
    int getJointIndex(const string& name) const
    {   assert(hasJoint(name)); return name2index.find(name)->second; }
    JointInfo& getJoint(const string& name) {return getJoint(getJointIndex(name));}
    JointInfo& getJoint(int index) {return jointByIndex[index];}
    int size() const {return (int)jointByIndex.size();}

    std::vector<JointInfo>  jointByIndex;
    std::map<string,int>    name2index;

    // These are the tree joints (mobilizers) ordered in levels from world
    // outwards.
    std::map<int,std::vector<int> > mobilizersByLevel;
};

//==============================================================================
//                                  MAIN
//==============================================================================
int main(int argc, const char* argv[]) {
  try {
    const std::string dir = "/Users/Sherm/Desktop/";
    //Xml::Document sdf(dir + "double_pendulum.sdf");
    //Xml::Document sdf(dir + "double_pendulum2.sdf");
    //Xml::Document sdf(dir + "ragdoll.sdf");
    Xml::Document sdf(dir + "ragdoll2.sdf");

    SimTK_ERRCHK1_ALWAYS(   sdf.getRootTag() == "sdf"
                         || sdf.getRootTag() == "gazebo", "Gazebo2Simbody",
        "Expected to see document tag <sdf> or <gazebo> but saw <%s> instead.",
        sdf.getRootTag().c_str());

    // This is a Gazebo sdf document.
    Xml::Element root = sdf.getRootElement();
    cout << "sdf version=" 
         << root.getOptionalAttributeValue("version", "unspecified") << endl;

    Xml::Element world = root.getRequiredElement("world");
    const Vec3 gravity = 
        world.getOptionalElementValueAs<Vec3>("gravity", Vec3(0,0,-9.81));

    Array_<Xml::Element> models = world.getAllElements("model");
    cout << "World '" << world.getOptionalAttributeValue("name","NONAME")
         << "' contains " << models.size() << " model(s).\n";
    cout << "Gravity=" << gravity << endl;

    if (models.empty()) {
        cout << "File contained no model -- nothing to do. Goodbye.\n";
        return 0;
    }

    MultibodySystem mbs;
    SimbodyMatterSubsystem matter(mbs);
    GeneralForceSubsystem forces(mbs);
    ContactTrackerSubsystem tracker(mbs);
    CompliantContactSubsystem contact(mbs, tracker);
    Force::UniformGravity(forces, matter, gravity);

    Xml::Element model = models[0]; // ignore others for now
    Array_<Xml::Element> linkElts = model.getAllElements("link");
    Array_<Xml::Element> jointElts = model.getAllElements("joint");
    printf("Building model '%s' with ground + %d links, %d joints.\n",
        model.getRequiredAttributeValue("name").c_str(), 
        linkElts.size(), jointElts.size());

    Links links;
    Joints joints;

    // Step 1: create the links collection including a "world" link and
    // all links that appeared in the model. An index number is assigned for
    // each link, with world as link 0. A map from name->index is built.
    LinkInfo worldLink("world");
    worldLink.level = 0;
    links.addLink(worldLink);
    for (unsigned i=0; i < linkElts.size(); ++i) {
        LinkInfo linfo(linkElts[i].getRequiredAttributeValue("name"));
        linfo.element = linkElts[i];
        links.addLink(linfo);
    }

    // Step 2: create the joints collection including at first only the 
    // joints specified in the model. We will also update the Links data
    // to note if and how the links appear in joints.

    for (unsigned i=0; i < jointElts.size(); ++i) {
        const string name = jointElts[i].getRequiredAttributeValue("name");
        const string type = jointElts[i].getRequiredAttributeValue("type");
        JointInfo jinfo(name,type);
        jinfo.element = jointElts[i];

        const string parent = jointElts[i].getRequiredElementValue("parent");
        const string child  = jointElts[i].getRequiredElementValue("child");
        assert(links.hasLink(parent) && links.hasLink(child));
        assert(child != parent);

        jinfo.parentIx = links.getLinkIndex(parent);
        jinfo.childIx  = links.getLinkIndex(child);

        const int jx = joints.addJoint(jinfo);
        links.getLink(jinfo.parentIx).jointsAsParent.push_back(jx);
        links.getLink(jinfo.childIx).jointsAsChild.push_back(jx);
    }

    // Step 3: add free joints to world for any link that didn't appear
    // in any joint (except world!).
    for (int lx=1; lx < links.size(); ++lx) { 
        LinkInfo& link = links.getLink(lx);
        if (link.jointsAsChild.empty() && link.jointsAsParent.empty()) 
            joints.addJointToWorld(links,lx);
    }

    // Step 4: repeat until done:
    //   - try to build the tree
    //   - add one missing connection to world

    while (true) {
        joints.makeTree(links);
        int newBaseBody = links.chooseNewBaseBody();
        if (newBaseBody < 0) 
            break; // all links are in use
        // Add joint to world.
        printf("Adding new base body %d\n", newBaseBody);
        joints.addJointToWorld(links, newBaseBody);
    }

   // Step 5: split the loops
   joints.breakLoops(links);

   // All joints should now be in the tree. 

    printf("links n=%d:\n", links.size());
    for (int lx=0; lx < links.size(); ++lx) { 
        const LinkInfo& link = links.getLink(lx);
        printf("%2d %2d: %15s mob=%2d", link.level, lx, link.name.c_str(),
            link.mobilizer); 
        printf(" childIn=");
        for (int cx=0; cx<(int)link.jointsAsChild.size(); ++cx)
            printf("%2d ", link.jointsAsChild[cx]);
        printf(" parentIn=");
        for (int px=0; px<(int)link.jointsAsParent.size(); ++px)
            printf("%2d ", link.jointsAsParent[px]);
        printf("\n");
    }   

    printf("mobilizers, max level=%d:\n", joints.mobilizersByLevel.size());
    std::map<int,std::vector<int> >::const_iterator p = 
        joints.mobilizersByLevel.begin();
    int nMobilizers = 0;
    for (; p != joints.mobilizersByLevel.end(); ++p) { 
        for (unsigned i=0; i < p->second.size(); ++i) {
            ++nMobilizers;
            const int jx = p->second[i];
            const JointInfo& joint = joints.getJoint(jx);
            assert(joint.level == p->first);
            printf("%2d %2d: %2d->%2d %10s %20s", joint.level,
                jx, joint.parentIx, joint.childIx,
                joint.type.c_str(), joint.name.c_str());
            printf(" rev=%d inb=%2d outb=%2d",
                joint.isReversed, joint.getInboardLink(), joint.getOutboardLink());
            if (joint.isLoopCloser())
                printf(" LOOP");
            printf("\n");
        }
    }  
    printf("%d mobilizers in tree\n", nMobilizers);

  } catch (const std::exception& e) {
    cout << "EXCEPTION: " << e.what() << "\n";
    cout << "Working directory: " << Pathname::getCurrentWorkingDirectory() 
         << endl;
  }

    return 0;
}

//------------------------------------------------------------------------------
//                          LINKS IMPLEMENTATION
//------------------------------------------------------------------------------
int Links::addLink(const LinkInfo& info) {
    int lx = linkByIndex.size();
    linkByIndex.push_back(info);
    SimTK_ERRCHK1_ALWAYS(name2index.find(info.name)==name2index.end(),
        "Gazebo2Simbody", "Link name %s was used more than once.", 
        info.name.c_str());
    name2index[info.name] = lx;
    return lx;
}

int Links::splitLink(int lx) {
    int sx = linkByIndex.size(); // next available
    LinkInfo& masterLink = linkByIndex[lx];
    const string& name = masterLink.name;
    // First shadow is number 1, shadow 0 is the master.
    char snum[8]; sprintf(snum, "%d_", masterLink.getNumShadows()+1);
    LinkInfo shadowInfo("#shadow" + string(snum) + name);
    shadowInfo.masterIx = lx;
    masterLink.shadowLinks.push_back(sx);
    linkByIndex.push_back(shadowInfo);
    // no name lookup for shadow links
    return sx;
}

int Links::chooseNewBaseBody() const {
    bool parentOnlyLinkSeen = false;
    int bestLink = -1; int nChildren=-1;
    for (int lx=0; lx < size(); ++lx) {
        const LinkInfo& linfo = linkByIndex[lx];
        if (linfo.isInTree()) continue;
        if (parentOnlyLinkSeen && !linfo.jointsAsChild.empty()) continue;
        if (!parentOnlyLinkSeen && linfo.jointsAsChild.empty()) {
            parentOnlyLinkSeen = true;
            bestLink = lx; nChildren = linfo.jointsAsParent.size();
        } else {
            if ((int)linfo.jointsAsParent.size() > nChildren) {
                bestLink  = lx; 
                nChildren = linfo.jointsAsParent.size();
            }
        }
    }
    return bestLink;
}

//------------------------------------------------------------------------------
//                          JOINTS IMPLEMENTATION
//------------------------------------------------------------------------------
int Joints::addJoint(const JointInfo& info) {
    unsigned index = jointByIndex.size();
    jointByIndex.push_back(info);
    SimTK_ERRCHK1_ALWAYS(name2index.find(info.name)==name2index.end(),
        "Gazebo2Simbody", "Joint name %s was used more than once.", 
        info.name.c_str());
    name2index[info.name] = index;
    return index;
}

int Joints::addJointToWorld(Links& links, int linkIndex) {
    LinkInfo& link = links.getLink(linkIndex);
    assert(!link.isInTree());
    JointInfo jinfo("#world_"+link.name, "free");
    jinfo.parentIx = 0; jinfo.childIx = linkIndex;
    const int jx = addJoint(jinfo);
    link.jointsAsChild.push_back(jx); // new link is child
    links.getLink(0).jointsAsParent.push_back(jx); // world is parent
    return jx;
}

void Joints::makeTree(Links& links) {
    for (int level=1; ;++level) { // level of outboard (mobilized) body
        bool anyJointAdded = false;
        int whichLevel = -1;
        for (int jx=0; jx<size(); ++jx) {
            JointInfo& jinfo = jointByIndex[jx];
            if (jinfo.isInTree()) continue; // already done
            LinkInfo&  pinfo = links.getLink(jinfo.parentIx);
            LinkInfo&  cinfo = links.getLink(jinfo.childIx);
            if (pinfo.isInTree()) {
                if (cinfo.isInTree()) continue; // a loop joint
                if (pinfo.level != level-1) continue; // not time yet
                jinfo.isReversed = false;
                jinfo.level=cinfo.level=level;
                cinfo.mobilizer = jx;
                mobilizersByLevel[level].push_back(jx);
                anyJointAdded=true;
            } else if (cinfo.isInTree()) {
                if (cinfo.level != level-1) continue; // not time yet
                jinfo.isReversed = true;
                jinfo.level=pinfo.level=level;
                pinfo.mobilizer = jx;
                mobilizersByLevel[level].push_back(jx);
                anyJointAdded=true;
            }
        }
        if (!anyJointAdded) 
            break;
    }
}

void Joints::breakLoops(Links& links) {
    for (int jx=0; jx<size(); ++jx) {
        JointInfo& jinfo = jointByIndex[jx];
        if (jinfo.isInTree()) continue; // already done
        const int px = jinfo.parentIx, cx = jinfo.childIx;
        assert(links.getLink(px).isInTree() && links.getLink(cx).isInTree());
        // Adds new link -- invalidates old references.
        int shadowIndex = links.splitLink(jinfo.childIx);
        jinfo.masterIx = jinfo.childIx;
        jinfo.childIx = shadowIndex; // replace child with its shadow
        // Add joint to tree.
        LinkInfo&  sinfo = links.getLink(shadowIndex);
        jinfo.level=sinfo.level=links.getLink(px).level+1;
        sinfo.mobilizer = jx;
        sinfo.jointsAsChild.push_back(jx);
        mobilizersByLevel[jinfo.level].push_back(jx);
    }
}

