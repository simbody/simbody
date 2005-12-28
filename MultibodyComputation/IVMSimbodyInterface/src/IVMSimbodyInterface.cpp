/* Copyright (c) 2005-6 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 * Implementation of the IVM Simbody interface.
 */

#include "simbody/internal/SimbodyCommon.h"
#include "simbody/IVMSimbodyInterface.h"

#include "IVMSimbodyInterfaceRep.h"

#include "RigidBodyTree.h"

#include <iostream>
using std::cout;
using std::endl;

#include <vector>

    // IVM SIMBODY INTERFACE //

IVMSimbodyInterface::IVMSimbodyInterface(const Multibody& m) : rep(0) {
    rep = new IVMSimbodyInterfaceRep(m);
    rep->setMyHandle(*this);
}

    // IVM SIMBODY INTERFACE REP //

IVMSimbodyInterfaceRep::IVMSimbodyInterfaceRep(const Multibody& m) 
  : handle(0), mbs(m)
{
    // First find a tree within the multibody system.

    std::vector<const Joint*> joints;
    for (int i=0; i<mbs.getNSubsystems(); ++i)
        if (Joint::isInstanceOf(mbs[i]))
            joints.push_back(&Joint::downcast(mbs[i]));


    size_t nxt=0;
    mbs2tree.push_back(TreeMap(&Body::downcast(mbs["Ground"]),0,0,0));
    while (nxt < mbs2tree.size()) {
        for (size_t i=0; i<joints.size(); ++i) {
            if (!Body::getPlacementBody(joints[i]->getReferenceFrame())
                     .isSameSubsystem(mbs2tree[nxt].getBody())) continue;
            mbs2tree.push_back(TreeMap(&Body::getPlacementBody(joints[i]->getMovingFrame()),
                                       &Body::getPlacementBody(joints[i]->getReferenceFrame()),
                                       joints[i],
                                       mbs2tree[nxt].getLevel() + 1));
        }
        ++nxt;
    }

    mbs.realize(Stage::Startup);


    cout << "**** TREE ****" << endl;
    for (size_t i=0; i<mbs2tree.size(); ++i) {
        cout << mbs2tree[i].getLevel() << ": " << mbs2tree[i].getBody().getFullName() << endl;
        if (!mbs2tree[i].getLevel()) continue;
        cout << " Joint: "  << mbs2tree[i].getJoint().getFullName() 
             << " Parent: " << mbs2tree[i].getParent().getFullName() << endl;

        const Real& mass = mbs2tree[i].getBody().getMass().getValue();
        const Vec3& com  = mbs2tree[i].getBody().getMassCenter().getValue();
        const MatInertia& iner = mbs2tree[i].getBody().getInertia().getValue();
        cout << "mass=" << mass << " com=" << com << " iner=" << iner << endl;

        const Frame& frame = mbs2tree[i].getJoint().getMovingFrame().getValue();
        cout << "frame=" << frame << endl;
    }


}

