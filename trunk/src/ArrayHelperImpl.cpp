/* -------------------------------------------------------------------------- *
 *                      SimTK Core: SimTKcommon                               *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK Core biosimulation toolkit originating from      *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org.               *
 *                                                                            *
 * Portions copyright (c) 2005-7 Stanford University and the Authors.         *
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

#include "SimTKcommon/internal/common.h"

#include "ArrayHelperImpl.h"

#include <cstddef>
#include <cstring>
#include <algorithm>

// Pre-allocating constructor without initialization.
SimTKimpl::ArrayData::ArrayData(const SimTKimpl::TypeManipulatorT& tmt, int n)
  : tmanip(tmt), nElts(n), nAlloc(n), data(0) 
{   assert(n>=0);
    data = tmanip.createArrayOfT(nAlloc,0);
}

// Pre-allocating and initializing constructor.
// If init is non zero, it either points to a single element or an array; repeat==true
// if it is just one element. Pass init=0 to leave elements default constructed in
// which case repeat is ignored.
SimTKimpl::ArrayData::ArrayData
   (const SimTKimpl::TypeManipulatorT& tmt, int n, const void* init, bool repeat)
  : tmanip(tmt), nElts(n), nAlloc(n), data(0) 
{   assert(n>=0);
    data = tmanip.createArrayOfT(nAlloc,0);
    if (!init) return;
    if (repeat) tmanip.setT(data, init, nElts);
    else tmanip.assignArrayOfT(data, init, nElts);
}

// Copy constructor
SimTKimpl::ArrayData::ArrayData(const ArrayData& ahi)
  : tmanip(ahi.tmanip), nElts(ahi.nElts), nAlloc(ahi.nElts), data(0)
{
	data = tmanip.createArrayOfT(nElts,0);
	tmanip.assignArrayOfT(data, ahi.data, nElts);
}

// Restricted copy constructor
SimTKimpl::ArrayData::ArrayData(const ArrayData& ahi, const ArrayMask& mask)
    : tmanip(ahi.tmanip), nElts(mask.size()), nAlloc(mask.size()), data(0)
{
    assert(mask[0]+mask.size() <= ahi.size());  // note that '=' is allowed if size is 0
    if (nElts) {
        data = tmanip.createArrayOfT(nElts,0);
        tmanip.assignArrayOfT(data, ahi.blobAddr(mask[0]), nElts);
    }
}

// View-mediated copy. Data must be the same type and lengths must match.
void
SimTKimpl::ArrayData::copyInThroughMasks
    (const ArrayMask& myMask, const ArrayData& src, const ArrayMask& srcMask)
{
    assert(myMask.size() == srcMask.size());
    assert(isSameType(src));
    // XXX Should optimize for the case when views are contiguous elements.
    for (int i=0; i < myMask.size(); ++i)
        tmanip.setT(blobAddr(myMask[i]), src.blobAddr(srcMask[i]), 1);     
}

void
SimTKimpl::ArrayData::reverseThroughMask(const ArrayMask& myMask) {
    void* temp = tmanip.createOneT(0);
    for (int i=0; i < myMask.size()/2; ++i) {
        const int i1 = myMask[i];
        const int i2 = myMask[myMask.size()-i-1];
        tmanip.setT(temp, blobAddr(i1), 1);
        tmanip.setT(blobAddr(i1), blobAddr(i2), 1);
        tmanip.setT(blobAddr(i2), temp, 1);
    }
    tmanip.destructOneT(temp);
}


void 
SimTKimpl::ArrayData::push_back(const void* blob)
{
	if (nAlloc == nElts)
		reserve(nElts+ChunkSize);
		
	tmanip.setT(blobAddr(nElts), blob, 1);
	++nElts;
}

void
SimTKimpl::ArrayData::reserve(int n)
{   assert(n>=0);
	if (nAlloc < n) {
		void* newData = tmanip.createArrayOfT(n,0); // construct, but don't initialize
		tmanip.assignArrayOfT(newData, data, nElts);	// copy the old stuff
		tmanip.destructArrayOfT(data);
		data = newData;
		nAlloc = n;	
	}
}

void 
SimTKimpl::ArrayData::resize(int n, const void* initBlob)
{   assert(n>=0);
	reserve(n);
	if (n > nElts && initBlob)
		tmanip.setT(blobAddr(nElts), initBlob, n-nElts); // init new stuff
	nElts = n;	
}

void
SimTKimpl::ArrayData::clear() {
    resize(0);
}

