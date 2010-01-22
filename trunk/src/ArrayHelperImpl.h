#ifndef SimTK_SimTKCOMMON_ARRAYHELPERIMPL_H_
#define SimTK_SimTKCOMMON_ARRAYHELPERIMPL_H_

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
#include "SimTKcommon/internal/Exception.h"
#include "SimTKcommon/internal/OldArray.h"

#include <cstddef>
#include <cassert>

using namespace SimTK;

namespace SimTKimpl {

/**
 * This class maps a logical view to the physical one. For now it is only
 * able to pick out a contiguous segment and record if this is a read-only
 * view; later it can be fancier if needed.
 */
class ArrayMask {
public:
    ArrayMask() : offset(0), length(0), readOnly(false) { }
    ArrayMask(int off, int len, bool ro=false) 
      : offset(off), length(len), readOnly(ro) {
        assert(off >= 0 && len >= 0);
    }
    explicit ArrayMask(int len, bool ro=false) 
      : offset(0), length(len), readOnly(ro) {
        assert(len >= 0);
    }
    
    // Create a subview of an exiting view. Readonly flag is set to more
    // restrictive of view.readOnly and ro.
    ArrayMask(const ArrayMask& mask, int off, int len, bool ro=false)
    {   assert(off >= 0 && len >= 0);
        assert(off+len <= mask.length);  // note "=" is OK with length 0
        offset = mask.offset+off; length = len;
        readOnly = (mask.readOnly || ro);
    }
            
    // default copy, assignment, constructor
    
    bool isReadOnly() const { return readOnly; }   
    int  size() const { return length; }  
      
    // Allow an index one past the end.
    int operator[](int ix) const { assert(ix <= length); return offset+ix; }  

private:
    int  offset;
    int  length;
    bool readOnly;
};

/*
 * This implements something like std:vector<blob> for a particular size blob.
 * We don't know the types of the blobs so we can only return their addresses. 
 */	
class ArrayData {
public:
	static const int ChunkSize = 10;
	
	explicit ArrayData(const TypeManipulatorT& tmt)
		: tmanip(tmt), nElts(0), nAlloc(0), data(0) { }
    ArrayData(const TypeManipulatorT& tmt, int n);        
	ArrayData(const TypeManipulatorT& dt, int n, const void* init, bool repeat);

    // Copy constructors for whole object or just a subset of its elements.
	ArrayData(const ArrayData&);
    ArrayData(const ArrayData&, const ArrayMask&);
    
    // View-mediated copy. Data must be the same type and lengths must match.
    void copyInThroughMasks
        (const ArrayMask& myMask, const ArrayData& src, const ArrayMask& srcMask);
    
	~ArrayData() { tmanip.destructArrayOfT(data); }
    
    bool isSameType(const ArrayData& ad) const 
    { return tmanip == ad.tmanip; }

    void reverseThroughMask(const ArrayMask& myMask);

	int  capacity() const { return nAlloc;   }
	int  size()     const { return nElts;    }
	bool empty()    const { return nElts==0; }
	
	void push_back(const void* blob);
	
	void pop_back()
	{
		assert(nElts);
		--nElts;
	}	

	void reserve(int n);
	void resize(int n, const void* initBlob = 0);
    void clear();

	const void* operator[](int i) const 
		{ assert(i<nElts); return blobAddr(i); }
	void* operator[](int i) 
		{ assert(i<nElts); return blobAddr(i); }
		
private:
	const TypeManipulatorT	tmanip;	// about the elements
	int                     nElts;
	int                     nAlloc;
	void*                   data;

// helpers
private:
	void* blobAddr(int i) { return tmanip.indexT(data,i); }
	const void* blobAddr(int i) const { return tmanip.indexConstT(data,i); }
	
// not allowed
private:
	ArrayData& operator=(const ArrayHelperImpl& ahi) 
		{ assert(false); return *this; }
};


/*
 * This is the class to which a SimTKimpl::ArrayBase<T> holds a handle. It is some kind
 * of restricted mask applied to an ArrayData object. An ArrayHelperImpl object may
 * own the underlying ArrayData, or may be sharing it. We're not reference
 * counting, so terrible things will happen if the ArrayData object is destructed
 * prior to some ArrayHelperImpl object that references it!
 */ 
class ArrayHelperImpl {
public:
    static const int ChunkSize = 10;
    
    explicit ArrayHelperImpl(const TypeManipulatorT& tmt, int n=0)
        : owner(true), arrayData(new ArrayData(tmt,n)), mask(n) { }

    ArrayHelperImpl(const TypeManipulatorT& tmt, int n, const void* init, bool repeat)
        : owner(true), arrayData(new ArrayData(tmt,n,init,repeat)), mask(n) { }
        
    // Slice constructor creates a non-owner, read only and writable versions.
    // Note that the writable one will still produce a read only view if the
    // original was read only.
    ArrayHelperImpl(const ArrayHelperImpl& ahi, int offset, int length)
      : owner(false), arrayData(const_cast<ArrayData*>(&ahi.getArrayData())), 
        mask(ahi.getMask(),offset,length,true) { 
        assert(offset>=0 && length>=0);
    }
    ArrayHelperImpl(ArrayHelperImpl& ahi, int offset, int length)
      : owner(false), arrayData(&ahi.updArrayData()), 
        mask(ahi.getMask(),offset,length,false) {
        assert(offset>=0 && length>=0);
    }
          
    // Copy constructor creates an owner even if the source was not.
    ArrayHelperImpl(const ArrayHelperImpl& ahi)
      : owner(true), arrayData(new ArrayData(ahi.getArrayData(), ahi.getMask())), 
        mask(ahi.size()) { }

    // If this is an owner, assignment behaves like the copy constructor although
    // we insist on matching underlying types.
    // If not an owner, then RHS dimensions must match ours. 
    ArrayHelperImpl& operator=(const ArrayHelperImpl& ahi) 
    { 
        assert(getArrayData().isSameType(ahi.getArrayData()));
        if (owner) {
            delete arrayData;
            arrayData = new ArrayData(ahi.getArrayData());
            resizeOwnerMask();
        } else {
           if (mask.isReadOnly())
                SimTK_THROW1(Exception::OperationNotAllowedOnNonconstReadOnlyView,"operator="); 
            updArrayData().copyInThroughMasks(mask,ahi.getArrayData(),ahi.getMask());
        }
        return *this; 
    }     
            
    ~ArrayHelperImpl() { if (owner) delete arrayData; }

    int  size()     const {return mask.size();}
    int  capacity() const {return size();     }
    bool empty()    const {return size()==0;  }
    
    void push_back(const void* blob)
    {   if (!owner) SimTK_THROW1(Exception::OperationNotAllowedOnView,"push_back"); 
        updArrayData().push_back(blob); resizeOwnerMask(); }
    
    void pop_back()
    {   if (!owner) SimTK_THROW1(Exception::OperationNotAllowedOnView,"pop_back"); 
        updArrayData().pop_back(); resizeOwnerMask(); }

    void reserve(int n)    
    {   assert(n>=0);
        if (n==size()) return;  // we'll allow this even on a view
        if (!owner) SimTK_THROW1(Exception::OperationNotAllowedOnView,"reserve"); 
        updArrayData().reserve(n);
    }
         
    void resize(int n, const void* initBlob = 0)
    {   assert(n>=0);
        if (n==size()) return;  // we'll allow this even on a view
        if (!owner) SimTK_THROW1(Exception::OperationNotAllowedOnView,"resize"); 
        updArrayData().resize(n,initBlob); resizeOwnerMask();
    }
        
    void clear() 
    {   if (0==size()) return;  // we'll allow this even on a view
        if (!owner) SimTK_THROW1(Exception::OperationNotAllowedOnView,"clear"); 
        updArrayData().clear(); resizeOwnerMask();
    }

    void reverse() {
        updArrayData().reverseThroughMask(getMask());
    }

    const void* operator[](int i) const 
      { assert(0<=i && i<size()); return (*arrayData)[mask[i]]; }
    void* operator[](int i) 
    {   assert(0<=i && i<size());
        if (mask.isReadOnly())
            SimTK_THROW1(Exception::OperationNotAllowedOnNonconstReadOnlyView,"operator[]"); 
        return (*arrayData)[mask[i]]; }
        
private:
    bool        owner;      // did we allocate this?
    ArrayData*  arrayData;
    
    // If we are the owner this *must* refer to the entire array, meaning
    // offset==0 and length=arrayData.size().
    ArrayMask   mask;
    
    ArrayData&       updArrayData()       {assert(arrayData); return *arrayData;}
    const ArrayData& getArrayData() const {assert(arrayData); return *arrayData;}
    const ArrayMask& getMask()      const {return mask;}
    void resizeOwnerMask() {assert(owner); mask=ArrayMask(getArrayData().size());}
};

} //namespace SimTKimpl

#endif // SimTK_SimTKCOMMON_ARRAYHELPERIMPL_H_
