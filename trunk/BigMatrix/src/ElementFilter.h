#ifndef SimTK_SIMMATRIX_ELEMENT_FILTER_H_
#define SimTK_SIMMATRIX_ELEMENT_FILTER_H_

/* Portions copyright (c) 2005-7 Stanford University and Michael Sherman.
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
 * IN NO EVENT SHALL THE AUTHORS, CONTRIBUTORS, OR COPYRIGHT HOLDERS BE
 * LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION
 * OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
 * WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */

/**@file
 * This declares the ElementFilter class (TODO: this should be an abstract
 * class with concrete derived classes implementing different filtering
 * strategies.)
 *
 * These classes are part of the hidden implementation of MatrixHandle
 * and are not visible to user programs.
 */

#include "SimTKcommon/internal/common.h"
#include "SimTKcommon/internal/Scalar.h"
#include "SimTKcommon/internal/SmallMatrix.h"
#include "SimTKcommon/internal/BigMatrix.h"


namespace SimTK {


/**
 * Describe a subset of the elements of a DataDescriptor object referenced by
 * the same MatrixHelper that references this ElementFilter object. Note that
 * this is done entirely in terms of logical "elements" composed of one or more
 * scalars -- this class is NOT templatized by element type.
 *
 * The DataDescriptor object contains only scalars. It is guaranteed
 * that the storage for an element is composed of consecutive scalars, but 
 * otherwise we don't know anything about them here. Note that different views of
 * the same data can claim elements of different sizes as long as the underlying
 * scalar types match. (That would be used, for example, to select the real or
 * imaginary submatrix of a complex matrix.)
 * 
 * A MatrixBase object does not have to contain an ElementFilter object if it
 * permits unfettered access to all the data elements. In that case all MatrixBase 
 * operations are forwarded directly to the DataDescriptor along with the
 * element size we think we're looking at; otherwise the ElementFilter
 * object must serve as an intermediary, slowing things down. 
 * 
 * This object expects to be intimately coupled with a DataDescriptor object
 * via the MatrixHelper. Don't move it around separately!
 */
class ElementFilter {
public:
    /**
     * This class handles the mundane details of mapping from logical, element-oriented
     * view indices to physical but still element-oriented indices of the 
     * DataDescriptor object, which handles the final mapping of indices to 
     * memory address.
     */        
    class Indexer {
    public:
        // All arguments refer to *elements*, not *scalars*.
        Indexer(int r, int c, int drdx_=1, int drdy_=0, int dcdx_=0, int dcdy_=1)
          : r0(r), c0(c), drdx(drdx_), drdy(drdy_), dcdx(dcdx_), dcdy(dcdy_) 
        { 
        }
 
        // Compose a physical indexer (old) with a relative one (that is, an
        // indexer relative to this view) to produce a new physical one 
        // suitable for use on the original data.          
        Indexer(const Indexer& old, const Indexer& ix)
          : r0(old.row(ix.r0,ix.c0)), c0(old.col(ix.r0,ix.c0)),
            drdx(old.drdx*ix.drdx + old.drdy*ix.dcdx), 
            drdy(old.drdy*ix.dcdy + old.drdx*ix.drdy),
            dcdx(old.dcdx*ix.drdx + old.dcdy*ix.dcdx), 
            dcdy(old.dcdy*ix.dcdy + old.dcdx*ix.drdy) 
        { 
        }

        // Given element index (x,y) relative to this view, return element
        // indices (row,col) from which the data object can retrieve the
        // desired element.        
        int row(int x, int y) const { return r0 + drdx*x + drdy*y; }
        int col(int x, int y) const { return c0 + dcdx*x + dcdy*y; }

        // Make an indexer just like this one but with the roles of x and y reversed.
        // Still has same (0,0) element.
        Indexer transpose() {
            return Indexer(r0,c0,drdy,drdx,dcdy,dcdx);
        }
        
    private:             
        int r0,c0;        // indices of (0,0) element
        int drdx, drdy;   // row selection
        int dcdx, dcdy;   // column selection
        
        // no default construction
        Indexer();
    }; 

    // Like copy constructor but can *reduce* writability. Can't create writability
    // where none was permitted before, however.    
    ElementFilter(ElementFilter& v, bool wrt)
      : writable(v.writable && wrt), nr(v.nr), nc(v.nc), indexer(v.indexer) 
    { 
    }
        
    // Copy constructor takes a const ElementFilter and loses writability.
    ElementFilter(const ElementFilter& v)
      : writable(false), nr(v.nr), nc(v.nc), indexer(v.indexer) 
    { 
    }

    // This is the basic constructor, giving the logical shape, writability
    // and an indexer to use to extract the elements from the orginal data.
    ElementFilter(bool wrt, int m, int n, const Indexer& ix)
      : writable(wrt), nr(m), nc(n), indexer(ix) 
    { 
    }
        
    // Offset constructor -- combine an old one and new instructions to get
    // another ElementFilter suitable for use with the original data.
    ElementFilter(const ElementFilter& old, bool wrt, int m, int n,
                  const Indexer& ix)
      : writable(wrt), nr(m), nc(n), indexer(old.indexer, ix) 
    { 
    } 
        
    int  nrow() const { return nr; }
    int  ncol() const { return nc; } 
    long size() const { return nr*nc; }    
        
    bool isViewWritable() const { return writable; }
    
    int r(int i, int j) const { assert(i<nr&&j<nc); return indexer.row(i,j); }
    int c(int i, int j) const { assert(i<nr&&j<nc); return indexer.col(i,j); }
    void  rc(int i, int j, int& r, int& c) const 
        { assert(i<nr&&j<nc); r=indexer.row(i,j); c=indexer.col(i,j); }   
                       
private:
    bool    writable;           // does this view allow writing?
    int     nr,nc;              // logical shape
    Indexer indexer;            // view->data mapping
};


} // namespace SimTK   


#endif // SimTK_SIMMATRIX_ELEMENT_FILTER_H_
