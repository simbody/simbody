#ifndef SimTK_SIMMATRIX_ELEMENT_FILTER_H_
#define SimTK_SIMMATRIX_ELEMENT_FILTER_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2005-12 Stanford University and the Authors.        *
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

/**@file
 * This declares the ElementFilter class (TODO: this should be an abstract
 * class with concrete derived classes implementing different filtering
 * strategies.)
 *
 * These classes are part of the hidden implementation of MatrixHandle
 * and are not visible to user programs.
 */

#include "SimTKcommon/Scalar.h"
#include "SimTKcommon/SmallMatrix.h"

#include "SimTKcommon/internal/BigMatrix.h"


namespace SimTK {



struct RowCol {
    RowCol(int r, int c) : row(r), col(c) {}

    RowCol transpose() const {return RowCol(col,row);}

    int row, col;

    RowCol& operator+=(const RowCol& rc) {row+=rc.row; col+=rc.col; return *this;}
    RowCol& operator-=(const RowCol& rc) {row-=rc.row; col-=rc.col; return *this;}
};

std::ostream& operator<<(std::ostream&, const RowCol&);

struct NRowCol {
    NRowCol() : nrow(0), ncol(0) {}
    NRowCol(int nr, int nc) : nrow(nr), ncol(nc) {}

    RowCol firstElt() const {return RowCol(0,0);}
    RowCol lastElt()  const {return RowCol(nrow-1,ncol-1);}

    NRowCol transpose() const {return NRowCol(ncol,nrow);}

    int nrow, ncol;
};

inline bool operator==(const NRowCol& left, const NRowCol& right) 
{   return left.nrow == right.nrow && left.ncol == right.ncol; }
inline bool operator!=(const NRowCol& left, const NRowCol& right) 
{   return left.nrow != right.nrow || left.ncol != right.ncol; }
inline bool operator<(const NRowCol& left, const NRowCol& right) 
{   return left.nrow < right.nrow && left.ncol < right.ncol; }
inline bool operator<=(const NRowCol& left, const NRowCol& right) 
{   return left.nrow <= right.nrow && left.ncol <= right.ncol; }
inline bool operator>(const NRowCol& left, const NRowCol& right) 
{   return left.nrow > right.nrow && left.ncol > right.ncol; }
inline bool operator>=(const NRowCol& left, const NRowCol& right) 
{   return left.nrow >= right.nrow && left.ncol >= right.ncol; }

std::ostream& operator<<(std::ostream&, const NRowCol&);

inline RowCol operator+(const RowCol& left, const RowCol& right) {
    return RowCol(left) += right;
}

inline RowCol operator-(const RowCol& left, const RowCol& right) {
    return RowCol(left) -= right;
}

inline bool operator==(const RowCol& left, const RowCol& right) 
{   return left.row == right.row && left.col == right.col; }
inline bool operator!=(const RowCol& left, const RowCol& right) 
{   return left.row != right.row || left.col != right.col; }
inline bool operator<(const RowCol& left, const RowCol& right) 
{   return left.row < right.row && left.col < right.col; }
inline bool operator<=(const RowCol& left, const RowCol& right) 
{   return left.row <= right.row && left.col <= right.col; }
inline bool operator>(const RowCol& left, const RowCol& right) 
{   return left.row > right.row && left.col > right.col; }
inline bool operator>=(const RowCol& left, const RowCol& right) 
{   return left.row >= right.row && left.col >= right.col; }

//------------------------------- EltIndexer -----------------------------------
//
// For a matrix whose elements are regularly spaced with respect to row and
// column indices (i,j), this class captures the spacing *in elements* between
// elements accessed by element index. We represent these like partial
// derivatives dr/di, dr/dj, dc/di, dc/dj where the row and column spacings
// (r,c) are in elements. Note that to be regular each of these values must
// be independent of the current values of i and j.
//------------------------------------------------------------------------------
class EltIndexer {
public:
    EltIndexer() : drdi(1), drdj(0), dcdi(0), dcdj(1) {} // no-op

    EltIndexer(int drdi, int drdj, int dcdi, int dcdj)
    :   drdi(drdi), drdj(drdj), dcdi(dcdi), dcdj(dcdj) {}

    // Return an indexer in which the dependencies on i and j are reversed.
    EltIndexer transpose() const
    {   return EltIndexer(drdj,drdi,dcdj,dcdi); }

    // Apply a new indexer to this one to produce a single indexer with
    // the composite effect (r,c) = post( this(i,j) ).
    EltIndexer postIndexBy(const EltIndexer& post) const {
        return EltIndexer(drdi*post.drdi + drdj*post.dcdi,
                          drdj*post.dcdj + drdi*post.drdj,
                          dcdi*post.drdi + dcdj*post.dcdi,
                          dcdj*post.dcdj + dcdi*post.drdj);
    }

    // Returns true if this indexer maps (i,j)->(i,j) and is thus a no-op.
    bool doesNothing() const
    {   return drdi==1 && drdj==0 && dcdi==0 && dcdj==1; }

    // Returns true if this indexer maps (i,j)->(j,i) and is thus a pure transpose.
    bool isTransposeOnly() const
    {   return drdi==0 && drdj==1 && dcdi==1 && dcdj==0; }

    // Given element index (i,j) relative to this view, return element
    // indices (row,col) from which the data object can retrieve the
    // desired element.        
    int row(int i, int j) const {return drdi*i + drdj*j;}
    int col(int i, int j) const {return dcdi*i + dcdj*j;}
    
private:             
    int drdi, drdj;   // row selection
    int dcdi, dcdj;   // column selection
};

//---------------------------------- EltBlock ----------------------------------
//
// This selects a sub-block of a matrix by giving the element index of its
// upper-left-hand element, and the size of the block.
//------------------------------------------------------------------------------
class EltBlock {
public:
    EltBlock(int nr, int nc)
    :   r0(0), c0(0), nr(nr), nc(nc)
    {   assert(nr>=0 && nc>=0); }
    EltBlock(int r0, int c0, int nr, int nc)
    :   r0(r0), c0(c0), nr(nr), nc(nc) 
    {   assert(r0>=0 && c0>=0 && nr>=0 && nc>=0); }

    // Returns true if, for a matrix of the supplied dimensions, this block
    // would show the whole thing.
    bool isWholeMatrix(int nrow, int ncol) const 
    {   return r0==0 && c0==0 && nr==nrow && nc==ncol; }

    int row0() const {return r0;}
    int col0() const {return c0;}
    int nrow() const {return nr;}
    int ncol() const {return nc;}
    ptrdiff_t nelt() const {return ptrdiff_t(nr)*nc;}

private:
    int r0, c0;
    int nr, nc;
};

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
        Indexer(int r, int c, int drdx=1, int drdy=0, int dcdx=0, int dcdy=1)
          : r0(r), c0(c), drdx(drdx), drdy(drdy), dcdx(dcdx), dcdy(dcdy) 
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
    // and an indexer to use to extract the elements from the original data.
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

    ElementFilter* clone() const {return new ElementFilter(*this);}
        
    int    nrow() const { return nr; }
    int    ncol() const { return nc; } 
    size_t nelt() const { return (size_t)nr*nc; }

    const Indexer& getIndexer() const {return indexer;}

    void setIndexer(const Indexer& newIndexer) {indexer=newIndexer;}
        
    bool isViewWritable() const { return writable; }
    
    int r(int i, int j) const { assert(i<nr&&j<nc); return indexer.row(i,j); }
    int c(int i, int j) const { assert(i<nr&&j<nc); return indexer.col(i,j); }
    void  rc(int i, int j, int& r, int& c) const 
        { assert(i<nr&&j<nc); r=indexer.row(i,j); c=indexer.col(i,j); }   

    RowCol rowcol(int i, int j) const {
        assert(i<nr&&j<nc);
        return RowCol(indexer.row(i,j), indexer.col(i,j));
    }
                       
private:
    bool    writable;           // does this view allow writing?
    int     nr,nc;              // logical shape
    Indexer indexer;            // view->data mapping
};


} // namespace SimTK   


#endif // SimTK_SIMMATRIX_ELEMENT_FILTER_H_
