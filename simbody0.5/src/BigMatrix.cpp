/* Copyright (c) 2005 Stanford University and Michael Sherman.
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

#include "SimTKcommon.h"
#include "simmatrix/internal/common.h"
#include "simmatrix/internal/Scalar.h"
#include "simmatrix/internal/SmallMatrix.h"
#include "simmatrix/internal/BigMatrix.h"

#include "BigMatrixImpl.h"

namespace SimTK {
    
using namespace simtkimpl;
    
/*===========================*
 * MatrixHelper definitions. *
 *===========================*/

template <class S> 
MatrixHelper<S>::MatrixHelper(unsigned int esz, int m, int n, bool lockNrow, bool lockNcol)
  : eltSize(esz), flags(0), view(0), data(new MatrixDataImpl<S>(esz,m,n,lockNrow,lockNcol))
{ }

// These are non-resizeable owners (of the data descriptor) but share pre-existing raw data.
// We have read-only and writable varieties.
template <class S> 
MatrixHelper<S>::MatrixHelper(unsigned int esz, int m, int n, int leadingDim, const S* s)
  : eltSize(esz), flags(0), view(0), data(new MatrixDataImpl<S>(esz,m,n,leadingDim,s)) { }
template <class S> 
MatrixHelper<S>::MatrixHelper(unsigned int esz, int m, int n, int leadingDim, S* s)
  : eltSize(esz), flags(0), view(0), data(new MatrixDataImpl<S>(esz,m,n,leadingDim,s)) { }

// copy constructor is suppressed; these are closest things but require
// specification of deep or shallow copy

// Create a read-only view of existing data.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixHelper& h, const ShallowCopy&)
  : eltSize(h.eltSize), flags(h.flags), view(0), data(0) 
{
    assert(h.data);
    view = h.view ? new MatrixViewImpl(*h.view) // loses writability
                  : new MatrixViewImpl(false,h.data->nrowElt(eltSize),h.data->ncolElt(eltSize),
                                       MatrixViewImpl::Indexer(eltSize,0,0));
    data = h.data;
}

// Create a (possibly) writable view of existing data. 
template <class S>
MatrixHelper<S>::MatrixHelper(MatrixHelper& h, const ShallowCopy&) 
  : eltSize(h.eltSize), flags(h.flags), view(0), data(0) 
{
    assert(h.data);
    view = h.view ? new MatrixViewImpl(*h.view, true)    // keep writability if possible
                  : new MatrixViewImpl(true,h.data->nrowElt(eltSize),h.data->ncolElt(eltSize),
                                       MatrixViewImpl::Indexer(eltSize,0,0));
    data = h.data;
}

// Deep copy. "This" will have no view, and the result is always writable, packed storage.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixHelper& h, const DeepCopy&)
  : eltSize(h.eltSize), flags(h.flags), view(0), data(0) 
{
    if (h.data)
        data = h.view ? new MatrixDataImpl<S>(eltSize,*h.view,*h.data)
                      : new MatrixDataImpl<S>(*h.data);
}

template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixHelper<typename CNT<S>::TNeg>& h, const DeepCopy&)
  : eltSize(h.eltSize), flags(h.flags), view(0), data(0) 
{
    if (h.data) {
        const MatrixDataImpl<S>* negSrc = reinterpret_cast<const MatrixDataImpl<S>*>(h.data);
        data = h.view ? new MatrixDataImpl<S>(eltSize,*h.view,*negSrc)
                      : new MatrixDataImpl<S>(*negSrc);
        scaleBy(typename CNT<S>::StdNumber(-1));
    }
}

// assignment depends on whether lhs (this) is a view
template <class S> MatrixHelper<S>&
MatrixHelper<S>::operator=(const MatrixHelper& h) {
    if (&h == this) return *this;

    assert(getElementSize()==h.getElementSize());
    assert(view==0 || (nrow()==h.nrow() && ncol()==h.ncol()));
    if (view == 0) {
        delete data; data = 0;
        if (h.data)
            data = h.view ? new MatrixDataImpl<S>(getElementSize(),*h.view,*h.data)
                          : new MatrixDataImpl<S>(*h.data);
    } else {
        // "this" has a view. TODO: this is really bad!
        for (int j=0; j<ncol(); ++j)
            for (int i=0; i<nrow(); ++i) 
                MatrixDataImpl<S>::copyElement(getElementSize(),
                                               updElt(i,j), h.getElt(i,j));
    }

    return *this;
}

template <class S> MatrixHelper<S>&
MatrixHelper<S>::operator=(const MatrixHelper<typename CNT<S>::TNeg>& h) {
    this->operator=(reinterpret_cast<const MatrixHelper<S>&>(h)); // negative of result
    this->scaleBy(typename CNT<S>::StdNumber(-1));
    return *this;
}

// construct read-only view
template <class S> 
MatrixHelper<S>::MatrixHelper(const MatrixHelper& h, int i, int j, int m, int n)
  : eltSize(h.eltSize), flags(h.flags), view(0), data(0) 
{
    const MatrixViewImpl::Indexer ix(i,j);
    if (h.view) view = new MatrixViewImpl(*h.view,false,m,n,ix);
    else view = new MatrixViewImpl(false,m,n,ix);
    data = h.data;
}

// construct writable view
template <class S> 
MatrixHelper<S>::MatrixHelper(MatrixHelper& h, int i, int j, int m, int n) 
  : eltSize(h.eltSize), flags(h.flags), view(0), data(0) 
{
    const MatrixViewImpl::Indexer ix(i,j);
    if (h.view) view = new MatrixViewImpl(*h.view,true,m,n,ix);
    else view = new MatrixViewImpl(true,m,n,ix);
    data = h.data;
} 

// Construct read only transposed view of passed-in helper.
template <class S>
MatrixHelper<S>::MatrixHelper(const MatrixHelper<typename CNT<S>::THerm>& h,
                              const TransposeView&)
  : eltSize(h.eltSize), flags(h.flags), view(0), data(0)
{
    const MatrixViewImpl::Indexer ix(MatrixViewImpl::Indexer(0,0).transpose());
    if (h.view) view = new MatrixViewImpl(*h.view,false,h.ncol(),h.nrow(),ix);
    else view = new MatrixViewImpl(false,h.ncol(),h.nrow(),ix);
    data = reinterpret_cast<MatrixDataImpl<S>*>(h.data);
}

// Construct writable transposed view of passed-in helper.
template <class S>
MatrixHelper<S>::MatrixHelper(MatrixHelper<typename CNT<S>::THerm>& h,
                              const TransposeView&)
  : eltSize(h.eltSize), flags(h.flags), view(0), data(0)
{    
    const MatrixViewImpl::Indexer ix(MatrixViewImpl::Indexer(0,0).transpose());
    if (h.view) view = new MatrixViewImpl(*h.view,true,h.ncol(),h.nrow(),ix);
    else view = new MatrixViewImpl(true,h.ncol(),h.nrow(),ix);
    data = reinterpret_cast<MatrixDataImpl<S>*>(h.data);
}

template <class S> void
MatrixHelper<S>::scaleBy(const typename CNT<S>::StdNumber& s) {
    if (!view) { data->scaleBy(s); return; }
    // XXX -- really, really bad! Optimize for contiguous data!
    const int sz = getElementSize();
    for (int j=0; j<view->ncol(); ++j)
        for (int i=0; i<view->nrow(); ++i) 
            MatrixDataImpl<S>::scaleElement(sz,updElt(i,j),s);
}  
     
template <class S> void
MatrixHelper<S>::addIn(const MatrixHelper& h) {
    assert(nrow()==h.nrow() && ncol()==h.ncol());
    assert(getElementSize()==h.getElementSize());
    // XXX -- really, really bad! Optimize for contiguous data, missing views, etc.!
    const int sz = getElementSize();
    for (int j=0; j<ncol(); ++j)
        for (int i=0; i<nrow(); ++i)
            MatrixDataImpl<S>::addToElement(sz,updElt(i,j),h.getElt(i,j));
} 

template <class S> void
MatrixHelper<S>::addIn(const MatrixHelper<typename CNT<S>::TNeg>& h) {
    subIn(reinterpret_cast<const MatrixHelper<S>&>(h));
}
     
template <class S> void
MatrixHelper<S>::subIn(const MatrixHelper& h) {
    assert(nrow()==h.nrow() && ncol()==h.ncol());
    assert(getElementSize()==h.getElementSize());
    // XXX -- really, really bad! Optimize for contiguous data, missing views, etc.!
    const int sz = getElementSize();
    for (int j=0; j<ncol(); ++j)
        for (int i=0; i<nrow(); ++i)
            MatrixDataImpl<S>::subFromElement(sz,updElt(i,j),h.getElt(i,j));
}  

template <class S> void
MatrixHelper<S>::subIn(const MatrixHelper<typename CNT<S>::TNeg>& h) {
    addIn(reinterpret_cast<const MatrixHelper<S>&>(h));
}

template <class S> void
MatrixHelper<S>::fillWith(const S* eltp) {
    // XXX -- really, really bad! Optimize for contiguous data, missing views, etc.!
    const int sz = getElementSize();
    for (int j=0; j<ncol(); ++j)
        for (int i=0; i<nrow(); ++i)
            MatrixDataImpl<S>::copyElement(sz,updElt(i,j),eltp);
} 

template <class S> void 
MatrixHelper<S>::copyInByRows(const S* elts) {
    // XXX -- really, really bad! Optimize for contiguous data, missing views, etc.!
    const int sz = getElementSize();
    for (int j=0; j<ncol(); ++j)
        for (int i=0; i<nrow(); ++i)
            MatrixDataImpl<S>::copyElement(sz,updElt(i,j),elts+i*ncol()+j);
}

template <class S> 
MatrixHelper<S>::~MatrixHelper()
{
    if (view==0) delete data;
    delete view;
}

template <class S> void
MatrixHelper<S>::invertInPlace() {
    assert(view==0 && eltSize == 1);
    if (data)
        data->invertInPlace();
}

template <class S>
template <class SA, class SB>
void MatrixHelper<S>::matmul(
    const StdNumber& beta,   // applied to 'this'
    const StdNumber& alpha, const MatrixHelper<SA>& A, const MatrixHelper<SB>& B)
{
    assert(false); //TODO
}

template <class S> const S* 
MatrixHelper<S>::getElt(int i, int j) const {
    assert(0<=i && i<nrow() && 0<=j && j<ncol());
    if (view) return data->getElt(getElementSize(),view->r(i,j), view->c(i,j));
    else return data->getElt(getElementSize(),i,j);
}       
    
template <class S> S* 
MatrixHelper<S>::updElt(int i, int j) {
    assert(0<=i && i<nrow() && 0<=j && j<ncol());
    if (view) {
        assert(view->isViewWritable());
        return data->updElt(getElementSize(), view->r(i,j), view->c(i,j));
    }
    return data->updElt(getElementSize(),i,j);
}

template <class S> const S& 
MatrixHelper<S>::getScalar(const int i, const int j) const {
    assert(0<=i && i<nrow() && 0<=j && j<ncol());
    if (view) return data->getScalar(view->r(i,j), view->c(i,j));
    else return data->getScalar(i,j);
} 

template <class S> S& 
MatrixHelper<S>::updScalar(const int i, const int j) {
    assert(0<=i && i<nrow() && 0<=j && j<ncol());
    if (view) {
        assert(view->isViewWritable());
        return data->updScalar(view->r(i,j), view->c(i,j));
    }
    return data->updScalar(i,j);
}

template <class S> int 
MatrixHelper<S>::nrow() const {return view?view->nrow():(data?data->nrowElt(eltSize):0);}

template <class S> int 
MatrixHelper<S>::ncol() const {return view?view->ncol():(data?data->ncolElt(eltSize):0);} 

template <class S> ptrdiff_t 
MatrixHelper<S>::size() const {return view?view->size():(data?data->sizeElt(eltSize):0);}

template <class S> bool
MatrixHelper<S>::isDataOwned() const {return view==0;} 

template <class S> void 
MatrixHelper<S>::resize(int m, int n) {
    if (nrow()==m && ncol()==n) return;
    assert(isDataOwned());
    if (data) data->resize(getElementSize(),m,n);
    else data=new MatrixDataImpl<S>(getElementSize(),m,n,false,false);          
} 

template <class S> void
MatrixHelper<S>::sum(S* const answer) const {
    const int sz = getElementSize();
    if (sz==1) *answer = scalarSum();
    else {
        S* csum = new S[sz];
        MatrixDataImpl<S>::fillElement(sz, answer, typename CNT<S>::StdNumber(0));
        for (int j=0; j<ncol(); ++j) {
            colSum(j, csum);
            MatrixDataImpl<S>::addToElement(sz, answer, csum); // answer+=csum
        }
        delete csum;
    }
}        

template <class S> void
MatrixHelper<S>::colSum(const int j, S* const answer) const {
    const int sz = getElementSize();
    if (sz==1) *answer = scalarColSum(j);
    else {
        MatrixDataImpl<S>::fillElement(sz, answer, typename CNT<S>::StdNumber(0));
        for (int i=0; i<nrow(); ++i)
            MatrixDataImpl<S>::addToElement(sz, answer, getElt(i,j));
    }
} 

template <class S> void
MatrixHelper<S>::rowSum(const int i, S* const answer) const {
    const int sz = getElementSize();
    if (sz==1) *answer = scalarRowSum(i);
    else {
        MatrixDataImpl<S>::fillElement(sz, answer, typename CNT<S>::StdNumber(0));
        for (int j=0; j<ncol(); ++j)
            MatrixDataImpl<S>::addToElement(sz, answer, getElt(i,j));
    }
} 

template <class S> void
MatrixHelper<S>::fillWithScalar(const typename CNT<S>::StdNumber& s) {
    // XXX -- really, really bad! Optimize for contiguous data, missing views, etc.!
    const int sz = getElementSize();
    for (int j=0; j<ncol(); ++j)
        for (int i=0; i<nrow(); ++i)
            MatrixDataImpl<S>::fillElement(sz,updElt(i,j),s);
} 

template <class S> S
MatrixHelper<S>::scalarColSum(const int j) const {
    assert(getElementSize()==1);
    S sum = S(0);
    for (int i=0; i<nrow(); ++i) sum += getScalar(i,j);
    return sum;
} 

template <class S> S
MatrixHelper<S>::scalarRowSum(const int i) const {
    assert(getElementSize()==1);
    S sum = S(0);
    for (int j=0; j<ncol(); ++j) sum += getScalar(i,j);
    return sum;
} 

template <class S> S
MatrixHelper<S>::scalarSum() const {
    assert(getElementSize()==1);
    S sum = S(0);
    for (int j=0; j<ncol(); ++j) sum += scalarColSum(j);
    return sum;
}     

// Instantiations for each of the 18 Scalar types.
template class MatrixHelper< float >;
template class MatrixHelper< double >;
template class MatrixHelper< long double >;
template class MatrixHelper< std::complex<float> >;
template class MatrixHelper< std::complex<double> >;
template class MatrixHelper< std::complex<long double> >;
template class MatrixHelper< conjugate<float> >;
template class MatrixHelper< conjugate<double> >;
template class MatrixHelper< conjugate<long double> >;

template class MatrixHelper<negator< float > >;
template class MatrixHelper<negator< double > >;
template class MatrixHelper<negator< long double > >;
template class MatrixHelper<negator< std::complex<float> > >;
template class MatrixHelper<negator< std::complex<double> > >;
template class MatrixHelper<negator< std::complex<long double> > >;
template class MatrixHelper<negator< conjugate<float> > >;
template class MatrixHelper<negator< conjugate<double> > >;
template class MatrixHelper<negator< conjugate<long double> > >;
   
} // namespace SimTK   
