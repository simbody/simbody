#ifndef __cdsListAutoPtr_hh__
#define __cdsListAutoPtr_hh__

//
// a list of auto_ptrs. The append method is specialized. This is necessary
// because auto_ptr::operator= modifies thr rhs.
//

#include <cdsList.h>
#include <cdsAuto_ptr.h>

#include <cassert>

template<class T>
class CDSListAutoPtr : public CDSList< CDS::auto_ptr<T> > {
    typedef CDSList< CDS::auto_ptr<T> > Base;
public:
    //specialize
    void append(const CDS::auto_ptr<T>& m) {
        assert("CDSListAutoPtr: this append should not be used.");
    }

    T* append(T* p) { 
        resize(Base::size()+1);
        CDS::auto_ptr<T> ap(p);
        (*this)[Base::size()-1] = ap;
        return p;
    }
}; 

#endif /* __cdsListAutoPtr_hh__ */
