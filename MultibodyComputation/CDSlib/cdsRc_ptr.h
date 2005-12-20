
#ifndef __cdsRc_ptr_h__
#define __cdsRc_ptr_h__

/*
	MODULE:			REFCOUNT.H
	DESCRIPTION:	Reference counting smart pointer class
*/

namespace CDS {

template<class T>
class rc_ptr
{
public:
  rc_ptr(T* realPtr = 0);
  rc_ptr(const rc_ptr& rhs);
  ~rc_ptr();
  
  rc_ptr& operator=(const rc_ptr& rhs);
  T* operator->() const;
  T& operator*() const;
  T* ptr() const;

#ifdef TESTING
  static int test();
#endif
private:
  struct CountHolder
  {
    CountHolder(T* realPtr = 0) {pointee=realPtr;count=0;};
    ~CountHolder() {delete pointee;};
    T *pointee;
    int count;
    int operator++() {return ++count;};
    int operator--() {
     if (--count) return count; 
     delete this; 
     return 0;
    }
  };
  CountHolder * counter;
};

template<class T>
rc_ptr<T>::rc_ptr(T* realPtr)
  : counter(new CountHolder(realPtr))
{
 ++(*counter);
}

template<class T>
rc_ptr<T>::rc_ptr(const rc_ptr& rhs)
  : counter(rhs.counter)
{
 ++(*counter);
}

template<class T>
rc_ptr<T>::~rc_ptr()
{
 if (counter) --(*counter);
}

template<class T>
rc_ptr<T>& 
rc_ptr<T>::operator=(const rc_ptr& rhs)
{
 if (counter != rhs.counter)
   {
     if (counter) --(*counter);
     counter = rhs.counter;
     ++(*counter);
   }
 return *this;
}

template<class T>
T* 
rc_ptr<T>::operator->() const
{
 return counter->pointee;
}

template<class T>
T& 
rc_ptr<T>::operator*() const
{
 return *(counter->pointee);
}

template<class T>
T* 
rc_ptr<T>::ptr() const
{
 return counter->pointee;
}

}; /* namespace CDS */

#endif /* __cdsRc_ptr_h__ */

