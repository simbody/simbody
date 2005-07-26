
#ifndef __cdsAuto_ptr_h__
#define __cdsAuto_ptr_h__ 1

#include <sthead.h>

CDS_BEGIN_NAMESPACE

template<class T>
class auto_ptr {
private:
  T *data;

public:
  explicit auto_ptr(T *p =0);  //default constructor
  auto_ptr(auto_ptr&);         //copy constructor
  void operator=(auto_ptr&);   //assignment operator
  template<class S>
  void operator=(auto_ptr<S>&);   //assignment operator
  ~auto_ptr();                 //destructor

  //member functions
  T& operator*() const;
  T* operator->() const;
  T* get() const;
  T* release();
  T* reset(T *p = 0);

#ifdef TESTING
  static int test();
#endif
};

template<class T>
auto_ptr<T>::auto_ptr(T* p)
{
 data = p;
}

template<class T>
auto_ptr<T>::auto_ptr(auto_ptr &a)
{
 data = a.release();
}

template<class T>
void auto_ptr<T>::operator=(auto_ptr &a)
{
 if (&a == this) return;  //check for self-assignment

 reset(a.release());
}

template<class T>
template<class S>
void auto_ptr<T>::operator=(auto_ptr<S> &a)
{
 reset(a.release());
}

template<class T>
auto_ptr<T>::~auto_ptr()
{
 delete data;
}

template<class T>
T& auto_ptr<T>::operator*() const
{
 return *get();
}

template<class T>
T* auto_ptr<T>::get() const
{
 return data;
}

template<class T>
T* auto_ptr<T>::operator->() const
{
 return get();
}

template<class T>
T* auto_ptr<T>::release()
{
 T* op = data;
 data = 0;
 return op;
}

template<class T>
T* auto_ptr<T>::reset(T* p)
{
 delete data;
 data = p;
 return data;
}


CDS_END_NAMESPACE





//  20.4.5  Template class auto_ptr                         [lib.auto.ptr]
//
//1 Template  auto_ptr  holds  onto a pointer obtained via new and deletes
//  that object when it itself is destroyed (such as  when  leaving  block
//  scope _stmt.dcl_).
//  namespace std {
//    template<class X> class auto_ptr {
//    public:
//    // _lib.auto.ptr.cons_ construct/copy/destroy:
//	explicit auto_ptr(X* p =0);
//	auto_ptr(auto_ptr&);
//	void operator=(auto_ptr&);
//     ~auto_ptr();
//    // _lib.auto.ptr.members_ members:
//	X& operator*() const;
//	X* operator->() const;
//	X* get() const;
//	X* release();
//	X* reset(X* p =0);
//    };
//  }
//
//2 The  auto_ptr provides a semantics of strict ownership.  An object may
//  be safely pointed to by only one  auto_ptr,  so  copying  an  auto_ptr
//  copies the pointer and transfers ownership to the destination.
//
//  20.4.5.1  auto_ptr constructors                    [lib.auto.ptr.cons]
//
//  explicit auto_ptr(X* p =0);
//
//  Requires:
//    p points to an object of class X or a class derived from X for which
//    delete p is defined and accessible, or else p is a null pointer.
//  Postcondition:
//    get() == p
//
//  auto_ptr(auto_ptr& a);
//
//  Effects:
//    copies the argument a to *this.
//    Calls a.release().
//  Postcondition:
//    get() == the value returned from a.release().5)
//
//  void operator=(auto_ptr& a);
//
//  Effects:
//    copies the argument a to *this.
//    Calls reset(a.release()).
//  Postcondition:
//    get() == the value returned from a.release().
//
//  ~auto_ptr();
//
//  Effects:
//    delete get()
//
//  20.4.5.2  auto_ptr members                      [lib.auto.ptr.members]
//
//  X& operator*() const;
//
//  Requires:
//    get() != 0
//  Returns:
//    *get()
//
//  X* get() const;
//
//  Returns:
//    The   pointer  p  specified  as  the  argument  to  the  constructor
//    auto_ptr(X* p) or as  the  argument  to  the  most  recent  call  to
//    reset(X* p).
//
//  _________________________
//  5) That is, the value returned by  a.get()  before  clearing  it  with
//  a.release().
//
//  X* operator->() const;
//
//  Returns:
//    get()->m
//
//  X* release();
//
//  Postcondition:
//    get() == 0
//
//  X* reset(X* p =0);
//
//  Requires:
//    p points to an object of class X or a class derived from X for which
//    delete p is defined and accessible, or else p is a null pointer
//  Postcondition:
//    get() == p
//
#endif /* __cdsAuto_ptr_h__ */
