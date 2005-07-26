
#ifndef __accessor_hh__
#define __accessor_hh__

#include <modified.h>

template<class T> 
class ModAccessor {
  T data;
  Modified& modified;
  //  Accessor(const Accessor<T>&);
public:
  ModAccessor(Modified& modified) : modified(modified) { }
  ModAccessor(const T&        i,
		    Modified& modified) : data(i), modified(modified)  { }
  ~ModAccessor() { }

  void operator=(const ModAccessor<T>& a)
    { if (&a == this) return; data = a.data; }

  void operator=(const T& x) { set(x); }

  T& raw() { return data; }

  const T& operator()() const   {return data;}
  T& get() { return data; }
  void set(const T& d) {data = d; modified.set();}
};

template<class T> 
class Accessor {
  T data;
public:
  Accessor(){ }
  Accessor(const T&        i) :  data(i) { }
  ~Accessor() { }

  void operator=(const Accessor<T>& a)
    { if (&a == this) return; data = a.data; }

  T& raw() { return data; }

  const T& operator()() const   {return data;}
  T& get() { return data; }
  void set(const T& d) {data = d;}
};

#define ACCESSOR(name,cap,type) \
Accessor<type> name; \
void set ## cap(const type& val) { name.set(val); }

#define MODACCESSOR(name,cap,type) \
ModAccessor<type> name; \
void set ## cap(const type& val) { name.set(val); }



// this version uses pointer so don't need to include class def
// in header files
template<class T>
class PointerAccessor {
  T* data;
  PointerAccessor(const PointerAccessor<T>&);
public:
  PointerAccessor() : data(new T) { }
  PointerAccessor(const T& i) : data(new T) { *data = i;}
  ~PointerAccessor() { delete data; }
  void operator=(const PointerAccessor<T>& a)
    { if (&a == this) return; *data = *a.data; }
  const T& operator()() const   {return *data;}
  const T& get() { return data; }
  void set(const T& d) {*data = d;}
};

#endif /* __accessor_hh__ */
