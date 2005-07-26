//auto_arr.h
//automatically deleted array - based on auto_ptr
//CDS 9/30/96

#ifndef __cdsAuto_arr_h__
#define __cdsAuto_arr_h__ 1

template<class T>
class auto_arr {
private:
  T *data;

public:
  explicit auto_arr(T *p =0);  //default constructor
  auto_arr(auto_arr<T>&);         //copy constructor
  void operator=(auto_arr<T>&);   //assignment operator
  ~auto_arr();                 //destructor

  //member functions
  T& operator*() const;
  operator T*() {return data;};
//  T* operator->() const;
  T* get() const;
  T* release();
  T* reset(T *p = 0);
};

template<class T>
auto_arr<T>::auto_arr(T* p)
{
 data = p;
}

template<class T>
auto_arr<T>::auto_arr(auto_arr<T> &a)
{
 data = a.release();
}

template<class T>
void auto_arr<T>::operator=(auto_arr<T> &a)
{
 if (&a == this) return;  //check for self-assignment

 reset(a.release());
}

template<class T>
auto_arr<T>::~auto_arr()
{
 delete [] data;
}

template<class T>
T& auto_arr<T>::operator*() const
{
 return *get();
}

template<class T>
T* auto_arr<T>::get() const
{
 return data;
}

//template<class T>
//T* auto_arr<T>::operator->() const
//{
// return get();
//}

template<class T>
T* auto_arr<T>::release()
{
 T* op = data;
 data = 0;
 return op;
}

template<class T>
T* auto_arr<T>::reset(T* p)
{
 data = p;
 return data;
}

int auto_arr_test();

#endif /*__cdsAuto_arr_h__*/

