
#include "circDLinkedList.h"

#include <cdsExcept.h>

namespace CDS {

template<class T>
CircDLinkedList<T>::CircDLinkedList() : _start(0) {}

template<class T>
CircDLinkedList<T>::~CircDLinkedList()
{
 while ( start().pos() )
   remove(start());
} /* destructor */

template<class T>
void
CircDLinkedList<T>::remove(const Iterator<CircDLinkedList<T> >& i)
{
 Letter* curr = i.pos();
 Letter* prev = curr->prev;
 Letter* next = curr->next;

 prev->next = next;
 next->prev = prev;

 if ( _start->next == _start ) {
   _start = 0;
 } else {
   if (_start == curr) _start = next;
 }
 
 delete curr;
 
} /* remove */

template<class T>
void
CircDLinkedList<T>::add(const T& elem)
{
 if (!_start) {
   _start = new Letter(0,0,elem);
   _start->prev = _start->next = _start;
 } else {
   new Letter(_start->prev,_start,elem);
 }
}  /* add */ 

template<class T>
void
CircDLinkedList<T>::setStart(const Iterator<CircDLinkedList<T> >& i)
{
 if ( &i.container() != this )
   throw CDS::exception("CircDLinkedList::setStart: container mismatch");

 _start = i.pos();
} /* setStart */


template<class CONTAINER>
Iterator<CONTAINER>::Iterator(CONTAINER&                    container,
			      typename CONTAINER::IndexType index) :
  _container(container), _pos(index?index:container.startPos())
{}

template<class CONTAINER>
Iterator<CONTAINER>::Iterator(const Iterator<CONTAINER>& i) :
  _container(i.container()), _pos(i.pos())
{}

template<class CONTAINER>
Iterator<CONTAINER>&
Iterator<CONTAINER>::operator=(const Iterator<CONTAINER>& i) 
{
 _container = i.container();
 _pos = i.pos();
 return *this;
} /* operator= */

template<class CONTAINER>
typename CONTAINER::ElemType&
Iterator<CONTAINER>::current() 
{
 return container().get(pos());
} /* current */

template<class CONTAINER>
Iterator<CONTAINER>& 
Iterator<CONTAINER>::operator++() 
{
 _pos = container().next(pos());
 return *this;
} /* current */

template<class CONTAINER>
Iterator<CONTAINER>& 
Iterator<CONTAINER>::operator--() 
{
 _pos = container().prev(pos());
 return *this;
} /* current */

template<class CONTAINER>
bool
Iterator<CONTAINER>::operator==(const Iterator<CONTAINER>& i) 
{
 return (pos() == i.pos());
} /* operator== */


}



#ifdef TESTING

#include <stdlib.h>
#include <string.h>
#include <cdsSStream.h>

namespace CDS {

 struct TC {
   TC(int a) : a(a) {}
   int a;
 };

template<class T>
int
CircDLinkedList<T>::test()
{
 int exit=0;
 cout << "testing CircDLinkedList...";

 {
   CircDLinkedList<int> l;
   // CircDLinkedList<int> l2(1,2);
   // CircDLinkedList<int> l3 = l2;

   l.add(1);
   l.add(2);
   l.add(3);
   
   typedef CircDLinkedList<int> IList;
   typedef CDS::Iterator<IList> IIterator;
   
   IIterator i(l);
   
   if (*i != 1) {
     cerr << "iterator value failed: " << *i << " != 1\n";
     exit=1;
   }

   ++i;
   if (*i != 2) {
     cerr << "iterator value failed: " << *i << " != 2\n";
     exit=1;
   }

   if (*l.end().prev() != 3) {
     cerr << "list end failed: " << *l.end().prev() << " != 3\n";
     exit=1;
   }
   
   ++i;
   if (l.end().prev() != i) {
     cerr << "iterator comparison: " << l.endPos() << " != " << i.pos();
     exit=1;
   }
   
   l.remove( l.end().prev() );
   if (*l.end().prev() != 2) {
     cerr << "list remove failed: " << *i << " != 2\n";
     exit=1;
   }
 }

 typedef CDS::Iterator<CDS::CircDLinkedList<TC> > TCIterator;
 CircDLinkedList<TC> l;

 for (int i=0 ; i<10 ; i++)
   l.add(TC(i));
 
 TCIterator i(l);
 
 ++i;++i;
 if (i->a != 2) {
   cerr << "iterator -> failed: " << i->a << " != 2\n";
   exit=1;
 }

 l.setStart(i);
 if ( l.start()->a != 2) {
   cerr << "list.start failed: " << l.start()->a << " != 2\n";
   exit=1;
 }


 TCIterator i2(i);

 ++i2;
 if (i2->a != 3) {
   cerr << "iterator copy constructor failed: " << i2->a << " != 3\n";
   exit=1;
 }

 i2 = i;
 if (i2->a != 2) {
   cerr << "iterator assignment operator failed: " << i2->a << " != 2\n";
   exit=1;
 }

 // TCIterator i2(i);

//
//#ifdef DEBUG_ALLOC
// { // two stress tests
//   for (int i=0 ; i<8 ; i++) {
//     // 8 x 1e7 x sizeof(int)
//     CDSList<int> l;
//     for (int j=0 ; j<10000000 ; j++ )
//       l.append(j);
//     cout << i << ' ';
//     cout.flush();
//   }
//   CDSList<int> l;
//   for (l_int i=0 ; i<80 ; i++) {
//     l.resize(100000000);
//     l[99999999] = i;
//     l[0] = 0;
//   }
// }
//#endif /* DEBUG_ALLOC */
//
//
 cout << (exit?"failed":"ok") << endl;
 return exit;
} /* test */

}

#endif /* TESTING */
