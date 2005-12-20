
#ifndef __circDLinkedList_hh__
#define __circDLinkedList_hh__

//
// Circular doubly-linked list data structure.
// should be traversed using the Iterator class.
// NOTE: the list's end() and start() methods return the same iterator!



namespace CDS {

template<class CONTAINER> class Iterator;

template<class T>
class CircDLinkedList {
public:

  struct Letter {
    Letter* prev;
    Letter* next;
    T data;

    Letter(Letter*  prev,
	   Letter*  next,
	   const T& elem) : prev(prev), next(next), data(elem) {
    if (prev) prev->next=this;
    if (next) next->prev=this;}
  };

  typedef T       ElemType;
  typedef Letter* IndexType;


  CircDLinkedList();
  CircDLinkedList(T&);
  CircDLinkedList(const CircDLinkedList&);
  ~CircDLinkedList();

  CircDLinkedList operator++();

  //operator=(CircDLinkedList( CircDLinkedList&));
  //operator==( CircDLinkedList);

  typedef Iterator<CircDLinkedList<T> > LIterator;

  IndexType startPos() { return _start; }
  IndexType endPos() { return _start; }

  LIterator start() { return LIterator(*this,_start); }
  LIterator end()   { return LIterator(*this,_start); }

  void setStart(const LIterator& i);


  T& get(IndexType pos) { return pos->data; }

  void add(const T&); //inserts at before end()
  void remove(const Iterator<CircDLinkedList<T> >&);

  Letter* next(Letter* p) {  return p->next;  }
  Letter* prev(Letter* p) { return p->prev; }

  Letter* _start;

#ifdef TESTING
  static int test();
#endif
};

template<class CONTAINER>
class Iterator {
  mutable CONTAINER& _container;
  typename CONTAINER::IndexType _pos;

  Iterator();
public:
  typedef typename CONTAINER::ElemType  ElemType;
  typedef typename CONTAINER::IndexType IndexType;

  Iterator(CONTAINER&,IndexType=0);
  Iterator(const Iterator&);
  Iterator& operator=(const Iterator&);

  Iterator& operator++();
  Iterator& operator--();

  ElemType& current();
  Iterator next() { Iterator<CONTAINER> r(*this); ++r; return r; }
  Iterator prev() { Iterator<CONTAINER> r(*this); --r; return r; }

  CONTAINER& container() { return _container; }
  IndexType& pos() { return _pos; }
  CONTAINER& container() const { return _container; }
  const IndexType& pos() const { return _pos; }


  //reference methods
  ElemType* operator->() { return &current(); }
  //  ElemType& operator()() { return current(); }
  ElemType& operator*() { return current(); }

  const ElemType& operator*() const { return current(); }


  bool operator==(const Iterator<CONTAINER>&);
  bool operator!=(const Iterator<CONTAINER>& i)
  { return !operator==(i); }
  //  bool operator==(const ElemType&);
};

}

#endif /* __circDLinkedList_hh__ */
