
#ifndef __cdsMap_hh__
#define __cdsMap_hh__

#include <cdsList.h>

template<class Key, class Atom>
class CDSMap;

template<class Key, class Atom>
struct CDSMapNode {
  Key key;
  Atom value;
  //CDSMapNode(const Key& key) : key(key) {}
  CDSMapNode(const Key& key,
	     const Atom& value) : key(key), value(value) {}
private:
  CDSMapNode() {}
};

template<class Key, class Atom>
class CDSMapRep {
  CDSList<CDSMapNode<Key,Atom> > list;
  int count;

  CDSMapRep(const CDSMapRep&);  //inaccessible
  CDSMapRep<Key,Atom>& operator=(const CDSMapRep&);
  CDSMapRep(const CDSMapRep* r) : 
    count(1)
  { 
   list = r->list;
  }
public:
  CDSMapRep() : count(1) {} 
  ~CDSMapRep() {assert( count==0 ); };
  friend class CDSMap<Key,Atom>;
};

template<class Key, class Atom>
class CDSMap {
  CDSMapRep<Key,Atom>* rep;
protected:
  inline void splitRep() 
  { if (rep->count>1) {rep->count--;rep= new CDSMapRep<Key,Atom>(rep);}}
  bool findSlot(const Key& key  ,
		      int& index) const;
public:
  CDSMap() : rep(new CDSMapRep<Key,Atom>) {}
  ~CDSMap();
  CDSMap(const  CDSMap<Key,Atom>&); //copy constructor
  CDSMap& operator=(const CDSMap<Key,Atom>&);

  void add(const Key&,
	   const Atom&);  // does not require default constructor
  Atom& operator[](const Key&); //requires Atom to have default constructor
  const Atom&  operator[](const Key&) const;
  CDSList<Key> keys() const;
  void remove(const Key&);
  int size()                  const { return rep->list.size(); }
  bool exists(const Key& key) const { int dum=0;return findSlot(key,dum); }

#ifdef TESTING
  static int test();
#endif
};

template<class Key, class Atom>
CDSMap<Key,Atom>::~CDSMap()
{
 assert( rep->count>0 );
 if (--rep->count <= 0) delete rep;
} /* CDSMapRep::~CDSMapRep */

#endif /* __cdsMap_hh__ */


