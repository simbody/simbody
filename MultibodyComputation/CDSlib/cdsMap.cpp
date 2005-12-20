

#include <sthead.h>
#include <cdsMap.h>
#include <cdsExcept.h>

template<class Key, class Atom> 
CDSMap<Key,Atom>::CDSMap(const CDSMap<Key,Atom>& m)
  //
  // copy constructor
  //
{
 rep = m.rep;
 rep->count++;
} /* CDSMap(const CDSMap&) */

template<class Key, class Atom> 
CDSMap<Key,Atom>& 
CDSMap<Key,Atom>::operator=(const CDSMap<Key,Atom> &m)
  //
  // assignment of list to list
  //
{
 if (&m==this)  //self-assignment
   return *this;

 m.rep->count++;
 if (--rep->count <= 0) delete rep;
 rep = m.rep;

 return *this;
} /* operator=(map) */

template<class Key, class Atom>
void
CDSMap<Key,Atom>::add(const Key& key,
		      const Atom& atom)
{
 splitRep();

 CDSList< CDSMapNode<Key,Atom> >& l = rep->list;
 int indx=0;
 
 if ( findSlot(key,indx) ) {
   l[indx] = CDSMapNode<Key,Atom>(key,atom);
 } else {
   indx++;
   if (indx<l.size() ) {
     l.append( l[l.size()-1] );
     for (int i=l.size()-2 ; i>indx ; i--)
       l[i] = l[i-1];
     l[indx] = CDSMapNode<Key,Atom>(key,atom);
   } else {
     indx = l.size();
     l.append( CDSMapNode<Key,Atom>(key,atom) );
   }
 }
} /* add */

template<class Key, class Atom>
Atom&
CDSMap<Key,Atom>::operator[](const Key& key)
{
 splitRep();

 CDSList< CDSMapNode<Key,Atom> >& l = rep->list;
 int indx=0;
 
 if ( !findSlot(key,indx) ) {
   indx++;
   if (indx<l.size() ) {
     l.append( l[l.size()-1] );
     for (int i=l.size()-2 ; i>indx ; i--)
       l[i] = l[i-1];
     l[indx] = CDSMapNode<Key,Atom>(key,Atom());
   } else {
     indx = l.size();
     l.append( CDSMapNode<Key,Atom>(key,Atom()) );
   }
 }
 return rep->list[indx].value;
} /* operator[] */

template<class Key, class Atom>
const Atom&
CDSMap<Key,Atom>::operator[](const Key& key) const
{
 int indx = 0;
 if ( findSlot(key,indx) )
   return rep->list[indx].value;
 else
   throw CDS::out_of_range("CDSMap::operator[] const");
} /* operator[] const */

template<class Key, class Atom> 
void
CDSMap<Key,Atom>::remove(const Key& key)
{
 splitRep();

 int indx = 0;
 if ( findSlot(key,indx) )
   rep->list.remove(indx);
} /* remove */

template<class Key, class Atom> 
bool
CDSMap<Key,Atom>::findSlot(const Key& key,
				 int& indx) const
  // find index for which l[indx] <= key
  // returns 0 for list of zero size
{
 CDSList< CDSMapNode<Key,Atom> >& l = rep->list;

 if ( l.size() == 0 ) {
   indx=0;
   return 0;
 }

 // use bisection
 int jl = 0;
 int jm = 0;
 int ju = l.size()+1;
 while ( ju-jl>1 ) {
   jm = (ju+jl)/2;
   if ( key < l[jm-1].key  )
     ju = jm;
   else 
     jl = jm;
 }
 indx = jl-1;
 if (jl>0 && l[indx].key==key)
   return 1;
 else
   return 0;
} /* findSlot */
 
template<class Key, class Atom> 
CDSList<Key> 
CDSMap<Key,Atom>::keys() const
{
 CDSList<Key> ret;
 for (int i=0 ; i<rep->list.size() ; i++)
   ret.append( rep->list[i].key );
 return ret;
} /* keys */
 

#ifdef TESTING

#include <cdsString.h>
#include <cdsSStream.h>

template<class Key, class Atom>
int
CDSMap<Key,Atom>::test()
{
 int exit=0;
 cout << "testing CDSMap...";

 CDSMap<String,int> m;
 m["abcd"] = 45;
 
 if ( m["abcd"] != 45 ) {
   cerr << "first test failed.\n";
   exit=1;
 }

 m["xyz"] = 2;
 if ( m["xyz"] != 2 ) {
   cerr << "second test failed.\n";
   exit=1;
 }

 CDSMap<String,int> m2 = m;
 if ( m2["abcd"] != 45 ) {
   cerr << "copy constructor failed.\n";
   exit=1;
 }

 if ( ! m.exists("xyz") || m.exists("wxyz") ) {
   cerr << "exists test failed.\n";
   exit=1;
 }

 CDSList<String> l = m.keys();
 if ( l[0] != "abcd" ||
      l[1] != "xyz"    ) {
   cerr << "keys test failed: " << m.keys() << "\n";
   exit=1;
 }

 m.remove("xyz");
 if ( ! m.exists("abcd") || m.exists("xyz") ) {
   cerr << "remove test failed.\n";
   exit=1;
 }

 m["bbb"] = 33;
 m["aaa"] = 39;
 if ( ! m.exists("abcd") || !m.exists("aaa") ) {
   cerr << "insert test failed: " << m.keys() << "\n";
   exit=1;
 }

 const CDSMap<String,int> cm(m);
 if ( cm["bbb"] != 33 ) {
   cerr << "constant accessor test failed: " << cm["bbb"] << "\n";
   exit=1;
 }
 
 { // test add
   CDSMap<String,int> m;
   m.add("b",41);
   m.add("a",42);
   m.add("b",43);
   if ( m["a"] != 42 ||
	m["b"] != 43   ) {
     cerr << "add test failed: " << m["a"] << ' ' << m["b"] << "\n";
     exit=1;
   }
 }
 

 CDSMap<double,double> mp;
 mp[0.1] = 0.99;
 mp[0.2] = 0.99;
 mp[0.0] = 0.99;
 mp[0.3] = 0.45;
 CDSList<double> keys = mp.keys();
 if ( keys[0] != 0.0 || keys[1] != 0.1 ||
      keys[2] != 0.2 || keys[3] != 0.3  ) {
   cerr << "double keys test failed: " << mp.keys() << "\n";
   exit=1;
 }

 {
   CDSMap<double,double> mp;
   mp[0.1] = 0.99;
   mp[0.2] = 0.99;
   mp[0.3] = 0.45;
   mp[0.0] = 0.99;
   CDSList<double> keys = mp.keys();
   if ( keys[0] != 0.0 || keys[1] != 0.1 ||
	keys[2] != 0.2 || keys[3] != 0.3  ) {
     cerr << "key order test failed: " << mp.keys() << "\n";
     exit=1;
   }
 }

#ifdef DEBUG_ALLOC
 { // stress tests
   CDSMap<int,int> map;
   for (int i=0 ; i<8 ; i++) {
     // 8 x 1e5 x sizeof(int)
     CDSList<int> l;
     for (int j=0 ; j<1000000 ; j++ )
       map[j] = j;
     cout << i << ' ';
     cout.flush();
   }
 }
#endif /* DEBUG_ALLOC */

 cout << (exit?"failed":"ok") << endl;
 return exit;
} /* test */

#endif /* TESTING */


