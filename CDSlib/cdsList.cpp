#include <sthead.h>
#include <cdsList.h>
#include <cdsMath.h>


#include <cstdlib>
#include <cstring>
#include <cdsSStream.h>

int
CDSList_test()
{
 int exit=0;
 cout << "testing CDSList...";

 CDSList<int> l;
 CDSList<int> l2(1,2);
 CDSList<int> l3 = l2;

 l.append(1);
 
 if (l[0] != 1) {
   cerr << "index failed: " << l[0] << " != 1\n";
   exit=1;
 }

 if (l3.size()!=1 || l3.getRep()->asize!=2) {
   cerr << "internal error 1\n";
   exit=1;
 }

 {
   CDSList<int> l2; l2.append(5); l2.append(7);
   CDSList<int> l3; l3 = l2.vector();
   if (l3.size()!=2 || l3[0]!=5 || l3[1] !=7 ) {
     cerr << "problem with assignment from generic vector:\n";
     cerr << '\t' << l2 << " != " << l3 << '\n';
     exit=1;
   }
 }

 l3.resize(3);
 if (l3.size()!=3 || l3.getRep()->asize!=5) {
   cerr << "internal error 2\n";
   exit=1;
 }
 
 l.append(4);
 if ( !l.contains(4) ) {
   cerr << "contains error\n";
   exit=1;
 }

 if ( l.getIndex(4)!=1 ) {
   cerr << "getIndex error\n";
   exit=1;
 }

 OStringStream os; os << l << ends;
 if ( strcmp(os.str(),"{ 1, 4 }") ) {
   cerr << "operator<< error: " << os.str() << " != { 1, 4 }\n";
   exit=1;
 }

 l.prepend(5);
 if ( l[0] != 5 || l[1] != 1 || l[2] != 4) {
   cerr << "prepend error: " << l << " != { 5, 1, 4 }\n";
   exit=1;
 }

 l.remove(0);
 if ( l[0] != 1 || l[1] != 4) {
   cerr << "remove error: " << l << " != { 1, 4 }\n";
   exit=1;
 }

 int t=l.pop();
 if ( t != 4 || l.size()!=1 ) {
   cerr << "pop error: " << t << " " << l.size() << '\n';
   exit=1;
 }

 l.resize(0);
 l.append(42);
 if ( l.size() != 1 || l[0]!=42 ) {
   cerr << "resize(0) error: " << l << " != { 42 }\n";
   exit=1;
 }

 {
   CDSList<int> l(0,0);
   l.append(42);
   if ( l.size()!= 1 || l[0] != 42 || l.getRep()->asize<1){
     cerr << "CDSList(0,0) error on append: asize= " << l.getRep()->asize << '\n';
     exit=1;
   }
 }

 {
   CDSList< CDSList<int> > ll(2);
   ll[0].append(3); ll[0].append(4); 
   ll[1].append(5); 
   CDSList< CDSList<int> > ll2 = ll;
   ll.resize(0);
   if ( ll.size() || (ll2[0][0] != 3) || (ll2[0][1] !=4) || (ll2[1][0] !=5) ) {
     cerr << "error: ll: " << ll 
	  << "\tll2: " << ll2 << endl;
     exit=1;
   }
 }


 {
   CDSList<int> l;
   l.append(3); l.append(5); l.append(2); l.append(1);
   l.sort_si(CDSList<int>::stdComparer);
   if ( l[0] != 1 ||
	l[1] != 2 ||
	l[2] != 3 ||
	l[3] != 5   ) {
     cerr << "sort_si error: " << l << " != { 1, 2, 3, 5 }\n";
     exit=1;
   }
 }

 {
   CDSList<int> l;
   l.append(3); l.append(5); l.append(2); l.append(1);
   l.sort_shell(CDSList<int>::stdComparer);
   if ( l[0] != 1 ||
	l[1] != 2 ||
	l[2] != 3 ||
	l[3] != 5   ) {
     cerr << "sort_shell error: " << l << " != { 1, 2, 3, 5 }\n";
     exit=1;
   }
 }

 {
   CDSList<int> l;
   l.append(3); l.append(5); l.append(2); l.append(1);
   l.sort_heap(CDSList<int>::stdComparer);
   if ( l[0] != 1 ||
	l[1] != 2 ||
	l[2] != 3 ||
	l[3] != 5   ) {
     cerr << "sort_heap error: " << l << " != { 1, 2, 3, 5 }\n";
     exit=1;
   }
 }

 {
   CDSList<int> l;
   l.append(3); l.append(5); l.append(2); l.append(1);
   l.reverse();
   if ( l[0] != 1 ||
	l[1] != 2 ||
	l[2] != 5 ||
	l[3] != 3   ) {
     cerr << "reverse error: " << l << " != { 1, 2, 5, 3 }\n";
     exit=1;
   }
 }

#ifdef DEBUG_ALLOC
 { // two stress tests
   for (int i=0 ; i<8 ; i++) {
     // 8 x 1e7 x sizeof(int)
     CDSList<int> l;
     for (int j=0 ; j<10000000 ; j++ )
       l.append(j);
     cout << i << ' ';
     cout.flush();
   }
   CDSList<int> l;
   for (l_int i=0 ; i<80 ; i++) {
     l.resize(100000000);
     l[99999999] = i;
     l[0] = 0;
   }
 }
#endif /* DEBUG_ALLOC */


 cout << (exit?"failed":"ok") << endl;
 return exit;
} /* test */
