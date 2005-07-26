
#ifndef __cdsPair_h__
#define __cdsPair_h__


//
// A pair is a set of two integers 
// in which member order is irrelevant
//

class Pair {  

public: 

  //
  // instance vbls
  //

  int a,b;

  //
  // constructors
  //

  Pair() : a(0), b(0) {}

  Pair(int a,
       int b) : a(a), b(b) {}

  //
  // equals operator
  //

  bool operator==(const Pair& p) {
    return ((a==p.a && b==p.b) ||
	    (a==p.b && b==p.a)  );
  }

};

#endif /* __cdsPair_h__ */
