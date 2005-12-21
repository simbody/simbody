
#ifndef __chull_hh__
#define __chull_hh__

#include <cdsList.h>
#include <fixedVector.h>

template<class VERTEX>
class CHullImpl;

// requirement for class VERTEX:
//   methods: methods of class CDSVec3
//   a typedef named Point3 which defines the Vertex type
//   an integer member named index

template<class VERTEX>
class CHull {
  CHullImpl<VERTEX>* chullp;
public:
  CHull();
  ~CHull();
  
  typedef FixedVector<int,3> Face;

  bool debug() const;
  bool check() const;

  void setDebug(bool);
  void setCheck(bool);

  void create();

  bool checkHull();

  void addVertex(const VERTEX&);

  CDSList<VERTEX> vertices() const;
  //returns list of index triplets
  //  these indices refer to the vertex list returned by vertices()
  CDSList<Face> faces() const;  

};

#ifdef TESTING
int testCHull();
#endif

#endif /* __chull_hh__ */
