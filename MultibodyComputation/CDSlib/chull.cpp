/*
This code is described in "Computational Geometry in C" (Second Edition),
Chapter 4.  It is not written to be comprehensible without the 
explanation in that book.

Input: 3n integer coordinates for the points.
Output: the 3D convex hull, in postscript with embedded comments
        showing the vertices and faces.

Compile: gcc -o chull chull.c (or simply: make)

Written by Joseph O'Rourke, with contributions by 
  Kristy Anderson, John Kutcher, Catherine Schevon, Susan Weller.
Last modified: May 2000
Questions to orourke@cs.smith.edu.

--------------------------------------------------------------------
This code is Copyright 2000 by Joseph O'Rourke.  It may be freely 
redistributed in its entirety provided that this copyright notice is 
not removed.
--------------------------------------------------------------------
*/

#include <sthead.h>
#include <chull.h>
#include <fixedVector.h>
#include <circDLinkedList.h>
#include <cdsString.h>
#include <cdsMap.h>
#include <cdsIomanip.h>

// constant flags
const int ONHULL   	=TRUE;
const int REMOVED  	=TRUE;
const int VISIBLE  	=TRUE;
const int PROCESSED	=TRUE;

template<class VERTEX>
class CHullImpl {
public:
  CHullImpl() : debug(0), check(0) {}

  typedef typename VERTEX::Point3 Point3;

  struct Edge;
  struct Face;

  struct Vertex {
    static int count;
  public:
    Vertex(const VERTEX& v) :
      v(v), vnum(count++), duplicate(0), onhull(!ONHULL), mark(!PROCESSED) {}

    
    VERTEX v;
    int	  vnum;
    Edge*   duplicate;	        /* pointer to incident cone edge (or NULL) */
    bool    onhull;		/* T iff point on hull. */
    bool	  mark;		/* T iff point already processed. */
  };

  struct Edge {
    Edge() :
      newface(NULL), remove(!REMOVED) 
    { adjface[0] = adjface[1] = NULL; endpts[0] = endpts[1] = NULL; }
    
    
    FixedVector<Face*,2>    adjface;
    FixedVector<Vertex*,2>  endpts;
    Face*   newface;            /* pointer to incident cone face. */
    bool     remove;		/* T iff edge should be remove. */
  };

  class Face {
  public:
    Face() : visible(!VISIBLE) 
    { 
     edge[0] = edge[1] = edge[2] = NULL;
     vertex[0] = vertex[1] = vertex[2] = NULL;
    } 
    
    FixedVector<Edge*,3>    edge;
    FixedVector<Vertex*,3>  vertex;
    bool                    visible;   // T iff face visible from new point.
  };

  typedef CDS::CircDLinkedList<Vertex> VertexList;
  typedef CDS::CircDLinkedList<Edge>   EdgeList;
  typedef CDS::CircDLinkedList<Face>   FaceList;

  typedef CDS::Iterator<VertexList> VIterator;
  typedef CDS::Iterator<EdgeList>   EIterator;
  typedef CDS::Iterator<FaceList>   FIterator;

  VertexList vertices;
  EdgeList   edges;
  FaceList   faces;

  bool debug;
  bool check;

  void doubleTriangle();
  void constructHull();
  bool addOne(Vertex& p);
  int  volumeSign(const Face&   f, 
		  const Vertex& p);
  float_type volume(const Face&   f,
		    const Vertex& p);
  Face*	makeConeFace(Edge&   e, 
		     Vertex& p );
  void    makeCcw(Face&   f, 
		  Edge&   e,
		  Vertex& p);
  Edge*  makeNullEdge();
  Face*  makeNullFace( void );
  Face*  makeFace(Vertex& v0, 
		  Vertex& v1, 
		  Vertex& v2, 
		  Face*   fp);
  void    cleanUp(VIterator& pvnext);
  void    cleanEdges();
  void    cleanFaces();
  void    cleanVertices(VIterator& pvnext);
  bool  checkEuler(int V, int E, int F );
  bool  consistency();
  bool  convexity();
  bool  checkEndpts();
  bool    checks();
  void	printOut(const Vertex& v);
  void	printVertices();
  void	printEdges();
  void	printFaces();
  void	edgeOrderOnFaces();

  friend class CHull<VERTEX>;
};





template<class VERTEX>
int CHullImpl<VERTEX>::Vertex::count=0;

template<class VERTEX>
CHull<VERTEX>::CHull() 
{
 chullp = new CHullImpl<VERTEX>();
}

template<class VERTEX>
CHull<VERTEX>::~CHull() 
{
 delete chullp;
}


template<class VERTEX>
static bool
collinear(const VERTEX& a,
	  const VERTEX& b, 
	  const VERTEX& c )
  //
  //Collinear checks to see if the three points given are collinear,
  //by checking to see if each element of the cross product is zero.
  //
{
 return  ( c.v.z() - a.v.z() ) * ( b.v.y() - a.v.y() ) -
         ( b.v.z() - a.v.z() ) * ( c.v.y() - a.v.y() ) == 0
      && ( b.v.z() - a.v.z() ) * ( c.v.x() - a.v.x() ) -
         ( b.v.x() - a.v.x() ) * ( c.v.z() - a.v.z() ) == 0
      && ( b.v.x() - a.v.x() ) * ( c.v.y() - a.v.y() ) -
         ( b.v.y() - a.v.y() ) * ( c.v.x() - a.v.x() ) == 0  ;
}

template<class VERTEX>
void
CHullImpl<VERTEX>::doubleTriangle()
  // doubleTriangle builds the initial double triangle.  It first finds 3 
  // noncollinear points and makes two faces out of them, in opposite order.
  // It then finds a fourth point that is not coplanar with that face.  The  
  // vertices are stored in the face structure in counterclockwise order so 
  // that the volume between the face and the point is negative. Lastly, the
  // 3 new faces to the fourth point are constructed and the data structures
  // are cleaned up. 
{
 /* Find 3 noncollinear points. */
 // VertexList v0 = vertices;
 VIterator v0 = vertices.start();
 while ( collinear( *v0, *v0.next(), *v0.next().next() ) )
   if ( ( ++v0 ) == vertices.end() ) {
     cout << "DoubleTriangle:  All points are Collinear!\n";
     throw CDS::exception("DoubleTriangle:  All points are Collinear!");
   }
 VIterator v1 = v0.next();
 VIterator v2 = v1.next();
 
 /* Mark the vertices as processed. */
 v0->mark = PROCESSED;
 v1->mark = PROCESSED;
 v2->mark = PROCESSED;
   
 /* Create the two "twin" faces. */
 Face* f1 = 0;
 Face* f0 = makeFace( *v0, *v1, *v2, f1 );
 f1 = makeFace( *v2, *v1, *v0, f0 );

 /* Link adjacent face fields. */
 f0->edge[0]->adjface[1] = f1;
 f0->edge[1]->adjface[1] = f1;
 f0->edge[2]->adjface[1] = f1;
 f1->edge[0]->adjface[1] = f0;
 f1->edge[1]->adjface[1] = f0;
 f1->edge[2]->adjface[1] = f0;
 
 /* Find a fourth, noncoplanar point to form tetrahedron. */
 VIterator v3 = v2.next();
 int vol = volumeSign( *f0, *v3 );
 while ( !vol )   {
   if ( ( v3 = v3.next() ) == v0 ) {
     cout << "DoubleTriangle:  All points are coplanar!\n";
     throw CDS::exception("DoubleTriangle:  All points are coplanar!");
   }
   vol = volumeSign( *f0, *v3 );
 }
	
 /* Insure that v3 will be the first added. */
 vertices.setStart(v3);
 if ( debug ) {
   cerr << "DoubleTriangle: finished. Head repositioned at v3.\n";
   printOut(*vertices.start() );
 }
} /* doubleTriangle */

  
template<class VERTEX>
void
CHullImpl<VERTEX>::constructHull()
  //
  // ConstructHull adds the vertices to the hull one at a time.  The hull
  // vertices are those in the list marked as onhull.
  //
{
 VIterator v = vertices.start();
 do {
   VIterator vnext = v.next();
   if ( !v->mark ) {
     v->mark = PROCESSED;
     bool changed = addOne(*v );  // true if addition changes hull; not used.
     cleanUp(vnext);              // Pass down vnext in case it gets removed.
     
     if ( check ) {
       cerr << "ConstructHull: After Add of " << v->vnum << " & Cleanup:\n";
       checks();
     }
     if ( debug )
       printOut(*v);
   }
   v = vnext;
 } while ( v != vertices.start() );
} /* constructHull */

template<class VERTEX>
bool
CHullImpl<VERTEX>::addOne(Vertex& p    )
  //
  // AddOne is passed a vertex.  It first determines all faces visible from 
  // that point.  If none are visible then the point is marked as not 
  // onhull.  Next is a loop over edges.  If both faces adjacent to an edge
  // are visible, then the edge is marked for deletion.  If just one of the
  // adjacent faces is visible then a new face is constructed.
  //
{
 if ( debug ) {
   cerr << "AddOne: starting to add v" << p.vnum << ".\n";
   printOut(*vertices.start());
 }

 /* Mark faces visible from p. */
 bool	  vis = FALSE;
 FIterator f = faces.start();
 do {
   int vol = volumeSign( *f, p );
   if (debug)
     cerr << "addOne: faddr: " << setw(6) << &f
	  << "   paddr: " << setw(6) << &p 
	  << "   Vol = " << vol << "\n";
   if ( vol < 0 ) {
     f->visible = VISIBLE;  
     vis = TRUE;                      
   }
 } while ( ++f != faces.end() );

 /* If no faces are visible from p, then p is inside the hull. */
 if ( !vis ) {
   p.onhull = !ONHULL;  
   return FALSE; 
 }

 /* Mark edges in interior of visible region for deletion.
    Erect a newface based on each border edge. */
 EIterator e = edges.start();
 do {
   EIterator temp = e.next();
   if ( e->adjface[0]->visible && e->adjface[1]->visible )
     // e interior: mark for deletion.
     e->remove = REMOVED;
   else if ( e->adjface[0]->visible || e->adjface[1]->visible ) 
     /* e border: make a new face. */
     e->newface = makeConeFace(*e, p );
   e = temp;
 } while ( e != edges.end() );
 return TRUE;
} /* addOne */

template<class VERTEX>
int  
CHullImpl<VERTEX>::volumeSign(const Face&   f, 
			      const Vertex& p )
  //
  // VolumeSign returns the sign of the volume of the tetrahedron determined by
  // f and p.  VolumeSign is +1 iff p is on the negative side of f,
  // where the positive side is determined by the rh-rule.  So the volume 
  // is positive if the ccw normal to f points outside the tetrahedron.
  // The final fewer-multiplications form is due to Bob Williamson.
  //
{
 Point3 a = Point3(f.vertex[0]->v) - Point3(p.v);
 Point3 b = Point3(f.vertex[1]->v) - Point3(p.v);
 Point3 c = Point3(f.vertex[2]->v) - Point3(p.v);

 float_type vol = a.x() * (b.y()*c.z() - b.z()*c.y()) +
		  a.y() * (b.z()*c.x() - b.x()*c.z()) +
		  a.z() * (b.x()*c.y() - b.y()*c.x());

 if ( debug ) {
   /* Compute the volume using integers for comparison. */
   float_type voli = volume( f, p );
   cerr << "Face=" << setw(8) << &f << "; Vertex=" << p.vnum
	<< ": vol = " << vol << "\n";
 }

 /* FIX: The volume should be an integer. */
 const float_type smallNum = 1e-10;
 if      ( vol >  smallNum )  return  1;
 else if ( vol < -smallNum )  return -1;
 else                         return  0;
} /* volumeSign */

/*---------------------------------------------------------------------
Same computation, but computes using ints, and returns the actual volume.
---------------------------------------------------------------------*/
template<class VERTEX>
float_type
CHullImpl<VERTEX>::volume(const Face&   f, 
			  const Vertex& p)
{
 Point3 a = Point3(f.vertex[0]->v) - Point3(p.v);
 Point3 b = Point3(f.vertex[1]->v) - Point3(p.v);
 Point3 c = Point3(f.vertex[2]->v) - Point3(p.v);

 float_type vol =  (  a.x() * (b.y()*c.z() - b.z()*c.y())
		    + a.y() * (b.z()*c.x() - b.x()*c.z())
		    + a.z() * (b.x()*c.y() - b.y()*c.x()));

 return vol;
} /* volume */


template<class VERTEX>
typename CHullImpl<VERTEX>::Face*
CHullImpl<VERTEX>::makeConeFace(Edge&   e,
				Vertex& p)
//
// MakeConeFace makes a new face and two new edges between the 
// edge and the point that are passed to it. It returns a pointer to
// the new face.
//
{
 Edge*  new_edge[2];

 /* Make two new edges (if don't already exist). */
 for (int i=0 ; i<2 ; ++i ) 
   /* If the edge exists, copy it into new_edge. */
   if ( !( new_edge[i] = e.endpts[i]->duplicate) ) {
     /* Otherwise (duplicate is NULL), MakeNullEdge. */
     new_edge[i] = makeNullEdge();
     new_edge[i]->endpts[0] = e.endpts[i];
     new_edge[i]->endpts[1] = &p;
     e.endpts[i]->duplicate = new_edge[i];
   }

 /* Make the new face. */
 Face* new_face = makeNullFace();   
 new_face->edge[0] = &e;
 new_face->edge[1] = new_edge[0];
 new_face->edge[2] = new_edge[1];
 makeCcw( *new_face, e, p ); 
 
 /* Set the adjacent face pointers. */
 for (int i=0 ; i<2 ; ++i )
   for (int j=0 ; j<2 ; ++j )  
     /* Only one NULL link should be set to new_face. */
     if ( !new_edge[i]->adjface[j] ) {
       new_edge[i]->adjface[j] = new_face;
       break;
     }
        
 return new_face;
} /* makeConeFace */



template<class VERTEX>
void	
CHullImpl<VERTEX>::makeCcw(Face&   f,
			   Edge&   e, 
			   Vertex& p)
//
// MakeCcw puts the vertices in the face structure in counterclock wise 
// order.  We want to store the vertices in the same 
// order as in the visible face.  The third vertex is always p.
// 
// Although no specific ordering of the edges of a face are used
// by the code, the following condition is maintained for each face f:
// one of the two endpoints of f->edge[i] matches f->vertex[i]. 
// But note that this does not imply that f->edge[i] is between
// f->vertex[i] and f->vertex[(i+1)%3].  (Thanks to Bob Williamson.)
//
{
 Face* fv = e.adjface[1];   // The visible face adjacent to e 
 if  ( e.adjface[0]->visible )      
   fv = e.adjface[0];
       
 /* Set vertex[0] & [1] of f to have the same orientation
    as do the corresponding vertices of fv. */ 
 int    i;    /* Index of e->endpoint[0] in fv. */
 for (i=0 ; fv->vertex[i] != e.endpts[0] ; ++i ) 
   if (i>2)
     throw CDS::exception("CHullImpl::makeCcw: face/edge inconsistency");

 /* Orient f the same as fv. */
 if ( fv->vertex[ (i+1) % 3 ] != e.endpts[1] ) {
   f.vertex[0] = e.endpts[1];  
   f.vertex[1] = e.endpts[0];    
 } else {                               
   f.vertex[0] = e.endpts[0];   
   f.vertex[1] = e.endpts[1];      
   CDS::swap( f.edge[1], f.edge[2] );
 }
 /* This swap is tricky. e is edge[0]. edge[1] is based on endpt[0],
    edge[2] on endpt[1].  So if e is oriented "forwards," we
    need to move edge[1] to follow [0], because it precedes. */
 
 f.vertex[2] = &p;
} /* makeCcw */
 
template<class VERTEX>
typename CHullImpl<VERTEX>::Edge*
CHullImpl<VERTEX>::makeNullEdge()
  //
  // MakeNullEdge creates a new cell and initializes all pointers to NULL
  // and sets all flags to off.  It returns a pointer to the empty cell.
  //
{
 Edge e;
 edges.add(e);
 return &(*edges.end().prev());
} /* makeNullEdge */

template<class VERTEX>
typename CHullImpl<VERTEX>::Face*
CHullImpl<VERTEX>::makeNullFace()
  //
  // MakeNullFace creates a new face structure and initializes all of its
  // flags to NULL and sets all the flags to off.  It returns a pointer
  // to the empty cell.
  //
{
 Face f;
 faces.add(f);
 return &(*faces.end().prev());
} /* makeNullFace */

template<class VERTEX>
typename CHullImpl<VERTEX>::Face*
CHullImpl<VERTEX>::makeFace(Vertex& v0, 
			    Vertex& v1, 
			    Vertex& v2, 
			    Face*   fold )
  //
  // MakeFace creates a new face structure from three vertices (in ccw
  // order).  It returns a pointer to the face.
  //
{
 /* Create edges of the initial triangle. */
 Edge  *e0, *e1, *e2;
 if( !fold ) {
   e0 = makeNullEdge();
   e1 = makeNullEdge();
   e2 = makeNullEdge();
 }
 else { /* Copy from fold, in reverse order. */
   e0 = fold->edge[2];
   e1 = fold->edge[1];
   e2 = fold->edge[0];
 }
 e0->endpts[0] = &v0;              e0->endpts[1] = &v1;
 e1->endpts[0] = &v1;              e1->endpts[1] = &v2;
 e2->endpts[0] = &v2;              e2->endpts[1] = &v0;

 /* Create face for triangle. */
 Face* f = makeNullFace();
 f->edge[0]   = e0;  f->edge[1]   = e1; f->edge[2]   = e2;
 f->vertex[0] = &v0;  f->vertex[1] = &v1; f->vertex[2] = &v2;
	
 /* Link edges to face. */
 e0->adjface[0] = e1->adjface[0] = e2->adjface[0] = f;
	
 return f;
} /* makeFace */

template<class VERTEX>
void
CHullImpl<VERTEX>::cleanUp(VIterator& pvnext)
  //
  // CleanUp goes through each data structure list and clears all
  // flags and NULLs out some pointers.  The order of processing
  // (edges, faces, vertices) is important.
  //
{
 cleanEdges();
 cleanFaces();
 cleanVertices(pvnext);
} /* cleanUp */

template<class VERTEX>
void
CHullImpl<VERTEX>::cleanEdges()
  //
  // CleanEdges runs through the edge list and cleans up the structure.
  // If there is a newface then it will put that face in place of the 
  // visible face and NULL out newface. It also removes so marked edges.
  //
{
 /* Integrate the newface's into the data structure. */
 /* Check every edge. */
 EIterator e = edges.start();
 do {
   if ( e->newface ) { 
     if ( e->adjface[0]->visible )
       e->adjface[0] = e->newface; 
     else	e->adjface[1] = e->newface;
     e->newface = NULL;
   }
 } while ( ++e != edges.end() );

 /* Remove any edges marked for deletion. */
 while ( edges.start().pos() && edges.start()->remove ) { 
   EIterator e = edges.start();
   edges.setStart( edges.start().next() );
   edges.remove(e);
 }
 e = edges.start().next();
 do {
   if ( e->remove ) {
     EIterator t = e;
     ++e;
     edges.remove(t);
   } else
     ++e;
 } while ( e != edges.end());

} /* cleanEdges */

template<class VERTEX>
void
CHullImpl<VERTEX>::cleanFaces()
  //
  // CleanFaces runs through the face list and removes any face marked visible.
  //
{
 while ( faces.start().pos() && faces.start()->visible ) { 
   FIterator f = faces.start();
   faces.setStart( f.next() );
   faces.remove(f);
 }
 FIterator f = faces.start().next();
 do {
   if ( f->visible ) {
     FIterator t = f;
     ++f;
     faces.remove(t);
   }
   else 
     ++f;
 } while ( f != faces.end() );
} /* cleanFaces */

template<class VERTEX>
void
CHullImpl<VERTEX>::cleanVertices(VIterator& pvnext)
  //
  // CleanVertices runs through the vertex list and deletes the 
  // vertices that are marked as processed but are not incident to any 
  // undeleted edges. 
  // The pointer to vnext, pvnext, is used to alter vnext in
  // ConstructHull() if we are about to delete vnext.
{
 /* Mark all vertices incident to some undeleted edge as on the hull. */
 EIterator e(edges);
 do {
   e->endpts[0]->onhull = e->endpts[1]->onhull = ONHULL;
 } while (++e != edges.end());
	
 /* Delete all vertices that have been processed but
    are not on the hull. */
 VIterator v(vertices);
 while ( v.pos() && v->mark && !v->onhull ) { 
   /* If about to remove vnext, advance it first. */
   if ( v == pvnext ) {
     pvnext = v.next();
   }
   vertices.setStart( v.next() );
   vertices.remove(v);
   v = vertices.start();
 }

 v = vertices.start().next();
 do {
   if ( v->mark && !v->onhull ) {    
     VIterator t = v; 
     ++v;
     vertices.remove(t);
   } else 
     ++v;
 } while ( v != vertices.end() );
	
 /* Reset flags. */
 v = vertices.start();
 do {
   v->duplicate = NULL; 
   v->onhull = !ONHULL; 
 } while ( ++v != vertices.end() );
} /* cleanVertices */


template<class VERTEX>
bool
CHullImpl<VERTEX>::consistency()
  //
  // Consistency runs through the edge list and checks that all
  // adjacent faces have their endpoints in opposite order.  This verifies
  // that the vertices are in counterclockwise order.
  //
{
 EIterator e = edges.start();
 do {
   /* find index of endpoint[0] in adjacent face[0] */
   int i;
   for (i = 0; e->adjface[0]->vertex[i] != e->endpts[0]; ++i )
     ;
   
   /* find index of endpoint[0] in adjacent face[1] */
   int j;
   for (j = 0; e->adjface[1]->vertex[j] != e->endpts[0]; ++j )
     ;

   /* check if the endpoints occur in opposite order */
   if ( !( e->adjface[0]->vertex[ (i+1) % 3 ] ==
	   e->adjface[1]->vertex[ (j+2) % 3 ] ||
	   e->adjface[0]->vertex[ (i+2) % 3 ] ==
	   e->adjface[1]->vertex[ (j+1) % 3 ] )  )
     break;
 } while ( ++e != edges.end() );

 bool ret = 1;
 if ( e != edges.start() ) {
   cerr << "Chull::Checks: edges are NOT consistent.\n";
   ret=0;
 }
 return ret;
} /* consistency */


template<class VERTEX>
bool
CHullImpl<VERTEX>::convexity()
  //
  // Convexity checks that the volume between every face and every
  // point is negative.  This shows that each point is inside every face
  // and therefore the hull is convex.
{
 FIterator f = faces.start();
 do {
   VIterator v = vertices.start();
   do {
     if ( v->mark ) {
       int vol = volumeSign( *f, *v );
       if ( vol < 0 )
	 break;
     }
   } while ( ++v != vertices.end() );
   
 } while ( ++f != faces.end() );

 bool ret=1;
 if ( f != faces.start() ) {
   cerr << "CHull::Checks: NOT convex.\n";
   ret=0;
 }
 return ret;
} /* convexity */

template<class VERTEX>
bool
CHullImpl<VERTEX>::checkEuler(int V, 
			      int E, 
			      int F )
  //
  // CheckEuler checks Euler's relation, as well as its implications when
  // all faces are known to be triangles.  Only prints positive information
  // when debug is true, but always prints negative information.
  //
{
 if ( check )
   cerr <<  "CHull::Checks: V, E, F = " << V << ' ' << E << ' '  
	<< F << "    ";

 bool ret=1;
 if ( (V - E + F) != 2 ) { 
   cerr << "Checks: V-E+F != 2\n";
   ret=0;
 } else if (check) cerr << "V-E+F = 2\t";

 if ( F != (2 * V - 4) ) {
   cerr << "Checks: F=" << F 
	<< " != 2V-4=" << (2*V-4) << "; V=" << V << "\n";
   ret=0;
 } else if (check) cerr << "F = 2V-4  ";
 
 if ( (2 * E) != (3 * F) ) {
   cerr << "Checks: 2E=" << (2*E) << " != 3F=" << (3*F) 
	<< "; E=" << E <<  ", F=" << F << '\n';
   ret=0;
 } else if (check) cerr << "2E = 3F\n";
 return ret;
} /* checkEuler */

template<class VERTEX>
bool
CHullImpl<VERTEX>::checks()
{
 bool ret=1;
 
 if ( !consistency() ) ret = 0;
 if ( !convexity()   ) ret = 0;
 
 int V=0;
 VIterator v = vertices.start();
 if ( v.pos() )
   do if (v->mark) V++; while ( ++v != vertices.end() );
 
 int E=0;
 EIterator e = edges.start();
 if ( e.pos() )
   do E++; while ( ++e != edges.end() );
 
 int F=0;
 FIterator f = faces.start();
 if ( f.pos() )
   do F++; while ( ++f  != faces.end() );
 
 if ( !checkEuler( V, E, F ) ) ret = 0;
 if ( !checkEndpts()         ) ret = 0;

 return ret;
} /* checks */


/*===================================================================
These functions are used whenever the debug flag is set.
They print out the entire contents of each data structure.  
Printing is to standard error.  To grab the output in a file in the csh, 
use this:
	chull < i.file >&! o.file
=====================================================================*/
/*-------------------------------------------------------------------*/
template<class VERTEX>
void
CHullImpl<VERTEX>::printOut(const Vertex& v)
{
 cerr << "\nHead vertex " << v.vnum << " = " << setw(8) << &v << ":\n";
 printVertices();
 printEdges();
 printFaces();
} /* printOut */

/*-------------------------------------------------------------------*/
template<class VERTEX>
void
CHullImpl<VERTEX>::printVertices()
{
 cerr << "Vertex List\n";
 VIterator v = vertices.start();
 if (v.pos()) do {
   cerr << "  addr " << setw(6) << &(v.current()) << "  "
	<< "  vnum " << setw(4) << v->vnum
	<< " (" 
	<< setw(6) << setprecision(2) << v->v.x() << ','
	<< setw(6) << setprecision(2) << v->v.x() << ','
	<< setw(6) << setprecision(2) << v->v.x() << ")"
	<< "  active:" << setw(3) << v->onhull
	<< "  dup:" << setw(5) << v->duplicate
	<< "  mark:" << setw(2) << v->mark <<'\n';
 } while ( ++v != vertices.end() );
} /* printVertices */

/*-------------------------------------------------------------------*/
template<class VERTEX>
void
CHullImpl<VERTEX>::printEdges()
{

 cerr << "Edge List\n";
 EIterator e = edges.start();
 if (e.pos()) do {
   cerr << "  addr " << setw(8) << &(*e)
	<< " adj " << setw(8) << &(e->adjface[0]) << ' ' << &(e->adjface[1]) 
	<< "   endpts:" 
	<< setw(5) << e->endpts[0]->vnum << ' ' << e->endpts[1]->vnum
	<< "  del:" << setw(3) << e->remove << '\n';
 } while ( ++e != edges.end() );
} /* printEdges */

template<class VERTEX>
void
CHullImpl<VERTEX>::printFaces()
{
 cerr << "Face List\n";
 FIterator f = faces.start();  
 if (f.pos()) do {
   cerr << "  addr: " << setw(8) << &(*f)
	<< "    edges:";
   for(int i=0 ; i<3 ; ++i )
     cerr << setw(8) << f->edge[i] << ' ' ;
   cerr << "  vert:";
   for (int i=0 ; i<3 ; ++i)
     cerr << setw(4) << f->vertex[i]->vnum;
   cerr << "  vis: " << f->visible << '\n';
 } while ( ++f != faces.end() );
} /* printFaces */

template<class VERTEX>
bool
CHullImpl<VERTEX>::checkEndpts()
  //
  // Checks that, for each face, for each i={0,1,2}, the [i]th vertex of
  // that face is either the [0]th or [1]st endpoint of the [ith] edge of
  // the face.
  //
{
 bool error = FALSE;
 FIterator f = faces.start();
 if (f.pos()) do {
   for(int i=0 ; i<3 ; ++i) {
     Vertex* v = f->vertex[i];
     Edge*   e = f->edge[i];
     if ( v != e->endpts[0] && v != e->endpts[1] ) {
       error = TRUE;
       cerr << "CheckEndpts: Error!\n";
       cerr << "  addr: " << setw(8) << &f
	    << "  edges:(" 
	    << setw(3) << e->endpts[0]->vnum << ',' << e->endpts[1]->vnum
	    << ")\n";
     }
   }
 } while ( ++f != faces.end() );

 if ( error ) {
   cerr << "Checks: ERROR found and reported above.\n";
   return 0;
 } else
   return 1;
} /* checkEndpts */

template<class VERTEX>
void
CHullImpl<VERTEX>::edgeOrderOnFaces()
  //
  //   EdgeOrderOnFaces: puts e0 between v0 and v1, e1 between v1 and v2,
  //   e2 between v2 and v0 on each face.  This should be unnecessary, alas.
  //   Not used in code, but useful for other purposes.
  //
{
  FIterator f(faces);
  do {
    for (int i = 0; i < 3; i++) {
      if (!(((f->edge[i]->endpts[0] == f->vertex[i]) &&
             (f->edge[i]->endpts[1] == f->vertex[(i+1)%3])) ||
            ((f->edge[i]->endpts[1] == f->vertex[i]) &&
             (f->edge[i]->endpts[0] == f->vertex[(i+1)%3])))) {
        /* Change the order of the edges on the face: */
        for (int j = 0; j < 3; j ++) {
          /* find the edge that should be there */
          if (((f->edge[j]->endpts[0] == f->vertex[i]) &&
               (f->edge[j]->endpts[1] == f->vertex[(i+1)%3])) ||
              ((f->edge[j]->endpts[1] == f->vertex[i]) &&
               (f->edge[j]->endpts[0] == f->vertex[(i+1)%3]))) {
            /* Swap it with the one erroneously put into its place: */
            if ( debug )
	      cerr << "Making a swap in EdgeOrderOnFaces: "
		   << "F(" 
		   << f->vertex[0]->vnum << ','
		   << f->vertex[1]->vnum << ','
		   << f->vertex[2]->vnum << "): e["
		   << i << "] and e[" << j << "]\n";
	    CDS::swap( f->edge[i] , f->edge[j] );
          }
        }
      }
    }
  } while (++f != faces.end());
} /* edgeOrderOnFaces */

template<class VERTEX>
bool
CHull<VERTEX>::debug() const
{
 return chullp->debug;
}

template<class VERTEX>
bool
CHull<VERTEX>::check() const
{
 return chullp->check;
}

template<class VERTEX>
void
CHull<VERTEX>::setDebug(bool b)
{
 chullp->debug = b;
}

template<class VERTEX>
void
CHull<VERTEX>::setCheck(bool b)
{
 chullp->check = b;
}

template<class VERTEX>
void
CHull<VERTEX>::addVertex(const VERTEX& v)
{
 chullp->vertices.add( CHullImpl<VERTEX>::Vertex(v) );
} /* addVertex */

template<class VERTEX>
void
CHull<VERTEX>::create()
{
 chullp->doubleTriangle();
 chullp->constructHull();
 chullp->edgeOrderOnFaces();
} /* create */

template<class VERTEX>
CDSList<VERTEX>
CHull<VERTEX>::vertices() const
{
 CDSList<VERTEX> ret;

 typename CHullImpl<VERTEX>::VIterator i(chullp->vertices);
 do {
   ret.append( i->v );
   ret[ ret.size()-1 ].index = i->vnum;
 } while ( ++i != chullp->vertices.end() );
 return ret;
} /* vertices */

template<class VERTEX>
CDSList<typename CHull<VERTEX>::Face> 
CHull<VERTEX>::faces() const
{
 CDSList<Face> ret;

 typename CHullImpl<VERTEX>::VIterator i(chullp->vertices);
 CDSMap<int,int> indexMap;
 int cnt=0;
 do {
   indexMap[i->vnum] = cnt;
   cnt++;
 } while ( ++i != chullp->vertices.end() );

 typename CHullImpl<VERTEX>::FIterator f(chullp->faces);
 do {
   Face face;
   face[0] = indexMap[f->vertex[0]->vnum];
   face[1] = indexMap[f->vertex[1]->vnum];
   face[2] = indexMap[f->vertex[2]->vnum];
   ret.append( face );
 } while ( ++f != chullp->faces.end() );
 return ret;
} /* faces */

template<class VERTEX>
bool
CHull<VERTEX>::checkHull()
{
 return chullp->checks();
} /* checkHull */

#ifdef TESTING

#include "cdsVec3.h"
struct Vec3Vertex: public CDSVec3 {
  int index;

  Vec3Vertex(const CDSVec3& v) : CDSVec3(v), index(-1) {}

  typedef CDSVec3 Point3;
};

int testCHull()
{  
 int exit = 0;
 cout << "testing CHull...";
 
 

 CHull<Vec3Vertex> chull;

 chull.addVertex( CDSVec3( -1.181, -14.208,  -2.087) );
 chull.addVertex( CDSVec3( -1.414, -15.034,  -1.432) );
 chull.addVertex( CDSVec3(  0.322, -13.916,  -2.033) );
 chull.addVertex( CDSVec3(  0.637, -13.844,  -1.003) );
 chull.addVertex( CDSVec3(  0.525, -12.983,  -2.537) );
 chull.addVertex( CDSVec3(  1.091, -15.044,  -2.722) );
 chull.addVertex( CDSVec3(  0.787, -15.108,  -3.756) );
 chull.addVertex( CDSVec3(  0.880, -15.980,  -2.226) );
 chull.addVertex( CDSVec3(  2.866, -14.700,  -2.633) );
 chull.addVertex( CDSVec3(  3.365, -15.692,  -4.063) );
 chull.addVertex( CDSVec3(  4.443, -15.730,  -4.117) );
 chull.addVertex( CDSVec3(  2.975, -15.245,  -4.965) );
 chull.addVertex( CDSVec3(  2.975, -16.694,  -3.960) );
 chull.addVertex( CDSVec3( -1.940, -12.965,  -1.634) );
 chull.addVertex( CDSVec3( -2.215, -12.073,  -2.438) );
 chull.addVertex( CDSVec3( -1.578, -14.575,  -3.441) );
 chull.addVertex( CDSVec3( -1.350, -13.797,  -4.093) );
 chull.addVertex( CDSVec3( -2.601, -14.758,  -3.465) );
 chull.addVertex( CDSVec3( -1.066, -15.432,  -3.732) );
 chull.addVertex( CDSVec3( -2.232, -12.888,  -0.349) );
 chull.addVertex( CDSVec3( -2.000, -13.630,   0.247) );
 chull.addVertex( CDSVec3( -2.896, -11.708,   0.190) );
 chull.addVertex( CDSVec3( -3.499, -11.257,  -0.585) );
 chull.addVertex( CDSVec3( -3.797, -12.076,   1.376) );
 chull.addVertex( CDSVec3( -4.547, -12.782,   1.050) );
 chull.addVertex( CDSVec3( -4.280, -11.185,   1.748) );
 chull.addVertex( CDSVec3( -2.960, -12.704,   2.494) );
 chull.addVertex( CDSVec3( -2.425, -11.928,   3.022) );
 chull.addVertex( CDSVec3( -2.254, -13.400,   2.067) );
 chull.addVertex( CDSVec3( -3.878, -13.444,   3.469) );
 chull.addVertex( CDSVec3( -3.434, -14.293,   4.215) );
 chull.addVertex( CDSVec3( -5.150, -13.155,   3.493) );
 chull.addVertex( CDSVec3( -5.506, -12.463,   2.898) );
 chull.addVertex( CDSVec3( -5.749, -13.630,   4.106) );
 chull.addVertex( CDSVec3( -1.805, -10.740,   0.629) );
 chull.addVertex( CDSVec3( -0.959, -11.089,   1.455) );
 chull.addVertex( CDSVec3( -1.867,  -9.526,   0.123) );
 chull.addVertex( CDSVec3( -2.569,  -9.310,  -0.526) );
 chull.addVertex( CDSVec3( -0.908,  -8.498,   0.512) );
 chull.addVertex( CDSVec3( -0.113,  -8.954,   1.084) );
 chull.addVertex( CDSVec3( -0.310,  -7.814,  -0.721) );
 chull.addVertex( CDSVec3(  0.300,  -6.981,  -0.404) );
 chull.addVertex( CDSVec3( -1.110,  -7.454,  -1.350) );
 chull.addVertex( CDSVec3(  0.542,  -8.783,  -1.507) );
 chull.addVertex( CDSVec3( -0.049,  -9.590,  -2.488) );
 chull.addVertex( CDSVec3( -1.112,  -9.528,  -2.669) );
 chull.addVertex( CDSVec3(  1.921,  -8.867,  -1.270) );
 chull.addVertex( CDSVec3(  2.378,  -8.244,  -0.516) );
 chull.addVertex( CDSVec3(  0.737, -10.477,  -3.234) );
 chull.addVertex( CDSVec3(  0.282, -11.093,  -3.995) );
 chull.addVertex( CDSVec3(  2.707,  -9.759,  -2.011) );
 chull.addVertex( CDSVec3(  3.769,  -9.826,  -1.826) );
 chull.addVertex( CDSVec3(  2.114, -10.564,  -2.992) );
 chull.addVertex( CDSVec3(  2.891, -11.432,  -3.730) );
 chull.addVertex( CDSVec3(  3.767, -11.453,  -3.338) );
 chull.addVertex( CDSVec3( -1.658,  -7.479,   1.353) );
 chull.addVertex( CDSVec3( -2.859,  -7.282,   1.169) );
 chull.addVertex( CDSVec3( -0.948,  -6.835,   2.253) );
 chull.addVertex( CDSVec3(  0.011,  -7.026,   2.317) );
 chull.addVertex( CDSVec3( -1.523,  -5.850,   3.160) );
 chull.addVertex( CDSVec3( -2.584,  -5.766,   2.977) );
 chull.addVertex( CDSVec3( -1.287,  -6.298,   4.604) );
 chull.addVertex( CDSVec3( -0.225,  -6.409,   4.767) );
 chull.addVertex( CDSVec3( -1.776,  -7.248,   4.760) );
 chull.addVertex( CDSVec3( -1.839,  -5.287,   5.609) );
 chull.addVertex( CDSVec3( -2.872,  -5.075,   5.378) );
 chull.addVertex( CDSVec3( -1.264,  -4.374,   5.556) );
 chull.addVertex( CDSVec3( -1.740,  -5.873,   7.021) );
 chull.addVertex( CDSVec3( -1.850,  -5.080,   7.745) );
 chull.addVertex( CDSVec3( -0.775,  -6.342,   7.146) );
 chull.addVertex( CDSVec3( -2.841,  -6.914,   7.239) );
 chull.addVertex( CDSVec3( -2.724,  -7.718,   6.528) );
 chull.addVertex( CDSVec3( -3.808,  -6.452,   7.104) );
 chull.addVertex( CDSVec3( -2.736,  -7.454,   8.624) );
 chull.addVertex( CDSVec3( -1.793,  -7.869,   8.764) );
 chull.addVertex( CDSVec3( -3.461,  -8.186,   8.767) );
 chull.addVertex( CDSVec3( -2.882,  -6.684,   9.308) );
 chull.addVertex( CDSVec3( -0.851,  -4.504,   2.928) );
 chull.addVertex( CDSVec3(  0.351,  -4.424,   2.685) );
 chull.addVertex( CDSVec3( -1.651,  -3.463,   3.076) );
 chull.addVertex( CDSVec3( -2.603,  -3.617,   3.247) );
 chull.addVertex( CDSVec3( -1.160,  -2.090,   2.992) );
 chull.addVertex( CDSVec3( -0.092,  -2.110,   2.831) );
 chull.addVertex( CDSVec3( -1.822,  -1.314,   1.842) );
 chull.addVertex( CDSVec3( -2.895,  -1.370,   1.952) );

 chull.create();

 if ( !chull.checkHull() ) {
   cout << "CHull:checkHull failed\n";
   exit=1;
 }

 CDSList<CDSVec3> overtices;
 overtices.append( CDSVec3( 0.88, -15.98, -2.226 ) );
 overtices.append( CDSVec3( 4.443, -15.73, -4.117 ) );
 overtices.append( CDSVec3( 2.975, -15.245, -4.965 ) );
 overtices.append( CDSVec3( 2.975, -16.694, -3.96 ) );
 overtices.append( CDSVec3( -1.35, -13.797, -4.093 ) );
 overtices.append( CDSVec3( -2.601, -14.758, -3.465 ) );
 overtices.append( CDSVec3( -1.066, -15.432, -3.732 ) );
 overtices.append( CDSVec3( -3.499, -11.257, -0.585 ) );
 overtices.append( CDSVec3( -3.434, -14.293, 4.215 ) );
 overtices.append( CDSVec3( -5.506, -12.463, 2.898 ) );
 overtices.append( CDSVec3( -5.749, -13.63, 4.106 ) );
 overtices.append( CDSVec3( -1.112, -9.528, -2.669 ) );
 overtices.append( CDSVec3( 0.282, -11.093, -3.995 ) );
 overtices.append( CDSVec3( 3.769, -9.826, -1.826 ) );
 overtices.append( CDSVec3( 2.891, -11.432, -3.73 ) );
 overtices.append( CDSVec3( 3.767, -11.453, -3.338 ) );
 overtices.append( CDSVec3( -1.85, -5.08, 7.745 ) );
 overtices.append( CDSVec3( -0.775, -6.342, 7.146 ) );
 overtices.append( CDSVec3( -3.808, -6.452, 7.104 ) );
 overtices.append( CDSVec3( -1.793, -7.869, 8.764 ) );
 overtices.append( CDSVec3( -3.461, -8.186, 8.767 ) );
 overtices.append( CDSVec3( -2.882, -6.684, 9.308 ) );
 overtices.append( CDSVec3( -1.16, -2.09, 2.992 ) );
 overtices.append( CDSVec3( -0.092, -2.11, 2.831 ) );
 overtices.append( CDSVec3( -1.822, -1.314, 1.842 ) );
 overtices.append( CDSVec3( -2.895, -1.37, 1.952 ) );

 CDSList<Vec3Vertex> vertices = chull.vertices();
 for (int i=0 ; i<vertices.size() ; i++)
   if ( norm(CDSVec3(vertices[i]) - overtices[i]) >1e-5 )
     cout << "vertex " << i << " " << vertices[i] 
	  << " != " << overtices[i] << '\n';
 
 cout << (exit?"failed":"ok") << endl;
 return exit;
} /* test */

#endif /* TESTING */
