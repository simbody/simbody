
#ifndef __subMatrixDecl_h__
#define __subMatrixDecl_h__

#include <assert.h>
#include <cdsVector.h>
#include <cdsExcept.h>

#include <iostream>
#include <iomanip>

using namespace std;

template<class MATRIX>
class SubMatrix {
  MATRIX&  m;
  int offset1_;
  int offset2_;
  int size1;
  int size2;
public:
  typedef typename MATRIX::ElementType ElementType;
  SubMatrix(MATRIX&  m,
	    int offset1,
	    int offset2,
	    int size1,
	    int size2) : 
    m(m), offset1_(offset1), offset2_(offset2), size1(size1), size2(size2) 
  {
   assert(offset1>=m.offset1());
   assert(offset2>=m.offset2());
   assert(offset1+size1<=m.offset1()+m.nrow());
   assert(offset2+size2<=m.offset2()+m.ncol());
  }
  int nrow()  	const { return size1;}
  int ncol()  	const { return size2;}
  int offset1() const { return 0;}    //indexing into m is zero-offset
  int offset2() const { return 0;}    //indexing into m is zero-offset

  SubMatrix<MATRIX>&
  operator=(const SubMatrix<MATRIX>& mright) {
   assert(mright.nrow() == nrow());
   assert(mright.ncol() == ncol());
   for (int i=0 ; i<mright.nrow() ; i++)
     for (int j=0 ; j<mright.ncol() ; j++)
       m(offset1_+i,offset2_+j) = mright(mright.offset1()+i,
					 mright.offset2()+j);
   return *this;
  }
  template<class Matrix2>
  SubMatrix<MATRIX>&
  operator=(const Matrix2& mright) {
   assert(mright.nrow() == nrow());
   assert(mright.ncol() == ncol());
   for (int i=0 ; i<mright.nrow() ; i++)
     for (int j=0 ; j<mright.ncol() ; j++)
       m(offset1_+i,offset2_+j) = mright(mright.offset1()+i,
					 mright.offset2()+j);
   return *this;
  }

  //assignment from Vector
  template<class VECTOR>
  SubMatrix<MATRIX>&
  assignFromVector(const VECTOR& vector);

  CDSVector<ElementType,0> convertToVector();
  
  // zero-offset index operator
  const ElementType& operator()(int i,
				int j) const {       
   assert(i>=0 && i<size1);
   assert(j>=0 && j<size2);
   return m(offset1_+i,offset2_+j);
  }
  ElementType& operator()(int i,        //non-const version
			  int j) {       
   assert(i>=0 && i<size1);
   assert(j>=0 && j<size2);
   return m(offset1_+i,offset2_+j);
  }
  template<class Matrix2>           //operator+=
  SubMatrix<MATRIX>&
  operator+=(const Matrix2& mr) {
   assert(mr.nrow() == nrow());
   assert(mr.ncol() == ncol());
   for (int i=0 ; i<mr.nrow() ; i++)
     for (int j=0 ; j<mr.ncol() ; j++)
       m(offset1_+i,offset2_+j) += mr(mr.offset1()+i,mr.offset2()+j);
   return *this;
  }

};

template<class MATRIX>
ostream& operator<<(      ostream&           s,
		    const SubMatrix<MATRIX>& m)
{
 int width = s.width(); // apply width to numeric fields only
 s << "{ ";
 for (int i=0 ; i<m.nrow() ; i++) {
   s << "{ ";
   s << setw(width) << m(i,0); //will fail for vectors of size 0!
   for (int j=1 ; j<m.ncol() ; j++) 
     s << ", " << setw(width) << m(i,j) ;
   s << " }";
   if (i<m.nrow()-1) s << ", ";
 }
 s << " }";
 return s;
} /* operator<< */

//only for square, FixedMatrix
template<class MatrixType, class MATRIX> 
MatrixType
operator*(const MatrixType&         m1,
	  const SubMatrix<MATRIX>&  m2);

//only for square, FixedMatrix
template<class MatrixType, class MATRIX>
MatrixType
operator*(const SubMatrix<MATRIX>&  m1,
	  const MatrixType&         m2);

template<class MatrixType, class MATRIX>
MatrixType
operator+(const MatrixType&         m1,
	  const SubMatrix<MATRIX>&  m2);

template<class MatrixType, class MATRIX>
MatrixType
operator+(const SubMatrix<MATRIX>&  m1,
	  const MatrixType&         m2)
{ return m2+m1; }

template<class MatrixType, class MATRIX>
MatrixType
operator-(const MatrixType&         m1,
	  const SubMatrix<MATRIX>&  m2);

template<class MatrixType, class MATRIX>
MatrixType
operator-(const SubMatrix<MATRIX>&  m1,
	  const MatrixType&         m2);

int SubMatrix_test();

#endif /* __subMatrixDecl_h__ */
