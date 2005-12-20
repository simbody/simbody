#ifndef __subMatrixDefs_h__
#define __subMatrixDefs_h__

template<class MATRIX>
template<class VECTOR>
SubMatrix<MATRIX>&
SubMatrix<MATRIX>::assignFromVector(const VECTOR& vector)
{
 if ( rows()==1 && cols()==vector.size() )
   for (int i=0 ; i<cols() ; i++)
     (*this)(0,i) = vector(i);
 else if ( cols()==1 && rows()==vector.size() )
   for (int i=0 ; i<rows() ; i++)
     (*this)(i,0) = vector(i);
 else
   throw CDS::exception("SubMatrix::assignFromVector: invalid dimensions.");
 return *this;
} /* assignFromVector */

template<class MATRIX>
CDSVector<typename MATRIX::ElementType>
SubMatrix<MATRIX>::convertToVector()
{
 CDSVector<ElementType> ret;
 if ( rows()==1 ) {
   ret.resize(cols());
   for (int i=0 ; i<cols() ; i++)
     ret(i) = (*this)(0,i);
 } else if ( cols()==1 ) {
   ret.resize(rows());
   for (int i=0 ; i<rows() ; i++)
     ret(i) = (*this)(i,0);
 } else
   throw CDS::exception("SubMatrix::convertToVector: invalid dimensions.");
 return ret;
} /* convertToVector */

template<class MatrixType, class MATRIX>
MatrixType
operator*(const MatrixType&         m1,
	  const SubMatrix<MATRIX>&  m2)
{
 assert( m1.cols() == m2.rows() );

 MatrixType r((typename MATRIX::ElementType)0); //works only for square m1,m2!!
 for (int i=0 ; i<m1.rows() ; i++)
   for (int k=0 ; k<m2.cols() ; k++)
     for (int j=0 ; j<m1.cols() ; j++) 
       r(i,k) += m1(i,j) * m2(j,k);

 return r;
} /* operator* (matrix,matrix) */

template<class MatrixType, class MATRIX>
MatrixType
operator*(const SubMatrix<MATRIX>&  m1,
	  const MatrixType&         m2)
{
 assert( m1.cols() == m2.rows() );

 MatrixType r((typename MATRIX::ElementType)0); //works only for square m1,m2!!
 for (int i=0 ; i<m1.rows() ; i++)
   for (int k=0 ; k<m2.cols() ; k++)
     for (int j=0 ; j<m1.cols() ; j++) 
       r(i,k) += m1(i,j) * m2(j,k);

 return r;
} /* operator* (matrix,matrix) */

template<class MatrixType, class MATRIX>
MatrixType
operator+(const MatrixType&         m1,
	  const SubMatrix<MATRIX>&  m2)
{
 assert( m1.cols() == m2.cols() );
 assert( m1.rows() == m2.rows() );
 MatrixType r = m1;
 for (int i=0 ; i<m1.rows() ; i++)
   for (int j=0 ; j<m1.cols() ; j++) 
     r(i,j) += m2(i,j);

 return r;
} /* operator* (matrix,matrix) */


template<class MatrixType, class MATRIX>
MatrixType
operator-(const MatrixType&         m1,
	  const SubMatrix<MATRIX>&  m2)
{
 assert( m1.cols() == m2.cols() );
 assert( m1.rows() == m2.rows() );
 MatrixType r = m1;
 for (int i=0 ; i<m1.rows() ; i++)
   for (int j=0 ; j<m1.cols() ; j++) 
     r(i,j) -= m2(i,j);

 return r;
} /* operator* (matrix,matrix) */

template<class MatrixType, class MATRIX>
MatrixType
operator-(const SubMatrix<MATRIX>&  m1,
	  const MatrixType&         m2)
{
 assert( m1.cols() == m2.cols() );
 assert( m1.rows() == m2.rows() );
 MatrixType r = m2;
 r.scale(-1);
 for (int i=0 ; i<m1.rows() ; i++)
   for (int j=0 ; j<m1.cols() ; j++) 
     r(i,j) += m1(i,j);

 return r;
} /* operator* (matrix,matrix) */


#endif /* __subMatrixDefs_h__ */
