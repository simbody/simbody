#ifndef _FACTOR_REP_H_ 
#define _FACTOR_REP_H_
/* Portions copyright (c) 2007 Stanford University and Jack Middleton.
 * Contributors:
 *
 * Permission is hereby granted, free of charge, to any person obtaining
 * a copy of this software and associated documentation files (the
 * "Software"), to deal in the Software without restriction, including
 * without limitation the rights to use, copy, modify, merge, publish,
 * distribute, sublicense, and/or sell copies of the Software, and to
 * permit persons to whom the Software is furnished to do so, subject
 * to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included
 * in all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
 * OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
 * MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
 * IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY
 * CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
 * TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE
 * SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
 */
#include "SimTKmath.h"
#include "WorkSpace.h"

namespace SimTK {

class FactorLURepBase {
    public:

    virtual ~FactorLURepBase(){};
    virtual FactorLURepBase* clone() const {};

   virtual void solve( const Vector_<float>& b, Vector_<float>& x ) const {
       checkIfFactored("solve");
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","solve",
       "solve called with rhs of type <float>  which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<double>& b, Vector_<double>& x ) const {
       checkIfFactored("solve");
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","solve",
       " solve called with rhs of type <double>  which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<std::complex<float> >& b, Vector_<std::complex<float> >& x ) const {
       checkIfFactored("solve");
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","solve",
       " solve called with rhs of type complex<float> which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<std::complex<double> >& b, Vector_<std::complex<double> >& x ) const {
       checkIfFactored("solve");
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","solve",
       " solve called with rhs of type complex<double>  which does not match type of original linear system \n");   
   }
     virtual void solve( const Matrix_<float>& b, Matrix_<float>& x ) const {
       checkIfFactored("solve");
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","solve",
       " solve called with rhs of type <float>  which does not match type of original linear system \n");
   }
   virtual void solve( const Matrix_<double>& b, Matrix_<double>& x ) const {
       checkIfFactored("solve");
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","solve",
       " solve called with rhs of type <double>  which does not match type of original linear system \n");
   }
   virtual void solve( const Matrix_<std::complex<float> >& b, Matrix_<std::complex<float> >& x ) const {
       checkIfFactored("solve");
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","solve",
       " solve called with rhs of type complex<float> which does not match type of original linear system \n");
   }
   virtual void solve( const Matrix_<std::complex<double> >& b, Matrix_<std::complex<double> >& x ) const {
       checkIfFactored("solve");
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","solve",
       " solve called with rhs of type complex<double>  which does not match type of original linear system \n");   
   }
     virtual void getL( Matrix_<float>& l) const{
       checkIfFactored( "getL" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getL",
       " getL called with L of type <float> which does not match type of original linear system \n");
   }
   virtual void getL( Matrix_<double>& l ) const{
       checkIfFactored( "getL" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getL",
       " getL called with L of type <double>  which does not match type of original linear system \n");
   }
   virtual void getL( Matrix_<std::complex<float> >& l) const{
       checkIfFactored( "getL" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getL",
       " getL called with L of type complex<float> which does not match type of original linear system \n");
   }
   virtual void getL( Matrix_<std::complex<double> >& l ) const{
       checkIfFactored( "getL" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getL",
       " getL called with L of type complex<double>  which does not match type of original linear system \n");   
   }
   virtual void getU( Matrix_<float>& u) const{
       checkIfFactored( "getU" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getU",
       " getU called with U of type <float> which does not match type of original linear system \n");
   }
   virtual void getU( Matrix_<double>& u ) const{
       checkIfFactored( "getU" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getU",
       " getU called with U of type <double>  which does not match type of original linear system \n");
   }
   virtual void getU( Matrix_<std::complex<float> >& u) const{
       checkIfFactored( "getU" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getU",
       " getU called with U of type complex<float> which does not match type of original linear system \n");
   }
   virtual void getU( Matrix_<std::complex<double> >& u ) const{
       checkIfFactored( "getU" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getU",
       " getU called with U of type complex<double>  which does not match type of original linear system \n");   
   }
   virtual void getD( Matrix_<float>& d) const{
       checkIfFactored( "getD" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getD",
       " getD called with D of type <float> which does not match type of original linear system \n");
   }
   virtual void getD( Matrix_<double>& d ) const{
       checkIfFactored( "getD" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getD",
       " getD called with D of type <double>  which does not match type of original linear system \n");
   }
   virtual void getD( Matrix_<std::complex<float> >& d) const{
       checkIfFactored( "getD" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getD",
       " getD called with D of type complex<float> which does not match type of original linear system \n");
   }
   virtual void getD( Matrix_<std::complex<double> >& d ) const{
       checkIfFactored( "getD" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getD",
       " getD called with D of type complex<double>  which does not match type of original linear system \n");   
   }
    virtual void inverse(  Matrix_<double>& inverse ) const{
        SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","inverse",         "inverse(  <double> ) called with type that is inconsistant with the original matrix  \n");
    }
    virtual void inverse(  Matrix_<float>& inverse ) const{
        SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","inverse",         "inverse(  <float> ) called with type that is inconsistant with the original matrix  \n");
    }
    virtual void inverse(  Matrix_<std::complex<float> >& inverse ) const{
        SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","inverse",         "inverse(  std::complex<float> ) called with type that is inconsistant with the original matrix  \n");
    }
    virtual void inverse(  Matrix_<std::complex<double> >& inverse ) const{
        SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","inverse",
        "inverse(  std::complex<double> ) called with type that is inconsistant with the original matrix  \n");
    }

   virtual bool isSingular() const{ return false;};
   virtual int getSingularIndex() const{ return 1; };
   virtual  Real getConditionNumber() const{ return 0.0;};
   
   virtual void getErrorBounds( const Vector_<float>& err, Vector_<float>& berr ){
       checkIfFactored( "getErrorBounds" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getErrorBounds",
       " getErrorBounds called with rhs arguments type <float>  which does not match type of original linear system \n");
   }
   virtual void getErrorBounds( const Vector_<double>& err, Vector_<float>& berr ){
       checkIfFactored( "getErrorBounds" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getErrorBounds",
       " getErrorBounds called with arguments of type <double>  which does not match type of original linear system \n");
   }
   virtual void getErrorBounds( const Vector_<std::complex<float> >& err, Vector_<float>& berr ){
       checkIfFactored( "getErrorBounds" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getErrorBounds",
       " getErrorBounds called with arguments of type complex<float> which does not match type of original linear system \n");
   }
   virtual void getErrorBounds( const Vector_<std::complex<double> >& err, Vector_<float>& berr ){
       checkIfFactored( "getErrorBounds" );
       SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","getErrorBounds",
       " getErrorBounds called with arguments of type complex<double>  which does not match type of original linear system \n");   
   } 
   bool isFactored;

   private:
   void checkIfFactored( const char* function_name)  const {
       if( !isFactored ) {
           SimTK_APIARGCHECK_ALWAYS(false,"FactorLU",function_name,
           "matrix has not been factored \n");
       }

       return;
   }



}; // class FactorLURepBase

class FactorLUDefault : public FactorLURepBase {
   public:
   FactorLUDefault();
   FactorLURepBase* clone() const;

};

template <typename T>
class FactorLURep : public FactorLURepBase {
   public:
   template <class ELT> FactorLURep( const Matrix_<ELT>&  );
   FactorLURep();

   ~FactorLURep();
   FactorLURepBase* clone() const;

   template < class ELT > void factor(const Matrix_<ELT>& ); 
   void solve( const Vector_<T>& b, Vector_<T>& x ) const;
   void solve( const Matrix_<T>& b, Matrix_<T>& x ) const;
   void inverse( Matrix_<T>& m ) const;

   void  getL( Matrix_<T>& l ) const; 
   void  getU( Matrix_<T>& u ) const;
   void  getD( Matrix_<T>& d ) const;
   void getErrorBounds( Vector_<T>& err, Vector_<T>& berr) const;
   Real getConditionNumber() const;
   bool isSingular() const;
   int getSingularIndex() const;
 
   private:

// factored matrix stored in LAPACK LU format
   template < class ELT> int getType(ELT*);   
   bool isLUinitialized;
   bool positiveDefinite;
   int nRow;
   int nCol;
   int mn;        // min(m,n)
   int LUtype;
   int singularIndex;
   int elementSize;
   int imagOffset;
   MatrixStructures::Structure structure;
   MatrixConditions::Condition condition; 
   MatrixShapes::Shape shape;  
   MatrixSparseFormats::Sparsity sparsity;
   MatrixStorageFormats::Storage storage;


   TypedWorkSpace<int>  pivots;
   TypedWorkSpace<T>    lu;
}; // end class FactorLURep
} // namespace SimTK
#endif   //  _FACTOR_REP_H_
