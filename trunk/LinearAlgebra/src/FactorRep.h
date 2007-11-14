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

// TODO change printf's to throw exceptions 
   virtual void solve( const Vector_<float>& b, Vector_<float>& x ){
       printf(" solve called with rhs of type <float>  which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<double>& b, Vector_<double>& x ){
       printf(" solve called with rhs of type <double>  which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<std::complex<float> >& b, Vector_<std::complex<float> >& x ){
       printf(" solve called with rhs of type complex<float> which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<std::complex<double> >& b, Vector_<std::complex<double> >& x ){
       printf(" solve called with rhs of type complex<double>  which does not match type of original linear system \n");   
   }
     virtual void solve( const Matrix_<float>& b, Matrix_<float>& x ){
       printf(" solve called with rhs of type <float>  which does not match type of original linear system \n");
   }
   virtual void solve( const Matrix_<double>& b, Matrix_<double>& x ){
       printf(" solve called with rhs of type <double>  which does not match type of original linear system \n");
   }
   virtual void solve( const Matrix_<std::complex<float> >& b, Matrix_<std::complex<float> >& x ){
       printf(" solve called with rhs of type complex<float> which does not match type of original linear system \n");
   }
   virtual void solve( const Matrix_<std::complex<double> >& b, Matrix_<std::complex<double> >& x ){
       printf(" solve called with rhs of type complex<double>  which does not match type of original linear system \n");   
   }
     virtual void getL( Matrix_<float>& l) const{
       printf(" getL called with L of type <float> which does not match type of original linear system \n");
   }
   virtual void getL( Matrix_<double>& l ) const{
       printf(" getL called with L of type <double>  which does not match type of original linear system \n");
   }
   virtual void getL( Matrix_<std::complex<float> >& l) const{
       printf(" getL called with L of type complex<float> which does not match type of original linear system \n");
   }
   virtual void getL( Matrix_<std::complex<double> >& l ) const{
       printf(" getL called with L of type complex<double>  which does not match type of original linear system \n");   
   }
   virtual void getU( Matrix_<float>& u) const{
       printf(" getU called with U of type <float> which does not match type of original linear system \n");
   }
   virtual void getU( Matrix_<double>& u ) const{
       printf(" getU called with U of type <double>  which does not match type of original linear system \n");
   }
   virtual void getU( Matrix_<std::complex<float> >& u) const{
       printf(" getU called with U of type complex<float> which does not match type of original linear system \n");
   }
   virtual void getU( Matrix_<std::complex<double> >& u ) const{
       printf(" getU called with U of type complex<double>  which does not match type of original linear system \n");   
   }
   virtual void getD( Matrix_<float>& d) const{
       printf(" getD called with D of type <float> which does not match type of original linear system \n");
   }
   virtual void getD( Matrix_<double>& d ) const{
       printf(" getD called with D of type <double>  which does not match type of original linear system \n");
   }
   virtual void getD( Matrix_<std::complex<float> >& d) const{
       printf(" getD called with D of type complex<float> which does not match type of original linear system \n");
   }
   virtual void getD( Matrix_<std::complex<double> >& d ) const{
       printf(" getD called with D of type complex<double>  which does not match type of original linear system \n");   
   }
   virtual bool isSingular() const{ return false;};
   virtual int getSingularIndex() const{ return 1; };
   virtual void display(int){};
   virtual  Real getConditionNumber() const{ return 0.0;};
   
   virtual void getErrorBounds( const Vector_<float>& err, Vector_<float>& berr ){
       printf(" getErrorBounds called with rhs arguments type <float>  which does not match type of original linear system \n");
   }
   virtual void getErrorBounds( const Vector_<double>& err, Vector_<float>& berr ){
       printf(" getErrorBounds called with arguments of type <double>  which does not match type of original linear system \n");
   }
   virtual void getErrorBounds( const Vector_<std::complex<float> >& err, Vector_<float>& berr ){
       printf(" getErrorBounds called with arguments of type complex<float> which does not match type of original linear system \n");
   }
   virtual void getErrorBounds( const Vector_<std::complex<double> >& err, Vector_<float>& berr ){
       printf(" getErrorBounds called with arguments of type complex<double>  which does not match type of original linear system \n");   
   }


}; // class FactorLURepBase

template <typename T>
class FactorLURep : public FactorLURepBase {
   public:
   template <class ELT> FactorLURep( const Matrix_<ELT>&  );

   ~FactorLURep();

   template < class ELT > void factor(const Matrix_<ELT>& ); 
   template < class ELT > void solve( const Vector_<ELT>& b, Vector_<ELT>& x );
   template < class ELT > void solve( const Matrix_<ELT>& b, Matrix_<ELT>& x );

   template < class ELT > void  getL( Matrix_<ELT>& l ) const; 
   template < class ELT > void  getU( Matrix_<ELT>& u ) const;
   template < class ELT > void  getD( Matrix_<ELT>& d ) const;
   template < class ELT > void getErrorBounds( Vector_<ELT>& err, Vector_<ELT>& berr) const;
   Real getConditionNumber() const;
   bool isSingular() const;
   void display(int);
   int getSingularIndex() const;
   void copyElement( const int i, const int j, T* ptr );
   template < class ELT > void copyElement( const int i, const int j, std::complex<ELT>* ptr );
 
   private:

// factored matrix stored in LAPACK LU format
   template < class ELT> int getType(ELT*);   
   void printElement( int i, int j);   
   template < class ELT> void initLU(const Matrix_<ELT>&);   
   template < class ELT> void initLU(const Matrix_<negator<ELT> >&);   
   template < class ELT> void initLU(const Matrix_<std::complex<ELT> >&);   
   template < class ELT> void initLU(const Matrix_<negator<std::complex<ELT> > >&);   
   template < class ELT> void initLU(const Matrix_<conjugate<ELT> >&);   
   template < class ELT> void initLU(const Matrix_<negator<conjugate<ELT> > >&);   
   bool isLUinitialized;
   int nRow;
   int nCol;
   int LUtype;
   int singularIndex;
   int elementSize;
   int imagOffset;
   TypedWorkSpace<int>  pivots;
   TypedWorkSpace<T>    lu;
}; // end class FactorLURep
} // namespace SimTK
#endif   //  _FACTOR_REP_H_
