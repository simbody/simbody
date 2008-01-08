#ifndef _FACTOR_QTZ_REP_H_ 
#define _FACTOR_QTZ_REP_H_
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

class FactorQTZRepBase {
    public:

    virtual ~FactorQTZRepBase(){};

   virtual void solve( const Vector_<float>& b, Vector_<float>& x ){
       SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
       "solve called with rhs of type <float>  which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<double>& b, Vector_<double>& x ){
       SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
       " solve called with rhs of type <double>  which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<std::complex<float> >& b, Vector_<std::complex<float> >& x ){
       SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
       " solve called with rhs of type complex<float> which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<std::complex<double> >& b, Vector_<std::complex<double> >& x ){
       SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
       " solve called with rhs of type complex<double>  which does not match type of original linear system \n");   
   }
     virtual void solve( const Matrix_<float>& b, Matrix_<float>& x ){
       SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
       " solve called with rhs of type <float>  which does not match type of original linear system \n");
   }
   virtual void solve( const Matrix_<double>& b, Matrix_<double>& x ){
       SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
       " solve called with rhs of type <double>  which does not match type of original linear system \n");
   }
   virtual void solve( const Matrix_<std::complex<float> >& b, Matrix_<std::complex<float> >& x ){
       SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
       " solve called with rhs of type complex<float> which does not match type of original linear system \n");
   }
   virtual void solve( const Matrix_<std::complex<double> >& b, Matrix_<std::complex<double> >& x ){
       SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
       " solve called with rhs of type complex<double>  which does not match type of original linear system \n");   
   }


}; // class FactorQTZRepBase

template <typename T>
class FactorQTZRep : public FactorQTZRepBase {
   public:
   template <class ELT> FactorQTZRep( const Matrix_<ELT>&  );
   template <class ELT> FactorQTZRep( const Matrix_<ELT>&, typename CNT<T>::TReal  );

   ~FactorQTZRep();

   template < class ELT > void factor(const Matrix_<ELT>& ); 
   void solve( const Vector_<T>& b, Vector_<T>& x );
   void solve( const Matrix_<T>& b, Matrix_<T>& x );
 
   private:
  
   void doSolve( Matrix_<T>& b, Matrix_<T>& x );

   int mn;           // min of number of rows or columns
   int maxmn;        // max of number of rows or columns
   int nRow;
   int nCol;
   int rank;
   bool scaleLinSys;
   bool scaleRHS;
   typename CNT<T>::TReal linSysScaleF;
   typename CNT<T>::TReal rhsScaleF;
   typename CNT<T>::TReal anrm;
   typename CNT<T>::TReal rcond;
   TypedWorkSpace<int>    pivots;
   TypedWorkSpace<T>      qtz;
   TypedWorkSpace<T>      tauGEQP3;
   TypedWorkSpace<T>      tauORMQR;

}; // end class FactorQTZRep
} // namespace SimTK
#endif   //  _FACTOR_QTZ_REP_H_
