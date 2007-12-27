#ifndef _EIGEN_REP_H_ 
#define _EIGEN_REP_H_
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

class EigenRepBase {
    public:

    virtual ~EigenRepBase(){};

   virtual void getValues( Vector_<float>& values ){
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getValues",
       "getValues called with vector of type <float>   \n");
   }
   virtual void getValues( Vector_<double>& values ){
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getValues",
       "  getValues called with illegal type <double>   \n");
   }
   virtual void getValues( Vector_<std::complex<float> >& values ){
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getValues",
       "  getValues called with illegal type complex<float>  \n");
   }
   virtual void getValues( Vector_<std::complex<double> >& values ){
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getValues",
       "  getValues called with illegal type complex<double>   \n");   
   }
   virtual void getVectors( Matrix_<float>& vectors ){
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getVectors",
       "getVectors called with vector of type <float>   \n");
   }
   virtual void getVectors( Matrix_<double>& vectors ){
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getVectors",
       "  getVectors called with illegal type <double>   \n");
   }
   virtual void getVectors( Matrix_<std::complex<float> >& vectors ){
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getVectors",
       "  getVectors called with illegal type complex<float>  \n");
   }
   virtual void getVectors( Matrix_<std::complex<double> >& vectors ){
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getVectors",
       "  getVectors called with illegal type complex<double>   \n");   
   }

   virtual void selectValues( const std::vector<bool>& selectedValues) { return; }
   virtual void selectVectors( const std::vector<bool>& selectedVectors) { return;  }

}; // class EigenRepBase

template <typename T>
class EigenRep : public EigenRepBase {
   public:
   template <class ELT> EigenRep( const Matrix_<ELT>&  );

   ~EigenRep();

   void getValues( Vector_<std::complex< typename CNT<T>::TReal> >& values );
   void getVectors( Matrix_<std::complex< typename CNT<T>::TReal> >& vectors );
   void getValues( Vector_< typename CNT<T>::TReal >& values );
   void getVectors( Matrix_< typename CNT<T>::TReal >& vectors );
   void selectValues( const std::vector<bool>& selectedValues);
   void selectVectors( const std::vector<bool>& selectedVectors);
   void computeValues();
 
   private:
  
   int n;
   bool needAllValues;
   bool needAllVectors;
   bool computedValues;     // true if eigen values  already computed
   bool computedVectors;    // true if eigen vectors already computed
   TypedWorkSpace<bool> selectedValues;
   TypedWorkSpace<bool> selectedVectors;
   TypedWorkSpace<T> inputMatrix;
   TypedWorkSpace<std::complex< typename CNT<T>::TReal > >eigenValues;
   Matrix_<std::complex< typename CNT<T>::TReal > > rightVectors;

}; // end class EigenRep
} // namespace SimTK
#endif   //  _EIGEN_REP_H_
