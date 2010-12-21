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

    virtual EigenRepBase* clone() const { return 0; };

   virtual void getAllEigenValuesAndVectors( Vector_<float>& values, Matrix_<float>& vectors ){
       checkIfFactored( "getAllEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValuesAndVectors",
       "getAllEigenValuesAndVectors called with vector of type <float, float>   \n");
   }
   virtual void getAllEigenValuesAndVectors( Vector_<double>& values, Matrix_<double>& vectors ){
       checkIfFactored( "getAllEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValuesAndVectors",
       "getAllEigenValuesAndVectors called with vector of type <double, double>   \n");
   }
   virtual void getAllEigenValuesAndVectors( Vector_<double>& values, Matrix_<std::complex<double> >& vectors ){
       checkIfFactored( "getAllEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValuesAndVectors",
       "getAllEigenValuesAndVectors called with vector of types  <double, std::complex<double>   \n");
   }
   virtual void getAllEigenValuesAndVectors( Vector_<std::complex<double> >& values, Matrix_<std::complex<double> >& vectors ){
       checkIfFactored( "getAllEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValuesAndVectors",
       "getAllEigenValuesAndVectors called with vector of type <std::complex<double>, std::complex<double> >  \n");
   }
   virtual void getAllEigenValuesAndVectors( Vector_<float>& values, Matrix_<std::complex<float> >& vectors ){
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValuesAndVectors",
       "getAllEigenValuesAndVectors called with vector of type <float, std::complex<float> >   \n");
   }
   virtual void getAllEigenValuesAndVectors( Vector_<std::complex<float> >& values, Matrix_<std::complex<float> >& vectors ){
       checkIfFactored( "getAllEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValuesAndVectors",
       "getAllEigenValuesAndVectors called with vector of type <std::complex<float, std::complex<float> >   \n");
   }

   virtual void getAllEigenValues( Vector_<float>& values ){
       checkIfFactored( "getAllEigenValues" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValues",
       "getAllEigenValues called with vector of type <float>   \n");
   }
   virtual void getAllEigenValues( Vector_<double>& values ){
       checkIfFactored( "getAllEigenValues" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValues",
       "getAllEigenValues called with vector of type <double>   \n");
   }
   virtual void getAllEigenValues( Vector_<std::complex<double> >& values){
       checkIfFactored( "getAllEigenValues" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValues",
       "getAllEigenValues called with vector of type <std::complex<double>   \n");
   }
   virtual void getAllEigenValues( Vector_<std::complex<float> >& values ){
       checkIfFactored( "getAllEigenValues" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getAllEigenValues",
       "getAllEigenValues called with vector of type <std::complex<float>   \n");
   }
   virtual void getFewEigenValuesAndVectors( Vector_<float>& values, Matrix_<float>& vectors, int ilow, int ihi ){
       checkIfFactored( "getFewEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
       "getFewEigenValuesAndVectors called with vector of type <float,float>   \n");
   }
   virtual void getFewEigenValuesAndVectors( Vector_<double>& values, Matrix_<double>& vectors, int ilow, int ihi ){
       checkIfFactored( "getFewEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
       "getFewEigenValuesAndVectors called with vector of type <double,double>   \n");
   }
   virtual void getFewEigenValuesAndVectors( Vector_<std::complex<float> >& values, Matrix_<std::complex<float> >& vectors, int ilow, int ihi ){
       checkIfFactored( "getFewEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
       "getFewEigenValuesAndVectors called with vector of type <std::complex<float>, std::complex<float> >   \n");
   }
   virtual void getFewEigenValuesAndVectors( Vector_<float>& values, Matrix_<std::complex<float> >& vectors, int ilow, int ihi ){
       checkIfFactored( "getFewEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
       "getFewEigenValuesAndVectors called with vector of type <float, std::complex<float> >   \n");
   }
   virtual void getFewEigenValuesAndVectors( Vector_<std::complex<double> >& values, Matrix_<std::complex<double> >& vectors, int ilow, int ihi ){
       checkIfFactored( "getFewEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
       "getFewEigenValuesAndVectors called with vector of type <std::complex<double>, std::complex<double> >   \n");
   }
   virtual void getFewEigenValuesAndVectors( Vector_<double>& values, Matrix_<std::complex<double> >& vectors, int ilow, int ihi ){
       checkIfFactored( "getFewEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
       "getFewEigenValuesAndVectors called with vector of type < double, std::complex<double> >   \n");
   }
   virtual void getFewEigenValuesAndVectors( Vector_<float>& values, Matrix_<float>& vectors, float rlow, float rhi ){
       checkIfFactored( "getFewEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
       "getFewEigenValuesAndVectors called with vector of type <float, float>   \n");
   }
   virtual void getFewEigenValuesAndVectors( Vector_<double>& values, Matrix_<double>& vectors, double rlow, double rhi ){
       checkIfFactored( "getFewEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
       "getFewEigenValuesAndVectors called with vector of type <double, double>   \n");
   }
   virtual void getFewEigenValuesAndVectors( Vector_<std::complex<float> >& values, Matrix_<std::complex<float> >& vectors, float rlow, float rhi ){
       checkIfFactored( "getFewEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
       "getFewEigenValuesAndVectors called with vector of type <std::complex<float>, std::complex<float> >   \n");
   }
   virtual void getFewEigenValuesAndVectors( Vector_<float>& values, Matrix_<std::complex<float> >& vectors, float rlow, float rhi ){
       checkIfFactored( "getFewEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
       "getFewEigenValuesAndVectors called with vector of type <float, std::complex<float> >   \n");
   }
   virtual void getFewEigenValuesAndVectors( Vector_<std::complex<double> >& values, Matrix_<std::complex<double> >& vectors, double rlow, double rhi ){
       checkIfFactored( "getFewEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
       "getFewEigenValuesAndVectors called with vector of type <std::complex<double>, std::complex<double> >   \n");
   }
   virtual void getFewEigenValuesAndVectors( Vector_<double>& values, Matrix_<std::complex<double> >& vectors, double rlow, double rhi ){
       checkIfFactored( "getFewEigenValuesAndVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValuesAndVectors",
       "getFewEigenValuesAndVectors called with vector of type <double, std::complex<double> >   \n");
   }
   virtual void getFewEigenValues( Vector_<float>& values, int ilow, int ihi ){
       checkIfFactored( "getFewEigenValues" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValues",
       "getFewEigenValues called with vector of type <float>   \n");
   }
   virtual void getFewEigenValues( Vector_<double>& values, int ilow, int ihi ){
       checkIfFactored( "getFewEigenValues" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValues",
       "getFewEigenValues called with vector of type <double>   \n");
   }
   virtual void getFewEigenValues( Vector_<std::complex<float> >& values, int ilow, int ihi ){
       checkIfFactored( "getFewEigenValues" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValues",
       "getFewEigenValues called with vector of type <std::complex<float> >   \n");
   }
   virtual void getFewEigenValues( Vector_<std::complex<double> >& values, int ilow, int ihi ){
       checkIfFactored( "getFewEigenValues" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValues",
       "getFewEigenValues called with vector of type <std::complex<double> >   \n");
   }
   virtual void getFewEigenValues( Vector_<float>& values, float rlow, float rhi ){
       checkIfFactored( "getFewEigenValues" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValues",
       "getFewEigenValues called with vector of type <float>   \n");
   }
   virtual void getFewEigenValues( Vector_<double>& values, double rlow, double rhi ){
       checkIfFactored( "getFewEigenValues" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValues",
       "getFewEigenValues called with vector of type <double>   \n");
   }
   virtual void getFewEigenValues( Vector_<std::complex<float> >& values, float rlow, float rhi ){
       checkIfFactored( "getFewEigenValues" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValues",
       "getFewEigenValues called with vector of type <std::complex<float> >   \n");
   }
   virtual void getFewEigenValues( Vector_<std::complex<double> >& values, double rlow, double rhi ){
       checkIfFactored( "getFewEigenValues" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenValues",
       "getFewEigenValues called with vector of type <std::complex<double> >   \n");
   }
   virtual void getFewEigenVectors( Matrix_<float>& vectors, int ilow, int ihi ){
       checkIfFactored( "getFewEigenVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenVectors",
       "getFewEigenVectors called with vector of type <float>   \n");
   }
   virtual void getFewEigenVectors(  Matrix_<double>& vectors, int ilow, int ihi ){
       checkIfFactored( "getFewEigenVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenVectors",
       "getFewEigenVectors called with vector of type <double>   \n");
   }
   virtual void getFewEigenVectors( Matrix_<std::complex<float> >& vectors, int ilow, int ihi ){
       checkIfFactored( "getFewEigenVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenVectors",
       "getFewEigenVectors called with vector of type <std::complex<float> >   \n");
   }
   virtual void getFewEigenVectors(  Matrix_<std::complex<double> >& vectors, int ilow, int ihi ){
       checkIfFactored( "getFewEigenVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenVectors",
       "getFewEigenVectors called with vector of type <std::complex<double> >   \n");
   }
   virtual void getFewEigenVectors(  Matrix_<float>& vectors, float rlow, float rhi ){
       checkIfFactored( "getFewEigenVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenVectors",
       "getFewEigenVectors called with vector of type <float>   \n");
   }
   virtual void getFewEigenVectors(  Matrix_<double>& vectors, double rlow, double rhi ){
       checkIfFactored( "getFewEigenVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenVectors",
       "getFewEigenVectors called with vector of type <double>   \n");
   }
   virtual void getFewEigenVectors( Matrix_<std::complex<float> >& vectors, float rlow, float rhi ){
       checkIfFactored( "getFewEigenVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenVectors",
       "getFewEigenVectors called with vector of type <std::complex<float> >   \n");
   }
   virtual void getFewEigenVectors( Matrix_<std::complex<double> >& vectors, double rlow, double rhi ){
       checkIfFactored( "getFewEigenVectors" );
       SimTK_APIARGCHECK_ALWAYS(false,"Eigen","getFewEigenVectors",
       "getFewEigenVectors called with vector of type <std::complex<double> >   \n");
   }
   bool isFactored;

   private:
   void checkIfFactored(const char* function_name)  const {
       if( !isFactored ) {
           SimTK_APIARGCHECK_ALWAYS(false,"Eigen",function_name,
           "matrix has not been factored\n");
       }

       return;
   }


}; // class EigenRepBase

class EigenDefault : public EigenRepBase {
   public:
   EigenDefault();
   EigenRepBase* clone() const;
};


template <typename T>
class EigenRep : public EigenRepBase {
   public:
   template <class ELT> EigenRep( const Matrix_<ELT>&  );

    enum EigenRange {
        AllValues   = 0,
        IndexRange  = 1,
        ValueRange  = 2 
    };

    ~EigenRep();
    EigenRepBase* clone() const;

    typedef typename CNT<T>::TReal RType;
//    template <class VAL, class VEC> void getAllEigenValuesAndVectors( Vector_<VAL>& values, Matrix_<VEC>& vectors);
    void getAllEigenValuesAndVectors( Vector_<std::complex<RType> >& values, Matrix_<std::complex<RType> >& vectors);
    void getAllEigenValuesAndVectors( Vector_<RType>& values, Matrix_<RType>& vectors);
    void getAllEigenValuesAndVectors( Vector_<RType>& values, Matrix_<std::complex<RType> >& vectors);
    void getAllEigenValues( Vector_<RType>& values);
    void getAllEigenValues( Vector_<std::complex<RType> >& values);

    void getFewEigenValuesAndVectors( Vector_<RType>& values, Matrix_<RType>& vectors, int ilow, int ihi);
    void getFewEigenValuesAndVectors( Vector_<std::complex<RType> >& values, Matrix_<std::complex<RType> >& vectors, int ilow, int ihi);
    void getFewEigenValuesAndVectors( Vector_<RType>& values, Matrix_<std::complex<RType> >& vectors, int ilow, int ihi);
    void getFewEigenVectors( Matrix_<RType>& vectors, int ilow, int ihi );
    void getFewEigenVectors( Matrix_<std::complex<RType> >& vectors, int ilow, int ihi );
    void getFewEigenValues( Vector_<RType>& values, int ilow, int ihi );
    void getFewEigenValues( Vector_<std::complex<RType> >& values, int ilow, int ihi );

    void getFewEigenValuesAndVectors( Vector_<RType>& values, Matrix_<RType>& vectors, RType rlow, RType ihi);
    void getFewEigenValuesAndVectors( Vector_<std::complex<RType> >& values, Matrix_<std::complex<RType> >& vectors, RType rlow, RType ihi);
    void getFewEigenValuesAndVectors( Vector_<RType>& values, Matrix_<std::complex<RType> >& vectors, RType rlow, RType ihi);
    void getFewEigenVectors( Matrix_<RType>& vectors, RType rlow, RType ihi );
    void getFewEigenVectors( Matrix_<std::complex<RType> >& vectors, RType rlow, RType ihi );
    void getFewEigenValues( Vector_<RType>& values, RType rlow, RType ihi );
    void getFewEigenValues( Vector_<std::complex<RType> >& values, RType rlow, RType ihi );

    private:
  
    void computeValues(bool);
   
    template <typename P>  void copyVectors( Matrix_<P>& vectors );
    void copyValues( Vector_<float>& values );
    void copyValues( Vector_<double>& values );
    void copyValues( Vector_<std::complex<float> >& values );
    void copyValues( Vector_<std::complex<double> >& values );

    int n;
    int lowIndex, hiIndex;   // min and max indexes for computing a few eigenvalues 
    RType lowValue, hiValue; // min and max values for computing a few eigenvalues
    RType abstol;            // convergence tolerance used by interative eigen routines
    bool needValues;         // true if eigen values  need to be computed
    bool needVectors;        // true if eigen vectors need to be computed
    bool vectorsInMatrix;    // true if eigen vectors stored in inputMatrix
    int valuesFound;         // number of eigen values found when computing few values
    EigenRange range;
    MatrixStructure structure;
    TypedWorkSpace<T> inputMatrix;
    TypedWorkSpace< RType> realEigenValues;
    TypedWorkSpace< std::complex<RType> > complexEigenValues;
    TypedWorkSpace< std::complex<RType> > complexEigenVectors;
    TypedWorkSpace<int> ifail;   // indexes of eigenvalues that did not converge
    TypedWorkSpace<T> symmetricEigenVectors;

}; // end class EigenRep
} // namespace SimTK
#endif   //  _EIGEN_REP_H_
