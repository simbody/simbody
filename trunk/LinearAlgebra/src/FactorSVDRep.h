#ifndef _FACTORSVD_REP_H_ 
#define _FACTORSVD_REP_H_
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

class FactorSVDRepBase {
    public:

    virtual ~FactorSVDRepBase(){};

    virtual FactorSVDRepBase* clone() const { return 0; };


    virtual void getSingularValues( Vector_<float>& values ){
        checkIfFactored( "getSingularValues" );
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","getSingularValues",
        "getSingularValues( float) called with types that are inconsistant with the original matrix  \n");
    }
    virtual void getSingularValues( Vector_<double>& values ){
        checkIfFactored( "getSingularValues" );
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","getSingularValues",
       "getSingularValues( double) called with types that are inconsistant with the original matrix  \n");
    }
    virtual void getSingularValuesAndVectors( Vector_<float>& values, Matrix_<float>& leftVectors, Matrix_<float>& rightVectors ){
        checkIfFactored( "getSingularValuesAndVectors" );
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","getSingularValuesAndVectors",
        "getSingularValuesAndVectors( float, float, float) called with types that are inconsistant with the original matrix  \n");
    }
    virtual void getSingularValuesAndVectors( Vector_<double>& values, Matrix_<double>& leftVectors, Matrix_<double>& rightVectors ){
        checkIfFactored( "getSingularValuesAndVectors" );
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","getSingularValuesAndVectors",
        "getSingularValuesAndVectors( double, double, double) called with types that are inconsistant with the original matrix  \n");
    }

    virtual void getSingularValuesAndVectors( Vector_<float>& values, Matrix_<std::complex<float> >& leftVectors, Matrix_<std::complex<float> >& rightVectors ){
        checkIfFactored( "getSingularValuesAndVectors" );
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","getSingularValuesAndVectors",
        "getSingularValuesAndVectors( float, std::complex<float> , std::complex<float> ) called with types that are inconsistant with the original matrix  \n");
    }

    virtual void getSingularValuesAndVectors( Vector_<double>& values, Matrix_<std::complex<double> >& leftVectors, Matrix_<std::complex<double> >& rightVectors ){
        checkIfFactored( "getSingularValuesAndVectors" );
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","getSingularValuesAndVectors",
        "getSingularValuesAndVectors( double, std::complex<double> , std::complex<double> ) called with types that are inconsistant with the original matrix  \n");
    }
    virtual void solve( const Vector_<float>& b, Vector_<float>& x ) {
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","solve",
        "solve called with rhs of type <float>  which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<double>& b, Vector_<double>& x ){
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","solve",
        "solve called with rhs of type <double>  which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<std::complex<float> >& b, Vector_<std::complex<float> >& x ) {
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","solve",
        "solve called with rhs of type complex<float> which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<std::complex<double> >& b, Vector_<std::complex<double> >& x ) {
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","solve",
        "solve called with rhs of type complex<double>  which does not match type of original linear system \n");   
   }
    virtual void solve( const Matrix_<float>& b, Matrix_<float>& x ) {
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","solve",
        "solve called with rhs of type <float>  which does not match type of original linear system \n");
   }
   virtual void solve( const Matrix_<double>& b, Matrix_<double>& x ) {
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","solve",
        "solve called with rhs of type <double>  which does not match type of original linear system \n");
   }
   virtual void solve( const Matrix_<std::complex<float> >& b, Matrix_<std::complex<float> >& x ) {
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","solve",
        "solve called with rhs of type complex<float> which does not match type of original linear system \n");
   }
   virtual void solve  ( const Matrix_<std::complex<double> >& b, Matrix_<std::complex<double> >& x ) {
       checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","solve",
        "solve called with rhs of type complex<double>  which does not match type of original linear system \n");   
   }
    virtual void inverse(  Matrix_<double>& inverse ){
        checkIfFactored( "inverse" );
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","inverse",
        "inverse(  <double> ) called with type that is inconsistant with the original matrix  \n");
    }
    virtual void inverse(  Matrix_<float>& inverse ){
        checkIfFactored( "inverse" );
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","inverse",
        "inverse(  <float> ) called with type that is inconsistant with the original matrix  \n");
    }
    virtual void inverse(  Matrix_<std::complex<float> >& inverse ){
        checkIfFactored( "inverse" );
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","inverse",
        "inverse(  std::complex<float> ) called with type that is inconsistant with the original matrix  \n");
    }
    virtual void inverse(  Matrix_<std::complex<double> >& inverse ){
        checkIfFactored( "inverse" );
        SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD","inverse",
        "inverse(  std::complex<double> ) called with type that is inconsistant with the original matrix  \n");
    }
    virtual int getRank() const {
       checkIfFactored( "getRank" );
	   return(0);
    }

    void checkIfFactored(const char* function_name)  const {
        if( !isFactored ) {
            SimTK_APIARGCHECK_ALWAYS(false,"FactorSVD",function_name,
            "matrix has not been specified\n");
        }
        return;
    }
    bool isFactored;




}; // class FactorSVDRepBase

class FactorSVDDefault : public FactorSVDRepBase {
   public:
   FactorSVDDefault();
    FactorSVDRepBase* clone() const;
};


template <typename T>
class FactorSVDRep : public FactorSVDRepBase {
   public:
   template <class ELT> FactorSVDRep( const Matrix_<ELT>&, typename CNT<T>::TReal  );

    ~FactorSVDRep();
    FactorSVDRepBase* clone() const;

    typedef typename CNT<T>::TReal RType;

    void getSingularValuesAndVectors( Vector_<RType>& values,   Matrix_<T>& leftVectors,  Matrix_<T>& rightVectors );
    void getSingularValues( Vector_<RType>& values );
    int getRank();
    void solve( const Vector_<T>& b, Vector_<T>& x );
    void solve( const Matrix_<T>& b, Matrix_<T>& x );


    private:

    void computeSVD( bool, RType*, T*, T* );
    void doSolve( Matrix_<T>& b, Matrix_<T>& x );
    void inverse( Matrix_<T>& b );

    int nCol;       // number of columns in original matrix
    int nRow;       // number of rows in original matrix
    int mn;      // min(m,n)
    int maxmn;   // max(m,n)
    int rank;
    RType rcond;   // reciprocol condition number
    RType abstol;
    MatrixStructure structure;
    TypedWorkSpace<T> inputMatrix;
    TypedWorkSpace<RType> singularValues;

}; // end class FactorSVDRep
} // namespace SimTK
#endif   //  _FACTORSVD_REP_H_
