#ifndef SimTK_SIMMATH_FACTOR_QTZ_REP_H_ 
#define SimTK_SIMMATH_FACTOR_QTZ_REP_H_

/* -------------------------------------------------------------------------- *
 *                        Simbody(tm): SimTKmath                              *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2007-12 Stanford University and the Authors.        *
 * Authors: Jack Middleton                                                    *
 * Contributors:                                                              *
 *                                                                            *
 * Licensed under the Apache License, Version 2.0 (the "License"); you may    *
 * not use this file except in compliance with the License. You may obtain a  *
 * copy of the License at http://www.apache.org/licenses/LICENSE-2.0.         *
 *                                                                            *
 * Unless required by applicable law or agreed to in writing, software        *
 * distributed under the License is distributed on an "AS IS" BASIS,          *
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.   *
 * See the License for the specific language governing permissions and        *
 * limitations under the License.                                             *
 * -------------------------------------------------------------------------- */

#include "SimTKmath.h"
#include "WorkSpace.h"

namespace SimTK {

class FactorQTZRepBase {
public:
    FactorQTZRepBase() : isFactored(false), rank(0), actualRCond(0) {}

    virtual ~FactorQTZRepBase(){};

    virtual FactorQTZRepBase* clone() const { return 0; };
    virtual void solve( const Vector_<float>& b, Vector_<float>& x ) const {
        checkIfFactored();
        SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
        "solve called with rhs of type <float>  which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<double>& b, Vector_<double>& x ) const {
        checkIfFactored();
        SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
        "solve called with rhs of type <double>  which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<std::complex<float> >& b, Vector_<std::complex<float> >& x ) const {
        checkIfFactored();
        SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
        "solve called with rhs of type complex<float> which does not match type of original linear system \n");
   }
   virtual void solve( const Vector_<std::complex<double> >& b, Vector_<std::complex<double> >& x ) const {
        checkIfFactored();
        SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
        "solve called with rhs of type complex<double>  which does not match type of original linear system \n");   
   }
    virtual void solve( const Matrix_<float>& b, Matrix_<float>& x ) const {
        checkIfFactored();
        SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
        "solve called with rhs of type <float>  which does not match type of original linear system \n");
   }
   virtual void solve( const Matrix_<double>& b, Matrix_<double>& x ) const {
        checkIfFactored();
        SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
        "solve called with rhs of type <double>  which does not match type of original linear system \n");
   }
   virtual void solve( const Matrix_<std::complex<float> >& b, Matrix_<std::complex<float> >& x ) const {
        checkIfFactored();
        SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
        "solve called with rhs of type complex<float> which does not match type of original linear system \n");
   }
   virtual void solve  ( const Matrix_<std::complex<double> >& b, Matrix_<std::complex<double> >& x ) const {
       checkIfFactored();
        SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
        "solve called with rhs of type complex<double>  which does not match type of original linear system \n");   
   }
    virtual void inverse(  Matrix_<double>& inverse ) const {
        SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","inverse",
        "inverse(  <double> ) called with type that is inconsistant with the original matrix  \n");
    }
    virtual void inverse(  Matrix_<float>& inverse ) const {
        SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","inverse",
        "inverse(  <float> ) called with type that is inconsistant with the original matrix  \n");
    }
    virtual void inverse(  Matrix_<std::complex<float> >& inverse ) const {
        SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","inverse",
        "inverse(  std::complex<float> ) called with type that is inconsistant with the original matrix  \n");
    }
    virtual void inverse(  Matrix_<std::complex<double> >& inverse ) const {
        SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","inverse",
        "inverse(  std::complex<double> ) called with type that is inconsistant with the original matrix  \n");
    }


   bool isFactored;
   int rank;     // esitmated rank computed during factorization
   double actualRCond; // 1 / condition number we actually got (est.)


   void checkIfFactored()  const {
       if( !isFactored ) {
           SimTK_APIARGCHECK_ALWAYS(false,"FactorQTZ","solve",
           "solve called before the matrix was factored \n");
       }
   }

}; // class FactorQTZRepBase

class FactorQTZDefault : public FactorQTZRepBase {
   public:
	   FactorQTZDefault();
	   FactorQTZRepBase* clone() const;
};

template <typename T>
class FactorQTZRep : public FactorQTZRepBase {
public:
   template <class ELT> FactorQTZRep( const Matrix_<ELT>&, typename CNT<T>::TReal );
   FactorQTZRep();

   ~FactorQTZRep();

   template < class ELT > void factor(const Matrix_<ELT>& ); 
   void inverse( Matrix_<T>& ) const; 
   void solve( const Vector_<T>& b, Vector_<T>& x ) const;
   void solve( const Matrix_<T>& b, Matrix_<T>& x ) const;

   FactorQTZRepBase* clone() const;
 
private:
   void doSolve( Matrix_<T>& b, Matrix_<T>& x ) const;

   int                      mn;           // min of number of rows or columns
   int                      maxmn;        // max of number of rows or columns
   int                      nRow;         // number of rows in original matrix
   int                      nCol;         // number of columns in original matrix
   bool                     scaleLinSys; // true if matrix was scaled during factorization
   typename CNT<T>::TReal   linSysScaleF; // scale factor applied to matrix 
   typename CNT<T>::TReal   anrm;
   typename CNT<T>::TReal   rcond;      // 1 / (largest acceptable condition number)
   TypedWorkSpace<int>      pivots;
   TypedWorkSpace<T>        qtz;     // factored matrix
   TypedWorkSpace<T>        tauGEQP3;
   TypedWorkSpace<T>        tauORMQR;

}; // end class FactorQTZRep

} // namespace SimTK

#endif   // SimTK_SIMMATH_FACTOR_QTZ_REP_H_
