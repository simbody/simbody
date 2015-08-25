#ifndef SimTK_SIMMATH_FACTOR_REP_H_ 
#define SimTK_SIMMATH_FACTOR_REP_H_

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

class FactorLURepBase {
    public:

    virtual ~FactorLURepBase(){};
    virtual FactorLURepBase* clone() const { return 0; };

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
        SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","inverse",         "inverse(  <double> ) called with type that is inconsistent with the original matrix  \n");
    }
    virtual void inverse(  Matrix_<float>& inverse ) const{
        SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","inverse",         "inverse(  <float> ) called with type that is inconsistent with the original matrix  \n");
    }
    virtual void inverse(  Matrix_<std::complex<float> >& inverse ) const{
        SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","inverse",         "inverse(  std::complex<float> ) called with type that is inconsistent with the original matrix  \n");
    }
    virtual void inverse(  Matrix_<std::complex<double> >& inverse ) const{
        SimTK_APIARGCHECK_ALWAYS(false,"FactorLU","inverse",
        "inverse(  std::complex<double> ) called with type that is inconsistent with the original matrix  \n");
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
   int nRow;
   int nCol;
   int mn;        // min(m,n)
   int singularIndex;

   TypedWorkSpace<int>  pivots;
   TypedWorkSpace<T>    lu;
}; // end class FactorLURep
} // namespace SimTK
#endif   //  SimTK_SIMMATH_FACTOR_REP_H_
