
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
/**@file
 *
 * Solves for eigen values and eigen vectors 
 */


#include <iostream> 
#include <malloc.h>
#include <math.h>
#include <complex>
#include "SimTKcommon.h"
#include "LapackInterface.h"
#include "LinearAlgebra.h"
#include "EigenRep.h"
#include "WorkSpace.h"
#include "LATraits.h"
#include "LapackConvert.h"


namespace SimTK {

   ///////////////
   // Eigen //
   ///////////////
Eigen::~Eigen() {
    delete rep;
}

template < class ELT >
Eigen::Eigen( const Matrix_<ELT>& m ) {
    rep = new EigenRep<typename CNT<ELT>::StdNumber>(m); 
}

void Eigen::selectVectors( const std::vector<bool>& selectedVectors ) {
    rep->selectVectors( selectedVectors );
}

void Eigen::selectValues( const std::vector<bool>& selectedValues ) {
    rep->selectValues( selectedValues );
}

template < typename T >
void Eigen::getValues( Vector_<T>& values ) {
    rep->getValues( values );
    return;
}
template < class T >
void Eigen::getVectors(  Matrix_<T>& vectors ) {
    rep->getVectors(  vectors );
    return;
}


   /////////////////
   // EigenRep //
   /////////////////
template <typename T >
    template < typename ELT >
EigenRep<T>::EigenRep( const Matrix_<ELT>& mat)  {
        
//	EigenRep<T>::computeValues( mat );
}

template <typename T >
EigenRep<T>::~EigenRep() {}

template < class T >
void EigenRep<T>::getValues( Vector_<std::complex<typename CNT<T>::TReal> >& values ) {
    values.resize(numValues);
    return;
}
template < class T >
void EigenRep<T>::getValues( Vector_<typename CNT<T>::TReal>& values ) {
    values.resize(numValues);
    return;
}

template < class T >
void EigenRep<T>::getVectors( Matrix_<std::complex<typename CNT<T>::TReal> >& vectors ) {
    vectors.resize(nRow,nCol);
    return;
}

template < class T >
void EigenRep<T>::getVectors( Matrix_<typename CNT<T>::TReal>& vectors ) {
    vectors.resize(nRow,nCol);
    return;
}

template < class T >
void EigenRep<T>::selectVectors( const std::vector<bool>& selectedVectors ) {
}

template < class T >
void EigenRep<T>::selectValues( const std::vector<bool>& selectedValues ) {
}



// instantiate
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<double>& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<float>& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<std::complex<float> >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<std::complex<double> >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<conjugate<float> >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<conjugate<double> >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<negator< double> >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<negator< float> >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<negator< std::complex<float> > >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<negator< std::complex<double> > >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<negator< conjugate<float> > >& m );
template SimTK_SIMMATH_EXPORT Eigen::Eigen( const Matrix_<negator< conjugate<double> > >& m );

template class EigenRep<double>;
template EigenRep<double>::EigenRep( const Matrix_<double>& m);
template EigenRep<double>::EigenRep( const Matrix_<negator<double> >& m);

template class EigenRep<float>;
template EigenRep<float>::EigenRep( const Matrix_<float>& m );
template EigenRep<float>::EigenRep( const Matrix_<negator<float> >& m );

template class EigenRep<std::complex<double> >;
template EigenRep<std::complex<double> >::EigenRep( const Matrix_<std::complex<double> >& m );
template EigenRep<std::complex<double> >::EigenRep( const Matrix_<negator<std::complex<double> > >& m );
template EigenRep<std::complex<double> >::EigenRep( const Matrix_<conjugate<double> >& m );
template EigenRep<std::complex<double> >::EigenRep( const Matrix_<negator<conjugate<double> > >& m );

template class EigenRep<std::complex<float> >;
template EigenRep<std::complex<float> >::EigenRep( const Matrix_<std::complex<float> >& m );
template EigenRep<std::complex<float> >::EigenRep( const Matrix_<negator<std::complex<float> > >& m );
template EigenRep<std::complex<float> >::EigenRep( const Matrix_<conjugate<float> >& m );
template EigenRep<std::complex<float> >::EigenRep( const Matrix_<negator<conjugate<float> > >& m );

template SimTK_SIMMATH_EXPORT void Eigen::getValues<float >(Vector_<float>& );
template SimTK_SIMMATH_EXPORT void Eigen::getValues<double>(Vector_<double>& );
template SimTK_SIMMATH_EXPORT void Eigen::getValues<std::complex<float> >(Vector_<std::complex<float> >& );
template SimTK_SIMMATH_EXPORT void Eigen::getValues<std::complex<double> >(Vector_<std::complex<double> >& );
template SimTK_SIMMATH_EXPORT void Eigen::getVectors<float>(Matrix_<float>&);
template SimTK_SIMMATH_EXPORT void Eigen::getVectors<double>(Matrix_<double>& );
template SimTK_SIMMATH_EXPORT void Eigen::getVectors<std::complex<float> >(Matrix_<std::complex<float> >&);
template SimTK_SIMMATH_EXPORT void Eigen::getVectors<std::complex<double> >(Matrix_<std::complex<double> >& );

} // namespace SimTK
