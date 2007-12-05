
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
 * Implementation numerical linear algebra operations.
 */


#include <iostream> 
#include <malloc.h>
#include <math.h>
#include "SimTKcommon.h"
#include "simmath/LinearAlgebra.h"
#include "simmath/internal/common.h"
#include "LapackInterface.h"


namespace SimTK {

/**
 *  Preliminary implementation!
 *  Wrapper for LAPACK to compute the eigen values and right eigen vectors of a 
 *  symmetric matrix.
 *  Current implemenation only returns the real values of the right eigen vectors. 
 *  Returns true if successful. Returns false if computation failed.
 */

template <class P>
bool calcEigenValuesRightEigenVectors( Matrix_<P> &m, Vector_< std::complex<P> > &eigenValues, Matrix_< std::complex<P> > &eigenVectors ) {

    int i,j;
    int info;
    int inputMatrixNumberOfRowsOrCols = m.ncol();     // dimension of matrix
    int inputMatrixLeadingDimension = m.ncol();       // maxtrix has contiguous data
    int leadingDimensionOfLeftEigenvectorsArray = m.ncol(); 
    int leadingDimensionOfRightEigenvectorsArray = m.ncol();
    char calcLeftEigenVectors  = 'N';     // tell LAPACK not to compute left eigen vectors
    char calcRightEigenVectors = 'V';     // tell LAPACK to compute right eigen vectors
    P size[1];              // dimension of opitimal workspace computed in 1st call to geev 
    P *inputMatrix = new P[inputMatrixNumberOfRowsOrCols*inputMatrixNumberOfRowsOrCols];  // alloc temporary matrix
    P *realPartOfEigenValues = new P[inputMatrixNumberOfRowsOrCols];     // alloc real portion of results
    P *imagPartOfEigenValues = new P[inputMatrixNumberOfRowsOrCols];     // alloc imaginary portion of results
    P *arrayForRightEigenVectors = new P[inputMatrixNumberOfRowsOrCols*inputMatrixNumberOfRowsOrCols]; // alloc right eigen vectors
    P *arrayForLeftEigenVectors = NULL;           // no left vectors are calculated

    /* check if the matrix is square and has no gaps between columns */
    assert( m.hasContiguousData() );
    assert( m.ncol() == m.nrow());


    /* make a copy of input matrix because LAPACK will overwrite it */
    for (j=0;j<inputMatrixNumberOfRowsOrCols;j++) {
        for (i=0;i<inputMatrixNumberOfRowsOrCols;i++) {
            inputMatrix[j*inputMatrixNumberOfRowsOrCols+i] = m(i,j);
        }
    }
    // compute and allocate optimial workspace
    int dimensionOfWorkSpace = -1;
    LapackInterface::geev<P>( calcLeftEigenVectors,     calcRightEigenVectors,
                              inputMatrixNumberOfRowsOrCols, inputMatrix, inputMatrixLeadingDimension,
                              realPartOfEigenValues,    imagPartOfEigenValues,
                              arrayForLeftEigenVectors, leadingDimensionOfLeftEigenvectorsArray,
                              arrayForRightEigenVectors,leadingDimensionOfRightEigenvectorsArray,
                              size, dimensionOfWorkSpace, 
                              info);

    dimensionOfWorkSpace = (int)size[0];
    P *workSpace = new P[dimensionOfWorkSpace]; // allocate workspace
   
    // compute eigen values and right eigen vectors 
    LapackInterface::geev<P>( calcLeftEigenVectors,      calcRightEigenVectors,
                              inputMatrixNumberOfRowsOrCols, inputMatrix, inputMatrixLeadingDimension,
                              realPartOfEigenValues,     imagPartOfEigenValues,
                              arrayForLeftEigenVectors,  leadingDimensionOfLeftEigenvectorsArray,
                              arrayForRightEigenVectors, leadingDimensionOfRightEigenvectorsArray,
                              workSpace,                 dimensionOfWorkSpace, 
                              info);

    // resize output vector and matrix  if necessary 
    eigenValues.resize( inputMatrixNumberOfRowsOrCols );
    eigenVectors.resize( inputMatrixNumberOfRowsOrCols, inputMatrixNumberOfRowsOrCols );

    // copy computed eigen values and eigen vectors into caller's arguements 
    if( info == 0 ) {
        for (j=0;j<inputMatrixNumberOfRowsOrCols;j++) {
            eigenValues(j) = std::complex<P>(realPartOfEigenValues[j], imagPartOfEigenValues[j]);
//            eigenValues(j).real = realPartOfEigenValues[j];
//            eigenValues(j).imag = imagPartOfEigenValues[j];
            for (i=0;i<inputMatrixNumberOfRowsOrCols;i++) {
//               eigenVectors(i,j).real = arrayForRightEigenVectors[j*inputMatrixNumberOfRowsOrCols+i];
               eigenVectors(i,j) = std::complex<P>(arrayForRightEigenVectors[j*inputMatrixNumberOfRowsOrCols+i], 0.0 );
            }
        }
    }
    // free all temporary memory we alloc'ed 
    delete [] inputMatrix;
    delete [] workSpace;
    delete [] imagPartOfEigenValues;
    delete [] realPartOfEigenValues;
    delete [] arrayForRightEigenVectors;
 
    if( info == 0 )   {
        return true;
    } else  {
        return false;
    }
}

// instantiate for only float and double 
template SimTK_SIMMATH_EXPORT bool calcEigenValuesRightEigenVectors( Matrix_<float>&, Vector_< std::complex<float> >&, Matrix_< std::complex<float> >&);
template SimTK_SIMMATH_EXPORT bool calcEigenValuesRightEigenVectors( Matrix_<double>&, Vector_< std::complex<double> >&, Matrix_< std::complex<double> >&);

} // end namespace SimTK
