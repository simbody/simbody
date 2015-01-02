#ifndef SimTK_SimTKCOMMON_VECTOR_MATH_H_
#define SimTK_SimTKCOMMON_VECTOR_MATH_H_

/* -------------------------------------------------------------------------- *
 *                       Simbody(tm): SimTKcommon                             *
 * -------------------------------------------------------------------------- *
 * This is part of the SimTK biosimulation toolkit originating from           *
 * Simbios, the NIH National Center for Physics-Based Simulation of           *
 * Biological Structures at Stanford, funded under the NIH Roadmap for        *
 * Medical Research, grant U54 GM072970. See https://simtk.org/home/simbody.  *
 *                                                                            *
 * Portions copyright (c) 2008-12 Stanford University and the Authors.        *
 * Authors: Peter Eastman                                                     *
 * Contributors: Michael Sherman                                              *
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

#include "SimTKcommon/basics.h"
#include "SimTKcommon/Simmatrix.h"

#include <cmath>     // for std:sin, sqrt, etc.
#include <algorithm> // for std:sort, nth_element, etc.

/** @file
 * This file defines a large number of standard math functions that can be 
 * applied to vectors and matrices (both the large matrix and small matrix 
 * classes).
 */

namespace SimTK {

// We can use a single definition for a number of functions that simply call a 
// function on each element, returning a value of the same type.
// Note that some of these intentionally copy their argument for use as a temp.

#define SimTK_ELEMENTWISE_FUNCTION(func)               \
template <class ELEM>                                  \
VectorBase<ELEM> func(const VectorBase<ELEM>& v) {     \
    const int size = v.size();                         \
    Vector_<ELEM> temp(size);                          \
    for (int i = 0; i < size; ++i)                     \
        temp[i] = std::func(v[i]);                     \
    return temp;                                       \
}                                                      \
template <class ELEM>                                  \
RowVectorBase<ELEM> func(const RowVectorBase<ELEM>& v){\
    const int size = v.size();                         \
    RowVector_<ELEM> temp(size);                       \
    for (int i = 0; i < size; ++i)                     \
        temp[i] = std::func(v[i]);                     \
    return temp;                                       \
}                                                      \
template <class ELEM>                                  \
MatrixBase<ELEM> func(const MatrixBase<ELEM>& v) {     \
    const int rows = v.nrow(), cols = v.ncol();        \
    Matrix_<ELEM> temp(rows, cols);                    \
    for (int i = 0; i < rows; ++i)                     \
        for (int j = 0; j < cols; ++j)                 \
            temp(i, j) = std::func(v(i, j));           \
    return temp;                                       \
}                                                      \
template <int N, class ELEM>                           \
Vec<N, ELEM> func(Vec<N, ELEM> v) {                    \
    for (int i = 0; i < N; ++i)                        \
        v[i] = std::func(v[i]);                        \
    return v;                                          \
}                                                      \
template <int N, class ELEM>                           \
Row<N, ELEM> func(Row<N, ELEM> v) {                    \
    for (int i = 0; i < N; ++i)                        \
        v[i] = std::func(v[i]);                        \
    return v;                                          \
}                                                      \
template <int M, int N, class ELEM>                    \
Mat<M, N, ELEM> func(Mat<M, N, ELEM> v) {              \
    for (int i = 0; i < M; ++i)                        \
        for (int j = 0; j < N; ++j)                    \
            v(i, j) = std::func(v(i, j));              \
    return v;                                          \
}                                                      \
template <int N, class ELEM>                           \
SymMat<N, ELEM> func(SymMat<N, ELEM> v) {              \
    for (int i = 0; i < N; ++i)                        \
        for (int j = 0; j <= i; ++j)                   \
            v(i, j) = std::func(v(i, j));              \
    return v;                                          \
}                                                      \

SimTK_ELEMENTWISE_FUNCTION(exp)
SimTK_ELEMENTWISE_FUNCTION(log)
SimTK_ELEMENTWISE_FUNCTION(sqrt)
SimTK_ELEMENTWISE_FUNCTION(sin)
SimTK_ELEMENTWISE_FUNCTION(cos)
SimTK_ELEMENTWISE_FUNCTION(tan)
SimTK_ELEMENTWISE_FUNCTION(asin)
SimTK_ELEMENTWISE_FUNCTION(acos)
SimTK_ELEMENTWISE_FUNCTION(atan)
SimTK_ELEMENTWISE_FUNCTION(sinh)
SimTK_ELEMENTWISE_FUNCTION(cosh)
SimTK_ELEMENTWISE_FUNCTION(tanh)

#undef SimTK_ELEMENTWISE_FUNCTION

// The abs() function.

template <class ELEM>
VectorBase<typename CNT<ELEM>::TAbs> abs(const VectorBase<ELEM>& v) {
    return v.abs();
}
template <class ELEM>
RowVectorBase<typename CNT<ELEM>::TAbs> abs(const RowVectorBase<ELEM>& v) {
    return v.abs();
}
template <class ELEM>
MatrixBase<typename CNT<ELEM>::TAbs> abs(const MatrixBase<ELEM>& v) {
    return v.abs();
}
template <int N, class ELEM>
Vec<N, typename CNT<ELEM>::TAbs> abs(const Vec<N, ELEM>& v) {
    return v.abs();
}
template <int N, class ELEM>
Row<N, typename CNT<ELEM>::TAbs> abs(const Row<N, ELEM>& v) {
    return v.abs();
}
template <int M, int N, class ELEM>
Mat<M, N, typename CNT<ELEM>::TAbs> abs(const Mat<M, N, ELEM>& v) {
    return v.abs();
}
template <int N, class ELEM>
SymMat<N, typename CNT<ELEM>::TAbs> abs(const SymMat<N, ELEM>& v) {
    return v.abs();
}

// The sum() function.

template <class ELEM>
ELEM sum(const VectorBase<ELEM>& v) {
    return v.sum();
}
template <class ELEM>
ELEM sum(const RowVectorBase<ELEM>& v) {
    return v.sum();
}
template <class ELEM>
RowVectorBase<ELEM> sum(const MatrixBase<ELEM>& v) {
    return v.sum();
}
template <int N, class ELEM>
ELEM sum(const Vec<N, ELEM>& v) {
    return v.sum();
}
template <int N, class ELEM>
ELEM sum(const Row<N, ELEM>& v) {
    return v.sum();
}
template <int M, int N, class ELEM>
Row<N, ELEM> sum(const Mat<M, N, ELEM>& v) {
    return v.sum();
}
template <int N, class ELEM>
Row<N, ELEM> sum(const SymMat<N, ELEM>& v) {
    return v.sum();
}

// The min() function.

template <class ELEM>
ELEM min(const VectorBase<ELEM>& v) {
    const int size = v.size();
    ELEM min = NTraits<ELEM>::getMostPositive();
    for (int i = 0; i < size; ++i) {
        ELEM val = v[i];
        if (val < min)
            min = val;
    }
    return min;
}
template <class ELEM>
ELEM min(const RowVectorBase<ELEM>& v) {
    const int size = v.size();
    ELEM min = NTraits<ELEM>::getMostPositive();
    for (int i = 0; i < size; ++i) {
        ELEM val = v[i];
        if (val < min)
            min = val;
    }
    return min;
}
template <class ELEM>
RowVectorBase<ELEM> min(const MatrixBase<ELEM>& v) {
    int cols = v.ncol();
    RowVectorBase<ELEM> temp(cols);
    for (int i = 0; i < cols; ++i)
        temp[i] = min(v(i));
    return temp;
}
template <int N, class ELEM>
ELEM min(const Vec<N, ELEM>& v) {
    ELEM min = NTraits<ELEM>::getMostPositive();
    for (int i = 0; i < N; ++i) {
        ELEM val = v[i];
        if (val < min)
            min = val;
    }
    return min;
}
template <int N, class ELEM>
ELEM min(const Row<N, ELEM>& v) {
    ELEM min = NTraits<ELEM>::getMostPositive();
    for (int i = 0; i < N; ++i) {
        ELEM val = v[i];
        if (val < min)
            min = val;
    }
    return min;
}
template <int M, int N, class ELEM>
Row<N, ELEM> min(const Mat<M, N, ELEM>& v) {
    Row<N, ELEM> temp;
    for (int i = 0; i < N; ++i)
        temp[i] = min(v(i));
    return temp;
}
template <int N, class ELEM>
Row<N, ELEM> min(const SymMat<N, ELEM>& v) {
    Row<N, ELEM> temp(~v.getDiag());
    for (int i = 1; i < N; ++i)
        for (int j = 0; j < i; ++j) {
            ELEM value = v.getEltLower(i, j);
            if (value < temp[i])
                temp[i] = value;
            if (value < temp[j])
                temp[j] = value;
        }
    return temp;
}

// The max() function.

template <class ELEM>
ELEM max(const VectorBase<ELEM>& v) {
    const int size = v.size();
    ELEM max = NTraits<ELEM>::getMostNegative();
    for (int i = 0; i < size; ++i) {
        ELEM val = v[i];
        if (val > max)
            max = val;
    }
    return max;
}
template <class ELEM>
ELEM max(const RowVectorBase<ELEM>& v) {
    const int size = v.size();
    ELEM max = NTraits<ELEM>::getMostNegative();
    for (int i = 0; i < size; ++i) {
        ELEM val = v[i];
        if (val > max)
            max = val;
    }
    return max;
}
template <class ELEM>
RowVectorBase<ELEM> max(const MatrixBase<ELEM>& v) {
    int cols = v.ncol();
    RowVectorBase<ELEM> temp(cols);
    for (int i = 0; i < cols; ++i)
        temp[i] = max(v(i));
    return temp;
}
template <int N, class ELEM>
ELEM max(const Vec<N, ELEM>& v) {
    ELEM max = NTraits<ELEM>::getMostNegative();
    for (int i = 0; i < N; ++i) {
        ELEM val = v[i];
        if (val > max)
            max = val;
    }
    return max;
}
template <int N, class ELEM>
ELEM max(const Row<N, ELEM>& v) {
    ELEM max = NTraits<ELEM>::getMostNegative();
    for (int i = 0; i < N; ++i) {
        ELEM val = v[i];
        if (val > max)
            max = val;
    }
    return max;
}
template <int M, int N, class ELEM>
Row<N, ELEM> max(const Mat<M, N, ELEM>& v) {
    Row<N, ELEM> temp;
    for (int i = 0; i < N; ++i)
        temp[i] = max(v(i));
    return temp;
}
template <int N, class ELEM>
Row<N, ELEM> max(const SymMat<N, ELEM>& v) {
    Row<N, ELEM> temp(~v.getDiag());
    for (int i = 1; i < N; ++i)
        for (int j = 0; j < i; ++j) {
            ELEM value = v.getEltLower(i, j);
            if (value > temp[i])
                temp[i] = value;
            if (value > temp[j])
                temp[j] = value;
        }
    return temp;
}

// The mean() function.

template <class ELEM>
ELEM mean(const VectorBase<ELEM>& v) {
    return sum(v)/v.size();
}
template <class ELEM>
ELEM mean(const RowVectorBase<ELEM>& v) {
    return sum(v)/v.size();
}
template <class ELEM>
RowVectorBase<ELEM> mean(const MatrixBase<ELEM>& v) {
    return sum(v)/v.nrow();
}
template <int N, class ELEM>
ELEM mean(const Vec<N, ELEM>& v) {
    return sum(v)/N;
}
template <int N, class ELEM>
ELEM mean(const Row<N, ELEM>& v) {
    return sum(v)/N;
}
template <int M, int N, class ELEM>
Row<N, ELEM> mean(const Mat<M, N, ELEM>& v) {
    return sum(v)/M;
}
template <int N, class ELEM>
Row<N, ELEM> mean(const SymMat<N, ELEM>& v) {
    return sum(v)/N;
}

// The sort() function.

template <class ELEM>
VectorBase<ELEM> sort(const VectorBase<ELEM>& v) {
    VectorBase<ELEM> temp(v);
    std::sort(temp.begin(), temp.end());
    return temp;
}
template <class ELEM>
RowVectorBase<ELEM> sort(const RowVectorBase<ELEM>& v) {
    RowVectorBase<ELEM> temp(v);
    std::sort(temp.begin(), temp.end());
    return temp;
}
template <class ELEM>
MatrixBase<ELEM> sort(const MatrixBase<ELEM>& v) {
    const int cols = v.ncol();
    MatrixBase<ELEM> temp(v);
    for (int i = 0; i < cols; ++i)
        temp.updCol(i) = sort(temp.col(i));
    return temp;
}
template <int N, class ELEM>
Vec<N, ELEM> sort(Vec<N, ELEM> v) { // intentional copy of argument
    ELEM* pointer = reinterpret_cast<ELEM*>(&v);
    std::sort(pointer, pointer+N);
    return v;
}
template <int N, class ELEM>
Row<N, ELEM> sort(Row<N, ELEM> v) { // intentional copy of argument
    ELEM* pointer = reinterpret_cast<ELEM*>(&v);
    std::sort(pointer, pointer+N);
    return v;
}
template <int M, int N, class ELEM>
Mat<M, N, ELEM> sort(Mat<M, N, ELEM> v) { // intentional copy of argument
    for (int i = 0; i < N; ++i)
        v.col(i) = sort(v.col(i));
    return v;
}
template <int N, class ELEM>
Mat<N, N, ELEM> sort(const SymMat<N, ELEM>& v) {
    return sort(Mat<N, N, ELEM>(v));
}

// The median() function.

template <class ELEM, class RandomAccessIterator>
ELEM median(RandomAccessIterator start, RandomAccessIterator end) {
    const ptrdiff_t size = (ptrdiff_t)(end-start);
    RandomAccessIterator mid = start+(size-1)/2;
    std::nth_element(start, mid, end);
    if (size%2 == 0 && mid+1 < end) {
        // An even number of element.  The median is the mean of the two middle elements.
        // nth_element has given us the first of them and partially sorted the list.
        // We need to scan through the rest of the list and find the next element in
        // sorted order.
        
        RandomAccessIterator min = mid+1;
        for (RandomAccessIterator iter = mid+1; iter < end; iter++) {
            if (*iter < *min)
                min = iter;
        }
        return (*mid+*min)/2;
    }
    return *mid;
}
template <class ELEM>
ELEM median(const VectorBase<ELEM>& v) {
    VectorBase<ELEM> temp(v);
    return median<ELEM>(temp.begin(), temp.end());
}
template <class ELEM>
ELEM median(const RowVectorBase<ELEM>& v) {
    RowVectorBase<ELEM> temp(v);
    return median<ELEM>(temp.begin(), temp.end());
}
template <class ELEM>
RowVectorBase<ELEM> median(const MatrixBase<ELEM>& v) {
    int cols = v.ncol(), rows = v.nrow();
    RowVectorBase<ELEM> temp(cols);
    VectorBase<ELEM> column(rows);
    for (int i = 0; i < cols; ++i) {
        column = v.col(i);
        temp[i] = median<ELEM>(column);
    }
    return temp;
}
template <int N, class ELEM>
ELEM median(Vec<N, ELEM> v) { // intentional copy of argument
    ELEM* pointer = reinterpret_cast<ELEM*>(&v);
    return  median<ELEM>(pointer, pointer+N);
}
template <int N, class ELEM>
ELEM median(Row<N, ELEM> v) { // intentional copy of argument
    ELEM* pointer = reinterpret_cast<ELEM*>(&v);
    return  median<ELEM>(pointer, pointer+N);
}
template <int M, int N, class ELEM>
Row<N, ELEM> median(const Mat<M, N, ELEM>& v) {
    Row<N, ELEM> temp;
    for (int i = 0; i < N; ++i)
        temp[i] = median(v(i));
    return temp;
}
template <int N, class ELEM>
Row<N, ELEM> median(const SymMat<N, ELEM>& v) {
    return median(Mat<N, N, ELEM>(v));
}

} // namespace SimTK

#endif // SimTK_SimTKCOMMON_VECTOR_MATH_H_
