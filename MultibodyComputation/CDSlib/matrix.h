
#ifndef __matrix__hh__
#define __matrix__hh__

template<class T>
class Matrix {
public:
  typedef T ElementType;
};

template<class T>
class SymmetricMatrix : public Matrix<T> {
};

template<class T>
class FullMatrix : public Matrix<T> {
};



#endif /* __matrix__hh__ */
