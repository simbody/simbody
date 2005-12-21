
#ifndef __matrix__hh__
#define __matrix__hh__

template<class T>
class GenericMatrix {
public:
  typedef T ElementType;
};

template<class T>
class SymmetricMatrix : public GenericMatrix<T> {
};

template<class T>
class FullMatrix : public GenericMatrix<T> {
};



#endif /* __matrix__hh__ */
