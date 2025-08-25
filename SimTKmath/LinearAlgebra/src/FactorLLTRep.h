#ifndef SimTK_SIMMATH_FACTORLLT_REP_H_
#define SimTK_SIMMATH_FACTORLLT_REP_H_

#include "SimTKmath.h"
#include "WorkSpace.h"

namespace SimTK {

class FactorLLTRepBase {
   public:
    virtual ~FactorLLTRepBase() {}

    virtual FactorLLTRepBase* clone() const { return 0; }

    virtual void solve(const Vector_<float>& b, Vector_<float>& x) const {
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "solve",
            "solve called with rhs of type <float>  which does not match type "
            "of original linear system \n");
    }
    virtual void solve(const Vector_<double>& b, Vector_<double>& x) const {
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "solve",
            " solve called with rhs of type <double>  which does not match "
            "type of original linear system \n");
    }
    virtual void solve(const Vector_<std::complex<float> >& b,
                       Vector_<std::complex<float> >& x) const {
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "solve",
            " solve called with rhs of type complex<float> which does not "
            "match type of original linear system \n");
    }
    virtual void solve(const Vector_<std::complex<double> >& b,
                       Vector_<std::complex<double> >& x) const {
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "solve",
            " solve called with rhs of type complex<double>  which does not "
            "match type of original linear system \n");
    }
    virtual void solve(const Matrix_<float>& b, Matrix_<float>& x) const {
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "solve",
            " solve called with rhs of type <float>  which does not match type "
            "of original linear system \n");
    }
    virtual void solve(const Matrix_<double>& b, Matrix_<double>& x) const {
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "solve",
            " solve called with rhs of type <double>  which does not match "
            "type of original linear system \n");
    }
    virtual void solve(const Matrix_<std::complex<float> >& b,
                       Matrix_<std::complex<float> >& x) const {
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "solve",
            " solve called with rhs of type complex<float> which does not "
            "match type of original linear system \n");
    }
    virtual void solve(const Matrix_<std::complex<double> >& b,
                       Matrix_<std::complex<double> >& x) const {
        checkIfFactored("solve");
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "solve",
            " solve called with rhs of type complex<double>  which does not "
            "match type of original linear system \n");
    }
    virtual void getL(Matrix_<float>& l) const {
        checkIfFactored("getL");
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "getL",
            " getL called with L of type <float> which does not match type of "
            "original linear system \n");
    }
    virtual void getL(Matrix_<double>& l) const {
        checkIfFactored("getL");
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "getL",
            " getL called with L of type <double>  which does not match type "
            "of original linear system \n");
    }
    virtual void getL(Matrix_<std::complex<float> >& l) const {
        checkIfFactored("getL");
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "getL",
            " getL called with L of type complex<float> which does not match "
            "type of original linear system \n");
    }
    virtual void getL(Matrix_<std::complex<double> >& l) const {
        checkIfFactored("getL");
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "getL",
            " getL called with L of type complex<double>  which does not match "
            "type of original linear system \n");
    }

    virtual void inverse(Matrix_<double>& inverse) const {
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "inverse",
            "inverse(  <double> ) called with type that is inconsistent with "
            "the original matrix  \n");
    }
    virtual void inverse(Matrix_<float>& inverse) const {
        SimTK_APIARGCHECK_ALWAYS(false, "FactorLLT", "inverse",
                                 "inverse(  <float> ) called with type that is "
                                 "inconsistent with the original matrix  \n");
    }
    virtual void inverse(Matrix_<std::complex<float> >& inverse) const {
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "inverse",
            "inverse(  std::complex<float> ) called with type that is "
            "inconsistent with the original matrix  \n");
    }
    virtual void inverse(Matrix_<std::complex<double> >& inverse) const {
        SimTK_APIARGCHECK_ALWAYS(
            false, "FactorLLT", "inverse",
            "inverse(  std::complex<double> ) called with type that is "
            "inconsistent with the original matrix  \n");
    }

   protected:
    void checkIfFactored(const char* fn) const {
        if (!isFactored) {
            SimTK_APIARGCHECK_ALWAYS(false, "FactorLLT", fn,
                                     "Matrix has not been factorized.");
        }
    }

    bool isFactored = false;
};  // class FactorLLTRepBase

class FactorLLTDefault : public FactorLLTRepBase {
   public:
    FactorLLTDefault();
    FactorLLTRepBase* clone() const override;
};

template <typename T>
class FactorLLTRep : public FactorLLTRepBase {
   public:
    template <class ELT>
    FactorLLTRep(const Matrix_<ELT>&);
    FactorLLTRep();

    ~FactorLLTRep();
    FactorLLTRepBase* clone() const override;

    template <class ELT>
    void factor(const Matrix_<ELT>&);
    void solve(const Vector_<T>& b, Vector_<T>& x) const override;
    void solve(const Matrix_<T>& b, Matrix_<T>& x) const override;
    void inverse(Matrix_<T>& m) const override;

    void getL(Matrix_<T>& l) const override;

   private:
    template <class ELT>
    int getType(ELT*);
    int nRow;
    int nCol;
    int mn;  // min(m,n)
    int singularIndex;

    TypedWorkSpace<T> l;
};  // class FactorLLTRep
}  // namespace SimTK

#endif  // SimTK_SIMMATH_FACTORLLT_REP_H_
