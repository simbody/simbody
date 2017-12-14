
#ifndef __IPLAPACKSOLVERINTERFACE_HPP__
#define __IPLAPACKSOLVERINTERFACE_HPP__

#include "IpSparseSymLinearSolverInterface.hpp"

namespace SimTKIpopt
{

  /** Interface to the linear solver Lapack, derived from
   *  SparseSymLinearSolverInterface.  For details, see description of
   *  SparseSymLinearSolverInterface base class.
   */
  class LapackSolverInterface: public SparseSymLinearSolverInterface
  {
  public:
    /** @name Constructor/Destructor */
    //@{
    /** Constructor */
    LapackSolverInterface();

    /** Destructor */
    virtual ~LapackSolverInterface();
    //@}

    /** overloaded from AlgorithmStrategyObject */
    bool InitializeImpl(const OptionsList& options,
                        const std::string& prefix) override;


    /** @name Methods for requesting solution of the linear system. */
    //@{
    /** Method for initializing internal stuctures. */
    virtual ESymSolverStatus InitializeStructure(Index dim, Index nonzeros, const Index *ia, const Index *ja) override;

    /** Method returning an internal array into which the nonzero
     *  elements are to be stored. */
    virtual Number* GetValuesArrayPtr() override;

    /** Solve operation for multiple right hand sides. */
    virtual ESymSolverStatus MultiSolve(bool new_matrix,
                                        const Index* ia,
                                        const Index* ja,
                                        Index nrhs,
                                        Number* rhs_vals,
                                        bool check_NegEVals,
                                        Index numberOfNegEVals) override;

    /** Number of negative eigenvalues detected during last
     *  factorization.
     */
    virtual Index NumberOfNegEVals() const override;
    //@}

    //* @name Options of Linear solver */
    //@{
    /** Request to increase quality of solution for next solve.
     */
    virtual bool IncreaseQuality() override;

    /** Query whether inertia is computed by linear solver.
     *  Returns true, if linear solver provides inertia.
     */
    virtual bool ProvidesInertia() const override
    {
      return true;
    }
    /** Query of requested matrix type that the linear solver
     *  understands.
     */
    EMatrixFormat MatrixFormat() const override
    {
      return Dense_Format;
    }
    //@}

    /** Methods for IpoptType */
    //@{
    static void RegisterOptions(SmartPtr<RegisteredOptions> roptions);
    //@}

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Copy Constructor */
    LapackSolverInterface(const LapackSolverInterface&);

    /** Overloaded Equals Operator */
    void operator=(const LapackSolverInterface&);
    //@}

    /** @name Information about the matrix */
    //@{
    /** @name Information about the matrix */
    //@{
    /** Number of rows and columns of the matrix */
    Index n;

    /** Number of nonzeros of the matrix in triplet representation. */
    Index nz;

    /** Array for storing the values of the matrix. */
    Number* a;

    /** Array for storing the values of the factored matrix. */
    Number* afact;

    /** Array for storing the row indices of the matrix */
    int* irn_;
    /** Array for storing the column indices of the matrix */
    int* jcn_;
    //@}
    /** @name Information about most recent factorization/solve */
    //@{
    /** Number of negative eigenvalues */
    Index negevals_;
    /** Array for storing the pivot order after factorization. */
    int *ipiv_;
    bool isFactored;

    //@}

    /** @name Solver specific options */
    //@{
    //@}

    /** @name Initialization flags */
    //@{
    //@}

    /** @name Solver specific information */
    //@{
    /**@name Some counters for debugging */
    //@{
    //@}

    /** @name Internal functions */
    //@{
    /** Call Lapack to do the analysis phase.
     */
    /** Call Lapack to factorize the Matrix.
     */
    ESymSolverStatus Factorization(const Index* ia,
                                   const Index* ja,
                                   bool check_NegEVals,
                                   Index numberOfNegEVals);

    /** Call Lapack to do the Solve.
     */
    ESymSolverStatus Solve(const Index* ia,
                           const Index* ja,
                           Index nrhs,
                           Number *rhs_vals);
    //@}
    //MUMPS data structure
//    DMUMPS_STRUC_C mumps_data;

  };

} // namespace Ipopt
#endif
