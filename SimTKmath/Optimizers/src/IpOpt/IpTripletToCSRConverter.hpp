// Copyright (C) 2005, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpTripletToCSRConverter.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2005-03-13

#ifndef __IPTRIPLETTOCSRCONVERTER_HPP__
#define __IPTRIPLETTOCSRCONVERTER_HPP__

#include "IpUtils.hpp"
#include "IpReferenced.hpp"
namespace SimTKIpopt
{

  /** Class for converting symmetric matrices given in triplet format
   *  to matrices in compressed sparse row (CSR) format of the upper
   *  triangual part (or, equivalently, compressed sparse column (CSC)
   *  format for the lower triangular part).  In the description for
   *  this class, we assume that we discuss the CSR format.
   */
  class TripletToCSRConverter: public ReferencedObject
  {
    /** Class for one triplet position entry. */
    class TripletEntry
    {
    public:
      /** @name Constructor/Destructor */
      //@{
      /** Constructor. */
      TripletEntry()
      {}

      /** Destructor */
      ~TripletEntry()
      {}

      /** Dummy copy constructor.  Note that nothing is really copied!
       *  This is just implemented to that the std::list can be
       *  created with uninitialized entries.  The values are
       *  afterwards set with the Set method. */
      TripletEntry(const TripletEntry&)
      {}
      //@}

      /** Set the values of an entry */
      void Set(Index i_row, Index j_col, Index i_pos_triplet)
      {
        if (i_row>j_col) {
          i_row_ = j_col;
          j_col_ = i_row;
        }
        else {
          i_row_ = i_row;
          j_col_ = j_col;
        }
        i_pos_triplet_ = i_pos_triplet;
      }

      /** @name Accessor methods. */
      //@{
      /** Row position. */
      Index IRow() const
      {
        return i_row_;
      }
      /** Column position. */
      Index JCol() const
      {
        return j_col_;
      }
      /** Index in original triplet matrix. */
      Index PosTriplet() const
      {
        return i_pos_triplet_;
      }
      //@}

      /** Comparison operator.  This is required for the sort function. */
      bool operator< (const TripletEntry& Tentry) const
      {
        return ((i_row_ < Tentry.i_row_) ||
                (i_row_==Tentry.i_row_ && j_col_<Tentry.j_col_));
      }

    private:
      /**@name Default Compiler Generated Methods
       * (Hidden to avoid implicit creation/calling).
       * These methods are not implemented and 
       * we do not want the compiler to implement
       * them for us, so we declare them private
       * and do not define them. This ensures that
       * they will not be implicitly created/called. */
      //@{
      /** Default Constructor */
      //TripletEntry();

      /** Copy Constructor */
      /*
      TripletEntry(const TripletEntry&);
      */

      /** Overloaded Equals Operator */
      void operator=(const TripletEntry&);
      //@}

      /** @name Entry content. */
      //@{
      Index i_row_;
      Index j_col_;
      Index i_pos_triplet_;
      //@}
    };

  public:
    /** @name Constructor/Destructor */
    //@{
    /* Constructor.  If offset is 0, then the counting of indices in
       the compressed format starts a 0 (C-style numbering); if offset
       is 1, then the counting starts at 1 (Fortran-type
       numbering). */
    TripletToCSRConverter(Index offset);

    /** Destructor */
    virtual ~TripletToCSRConverter();
    //@}

    /** Initialize the converter, given the fixed structure of the
     *  matrix.  There, ndim gives the number of rows and columns of
     *  the matrix, nonzeros give the number of nonzero elements, and
     *  airn and acjn give the positions of the nonzero elements.  The
     *  return value is the number of nonzeros in the condensed
     *  matrix.  (Since nonzero elements can be listed several times
     *  in the triplet format, it is possible that this value is
     *  different from the input value nonzeros.)  This method must be
     *  called before the GetIA, GetJA, GetValues methods are called.
     */
    Index InitializeConverter(Index dim, Index nonzeros,
                              const Index* airn,
                              const Index* ajcn);

    /** Return the IA array for the condensed format. */
    const Index* IA() const
    {
      DBG_ASSERT(initialized_);
      return ia_;
    }

    /** Return the JA array for the condensed format. */
    const Index* JA() const
    {
      DBG_ASSERT(initialized_);
      return ja_;
    }

    /** Convert the values of the nonzero elements.  Given the values
     *  a_triplet for the triplet format, return the array of values
     *  for the condensed format in a_condensed. nonzeros_condensed is
     *  the length of the array a_condensed and must be identical to
     *  the return value of InitializeConverter. */
    void ConvertValues(Index nonzeros_triplet, const Number* a_triplet,
                       Index nonzeros_compressed, Number* a_compressed);

  private:
    /**@name Default Compiler Generated Methods
     * (Hidden to avoid implicit creation/calling).
     * These methods are not implemented and 
     * we do not want the compiler to implement
     * them for us, so we declare them private
     * and do not define them. This ensures that
     * they will not be implicitly created/called. */
    //@{
    /** Default Constructor */
    TripletToCSRConverter();

    /** Copy Constructor */
    TripletToCSRConverter(const TripletToCSRConverter&);

    /** Overloaded Equals Operator */
    void operator=(const TripletToCSRConverter&);
    //@}

    /** Offset for CSR numbering. */
    Index offset_;

    /** Array storing the values for IA in the condensed format */
    Index* ia_;

    /** Array storing the values for JA in the condensed format */
    Index* ja_;

    /** Dimension of the matrix. */
    Index dim_;

    /** Number of nonzeros in the triplet format. */
    Index nonzeros_triplet_;

    /** Number of nonzeros in the compressed format. */
    Index nonzeros_compressed_;

    /** Flag indicating if initialize method had been called. */
    bool initialized_;

    /** @name Arrays for cross-positions for the conversion of values. */
    //@{
    /** First elements assignement. For i with 0 <= i <=
     *  nonzeros_compressed-1, the i-th element in the compressed
     *  format is obtained from copying the ipos_filter_[i]-th element
     *  from the triplet format.  */
    Index* ipos_first_;
    /** Position of multiple elements in triplet matrix.  For i =
     *  0,..,nonzeros_triplet_-nonzeros_compressed_, the
     *  ipos_double_triplet_[i]-th element in the triplet matrix has
     *  to be added to the ipos_double_compressed_[i]-th element in
     *  the compressed matrix. */
    Index* ipos_double_triplet_;
    /** Position of multiple elements in compressed matrix. */
    Index* ipos_double_compressed_;
    //@}
  };


} // namespace Ipopt

#endif
