// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpFilter.hpp 759 2006-07-07 03:07:08Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#ifndef __IPFILTER_HPP__
#define __IPFILTER_HPP__

#include "IpJournalist.hpp"
#include "IpDebug.hpp"
#include <list>
#include <vector>

namespace Ipopt
{

  /** Class for one filter entry. */
  class FilterEntry
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Constructor with the two components and the current iteration count */
    FilterEntry(std::vector<Number> vals, Index iter);

    /** Default Destructor */
    ~FilterEntry();
    //@}

    /** Check acceptability of pair (phi,theta) with respect
     *  to this filter entry.  Returns true, if pair is acceptable.
     */
    bool Acceptable(std::vector<Number> vals) const
    {
      Index ncoor = (Index)vals_.size();
      DBG_ASSERT((Index)vals.size() == ncoor);

      // ToDo decide if we need Compare_le
      bool retval = false;
      for (Index i=0; i<ncoor; i++) {
        if (vals[i] <= vals_[i]) {
          retval = true;
          break;
        }
      }

      return retval;
    }

    /** Check if this entry is dominated by given coordinates.
     *  Returns true, if this entry is dominated.
     */
    bool Dominated(std::vector<Number> vals) const
    {
      Index ncoor = (Index)vals_.size();
      DBG_ASSERT((Index)vals.size() == ncoor);

      bool retval = true;
      for (Index i=0; i<ncoor; i++) {
        if (vals[i] > vals_[i]) {
          retval = false;
          break;
        }
      }

      return retval;
    }

    /** @name Accessor functions */
    //@{
    Number val(Index i) const
    {
      return vals_[i];
    }
    Index iter() const
    {
      return iter_;
    }
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
    /** Default Constructor */
    FilterEntry();
    /** Copy Constructor */
    FilterEntry(const FilterEntry&);

    /** Overloaded Equals Operator */
    void operator=(const FilterEntry&);
    //@}

    /** values defining the coordinates of the entry */
    std::vector<Number> vals_;
    /** iteration number in which this entry was added to filter */
    const Index iter_;
  };

  /** Class for the filter.  This class contains all filter entries.
   *  The entries are stored as the corner point, including the
   *  margin. */
  class Filter
  {
  public:
    /**@name Constructors/Destructors */
    //@{
    /** Default Constructor */
    Filter(Index dim);
    /** Default Destructor */
    ~Filter()
    {
      //ToDo figure out if that here is necessary
      Clear();
    }
    //@}

    /** Check acceptability of given coordinates with respect
     *  to the filter.  Returns true, if pair is acceptable
     */
    bool Acceptable(std::vector<Number> vals) const;

    /** Add filter entry for given coordinates.  This will also
     *  delete all dominated entries in the current filter. */
    void AddEntry(std::vector<Number> vals, Index iteration);

    /** @name Wrappers for 2-dimensional filter. */
    //@{
    bool Acceptable(Number val1, Number val2) const
    {
      std::vector<Number> vals(2);
      vals[0] = val1;
      vals[1] = val2;

      return Acceptable(vals);
    }

    void AddEntry(Number val1, Number val2, Index iteration)
    {
      std::vector<Number> vals(2);
      vals[0] = val1;
      vals[1] = val2;

      AddEntry(vals, iteration);
    }
    //@}

    /** Delete all filter entries */
    void Clear();

    /** Print current filter entries */
    void Print(const Journalist& jnlst);

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
    Filter();
    /** Copy Constructor */
    Filter(const Filter&);

    /** Overloaded Equals Operator */
    void operator=(const Filter&);
    //@}

    /** Dimension of the filter (number of coordinates per entry) */
    Index dim_;

    /** List storing the filter entries */
    mutable std::list<FilterEntry*> filter_list_;
  };

} // namespace Ipopt

#endif
