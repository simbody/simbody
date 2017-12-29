// Copyright (C) 2004, 2006 International Business Machines and others.
// All Rights Reserved.
// This code is published under the Common Public License.
//
// $Id: IpFilter.cpp 765 2006-07-14 18:03:23Z andreasw $
//
// Authors:  Carl Laird, Andreas Waechter     IBM    2004-08-13

#include "IpFilter.hpp"
#include "IpJournalist.hpp"

namespace SimTKIpopt
{

#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  ///////////////////////////////////////////////////////////////////////////
  //                            Filter entries                             //
  ///////////////////////////////////////////////////////////////////////////

  FilterEntry::FilterEntry(std::vector<Number> vals, Index iter)
      :
      vals_(vals),
      iter_(iter)
  {}

  FilterEntry::~FilterEntry()
  {}

  ///////////////////////////////////////////////////////////////////////////
  //                                 Filter                                //
  ///////////////////////////////////////////////////////////////////////////

  Filter::Filter(Index dim)
      :
      dim_(dim)
  {}

  bool Filter::Acceptable(std::vector<Number> vals) const
  {
    DBG_START_METH("FilterLineSearch::Filter::Acceptable", dbg_verbosity);
    DBG_ASSERT((Index)vals.size()==dim_);
    bool acceptable = true;
    std::list<FilterEntry*>::iterator iter;
    for (iter = filter_list_.begin(); iter != filter_list_.end();
         ++iter) {
      if (!(*iter)->Acceptable(vals)) {
        acceptable = false;
        break;
      }
    }
    return acceptable;
  }

  void Filter::AddEntry(std::vector<Number> vals, Index iteration)
  {
    DBG_START_METH("FilterLineSearch::Filter::AddEntry", dbg_verbosity);
    DBG_ASSERT((Index)vals.size()==dim_);
    std::list<FilterEntry*>::iterator iter;
    iter = filter_list_.begin();
    while (iter != filter_list_.end()) {
      if ((*iter)->Dominated(vals)) {
        std::list<FilterEntry*>::iterator iter_to_remove = iter;
        ++iter;
        FilterEntry* entry_to_remove = *iter_to_remove;
        filter_list_.erase(iter_to_remove);
        delete entry_to_remove;
      }
      else {
        ++iter;
      }
    }
    FilterEntry* new_entry = new FilterEntry(vals, iteration);
    filter_list_.push_back(new_entry);
  }

  void Filter::Clear()
  {
    DBG_START_METH("FilterLineSearch::Filter::Clear", dbg_verbosity);
    while (!filter_list_.empty()) {
      FilterEntry* entry = filter_list_.back();
      filter_list_.pop_back();
      delete entry;
    }
  }

  void Filter::Print(const Journalist& jnlst)
  {
    DBG_START_METH("FilterLineSearch::Filter::Print", dbg_verbosity);
    jnlst.Printf(J_DETAILED, J_LINE_SEARCH,
                 "The current filter has %d entries.\n", filter_list_.size());
    if (!jnlst.ProduceOutput(J_VECTOR, J_LINE_SEARCH)) {
      return;
    }
    std::list<FilterEntry*>::iterator iter;
    Index count = 0;
    for (iter = filter_list_.begin(); iter != filter_list_.end();
         ++iter) {
      if (count % 10 == 0) {
        jnlst.Printf(J_VECTOR, J_LINE_SEARCH,
                     "                phi                    theta            iter\n");
      }
      count++;
      jnlst.Printf(J_VECTOR, J_LINE_SEARCH, "%5d ", count);
      for (Index i=0; i<dim_; i++) {
        jnlst.Printf(J_VECTOR, J_LINE_SEARCH, "%23.16e ", (*iter)->val(i));
      }
      jnlst.Printf(J_VECTOR, J_LINE_SEARCH, "%5d\n",(*iter)->iter());
    }
  }

} // namespace Ipopt
