
#include "IpTripletToCSRConverter.hpp"
#include <list>

namespace SimTKIpopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  TripletToCSRConverter::TripletToCSRConverter(Index offset)
      :
      offset_(offset),
      ia_(NULL),
      ja_(NULL),
      dim_(0),
      nonzeros_triplet_(0),
      nonzeros_compressed_(0),
      initialized_(false),
      ipos_first_(NULL),
      ipos_double_triplet_(NULL),
      ipos_double_compressed_(NULL)
  {
    DBG_ASSERT(offset==0|| offset==1);
  }

  TripletToCSRConverter::~TripletToCSRConverter()
  {
    delete[] ia_;
    delete[] ja_;
    delete[] ipos_first_;
    delete[] ipos_double_triplet_;
    delete[] ipos_double_compressed_;
  }

  Index TripletToCSRConverter::InitializeConverter(Index dim, Index nonzeros,
      const Index* airn,
      const Index* ajcn)
  {
    DBG_START_METH("TSymLinearSolver::InitializeStructure",
                   dbg_verbosity);

    DBG_ASSERT(dim>0);
    DBG_ASSERT(nonzeros>0);

    delete[] ia_;
    delete[] ja_;
    delete[] ipos_first_;
    delete[] ipos_double_triplet_;
    delete[] ipos_double_compressed_;

    dim_ = dim;
    nonzeros_triplet_ = nonzeros;

    // Create a list with all triplet entries
    std::list<TripletEntry> entry_list(nonzeros);
    std::list<TripletEntry>::iterator list_iterator = entry_list.begin();
    for (Index i=0; i<nonzeros; i++) {
      list_iterator->Set(airn[i], ajcn[i], i);
      ++list_iterator;
    }
    DBG_ASSERT(list_iterator == entry_list.end());

    if (DBG_VERBOSITY()>=2) {
      for (Index i=0; i<nonzeros; i++) {
        DBG_PRINT((2, "airn[%5d] = %5d acjn[%5d] = %5d\n", i, airn[i], i, ajcn[i]));
      }
    }

    // sort the list
    entry_list.sort();

    // Now got through the list and compute ipos_ arrays and the
    // number of elements in the compressed format
    Index* ja_tmp = new Index[nonzeros];    // overestimate memory requirement
    ia_ = new Index[dim_+1];
    Index* ipos_first_tmp = new Index[nonzeros];  // overestimate memory requirement
    Index* ipos_double_triplet_tmp = new Index[nonzeros];  // overestimate memory requirement
    Index* ipos_double_compressed_tmp = new Index[nonzeros];  // overestimate memory requirement

    nonzeros_compressed_ = 0;
    Index cur_row = 1;

    // The first element must be the first diagonal element
    list_iterator = entry_list.begin();
    DBG_ASSERT(list_iterator->IRow()==1);
    DBG_ASSERT(list_iterator->JCol()==1);
    ia_[0] = 0;
    ja_tmp[0] = 1;
    ipos_first_tmp[0] = list_iterator->PosTriplet();

    ++list_iterator;
    Index idouble = 0;
    while (list_iterator != entry_list.end()) {
      Index irow = list_iterator->IRow();
      Index jcol = list_iterator->JCol();
      if (cur_row == irow && ja_tmp[nonzeros_compressed_] == jcol) {
        // This element appears repeatedly, add to the double list
        ipos_double_triplet_tmp[idouble] = list_iterator->PosTriplet();
        ipos_double_compressed_tmp[idouble] = nonzeros_compressed_;
        idouble++;
      }
      else {
        // This is a new element
        nonzeros_compressed_++;
        ja_tmp[nonzeros_compressed_] = jcol;
        ipos_first_tmp[nonzeros_compressed_] = list_iterator->PosTriplet();
        if (cur_row != irow) {
          // this is in a new row

          // make sure that the diagonal element is given and that no
          // row is omitted
          DBG_ASSERT(irow==jcol);
          DBG_ASSERT(cur_row+1==irow);
          ia_[cur_row] = nonzeros_compressed_;
          cur_row++;
        }
      }

      ++list_iterator;
    }
    nonzeros_compressed_++;
    ia_[dim_] = nonzeros_compressed_;

    DBG_ASSERT(cur_row == dim_);
    DBG_ASSERT(idouble == nonzeros_triplet_-nonzeros_compressed_);

    // Now copy the ja_tmp array to the (shorter) final one and make
    // sure that the correct offset is used
    ja_ = new Index[nonzeros_compressed_];
    if (offset_==0) {
      for (Index i=0; i<nonzeros_compressed_; i++) {
        ja_[i] = ja_tmp[i] - 1;
      }
    }
    else {
      for (Index i=0; i<nonzeros_compressed_; i++) {
        ja_[i] = ja_tmp[i];
      }
      for (Index i=0; i<=dim_; i++) {
        ia_[i] = ia_[i] + 1;
      }
    }
    delete[] ja_tmp;

    // Reallocate memory for the "first" array
    ipos_first_ = new Index[nonzeros_compressed_];
    for (Index i=0; i<nonzeros_compressed_; i++) {
      ipos_first_[i] = ipos_first_tmp[i];
    }
    delete[] ipos_first_tmp;

    // Reallocate memory for the "double" arrays
    ipos_double_triplet_ = new Index[idouble];
    ipos_double_compressed_ = new Index[idouble];
    for (Index i=0; i<idouble; i++) {
      ipos_double_triplet_[i] = ipos_double_triplet_tmp[i];
      ipos_double_compressed_[i] = ipos_double_compressed_tmp[i];
    }
    delete[] ipos_double_triplet_tmp;
    delete[] ipos_double_compressed_tmp;

    initialized_ = true;

    if (DBG_VERBOSITY()>=2) {
      for (Index i=0; i<=dim_; i++) {
        DBG_PRINT((2, "ia[%5d] = %5d\n", i, ia_[i]));
      }
      for (Index i=0; i<nonzeros_compressed_; i++) {
        DBG_PRINT((2, "ja[%5d] = %5d ipos_first[%5d] = %5d\n", i, ja_[i], i, ipos_first_[i]));
      }
      for (Index i=0; i<nonzeros_triplet_-nonzeros_compressed_; i++) {
        DBG_PRINT((2, "ipos_double_triplet[%5d] = %5d ipos_double_compressed[%5d] = %5d\n", i, ipos_double_triplet_[i], i, ipos_double_compressed_[i]));
      }
    }

    return nonzeros_compressed_;
  }

  void TripletToCSRConverter::ConvertValues(Index nonzeros_triplet,
      const Number* a_triplet,
      Index nonzeros_compressed,
      Number* a_compressed)
  {
    DBG_START_METH("TSymLinearSolver::ConvertValues",
                   dbg_verbosity);

    DBG_ASSERT(initialized_);

    DBG_ASSERT(nonzeros_triplet_==nonzeros_triplet);
    DBG_ASSERT(nonzeros_compressed_==nonzeros_compressed);

    for (Index i=0; i<nonzeros_compressed_; i++) {
      a_compressed[i] = a_triplet[ipos_first_[i]];
    }
    for (Index i=0; i<nonzeros_triplet_-nonzeros_compressed_; i++) {
      a_compressed[ipos_double_compressed_[i]] +=
        a_triplet[ipos_double_triplet_[i]];
    }

    if (DBG_VERBOSITY()>=2) {
      for (Index i=0; i<nonzeros_triplet; i++) {
        DBG_PRINT((2, "atriplet[%5d] = %24.16e\n", i, a_triplet[i]));
      }
      for (Index i=0; i<nonzeros_compressed; i++) {
        DBG_PRINT((2, "acompre[%5d] = %24.16e\n", i, a_compressed[i]));
      }
    }
  }

} // namespace Ipopt
