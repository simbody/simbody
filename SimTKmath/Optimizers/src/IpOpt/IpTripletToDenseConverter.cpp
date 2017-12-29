
#include "IpTripletToDenseConverter.hpp"
#include <list>

namespace SimTKIpopt
{
#ifdef IP_DEBUG
  static const Index dbg_verbosity = 0;
#endif

  TripletToDenseConverter::TripletToDenseConverter(Index offset)
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

  TripletToDenseConverter::~TripletToDenseConverter()
  {
    delete[] ia_;
    delete[] ja_;
    delete[] ipos_first_;
    delete[] ipos_double_triplet_;
    delete[] ipos_double_compressed_;
  }

  Index TripletToDenseConverter::InitializeConverter(Index dim, Index nonzeros,
      const Index* airn,
      const Index* ajcn)
  {
    int i;
    DBG_START_METH("TripletToDenseConverter::InitializeConverter",
                   dbg_verbosity);

    DBG_ASSERT(dim>0);
    DBG_ASSERT(nonzeros>0);

//printf("TripletToDenseConverter::InitializeConverter");
    ia_ = new Index[nonzeros];
    ja_ = new Index[nonzeros];

//printf("TSymLinearSolver::InitializeStructure indexes =\n");
    for(i=0;i<nonzeros;i++ ) {
       ia_[i] = airn[i] - 1;
       ja_[i] = ajcn[i] - 1;
//    printf("%d %d \n",ia_[i],ja_[i]);
    }

    return nonzeros;
  }

  void TripletToDenseConverter::ConvertValues(Index nonzeros_triplet,
      const Number* a_triplet,
      Index dim,
      Number* a)
  {
    int i,j;
    DBG_START_METH("TSymLinearSolver::ConvertValues",
                   dbg_verbosity);
//printf(" TripletToDenseConverter:::ConvertValues\n");

    DBG_ASSERT(initialized_);

    DBG_ASSERT(nonzeros_triplet_==nonzeros_triplet);
    DBG_ASSERT(nonzeros_compressed_==nonzeros_compressed);

    /* initialize matrix to all zeros */
    for (i=0;i<dim;i++ ) {
       for (j=0;j<dim;j++ ) {
         a[i*dim+j] = 0.0;
       }
    }
/*
printf("TripletToDenseConverter::ConvertValues input = \n");     
 for (Index i=0; i<nonzeros_triplet; i++) {
  printf(" %d %d %f \n",ia_[i],ja_[i],a_triplet[i]);
}
*/

    /* fill in non zeros assuming symmetric input matrix */
    for (Index i=0; i<nonzeros_triplet; i++) {
//       printf("dim = %d value = %f index = %d %d \n",dim, a_triplet[i], ia_[i]*dim+ja_[i], ja_[i]*dim+ia_[i]);
       a[ja_[i]*dim+ia_[i]] += a_triplet[i]; // lower triangle
       if(ia_[i] != ja_[i] ) { // only fill in off diagonal elements 
          a[ia_[i]*dim+ja_[i]] += a_triplet[i]; // upper triangle
       }
    }
/*
printf("TripletToDenseConverter::ConvertValues  a matrix : \n");
    for (i=0;i<dim;i++ ) {
       for (j=0;j<dim;j++ ) {
          printf("%f ",a[i*dim+j] );
       }
       printf("\n");
    }
*/


  }

} // namespace Ipopt
