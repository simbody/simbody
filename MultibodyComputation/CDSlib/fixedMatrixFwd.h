// This file is necessary because FixedMatrix is a template class
// with default template arguments. The defaults must appear first
// in a compilation unit (sherm 050528)


#ifndef __fixedMatrixFwd_h__
#define __fixedMatrixFwd_h__ 1

template<class T,int size1,int size2=size1,int OFFSET1=0,int OFFSET2=OFFSET1>
class FixedMatrix;


#endif /* __fixedMatrixFwd_h__ */
