#ifdef NOTDEF

//
// my definitions of global new and delete
//

//
// these print debug messages using iostreams: if iostreams use ::new
//  there will be an infinite loop....
//

//
// things to know: 
//
//  new must always return a unique pointer, even when given a size-zero
//  argument
//
//  delete must work when given a NULL-valued pointer (i.e. do nothing)
//
//

#include <sthead.h>

#include <stdio.h>



#if DEC_COMPILERS

#include <new>                                          
                                                        
// Redefine global operator new(),                      
// entry point into C++ Standard Library.               
// Also override user calls to operator new()           
                                                        
//void *operator new(size_t size) throw(std::bad_alloc) { 
//printf("in my global new\n");                           
//...                                                     
//}                                                       
void *operator new(size_t size) throw(std::bad_alloc)
     //global allocator
{
 void *p;

 if (size==0)
   size=1; // new must return a unique pointer

 p = malloc(size);
 if (p == NULL) {
   printf("new: ERROR: Malloc failed on allocation size of %d\n",size);
   
   abort();
 }

 printf("::new: allocation at %x of size %d\n",p,size);
 return p;
}
                                                        
// Entry point into C++ Class Library                   
                                                        
#ifdef __MODEL_ANSI                                     
extern "C" void *__7__nw__FUl(size_t size) {            
#else  //__MODEL_ARM                                    
extern "C" void *__nw__XUl(size_t size) {               
#endif                                                  
     return ::operator new(size);                       
}                                                       

void operator delete(void *p)
     //global deallocator
{
 printf("::delete: %x\n",p);
 if (p!=NULL)
   free(p);
 p = NULL;
}

void *operator new[](size_t size)
     //global array allocator
{
 void *p;

 if (size==0) 
   size = 1;

 p = malloc(size);
 if (p == NULL) {
   printf("new[]: ERROR: Malloc failed on allocation size of %d\n",size);


   abort();
 }
 printf("::new[]: allocation at %x of size %d\n",p,size);
 return p;
}
   
void operator delete[](void *p)
     //global deallocator
{
 printf("::delete[]: %x\n",p);
 if (p!=NULL)
   free(p);
 p = NULL;
}

#else  /* not DEC compilers */

#if MIPSEB && USE_TRACEBACK
#     include <libexc.h>     // for traceback
#endif

void *operator new(size_t size)
     //global allocator
{
 void *p;

 if (size==0)
   size=1; // new must return a unique pointer

 p = malloc(size);
 if (p == NULL) {
   printf("new: ERROR: Malloc failed on allocation size of %d\n",size);
   
#if MIPSEB && USE_TRACEBACK
   trace_back_stack_and_print();
#endif

   abort();
 }

 printf("::new: allocation at %x of size %d\n",p,size);
 return p;
}
   
void operator delete(void *p)
     //global deallocator
{
 printf("::delete: %x\n",p);
 if (p!=NULL)
   free(p);
 p = NULL;
}





#ifndef AIX
//something wierd here with AIX - doesn't seem to understand operator new[]
//SGI can't redefined the global allocator???


void *operator new[](size_t size)
     //global array allocator
{
 void *p;

 if (size==0) 
   size = 1;

 p = malloc(size);
 if (p == NULL) {
   printf("new[]: ERROR: Malloc failed on allocation size of %d\n",size);

#if MIPSEB && USE_TRACEBACK
   trace_back_stack_and_print();
#endif

   abort();
 }
 printf("::new[]: allocation at %x of size %d\n",p,size);
 return p;
}
   
void operator delete[](void *p)
     //global deallocator
{
 printf("::delete[]: %x\n",p);
 if (p!=NULL)
   free(p);
 p = NULL;
}
#endif /*not AIX*/

#endif /* not DEC compilers */

#endif
