#ifdef NOTDEF


#include "mmapAlloc.h"

#include <cdsMath.h>
#include <cdsSStream.h>
#include <cdsExcept.h>
//#include <unistd.h>
//#include <sys/mman.h> 
#include <sys/types.h>
#include <errno.h>
#include <fcntl.h>
#include <string.h>


static bool debugOutput=0;
static int totalAllocated=0;
static int callCnt=0;
static int offsetSize=0;

#ifndef MAP_ANONYMOUS
#  define MAP_ANONYMOUS 0
#endif

void*
MMapAlloc::alloc(size_t size,
		 bool   deleteFile,
		 int    mapID)
{
 static int count=0;

 totalAllocated += size;
 callCnt++;
 if ( debugOutput )
   cout << "MMapAlloc::alloc: size requested: " << size << endl;
 if ( debugOutput )
   cout << "MMapAlloc::alloc: total requests: " << callCnt
	<< " total amount: " << totalAllocated << endl;
 using CDS::exception;
 String routineName = "MMapAlloc::alloc";

 OStringStream mmapFilename;
 char* logName = getlogin();

 // this should ensure that the returned pointer is properly aligned.
 // it might be better to make this platform-specific.
 if ( offsetSize==0 )
   offsetSize = CDSMath::max(sizeof(double) , sizeof(size_t));
 size += offsetSize;
 

 if ( mapID>=0 )
   // NOTE: this assumes that everyone calling this is in the same
   // process group.
   mmapFilename << "/var/tmp/mmapAlloc-" << (logName?logName:"unknown") << '-'
		<< getpgrp() << "-manual-" << mapID << ends;
 else {
   mmapFilename << "/var/tmp/mmapAlloc-" << (logName?logName:"unknown") << '-'
		<< getpid() << '-' << count << ends;
   count++;
 }

 if ( debugOutput )
   cout << "MMapAlloc::alloc: deleteFile: " << deleteFile
	<< "\n\tfilename: " << mmapFilename.str() << '\n';

 bool alreadyExists=0;
 int mmapFlags = MAP_SHARED | X_MMAP_FLAGS;
 if ( mapID>=0 )
   // don't know how to use anonymous map between processes
   mmapFlags &= ~(int)MAP_ANONYMOUS; 
		 
 int fd = -1;
 if ( ! (mmapFlags & MAP_ANONYMOUS) ) {
   if ( (fd = open(mmapFilename.str(),O_RDWR,0600)) != -1 )
     alreadyExists=1;
   else 
     fd = open(mmapFilename.str(),O_CREAT|O_RDWR,0600);

   if ( fd==-1 )
     throw exception(routineName + ": error in open():\n\t" +
		     strerror(errno));


   if ( deleteFile )
     unlink(mmapFilename.str());

   if ( !alreadyExists && 
	write(fd,&size,sizeof(size_t)) != sizeof(size_t) )
     throw exception(routineName + ": error in write():\n\t" +
		     strerror(errno));
 }

 //  //FIX: shouldn't be necessary
 //  for (int i=0 ; i<size ; i++)
 //    write(fd,"\0",1);


 if ( !alreadyExists &&
      fd>=0 && 
      (lseek(fd,size-1,SEEK_SET) != (off_t)(size-1)) )
   throw exception(routineName + ": error in lseek():\n\t" +
		   strerror(errno));

 if ( !alreadyExists &&
      fd>=0 && (write(fd,"\0",1) != 1) )
   throw exception(routineName + ": error in second write():\n\t" +
		   strerror(errno));

 if ( debugOutput )
   cout << "file descriptor: " << fd << endl;
 char* ret = (char*)mmap(NULL,size,PROT_READ|PROT_WRITE,
			 mmapFlags,fd,(off_t)0);
 
 if ( debugOutput )
   cout << "MMapAlloc::alloc: mmap returns: " << (void*)ret << endl;

 if ( ret == MAP_FAILED ) {
    if ( debugOutput )
      cout << "MMapAlloc::alloc: errno: " << errno << endl;
    throw exception(routineName + ": error in mmap():\n\t" +
		    strerror(errno));
 }

 if ( alreadyExists ) {
   if ( *((size_t*)ret) != size ) {
     cout << "MMapAlloc::alloc: allocation mismatch: " 
	  << *((size_t*)ret) << " != " << size <<endl;
     throw exception(routineName + ": size allocation mismatch.");
   } 
 } else
   *((size_t*)ret) = size;

 if ( debugOutput )
   cout << "MMapAlloc::alloc: first byte: " << *ret << "\n\t"
	<< "file descriptor: " << fd << endl;

 if ( fd>=0 && (close(fd) == -1) ) 
   throw exception(routineName + ": error in close():\n\t" +
		     strerror(errno));

 ret += offsetSize;

 if ( debugOutput )
   cout << "MMapAlloc::alloc: returns: " << (void*)ret << endl;
 return ret;
} /* alloc */

void
MMapAlloc::free(void* p)
{
 size_t* actualP = (size_t*)((char*)p - offsetSize);
 size_t size = *actualP;
 if ( debugOutput )
   cout << "MMapAlloc::free: unmapping " << size << " bytes at "
	<< actualP << endl;
 if ( munmap((char*)actualP,size)<0 )
   throw CDS::exception(String("MMapAlloc::free: error in munmap():\n\t") +
			strerror(errno));
} /* free */

#endif
