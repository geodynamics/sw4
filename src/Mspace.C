#include "Mspace.h"
#include <unordered_map>
#include "caliper.h"
struct global_variable_holder_struct global_variables = { .gpu_memory_hwm=0 , .curr_mem=0, .max_mem = 0 };
using namespace std;

typedef struct {
  const char *file;
  int line;
  size_t size;
  Space type;} pattr_t;

pattr_t *patpush(void *ptr, pattr_t *ss){
  static std::unordered_map<void*,pattr_t *> map;
  if (ss!=NULL) {
    map[ptr]=ss;
  } else {
    std::unordered_map<void*,pattr_t*>::const_iterator got = map.find (ptr);
    if (got==map.end()){
      //std:cerr<<"ELEMENT NOT FOUND IN MAP\n";
      return NULL;
    } else
      return got->second;
  }
  return NULL;
}
  

void check_mem(){
  size_t mfree,mtotal;
  SW4_CheckDeviceError(cudaMemGetInfo(&mfree,&mtotal));
  global_variables.gpu_memory_hwm=std::max((mtotal-mfree), global_variables.gpu_memory_hwm);
}
  
void * operator new(std::size_t size,Space loc) throw(std::bad_alloc){
#ifdef ENABLE_CUDA
if (loc==Managed){
  //std::cout<<"Managed allocation \n";
    if (size==0) size=1; // new has to return an valid pointer for 0 size.
    void *ptr;
    if (cudaMallocManaged(&ptr,size)!=cudaSuccess){
      std::cerr<<"Mananged memory allocation failed "<<size<<"\n";
      throw std::bad_alloc();
    } else {
      check_mem();
      global_variables.curr_mem+=size;
      global_variables.max_mem=std::max(global_variables.max_mem,global_variables.curr_mem);
      return ptr;
    }
    
  } else if (loc==Host){
  //std::cout<<"Calling my placement new \n";
    return ::operator new(size);
 } else if (loc==Device){
  //std::cout<<"Managed allocation \n";
    if (size==0) size=1; // new has to return an valid pointer for 0 size.
    void *ptr;
    if (cudaMalloc(&ptr,size)!=cudaSuccess){
      std::cerr<<"Device memory allocation failed "<<size<<"\n";
      throw std::bad_alloc();
    } else return ptr;
 } else if (loc==Pinned){ 
  if (size==0) size=1; // new has to return an valid pointer for 0 size.
  void *ptr;
  SW4_CheckDeviceError(cudaHostAlloc(&ptr,size,cudaHostAllocMapped));
  return ptr;
 }else {
  std::cerr<<"Unknown memory space for allocation request "<<loc<<"\n";
    throw std::bad_alloc();
  }
#else
 if ((loc==Managed)||(loc==Device)||(loc==Pinned)){
    //std::cout<<"Managed location not available yet \n";
    return ::operator new(size);
  } else if (loc==Host){
    //std::cout<<"Calling my placement new \n";
    return ::operator new(size);
  } else {
    std::cerr<<"Unknown memory space for allocation request\n";
    throw std::bad_alloc();
  }
#endif
}
void * operator new(std::size_t size,Space loc,char *file, int line) throw(std::bad_alloc){
  std::cout<<"Calling tracking new from "<<line<<" of "<<file<<"\n";
  pattr_t *ss=new pattr_t;
  ss->file=file;
  ss->line=line;
  ss->type=loc;
  ss->size=size;
  void *ret= ::operator new(size,loc);
  patpush(ret,ss);
  return ret;
}

void * operator new[](std::size_t size,Space loc) throw(std::bad_alloc){
#ifdef ENABLE_CUDA
if (loc==Managed){
  //std::cout<<"Managed [] allocation \n";
    if (size==0) size=1; // new has to return an valid pointer for 0 size.
    void *ptr;
    if (cudaMallocManaged(&ptr,size)!=cudaSuccess){
      std::cerr<<"Managed memory allocation failed "<<size<<"\n";
      throw std::bad_alloc();
    } else {
      check_mem();
      global_variables.curr_mem+=size;
      global_variables.max_mem=std::max(global_variables.max_mem,global_variables.curr_mem);
      return ptr;
    }
    
  } else if (loc==Host){
  // std::cout<<"Calling my placement new \n";
    return ::operator new(size);
  } else if (loc==Device){
  //std::cout<<"Managed allocation \n";
    if (size==0) size=1; // new has to return an valid pointer for 0 size.
    void *ptr;
    if (cudaMalloc(&ptr,size)!=cudaSuccess){
      std::cerr<<"Device memory allocation failed "<<size<<"\n";
      throw std::bad_alloc();
    } else return ptr;
 }else if (loc==Pinned){ 
  if (size==0) size=1; // new has to return an valid pointer for 0 size.
  void *ptr;
  SW4_CheckDeviceError(cudaHostAlloc(&ptr,size,cudaHostAllocMapped));
  return ptr;
 }  else {
  //cudaHostAlloc(&ptr,size+sizeof(size_t)*MEM_PAD_LEN,cudaHostAllocMapped));
  std::cerr<<"Unknown memory space for allocation request "<<loc<<"\n";
    throw std::bad_alloc();
  }
#else
 if ((loc==Managed)||(loc==Device)||(loc==Pinned)){
    //std::cout<<"Managed location not available yet \n";
    return ::operator new(size);
  } else if (loc==Host){
    //std::cout<<"Calling my placement new \n";
    return ::operator new(size);
  } else {
  std::cerr<<"Unknown memory space for allocation request "<<loc<<"\n";
    throw std::bad_alloc();
  }
#endif
}
void * operator new[](std::size_t size,Space loc,const char *file,int line){
  //std::cout<<"Calling tracking new from "<<line<<" of "<<file<<"\n";
  pattr_t *ss=new pattr_t;
  ss->file=file;
  ss->line=line;
  ss->type=loc;
  ss->size=size;
  void *ret= ::operator new(size,loc);
  patpush(ret,ss);
  return ret;
}

void operator delete(void *ptr, Space loc) throw(){
#ifdef ENABLE_CUDA
  if ((loc==Managed)||(loc==Device)){
    //std::cout<<"Managed delete\n";
    pattr_t *ss = patpush(ptr,NULL);
    if (ss!=NULL){
      global_variables.curr_mem-=ss->size;
      //global_variables.max_mem=std::max(global_variables.max_mem,global_variables.curr_mem);
    }
    SW4_CheckDeviceError(cudaFree(ptr));
  } else if (loc==Pinned)
    SW4_CheckDeviceError(cudaFreeHost(ptr)); 
  else if (loc==Host){
    //std:cout<<"Calling my placement delete\n";
    ::operator delete(ptr);
  } else {
    std::cerr<<"Unknown memory space for de-allocation request\n";
  }
#else
  if ((loc==Managed)||(loc==Device)){
    //std::cout<<"Managed delete not available yet \n";
    ::operator delete(ptr);
  } else if (loc==Host){
    //std:cout<<"Calling my placement delete\n";
    ::operator delete(ptr);
  } else {
    std::cerr<<"Unknown memory space for de-allocation request "<<loc<<"\n";
  }
#endif
}

void operator delete[](void *ptr, Space loc) throw(){
#ifdef ENABLE_CUDA
  if ((loc==Managed)||(loc==Device)){
    //std::cout<<"Managed [] delete\n";
    pattr_t *ss = patpush(ptr,NULL);
    if (ss!=NULL){
      global_variables.curr_mem-=ss->size;
      //global_variables.max_mem=std::max(global_variables.max_mem,curr_mem);
    }
    SW4_CheckDeviceError(cudaFree(ptr));
    
  }else if (loc==Pinned)
    SW4_CheckDeviceError(cudaFreeHost(ptr));
  else if (loc==Host){
    //std:cout<<"Calling my placement delete\n";
    ::operator delete(ptr);
  } else {
    std::cerr<<"Unknown memory space for de-allocation request "<<loc<<"\n";
  }
#else
    if ((loc==Managed)||(loc==Device)||(loc==Pinned)){
    //std::cout<<"Managed delete not available yet \n";
    ::operator delete(ptr);
  } else if (loc==Host){
    //std:cout<<"Calling my placement delete\n";
    ::operator delete(ptr);
  } else {
      std::cerr<<"Unknown memory space for de-allocation request "<<loc<<"\n";
  }
#endif
}

#if defined(SW4_TRACK_MEMORY_ALLOCATIONS)
void assert_check_managed(void *ptr, const char *file, int line){
  if (ptr==NULL) return;
  pattr_t *ss = patpush(ptr,NULL);
  if (ss!=NULL)
    {
      if (ss->type!=Managed)
	{
	  std::cerr<<"ASSERT_MANAGED FAILURE in line"<<line<<" of file "<<file<<"type = "<<ss->type<<"pointer  ="<<ptr<<"\n";
	  std::cerr<<"ASSERT_MANAGED failed on allocation from line"<<ss->line<<"of "<<ss->file<<"\n";
	}
    }
  else
    {
      std::cerr<<"ASSERT_MANAGED FAILURE in line"<<line<<" of file "<<file<<"type = "<<ss->type<<"pointer  ="<<ptr<<"\n";
      std::cerr<<"No info in map\n Use PointerAttributes\n";
      //printf("Address not in map\n Calling PrintPointerAttributes\n");
      //if ( PointerAttributes(ptr)!=HYPRE_MANAGED_POINTER){
	//fprintf(stderr,"ASSERT_MANAGED FAILURE in line %d of file %s \n NO ALLOCATION INFO\n",line,file);
      //}
    }
  
}
void assert_check_host(void *ptr, const char *file, int line){
  if (ptr==NULL) return;
  pattr_t *ss = patpush(ptr,NULL);
  if (ss!=NULL)
    {
      if (ss->type!=Host)
	{
	  std::cerr<<"ASSERT_HOST FAILURE in line"<<line<<" of file "<<file<<"type = "<<ss->type<<"pointer  ="<<ptr<<"\n";
	  std::cerr<<"ASSERT_HOST failed on allocation from line"<<ss->line<<"of "<<ss->file<<"\n";
	  
	}
    }
  else
    {
      std::cerr<<"ASSERT_HOST FAILURE in line"<<line<<" of file "<<file<<"type = "<<ss->type<<"pointer  ="<<ptr<<"\n";
      std::cerr<<"Ptr not in map\n Use PointerAttribs or somesuch\n";
      //printf("Address not in map\n Calling PrintPointerAttributes\n");
      // if ( PointerAttributes(ptr)!=HYPRE_HOST_POINTER){
      // 	fprintf(stderr,"ASSERT_HOST FAILURE in line %d of file %s \n NO ALLOCATION INFO\n",line,file);
      // }
    }
  
}
void ptr_push(void *ptr,Space type, size_t size,const char *file,int line){
  pattr_t *ss=new pattr_t;
  ss->file=file;
  ss->line=line;
  ss->type=type;
  ss->size=size;
  patpush(ptr,ss);
  return;
}

size_t getsize(void *ptr){
  if (ptr==NULL) return -1;
  pattr_t *ss = patpush(ptr,NULL);
  if (ss!=NULL) return ss->size;
  else return -1;
}
#endif

#if defined(ENABLE_CUDA)

void prefetch_to_device(const float_sw4 *ptr){
  if (ptr==NULL) return;
  pattr_t *ss = patpush((void*)ptr,NULL);
  if (ss!=NULL) {
    if (ss->size>0){
      SW4_MARK_BEGIN(" prefetch_to_device");
      SW4_CheckDeviceError(cudaMemPrefetchAsync(ptr,
						ss->size,
						0,
						0));
      SW4_MARK_END(" prefetch_to_device");
    } //else std::cerr<<"Zero size prefetch \n";
  } else std::cerr<<"NO prefetch due to unknown address\n";
}

void CheckError(cudaError_t const err, const char* file, char const* const fun, const int line)
{
    if (err)
    {
      std::cerr<<"CUDA Error Code["<<err<<"]: "<<cudaGetErrorString(err)<<" "<<file<<" "<<fun<<" Line number:  "<<line<<"\n";
    }
}
#endif

// AUTOPEEL CODE
#if defined(SW4_TRACK_MEMORY_ALLOCATIONS)
// AUTOPEEL uses getsize which only works when SW4_TRACK_MEMORY_ALLOCATIONS is defined
std::string line(int n,int C){ 
  std::ostringstream buf;
  buf<<"int arg"<<C<<" = "<<n<<";\n";
  return buf.str();
}
std::string line(int *n,int C){ 
  std::ostringstream buf;
  buf<<"int arg"<<C<<"[6]={";
  for(int i=0;i<5;i++)buf<<n[i]<<",";
  buf<<n[5]<<"};\n";
  return buf.str();
}
std::string line(double n,int C){ 
  std::ostringstream buf;
  buf<<"double arg"<<C<<" = "<<n<<";\n";
  return buf.str();
}
std::string line(double *n,int C){ 
  std::ostringstream buf;
  buf<<"double *arg"<<C<<";\n cudaMallocManaged((void*)&arg"<<C<<","<<getsize((void**)n)<<");\n";
  return buf.str();
}
std::string line(char n,int C){ 
  std::ostringstream buf;
  buf<<"char arg"<<C<<" = \""<<n<<"\";\n";
  return buf.str();
}



Apc::Apc(char *s){
  counter=0;
  ofile.open(s,std::ios::out);
  ofile<<"int main(int argc, char *argv[]){\n";
}
Apc::~Apc(){
  ofile<<"FUNCTION(";
  for(int i=0;i<counter;i++)ofile<<"arg"<<i<<",";
  ofile<<"arg"<<counter<<");\n";
  ofile<<"}\n";
  ofile.close();
}
#endif

// END AUTOPEEL CODE
