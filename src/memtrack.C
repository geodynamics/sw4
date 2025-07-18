#include <stddef.h>  // size_t
#include <stdio.h>  // fopen, fscanf
#include <errno.h>  // errno, not shown but the "file" functions set "errno" if error
#include <unistd.h>  // sysconf
#include <hip/hip_runtime.h>  // hipMemGetInfo
#include <iostream>
#include <iomanip>


#define HIP_CHECK(expression)                  \
{                                              \
    const hipError_t status = expression;      \
    if(status != hipSuccess){                  \
        std::cerr << "HIP error "              \
                  << status << ": "            \
                  << hipGetErrorString(status) \
                  << " at " << __FILE__ << ":" \
                  << __LINE__ << std::endl;    \
    }                                          \
}
// Get host-allocated memory by the calling process (units in bytes)
size_t get_process_allocated_memory() {
    // These values are measured in pages (default page size)
    // https://www.man7.org/linux/man-pages/man5/proc_pid_statm.5.html
    FILE * const fd = fopen("/proc/self/statm", "rb");
    unsigned long int used_mem = 0;
    // heap data + stack
    fscanf(fd, "%*lu %*lu %*lu %*lu %*lu %lu %*lu", &used_mem);
    return (size_t)used_mem * (size_t)sysconf(_SC_PAGESIZE);
}


// Get hipMalloc'ed memory of the logical HIP device currently set for the calling process (units in bytes)
size_t get_hip_allocated_memory() {
    size_t total_mem = 0, free_mem = 0;
    HIP_CHECK(hipMemGetInfo(&free_mem, &total_mem));
    return total_mem - free_mem;
}

void print_mem_usage(){
auto host = get_process_allocated_memory();
auto device = get_hip_allocated_memory();
auto total = host + device;
std::cout<<std::fixed<<std::setprecision(2)<<"Memory used host: "<<host/1024/1024/1024.0<<" GB, device "<<device/1024/1024/1024.0<<"GB, total "<<total/1024/1024/1024.0<<" GB\n"<<std::defaultfloat;
}
