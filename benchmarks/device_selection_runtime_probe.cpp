// Host-only measurement shim for independent CUDA API primitives. No kernels.
#include <cuda_runtime_api.h>
#include <dlfcn.h>
#include <time.h>
#include <cstdint>
#include <cstdio>
#include <cstdlib>

struct Event {
    uint64_t kind, bytes, direction, begin_wall, end_wall, begin_cpu, end_cpu, result;
};
struct Record {
    uint64_t begin_wall, end_wall, begin_cpu, end_cpu, return_wall, return_cpu, count, overflow;
    Event events[64];
};
static thread_local Record record{};
static thread_local bool enabled=false;
static uint64_t stamp(clockid_t clock) {
    timespec value; if(clock_gettime(clock,&value)) std::abort();
    return uint64_t(value.tv_sec)*1000000000ULL+uint64_t(value.tv_nsec);
}
template<class T> T next(const char* name) {
    void* value=dlsym(RTLD_NEXT,name);
    if(!value) { std::fprintf(stderr,"Unresolved CUDA runtime probe symbol %s: %s\n",name,dlerror()); std::abort(); }
    return reinterpret_cast<T>(value);
}
static void begin_record(bool capture) {
    record=Record{};
    record.begin_wall=stamp(CLOCK_MONOTONIC_RAW);
    record.begin_cpu=stamp(CLOCK_THREAD_CPUTIME_ID);
    enabled=capture;
}
extern "C" void selection_probe_begin() { begin_record(true); }
extern "C" void selection_probe_begin_uninstrumented() { begin_record(false); }
extern "C" void selection_probe_return() {
    record.return_cpu=stamp(CLOCK_THREAD_CPUTIME_ID);
    record.return_wall=stamp(CLOCK_MONOTONIC_RAW);
}
extern "C" void selection_probe_end(Record* output) {
    enabled=false;
    record.end_cpu=stamp(CLOCK_THREAD_CPUTIME_ID);
    record.end_wall=stamp(CLOCK_MONOTONIC_RAW);
    *output=record;
}
extern "C" uint64_t selection_probe_record_bytes() { return sizeof(Record); }
template<class F> cudaError_t observe(uint64_t kind,uint64_t bytes,uint64_t direction,F function) {
    if(!enabled) return function();
    const uint64_t index=record.count++;
    if(index>=64) { ++record.overflow; return function(); }
    Event& event=record.events[index];
    event.kind=kind;event.bytes=bytes;event.direction=direction;
    event.begin_wall=stamp(CLOCK_MONOTONIC_RAW);event.begin_cpu=stamp(CLOCK_THREAD_CPUTIME_ID);
    auto result=function();
    event.end_cpu=stamp(CLOCK_THREAD_CPUTIME_ID);event.end_wall=stamp(CLOCK_MONOTONIC_RAW);
    event.result=uint64_t(result);
    return result;
}
extern "C" cudaError_t cudaMemcpyAsync(void* to,const void* from,size_t bytes,cudaMemcpyKind kind,cudaStream_t stream) {
    static auto real=next<decltype(&cudaMemcpyAsync)>("cudaMemcpyAsync");
    return observe(1,bytes,uint64_t(kind),[&] {return real(to,from,bytes,kind,stream);});
}
extern "C" cudaError_t cudaStreamSynchronize(cudaStream_t stream) {
    static auto real=next<decltype(&cudaStreamSynchronize)>("cudaStreamSynchronize");
    return observe(2,0,0,[&] {return real(stream);});
}
extern "C" cudaError_t cudaDeviceSynchronize() {
    static auto real=next<decltype(&cudaDeviceSynchronize)>("cudaDeviceSynchronize");
    return observe(3,0,0,[&] {return real();});
}
extern "C" cudaError_t cudaEventSynchronize(cudaEvent_t event) {
    static auto real=next<decltype(&cudaEventSynchronize)>("cudaEventSynchronize");
    return observe(4,0,0,[&] {return real(event);});
}
extern "C" cudaError_t cudaMemcpy(void* to,const void* from,size_t bytes,cudaMemcpyKind kind) {
    static auto real=next<decltype(&cudaMemcpy)>("cudaMemcpy");
    return observe(5,bytes,uint64_t(kind),[&] {return real(to,from,bytes,kind);});
}
