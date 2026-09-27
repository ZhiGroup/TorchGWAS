/* Experiment-only NumPy data allocator, one owning thread per process.
 * No production import. Resident mode changes only allocation, not NumPy
 * numerical kernels. A bounded pre-touched arena is reset only after every
 * array allocated from it has been released. Outputs remain valid after the
 * default handler is restored. Realloc is deliberately unsupported in arena
 * mode and must be absent from accepted observations.
 */
#define _GNU_SOURCE
#define NPY_TARGET_VERSION NPY_1_22_API_VERSION
#define NPY_NO_DEPRECATED_API NPY_1_7_API_VERSION
#include <Python.h>
#include <numpy/arrayobject.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <stdatomic.h>
#include <sys/resource.h>
#include <unistd.h>

typedef struct { size_t size; uint64_t magic; } Header;
static const uint64_t MAGIC = UINT64_C(0x5ac173eae25107bb);
static PyDataMem_Handler parent;
static PyObject *saved_policy, *policy;
static unsigned char *arena;
static size_t capacity, offset, live, malloc_calls, calloc_calls, free_calls;
static size_t allocate_bytes, free_bytes, realloc_calls, failures, invalid_free;
static uint64_t allocate_ns, free_ns;
static uint64_t prefault_ns;
static size_t prefault_minor, prefault_major, page_bytes;
static unsigned long owner;
static int active, mode;
static _Atomic int foreign_thread;

static uint64_t now(void) {
    struct timespec t;
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID,&t)) abort();
    return (uint64_t)t.tv_sec*UINT64_C(1000000000)+(uint64_t)t.tv_nsec;
}
static int owns_thread(void) {
    if (PyThread_get_thread_ident()==owner) return 1;
    atomic_store(&foreign_thread,1); return 0;
}
static void *resident_alloc(size_t size) {
    size_t bytes=size ? size : 1;
    if (bytes>SIZE_MAX-127) return NULL;
    size_t extent=64+((bytes+63)&~(size_t)63);
    if (extent>capacity-offset) return NULL;
    Header *header=(Header *)(arena+offset);
    header->size=size; header->magic=MAGIC;
    void *out=arena+offset+64; offset+=extent;
    return out;
}
static void prefault(void *ptr,size_t size) {
    if (!ptr || mode!=3 || size<page_bytes) return;
    struct rusage before,after;
    if (getrusage(RUSAGE_THREAD,&before)) abort();
    uint64_t start=now();
    /* Only uninitialized malloc or already-zero calloc storage is touched.
     * Include an unaligned final page without accessing past the allocation. */
    volatile unsigned char *bytes=ptr;
    bytes[0]=0;
    size_t first=page_bytes-((uintptr_t)ptr%page_bytes);
    for(size_t index=first;index<size;index+=page_bytes) bytes[index]=0;
    prefault_ns+=now()-start;
    if (getrusage(RUSAGE_THREAD,&after)) abort();
    prefault_minor+=(size_t)(after.ru_minflt-before.ru_minflt);
    prefault_major+=(size_t)(after.ru_majflt-before.ru_majflt);
}
static void *allocate(void *ctx,size_t size) {
    (void)ctx;
    if (!owns_thread()) return NULL;
    uint64_t start=now();
    void *out=mode==2 ? resident_alloc(size) : parent.allocator.malloc(parent.allocator.ctx,size);
    prefault(out,size);
    allocate_ns+=now()-start; ++malloc_calls; allocate_bytes+=size;
    if (out) ++live; else ++failures;
    return out;
}
static void *allocate_zero(void *ctx,size_t count,size_t size) {
    (void)ctx;
    if (!owns_thread()) return NULL;
    if (size && count>SIZE_MAX/size) { ++failures; return NULL; }
    size_t bytes=count*size;uint64_t start=now();void *out;
    if (mode==2) { out=resident_alloc(bytes); if(out) memset(out,0,bytes); }
    else out=parent.allocator.calloc(parent.allocator.ctx,count,size);
    prefault(out,bytes);
    allocate_ns+=now()-start; ++calloc_calls; allocate_bytes+=bytes;
    if (out) ++live; else ++failures;
    return out;
}
static void release(void *ctx,void *ptr,size_t size) {
    (void)ctx;
    if (!ptr || !owns_thread()) return;
    uint64_t start=now();
    if (mode==2) {
        uintptr_t address=(uintptr_t)ptr,base=(uintptr_t)arena;
        if (address<base+64 || address>=base+offset || (address-base)%64) {++invalid_free;return;}
        Header *header=(Header *)((unsigned char *)ptr-64);
        if (header->magic!=MAGIC) {++invalid_free;return;}
        header->magic=0;
    } else parent.allocator.free(parent.allocator.ctx,ptr,size);
    free_ns+=now()-start; ++free_calls; free_bytes+=size;
    if (live) --live; else ++invalid_free;
}
static void *resize(void *ctx,void *ptr,size_t size) {
    (void)ctx;
    if (!owns_thread()) return NULL;
    ++realloc_calls;
    /* Reject in both modes: a probe must use matching allocation semantics. */
    (void)ptr;(void)size; ++failures; return NULL;
}
static PyDataMem_Handler handler={"selector_allocator_control",1,{NULL,allocate,allocate_zero,resize,release}};

static PyObject *configure(PyObject *self,PyObject *arg) {
    (void)self;size_t size=PyLong_AsSize_t(arg);
    if (PyErr_Occurred()) return NULL;
    if (active || live || arena || size<4096 || size>((size_t)512<<20)) {
        PyErr_SetString(PyExc_ValueError,"One bounded arena, configured before use, required");return NULL;
    }
    if (posix_memalign((void **)&arena,64,size)) return PyErr_NoMemory();
    capacity=size;memset(arena,0,size);Py_RETURN_NONE;
}
static PyObject *begin(PyObject *self,PyObject *arg) {
    (void)self;long requested=PyLong_AsLong(arg);
    if (PyErr_Occurred()) return NULL;
    if (active || live || requested<0 || requested>3 || (requested==2 && !arena)) {
        PyErr_SetString(PyExc_ValueError,"Restore handler and release all live arrays before reset");return NULL;
    }
    PyObject *current=PyDataMem_GetHandler();if(!current)return NULL;
    PyDataMem_Handler *p=PyCapsule_GetPointer(current,"mem_handler");
    if (!p) {Py_DECREF(current);return NULL;}
    if (strcmp(p->name,"default_allocator")) {
        Py_DECREF(current);PyErr_SetString(PyExc_ValueError,"Default NumPy allocator required");return NULL;
    }
    parent=*p;mode=(int)requested;owner=PyThread_get_thread_ident();offset=0;
    malloc_calls=calloc_calls=free_calls=allocate_bytes=free_bytes=realloc_calls=failures=invalid_free=0;
    allocate_ns=free_ns=0;atomic_store(&foreign_thread,0);
    prefault_ns=0;prefault_minor=prefault_major=0;
    if (mode) {
        PyObject *old=PyDataMem_SetHandler(policy);
        if(!old) {Py_DECREF(current);return NULL;}
        Py_DECREF(old);
    }
    saved_policy=current;active=1;Py_RETURN_NONE;
}
static PyObject *restore(PyObject *self,PyObject *arg) {
    (void)self;(void)arg;
    if (!active || !owns_thread()) {PyErr_SetString(PyExc_ValueError,"No owning active scope");return NULL;}
    PyObject *old=PyDataMem_SetHandler(saved_policy);if(!old)return NULL;
    Py_DECREF(old);Py_CLEAR(saved_policy);active=0;Py_RETURN_NONE;
}
static PyObject *snapshot(PyObject *self,PyObject *arg) {
    (void)self;(void)arg;PyObject *result=PyDict_New();if(!result)return NULL;
    #define PUT(name,value) do { PyObject *v=PyLong_FromUnsignedLongLong(value); \
        if(!v || PyDict_SetItemString(result,name,v)) {Py_XDECREF(v);Py_DECREF(result);return NULL;} Py_DECREF(v); } while(0)
    PUT("malloc_calls",malloc_calls);PUT("calloc_calls",calloc_calls);PUT("free_calls",free_calls);
    PUT("allocate_bytes",allocate_bytes);PUT("free_bytes",free_bytes);PUT("realloc_calls",realloc_calls);
    PUT("failures",failures);PUT("invalid_free",invalid_free);PUT("foreign_thread",atomic_load(&foreign_thread));
    PUT("live",live);PUT("arena_used_bytes",offset);PUT("allocate_cpu_ns",allocate_ns);PUT("free_cpu_ns",free_ns);
    PUT("prefault_cpu_ns",prefault_ns);PUT("prefault_minor_faults",prefault_minor);PUT("prefault_major_faults",prefault_major);
    #undef PUT
    return result;
}
static PyMethodDef methods[]={
    {"configure",configure,METH_O,"Allocate and pre-touch a bounded arena once."},
    {"begin",begin,METH_O,"Begin default (0), passthrough (1), resident (2), or default-prefault (3) scope."},
    {"restore",restore,METH_NOARGS,"Restore default; live outputs keep their recorded handler."},
    {"snapshot",snapshot,METH_NOARGS,"Allocator counters for the last scope, including later frees."},
    {NULL,NULL,0,NULL}};
static struct PyModuleDef module={PyModuleDef_HEAD_INIT,"selector_allocator_control",NULL,-1,methods};
PyMODINIT_FUNC PyInit_selector_allocator_control(void) {
    import_array();
    long page=sysconf(_SC_PAGESIZE);if(page<=0)return NULL;page_bytes=(size_t)page;
    policy=PyCapsule_New(&handler,"mem_handler",NULL);if(!policy)return NULL;
    return PyModule_Create(&module);
}
