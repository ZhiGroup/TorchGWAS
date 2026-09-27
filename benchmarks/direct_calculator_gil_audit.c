/* Independent CPU-service probe. No CUDA or association kernels.
 * LD_AUDIT observes extension calls even when CPython is in the executable.
 * Control entry points are redirected into this audit namespace, keeping
 * counters on the same thread-local instance as the intercepted calls.
 */
#define _GNU_SOURCE
#include <link.h>
#include <stdint.h>
#include <string.h>
#include <time.h>
#include <stdlib.h>
#include <stdatomic.h>

typedef void *(*save_function)(void);
typedef void (*restore_function)(void *);
typedef int (*ensure_function)(void);
typedef void (*release_function)(int);
static _Atomic(save_function) real_save;
static _Atomic(restore_function) real_restore;
static _Atomic(ensure_function) real_ensure;
static _Atomic(release_function) real_release;
static __thread int active, detached, errors;
static __thread uint64_t started, saved, released, intervals;

static uint64_t cpu_ns(void) {
    struct timespec t;
    if (clock_gettime(CLOCK_THREAD_CPUTIME_ID, &t)) abort();
    return (uint64_t)t.tv_sec * 1000000000ULL + (uint64_t)t.tv_nsec;
}

static void start_detached_meter(void) {
    if (active) {
        if (detached) ++errors;
        detached = 1;
        saved = cpu_ns();
    }
}

static void stop_detached_meter(void) {
    if (active) {
        if (!detached) ++errors;
        else { released += cpu_ns() - saved; ++intervals; }
        detached = 0;
    }
}

static void *observed_save(void) {
    void *state = atomic_load(&real_save)();
    start_detached_meter();
    return state;
}

static void observed_restore(void *state) {
    stop_detached_meter();
    /* Reacquisition CPU is conservatively left in the held/unknown remainder.
     * Wall-clock waiting is never charged as CPU service. */
    atomic_load(&real_restore)(state);
}

static int observed_ensure(void) {
    if (active && detached) {
        released += cpu_ns() - saved; ++intervals; detached = 0;
    }
    return atomic_load(&real_ensure)();
}

static void observed_release(int state) {
    atomic_load(&real_release)(state);
    /* PyGILState_UNLOCKED is 1; LOCKED is 0. */
    if (active && state == 1) {
        if (detached) ++errors;
        detached = 1; saved = cpu_ns();
    }
}

void gil_probe_begin_mode(int record) {
    active = record; detached = errors = 0; released = intervals = 0;
    started = cpu_ns();
}

void gil_probe_begin(void) { gil_probe_begin_mode(1); }

/* Meter-only control: same two CPU-clock/counter sections, with no native
 * work, GIL transitions or Python loop between them. Used only to estimate
 * measurement overhead; this deliberately does not label real GIL state. */
void gil_probe_meter_batch(unsigned count) {
    for (unsigned i = 0; i < count; ++i) {
        __asm__ volatile ("" ::: "memory");
        start_detached_meter();
        stop_detached_meter();
    }
}

void gil_probe_end(double *out) {
    uint64_t ended = cpu_ns();
    active = 0;
    out[0] = (ended - started) * 1e-9;
    out[1] = released * 1e-9;
    out[2] = (double)intervals;
    out[3] = (double)(errors + detached);
}

/* Fixed generic native work for the audit's positive control. */
double gil_probe_native_work(unsigned iterations) {
    volatile double value = 1.;
    for (unsigned i = 0; i < iterations; ++i) value = value * 1.00000001 + .00000001;
    return value;
}

/* Exercise the same auditable PLT calls as libtorch_python. _ctypes itself
 * uses GLOB_DAT relocations here, which GNU auditing does not intercept. */
extern void *PyEval_SaveThread(void);
extern void PyEval_RestoreThread(void *);
extern int PyGILState_Ensure(void);
extern void PyGILState_Release(int);
double gil_probe_detached_work(unsigned iterations) {
    void *state = PyEval_SaveThread();
    double result = gil_probe_native_work(iterations);
    PyEval_RestoreThread(state);
    return result;
}

double gil_probe_nested_work(unsigned iterations) {
    void *outer = PyEval_SaveThread();
    double result = gil_probe_native_work(iterations);
    int state = PyGILState_Ensure();
    result += gil_probe_detached_work(iterations);
    PyGILState_Release(state);
    result += gil_probe_native_work(iterations);
    PyEval_RestoreThread(outer);
    return result;
}

unsigned int la_version(unsigned int version) {
    return version < LAV_CURRENT ? 0 : LAV_CURRENT;
}

unsigned int la_objopen(struct link_map *map, Lmid_t lmid, uintptr_t *cookie) {
    (void)map; (void)lmid; (void)cookie;
    return LA_FLG_BINDTO | LA_FLG_BINDFROM;
}

uintptr_t la_symbind64(Elf64_Sym *symbol, unsigned int index, uintptr_t *ref,
                     uintptr_t *def, unsigned int *flags, const char *name) {
    (void)index; (void)ref; (void)def; (void)flags;
    if (!strcmp(name, "PyEval_SaveThread")) {
        atomic_store(&real_save, (save_function)symbol->st_value);
        return (uintptr_t)observed_save;
    }
    if (!strcmp(name, "PyEval_RestoreThread")) {
        atomic_store(&real_restore, (restore_function)symbol->st_value);
        return (uintptr_t)observed_restore;
    }
    if (!strcmp(name, "PyGILState_Ensure")) {
        atomic_store(&real_ensure, (ensure_function)symbol->st_value);
        return (uintptr_t)observed_ensure;
    }
    if (!strcmp(name, "PyGILState_Release")) {
        atomic_store(&real_release, (release_function)symbol->st_value);
        return (uintptr_t)observed_release;
    }
    if (!strcmp(name, "gil_probe_begin_mode")) return (uintptr_t)gil_probe_begin_mode;
    if (!strcmp(name, "gil_probe_meter_batch")) return (uintptr_t)gil_probe_meter_batch;
    if (!strcmp(name, "gil_probe_begin")) return (uintptr_t)gil_probe_begin;
    if (!strcmp(name, "gil_probe_end")) return (uintptr_t)gil_probe_end;
    return symbol->st_value;
}
