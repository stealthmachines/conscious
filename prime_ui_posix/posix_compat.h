/* posix_compat.h — Windows → POSIX/Linux compatibility shim for prime_ui
 *
 * Include this INSTEAD of <windows.h> and <bcrypt.h> on non-Windows targets.
 * Provides: Windows types, console API stubs, timing API, BCryptGenRandom,
 * file/path helpers, process priority + CPU affinity wrappers.
 *
 * Targets: Linux x86-64 with GCC ≥ 7 or Clang ≥ 6.
 * Build flags: -O2 -mavx2 -mfma -msse4.1 -maes -mpclmul -mrdseed -mrdrnd
 *              -D_GNU_SOURCE -lm
 */
#ifndef POSIX_COMPAT_H
#define POSIX_COMPAT_H

/* ── POSIX feature test ───────────────────────────────────────────────────── */
#ifndef _GNU_SOURCE
#  define _GNU_SOURCE
#endif

/* ── Standard headers ─────────────────────────────────────────────────────── */
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <time.h>
#include <termios.h>
#include <sys/stat.h>
#include <sys/time.h>
#include <sys/resource.h>
#include <sched.h>
#include <limits.h>
#include <errno.h>

/* getrandom() — Linux 3.17+ */
#if defined(__linux__)
#  include <sys/random.h>
#endif

/* x86 intrinsics — GCC needs x86intrin.h; Clang accepts intrin.h on Linux
 * but x86intrin.h is canonical and works on both.  Included before the rest
 * of prime_ui_posix.c uses <immintrin.h> etc. */
#if defined(__GNUC__) || defined(__clang__)
#  include <x86intrin.h>   /* __rdtsc, _rdseed64_step, _rdrand64_step */
#  include <immintrin.h>   /* AVX2/FMA/SSE4.1 */
#endif

/* ── Basic Windows types ──────────────────────────────────────────────────── */
typedef uint32_t   DWORD;
typedef uint64_t   DWORD_PTR;
typedef uintptr_t  UINT_PTR;
typedef int        BOOL;
typedef uint8_t    BYTE;
typedef uint16_t   WORD;
typedef int        HANDLE;   /* fd or dummy; unused on POSIX */
typedef uint16_t   WCHAR;    /* UTF-16LE; only ASCII range is used */

#define INVALID_HANDLE_VALUE  (-1)
#ifndef TRUE
#  define TRUE   1
#  define FALSE  0
#endif

/* ── MAX_PATH ─────────────────────────────────────────────────────────────── */
#ifndef PATH_MAX
#  define PATH_MAX  4096
#endif
#ifndef MAX_PATH
#  define MAX_PATH  PATH_MAX
#endif

/* ── LARGE_INTEGER ────────────────────────────────────────────────────────── */
typedef union {
    struct { uint32_t LowPart; int32_t HighPart; };
    int64_t QuadPart;
} LARGE_INTEGER;

/* ── FILETIME ─────────────────────────────────────────────────────────────── */
typedef struct {
    uint32_t dwLowDateTime;
    uint32_t dwHighDateTime;
} FILETIME;

/* ── Console mode constants (referenced but no-op on POSIX) ──────────────── */
/* Note: g_hout and g_hin are declared in prime_ui_posix.c via console_init() */
#define ENABLE_PROCESSED_INPUT             0x0001
#define ENABLE_LINE_INPUT                  0x0002
#define ENABLE_ECHO_INPUT                  0x0004
#define ENABLE_VIRTUAL_TERMINAL_PROCESSING 0x0004
#define ENABLE_PROCESSED_OUTPUT            0x0001
#define STD_OUTPUT_HANDLE  ((DWORD)(-11))
#define STD_INPUT_HANDLE   ((DWORD)(-10))
#define STD_ERROR_HANDLE   ((DWORD)(-12))
#define CP_UTF8            65001

/* GetConsoleMode / SetConsoleMode — no-op; Linux terminals speak VT by default */
static inline int GetConsoleMode(HANDLE h, DWORD *mode)    { (void)h; *mode = 0; return 1; }
static inline int SetConsoleMode(HANDLE h, DWORD mode)     { (void)h; (void)mode; return 1; }
static inline int SetConsoleOutputCP(DWORD cp)             { (void)cp; return 1; }
static inline int SetConsoleCP(DWORD cp)                   { (void)cp; return 1; }
static inline HANDLE GetStdHandle(DWORD n) {
    if (n == STD_INPUT_HANDLE)  return 0;
    if (n == STD_OUTPUT_HANDLE) return 1;
    return 2;
}

/* ── ReadConsoleW — reads a line from stdin, stores as WCHAR (ASCII subset) ─ */
/* The Windows ReadConsoleW includes the trailing \r\n in the buffer.          */
/* This implementation does the same so existing trim loops are unaffected.    */
static inline int ReadConsoleW(HANDLE h, WCHAR *buf, DWORD maxchars,
                               DWORD *nr, void *reserved)
{
    (void)h; (void)reserved;
    char nbuf[1024] = {0};
    if (!fgets(nbuf, (int)((maxchars < 1022u) ? maxchars : 1022u), stdin)) {
        *nr = 0;
        return 0;
    }
    DWORD n = 0;
    while (nbuf[n] && n < maxchars - 1) {
        buf[n] = (WCHAR)(unsigned char)nbuf[n];
        n++;
    }
    buf[n] = 0;
    *nr = n;
    return 1;
}

/* ── WideCharToMultiByte — trivial ASCII passthrough ─────────────────────── */
static inline int WideCharToMultiByte(DWORD cp, DWORD flags,
                                      const WCHAR *ws, int wlen,
                                      char *mb, int mbsz,
                                      const char *def, int *used)
{
    (void)cp; (void)flags; (void)def; (void)used;
    if (wlen < 0) { wlen = 0; while (ws[wlen]) wlen++; }
    int out = 0;
    while (out < wlen && out < mbsz - 1) {
        mb[out] = (char)(ws[out] & 0x7F);   /* ASCII passthrough */
        out++;
    }
    if (mbsz > 0) mb[out] = '\0';
    return out + 1;
}

/* ── QueryPerformanceFrequency / Counter — map to CLOCK_MONOTONIC ───────── */
static inline int QueryPerformanceFrequency(LARGE_INTEGER *f) {
    f->QuadPart = 1000000000LL;   /* nanosecond resolution */
    return 1;
}
static inline int QueryPerformanceCounter(LARGE_INTEGER *c) {
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    c->QuadPart = (int64_t)ts.tv_sec * 1000000000LL + (int64_t)ts.tv_nsec;
    return 1;
}

/* ── GetSystemTimeAsFileTime — FILETIME = 100-ns ticks since 1601-01-01 ──── */
static inline void GetSystemTimeAsFileTime(FILETIME *ft) {
    struct timespec ts;
    clock_gettime(CLOCK_REALTIME, &ts);
    /* Unix epoch to Windows FILETIME epoch: +11644473600 seconds */
    uint64_t v = ((uint64_t)ts.tv_sec + 11644473600ULL) * 10000000ULL
                 + (uint64_t)((unsigned long)ts.tv_nsec / 100u);
    ft->dwLowDateTime  = (uint32_t)(v & 0xFFFFFFFFu);
    ft->dwHighDateTime = (uint32_t)(v >> 32);
}

/* ── BCryptGenRandom → getrandom() (Linux 3.17+) / /dev/urandom fallback ── */
#define BCRYPT_USE_SYSTEM_PREFERRED_RNG 0x00000002u
static inline long BCryptGenRandom(void *alg, uint8_t *buf, size_t len, DWORD flags)
{
    (void)alg; (void)flags;
#if defined(__linux__) && defined(SYS_getrandom)
    {
        ssize_t got = getrandom(buf, len, 0);
        if (got == (ssize_t)len) return 0;   /* 0 = STATUS_SUCCESS */
    }
#endif
    /* fallback: /dev/urandom */
    FILE *f = fopen("/dev/urandom", "rb");
    if (!f) return -1L;
    size_t got = fread(buf, 1, len, f);
    fclose(f);
    return ((got == len) ? 0L : -1L);
}

/* ── GetCurrentDirectoryA → getcwd ──────────────────────────────────────── */
static inline int GetCurrentDirectoryA(DWORD sz, char *buf) {
    char *r = getcwd(buf, (size_t)sz);
    return r ? (int)strlen(buf) : 0;
}

/* ── GetFileAttributesA → stat ───────────────────────────────────────────── */
#define INVALID_FILE_ATTRIBUTES ((DWORD)-1u)
static inline DWORD GetFileAttributesA(const char *path) {
    struct stat st;
    return (stat(path, &st) == 0) ? 0x80u : INVALID_FILE_ATTRIBUTES;
}

/* ── Sleep(ms) → usleep ──────────────────────────────────────────────────── */
static inline void Sleep(DWORD ms) {
    usleep((useconds_t)ms * 1000u);
}

/* ── Process priority classes ────────────────────────────────────────────── */
#define REALTIME_PRIORITY_CLASS       0x00000100u
#define HIGH_PRIORITY_CLASS           0x00000080u
#define ABOVE_NORMAL_PRIORITY_CLASS   0x00008000u
#define NORMAL_PRIORITY_CLASS         0x00000020u
#define BELOW_NORMAL_PRIORITY_CLASS   0x00004000u
#define IDLE_PRIORITY_CLASS           0x00000040u

static inline HANDLE GetCurrentProcess(void) { return (HANDLE)0; }

static inline int SetPriorityClass(HANDLE hp, DWORD wclass) {
    (void)hp;
    int nice_val = 0;
    if      (wclass == REALTIME_PRIORITY_CLASS)      nice_val = -20;
    else if (wclass == HIGH_PRIORITY_CLASS)           nice_val = -10;
    else if (wclass == ABOVE_NORMAL_PRIORITY_CLASS)   nice_val = -5;
    else if (wclass == BELOW_NORMAL_PRIORITY_CLASS)   nice_val =  5;
    else if (wclass == IDLE_PRIORITY_CLASS)           nice_val =  19;
    /* setpriority may require CAP_SYS_NICE for negative values; ignore EPERM */
    return setpriority(PRIO_PROCESS, 0, nice_val) == 0 || errno == EPERM;
}

static inline int GetProcessAffinityMask(HANDLE hp,
                                         DWORD_PTR *proc_mask,
                                         DWORD_PTR *sys_mask)
{
    (void)hp;
    cpu_set_t cs;
    CPU_ZERO(&cs);
    if (sched_getaffinity(0, sizeof(cs), &cs) == 0) {
        DWORD_PTR mask = 0;
        for (int i = 0; i < 64 && i < CPU_SETSIZE; i++)
            if (CPU_ISSET(i, &cs)) mask |= (DWORD_PTR)1u << i;
        *proc_mask = mask;
        *sys_mask  = mask;
        return 1;
    }
    *proc_mask = 1;
    *sys_mask  = 1;
    return 0;
}

static inline int SetProcessAffinityMask(HANDLE hp, DWORD_PTR mask) {
    (void)hp;
    cpu_set_t cs;
    CPU_ZERO(&cs);
    for (int i = 0; i < 64 && i < CPU_SETSIZE; i++)
        if (mask & ((DWORD_PTR)1u << i)) CPU_SET(i, &cs);
    return sched_setaffinity(0, sizeof(cs), &cs) == 0;
}

/* ── _popen / _pclose → popen / pclose ──────────────────────────────────── */
#define _popen   popen
#define _pclose  pclose

/* ── _umul128 — Windows/MSVC intrinsic: use __uint128_t on GCC/Clang ─────── */
/* This macro covers the #else branch in mulmod64 when compiled with Clang    */
/* on Linux where __clang__ is defined but _WIN32 is not.                     */
#if !defined(_WIN32)
static inline uint64_t _umul128_posix(uint64_t a, uint64_t b, uint64_t *hi_out) {
    __uint128_t r = (__uint128_t)a * b;
    *hi_out = (uint64_t)(r >> 64);
    return (uint64_t)(r & (uint64_t)-1);
}
#  define _umul128(a, b, hi)  _umul128_posix((a), (b), (hi))
#endif

#endif /* POSIX_COMPAT_H */
