/* prime_ui.c — Windows TUI for the Quantum-Prime Library v1.0
 *
 *  [1] Prime Pipeline     sieve -> phi-filter -> Dn-rank (Mersenne candidates)
 *  [2] Number Analyzer    Miller-Rabin, phi-lattice, Dn score, psi-score
 *  [3] Mersenne Explorer  all M1..M51 + next candidate predictions
 *  [4] Zeta Zeros         zeta(1/2+it) hardcoded table + Gram approximation
 *  [5] Benchmark          time all 13 prime library functions
 *  [6] Alpine Install      boot lattice + GPU resonance hook
 *  [7] Lattice Shell       interactive Slot4096 REPL
 *  [8] Alpine OS Shell     spawn lattice-powered Alpine Linux (Docker/WSL)
 *  [Q] Quit
 *
 * Build:  build_prime_ui.bat
 *   -or-  clang -O2 -D_CRT_SECURE_NO_WARNINGS -D_USE_MATH_DEFINES ^
 *               prime_ui.c -o prime_ui.exe
 *
 * Self-contained.  All math inlined from bench_prime_funcs.c,
 * prime_pipeline.c, phi_mersenne_predictor.c.
 * No external dependencies beyond MSVCRT + kernel32.
 */

#ifndef _CRT_SECURE_NO_WARNINGS
#  define _CRT_SECURE_NO_WARNINGS
#endif
#ifndef _USE_MATH_DEFINES
#  define _USE_MATH_DEFINES
#endif
#define WIN32_LEAN_AND_MEAN
#include <windows.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdint.h>
#include <ctype.h>

/* ══════════════════════════ A. Constants ════════════════════════════════════ */

#define PHI      1.6180339887498948482
#define LN_PHI   0.4812118250596034748
#define LOG10PHI 0.2090150076824960
#define SQRT5    2.2360679774997896
#define PI       3.14159265358979323846
#define INV_E    0.36787944117144232159
#define TWO_PI   6.28318530717958647693
#define M_LN2_V  0.6931471805599453094
#define LOG10_2  0.3010299957316877

static const int PRIMES50[50] = {
     2,  3,  5,  7, 11, 13, 17, 19, 23, 29,
    31, 37, 41, 43, 47, 53, 59, 61, 67, 71,
    73, 79, 83, 89, 97,101,103,107,109,113,
   127,131,137,139,149,151,157,163,167,173,
   179,181,191,193,197,199,211,223,227,229
};

static const uint64_t MERSENNE_EXP[51] = {
    2ULL,3ULL,5ULL,7ULL,13ULL,17ULL,19ULL,31ULL,61ULL,89ULL,
    107ULL,127ULL,521ULL,607ULL,1279ULL,2203ULL,2281ULL,3217ULL,
    4253ULL,4423ULL,9689ULL,9941ULL,11213ULL,19937ULL,21701ULL,
    23209ULL,44497ULL,86243ULL,110503ULL,132049ULL,216091ULL,
    756839ULL,859433ULL,1257787ULL,1398269ULL,2976221ULL,3021377ULL,
    6972593ULL,13466917ULL,20996011ULL,24036583ULL,25964951ULL,
    30402457ULL,32582657ULL,37156667ULL,42643801ULL,43112609ULL,
    57885161ULL,74207281ULL,77232917ULL,136279841ULL
};
#define N_MERSENNE 51

static const double ZETA_ZEROS_80[80] = {
    14.134725141734693,  21.022039638771555,  25.010857580145688,
    30.424876125859513,  32.935061587739189,  37.586178158825671,
    40.918719012147495,  43.327073280914999,  48.005150881167159,
    49.773832477672302,  52.970321477714460,  56.446247697063246,
    59.347044002602352,  60.831778524609809,  65.112544048081560,
    67.079810529494173,  69.546401711173979,  72.067157674481907,
    75.704690699083933,  77.144840068874805,  79.337375020249367,
    82.910380854160462,  84.735492981074628,  87.425274613125229,
    88.809111207634465,  92.491899270593585,  94.651344040519681,
    95.870634228245332,  98.831194218193159, 101.317851006956152,
   103.725538040478419, 105.446623052947866, 107.168611184276793,
   111.029535543169970, 111.874659177229233, 114.320220915452460,
   116.226680321519019, 118.790782866217474, 121.370125002980428,
   122.946829294236573, 124.256818554513985, 127.516683879564406,
   129.578704200821853, 131.087688531430975, 133.497737202990660,
   134.756510050820649, 138.116042054533808, 139.736208952121808,
   141.123707404415728, 143.111845808910186, 146.000982487395827,
   147.422765343849989, 150.053520421293562, 150.925257612895526,
   153.024693791188948, 156.112909294982618, 157.597591818986345,
   158.849988365204885, 161.188964138954152, 163.030709687408168,
   165.537069188392498, 167.184439971994828, 169.094515416791259,
   169.911976498590630, 173.411536520135680, 174.754191523438771,
   176.441434188575954, 178.377407776468757, 179.916484018400656,
   182.207078484665730, 184.874467848130730, 185.598783678433914,
   187.228922291882030, 189.415759393773366, 192.026656325978780,
   193.079726604550355, 195.265396680495222, 196.876481841084053,
   198.015309585175508, 201.264751178782752
};

/* ══════════════════════════ B. Math functions ════════════════════════════════ */

static double fibonacci_real(double n) {
    return pow(PHI, n) / SQRT5 - pow(1.0/PHI, n) * cos(PI * n);
}

static double prime_product_index(double n, double beta) {
    int idx = ((int)floor(n + beta) + 50) % 50;
    return (double)PRIMES50[idx];
}

static double D_n(double n, double beta, double r, double k,
                  double Omega, double base) {
    double nb  = n + beta;
    double Fn  = fibonacci_real(nb);
    double Pn  = prime_product_index(n, beta);
    double dy  = pow(base, nb);
    double val = PHI * fmax(Fn, 1e-15) * dy * Pn * Omega;
    return sqrt(fmax(val, 1e-15)) * pow(fabs(r), k);
}

static double n_of_2p(uint64_t p) {
    double lx  = (double)p * M_LN2_V;
    double llx = log(lx / LN_PHI);
    if (llx <= 0.0) return -1.0;
    return llx / LN_PHI - 0.5 / PHI;
}

/* log10 of the Mersenne VALUE x = 2^p for lattice coordinate n */
static double log10_x_of_n(double n) {
    double inner     = n + 0.5 / PHI;
    double phi_inner = exp(inner * LN_PHI);   /* phi^(n + 1/(2*phi)) */
    return phi_inner * LOG10PHI;              /* log10(phi^phi_inner) */
}

static int phi_filter(uint64_t p) {
    double n = n_of_2p(p);
    if (n < 0.0) return 0;
    return (n - floor(n)) < 0.5;
}

static double gram_zero_k(int k) {
    if (k <= 0) return TWO_PI * (double)k;
    double z = (double)k * INV_E;
    double w = log(z + 1.0), ew;
    ew=exp(w); w-=(w*ew-z)/(ew*(w+1.0));
    ew=exp(w); w-=(w*ew-z)/(ew*(w+1.0));
    ew=exp(w); w-=(w*ew-z)/(ew*(w+1.0));
    ew=exp(w); w-=(w*ew-z)/(ew*(w+1.0));
    ew=exp(w); w-=(w*ew-z)/(ew*(w+1.0));
    ew=exp(w); w-=(w*ew-z)/(ew*(w+1.0));
    return TWO_PI * (double)k / w;
}

static double zeta_zero_cpu(int k) {
    return (k >= 0 && k < 80) ? ZETA_ZEROS_80[k] : gram_zero_k(k);
}

static double psi_score_cpu(double x, int B) {
    if (x <= 1.5) return 0.0;
    double lx  = log(x),  lxm = log(x - 1.0);
    double mx  = exp(0.5*lx), mm = exp(0.5*lxm);
    double px  = x, pm = x - 1.0;
    for (int k = 0; k < B; k++) {
        double t = zeta_zero_cpu(k);
        double d = 0.25 + t*t;
        px -= 2.0*mx*(0.5*cos(t*lx)  + t*sin(t*lx))  / d;
        pm -= 2.0*mm*(0.5*cos(t*lxm) + t*sin(t*lxm)) / d;
    }
    return px - pm;
}

/* ── Miller-Rabin 64-bit ─────────────────────────────────────────────────── */
#if defined(__GNUC__) && !defined(__clang__)
static uint64_t mulmod64(uint64_t a, uint64_t b, uint64_t m) {
    return (uint64_t)((__uint128_t)a * b % m);
}
#else
static uint64_t mulmod64(uint64_t a, uint64_t b, uint64_t m) {
    uint64_t r = 0; a %= m;
    while (b > 0) {
        if (b & 1) { r += a; if (r >= m) r -= m; }
        a = (a >= m - a) ? (a + a - m) : (a + a);
        b >>= 1;
    }
    return r;
}
#endif

static uint64_t powmod64(uint64_t base, uint64_t exp, uint64_t mod) {
    uint64_t r = 1; base %= mod;
    while (exp > 0) {
        if (exp & 1) r = mulmod64(r, base, mod);
        base = mulmod64(base, base, mod); exp >>= 1;
    }
    return r;
}

static int miller_rabin_witness(uint64_t n, uint64_t a) {
    if (n % a == 0) return (int)(n == a);
    uint64_t d = n - 1; int r = 0;
    while (!(d & 1)) { d >>= 1; r++; }
    uint64_t x = powmod64(a, d, n);
    if (x == 1 || x == n-1) return 1;
    for (int i = 0; i < r-1; i++) {
        x = mulmod64(x, x, n);
        if (x == n-1) return 1;
    }
    return 0;
}

static int is_prime_64(uint64_t n) {
    if (n < 2) return 0; if (n < 4) return 1;
    if (!(n & 1) || n % 3 == 0) return 0;
    static const uint64_t W[] = {2,3,5,7,11,13,17,19,23,29,31,37};
    for (int i = 0; i < 12; i++)
        if (!miller_rabin_witness(n, W[i])) return 0;
    return 1;
}

/* ── Segmented sieve ─────────────────────────────────────────────────────── */
static uint64_t *sieve_range(uint64_t lo, uint64_t hi, int *count) {
    *count = 0;
    if (hi < lo || (hi - lo) > 5000000ULL) return NULL;
    uint64_t cap  = 8*((hi-lo)/50 + 64);
    uint64_t *out = (uint64_t*)malloc(cap * sizeof(uint64_t));
    if (!out) return NULL;
    uint64_t span = hi - lo + 1;
    uint8_t *sv   = (uint8_t*)calloc((size_t)span, 1);
    if (!sv) { free(out); return NULL; }
    /* mark even composites */
    uint64_t s2 = (lo % 2 == 0) ? lo : lo + 1;
    for (uint64_t n = s2; n <= hi; n += 2) if (n > 2) sv[(size_t)(n-lo)] = 1;
    /* sieve small odd primes */
    uint64_t sq = (uint64_t)ceil(sqrt((double)hi)) + 1;
    for (uint64_t p = 3; p <= sq; p += 2) {
        int ok = 1;
        for (uint64_t d = 3; d*d <= p; d += 2) if (p%d==0) { ok=0; break; }
        if (!ok) continue;
        uint64_t first = ((lo + p - 1) / p) * p;
        if (first < p*p) first = p*p;
        if (first > hi) continue;
        for (uint64_t n = first; n <= hi; n += p) sv[(size_t)(n-lo)] = 1;
    }
    for (uint64_t n = lo; n <= hi; n++) {
        if (sv[(size_t)(n-lo)] || n < 2) continue;
        if ((size_t)*count >= (size_t)(cap-1)) {
            cap *= 2;
            uint64_t *tmp = (uint64_t*)realloc(out, cap*sizeof(uint64_t));
            if (!tmp) { free(sv); free(out); return NULL; }
            out = tmp;
        }
        out[(*count)++] = n;
    }
    free(sv);
    return out;
}

/* ── Dn ranking score for Mersenne candidate p ───────────────────────────── */
static double dn_score_mersenne(uint64_t p) {
    double n = n_of_2p(p);
    if (n < 0.0) return 0.0;
    double frac = n - floor(n);
    double k    = (n + 1.0) / 8.0;
    return D_n(n, 0.0, frac, k, 1.0, 2.0);
}

/* ══════════════════════════ C. Windows console helpers ══════════════════════ */

static HANDLE g_hout, g_hin;

static void console_init(void) {
    g_hout = GetStdHandle(STD_OUTPUT_HANDLE);
    g_hin  = GetStdHandle(STD_INPUT_HANDLE);
    DWORD mode;
    if (GetConsoleMode(g_hout, &mode))
        SetConsoleMode(g_hout, mode
            | ENABLE_VIRTUAL_TERMINAL_PROCESSING
            | ENABLE_PROCESSED_OUTPUT);
    SetConsoleOutputCP(CP_UTF8);
    SetConsoleCP(CP_UTF8);
}

/* ANSI color macros */
#define CR    "\033[0m"
#define BOLD  "\033[1m"
#define DIM   "\033[2m"
#define CYAN  "\033[1;36m"
#define YEL   "\033[1;33m"
#define GRN   "\033[1;32m"
#define RED   "\033[1;31m"
#define MAG   "\033[1;35m"
#define WHT   "\033[1;37m"

static double now_s(void) {
    LARGE_INTEGER f, c;
    QueryPerformanceFrequency(&f);
    QueryPerformanceCounter(&c);
    return (double)c.QuadPart / (double)f.QuadPart;
}

/* Readline using ReadConsoleW (cooked; handles echo + backspace) */
static int readline_prompt(const char *col, const char *prompt,
                           char *buf, int maxlen) {
    printf("%s%s" CR " ", col, prompt); fflush(stdout);
    DWORD old; GetConsoleMode(g_hin, &old);
    SetConsoleMode(g_hin, ENABLE_PROCESSED_INPUT|ENABLE_ECHO_INPUT|ENABLE_LINE_INPUT);
    WCHAR wb[512] = {0}; DWORD nr = 0;
    ReadConsoleW(g_hin, wb, 511, &nr, NULL);
    SetConsoleMode(g_hin, old);
    while (nr > 0 && (wb[nr-1]==L'\r'||wb[nr-1]==L'\n')) nr--;
    wb[nr] = 0;
    int n = WideCharToMultiByte(CP_UTF8, 0, wb, -1, buf, maxlen-1, NULL, NULL);
    if (n > 0) buf[n-1] = 0; else buf[0] = 0;
    return (int)strlen(buf);
}

static void wait_enter(void) {
    printf("\n  " DIM "[Enter to return to menu]" CR " "); fflush(stdout);
    DWORD old; GetConsoleMode(g_hin, &old);
    SetConsoleMode(g_hin, ENABLE_PROCESSED_INPUT|ENABLE_ECHO_INPUT|ENABLE_LINE_INPUT);
    WCHAR wb[4] = {0}; DWORD nr = 0;
    ReadConsoleW(g_hin, wb, 3, &nr, NULL);
    SetConsoleMode(g_hin, old);
}

/* ══════════════════════════ MODULE 1: Prime Pipeline ════════════════════════ */

typedef struct { uint64_t p; double nv; double frac; double dn; int pass; } Cand;

static int cmp_dn_desc(const void *a, const void *b) {
    double da = ((const Cand*)a)->dn, db = ((const Cand*)b)->dn;
    return (da > db) ? -1 : (da < db) ? 1 : 0;
}

static void module_pipeline(void) {
    printf("\n" CYAN "-- Prime Pipeline "
           "---------------------------------------------------\n" CR);
    printf("  Sieve prime exponents in [p_lo, p_hi], apply phi-filter,\n"
           "  compute Dn resonance score, rank Mersenne candidates.\n"
           "  Range limit: 5,000,000.\n\n");

    char buf[64]; uint64_t p_lo, p_hi;
    readline_prompt(YEL, "  p_lo", buf, sizeof(buf));
    if (!buf[0]) { printf("  Cancelled.\n"); return; }
    p_lo = strtoull(buf, NULL, 10);
    if (p_lo < 2) { printf(RED "  p_lo must be >= 2\n" CR); return; }

    readline_prompt(YEL, "  p_hi", buf, sizeof(buf));
    if (!buf[0]) { printf("  Cancelled.\n"); return; }
    p_hi = strtoull(buf, NULL, 10);

    if (p_hi < p_lo) { printf(RED "  p_hi < p_lo\n" CR); return; }
    if (p_hi - p_lo > 5000000ULL) {
        printf(RED "  Range > 5M — capped at p_lo + 5000000\n" CR);
        p_hi = p_lo + 5000000ULL;
    }

    printf("\n  Sieving [%llu, %llu]...\n",
           (unsigned long long)p_lo, (unsigned long long)p_hi);
    double t0 = now_s();
    int cnt = 0;
    uint64_t *primes = sieve_range(p_lo, p_hi, &cnt);
    if (!primes) { printf(RED "  Sieve failed.\n" CR); return; }
    printf("  Found %d primes in %.2f ms\n", cnt, (now_s()-t0)*1e3);
    if (cnt == 0) { free(primes); return; }

    Cand *cands = (Cand*)malloc(cnt * sizeof(Cand));
    if (!cands) { free(primes); return; }

    int phi_cnt = 0;
    for (int i = 0; i < cnt; i++) {
        uint64_t p = primes[i];
        double nv   = n_of_2p(p);
        double frac = (nv >= 0.0) ? nv - floor(nv) : 0.0;
        double dn   = dn_score_mersenne(p);
        int    pass = phi_filter(p);
        cands[i].p    = p;
        cands[i].nv   = nv;
        cands[i].frac = frac;
        cands[i].dn   = pass ? dn * 1.5 : dn;  /* lower-half bonus */
        cands[i].pass = pass;
        if (pass) phi_cnt++;
    }
    free(primes);
    qsort(cands, cnt, sizeof(Cand), cmp_dn_desc);

    readline_prompt(YEL, "  Top N to show (Enter = 20)", buf, sizeof(buf));
    int top = (buf[0] && atoi(buf) > 0) ? atoi(buf) : 20;
    if (top > cnt) top = cnt;

    printf("\n" CYAN
        "  +------+-----------+------------+-----------+-----------+------+\n"
        "  | Rank |     p     |   n(2^p)   |  frac(n)  |   Dn*     | phi  |\n"
        "  +------+-----------+------------+-----------+-----------+------+\n" CR);

    for (int i = 0; i < top; i++) {
        Cand *c = &cands[i];
        const char *pc = c->pass ? GRN : DIM;
        printf("  | %4d | %9llu | %10.6f | %9.6f | %9.3e | %s%s%s  |\n",
               i+1, (unsigned long long)c->p,
               c->nv, c->frac, c->dn,
               pc, c->pass ? "Y" : "n", CR);
    }
    printf(CYAN
        "  +------+-----------+------------+-----------+-----------+------+\n"
        CR);
    printf("  phi-filter: %d/%d passed (%.1f%%)\n"
           "  * Dn score x1.5 bonus for lower-half candidates.\n",
           phi_cnt, cnt, 100.0*phi_cnt/cnt);
    free(cands);
}

/* ══════════════════════════ MODULE 2: Number Analyzer ═══════════════════════ */

static void module_analyzer(void) {
    printf("\n" CYAN "-- Number Analyzer "
           "--------------------------------------------------\n" CR);
    printf("  Miller-Rabin primality, phi-lattice coordinate,\n"
           "  Dn resonance score, psi prime-counting score.\n\n");

    char buf[64];
    readline_prompt(YEL, "  n", buf, sizeof(buf));
    if (!buf[0]) { printf("  Cancelled.\n"); return; }
    uint64_t n = strtoull(buf, NULL, 10);
    if (n < 2) { printf(RED "  n must be >= 2\n" CR); return; }

    printf("\n  " WHT "n = %llu" CR "\n\n", (unsigned long long)n);

    /* Primality */
    int prime = is_prime_64(n);
    printf("  12-witness Miller-Rabin:  %s\n",
           prime ? GRN "PRIME" CR : RED "COMPOSITE" CR);

    /* phi-lattice */
    double nv = n_of_2p(n);
    if (nv >= 0.0) {
        double frac = nv - floor(nv);
        int pass = phi_filter(n);
        printf("  n(2^p) [p=%llu]:  %.8f\n", (unsigned long long)n, nv);
        printf("  frac(n):         %.8f  ->  phi-filter: %s\n",
               frac, pass ? GRN "PASS (lower half)" CR : RED "FAIL (upper half)" CR);
    }

    /* Dn score */
    double dn = dn_score_mersenne(n);
    printf("  Dn resonance score:  " YEL "%.6e" CR "\n", dn);

    /* psi-score */
    if (n <= 1000000000ULL) {
        double t0 = now_s();
        double psi = psi_score_cpu((double)n, 80);
        printf("  psi-score (B=80):    " YEL "%.6f" CR "  (%.1f ms)\n",
               psi, (now_s()-t0)*1e3);
        printf("  " DIM "(psi > 0.5 suggests prime; psi ~ 1.0 = strong signal)" CR "\n");
    } else {
        printf("  psi-score: skipped (n > 1e9)\n");
    }

    /* Factorization for small composites */
    if (!prime && n < 10000000ULL) {
        printf("  Factors: ");
        uint64_t m = n;
        for (uint64_t d = 2; d*d <= m && d < 100000; d++) {
            while (m % d == 0) {
                printf(MAG "%llu" CR " ", (unsigned long long)d); m /= d;
            }
        }
        if (m > 1) printf(MAG "%llu" CR, (unsigned long long)m);
        printf("\n");
    }

    /* Mersenne check */
    if (prime) {
        uint64_t mp = (n <= 62) ? ((1ULL << n) - 1ULL) : 0;
        if (mp > 0) {
            int mp_is = is_prime_64(mp);
            printf("  M_%llu = 2^%llu-1 = %llu  ->  %s\n",
                   (unsigned long long)n, (unsigned long long)n,
                   (unsigned long long)mp,
                   mp_is ? GRN "MERSENNE PRIME" CR : RED "composite" CR);
        } else {
            printf("  M_%llu = 2^%llu-1  (too large for 64-bit; use ll_analog.exe)\n",
                   (unsigned long long)n, (unsigned long long)n);
        }
    }
}

/* ══════════════════════════ MODULE 3: Mersenne Explorer ════════════════════ */

static void module_mersenne(void) {
    printf("\n" CYAN "-- Mersenne Explorer "
           "-------------------------------------------------\n" CR);
    printf("  All 51 known Mersenne primes M_p = 2^p - 1.\n\n");

    printf(CYAN
        "  +----+------------+------------+-----------+--------+-----------+\n"
        "  |  # |     p      |   n(2^p)   |  frac(n)  | phi-ok |  Dn score |\n"
        "  +----+------------+------------+-----------+--------+-----------+\n" CR);

    int phi_pass = 0;
    for (int i = 0; i < N_MERSENNE; i++) {
        uint64_t p = MERSENNE_EXP[i];
        double nv  = n_of_2p(p);
        double frac = (nv >= 0.0) ? nv - floor(nv) : -1.0;
        int    pass = phi_filter(p);
        double dn   = dn_score_mersenne(p);
        if (pass) phi_pass++;
        const char *pc = pass ? GRN : DIM;
        printf("  | %2d | %10llu | %10.6f | %9.6f |  %s%-3s%s   | %9.3e |\n",
               i+1, (unsigned long long)p, nv, frac,
               pc, pass ? "YES" : "no", CR, dn);
        /* page break every 25 rows */
        if (i == 24) {
            printf(CYAN "  ..." CR "\n");
            printf("  " DIM "[Enter for next page]" CR " "); fflush(stdout);
            DWORD old; GetConsoleMode(g_hin, &old);
            SetConsoleMode(g_hin, ENABLE_PROCESSED_INPUT|ENABLE_ECHO_INPUT|ENABLE_LINE_INPUT);
            WCHAR wb[4]={0}; DWORD nr=0;
            ReadConsoleW(g_hin, wb, 3, &nr, NULL);
            SetConsoleMode(g_hin, old);
            printf(CYAN
                "  +----+------------+------------+-----------+--------+-----------+\n"
                "  |  # |     p      |   n(2^p)   |  frac(n)  | phi-ok |  Dn score |\n"
                "  +----+------------+------------+-----------+--------+-----------+\n" CR);
        }
    }
    printf(CYAN
        "  +----+------------+------------+-----------+--------+-----------+\n" CR);
    printf("  phi-filter pass rate: " YEL "%d/%d = %.1f%%" CR
           "  (expected 50%% random; actual 67%%)\n\n",
           phi_pass, N_MERSENNE, 100.0*phi_pass/N_MERSENNE);

    /* Next-candidate predictions beyond M51 */
    double n51 = n_of_2p(136279841ULL);
    printf(CYAN "  -- Next Candidate Predictions (beyond M51, p=136279841) --\n" CR);
    printf("  M51: n(2^136279841) = %.6f   floor(n51) = %.0f\n\n", n51, floor(n51));
    printf("  Method: phi-lattice inverse  x(n) = phi^(phi^(n + 1/(2*phi)))\n"
           "          predicted p = 10^(log10_x) / log10(2)\n\n");
    printf(CYAN
        "  +-----+-----------+--------------------+------------------------+\n"
        "  |  +n |  lattice  |  predicted p       |  digits(M_p = 2^p-1)   |\n"
        "  +-----+-----------+--------------------+------------------------+\n" CR);

    for (int step = 1; step <= 14; step++) {
        double target_n  = floor(n51) + (double)step;
        double log10_x   = log10_x_of_n(target_n); /* log10(2^p) = p*log10(2) */
        double p_est     = log10_x / LOG10_2;       /* predicted exponent */
        double digits    = log10_x;                  /* digits of M_p */
        double log10_p   = log10(p_est);

        char p_str[32], d_str[32];
        if (log10_p < 9.0)
            snprintf(p_str, sizeof(p_str), "%llu", (unsigned long long)p_est);
        else
            snprintf(p_str, sizeof(p_str), "~10^%.3f", log10_p);
        snprintf(d_str, sizeof(d_str), "~%.4e", digits);

        printf("  | %+3d | %9.4f | %-18s | %-22s |\n",
               step, target_n, p_str, d_str);
    }
    printf(CYAN
        "  +-----+-----------+--------------------+------------------------+\n" CR);
    printf("  " DIM "Predictions are phi-lattice resonance points, not proofs.\n"
           "  Use prime_pipeline.exe + ll_cuda.exe to test specific ranges." CR "\n");
}

/* ══════════════════════════ MODULE 4: Zeta Zeros ════════════════════════════ */

static void module_zeta(void) {
    printf("\n" CYAN "-- Zeta Zeros "
           "--------------------------------------------------------\n" CR);
    printf("  Non-trivial zeros of zeta(s) on critical line s = 1/2 + it.\n"
           "  k=0..79 exact; k>=80 Gram approximation (6-iter Newton/W).\n\n");

    char buf[16];
    readline_prompt(YEL, "  Show k=0..K (Enter = 19)", buf, sizeof(buf));
    int K = (buf[0] && atoi(buf) >= 0) ? atoi(buf) : 19;
    if (K > 499) K = 499;

    printf("\n" CYAN
        "  +------+------------------------------+-----------------+\n"
        "  |   k  |  t_k  (imaginary part)       |  source         |\n"
        "  +------+------------------------------+-----------------+\n" CR);

    for (int k = 0; k <= K; k++) {
        double t = zeta_zero_cpu(k);
        const char *src = (k < 80) ? "exact" : "Gram approx";
        const char *col = (k < 80) ? GRN    : YEL;
        printf("  | %4d | %s%28.12f" CR " | %-15s |\n", k, col, t, src);
    }
    printf(CYAN
        "  +------+------------------------------+-----------------+\n" CR);
    printf("  psi_score_cpu(x, B) uses zeros k=0..B-1.  B=80 fast, B=500 accurate.\n");
}

/* ══════════════════════════ MODULE 5: Benchmark ═════════════════════════════ */

static void print_bench_row(const char *name, double us_total, long N) {
    double per = us_total / (double)N;
    const char *unit; double val;
    if      (per < 1.0)    { unit = "ns"; val = per * 1000.0; }
    else if (per < 1000.0) { unit = "us"; val = per; }
    else                   { unit = "ms"; val = per / 1000.0; }
    printf("  %-42s" YEL "%8.3f %s" CR "   (%ldK)\n",
           name, val, unit, N/1000);
}

static void module_benchmark(void) {
    printf("\n" CYAN "-- Benchmark "
           "---------------------------------------------------------\n" CR);
    printf("  Quick timing of all 13 prime library functions.\n\n");
    printf("  Running " YEL "13" CR " benchmarks...\n\n");

    volatile double acc = 0.0;
    double t0, t1; long N;

    N = 2000000; t0 = now_s();
    for (long i=0;i<N;i++) acc += fibonacci_real((double)(i%64));
    t1 = now_s(); print_bench_row("1.  fibonacci_real(n)", (t1-t0)*1e6, N);

    N = 2000000; t0 = now_s();
    for (long i=0;i<N;i++) acc += prime_product_index((double)(i%32), 0.5);
    t1 = now_s(); print_bench_row("2.  prime_product_index(n, beta)", (t1-t0)*1e6, N);

    N = 1000000; t0 = now_s();
    for (long i=0;i<N;i++) acc += D_n((double)(i%16),0.0,0.618,0.375,1.0,2.0);
    t1 = now_s(); print_bench_row("3.  D_n(n,beta,r,k,Omega,base)", (t1-t0)*1e6, N);

    N = 2000000; t0 = now_s();
    for (long i=0;i<N;i++) acc += n_of_2p((uint64_t)(i%10000+2));
    t1 = now_s(); print_bench_row("4.  n_of_2p(p)", (t1-t0)*1e6, N);

    N = 2000000; t0 = now_s();
    for (long i=0;i<N;i++) acc += phi_filter((uint64_t)(i%10000+2));
    t1 = now_s(); print_bench_row("5.  phi_filter(p)", (t1-t0)*1e6, N);

    N = 500000; t0 = now_s();
    for (long i=0;i<N;i++) acc += gram_zero_k(i%10000+1);
    t1 = now_s(); print_bench_row("6.  gram_zero_k(k)  [6-iter Newton]", (t1-t0)*1e6, N);

    N = 500000; t0 = now_s();
    for (long i=0;i<N;i++) acc += zeta_zero_cpu(i%10000);
    t1 = now_s(); print_bench_row("7.  zeta_zero_cpu(k)  [table+Gram]", (t1-t0)*1e6, N);

    N = 2000; t0 = now_s();
    for (long i=0;i<N;i++) acc += psi_score_cpu(1000000.0+(double)i, 80);
    t1 = now_s(); print_bench_row("8.  psi_score_cpu(x, B=80)", (t1-t0)*1e6, N);

    N = 200000; t0 = now_s();
    for (long i=0;i<N;i++) acc += is_prime_64((uint64_t)(i*6+5));
    t1 = now_s(); print_bench_row("9.  is_prime_64  [12-witness MR]", (t1-t0)*1e6, N);

    N = 2000000; t0 = now_s();
    for (long i=0;i<N;i++) acc += dn_score_mersenne(MERSENNE_EXP[i%N_MERSENNE]);
    t1 = now_s(); print_bench_row("10. dn_score_mersenne(p)", (t1-t0)*1e6, N);

    /* sieve */
    t0 = now_s();
    int sc = 0;
    uint64_t *sv = sieve_range(10000000ULL, 10100000ULL, &sc);
    if (sv) free(sv);
    t1 = now_s();
    printf("  %-42s" YEL "%8.3f ms" CR "  (%d primes)\n",
           "11. sieve_range([1e7, 1e7+100K])", (t1-t0)*1e3, sc);

    N = 100000; t0 = now_s();
    for (long i=0;i<N;i++)
        for (int j=0;j<N_MERSENNE;j++) acc += phi_filter(MERSENNE_EXP[j]);
    t1 = now_s(); print_bench_row("12. phi_filter_batch (51 exp x100K)", (t1-t0)*1e6, N*N_MERSENNE);

    N = 100000; t0 = now_s();
    for (long i=0;i<N;i++) {
        acc += zeta_zero_cpu(i%200);
        acc += dn_score_mersenne((uint64_t)(i*7+13));
    }
    t1 = now_s(); print_bench_row("13. Gram + Dn pipeline (100K)", (t1-t0)*1e6, N);

    printf("\n  " DIM "volatile acc = %.3f  (prevents dead-code elimination)" CR "\n", acc);
}

/* ── Forward declarations for lattice state (defined in MODULE 7) ─────────── */
#define LATTICE_MAX 4096
static double lattice[LATTICE_MAX];
static int    lattice_N                = 4096;
static int    lattice_alpine_installed = 0;
static int    lattice_seed_steps_done  = 0;
static void   lattice_seed_phi(int N, int steps);   /* defined in MODULE 7 */

/* ══════════════════════════ MODULE 6: Alpine Install ════════════════════════
 *
 *  Alpine Install: two-stage hook that marries prime_ui.exe with the
 *  bootloader + conscious pipeline.
 *
 *  Stage 1 -- bootloaderZ.exe  : initialize Slot4096 APA lattice on-disk
 *  Stage 2 -- conscious.exe    : GPU phi-resonance classifier on the lattice
 *
 *  Both are launched as child processes; stdout is captured line-by-line
 *  and rendered inside the TUI.
 * ══════════════════════════════════════════════════════════════════════════ */

/* Run a command, stream its output into the TUI indented, return exit code. */
static int run_capture(const char *cmd) {
    FILE *fp = _popen(cmd, "r");
    if (!fp) {
        printf("  " RED "[error] _popen failed: %s\n" CR, cmd);
        return -1;
    }
    char line[512];
    while (fgets(line, (int)sizeof(line), fp)) {
        int len = (int)strlen(line);
        while (len > 0 && (line[len-1] == '\n' || line[len-1] == '\r'))
            line[--len] = '\0';
        printf("    %s\n", line);
    }
    fflush(stdout);
    return _pclose(fp);
}

static void module_alpine(void) {
    printf("\n" CYAN
        "+--------------------------------------------------------------+\n"
        "|  Alpine Install  --  Slot4096 Lattice Boot + GPU Resonance   |\n"
        "+--------------------------------------------------------------+\n"
        CR "\n");
    printf("  Stage 1: " YEL "bootloaderZ.exe" CR
           "  -- initialize APA phi-lattice (Slot4096)\n"
           "  Stage 2: " YEL "conscious.exe" CR
           "   -- GPU dual-slot lambda-sigma resonance classifier\n\n");

    /* probe binaries next to the exe */
    int have_boot = (GetFileAttributesA("bootloaderZ.exe") != INVALID_FILE_ATTRIBUTES);
    int have_con  = (GetFileAttributesA("conscious.exe")   != INVALID_FILE_ATTRIBUTES);
    printf("  bootloaderZ.exe : %s\n",
           have_boot ? GRN "found"     CR : YEL "not found (self-seed fallback)" CR);
    printf("  conscious.exe   : %s\n\n",
           have_con  ? GRN "found"     CR : RED "NOT FOUND"                      CR);

    if (!have_con) {
        printf("  " RED "Cannot proceed -- run build_conscious.bat first.\n" CR);
        return;
    }

    /* parameters */
    char buf[64];
    uint64_t N  = 8192;
    int   steps = 1024;
    char  seed_arg[40] = {0};

    readline_prompt(YEL, "  N slots   [default 8192]", buf, (int)sizeof(buf));
    if (buf[0]) N = strtoull(buf, NULL, 10);
    readline_prompt(YEL, "  steps     [default 1024]", buf, (int)sizeof(buf));
    if (buf[0]) steps = (int)strtol(buf, NULL, 10);
    readline_prompt(YEL, "  seed hex  [blank = random]", buf, (int)sizeof(buf));
    if (buf[0]) snprintf(seed_arg, sizeof(seed_arg), "--seed %s", buf);

    /* ── Stage 1: bootloader ─────────────────────────────────────────── */
    if (have_boot) {
        printf("\n" CYAN "  ---- Stage 1: bootloaderZ ----\n" CR);
        fflush(stdout);
        int rc1 = run_capture("bootloaderZ.exe 2>&1");
        if (rc1 == 0)
            printf("  " GRN "[OK] Lattice initialized.\n" CR);
        else
            printf("  " YEL "[warn] bootloaderZ exited %d -- "
                   "conscious will self-seed.\n" CR, rc1);
    } else {
        printf("\n  " YEL "[skip] bootloaderZ.exe absent -- "
               "conscious will self-seed Slot4096.\n" CR);
    }

    /* ── Stage 2: conscious ──────────────────────────────────────────── */
    printf("\n" CYAN "  ---- Stage 2: conscious (GPU resonance) ----\n" CR);
    fflush(stdout);

    char cmd[256];
    snprintf(cmd, sizeof(cmd),
             "conscious.exe --N %llu --steps %d %s 2>&1",
             (unsigned long long)N, steps, seed_arg);
    printf("  " DIM "%s\n" CR "\n", cmd);

    int rc2 = run_capture(cmd);
    printf("\n");
    if (rc2 == 0) {
        printf("  " GRN "[OK] Alpine install complete -- "
               "lattice live, resonance verdict above.\n" CR);
        /* Seed the in-memory Slot4096 lattice so the Lattice Shell [7]
         * reflects the post-Alpine state (matches bootloader_init_lattice). */
        lattice_seed_phi(4096, 50);
        printf("  " GRN "[OK] In-memory lattice seeded: %d slots, 50 resonance steps.\n" CR,
               lattice_N);
    } else {
        printf("  " RED "[FAILED] conscious.exe exited %d.\n" CR, rc2);
        /* Even if GPU stage failed, bootloaderZ may have initialised the
         * lattice structure -- seed in-memory state from it regardless. */
        if (have_boot) {
            lattice_seed_phi(4096, 50);
            printf("  " YEL "[warn] GPU stage failed; in-memory lattice seeded from"
                   " bootloaderZ state only.\n" CR);
        }
    }
}

/* ══════════════════════════ MODULE 7: Lattice Shell ════════════════════════ */

static void lattice_reset(int N) {
    if (N > LATTICE_MAX) N = LATTICE_MAX;
    lattice_N = N;
    for (int i = 0; i < N; ++i)
        lattice[i] = 0.0;
    lattice_alpine_installed = 0;
    lattice_seed_steps_done  = 0;
}

static void lattice_step(void) {
    /* Simple resonance step: phi-lattice update */
    for (int i = 0; i < lattice_N; ++i)
        lattice[i] = PHI * (lattice[i] + 1.0) - floor(PHI * (lattice[i] + 1.0));
}

/*
 * lattice_seed_phi -- mirror of bootloader_init_lattice():
 *   slot[i] = frac(i * phi)  then advance `steps` resonance ticks.
 * Called automatically after a successful Alpine Install.
 */
static void lattice_seed_phi(int N, int steps) {
    if (N < 1) N = 4096;
    if (N > LATTICE_MAX) N = LATTICE_MAX;
    lattice_N = N;
    /* phi-irrational uniform spacing (van-der-Corput / Weyl sequence) */
    for (int i = 0; i < N; ++i)
        lattice[i] = (i * PHI) - floor(i * PHI);
    /* seed resonance steps */
    for (int s = 0; s < steps; ++s)
        lattice_step();
    lattice_alpine_installed = 1;
    lattice_seed_steps_done  = steps;
}

static void lattice_stats(void) {
    double min = lattice[0], max = lattice[0], sum = 0.0;
    for (int i = 0; i < lattice_N; ++i) {
        if (lattice[i] < min) min = lattice[i];
        if (lattice[i] > max) max = lattice[i];
        sum += lattice[i];
    }
    printf("  Lattice size: %d\n", lattice_N);
    printf("  Min: %.6f  Max: %.6f  Mean: %.6f\n", min, max, sum/lattice_N);
}

static void lattice_query(int idx) {
    if (idx < 0 || idx >= lattice_N) {
        printf("  [error] Index out of range (0-%d)\n", lattice_N-1);
        return;
    }
    printf("  lattice[%d] = %.8f\n", idx, lattice[idx]);
}

static void module_lattice_shell(void) {
    printf("\n" CYAN
        "+--------------------------------------------------------------+\n"
        "|  Lattice Shell  --  Interactive Slot4096 Lattice REPL       |\n"
        "+--------------------------------------------------------------+\n"
        CR "\n");

    /* Show Alpine status instead of always wiping the lattice */
    if (lattice_alpine_installed) {
        printf("  " GRN "[Alpine] Lattice live: %d slots, seeded with %d resonance steps.\n"
               CR "\n", lattice_N, lattice_seed_steps_done);
    } else {
        /* First entry: blank initialise */
        if (lattice_N == 0) lattice_reset(4096);
        printf("  " YEL "[Note] Lattice not Alpine-seeded. "
               "Run Alpine Install [6] first, or type 'install' / 'seed [steps]'.\n"
               CR "\n");
    }
    printf("  Commands: help, status, stats, step [n], reset [N], query <i>,\n"
           "            seed [steps], install, exit\n\n");

    char buf[128];
    for (;;) {
        printf(WHT "lattice> " CR); fflush(stdout);
        readline_prompt(WHT, "", buf, sizeof(buf));
        char *cmd = strtok(buf, " \t\r\n");
        if (!cmd) continue;

        if (!strcmp(cmd, "help")) {
            printf("  help           Show this help message\n");
            printf("  status         Show Alpine install / seeding status\n");
            printf("  stats          Show lattice statistics\n");
            printf("  step [n]       Advance lattice n resonance steps (default 1)\n");
            printf("  reset [N]      Reset lattice to zero (default N=4096)\n");
            printf("  query <i>      Show value at slot i\n");
            printf("  seed [steps]   Phi-seed lattice (default 50 steps, like Alpine)\n");
            printf("  install        Run Alpine Install pipeline + seed lattice\n");
            printf("  exit           Leave the shell\n");

        } else if (!strcmp(cmd, "status")) {
            if (lattice_alpine_installed)
                printf("  " GRN "[Alpine] Installed -- %d slots, %d seed steps.\n" CR,
                       lattice_N, lattice_seed_steps_done);
            else
                printf("  " YEL "[Not installed] Lattice is blank or manually reset.\n" CR);

        } else if (!strcmp(cmd, "stats")) {
            lattice_stats();

        } else if (!strcmp(cmd, "step")) {
            char *nstr = strtok(NULL, " \t\r\n");
            int n = nstr ? atoi(nstr) : 1;
            if (n < 1) n = 1;
            for (int s = 0; s < n; ++s) lattice_step();
            printf("  [OK] Lattice advanced %d step(s).\n", n);

        } else if (!strcmp(cmd, "reset")) {
            char *nstr = strtok(NULL, " \t\r\n");
            int N = nstr ? atoi(nstr) : 4096;
            lattice_reset(N);
            printf("  [OK] Lattice reset to zero, N=%d.\n", lattice_N);

        } else if (!strcmp(cmd, "query")) {
            char *istr = strtok(NULL, " \t\r\n");
            if (!istr) { printf("  Usage: query <i>\n"); continue; }
            int idx = atoi(istr);
            lattice_query(idx);

        } else if (!strcmp(cmd, "seed")) {
            char *sstr = strtok(NULL, " \t\r\n");
            int steps = sstr ? atoi(sstr) : 50;
            if (steps < 1) steps = 1;
            lattice_seed_phi(lattice_N, steps);
            printf("  " GRN "[OK] Lattice phi-seeded: %d slots, %d steps.\n" CR,
                   lattice_N, steps);

        } else if (!strcmp(cmd, "install")) {
            /* Run the full Alpine install pipeline inline, then seed lattice */
            printf("\n" CYAN "  ---- Alpine Install (inline) ----\n" CR);
            int have_boot = (GetFileAttributesA("bootloaderZ.exe") != INVALID_FILE_ATTRIBUTES);
            int have_con  = (GetFileAttributesA("conscious.exe")   != INVALID_FILE_ATTRIBUTES);
            printf("  bootloaderZ.exe : %s\n",
                   have_boot ? GRN "found" CR : YEL "not found" CR);
            printf("  conscious.exe   : %s\n",
                   have_con  ? GRN "found" CR : RED "NOT FOUND" CR);
            if (have_boot) {
                printf("  " CYAN "-- Stage 1: bootloaderZ --\n" CR);
                int rc1 = run_capture("bootloaderZ.exe 2>&1");
                if (rc1 != 0)
                    printf("  " YEL "[warn] bootloaderZ exited %d\n" CR, rc1);
            }
            if (have_con) {
                printf("  " CYAN "-- Stage 2: conscious --\n" CR);
                int rc2 = run_capture("conscious.exe --N 4096 --steps 1024 2>&1");
                if (rc2 != 0)
                    printf("  " YEL "[warn] conscious exited %d\n" CR, rc2);
            } else {
                printf("  " YEL "[skip] conscious.exe not found -- seeding from bootloaderZ only.\n" CR);
            }
            lattice_seed_phi(4096, 50);
            printf("  " GRN "[OK] In-memory lattice seeded: %d slots, 50 resonance steps.\n" CR,
                   lattice_N);

        } else if (!strcmp(cmd, "exit")) {
            printf("  Leaving lattice shell.\n");
            break;

        } else {
            printf("  [error] Unknown command: %s  (type 'help')\n", cmd);
        }
    }
}

/* ══════════════════════════ MODULE 8: Alpine OS Shell ══════════════════════
 *
 * The Slot4096 lattice IS the OS substrate.  The lattice determines:
 *   - Entropy  : all N slot IEEE-754 values fed raw into /dev/urandom
 *   - Hostname : phi-fold of slots[0..3] → hex string
 *   - APK mirror: slot[4] selects which CDN serves packages
 *   - Packages  : slots[5..15] bit-select which tools get installed
 *   - Timezone  : slot[16] maps to UTC offset
 *   - CPU nice  : slot[17] → process priority for all children
 *   - PS1 prompt: resonance seed embedded in shell prompt
 *   - MOTD      : full resonance fingerprint shown on login
 *   - /etc/profile.d/lattice.sh : LATTICE_* exported to every shell
 *   - /lattice/slots/ : individual slot values as readable files
 * ══════════════════════════════════════════════════════════════════════════ */

/* ── Derivation helpers ─────────────────────────────────────────────────── */

static uint64_t lattice_derive_seed(void) {
    uint64_t seed = 0xC0FFEE00DEAD1234ULL;
    int n = lattice_N < 64 ? lattice_N : 64;
    for (int i = 0; i < n; ++i) {
        union { double d; uint64_t u; } cv;
        cv.d = lattice[i];
        seed ^= cv.u ^ ((uint64_t)(i + 1) * 6364136223846793005ULL);
        seed  = (seed << 31) | (seed >> 33);
        seed *= 0x9e3779b97f4a7c15ULL;
    }
    return seed;
}

static void lattice_derive_hostname(char *buf, int bufsz) {
    /* XOR-fold first 4 slots into 32 bits → phi4096-XXXXXXXX */
    uint32_t h = 0;
    int n = lattice_N < 4 ? lattice_N : 4;
    for (int i = 0; i < n; ++i) {
        union { double d; uint64_t u; } cv; cv.d = lattice[i];
        h ^= (uint32_t)(cv.u >> 32) ^ (uint32_t)(cv.u & 0xFFFFFFFF);
        h  = (h << 13) | (h >> 19);
        h *= 0x9e3779b9u;
    }
    snprintf(buf, (size_t)bufsz, "phi4096-%08x", h);
}

static const char *lattice_derive_mirror(void) {
    static const char *mirrors[] = {
        "dl-cdn.alpinelinux.org",
        "mirror.leaseweb.net",
        "mirrors.dotsrc.org",
        "mirrors.edge.kernel.org",
        "alpine.global.ssl.fastly.net",
        "ftp.halifax.rwth-aachen.de",
        "mirror.yandex.ru",
        "mirror.internode.on.net"
    };
    double v = (lattice_N > 4) ? lattice[4] : 0.0;
    if (v < 0.0) v = 0.0; if (v >= 1.0) v = 0.999;
    return mirrors[(int)(v * 8)];
}

/* Returns space-separated package names derived from slots[5..15] */
static void lattice_derive_packages(char *buf, int bufsz) {
    static const char *pkgs[] = {
        "python3", "curl", "htop", "bash", "vim",
        "jq",      "git",  "openssl", "bc",  "ncurses", "file"
    };
    buf[0] = '\0';
    for (int i = 0; i < 11; ++i) {
        double v = (lattice_N > 5 + i) ? lattice[5 + i] : 0.0;
        if (v > 0.5) {
            if (buf[0]) strncat(buf, " ", (size_t)(bufsz - (int)strlen(buf) - 1));
            strncat(buf, pkgs[i],        (size_t)(bufsz - (int)strlen(buf) - 1));
        }
    }
    /* always install at least these baseline tools */
    if (!strstr(buf, "curl")) strncat(buf, " curl",     (size_t)(bufsz - (int)strlen(buf) - 1));
    if (!strstr(buf, "bash")) strncat(buf, " bash",     (size_t)(bufsz - (int)strlen(buf) - 1));
}

/* UTC timezone string derived from slot[16]: maps [0,1) → UTC-11..UTC+14 */
static void lattice_derive_tz(char *buf, int bufsz) {
    double v = (lattice_N > 16) ? lattice[16] : 0.5;
    if (v < 0.0) v = 0.0; if (v >= 1.0) v = 0.999;
    int offset = (int)(v * 26) - 11;   /* -11 .. +14 */
    if (offset == 0)
        snprintf(buf, (size_t)bufsz, "UTC");
    else if (offset > 0)
        snprintf(buf, (size_t)bufsz, "Etc/GMT-%d", offset);
    else
        snprintf(buf, (size_t)bufsz, "Etc/GMT+%d", -offset);
}

/* nice value [−20..19] derived from slot[17] */
static int lattice_derive_nice(void) {
    double v = (lattice_N > 17) ? lattice[17] : 0.5;
    if (v < 0.0) v = 0.0; if (v >= 1.0) v = 0.999;
    return (int)(v * 40) - 20;
}

/* UID for slot4096 user: slot[18] → 1000..9999 */
static int lattice_derive_uid(void) {
    double v = (lattice_N > 18) ? lattice[18] : 0.5;
    if (v < 0.0) v = 0.0; if (v >= 1.0) v = 0.999;
    return 1000 + (int)(v * 8999);
}

/* ulimit -n open files: slot[19] → 1024..65536, rounded to 1024 */
static int lattice_derive_nofile(void) {
    double v = (lattice_N > 19) ? lattice[19] : 0.5;
    if (v < 0.0) v = 0.0; if (v >= 1.0) v = 0.999;
    int n = 1024 + (int)(v * 64512);
    return (n / 1024) * 1024;
}

/* umask: slot[20] → one of 002, 007, 022, 027 */
static int lattice_derive_umask(void) {
    static const int masks[] = { 002, 007, 022, 027 };
    double v = (lattice_N > 20) ? lattice[20] : 0.5;
    if (v < 0.0) v = 0.0; if (v >= 1.0) v = 0.999;
    return masks[(int)(v * 4)];
}

/* HISTSIZE: slot[21] → 500..9500 */
static int lattice_derive_histsize(void) {
    double v = (lattice_N > 21) ? lattice[21] : 0.5;
    if (v < 0.0) v = 0.0; if (v >= 1.0) v = 0.999;
    return 500 + (int)(v * 9000);
}

/* Shell inactivity timeout TMOUT: slot[22] → 300..3600 s */
static int lattice_derive_tmout(void) {
    double v = (lattice_N > 22) ? lattice[22] : 0.5;
    if (v < 0.0) v = 0.0; if (v >= 1.0) v = 0.999;
    return 300 + (int)(v * 3300);
}

/* ── File writers ───────────────────────────────────────────────────────── */

/* Write all N slot double values as raw IEEE-754 bytes → entropy pool */
static void lattice_write_entropy_blob(const char *path) {
    FILE *f = fopen(path, "wb");
    if (!f) return;
    fwrite(lattice, sizeof(double), (size_t)lattice_N, f);
    fclose(f);
}

/* Write full lattice state as shell-sourceable env file */
static void lattice_write_state_env(const char *path) {
    FILE *f = fopen(path, "wb");  /* binary mode = Unix LF on Windows */
    if (!f) return;
    uint64_t seed = lattice_derive_seed();
    char hostname[64]; lattice_derive_hostname(hostname, sizeof(hostname));
    char tz[32];       lattice_derive_tz(tz, sizeof(tz));
    fprintf(f, "LATTICE_N=%d\n",               lattice_N);
    fprintf(f, "LATTICE_STEPS=%d\n",            lattice_seed_steps_done);
    fprintf(f, "LATTICE_INSTALLED=%d\n",        lattice_alpine_installed);
    fprintf(f, "LATTICE_SEED=0x%016llx\n",      (unsigned long long)seed);
    fprintf(f, "LATTICE_HOSTNAME=%s\n",          hostname);
    fprintf(f, "LATTICE_MIRROR=%s\n",            lattice_derive_mirror());
    fprintf(f, "LATTICE_TZ=%s\n",               tz);
    fprintf(f, "LATTICE_NICE=%d\n",             lattice_derive_nice());
    fprintf(f, "LATTICE_UID=%d\n",              lattice_derive_uid());
    fprintf(f, "LATTICE_NOFILE=%d\n",           lattice_derive_nofile());
    fprintf(f, "LATTICE_UMASK=%04o\n",          lattice_derive_umask());
    fprintf(f, "LATTICE_HISTSIZE=%d\n",         lattice_derive_histsize());
    fprintf(f, "LATTICE_TMOUT=%d\n",            lattice_derive_tmout());
    /* all slot values */
    for (int i = 0; i < lattice_N; ++i)
        fprintf(f, "SLOT_%d=%.16f\n", i, lattice[i]);
    fclose(f);
}

/*
 * Write the full lattice_init.sh that IS the Alpine boot sequence.
 * Every OS parameter comes from the lattice — nothing is default.
 */
static void lattice_write_init_sh(const char *sh_path,
                                   const char *ent_linux_path) {
    FILE *f = fopen(sh_path, "wb");  /* binary mode = Unix LF so Alpine sh parses it */
    if (!f) return;

    uint64_t seed     = lattice_derive_seed();
    char hostname[64]; lattice_derive_hostname(hostname, sizeof(hostname));
    const char *mirror = lattice_derive_mirror();
    char pkgs[256];    lattice_derive_packages(pkgs, sizeof(pkgs));
    char tz[32];       lattice_derive_tz(tz, sizeof(tz));
    int  nice_val      = lattice_derive_nice();
    int  uid_val       = lattice_derive_uid();
    int  nofile_val    = lattice_derive_nofile();
    int  umask_val     = lattice_derive_umask();
    int  histsize_val  = lattice_derive_histsize();
    int  tmout_val     = lattice_derive_tmout();

    fputs("#!/bin/sh\n"
          "# lattice_init.sh -- generated by Slot4096 resonance engine\n"
          "# Every parameter below was derived from the live phi-lattice state.\n"
          "# This script IS the OS boot sequence.\n\n", f);

    /* 1. Entropy -- hardcoded bind-mount path, not $LATTICE_ENT */
    fprintf(f, "# -- 1. Entropy: feed all %d slot values into /dev/urandom --\n", lattice_N);
    fprintf(f, "if [ -f '%s' ]; then\n"
               "  dd if='%s' of=/dev/urandom bs=4096 2>/dev/null\n"
               "  echo \"[lattice] Entropy: $(wc -c < '%s') bytes -> /dev/urandom\"\n"
               "fi\n\n", ent_linux_path, ent_linux_path, ent_linux_path);

    /* 2. Hostname + /etc/hostname file */
    fprintf(f, "# -- 2. Hostname (phi-fold of slots[0..3]) --\n"
               "hostname '%s'\n"
               "echo '%s' > /etc/hostname\n"
               "echo '[lattice] Hostname:  %s'\n\n", hostname, hostname, hostname);

    /* 3. Timezone */
    fprintf(f, "# -- 3. Timezone (slot[16] -> %s) --\n"
               "ln -sf /usr/share/zoneinfo/%s /etc/localtime 2>/dev/null\n"
               "echo '%s' > /etc/timezone\n"
               "echo '[lattice] Timezone:  %s'\n\n", tz, tz, tz, tz);

    /* 4. APK mirror */
    fprintf(f, "# -- 4. APK mirror (slot[4] -> %s) --\n"
               "cat > /etc/apk/repositories << 'EOAPK'\n"
               "https://%s/alpine/latest-stable/main\n"
               "https://%s/alpine/latest-stable/community\n"
               "EOAPK\n"
               "echo '[lattice] APK mirror: %s'\n\n",
               mirror, mirror, mirror, mirror);

    /* 5. Packages */
    fprintf(f, "# -- 5. Packages (slots[5..15] bit-select) --\n"
               "echo '[lattice] Installing: %s'\n"
               "apk update -q 2>/dev/null && apk add -q %s 2>/dev/null\n"
               "echo '[lattice] Packages ready.'\n\n", pkgs, pkgs);

    /* 6. CPU nice */
    fprintf(f, "# -- 6. Process priority (slot[17] -> nice %+d) --\n"
               "renice %d $$ >/dev/null 2>&1\n"
               "echo '[lattice] Nice: %+d'\n\n", nice_val, nice_val, nice_val);

    /* 7. slot4096 user (slot[18] -> UID %d) */
    fprintf(f, "# -- 7. slot4096 user (slot[18] -> UID %d) --\n"
               "addgroup -g %d slot4096 2>/dev/null\n"
               "adduser -u %d -G slot4096 -h /home/slot4096 -s /bin/bash -D slot4096 2>/dev/null\n"
               "mkdir -p /home/slot4096\n"
               "chown slot4096:slot4096 /home/slot4096\n"
               "echo '[lattice] User: slot4096 uid=%d'\n\n",
               uid_val, uid_val, uid_val, uid_val);

    /* 8. /etc/profile.d/lattice.sh — exported to every login shell */
    fprintf(f, "# -- 8. Profile: export LATTICE_* + ulimit/umask to all shells --\n"
               "mkdir -p /etc/profile.d\n"
               "cat > /etc/profile.d/lattice.sh << 'EOPROFILE'\n"
               "export LATTICE_N=%d\n"
               "export LATTICE_STEPS=%d\n"
               "export LATTICE_SEED=0x%016llx\n"
               "export LATTICE_HOSTNAME=%s\n"
               "export LATTICE_MIRROR=%s\n"
               "export LATTICE_TZ=%s\n"
               "export LATTICE_NICE=%d\n"
               "export LATTICE_UID=%d\n"
               "export LATTICE_NOFILE=%d\n"
               "export LATTICE_UMASK=%04o\n"
               "export LATTICE_HISTSIZE=%d\n"
               "export LATTICE_TMOUT=%d\n"
               "umask %04o\n"
               "ulimit -n %d 2>/dev/null\n"
               "export HISTSIZE=%d\n"
               "export HISTFILESIZE=%d\n"
               "export TMOUT=%d\n"
               "export PS1='\\[\\033[1;35m\\]slot4096@%s:\\[\\033[1;32m\\]\\w\\[\\033[0m\\]$ '\n"
               "EOPROFILE\n\n",
               lattice_N, lattice_seed_steps_done,
               (unsigned long long)seed,
               hostname, mirror, tz,
               nice_val, uid_val, nofile_val,
               umask_val, histsize_val, tmout_val,
               umask_val, nofile_val,
               histsize_val, histsize_val * 2, tmout_val,
               hostname);

    /* 9. ~/.bashrc for slot4096 — full lattice env on interactive login */
    fprintf(f, "# -- 9. ~/.bashrc for slot4096 --\n"
               "cat > /home/slot4096/.bashrc << 'EOBASHRC'\n"
               ". /etc/profile.d/lattice.sh\n"
               "alias ls='ls --color=auto'\n"
               "alias ll='ls -la'\n"
               "alias lattice='cat /lattice/state.env'\n"
               "alias slots='ls /lattice/slots/ | head -32'\n"
               "export HISTFILE=/home/slot4096/.bash_history\n"
               "EOBASHRC\n"
               "cp /home/slot4096/.bashrc /home/slot4096/.bash_profile\n"
               "chown slot4096:slot4096 /home/slot4096/.bashrc /home/slot4096/.bash_profile\n"
               "echo '[lattice] ~/.bashrc written for slot4096'\n\n");

    /* 10. /etc/issue — lattice fingerprint shown before login prompt */
    fprintf(f, "# -- 10. /etc/issue --\n"
               "printf '\\nSlot4096 Lattice-Powered Alpine Linux\\n"
               "  host: %s  seed: 0x%016llx\\n"
               "  slots: %d  steps: %d  nice: %+d  uid: %d\\n\\n' "
               "> /etc/issue\n"
               "echo '[lattice] /etc/issue written'\n\n",
               hostname, (unsigned long long)seed,
               lattice_N, lattice_seed_steps_done, nice_val, uid_val);

    /* 11. /etc/os-release — lattice-branded OS identity */
    fprintf(f, "# -- 11. /etc/os-release --\n"
               "cat > /etc/os-release << 'EOOSREL'\n"
               "NAME=\"Slot4096 Alpine Linux\"\n"
               "ID=alpine\n"
               "PRETTY_NAME=\"Slot4096 Lattice-Powered Alpine Linux\"\n"
               "LATTICE_SEED=0x%016llx\n"
               "LATTICE_HOST=%s\n"
               "LATTICE_SLOTS=%d\n"
               "LATTICE_STEPS=%d\n"
               "EOOSREL\n"
               "echo '[lattice] /etc/os-release written'\n\n",
               (unsigned long long)seed, hostname,
               lattice_N, lattice_seed_steps_done);

    /* 12. /lattice/slots/ — every slot value as a file */
    fputs("# -- 12. /lattice/slots/ : slot values as files --\n"
          "mkdir -p /lattice/slots\n", f);
    int fp_n = lattice_N < 256 ? lattice_N : 256;
    for (int i = 0; i < fp_n; ++i)
        fprintf(f, "echo '%.16f' > /lattice/slots/%d\n", lattice[i], i);
    fputs("chown -R slot4096:slot4096 /lattice/slots\n"
          "echo '[lattice] /lattice/slots/ written.'\n\n", f);

    /* 13. MOTD -- 64-char wide box */
    fprintf(f, "# -- 13. MOTD --\n"
               "cat > /etc/motd << 'EOMOTD'\n"
               "+--------------------------------------------------------------+\n"
               "| Slot4096 Lattice-Powered Alpine Linux                        |\n"
               "| seed : 0x%016llx                               |\n"
               "| slots: %-5d  steps: %-5d  nice: %-+4d                       |\n"
               "| host : %-30s                   |\n"
               "| tz   : %-20s                             |\n"
               "| user : slot4096  uid: %-5d  umask: %04o                     |\n"
               "| lim  : nofile=%-6d  hist=%-5d  tmout=%-4ds               |\n"
               "| pkgs : %-52s |\n"
               "+--------------------------------------------------------------+\n"
               "EOMOTD\n\n",
               (unsigned long long)seed,
               lattice_N, lattice_seed_steps_done, nice_val,
               hostname, tz,
               uid_val, umask_val,
               nofile_val, histsize_val, tmout_val,
               pkgs);

    /* 14. Activate and drop into slot4096 login shell */
    fputs("# -- 14. Launch: su into slot4096 lattice user --\n"
          "cat /etc/motd\n"
          "echo '[lattice] OS fully initialised -- resonance substrate active'\n"
          "exec su - slot4096 -s /bin/bash 2>/dev/null || exec su - slot4096\n", f);

    fclose(f);
}

static void module_alpine_os(void) {
    printf("\n" CYAN
        "+--------------------------------------------------------------+\n"
        "|  Alpine OS Shell  --  Slot4096 Lattice-Powered OS           |\n"
        "+--------------------------------------------------------------+\n"
        CR "\n");

    if (!lattice_alpine_installed) {
        printf("  " YEL "[warn] Lattice not yet seeded.\n"
               "         Run Alpine Install [6] first.\n"
               "         Proceeding with current (blank) lattice state.\n" CR "\n");
    } else {
        printf("  " GRN "[OK] Lattice live: %d slots, %d seed steps.\n" CR "\n",
               lattice_N, lattice_seed_steps_done);
    }

    /* Derive all OS parameters and display them */
    uint64_t seed = lattice_derive_seed();
    char hostname[64]; lattice_derive_hostname(hostname, sizeof(hostname));
    const char *mirror = lattice_derive_mirror();
    char pkgs[256];    lattice_derive_packages(pkgs, sizeof(pkgs));
    char tz[32];       lattice_derive_tz(tz, sizeof(tz));
    int  nice_val      = lattice_derive_nice();
    int  uid_val       = lattice_derive_uid();
    int  nofile_val    = lattice_derive_nofile();
    int  umask_val     = lattice_derive_umask();
    int  histsize_val  = lattice_derive_histsize();
    int  tmout_val     = lattice_derive_tmout();

    printf("  " BOLD "Lattice-derived OS configuration:" CR "\n");
    printf("    Hostname  : " GRN "%s" CR "\n",     hostname);
    printf("    APK mirror: " GRN "%s" CR "\n",     mirror);
    printf("    Packages  : " GRN "%s" CR "\n",     pkgs);
    printf("    Timezone  : " GRN "%s" CR "\n",     tz);
    printf("    Nice      : " GRN "%+d" CR "\n",    nice_val);
    printf("    User      : " GRN "slot4096  uid=%d" CR "\n",  uid_val);
    printf("    umask     : " GRN "%04o" CR "\n",   umask_val);
    printf("    ulimit -n : " GRN "%d" CR "\n",     nofile_val);
    printf("    HISTSIZE  : " GRN "%d" CR "\n",     histsize_val);
    printf("    TMOUT     : " GRN "%ds" CR "\n",    tmout_val);
    printf("    Entropy   : " GRN "%d slots × 8 bytes = %d bytes" CR "\n",
           lattice_N, lattice_N * 8);
    printf("    Seed      : " GRN "0x%016llx" CR "\n\n", (unsigned long long)seed);

    /* Write all artefacts */
    char cwd[MAX_PATH];
    GetCurrentDirectoryA(MAX_PATH, cwd);

    char env_path[MAX_PATH], ent_path[MAX_PATH], sh_path[MAX_PATH];
    snprintf(env_path, sizeof(env_path), "%s\\lattice_state.env", cwd);
    snprintf(ent_path, sizeof(ent_path), "%s\\lattice_entropy.bin", cwd);
    snprintf(sh_path,  sizeof(sh_path),  "%s\\lattice_init.sh",    cwd);

    lattice_write_state_env(env_path);
    lattice_write_entropy_blob(ent_path);

    /* Build Linux-style paths for use inside Docker/WSL */
    char ent_linux[MAX_PATH], sh_linux[MAX_PATH];
    snprintf(ent_linux, sizeof(ent_linux), "%s", ent_path);
    snprintf(sh_linux,  sizeof(sh_linux),  "%s", sh_path);
    for (char *p = ent_linux; *p; p++) if (*p == '\\') *p = '/';
    for (char *p = sh_linux;  *p; p++) if (*p == '\\') *p = '/';

    /* Write init script (uses linux paths for internal references) */
    char ent_docker[MAX_PATH];
    snprintf(ent_docker, sizeof(ent_docker), "/lattice/entropy.bin");
    lattice_write_init_sh(sh_path, ent_docker);

    printf("  " DIM "lattice_state.env    -> %s\n" CR, env_path);
    printf("  " DIM "lattice_entropy.bin  -> %d bytes raw slot data\n" CR, lattice_N * 8);
    printf("  " DIM "lattice_init.sh      -> full OS boot script\n" CR "\n");

    /* ── Probe Docker ─────────────────────────────────────────────────── */
    int have_docker = 0, have_docker_desktop = 0;
    {
        const char *dd = "C:\\Program Files\\Docker\\Docker\\Docker Desktop.exe";
        if (GetFileAttributesA(dd) != INVALID_FILE_ATTRIBUTES)
            have_docker_desktop = 1;
        FILE *tp = _popen("docker info >NUL 2>&1 && echo YES", "r");
        if (tp) {
            char tmp[8] = {0};
            if (fgets(tmp, sizeof(tmp), tp) && tmp[0] == 'Y') have_docker = 1;
            _pclose(tp);
        }
    }

    /* ── Probe WSL distros ────────────────────────────────────────────── */
    char wsl_distro[128] = {0};
    int  wsl_is_alpine   = 0;
    {
        FILE *wp = _popen("wsl --list --quiet 2>NUL", "r");
        if (wp) {
            char line[256];
            while (fgets(line, sizeof(line), wp)) {
                char clean[256]; int ci = 0;
                for (int i = 0; line[i] && ci < 255; i++)
                    if (line[i] != '\0' && line[i] != '\r' && line[i] != '\n')
                        clean[ci++] = line[i];
                clean[ci] = '\0';
                if (!clean[0]) continue;
                char lower[256];
                for (int i = 0; clean[i]; i++)
                    lower[i] = (char)tolower((unsigned char)clean[i]);
                lower[ci] = '\0';
                if (strstr(lower, "docker-desktop")) continue;
                if (strstr(lower, "alpine")) {
                    strncpy(wsl_distro, clean, sizeof(wsl_distro)-1);
                    wsl_is_alpine = 1;
                    break;
                }
                if (!wsl_distro[0])
                    strncpy(wsl_distro, clean, sizeof(wsl_distro)-1);
            }
            _pclose(wp);
        }
    }
    int have_wsl = (wsl_distro[0] != '\0');

    /* ── Status display ───────────────────────────────────────────────── */
    if (have_docker)
        printf("  Docker     : " GRN "running" CR "\n");
    else if (have_docker_desktop)
        printf("  Docker     : " YEL "installed (not running)" CR "\n");
    else
        printf("  Docker     : " YEL "not found" CR "\n");

    if (have_wsl)
        printf("  WSL        : " GRN "%s%s" CR "\n", wsl_distro,
               wsl_is_alpine ? " (Alpine)" : "");
    else
        printf("  WSL        : " YEL "no distros found" CR "\n");
    printf("\n");

    /* ── Wake Docker Desktop if sleeping ─────────────────────────────── */
    if (!have_docker && have_docker_desktop) {
        printf("  " YEL "Attempting to start Docker Desktop...\n" CR);
        fflush(stdout);
        system("start \"\" \"C:\\Program Files\\Docker\\Docker\\Docker Desktop.exe\"");
        for (int w = 0; w < 15 && !have_docker; w++) {
            Sleep(1000);
            printf("  waiting for Docker daemon (%d/15)...\r", w+1); fflush(stdout);
            FILE *tp2 = _popen("docker info >NUL 2>&1 && echo YES", "r");
            if (tp2) {
                char tmp[8] = {0};
                if (fgets(tmp, sizeof(tmp), tp2) && tmp[0] == 'Y') have_docker = 1;
                _pclose(tp2);
            }
        }
        printf("\n");
        if (have_docker) printf("  " GRN "[OK] Docker daemon up.\n" CR "\n");
        else             printf("  " YEL "[warn] Docker still down -- using WSL.\n" CR "\n");
    }

    /* ── Spawn ────────────────────────────────────────────────────────── */
    if (have_docker) {
        /* Check if container is already running */
        int container_running = 0, container_exists = 0;
        {
            FILE *cp = _popen("docker inspect --format={{.State.Status}} phi4096-lattice 2>NUL", "r");
            if (cp) {
                char st[32] = {0};
                if (fgets(st, sizeof(st), cp)) {
                    if (strstr(st, "running"))  container_running = 1;
                    if (st[0] && st[0] != '\n') container_exists  = 1;
                }
                _pclose(cp);
            }
        }

        if (container_running) {
            /* ── Re-attach to running container as slot4096 ─────────── */
            printf("  " GRN "[live] Container phi4096-lattice is running.\n"
                   "  " YEL "Attaching as slot4096...\n"
                   "  " DIM "(Detach: Ctrl-P Ctrl-Q   |   New shell: docker exec -it phi4096-lattice bash)\n"
                   CR "\n");
            fflush(stdout);
            system("docker exec -it phi4096-lattice su - slot4096 -s /bin/bash 2>/dev/null "
                   "|| docker exec -it phi4096-lattice bash");

        } else {
            /* ── Fresh boot ─────────────────────────────────────────── */
            if (container_exists)
                system("docker rm -f phi4096-lattice >nul 2>nul");

            char docker_sh[MAX_PATH], docker_ent[MAX_PATH], docker_env[MAX_PATH];
            snprintf(docker_sh,  sizeof(docker_sh),  "%s", sh_path);
            snprintf(docker_ent, sizeof(docker_ent), "%s", ent_path);
            snprintf(docker_env, sizeof(docker_env), "%s", env_path);
            for (char *p = docker_sh;  *p; p++) if (*p == '\\') *p = '/';
            for (char *p = docker_ent; *p; p++) if (*p == '\\') *p = '/';
            for (char *p = docker_env; *p; p++) if (*p == '\\') *p = '/';

            char dcmd[2048];
            snprintf(dcmd, sizeof(dcmd),
                "docker run -it --name phi4096-lattice "
                "--cap-add SYS_ADMIN "
                "-e LATTICE_N=%d -e LATTICE_STEPS=%d "
                "-e LATTICE_SEED=0x%016llx -e LATTICE_INSTALLED=%d "
                "-v \"%s:/lattice/init.sh:ro\" "
                "-v \"%s:/lattice/entropy.bin:ro\" "
                "-v \"%s:/lattice/state.env:ro\" "
                "-w /lattice alpine sh /lattice/init.sh",
                lattice_N, lattice_seed_steps_done,
                (unsigned long long)seed, lattice_alpine_installed,
                docker_sh, docker_ent, docker_env);

            printf("  " CYAN "Launching Docker Alpine (full lattice substrate)...\n"
                   "  " YEL  "Container name: phi4096-lattice\n"
                   "  " YEL  "Second shell:   docker exec -it phi4096-lattice bash\n"
                   "  " DIM  "Detach without stopping: Ctrl-P Ctrl-Q\n"
                   CR "\n");
            fflush(stdout);
            system(dcmd);
            /* Clean up stopped container so next launch boots fresh */
            system("docker rm phi4096-lattice >nul 2>nul");
        }

    } else if (have_wsl) {
        /*
         * Convert Windows path to WSL path via wslpath, then run init.sh.
         * The entropy blob and state env are readable via the WSL mount at
         * /mnt/c/... automatically.
         */
        char wsl_sh_cmd[1024];
        /* Get the WSL-translated path for the init script */
        char wpath_cmd[512];
        snprintf(wpath_cmd, sizeof(wpath_cmd),
            "wsl --distribution %s -- wslpath -u '%s'",
            wsl_distro, sh_path);
        char wsl_sh_path[512] = {0};
        FILE *wp2 = _popen(wpath_cmd, "r");
        if (wp2) {
            if (fgets(wsl_sh_path, sizeof(wsl_sh_path), wp2)) {
                /* trim trailing newline */
                int l = (int)strlen(wsl_sh_path);
                while (l > 0 && (wsl_sh_path[l-1]=='\n'||wsl_sh_path[l-1]=='\r'))
                    wsl_sh_path[--l] = '\0';
            }
            _pclose(wp2);
        }
        /* Also get wslpath for entropy bin */
        char wpath_ent[512];
        snprintf(wpath_ent, sizeof(wpath_ent),
            "wsl --distribution %s -- wslpath -u '%s'",
            wsl_distro, ent_path);
        char wsl_ent_path[512] = {0};
        FILE *wp3 = _popen(wpath_ent, "r");
        if (wp3) {
            if (fgets(wsl_ent_path, sizeof(wsl_ent_path), wp3)) {
                int l = (int)strlen(wsl_ent_path);
                while (l > 0 && (wsl_ent_path[l-1]=='\n'||wsl_ent_path[l-1]=='\r'))
                    wsl_ent_path[--l] = '\0';
            }
            _pclose(wp3);
        }

        if (wsl_sh_path[0]) {
            snprintf(wsl_sh_cmd, sizeof(wsl_sh_cmd),
                "wsl --distribution %s -- "
                "env LATTICE_N=%d LATTICE_STEPS=%d "
                "LATTICE_SEED=0x%016llx LATTICE_INSTALLED=%d "
                "LATTICE_ENT='%s' "
                "sh '%s'",
                wsl_distro,
                lattice_N, lattice_seed_steps_done,
                (unsigned long long)seed, lattice_alpine_installed,
                wsl_ent_path[0] ? wsl_ent_path : "/dev/null",
                wsl_sh_path);
            printf("  " CYAN "Launching WSL %s (full lattice substrate)...\n" CR "\n",
                   wsl_distro);
        } else {
            /* fallback: inline export */
            snprintf(wsl_sh_cmd, sizeof(wsl_sh_cmd),
                "wsl --distribution %s -- sh -c \""
                "export LATTICE_N=%d LATTICE_STEPS=%d "
                "LATTICE_SEED=0x%016llx LATTICE_INSTALLED=%d; "
                "echo '=== Slot4096 Lattice-Powered OS (%s) ==='; "
                "echo Seed: $LATTICE_SEED; sh\"",
                wsl_distro,
                lattice_N, lattice_seed_steps_done,
                (unsigned long long)seed, lattice_alpine_installed,
                wsl_distro);
            printf("  " YEL "Launching WSL %s (env-only fallback)...\n" CR "\n", wsl_distro);
        }
        fflush(stdout);
        system(wsl_sh_cmd);

    } else {
        printf("  " RED "No Linux runtime available.\n" CR "\n");
        printf("  Options:\n");
        printf("    1. Start Docker Desktop, then re-run [8]\n");
        printf("    2. wsl --install -d Alpine\n");
        printf("    3. wsl --install -d Ubuntu\n\n");
        printf("  " YEL "Generated files:\n" CR);
        printf("    %s\n    %s\n    %s\n", env_path, ent_path, sh_path);
        printf("  Source in any shell:  . lattice_state.env\n");
    }
}

/* ══════════════════════════ MAIN ════════════════════════════════════════════ */

static void print_banner(void) {
    printf(CYAN
        "+============================================================+\n"
        "|   Quantum Prime Library   v1.0  (Windows TUI)             |\n"
        "|   phi-lattice  *  Dn(r)  *  psi-score  *  Miller-Rabin    |\n"
        "+============================================================+\n"
        CR "\n");
}

static void print_menu(void) {
    printf(WHT "  Main Menu\n" CR);
    printf("  " YEL "[1]" CR " Prime Pipeline       sieve -> phi-filter -> Dn-rank\n");
    printf("  " YEL "[2]" CR " Number Analyzer      Miller-Rabin, phi-lattice, Dn, psi\n");
    printf("  " YEL "[3]" CR " Mersenne Explorer    M1..M51 + next candidate predictions\n");
    printf("  " YEL "[4]" CR " Zeta Zeros           zeta(1/2+it) table + Gram approximation\n");
    printf("  " YEL "[5]" CR " Benchmark            time all 13 prime functions\n");
    printf("  " YEL "[6]" CR " Alpine Install       boot lattice + GPU resonance hook\n");
    printf("  " YEL "[7]" CR " Lattice Shell        interactive Slot4096 REPL\n");
    printf("  " YEL "[8]" CR " Alpine OS Shell      spawn lattice-powered Alpine Linux\n");
    printf("  " YEL "[Q]" CR " Quit\n\n");
    printf("  " BOLD ">" CR " "); fflush(stdout);
}

int main(void) {
    console_init();
    print_banner();

    for (;;) {
        print_menu();
        DWORD old; GetConsoleMode(g_hin, &old);
        SetConsoleMode(g_hin, ENABLE_PROCESSED_INPUT|ENABLE_ECHO_INPUT|ENABLE_LINE_INPUT);
        WCHAR wb[8] = {0}; DWORD nr = 0;
        ReadConsoleW(g_hin, wb, 7, &nr, NULL);
        SetConsoleMode(g_hin, old);

        char ch = (wb[0] > 0 && wb[0] <= 127) ? (char)wb[0] : '?';
        char lc = (ch >= 'A' && ch <= 'Z') ? ch + 32 : ch;

        switch (lc) {
        case '1': module_pipeline();  break;
        case '2': module_analyzer();  break;
        case '3': module_mersenne();  break;
        case '4': module_zeta();      break;
        case '5': module_benchmark(); break;
        case '6': module_alpine();    break;
        case '7': module_lattice_shell(); break;
        case '8': module_alpine_os();     break;
        case 'q':
            printf("\n" CYAN "  Goodbye.\n" CR "\n");
            return 0;
        default:
            printf("\n  " RED "Unknown option '%c'\n" CR, ch);
            break;
        }
        wait_enter();
        printf("\n");
    }
}
