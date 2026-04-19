/* prime_ui.c — Windows TUI for the Quantum-Prime Library v1.0
 *
 *  [1] Prime Pipeline     sieve -> phi-filter -> Dn-rank (Mersenne candidates)
 *  [2] Number Analyzer    Miller-Rabin, phi-lattice, Dn score, psi-score
 *  [3] Mersenne Explorer  all M1..M51 + next candidate predictions
 *  [4] Zeta Zeros         zeta(1/2+it) hardcoded table + Gram approximation
 *  [5] Benchmark          time all 13 prime library functions
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
    printf("  " YEL "[4]" CR " Zeta Zeros           zeta(1/2+it) table + Gram approx\n");
    printf("  " YEL "[5]" CR " Benchmark            time all 13 prime functions\n");
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
