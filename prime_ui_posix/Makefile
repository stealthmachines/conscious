# Makefile — POSIX/Linux build for prime_ui (Quantum-Prime Library v1.0)
#
# Requires: clang (or gcc) with AVX2/FMA/SSE4.1/AES-NI/PCLMUL/RDSEED/RDRAND
# Tested on:  Ubuntu 22.04 / Debian 12, x86-64, clang 14+ or gcc 11+
#
# Usage:
#   make              build prime_ui binary
#   make CC=gcc       build with GCC instead of Clang
#   make clean        remove build artefacts
#   make install      install to /usr/local/bin
#   make ASAN=1       build with AddressSanitizer + LeakSanitizer
#
# Runtime requirements:
#   - Linux 3.17+ (getrandom syscall for BCryptGenRandom replacement)
#   - Docker (optional, for Alpine OS shell modules)

CC      ?= clang

# ── CPU feature flags ────────────────────────────────────────────────────────
# All of these are required for the math/crypto inlines:
#   -mavx2       AVX2 256-bit SIMD (phi-lattice, KDF)
#   -mfma        Fused multiply-add (phi_fold_hash, phi_stream)
#   -msse4.1     SSE 4.1 (_mm_floor_pd, _mm_blendv_pd)
#   -maes        AES-NI (_mm_aesenc_si128 etc.)
#   -mpclmul     Carry-less multiply (ghash_update in full-crypto module)
#   -mrdseed     RDSEED instruction (phi_hw_rng64)
#   -mrdrnd      RDRAND instruction (phi_hw_rng64 fallback)
CPU_FLAGS = -mavx2 -mfma -msse4.1 -maes -mpclmul -mrdseed -mrdrnd

# ── Common flags ─────────────────────────────────────────────────────────────
CFLAGS  = -O2 $(CPU_FLAGS) -D_GNU_SOURCE \
          -Wall -Wextra \
          -Wno-unused-parameter \
          -Wno-unused-function \
          -Wno-sign-compare \
          -Wno-implicit-fallthrough \
          -lm

# ── ASAN build (make ASAN=1) ──────────────────────────────────────────────────
ifdef ASAN
CFLAGS  += -fsanitize=address,undefined -fno-omit-frame-pointer -g
endif

TARGET  = prime_ui
SRC     = prime_ui_posix.c
HDR     = posix_compat.h

.PHONY: all clean install

all: $(TARGET)

$(TARGET): $(SRC) $(HDR)
	$(CC) $(CFLAGS) -o $@ $(SRC)

install: $(TARGET)
	install -m 755 $(TARGET) /usr/local/bin/$(TARGET)

clean:
	rm -f $(TARGET) *.o
