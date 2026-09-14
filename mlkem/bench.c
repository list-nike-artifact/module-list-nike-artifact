/* ML-KEM-768 from the pq-crystals reference C implementation, timed with the
   same harness as every other row: performance core, warm-up, median of 1000. */
#include "kem.h"
#include "params.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

#define RUNS 1000
static double GHZ;

extern int pthread_set_qos_class_self_np(unsigned int, int);
static void pin_to_p_core(void) { pthread_set_qos_class_self_np(0x21, 0); }

static void warm_up(void) {
  volatile unsigned long x = 0;
  struct timespec t0; clock_gettime(CLOCK_MONOTONIC_RAW, &t0);
  for (;;) {
    for (int i = 0; i < 100000; i++) x += i;
    struct timespec t1; clock_gettime(CLOCK_MONOTONIC_RAW, &t1);
    if ((t1.tv_sec - t0.tv_sec) * 1e9 + (t1.tv_nsec - t0.tv_nsec) > 3e8) break;
  }
}
static double now_ns(void) {
  struct timespec t; clock_gettime(CLOCK_MONOTONIC_RAW, &t);
  return t.tv_sec * 1e9 + t.tv_nsec;
}
static int cmp(const void *a, const void *b) {
  double x = *(const double *)a, y = *(const double *)b;
  return x < y ? -1 : x > y;
}
static double s[RUNS];
static void report(const char *label) {
  qsort(s, RUNS, sizeof(double), cmp);
  double m = s[RUNS / 2];
  printf("%-28s %8.2f us %9.0f cycles\n", label, m / 1e3, m * GHZ);
}

int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: %s <clock in GHz>\n", argv[0]); return 1; }
  GHZ = atof(argv[1]);
  pin_to_p_core();
  warm_up();

  static uint8_t pk[CRYPTO_PUBLICKEYBYTES], sk[CRYPTO_SECRETKEYBYTES];
  static uint8_t ct[CRYPTO_CIPHERTEXTBYTES], ss[CRYPTO_BYTES], ss2[CRYPTO_BYTES];

  for (int i = 0; i < RUNS; i++) {
    double t = now_ns();
    crypto_kem_keypair(pk, sk);
    s[i] = now_ns() - t;
  }
  report("ML-KEM-768 ref keygen");

  for (int i = 0; i < RUNS; i++) {
    double t = now_ns();
    crypto_kem_enc(ct, ss, pk);
    s[i] = now_ns() - t;
  }
  report("ML-KEM-768 ref encaps");

  for (int i = 0; i < RUNS; i++) {
    double t = now_ns();
    crypto_kem_dec(ss2, ct, sk);
    s[i] = now_ns() - t;
  }
  report("ML-KEM-768 ref decaps");

  if (memcmp(ss, ss2, CRYPTO_BYTES) != 0) { printf("SHARED SECRET MISMATCH\n"); return 1; }
  printf("shared secrets agree; |ct| %d  |pk| %d  |sk| %d\n",
         CRYPTO_CIPHERTEXTBYTES, CRYPTO_PUBLICKEYBYTES, CRYPTO_SECRETKEYBYTES);
  return 0;
}
