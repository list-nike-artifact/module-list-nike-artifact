/* X25519 from the SUPERCOP reference implementation `ref10`, on this machine.
   Checks itself against RFC 7748 before reporting anything. */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "crypto_scalarmult.h"

#define RUNS 1000
static double GHZ;

/* Same scheduling request the Rust harness makes: a performance core at full
   frequency.  Without it a short run can be measured on an efficiency core. */
extern int pthread_set_qos_class_self_np(unsigned int, int);

/* Busy-wait so the core is already at its top frequency when timing starts. */
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
  printf("%-16s %9.2f us %10.0f cycles\n", label, m / 1e3, m * GHZ);
}

static void unhex(unsigned char *out, const char *h) {
  for (int i = 0; i < 32; i++) sscanf(h + 2 * i, "%2hhx", &out[i]);
}

/* RFC 7748, Section 6.1. */
static int self_test(void) {
  unsigned char ask[32], apk[32], bsk[32], bpk[32], want[32], got[32], t[32];
  unhex(ask, "77076d0a7318a57d3c16c17251b26645df4c2f87ebc0992ab177fba51db92c2a");
  unhex(apk, "8520f0098930a754748b7ddcb43ef75a0dbf3a0d26381af4eba4a98eaa9b4e6a");
  unhex(bsk, "5dab087e624a8a4b79e17f8b83800ee66f3bb1292618b6fd1c2f8b27ff88e0eb");
  unhex(bpk, "de9edb7d7b7dc1b4d35b61c2ece435373f8343c85b78674dadfc7e146f882b4f");
  unhex(want, "4a5d9d5ba4ce2de1728e3bf480350f25e07e21c947d19e3376f09b3c1e161742");
  crypto_scalarmult_base(t, ask);
  if (memcmp(t, apk, 32)) return 1;
  crypto_scalarmult_base(t, bsk);
  if (memcmp(t, bpk, 32)) return 2;
  crypto_scalarmult(got, ask, bpk);
  if (memcmp(got, want, 32)) return 3;
  crypto_scalarmult(got, bsk, apk);
  if (memcmp(got, want, 32)) return 4;
  return 0;
}

int main(int argc, char **argv) {
  if (argc < 2) { fprintf(stderr, "usage: %s <GHz>\n", argv[0]); return 2; }
  GHZ = atof(argv[1]);
  int e = self_test();
  if (e) { fprintf(stderr, "RFC 7748 self-test failed at step %d\n", e); return 1; }
  printf("RFC 7748 test vectors: ok\n");

  pthread_set_qos_class_self_np(0x21, 0);
  warm_up();

  unsigned char sk[RUNS][32], pk[32], peer[32], shared[32];
  for (int i = 0; i < RUNS; i++)
    for (int j = 0; j < 32; j++) sk[i][j] = (unsigned char)(i * 31 + j * 7 + 1);
  crypto_scalarmult_base(peer, sk[0]);

  for (int i = 0; i < RUNS; i++) {
    double t = now_ns();
    crypto_scalarmult_base(pk, sk[i]);
    s[i] = now_ns() - t;
  }
  report("X25519 keygen");
  for (int i = 0; i < RUNS; i++) {
    double t = now_ns();
    crypto_scalarmult(shared, sk[i], peer);
    s[i] = now_ns() - t;
  }
  report("X25519 derive");
  printf("  pk %d bytes\n", 32);
  return 0;
}
