#include <stdint.h>
#define crypto_int8 int8_t

/* Supplied here: `scalarmult.c` calls this helper, which upstream generates into
   its own `crypto_int8.h` and which is therefore not part of the `ref10`
   directory.  It returns an all-ones mask when bit `pos mod 8` of `x` is set and
   zero otherwise; `scalarmult.c` negates it to recover the bit, which is what
   `fe_cswap` wants.  The shift amount is the loop counter and so is public. */
#define crypto_int8_bitmod_mask(x, pos) \
  ((crypto_int8) - (crypto_int8)(((uint8_t)(x) >> ((pos) & 7)) & 1))
