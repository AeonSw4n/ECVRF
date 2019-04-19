#include "../lib/elligator2.h"

static void elligator_fe_A(fe out)
{
    out[0] = 486662;
    out[1] = 0;
    out[2] = 0;
    out[3] = 0;
    out[4] = 0;
    out[5] = 0;
    out[6] = 0;
    out[7] = 0;
    out[8] = 0;
    out[9] = 0;

}

static void fe_sqAp2(fe out)
{
  out[0] = 8930344;
  out[1] = 9583591;
  out[2] = -26444492;
  out[3] = 3752533;
  out[4] = 26044487;
  out[5] = -743697;
  out[6] = -2900628;
  out[7] = 5634116;
  out[8] = 25139868;
  out[9] = -5270574;
}

static void fe_u252m2(fe out)
{
  out[0] = -32595791;
  out[1] = -7943725;
  out[2] = 9377950;
  out[3] = 3500415;
  out[4] = 12389472;
  out[5] = -272473;
  out[6] = -25146209;
  out[7] = -2005654;
  out[8] = 326686;
  out[9] = 11406482;
}

static void fe_sqr_legendre(fe sqr, fe e, const fe w){
  fe ww;
  fe l0;
  fe l1;
  fe l2;
  fe l3;
  int i;

  /*
   * Compute w ** ((p-1)/2) = w ** (2 ** 254 - 10) with the exponent as
   * 2 ** 254 - 10 = (2 ** 4) * (2 ** 250 - 1) + 6.
   */

  /* sw = w ** 2 */
  fe_sq(ww, w);

  /* l0 = w ** 2 */
  fe_copy(l0, ww);

  /* l3 = l0 * w = w ** 3 */
  fe_mul(l3, l0, w);

  /* l0 = l3 ** 2 = w ** 6 -- stash l0 away for the end. */
  fe_sq(l0, l3);

  /* l1 = l0 * w = w ** 7 */
  fe_mul(l1, l0, w);

  /* l2 = l1 ** 2 = w ** 14 */
  fe_sq(l2, l1);

  /* l2 = l2 ** 2 = w ** 28 */
  fe_sq(l2, l2);

  /* l1 = l3 * l2 = w ** 31 =  w ** (2 ** 5 - 1) */
  fe_mul(l1, l3, l2);

  /* l2 = l1 ** (2 ** 5) = w ** ((2 ** 5) * (2 ** 5 - 1)) */
  fe_sq(l2, l1);
  for (i = 1; i < 5; ++i) {
      fe_sq(l2, l2);
  }

  /* l1 = l1 * l2 = w ** ((2 ** 5 + 1) * (2 ** 5 - 1)) = w ** (2 ** 10 - 1) */
  fe_mul(l1, l2, l1);

  /* Continuing similarly... */

  /* l2 = w ** (2 ** 20 - 1) */
  fe_sq(l2, l1);
  for (i = 1; i < 10; ++i) {
      fe_sq(l2, l2);
  }
  fe_mul(l2, l2, l1);

  /* l2 = w ** (2 ** 40 - 1) */
  fe_sq(l3, l2);
  for (i = 1; i < 20; ++i) {
      fe_sq(l3, l3);
  }
  fe_mul(l2, l3, l2);

  /* l2 = w ** (2 ** 10) * (2 ** 40 - 1) */
  for (i = 0; i < 10; ++i) {
      fe_sq(l2, l2);
  }
  /* l1 = w ** (2 ** 50 - 1) */
  fe_mul(l1, l2, l1);

  /* l2 = w ** (2 ** 100 - 1) */
  fe_sq(l2, l1);
  for (i = 1; i < 50; ++i) {
      fe_sq(l2, l2);
  }
  fe_mul(l2, l2, l1);

  /* l2 = w ** (2 ** 200 - 1) */
  fe_sq(l3, l2);
  for (i = 1; i < 100; ++i) {
      fe_sq(l3, l3);
  }
  fe_mul(l2, l3, l2);

  /* l2 = w ** ((2 ** 50) * (2 ** 200 - 1) */
  fe_sq(l2, l2);
  for (i = 1; i < 50; ++i) {
      fe_sq(l2, l2);
  }

  /* l1 = w ** (2 ** 250 - 1) */
  fe_mul(l1, l2, l1);

  /* l1 = w ** (2 ** 251 - 2) */
  fe_sq(l1, l1);

  /* l1 = w ** (2 ** 252 - 4) */
  fe_sq(l1, l1);

  /* sqrt = w ** (2 ** 252 - 2) */
  fe_mul(sqr, l1, ww);

  /* li = w ** (2 ** 254 - 16) */
  fe_sq(l1, l1);
  fe_sq(l1, l1);

  /* Recall l0 = w ** 6; out = w ** (2 ** 254 - 10) */
  fe_mul(e, l1, l0);
}

static void fe_legendre(fe out, const fe z)
{
    fe t0;
    fe t1;
    fe t2;
    fe t3;
    int i;

    /*
     * Compute z ** ((p-1)/2) = z ** (2 ** 254 - 10) with the exponent as
     * 2 ** 254 - 10 = (2 ** 4) * (2 ** 250 - 1) + 6.
     */

    /* t0 = z ** 2 */
    fe_sq(t0, z);
    //return;
    /* t3 = t0 * z = z ** 3 */
    fe_mul(t3, t0, z);

    /* t0 = t3 ** 2 = z ** 6 -- stash t0 away for the end. */
    fe_sq(t0, t3);

    /* t1 = t0 * z = z ** 7 */
    fe_mul(t1, t0, z);

    /* t2 = t1 ** 2 = z ** 14 */
    fe_sq(t2, t1);

    /* t2 = t2 ** 2 = z ** 28 */
    fe_sq(t2, t2);

    /* t1 = t3 * t2 = z ** 31 =  z ** (2 ** 5 - 1) */
    fe_mul(t1, t3, t2);

    /* t2 = t1 ** (2 ** 5) = z ** ((2 ** 5) * (2 ** 5 - 1)) */
    fe_sq(t2, t1);
    for (i = 1; i < 5; ++i) {
        fe_sq(t2, t2);
    }

    /* t1 = t1 * t2 = z ** ((2 ** 5 + 1) * (2 ** 5 - 1)) = z ** (2 ** 10 - 1) */
    fe_mul(t1, t2, t1);

    /* Continuing similarly... */

    /* t2 = z ** (2 ** 20 - 1) */
    fe_sq(t2, t1);
    for (i = 1; i < 10; ++i) {
        fe_sq(t2, t2);
    }
    fe_mul(t2, t2, t1);

    /* t2 = z ** (2 ** 40 - 1) */
    fe_sq(t3, t2);
    for (i = 1; i < 20; ++i) {
        fe_sq(t3, t3);
    }
    fe_mul(t2, t3, t2);

    /* t2 = z ** (2 ** 10) * (2 ** 40 - 1) */
    for (i = 0; i < 10; ++i) {
        fe_sq(t2, t2);
    }
    /* t1 = z ** (2 ** 50 - 1) */
    fe_mul(t1, t2, t1);

    /* t2 = z ** (2 ** 100 - 1) */
    fe_sq(t2, t1);
    for (i = 1; i < 50; ++i) {
        fe_sq(t2, t2);
    }
    fe_mul(t2, t2, t1);

    /* t2 = z ** (2 ** 200 - 1) */
    fe_sq(t3, t2);
    for (i = 1; i < 100; ++i) {
        fe_sq(t3, t3);
    }
    fe_mul(t2, t3, t2);

    /* t2 = z ** ((2 ** 50) * (2 ** 200 - 1) */
    fe_sq(t2, t2);
    for (i = 1; i < 50; ++i) {
        fe_sq(t2, t2);
    }

    /* t1 = z ** (2 ** 250 - 1) */
    fe_mul(t1, t2, t1);

    /* t1 = z ** ((2 ** 4) * (2 ** 250 - 1)) */
    fe_sq(t1, t1);
    for (i = 1; i < 4; ++i) {
        fe_sq(t1, t1);
    }

    /* Recall t0 = z ** 6; out = z ** (2 ** 254 - 10) */
    fe_mul(out, t1, t0);
}

static void fe_legendre_fromby(fe a, by b)
{
  fe_frombytes(a, b);
  fe_legendre(a, a);
}

static uint32_t fe_legendre_ifromby(by b)
{
  fe fen, lfen;
  fe_frombytes(fen, b);
  fe_legendre(lfen, fen);
  return lfen[0];
}

static void elligator2_ed25519(const uint8_t *data, size_t size,
                       const uint8_t public_key[32],
                       uint8_t out_point[32])
{
    static const uint8_t SUITE[1] = {0x01};
    static const uint8_t ONE[1] = {0x01};

    // 3. hash(suite || one || pk || alpha)

    uint8_t hash[SHA512_DIGEST_LENGTH] = {0};
    SHA512_CTX hash_ctx;
    SHA512_Init(&hash_ctx);
    SHA512_Update(&hash_ctx, SUITE, sizeof(SUITE));
    SHA512_Update(&hash_ctx, ONE, sizeof(ONE));
    SHA512_Update(&hash_ctx, public_key, 32);
    SHA512_Update(&hash_ctx, data, size);
    SHA512_Final(hash, &hash_ctx);

    // 4. take first 32 bytes of the hash

    uint8_t truncatedHash[32] = {0};
    memcpy(truncatedHash, hash, 32);

    // 7. take highest order bit of truncated hash
    // 8. clear the bit in the source

    uint8_t x0 = truncatedHash[31] & 0x80;
    truncatedHash[31] &= 0x7f;

    // 9. convert to integer

    fe r = {0};
    fe_frombytes(r, truncatedHash);

    // 10. u = - A / (1 + 2*(r^2) ) mod p

    fe t0 = {0};
    fe t1 = {0};

    fe one = {0};
    fe_1(one);

    fe A = {0};
    elligator_fe_A(A);


    fe_sq(t0, r);       // r ** 2
    fe_add(t0, t0, t0);   // 2 * (r ** 2)
    fe_add(t0, one, t0);  // 1 + 2 * (r ** 2)
    fe_invert(t0, t0);    // 1 / (1 + 2 * (r ** 2))

    fe u = {0};
    fe_0(u);
    fe_sub(u, u, A);      // -A
    fe_mul(u, u, t0);     // -A / (1 + 2 * (r ** 2))

    // 11. w = u * (u^2 + A*u + 1) mod p

    fe_sq(t0, u);       // u ** 2
    fe_mul(t1, A, u);    // A * u
    fe_add(t0, t0, t1);  // u**2 + A*u
    fe_add(t0, t0, one); // u**2 + A*u + 1

    fe w = {0};
    fe_mul(w, u, t0);    // u * (u**2 + A*u + 1)

    // 12. e = Legendre symbol of w and p

    fe e = {0};
    fe_legendre(e, w);   // w ** ((p-1)/2)

    fe u2 = {0};
    fe_0(u2);
    fe_sub(u2, u2, A);   // -A
    fe_sub(u2, u2, u);   // -A - u

    e[0] = (e[0] + 1) & 2;
    unsigned int b = e[0] >> 1;
    fe_swap(u, u2, b);  // swaps if b == 1

    fe uf = {0};
    fe_copy(uf, u2);

    // 14. y coordinate

    fe_sub(t0, uf, one); // t0 = uf - 1

    fe_add(t1, uf, one); // t1 = uf + 1
    fe_invert(t1, t1);   // t1 = 1 / (uf + 1)

    fe y = {0};
    fe_mul(y, t0, t1);   // y = (uf - 1) / (uf + 1)


      /*
      by h;
      fe_tobytes(h, y);
      h[31] |= x0;
      */

    // 15. encode point

    ge_p1p1 hc;
    ge_p2 hp;

    fe_copy(hp.X, uf);
    fe_copy(hp.Y, y);
    fe_1(hp.Z);


    // 16. out_point = (uf, y) ^ 8
    ge_p2_dbl(&hc, &hp);
    ge_p1p1_to_p2(&hp, &hc);
    ge_p2_dbl(&hc, &hp);
    ge_p1p1_to_p2(&hp, &hc);
    ge_p2_dbl(&hc, &hp);
    ge_p1p1_to_p2(&hp, &hc);

    ge_tobytes(out_point, &hp);

/*
    uint8_t cofactor[32] = {0};
    cofactor[0] = 8;
    x25519_scalar_mult(h, cofactor, out_point);
*/
}

static void elligator2_cofactor_mult(by out, by in)
{

}
