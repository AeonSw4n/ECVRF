#include "../lib/elligator2.h"

//#define DEBUG

/**
 * Returns, via parameter out, a field element equal to
 * A, where A = 48662 is the coefficient in the Mongtomery Curve 25519 equation v^2 = u(u^2 + Au + 1)
 */
static void fe_A(fe out)
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

/**
 * Returns, via parameter out, a field element equal to
 * A^3, where A = 48662 is the coefficient in the Mongtomery Curve 25519 equation v^2 = u(u^2 + Au + 1)
 */
static void fe_A3(fe out)
{
  out[0] = -8127272;
  out[1] = 6246418;
  out[2] = 51;
  out[3] = 0;
  out[4] = 0;
  out[5] = 0;
  out[6] = 0;
  out[7] = 0;
  out[8] = 0;
  out[9] = 0;
}

/**
 * Returns, via parameter out, a field element equal to the
 * square root of A + 2 (= 486664), i.e. (A + 2)^{2^{252}-2} mod p (where p = 2^{255}-19)
 * A = 48662 is the value from the Mongtomery Curve 25519 equation v^2 = u(u^2 + Au + 1)
 * // TODO: describe which square root -- positive or negative? And what is positive and what is negative?
 */
static void fe_sqrt_A_plus_2(fe out)
{
  out[0] = -8930344;
  out[1] = -9583591;
  out[2] = 26444492;
  out[3] = -3752533;
  out[4] = -26044487;
  out[5] = 743697;
  out[6] = 2900628;
  out[7] = -5634116;
  out[8] = -25139868;
  out[9] = 5270574;
}

/**
 * Returns, via parameter out, a field element equal to the
 * square root of 2, i.e. 2^{2^{252}-2} mod p (where p = 2^{255}-19)
 * // TODO: confirm
 */
static void fe_sqrt_2(fe out)
{
  out[0] = 32595791;
  out[1] = 7943725;
  out[2] = -9377950;
  out[3] = -3500415;
  out[4] = -12389472;
  out[5] = 272473;
  out[6] = 25146209;
  out[7] = 2005654;
  out[8] = -326686;
  out[9] = -11406482;
}

static unsigned int fe_parity(const fe h)
{
  int32_t h0 = h[0];
  int32_t h1 = h[1];
  int32_t h2 = h[2];
  int32_t h3 = h[3];
  int32_t h4 = h[4];
  int32_t h5 = h[5];
  int32_t h6 = h[6];
  int32_t h7 = h[7];
  int32_t h8 = h[8];
  int32_t h9 = h[9];
  int32_t q;

  q = (19 * h9 + (((int32_t) 1) << 24)) >> 25;
  q = (h0 + q) >> 26;
  q = (h1 + q) >> 25;
  q = (h2 + q) >> 26;
  q = (h3 + q) >> 25;
  q = (h4 + q) >> 26;
  q = (h5 + q) >> 25;
  q = (h6 + q) >> 26;
  q = (h7 + q) >> 25;
  q = (h8 + q) >> 26;
  q = (h9 + q) >> 25;

  /* Goal: Output h-(2^255-19)q, which is between 0 and 2^255-20. */
  h0 += 19 * q;
  h0 &= kBottom26Bits;
  return (unsigned int) ((uint8_t)(h0) & 1);
}

static unsigned int fe_ispositive(const fe h)
{
  int32_t h0 = h[0];
  int32_t h1 = h[1];
  int32_t h2 = h[2];
  int32_t h3 = h[3];
  int32_t h4 = h[4];
  int32_t h5 = h[5];
  int32_t h6 = h[6];
  int32_t h7 = h[7];
  int32_t h8 = h[8];
  int32_t h9 = h[9];
  int32_t q;

  q = (19 * h9 + (((int32_t) 1) << 24)) >> 25;
  q = (h0 + q) >> 26;
  q = (h1 + q) >> 25;
  q = (h2 + q) >> 26;
  q = (h3 + q) >> 25;
  q = (h4 + q) >> 26;
  q = (h5 + q) >> 25;
  q = (h6 + q) >> 26;
  q = (h7 + q) >> 25;
  q = (h8 + q) >> 26;
  q = (h9 + q) >> 25;

  h0 += 19 * q;
  h0 &= kBottom26Bits;
  return !(unsigned int) ((uint8_t)(h0) & 1);
}

static void fe_sqrt_and_legendre(fe sqrt, fe e, const fe w)
{
  fe w_squared;
  fe l0;
  fe l1;
  fe l2;
  fe l3;
  int i;

  /*
   * Compute the Legendre symbol e of w as
   * w ** ((p-1)/2) = w ** (2 ** 254 - 10) with the exponent as
   * 2 ** 254 - 10 = (2 ** 4) * (2 ** 250 - 1) + 6.
   * We will do so using 253 squarings and 11 multiplications.
   * In the process, also compute w ** (p+3)/8 = (2 ** 252 - 2),
   * which costs us one extra multiplication.
   */

  /* ww = w ** 2 */
  fe_sq(w_squared, w);

  /* l0 = w ** 2 */
  fe_copy(l0, w_squared);

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

  /* l2 = l1 ** (2 ** 5) = w ** ((2 ** 5) * (2 ** 5 - 1)) = w ** (2 ** 10 - 2 ** 5)*/
  fe_sq(l2, l1);
  for (i = 1; i < 5; ++i) {
      fe_sq(l2, l2);
  }

  /* l1 = l2 * l1 = w ** (2 ** 10 - 2 ** 5 + 2 ** 5 - 1) = w ** (2 ** 10 - 1) */
  fe_mul(l1, l2, l1);

  /* l2 = l1 ** (2 ** 10) = w ** ((2 ** 10) * (2 ** 10 - 1)) = w ** (2 ** 20 - 2 ** 10) */
  fe_sq(l2, l1);
  for (i = 1; i < 10; ++i) {
      fe_sq(l2, l2);
  }

  /* l2 = l2 * l1 = w ** (2 ** 20 - 2 ** 10 + 2 ** 10 - 1) = w ** (2 ** 20 - 1) */
  fe_mul(l2, l2, l1);

  /* l3 = l2 ** (2 ** 20) = w ** (2 ** 40 - 2 ** 20)) */
  fe_sq(l3, l2);
  for (i = 1; i < 20; ++i) {
      fe_sq(l3, l3);
  }
  /* l2 = l3 * l2 = w ** (2 ** 40 - 2 ** 20 + 2 ** 20 - 1) = w ** (2 ** 40 - 1) */
  fe_mul(l2, l3, l2);

  /* l2 = l2 ** (2 ** 10) = w ** ((2 ** 10) * (2 ** 40 - 1)) = w ** (2 ** 50 - 2 ** 10) */
  for (i = 0; i < 10; ++i) {
      fe_sq(l2, l2);
  }

  /* l1 = l2 * l1 = w ** (2 ** 50 - 2 ** 10 + 2 ** 10 - 1) = w ** (2 ** 50 - 1) */
  fe_mul(l1, l2, l1);

  /* l2 = l1 ** (2 ** 50) = w ** (2 ** 100 - 2 ** 50) */
  fe_sq(l2, l1);
  for (i = 1; i < 50; ++i) {
      fe_sq(l2, l2);
  }
  /* l2 = l2 ** l1 = w ** (2 ** 100 - 2 ** 50 + 2 ** 50 - 1)) = w ** (2 ** 100 - 1) */
  fe_mul(l2, l2, l1);

  /* l3 = l2 ** (2 ** 100) = w ** (2 ** 200 - 2 ** 100) */
  fe_sq(l3, l2);
  for (i = 1; i < 100; ++i) {
      fe_sq(l3, l3);
  }
  /* l2 = l3 * l2 = w ** (2 ** 200 - 2 ** 100 + 2 ** 100 - 1) = w ** (2 ** 200 - 1) */
  fe_mul(l2, l3, l2);

  /* l2 = l2 ** (2 ** 50) = w ** ((2 ** 50) * (2 ** 200 - 1) = w ** (2 ** 250 - 2 ** 50) */
  fe_sq(l2, l2);
  for (i = 1; i < 50; ++i) {
      fe_sq(l2, l2);
  }

  /* l1 = l2 * l1 = w ** (2 ** 250 - 2 ** 50 + 2 ** 50 - 1) = w ** (2 ** 250 - 1) */
  fe_mul(l1, l2, l1);

  /* l1 = l1 ** 4 = w ** (2 ** 252 - 4) */
  fe_sq(l1, l1);
  fe_sq(l1, l1);

  /* sqrt = l1 * w_squared = w ** (2 ** 252 - 4 + 2) = w ** (2 ** 252 - 2)*/
  fe_mul(sqrt, l1, w_squared);

  /* l1 = l1 ** 4 = w ** (2 ** 254 - 16) */
  fe_sq(l1, l1);
  fe_sq(l1, l1);

  /* Recall l0 = w ** 6; out = w ** (2 ** 254 - 10) */
  fe_mul(e, l1, l0);
}

/*
 * Compute the Legendre symbol e of w as
 * w ** ((p-1)/2) = w ** (2 ** 254 - 10) with the exponent as
 * 2 ** 254 - 10 = (2 ** 4) * (2 ** 250 - 1) + 6.
 * We will do so using 253 squarings and 11 multiplications.
 * In the process, also compute w ** (p+3)/8 = (2 ** 252 - 2),
 * which costs us one extra multiplication.
 */
static void fe_fraction_sqrt_and_legendre(fe sqrt, fe e, const fe w)
{
  fe w_squared;
  fe l0;
  fe l1;
  fe l2;
  fe l3;
  int i;

  /* ww = w ** 2 */
  fe_sq(w_squared, w);

  /* l0 = w ** 2 */
  fe_copy(l0, w_squared);

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

  /* l2 = l1 ** (2 ** 5) = w ** ((2 ** 5) * (2 ** 5 - 1)) = w ** (2 ** 10 - 2 ** 5)*/
  fe_sq(l2, l1);
  for (i = 1; i < 5; ++i) {
      fe_sq(l2, l2);
  }

  /* l1 = l2 * l1 = w ** (2 ** 10 - 2 ** 5 + 2 ** 5 - 1) = w ** (2 ** 10 - 1) */
  fe_mul(l1, l2, l1);

  /* l2 = l1 ** (2 ** 10) = w ** ((2 ** 10) * (2 ** 10 - 1)) = w ** (2 ** 20 - 2 ** 10) */
  fe_sq(l2, l1);
  for (i = 1; i < 10; ++i) {
      fe_sq(l2, l2);
  }

  /* l2 = l2 * l1 = w ** (2 ** 20 - 2 ** 10 + 2 ** 10 - 1) = w ** (2 ** 20 - 1) */
  fe_mul(l2, l2, l1);

  /* l3 = l2 ** (2 ** 20) = w ** (2 ** 40 - 2 ** 20)) */
  fe_sq(l3, l2);
  for (i = 1; i < 20; ++i) {
      fe_sq(l3, l3);
  }
  /* l2 = l3 * l2 = w ** (2 ** 40 - 2 ** 20 + 2 ** 20 - 1) = w ** (2 ** 40 - 1) */
  fe_mul(l2, l3, l2);

  /* l2 = l2 ** (2 ** 10) = w ** ((2 ** 10) * (2 ** 40 - 1)) = w ** (2 ** 50 - 2 ** 10) */
  for (i = 0; i < 10; ++i) {
      fe_sq(l2, l2);
  }

  /* l1 = l2 * l1 = w ** (2 ** 50 - 2 ** 10 + 2 ** 10 - 1) = w ** (2 ** 50 - 1) */
  fe_mul(l1, l2, l1);

  /* l2 = l1 ** (2 ** 50) = w ** (2 ** 100 - 2 ** 50) */
  fe_sq(l2, l1);
  for (i = 1; i < 50; ++i) {
      fe_sq(l2, l2);
  }
  /* l2 = l2 ** l1 = w ** (2 ** 100 - 2 ** 50 + 2 ** 50 - 1)) = w ** (2 ** 100 - 1) */
  fe_mul(l2, l2, l1);

  /* l3 = l2 ** (2 ** 100) = w ** (2 ** 200 - 2 ** 100) */
  fe_sq(l3, l2);
  for (i = 1; i < 100; ++i) {
      fe_sq(l3, l3);
  }
  /* l2 = l3 * l2 = w ** (2 ** 200 - 2 ** 100 + 2 ** 100 - 1) = w ** (2 ** 200 - 1) */
  fe_mul(l2, l3, l2);

  /* l2 = l2 ** (2 ** 50) = w ** ((2 ** 50) * (2 ** 200 - 1) = w ** (2 ** 250 - 2 ** 50) */
  fe_sq(l2, l2);
  for (i = 1; i < 50; ++i) {
      fe_sq(l2, l2);
  }

  /* l1 = l2 * l1 = w ** (2 ** 250 - 2 ** 50 + 2 ** 50 - 1) = w ** (2 ** 250 - 1) */
  fe_mul(l1, l2, l1);

  /* l1 = l1 ** 4 = w ** (2 ** 252 - 4) */
  fe_sq(l1, l1);
  fe_sq(l1, l1);
  fe_mul(sqrt, l1, w);

  //fe_copy(sqrt, l1);

  /* l1 = l1 ** 4 = w ** (2 ** 254 - 16) */
  fe_sq(l1, l1);
  fe_sq(l1, l1);

  /* Recall l0 = w ** 6; out = w ** (2 ** 254 - 10) */
  fe_mul(e, l1, l0);
}

static void elligator2_ed25519(ge_p3 *out_point, by out_point_bytes, by out_point_mont, const uint8_t *data,
                  size_t size, by public_key)
{
    static const uint8_t SUITE  = 0x04;
    static const uint8_t ONE    = 0x01;

    // 3. hash(suite || one || pk || alpha)
    uint8_t hash[SHA512_DIGEST_LENGTH] = {0};
    SHA512_CTX hash_ctx;
    SHA512_Init(&hash_ctx);
    SHA512_Update(&hash_ctx, &SUITE, 1);
    SHA512_Update(&hash_ctx, &ONE, 1);
    SHA512_Update(&hash_ctx, public_key, 32);
    SHA512_Update(&hash_ctx, data, size);
    SHA512_Final(hash, &hash_ctx);

    // 4. take first 32 bytes of the hash
    by truncatedHash;
    memcpy(truncatedHash, hash, 32);

    // 7. take highest order bit of truncated hash
    // 8. clear the bit in the source

    //uint8_t x0 = truncatedHash[31] & 0x80;
    truncatedHash[31] &= 0x7f;
    //by_print(truncatedHash);

    fe r;
    fe_frombytes(r, truncatedHash);

    // Montgomery constant A
    fe A;
    fe_A(A);

    // Legendre symbols
    fe e1, e2;

    // Conditional value
    unsigned int b;

    // t0/t1 is y-Edwards
    fe t0, t1;

    // 1 constants for arithmetic / cmoves
    fe o0, o1;
    fe_1(o0);
    fe_1(o1);

    // Square root of u=2
    fe u252m2;
    fe_sqrt_2(u252m2);

    // Square root of A
    fe sqA;
    fe_sqrt_A_plus_2(sqA);

    // u is x-Montgomery
    fe u0, u1, uf;
    fe_0(u0);
    fe_0(u1);

    // w=v^2, where v is y-Montgomery
    fe w, w2;

    // y-Montgomery
    fe sqr, sqr2;

    // y-Edwards
    fe y;

    // (x,y) Edwards complete
    ge_p1p1 p1p1;

    // (x,y) Edwards projective
    ge_p2 p2;

    // 8*p out point
    by out;

    // 1. u   = -A / (1 + 2 * (r ** 2))
    fe_sq2(t0, r);                // t0        = 2 * (r ** 2)
    fe_add(t0, o0, t0);           // t0        = 1 + 2 * (r ** 2)
    fe_invert(t0, t0);            // t0        = 1 / (1 + 2 * (r ** 2))
    fe_sub(u0, u0, A);            // u0        = -A
    fe_mul(u0, u0, t0);           // u0        = -A / (1 + 2 * (r ** 2))

    // 2. v^2 = u * (u**2 + A*u + 1)
    fe_sq(t0, u0);                // t0        = u ** 2
    fe_mul(t1, A, u0);            // t1        = A * u
    fe_add(t0, t0, t1);           // t0        = u**2 + A*u
    fe_add(t0, t0, o0);           // t0        = u**2 + A*u + 1
    fe_mul(w, u0, t0);            // w         = u * (u**2 + A*u + 1)

    // 3. (sqr, e) = (v, legendre(w))
    fe_sqrt_and_legendre(sqr, e1, w);  // (sqr, e1) = (v, lengendre(w))
    fe_copy(e2, e1);              // e2        = e1

    // 4. v   = u             e == 1
    //    w2  = w             e == 1
    //    v   = -A - u        e == -1
    //    w2  = ur^2w         e == -1
    fe_sub(u1, u1, A);            // u1        = -A
    fe_sub(u1, u1, u0);           // u1        = -A - u

    e1[0] = (e1[0] + 1) & 2;
    b = e1[0] >> 1;  // e         = {0,1} for cmove
    fe_cmov(u1, u0, b);
    fe_copy(uf, u1);              // v         = uf

    fe_sq2(w2, r);                // w2        = 2 * r ** 2
    fe_mul(w2, w2, w);            // w2        = 2 * r ** 2 * (u ** 3 + A * u ** 2 + u)
    fe_cmov(w2, w, b);            // w2

    // 5. sqr = sqr           e == 1
    //    sqr = sqr*r2^(1/2)  e == -1
    fe_mul(u252m2, u252m2, r);
    fe_cmov(u252m2, o1, b);
    fe_mul(sqr, sqr, u252m2);     // sqr

    // 6. sqr = sqr             sqr ** 2 == w2
    //    sqr = sqr * 1^(1/2)   sqr ** 2 != w2
    fe_sq(sqr2, sqr);
    b = !(fe_parity(sqr2) ^ fe_parity(w2));    // Inverted for some reason. 0 if they differ (good). 1 if they are the same (bad)
    fe_cmov(o1, sqrtm1, b);
    fe_mul(sqr, sqr, o1);         // sqr

    // 7. Y   = uf - 1
    //    T   = uf + 1
    fe_sub(t0, uf, o0);           // t0        = uf - 1
    fe_add(t1, uf, o0);           // t1        = uf + 1
    fe_copy(p1p1.Y, t0);          // Y         = t0
    fe_copy(p1p1.T, t1);          // T         = t1

    // 8. X   = uf * A^(1/2)
    //    Z   = sqr
    fe_mul(sqA, sqA, uf);         // sqA       = uf * A^(1/2)

    // absolute value of (X/Z)

    fe_copy(p1p1.X, sqA);         // X         = sqA
    fe_copy(p1p1.Z, sqr);         // Z         = sqr

    // 9. pio2 = p1p1
    ge_p1p1_to_p2(&p2, &p1p1);

    // 10. pio2 = 8 * pio2
    ge_p2_dbl(&p1p1, &p2);
    ge_p1p1_to_p2(&p2, &p1p1);
    ge_p2_dbl(&p1p1, &p2);
    ge_p1p1_to_p2(&p2, &p1p1);
    ge_p2_dbl(&p1p1, &p2);
    //ge_p1p1_to_p2(&p2, &p1p1);

    // 11. piopoint = pio2
    //ge_tobytes(out_point, &p2);
    //by_print(out_point);
    ge_p1p1_to_p3(out_point, &p1p1);

    fe minvert;
    fe temp, temp2;
    fe_copy(minvert, sqr);                    // minvert = sqr
    fe_mul(minvert, minvert, out_point->Z);    // minvert = sqr*Z
    fe_copy(temp, out_point->Z);               // temp = Z
    fe_sub(temp, temp, out_point->Y);          // temp = Z-Y
    fe_mul(minvert, minvert, temp);           // minvert = sqr*Z*(Z-Y)
    fe_invert(minvert, minvert);              // minvert = (sqr * Z * (Z-Y))^-1
    fe sign;
    fe_mul(t0, minvert, temp);                // t0 = (sqr * Z)^-1
    fe_mul(sign, t0, out_point->Z);            // sign = (sqr)^-1
    fe_mul(sign, sign, sqA);                  // sign = sqA / sqr
    fe_neg(t1, out_point->X);
    fe_neg(temp2, out_point->T);
    b = fe_parity(sign);
    fe_cmov(out_point->X, t1, b);              // change sign of 8P.X if sign of P.X was wrong
    fe_cmov(out_point->T, temp2, b);
    fe_mul(sign, t0, sqr);                     // sign = (Z)^-1
    fe_mul(t1, sign, out_point->Y);            // t1 = Y/Z
    fe_mul(sign, sign, out_point->X);          // sign = X/Z

    fe_tobytes(out_point_bytes, t1);
    out_point_bytes[31] ^= fe_isnegative(sign) << 7;
    fe_mul(minvert, minvert, sqr);            // minvert = (Z*(Z-Y))^-1
    fe_mul(minvert, minvert, out_point->Z);    // minvert = (Z-Y)^-1
    fe_copy(temp, out_point->Z);               // temp = Z
    fe_add(temp, temp, out_point->Y);          // temp = Z+Y
    fe_mul(minvert, minvert, temp);           // minvert = (Z+Y)/(Z-Y)
    fe_tobytes(out_point_mont, minvert);


}

/*
    Implements elligator2 on curve 25519
    Point returned in projective system through out_point and in
    bytes format through out_point_bytes. Point derived from
    field element r.

    Variables:
    A           - Montgomery constant A
    A3          - A ** 3
    e           - Legendre symbol
    b           - Conditional value
    one         - 1 constant
    sqrt2       - Square root of u=2
    sqrtAp2     - Square root of A+2
    v, v2       - Montgomery y-coodrinates
    w, w2       - Square of v
    y           - Edwards y-coodrinate
    p1p1        - Edwards point, completed system
    p2          - Edwards point, projective system

    cost: 2S (exponentiations have similar cost to fe_sq per bit)
*/
static void ECVRF_hash_to_curve_elligator2_25519(ge_p3 *out_point, uint8_t out_point_bytes[32],
                                                    const fe r)
{
    fe A, A3, e, one, sqrt2, sqrtAp2, v, v2, w, w2, y;
    fe z0, z1, t0, t1, t2, t3, recip, recip2, recip3, recip9, recip21, minvert;
    unsigned int b;
    ge_p1p1 p1p1;
    ge_p2 p2;

    fe_A(A), fe_A3(A3), fe_1(one), fe_sqrt_2(sqrt2), fe_sqrt_A_plus_2(sqrtAp2);

    // recip = 1 + 2 * (r ** 2)
    fe_sq2(recip, r);
    fe_add(recip, one, recip);

    // recip2 = recip ** 2
    fe_sq(recip2, recip);

    // z0 = [A * recip ** 2] - [(recip - 1) * A ** 3]
    // z1 = z0
    fe_sub(z0, recip, one);
    fe_mul(z0, z0, A3);
    fe_mul(z1, recip2, A);
    fe_sub(z0, z1, z0);
    fe_copy(z1, z0);

    // recip3 = recip ** 3
    // recip9 = recip ** 9
    // recip21 = recip ** 21
    fe_mul(recip3, recip2, recip);
    fe_sq(recip9, recip3);
    fe_mul(recip9, recip9, recip3);
    fe_sq(recip21, recip9);
    fe_mul(recip21, recip21, recip3);

    // z0 = z0 * recip ** 21
    fe_mul(z0, z0, recip21);

    // (v, e) = (sqrt(z0), legendre(z0))
    fe_fraction_sqrt_and_legendre(v, e, z0);
    fe_mul(v, v, z1);
    fe_mul(v, v, recip9);

    // e ==  1:  (t2, z0) = (-A, z1)
    // e == -1:  (t2, z0) = (A(1 - recip), (2 * r ** 2) * z1)
    b = ((e[0] + 1) & 2) >> 1;
    fe_copy(t3, A);
    fe_neg(t3, t3);
    fe_sub(t2, one, recip);
    fe_mul(t2, A, t2);
    fe_sq2(z0, r);
    fe_mul(z0, z0, z1);
    fe_cmov(t2, t3, b);
    fe_cmov(z0, z1, b);

    // e ==  1:  v = v
    // e == -1:  v = v * r * sqrt(2)
    fe_mul(sqrt2, sqrt2, r);
    fe_cmov(sqrt2, one, b);
    fe_mul(v, v, sqrt2);

    // (recip ** 3 * v ** 2) == z0:  v = v
    // (recip ** 3 * v ** 2) != z0:  v = v * sqrt(-1)
    fe_sq(v2, v);
    fe_mul(v2, v2, recip3);
    b = fe_ispositive(v2) ^ fe_ispositive(z0);
    fe_cmov(one, sqrtm1, b);
    fe_mul(v, v, one);

    // (x, y) = (|X/Z|, Y/T) is the +/- Edwards curve point before cofactor clearing
    // x = (sqrt(A + 2) * t2) / (v * recip)
    // y = (t2 - recip) / (t2 + recip)
    fe_mul(sqrtAp2, sqrtAp2, t2);
    fe_mul(v, v, recip);
    fe_sub(p1p1.Y, t2, recip);
    fe_add(p1p1.T, t2, recip);
    fe_copy(p1p1.X, sqrtAp2);
    fe_copy(p1p1.Z, v);

    // out_point = 8(x, y)
    ge_p1p1_to_p2(&p2, &p1p1);
    ge_p2_dbl(&p1p1, &p2);
    ge_p1p1_to_p2(&p2, &p1p1);
    ge_p2_dbl(&p1p1, &p2);
    ge_p1p1_to_p2(&p2, &p1p1);
    ge_p2_dbl(&p1p1, &p2);
    ge_p1p1_to_p3(out_point, &p1p1);

    // minvert = 1 / (v * Z * (Z - Y))
    fe_copy(minvert, v);
    fe_mul(minvert, minvert, out_point->Z);
    fe_copy(t0, out_point->Z);
    fe_sub(t0, t0, out_point->Y);
    fe_mul(minvert, minvert, t0);
    fe_invert(minvert, minvert);

    // fix sign of 8(x, y) if (x, y) was negative
    fe_mul(t1, minvert, t0);
    fe_mul(t2, t1, out_point->Z);
    fe_mul(t2, t2, sqrtAp2);
    fe_neg(t3, out_point->X);
    fe_neg(t0, out_point->T);
    b = !fe_ispositive(t2);
    fe_cmov(out_point->X, t3, b);
    fe_cmov(out_point->T, t0, b);
    fe_mul(t2, t1, v);
    fe_mul(t3, t2, out_point->Y);
    fe_mul(t2, t2, out_point->X);
    fe_tobytes(out_point_bytes, t3);
    out_point_bytes[31] ^= fe_isnegative(t2) << 7;
}

/*
    Controler for ECVRF_hash_to_curve_elligator2_25519
    Derives input field element for elligator2 from alpha and the public key
*/
static void ECVRF_alpha_to_curve(ge_p3 *out_point, uint8_t out_point_bytes[32],
                  uint8_t public_key[32], const uint8_t *alpha, size_t alpha_len)
{

  static const uint8_t SUITE  = 0x04;
  static const uint8_t ONE    = 0x01;

  // hash(suite || one || pk || alpha)
  uint8_t hash[SHA512_DIGEST_LENGTH] = {0};
  SHA512_CTX hash_ctx;
  SHA512_Init(&hash_ctx);
  SHA512_Update(&hash_ctx, &SUITE, 1);
  SHA512_Update(&hash_ctx, &ONE, 1);
  SHA512_Update(&hash_ctx, public_key, 32);
  SHA512_Update(&hash_ctx, alpha, alpha_len);
  SHA512_Final(hash, &hash_ctx);

  // take first 32 bytes of the hash
  uint8_t truncatedHash[32];
  memcpy(truncatedHash, hash, 32);

  // take highest order bit of truncated hash
  // clear the bit in the source
  truncatedHash[31] &= 0x7f;

  // field element from the hash
  fe r;
  fe_frombytes(r, truncatedHash);

  // elligator2
  ECVRF_hash_to_curve_elligator2_25519(out_point, out_point_bytes, r);
}
