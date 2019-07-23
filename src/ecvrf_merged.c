#include "../lib/ecvrf_merged.h"
#include <stdio.h>
#include <time.h>

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
 * 2*A, where A = 48662 is the coefficient in the Mongtomery Curve 25519 equation v^2 = u(u^2 + Au + 1)
 */
static void fe_2A(fe out)
{
  out[0] = 973324;
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
static void fe_A_cubed(fe out)
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
 * LC: It seems this function is not needed
 * Returns, via parameter out, a field element equal to the
 * positive square root of A + 2 (= 486664), i.e. (A + 2)^{2^{252}-2} mod p (where p = 2^{255}-19)
 * A = 48662 is the value from the Mongtomery Curve 25519 equation v^2 = u(u^2 + Au + 1)
 * LC: (a) how is "positive" defined and (b) are we sure (A + 2)^{2^{252}-2} mod p is positive and not negative?
 * LC: if we are not sure about (b), just omit that part of the comment, and simply leave "positive square root" and
 * LC: an explanation of how "positive" is defined
 * LC: in fact, we probably don't care if it's positive or negative
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
 * LC: this needs a header like fe_sqrt_A_plus_2
 */
static void fe_sqrt_neg_A_plus_2(fe out)
{
  out[0] = -12222970;
  out[1] = -8312128;
  out[2] = -11511410;
  out[3] = 9067497;
  out[4] = -15300785;
  out[5] = -241793;
  out[6] = 25456130;
  out[7] = 14121551;
  out[8] = -12187136;
  out[9] = 3972024;
}

/**
 * Returns, via parameter out, a field element equal to the
 * fourth root of -4; specifically 2^{2^{252}-2} mod p (where p = 2^{255}-19)
 * LC: could you double check that this really is computed as 2^{2^{252}-2} mod p?
 */
static void fe_fourth_root_of_neg_4(fe out)
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

/**
 * Returns 0 if field element h is negative
 * Returns 1 if field element h is positive
 * LC: how is positive and negative defined?
 * Sign of a field element defned as in
 * "High-speed high-security signatures", D. Bernstein et al., 2007, Section 2.
 */
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

// this will be changed when dif bases added
static void fe_mul121666(fe h, fe f)
{
    int32_t f0 = f[0];
    int32_t f1 = f[1];
    int32_t f2 = f[2];
    int32_t f3 = f[3];
    int32_t f4 = f[4];
    int32_t f5 = f[5];
    int32_t f6 = f[6];
    int32_t f7 = f[7];
    int32_t f8 = f[8];
    int32_t f9 = f[9];
    int64_t h0 = f0 * (int64_t) 121666;
    int64_t h1 = f1 * (int64_t) 121666;
    int64_t h2 = f2 * (int64_t) 121666;
    int64_t h3 = f3 * (int64_t) 121666;
    int64_t h4 = f4 * (int64_t) 121666;
    int64_t h5 = f5 * (int64_t) 121666;
    int64_t h6 = f6 * (int64_t) 121666;
    int64_t h7 = f7 * (int64_t) 121666;
    int64_t h8 = f8 * (int64_t) 121666;
    int64_t h9 = f9 * (int64_t) 121666;
    int64_t carry0;
    int64_t carry1;
    int64_t carry2;
    int64_t carry3;
    int64_t carry4;
    int64_t carry5;
    int64_t carry6;
    int64_t carry7;
    int64_t carry8;
    int64_t carry9;

    carry9 = h9 + (1 << 24); h0 += (carry9 >> 25) * 19; h9 -= carry9 & kTop39Bits;
    carry1 = h1 + (1 << 24); h2 += carry1 >> 25; h1 -= carry1 & kTop39Bits;
    carry3 = h3 + (1 << 24); h4 += carry3 >> 25; h3 -= carry3 & kTop39Bits;
    carry5 = h5 + (1 << 24); h6 += carry5 >> 25; h5 -= carry5 & kTop39Bits;
    carry7 = h7 + (1 << 24); h8 += carry7 >> 25; h7 -= carry7 & kTop39Bits;

    carry0 = h0 + (1 << 25); h1 += carry0 >> 26; h0 -= carry0 & kTop38Bits;
    carry2 = h2 + (1 << 25); h3 += carry2 >> 26; h2 -= carry2 & kTop38Bits;
    carry4 = h4 + (1 << 25); h5 += carry4 >> 26; h4 -= carry4 & kTop38Bits;
    carry6 = h6 + (1 << 25); h7 += carry6 >> 26; h6 -= carry6 & kTop38Bits;
    carry8 = h8 + (1 << 25); h9 += carry8 >> 26; h8 -= carry8 & kTop38Bits;

    h[0] = (int32_t)h0;
    h[1] = (int32_t)h1;
    h[2] = (int32_t)h2;
    h[3] = (int32_t)h3;
    h[4] = (int32_t)h4;
    h[5] = (int32_t)h5;
    h[6] = (int32_t)h6;
    h[7] = (int32_t)h7;
    h[8] = (int32_t)h8;
    h[9] = (int32_t)h9;
}

// this will be changed when dif bases added
static void fe_cswap(fe f, fe g, unsigned int b)
{
    size_t i;

    b = 0-b;
    for (i = 0; i < 10; i++) {
        int32_t x = f[i] ^ g[i];
        x &= b;
        f[i] ^= x;
        g[i] ^= x;
    }
}

/**
 * Compute the Legendre symbol e of w as
 * w ** ((p-1)/2) = w ** (2 ** 254 - 10) with the exponent as
 * 2 ** 254 - 10 = (2 ** 4) * (2 ** 250 - 1) + 6.
 * We will do so using 253 squarings and 11 multiplications.
 * In the process, also compute sqrt_helper = w ** (p-5)/8 = (2 ** 252 - 3),
 * which costs us one extra multiplication (this value will be useful
 * for computing the square root of w, or square root of related values -- for example,
 * if w = n * (d**7), this value helps us compute the square root of
 * of n/d, because sqrt_helper * n * (d ** 3) = (n/d)**((p+3)/8),
 * (because 7(p-5)/8+3 = (p-1)-(p+3)/8, and d ** (p-1) = 1).
 */
static void fe_fraction_sqrt_and_legendre(fe sqrt_helper, fe e, const fe w)
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

  /* sqrt_helper = w ** (2**252 - 3) */
  fe_mul(sqrt_helper, l1, w);


  /* l1 = l1 ** 4 = w ** (2 ** 254 - 16) */
  fe_sq(l1, l1);
  fe_sq(l1, l1);

  /* Recall l0 = w ** 6; e = w ** (2 ** 254 - 10) */
  fe_mul(e, l1, l0);
}

/**
 * Conditional swap.
 * (t0, t1) := (t0, t1)  if b == 0
 * (t0, t1) := (t1, t0)  if b == 1
 *
 * Preconditions: b in {0,1}.
 */
static void cswap(uint8_t *t0, uint8_t *t1, unsigned int b)
{
  b = 0-b;
  uint8_t x = (*t0) ^ (*t1);
  x &= b;
  *t0 ^= x;
  *t1 ^= x;
}

static void ge_p2_to_p3(ge_p3 *r, const ge_p2 *p)
{
    fe_mul(r->X, p->X, p->Z);
    fe_mul(r->Y, p->Y, p->Z);
    fe_sq(r->Z, p->Z);
    fe_mul(r->T, r->X, r->Y);
}

/* r = p */
static void ge_p3_copy(ge_p3 *r, const ge_p3 *p)
{
  fe_copy(r->X, p->X);
  fe_copy(r->Y, p->Y);
  fe_copy(r->Z, p->Z);
  fe_copy(r->T, p->T);
}

/**
 * Conditional move.
 * p := p  if b == 0
 * p := q  if b == 1
 *
 * Preconditions: b in {0,1}.
 */
static void ge_p3_cmov(ge_p3 *p, const ge_p3 *q, unsigned int b)
{
  fe_cmov(p->X, q->X, b);
  fe_cmov(p->Y, q->Y, b);
  fe_cmov(p->Z, q->Z, b);
  fe_cmov(p->T, q->T, b);
}

/**
 * Conditional swap.
 * (p, q) := (p, q)  if b == 0
 * (p, q) := (q, p)  if b == 1
 *
 * Preconditions: b in {0,1}.
 */
static void ge_p3_cswap(ge_p3 *p, ge_p3 *q, unsigned int b)
{
  fe_cswap(p->X, q->X, b);
  fe_cswap(p->Y, q->Y, b);
  fe_cswap(p->Z, q->Z, b);
  fe_cswap(p->T, q->T, b);
}

static void ge_p3_neg(ge_p3 *r, const ge_p3 *p)
{
  fe_neg(r->X, p->X);
  fe_neg(r->T, p->T);
}

/**
 * Retrieve extended Edwards point from bytes in constant time
 *
 * Approcimate cost: 1S per bit (where S is the cost of fe squaring)
 */
static int ge_p3_frombytes(ge_p3 *h, const uint8_t s[32])
{
    fe u;
    fe v;
    fe v3;
    fe vxx;
    fe check;

    fe_frombytes(h->Y, s);
    fe_1(h->Z);
    fe_sq(u, h->Y);
    fe_mul(v, u, d);
    fe_sub(u, u, h->Z); /* u = y^2-1 */
    fe_add(v, v, h->Z); /* v = dy^2+1 */

    fe_sq(v3, v);
    fe_mul(v3, v3, v); /* v3 = v^3 */
    fe_sq(h->X, v3);
    fe_mul(h->X, h->X, v);
    fe_mul(h->X, h->X, u); /* x = uv^7 */

    fe_pow22523(h->X, h->X); /* x = (uv^7)^((q-5)/8) */
    fe_mul(h->X, h->X, v3);
    fe_mul(h->X, h->X, u); /* x = uv^3(uv^7)^((q-5)/8) */

    fe_sq(vxx, h->X);
    fe_mul(vxx, vxx, v);
    fe_sub(check, vxx, u); /* vx^2-u */

    unsigned int b;
    fe status, correct, error;
    fe_1(status);
    fe_1(correct);
    fe_neg(error, status);

    b = fe_isnonzero(check);

    fe check2;
    fe_add(check2, vxx, u);
    fe hX;
    fe_mul(hX, h->X, sqrtm1);

    fe_cmov(correct, error, b);
    fe_cmov(check, check2, b);
    fe_cmov(h->X, hX, b);


    b = fe_isnonzero(check);
    fe_cmov(status, correct, b);

    b = fe_isnegative(h->X) != (s[31] >> 7);
    fe_neg(hX, h->X);
    fe_cmov(h->X, hX, b);

    fe_mul(h->T, h->X, h->Y);
    return fe_isnegative(status);
}

/**
 * Merge two tobytes using the multiplication trick
 * out1 <-- p
 * out2 <-- q
 *
 * Cost: 1S
 */
static void ge_p3_merge_two_tobytes(uint8_t out1[32], uint8_t out2[32],
                                const ge_p3 *p, const ge_p3 *q)
{
  fe x, y, recip;

  fe_mul(recip, p->Z, q->Z);
  fe_invert(recip, recip);

  fe_mul(x, recip, q->Z);
  fe_mul(y, p->Y, x);
  fe_mul(x, p->X, x);
  fe_tobytes(out1, y);
  out1[31] ^= fe_isnegative(x) << 7;

  fe_mul(x, recip, p->Z);
  fe_mul(y, q->Y, x);
  fe_mul(x, q->X, x);
  fe_tobytes(out2, y);
  out2[31] ^= fe_isnegative(x) << 7;
}

/* LC: Either explain this algorithm or point to a paper where it is explained */
static void ge_p3_merge_two_to_montgomery(uint8_t u1[32], uint8_t v1[32], uint8_t u2[32], uint8_t v2[32],
                                              const ge_p3 *p, const ge_p3 *q)
{
    fe sqrtnAp2, up, uq, t0, t1, t2, t3, t4, t5;
    fe_sqrt_neg_A_plus_2(sqrtnAp2);
    fe_sub(t0, p->Z, p->Y);
    fe_sub(t1, q->Z, q->Y);

    fe_mul(t2, t0, t1);
    fe_mul(t3, p->X, q->X);
    fe_mul(t4, t2, t3);
    fe_invert(t4, t4);

    fe_add(t5, p->Z, p->Y);
    fe_mul(t5, t5, t4);
    fe_mul(t5, t5, t3);
    fe_mul(up, t5, t1);
    fe_tobytes(u1, up);

    fe_add(t5, q->Z, q->Y);
    fe_mul(t5, t5, t4);
    fe_mul(t5, t5, t3);
    fe_mul(uq, t5, t0);
    fe_tobytes(u2, uq);

    fe_mul(t5, up, t4);
    fe_mul(t5, t5, t2);
    fe_mul(t5, t5, q->X);
    fe_mul(t5, t5, p->Z);
    fe_mul(t5, t5, sqrtnAp2);
    fe_tobytes(v1, t5);

    fe_mul(t5, uq, t4);
    fe_mul(t5, t5, t2);
    fe_mul(t5, t5, p->X);
    fe_mul(t5, t5, q->Z);
    fe_mul(t5, t5, sqrtnAp2);
    fe_tobytes(v2, t5);
}

/**
 * dedicated / dual Edwards addition
 * r = p + q
 * (X3:Y3:Z3:T3) = (X1:Y1:Z1:T1) + (X2:Y2:Z2:T2)
 *
 * Cost: 8M total
 *
 * "Twisted Edwards Curves Revisited", H. Hisil et al., Section 3.2
 */
static void ge_dedicated_add(ge_p3 *r, const ge_p3 *p, const ge_p3 *q)
{
    fe A, B, C, D, E, F, G, H;

    fe_sub(A, p->Y, p->X);      //  1.  A = (Y1 - X1)
    fe_add(B, q->Y, q->X);      //  2.  B = (Y2 + X2)
    fe_mul(A, A, B);            //  3.  A = (Y1 - X1) * (Y2 + X2)
    fe_add(B, p->Y, p->X);      //  4.  B = (Y1 + X1)
    fe_sub(C, q->Y, q->X);      //  5.  C = (Y2 - X2)
    fe_mul(B, B, C);            //  6.  B = (Y1 + X1) * (Y2 - X2)
    fe_add(C, p->Z, p->Z);      //  7.  C = 2 Z1
    fe_mul(C, C, q->T);         //  8.  C = 2 Z1 * T2
    fe_add(D, p->T, p->T);      //  9.  D = 2 T1
    fe_mul(D, D, q->Z);         // 10.  D = 2 T1 * Z2
    fe_add(E, D, C);            // 11.  E = D + C
    fe_sub(F, B, A);            // 12.  F = B - A
    fe_add(G, B, A);            // 13.  G = B + A
    fe_sub(H, D, C);            // 14.  H = D - C
    fe_mul(r->X, E, F);         // 15.  X3 = E * F
    fe_mul(r->Y, G, H);         // 16.  Y3 = G * H
    fe_mul(r->T, E, H);         // 17.  T3 = E * H
    fe_mul(r->Z, F, G);         // 18.  Z3 = F * G
}

/**
 * Calculate (out1, out2) = (scalar1 * H, scalar2 * H) in constant time
 *
 * Cost: 12M + 4S + 1C per bit of scalar
 */
static void double_scalar_fixed_point_mult(uint8_t out1[32], uint8_t out2[32],
                        const ge_p3 *H, const uint8_t scalar1[32], const uint8_t scalar2[32])
{
  ge_p3 sum[4];
  ge_p1p1 tp1;
  ge_cached tca;

  uint8_t index[4];

  int pos;
  for(pos = 0; pos < 4; pos++){
    index[pos] = pos;
    ge_p3_0(&sum[pos]);
  }

  ge_p3 P;
  ge_p3_copy(&P, H);


  for(pos = 0; pos < 255; ++pos){
    unsigned int b1 = 1 & (scalar1[pos / 8] >> (pos & 7));
    unsigned int b2 = 1 & (scalar2[pos / 8] >> (pos & 7));

    for(int i=2; i>=0; i--){
      unsigned int b = (index[i+1] == (b1 + 2*b2));
      ge_p3_cswap(&sum[i], &sum[i+1], b);
      cswap(&index[i], &index[i+1], b);
    }

    ge_dedicated_add(&sum[0], &sum[0], &P);

    ge_p3_dbl(&tp1, &P);
    ge_p1p1_to_p3(&P, &tp1);

  }

  for(int i=0; i<3; i++){
    for(int j=2; j>=i; j--){
      unsigned int b;
      b = (index[j+1] == i);
      ge_p3_cswap(&sum[j], &sum[j+1], b);
      cswap(&index[j], &index[j+1], b);
    }
  }

  ge_dedicated_add(&sum[1], &sum[1], &sum[3]);
  ge_dedicated_add(&sum[2], &sum[2], &sum[3]);

  ge_p3_merge_two_tobytes(out1, out2, &sum[1], &sum[2]);

  OPENSSL_cleanse(index, sizeof(index));
}

static void montgomery_p2_to_ge_p3(ge_p3 *r, const fe U, const fe V, const fe Z)
{
    fe t0, sqrtnAp2;
    ge_p1p1 p;
    fe_sqrt_neg_A_plus_2(sqrtnAp2);
    fe_add(t0, U, Z);
    fe_mul(p.X, U, sqrtnAp2);
    fe_copy(p.Z, V);
    fe_sub(p.Y, U, Z);
    fe_add(p.T, U, Z);
    ge_p1p1_to_p3(r, &p);
}

/* LC: Either explain this algorithm or point to a paper where it is explained */
/**
 * Recovery formula for Montgomery y-Coordinate
 *
 * "Efficient Elliptic Curve Cryptosystems from a Scalar Multiplication Algorithms with Recovery
 *  of the y-Coordinate on a Montgomery-Form Elliptic Curve", K. Okeya, K. Sakurai, 2001, Section 3
 */
static void montgomery_recover_point(fe U, fe V, fe Z, const uint8_t ub[32], const uint8_t vb[32],
                                const fe U1, const fe Z1, const fe U2, const fe Z2)
{
    fe u, v, AA, T1, T2, T3, T4;
    fe_frombytes(u, ub);
    fe_frombytes(v, vb);
    fe_2A(AA);
    fe_mul(T1, u, Z1);  //  1.  T1 <--  u * Z1
    fe_add(T2, U1, T1); //  2.  T2 <-- U1 + T1
    fe_sub(T3, U1, T1); //  3.  T3 <-- U1 - T1
    fe_mul(T3, T3, T3); //  4.  T3 <-- T3 * T3
    fe_mul(T3, T3, U2); //  5.  T3 <-- T3 * U2
    fe_mul(T1, AA, Z1); //  6.  T1 <-- 2A * Z1
    fe_add(T2, T2, T1); //  7.  T2 <-- T2 + T1
    fe_mul(T4, u, U1);  //  8.  T4 <--  u * U1
    fe_add(T4, T4, Z1); //  9.  T4 <-- T4 + Z1
    fe_mul(T2, T2, T4); // 10.  T2 <-- T2 * T4
    fe_mul(T1, T1, Z1); // 11.  T1 <-- T1 * Z1
    fe_sub(T2, T2, T1); // 12.  T2 <-- T2 - T1
    fe_mul(T2, T2, Z2); // 13.  T2 <-- T2 * Z2
    fe_sub(V, T2, T3);  // 14.   V <-- T2 - T3
    fe_add(T1, v, v);   // 15.  T1 <--  2 * v
    fe_mul(T1, T1, Z1); // 16.  T1 <-- T1 * Z1
    fe_mul(T1, T1, Z2); // 17.  T1 <-- T1 * Z2
    fe_mul(U, T1, U1);  // 18.   U <-- T1 * U1
    fe_mul(Z, T1, Z1);  // 19.   Z <-- T1 * Z1
}

/* LC: Either explain this algorithm or point to a paper where it is explained */
/**
 * Montgomery Ladder
 *
 * Cost: 5M + 4S + 1C
 *
 * "Speeding the Pollard and Elliptic Curve Methods of Factorization", P. Montgomery, Section 10.
 * "Montgomery curves and the Montgomery ladder", D. Bernstein, T. Lange, 2017. Section 4.6
 */
static void montgomery_ladder(fe out_u2, fe out_z2, fe out_u3, fe out_z3,
                                       const uint8_t scalar[32], const size_t scalar_len, const uint8_t in_u[32])
{
    fe u1, u2, z2, u3, z3, tmp0, tmp1;
    uint8_t e[32];
    unsigned swap = 0;
    int pos;
    memcpy(e, scalar, scalar_len);
    fe_frombytes(u1, in_u);
    fe_1(u2);
    fe_0(z2);
    fe_copy(u3, u1);
    fe_1(z3);

    for (pos = 8*scalar_len-1; pos >= 0; --pos) {
        unsigned b = 1 & (e[pos / 8] >> (pos & 7));
        swap ^= b;
        fe_cswap(u2, u3, swap);
        fe_cswap(z2, z3, swap);
        swap = b;
        fe_sub(tmp0, u3, z3);
        fe_sub(tmp1, u2, z2);
        fe_add(u2, u2, z2);
        fe_add(z2, u3, z3);
        fe_mul(z3, tmp0, u2);
        fe_mul(z2, z2, tmp1);
        fe_sq(tmp0, tmp1);
        fe_sq(tmp1, u2);
        fe_add(u3, z3, z2);
        fe_sub(z2, z3, z2);
        fe_mul(u2, tmp1, tmp0);
        fe_sub(tmp1, tmp1, tmp0);
        fe_sq(z2, z2);
        fe_mul121666(z3, tmp1);
        fe_sq(u3, u3);
        fe_add(tmp0, tmp0, z3);
        fe_mul(z3, u1, z2);
        fe_mul(z2, tmp1, tmp0);
    }
    fe_cswap(u2, u3, swap);
    fe_cswap(z2, z3, swap);

    fe_copy(out_u2, u2);
    fe_copy(out_z2, z2);
    fe_copy(out_u3, u3);
    fe_copy(out_z3, z3);


    OPENSSL_cleanse(u1, sizeof(u1));
    OPENSSL_cleanse(e, sizeof(e)); /* LC: How do we decide when cleanse is needed? */
}

static void montgomery_ladder_scalar_mult(ge_p3 *r, const uint8_t *scalar, const size_t scalar_len,
                                                  const uint8_t u[32], const uint8_t v[32])
{
  fe u1, z1, u2, z2;
  fe u_rec, v_rec, z_rec;
  montgomery_ladder(u1, z1, u2, z2, scalar, scalar_len, u);
  montgomery_recover_point(u_rec, v_rec, z_rec, u, v, u1, z1, u2, z2);
  montgomery_p2_to_ge_p3(r, u_rec, v_rec, z_rec);
}

static void montgomery_double_scalar_mult_difference(ge_p3 *r, const ge_p3 *p, const ge_p3 *q,
                                                  const uint8_t *s1, const size_t s1_len,
                                                  const uint8_t *s2, const size_t s2_len)
{
  uint8_t p_u[32], p_v[32], q_u[32], q_v[32];
  ge_cached s2_q_cached;
  ge_p1p1 r_p1p1;
  ge_p3 s1_p, s2_q;
  ge_p3_merge_two_to_montgomery(p_u, p_v, q_u, q_v, p, q);
  montgomery_ladder_scalar_mult(&s1_p, s1, s1_len, p_u, p_v);
  montgomery_ladder_scalar_mult(&s2_q, s2, s2_len, q_u, q_v);
  ge_p3_to_cached(&s2_q_cached, &s2_q);
  ge_sub(&r_p1p1, &s1_p, &s2_q_cached);
  ge_p1p1_to_p3(r, &r_p1p1);
}

/**
 * Implements elligator2 on curve 25519
 * Point returned in projective system through out_point and in
 * bytes format through out_point_bytes. Point is derived from
 * field element r.
 *
 * Cost per bit of r: approximately 2S + 0.2M (S =  squaring, M = fe multiplication)
 * (more precisely, 514 squarings and 45 multplications, with fe_inverse and fe_fraction_sqrt_and_legendre taking up most of the time) */
static void ECVRF_hash_to_curve_elligator2_25519(ge_p3 *out_point, uint8_t out_point_bytes[32], const fe r)
{
    fe A, A_cubed, e, one, rt4_neg_4, sqrtnAp2, v, v_squared, y;
    fe new_w_numerator, w_numerator, t0, t1, t2, t3, u_numerator, u_denom, recip_squared, w_denom, w_denom_cubed, w_denom_to_7;
    unsigned int b;
    ge_p1p1 pre_cofactor_p1p1, p1p1;
    ge_p2 p2;

    fe_A(A), fe_A_cubed(A_cubed), fe_1(one), fe_fourth_root_of_neg_4(rt4_neg_4), fe_sqrt_neg_A_plus_2(sqrtnAp2);

    /* t0 = 2 * (r ** 2) */
    fe_sq2(t0, r);
    /* u_denom = 1 + 2 * (r ** 2) */
    fe_add(u_denom, one, t0);

    /* We want to evaluate the montgomery curve equation w = u (u**2 + Au + 1) for u = -A/u_denom
       Expressing this equation as a fraction, via a common denominator, we have
       w = ([(u_denom - 1) * A ** 3] - [A * u_denom ** 2]) / u_denom ** 3 */

    /* w_numerator = [(u_denom - 1) * A ** 3] - [A * u_denom ** 2] */
    fe_mul(w_numerator, t0, A_cubed);
    fe_sq(recip_squared, u_denom);
    fe_mul(t0, recip_squared, A);
    fe_sub(w_numerator, w_numerator, t0);

    /* w_denom = u_denom ** 3, i.e., the denominator of w */
    fe_mul(w_denom, recip_squared, u_denom);

    /* We now have the numerator and the denominator of w. We need to evaluate the Legendre of w and takes its square root or the square root
       of 2rw (depending on whether is a square) to get v.
      To avoid inversion, we will evaluate Legendre of w_numerator * w_denom**7 (which is the same as Legendre of w); this will also help us with square root.
      This trick is described by Bernstein et al. in "High-speed high-security signatures" in Sec. 5. LC: get a better reference and double check It is easy to see that the Legendre can also be acquired from that routine.
    */

    /* w_denom_cubed = w_denom ** 3 */
    /* w_denom_to_7 = w_denom ** 7 */
    fe_sq(w_denom_cubed, w_denom);
    fe_mul(w_denom_cubed, w_denom_cubed, w_denom);
    fe_sq(w_denom_to_7, w_denom_cubed);
    fe_mul(w_denom_to_7, w_denom_to_7, w_denom);

    /* t0 = w_numerator * w_denom ** 7 */
    fe_mul(t0, w_numerator, w_denom_to_7);

    /* (t1, e) = (t0 ** ((p-5)/8), legendre(t0) = legendre(w)) */
    fe_fraction_sqrt_and_legendre(t1, e, t0);

    /* t1 = t1 * w_numerator * w_denom ** 3 = w ** ((p+3)/8) */
    /* (This equation works out because t0 = [w_numerator ** ((p-1)/8)] * [w_denom ** (7(p-1)/8)], and so
       t1 = [w_numerator ** ((p-1)/8+1)] * [w_denom ** (7(p-1)/8+3)]. Note that 7(p-1)/8+3 = (p-1)-(p+3)/8.
       So w_denom ** (7(p-1)/8+3) = [w_denom ** (p-1)] * [w_denom ** (-(p+3)/8)] = 1 / [w_denom ** ((p+3)/8)].
     */
    fe_mul(t1, t1, w_numerator);
    fe_mul(t1, t1, w_denom_cubed);

    /* Note that now t1 ** 4 = w ** ((p+3)/2) = [w ** ((p-1)/2)] * [w ** 2] = e * [w ** 2] */

    /* If e == 1, we will keep the same u and w.
       If e == -1, then we will replace u by 2 * r**2 * u (== -A-u) and replace w by 2 * r**2 * w.
       Either way, we need the square root of w to get v (also known as "montgomery y coordinate"), which we will obtain from t1.
     */

    /* e ==  1:  (u_numerator, new_w_numerator) = (-A, w_numerator) */
    /* e == -1:  (u_numerator, new_w_numerator) = (A(1 - u_denom), (2 * r ** 2) * w_numerator) */
    b = ((e[0] + 1) & 2) >> 1; /* b == 0 iff e == -1 */
    fe_copy(t3, A);
    fe_neg(t3, t3);
    fe_sub(u_numerator, one, u_denom);
    fe_mul(u_numerator, A, u_numerator);
    fe_sq2(new_w_numerator, r);
    fe_mul(new_w_numerator, new_w_numerator, w_numerator);
    fe_cmov(u_numerator, t3, b);
    fe_cmov(new_w_numerator, w_numerator, b);

    /* e ==  1:  v = t1 */
    /* e == -1:  v = t1 * r * fourth_root_of_negative_4 */
    fe_mul(rt4_neg_4, rt4_neg_4, r);
    fe_cmov(rt4_neg_4, one, b);
    fe_mul(v, t1, rt4_neg_4);

    /* Now we have the following: letting new_w = new_w_numerator/w_denom:
       If e = 1 (i.e., w is a square), then v ** 4 = t1 ** 4 = w ** 2 = new_w ** 2
       If e = -1 (i.e., w is a nonsquare), then v ** 4 = [(t1 * r) ** 4] * [-4]  = -[w ** 2] * [r ** 4] * [-4] = [2 * r**2 * w]^2 = new_w ** 2.
       Thus, in either case, so v ** 2 = plusminus new_w = plusminus new_w_numerator / w_denom.
       We will now fix this plusminus issue by multiplying v  sqrt(-1) if needed.
     */

    /* If v ** 2 has a different sign from new_w (= new_w_numerator / w_denom), then  v = v * sqrt(-1) */
    fe_sq(v_squared, v);
    fe_mul(v_squared, v_squared, w_denom);
    b = fe_ispositive(v_squared) ^ fe_ispositive(new_w_numerator);
    fe_cmov(one, sqrtm1, b);
    fe_mul(v, v, one);

    /* The Montgomery point, up to sign of v, is now (u = u_numerator/u_denom, v) */
    /* Now we convert it to Edwards via
       x = sqrt(-(A+2)) * u / v = (sqrt(-(A + 2)) * u_numerator) / (v * u_denom)
       y = (u-1)/u+1 = (u_numerator - u_denom) / (u_numerator + u_denom) */
    /* We will use the complete coordinate system to represent x as X/Z and y as Y/T thus avoiding inversions */
    fe_mul(pre_cofactor_p1p1.X, sqrtnAp2, u_numerator);
    fe_mul(pre_cofactor_p1p1.Z, v, u_denom);
    fe_sub(pre_cofactor_p1p1.Y, u_numerator, u_denom);
    fe_add(pre_cofactor_p1p1.T, u_numerator, u_denom);

    /* Clear the cofactor: out_point = 8(x, y) */
    ge_p1p1_to_p2(&p2, &pre_cofactor_p1p1);
    ge_p2_dbl(&p1p1, &p2);
    ge_p1p1_to_p2(&p2, &p1p1);
    ge_p2_dbl(&p1p1, &p2);
    ge_p1p1_to_p2(&p2, &p1p1);
    ge_p2_dbl(&p1p1, &p2);
    ge_p1p1_to_p3(out_point, &p1p1);

    /* We will need inverses of both pre_cofactor_p1p1.Z and out_point->Z. To avoid two inversions, we will multiply them and invert the product */
    fe_mul(t0, pre_cofactor_p1p1.Z, out_point->Z);
    fe_invert(t0, t0);

    /* t2 = 1/pre_cofactor_p1p1.Z */
    fe_mul(t2, t0, out_point->Z);

    /* t1 = 1/out_point->Z */
    fe_mul(t1, t0, pre_cofactor_p1p1.Z);

    /* negate 8(x, y) if x in the point (x, y) before cofactor clearing was negative */
    fe_mul(t2, t2, pre_cofactor_p1p1.X);
    fe_neg(t3, out_point->X);
    fe_neg(t0, out_point->T);
    b = !fe_ispositive(t2);
    fe_cmov(out_point->X, t3, b);
    fe_cmov(out_point->T, t0, b);

    /* Convert out_point from p3 representation to bytes, without another inversion */
    fe_mul(t3, t1, out_point->Y);
    fe_mul(t2, t1, out_point->X);
    fe_tobytes(out_point_bytes, t3);
    out_point_bytes[31] ^= fe_isnegative(t2) << 7;
}

/**
 * Controler for ECVRF_hash_to_curve_elligator2_25519
 * Derives input field element for elligator2 from alpha and the public key
 */
static void ECVRF_alpha_to_curve(ge_p3 *out_point, uint8_t out_point_bytes[32],
                  const uint8_t public_key[32], const uint8_t *alpha, size_t alpha_len)
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

static int ECVRF_decode_proof(ge_p3 *Gamma, uint8_t *c, uint8_t *s, const uint8_t* pi)
{
  /*


    2.  let c_string = pi_string[ptLen]...pi_string[ptLen+n-1]

    3.  let s_string =pi_string[ptLen+n]...pi_string[ptLen+n+qLen-1]

    4.  Gamma = string_to_point(gamma_string)

    5.  if Gamma = "INVALID" output "INVALID" and stop.

    6.  c = string_to_int(c_string)

    7.  s = string_to_int(s_string)

    8.  Output Gamma, c, and s
  */

  // 1.  let gamma_string = pi_string[0]...p_string[ptLen-1]
  uint8_t gamma_string[32];
  memcpy(gamma_string, pi, 32);

  int status = ge_p3_frombytes(Gamma, gamma_string);

  memset(c, 0, 32);

  memcpy(c, pi+32, 16);
  memcpy(s, pi+48, 32);

  return status;

}

static int ECVRF_decode_proof_vartime(ge_p3 *Gamma, uint8_t *c, uint8_t *s, const uint8_t* pi)
{
  /*


    2.  let c_string = pi_string[ptLen]...pi_string[ptLen+n-1]

    3.  let s_string =pi_string[ptLen+n]...pi_string[ptLen+n+qLen-1]

    4.  Gamma = string_to_point(gamma_string)

    5.  if Gamma = "INVALID" output "INVALID" and stop.

    6.  c = string_to_int(c_string)

    7.  s = string_to_int(s_string)

    8.  Output Gamma, c, and s
  */

  // 1.  let gamma_string = pi_string[0]...p_string[ptLen-1]
  uint8_t gamma_string[32];
  memcpy(gamma_string, pi, 32);

  if(ge_frombytes_vartime(Gamma, gamma_string) == -1)
    return -1;

  memset(c, 0, 32);

  memcpy(c, pi+32, 16);
  memcpy(s, pi+48, 32);

  return 1;

}

static int ECVRF_proof_to_hash(uint8_t *beta, const uint8_t *pi)
{
  /*
    1.  D = ECVRF_decode_proof(pi_string)

    2.  If D is "INVALID", output "INVALID" and stop

    3.  (Gamma, c, s) = D

    4.  three_string = 0x03 = int_to_string(3, 1), a single octet with
       value 3

    5.  beta_string = Hash(suite_string || three_string ||
       point_to_string(cofactor * Gamma))

    6.  Output beta_string
  */


  static const uint8_t SUITE  = 0x04;
  static const uint8_t THREE  = 0x03;

  uint8_t gamma_cofactor[32];
  uint8_t c[32], s[32];
  ge_p3 gamma_p3;
  ge_p2 gamma_p2;
  ge_p1p1 gamma_p1p1;
  int status = ECVRF_decode_proof(&gamma_p3, c, s, pi);

  ge_p3_to_p2(&gamma_p2, &gamma_p3);

  ge_p2_dbl(&gamma_p1p1, &gamma_p2);
  ge_p1p1_to_p2(&gamma_p2, &gamma_p1p1);
  ge_p2_dbl(&gamma_p1p1, &gamma_p2);
  ge_p1p1_to_p2(&gamma_p2, &gamma_p1p1);
  ge_p2_dbl(&gamma_p1p1, &gamma_p2);
  ge_p1p1_to_p2(&gamma_p2, &gamma_p1p1);

  ge_tobytes(gamma_cofactor, &gamma_p2);

  uint8_t hash[SHA512_DIGEST_LENGTH] = {0};
  SHA512_CTX hash_ctx;
  SHA512_Init(&hash_ctx);
  SHA512_Update(&hash_ctx, &SUITE, 1);
  SHA512_Update(&hash_ctx, &THREE, 1);
  SHA512_Update(&hash_ctx, gamma_cofactor, 32);
  SHA512_Final(hash, &hash_ctx);

  memcpy(beta, hash, 64);

  return status;
}

static int ECVRF_proof_to_hash_vartime(uint8_t *beta, const uint8_t *pi)
{
  /*
    1.  D = ECVRF_decode_proof(pi_string)

    2.  If D is "INVALID", output "INVALID" and stop

    3.  (Gamma, c, s) = D

    4.  three_string = 0x03 = int_to_string(3, 1), a single octet with
       value 3

    5.  beta_string = Hash(suite_string || three_string ||
       point_to_string(cofactor * Gamma))

    6.  Output beta_string
  */

  uint8_t c[16], s[32];

  static const uint8_t SUITE  = 0x04;
  static const uint8_t THREE  = 0x03;

  uint8_t gamma_cofactor[32];
  ge_p3 gamma_p3;
  ge_p2 gamma_p2;
  ge_p1p1 gamma_p1p1;
  if(ECVRF_decode_proof_vartime(&gamma_p3, c, s, pi) == -1)
    return -1;


  ge_p3_to_p2(&gamma_p2, &gamma_p3);

  ge_p2_dbl(&gamma_p1p1, &gamma_p2);
  ge_p1p1_to_p2(&gamma_p2, &gamma_p1p1);
  ge_p2_dbl(&gamma_p1p1, &gamma_p2);
  ge_p1p1_to_p2(&gamma_p2, &gamma_p1p1);
  ge_p2_dbl(&gamma_p1p1, &gamma_p2);
  ge_p1p1_to_p2(&gamma_p2, &gamma_p1p1);

  ge_tobytes(gamma_cofactor, &gamma_p2);

  SHA512_CTX hash_ctx;
  SHA512_Init(&hash_ctx);
  SHA512_Update(&hash_ctx, &SUITE, 1);
  SHA512_Update(&hash_ctx, &THREE, 1);
  SHA512_Update(&hash_ctx, gamma_cofactor, 32);
  SHA512_Final(beta, &hash_ctx);

  return 1;
}

static void ECVRF_prove(uint8_t pi[80], const uint8_t SK_bytes[32],
                    const uint8_t* alpha, const size_t alpha_len)
{

  //1.  Use SK to derive the VRF secret scalar x and the VRF public key Y = x*B
  SHA512_CTX hash_ctx, nonce_ctx, c_ctx;
  const uint8_t SUITE  = 0x04;
  const uint8_t TWO    = 0x02;
  uint8_t hash[SHA512_DIGEST_LENGTH], nonce[SHA512_DIGEST_LENGTH], c_string[SHA512_DIGEST_LENGTH];
  uint8_t Y_bytes[32], H_bytes[32], G_bytes[32], kB_bytes[32], kH_bytes[32];
  uint8_t c[32], s[32], truncatedHash[32];
  ge_p3 Y, H;

  SHA512_Init(&hash_ctx);
  SHA512_Update(&hash_ctx, SK_bytes, 32);
  SHA512_Final(hash, &hash_ctx);

  memcpy(truncatedHash, hash, 32);

  truncatedHash[0]  &= 0xF8;
  truncatedHash[31] &= 0x7F;
  truncatedHash[31] |= 0x40;

  ge_scalarmult_base(&Y, truncatedHash);
  ge_p3_tobytes(Y_bytes, &Y);

  ECVRF_alpha_to_curve(&H, H_bytes, Y_bytes, alpha, alpha_len);

  SHA512_Init(&nonce_ctx);
  SHA512_Update(&nonce_ctx, hash + 32, 32);
  SHA512_Update(&nonce_ctx, H_bytes, 32);
  SHA512_Final(nonce, &nonce_ctx);


  x25519_sc_reduce(nonce);
  ge_scalarmult_base(&Y, nonce);
  ge_p3_tobytes(kB_bytes, &Y);

  double_scalar_fixed_point_mult(kH_bytes, G_bytes, &H, nonce, truncatedHash);
  //6.  c = ECVRF_hash_points(H, Gamma, k*B, k*H)

  SHA512_Init(&c_ctx);
  SHA512_Update(&c_ctx, &SUITE, 1);
  SHA512_Update(&c_ctx, &TWO, 1);
  SHA512_Update(&c_ctx, H_bytes, 32);
  SHA512_Update(&c_ctx, G_bytes, 32);
  SHA512_Update(&c_ctx, kB_bytes, 32);
  SHA512_Update(&c_ctx, kH_bytes, 32);
  SHA512_Final(c_string, &c_ctx);

  memset(c, 0, 32);
  memcpy(c, c_string, 16);

  //7.  s = (k + c*x) mod q
  sc_muladd(s, truncatedHash, c, nonce);

  //8.  pi_string = point_to_string(Gamma) || int_to_string(c, n) || int_to_string(s, qLen)
  //9.  Output pi_string
  memcpy(pi, G_bytes, 32);
  memcpy(pi+32, c, 16);
  memcpy(pi+48, s, 32);

  OPENSSL_cleanse(hash, sizeof(hash));
  OPENSSL_cleanse(nonce, sizeof(nonce));
  OPENSSL_cleanse(c_string, sizeof(c_string));
  OPENSSL_cleanse(truncatedHash, sizeof(truncatedHash));
  OPENSSL_cleanse(Y_bytes, sizeof(Y_bytes));
  OPENSSL_cleanse(H_bytes, sizeof(H_bytes));
  OPENSSL_cleanse(kB_bytes, sizeof(kB_bytes));
  OPENSSL_cleanse(kH_bytes, sizeof(kH_bytes));
  ge_p3_0(&Y);
  ge_p3_0(&H);

}

static int ECVRF_verify(const uint8_t *Y_bytes, const uint8_t *pi,
                  const uint8_t *alpha, const size_t alpha_len)
{
  /*
    1.  D = ECVRF_decode_proof(pi_string)

    2.  If D is "INVALID", output "INVALID" and stop

    3.  (Gamma, c, s) = D

    4.  H = ECVRF_hash_to_curve(suite_string, Y, alpha_string)

    5.  U = s*B - c*Y

    6.  V = s*H - c*Gamma

    7.  c' = ECVRF_hash_points(H, Gamma, U, V)

    8.  If c and c' are equal, output ("VALID",
        ECVRF_proof_to_hash(pi_string)); else output "INVALID"
  */
  static const uint8_t SUITE  = 0x04;
  static const uint8_t TWO    = 0x02;
  uint8_t cp_string[SHA512_DIGEST_LENGTH];
  uint8_t c[32], s[32], H_bytes[32], U_bytes[32], V_bytes[32];
  SHA512_CTX cp_ctx;
  ge_p3 H, G, U, V, Y;
  ge_p2 U_p2;
    // LC: why not ge_frombytes_vartime -- why waste time on constant-time?
  if(ECVRF_decode_proof_vartime(&G, c, s, pi) == -1)
    return 0;

  ECVRF_alpha_to_curve(&H, H_bytes, Y_bytes, alpha, alpha_len);

  ge_frombytes_vartime(&Y, Y_bytes);
  ge_p3_neg(&Y, &Y);
  ge_double_scalarmult_vartime(&U_p2, c, &Y, s);
  ge_p2_to_p3(&U, &U_p2);

    /* LC: explain why we choose montgomery algorithm here but edwards algorithm in the previous line */
  montgomery_double_scalar_mult_difference(&V, &H, &G, s, 32, c, 16);

  ge_p3_merge_two_tobytes(U_bytes, V_bytes, &U, &V); /* LC: This could be further sped up if we used vartime inversion, but not clear if worth implementing */

  SHA512_Init(&cp_ctx);
  SHA512_Update(&cp_ctx, &SUITE, 1);
  SHA512_Update(&cp_ctx, &TWO, 1);
  SHA512_Update(&cp_ctx, H_bytes, 32);
  SHA512_Update(&cp_ctx, pi, 32);
  SHA512_Update(&cp_ctx, U_bytes, 32);
  SHA512_Update(&cp_ctx, V_bytes, 32);
  SHA512_Final(cp_string, &cp_ctx);

  uint8_t cp[32];
  memset(cp, 0, 32);
  memcpy(cp, cp_string, 16);

  for(uint32_t i=0; i<32; i++)
    if(c[i] != cp[i])
      return 0;
  return 1;
}
