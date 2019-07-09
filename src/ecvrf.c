#include <stdio.h>
#include <time.h>
#include "../lib/ecvrf.h"

//#define DEBUG

/*
      F - finite field

      2n = 32 - length, in octets, of a field element in F; must be even

      E - elliptic curve (EC) defined over F

      m = 32 - length, in octets, of an EC point encoded as an octet string

      G - subgroup of E of large prime order

      q - prime order of group G; must be less than 2^{2n}

      cofactor = 8 - number of points on E divided by q

      g - generator of group G

      Hash - cryptographic hash function

      hLen - output length in octets of Hash; must be at least 2n

      suite = 0x03 - a single nonzero octet specifying the ECVRF ciphersuite,
      which determines the above options

*/

static void int_cmov(int *t0, int t1, unsigned int b)
{
  b = 0-b;
  int x = (*t0) ^ t1;
  x &= b;
  *t0 = (*t0)^x;
}

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

static void ge_p3_copy(ge_p3 *r, ge_p3 *p)
{
  fe_copy(r->X, p->X);
  fe_copy(r->Y, p->Y);
  fe_copy(r->Z, p->Z);
  fe_copy(r->T, p->T);
}

static void ge_p3_cmov(ge_p3 *p, ge_p3 *q, unsigned int b)
{
  fe_cmov(p->X, q->X, b);
  fe_cmov(p->Y, q->Y, b);
  fe_cmov(p->Z, q->Z, b);
  fe_cmov(p->T, q->T, b);
}

static void ge_p3_cswap(ge_p3 *p, ge_p3 *q, unsigned int b)
{
  fe_cswap(p->X, q->X, b);
  fe_cswap(p->Y, q->Y, b);
  fe_cswap(p->Z, q->Z, b);
  fe_cswap(p->T, q->T, b);
}

/*
  dedicated / dual Edwards addition
  r = p + q
  (X3:Y3:Z3:T3) = (X1:Y1:Z1:T1) + (X2:Y2:Z2:T2)
  cost: 8M
*/
static void ge_dedicated_add(ge_p3 *r, ge_p3 *p, ge_p3 *q)
{
  fe A, B, C, D, E, F, G, H;

  fe_sub(A, p->Y, p->X);      // A = (Y1 - X1)
  fe_add(B, q->Y, q->X);      // B = (Y2 + X2)
  fe_mul(A, A, B);            // A = (Y1 - X1) * (Y2 + X2)
  fe_add(B, p->Y, p->X);      // B = (Y1 + X1)
  fe_sub(C, q->Y, q->X);      // C = (Y2 - X2)
  fe_mul(B, B, C);            // B = (Y1 + X1) * (Y2 - X2)
  fe_add(C, p->Z, p->Z);      // C = 2 Z1
  fe_mul(C, C, q->T);         // C = 2 Z1 * T2
  fe_add(D, p->T, p->T);      // D = 2 T1
  fe_mul(D, D, q->Z);         // D = 2 T1 * Z2
  fe_add(E, D, C);            // E = D + C
  fe_sub(F, B, A);            // F = B - A
  fe_add(G, B, A);            // G = B + A
  fe_sub(H, D, C);            // H = D - C
  fe_mul(r->X, E, F);         // X3 = E * F
  fe_mul(r->Y, G, H);         // Y3 = G * H
  fe_mul(r->T, E, H);         // T3 = E * H
  fe_mul(r->Z, F, G);         // Z3 = F * G
}

/*
  merging two tobytes inversions using the multiplication trick
  uint8_t[32] out1 <-- p
  uint8_t[32] out2 <-- q
  cost: 1S
*/
static void ge_p3_merge_two_tobytes(by out1, by out2, ge_p3 *p, ge_p3 *q)
{
  fe recip;
  fe x;
  fe y;

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

static void double_scalar_fixed_point_mult(uint8_t out1[32], uint8_t out2[32], ge_p3 *H,
                                              const uint8_t scalar1[32], const uint8_t scalar2[32])
{
  ge_p3 sum[4];
  ge_p1p1 tp1;
  ge_cached tca;

  int pos;
  for(pos = 0; pos < 4; pos++)
    ge_p3_0(&sum[pos]);

  ge_p3 P;
  ge_p3_copy(&P, H);

  unsigned int swap1 = 0;
  unsigned int swap2 = 0;

  for(pos = 0; pos < 255; ++pos){
    unsigned int b1 = 1 & (scalar1[pos / 8] >> (pos & 7));
    unsigned int b2 = 1 & (scalar2[pos / 8] >> (pos & 7));

    ge_p3_cswap(&sum[0], &sum[3], (b1 && b2) ^ (swap1 && swap2));
    ge_p3_cswap(&sum[1], &sum[2], swap2 ^ b2);
    swap1 = b1 ^ b2;
    ge_dedicated_add(&sum[swap1], &sum[swap1], &P);
    swap1 = b1;
    swap2 = b2;

    ge_p3_dbl(&tp1, &P);
    ge_p1p1_to_p3(&P, &tp1);


  }
  ge_p3_cswap(&sum[1], &sum[2], swap1 && swap2);
  ge_p3_cswap(&sum[1], &sum[2], swap2);

  ge_dedicated_add(&sum[1], &sum[1], &sum[3]);
  ge_dedicated_add(&sum[2], &sum[2], &sum[3]);

  ge_p3_merge_two_tobytes(out1, out2, &sum[1], &sum[2]);
}

static void double_scalar_fixed_point_mult_secure(uint8_t out1[32], uint8_t out2[32], ge_p3 *H,
                                              const uint8_t scalar1[32], const uint8_t scalar2[32])
{
  ge_p3 sum[4];
  ge_p1p1 tp1;
  ge_cached tca;

  int pos;
  for(pos = 0; pos < 4; pos++)
    ge_p3_0(&sum[pos]);

  ge_p3 P;
  ge_p3_copy(&P, H);

  unsigned int swap1 = 0;

  for(pos = 0; pos < 255; ++pos){
    unsigned int b1 = 1 & (scalar1[pos / 8] >> (pos & 7));
    unsigned int b2 = 1 & (scalar2[pos / 8] >> (pos & 7));

    swap1 ^= b1;
    ge_p3_cswap(&sum[2], &sum[3], swap1);
    ge_p3_cswap(&sum[1], &sum[2], b2);
    ge_p3_cswap(&sum[0], &sum[1], b1 || b2);
    swap1 = b1;

    ge_dedicated_add(&sum[0], &sum[0], &P);

    ge_p3_cswap(&sum[0], &sum[1], b1 || b2);
    ge_p3_cswap(&sum[1], &sum[2], b2);

    ge_p3_dbl(&tp1, &P);
    ge_p1p1_to_p3(&P, &tp1);


  }

  ge_dedicated_add(&sum[1], &sum[1], &sum[3]);
  ge_dedicated_add(&sum[2], &sum[2], &sum[3]);

  ge_p3_merge_two_tobytes(out1, out2, &sum[1], &sum[2]);
}

static void double_scalar_fixed_point_mult_plain(uint8_t out1[32], uint8_t out2[32], ge_p3 *H,
                                              const uint8_t scalar1[32], const uint8_t scalar2[32])
{
  ge_p3 sum[4];
  ge_p1p1 tp1;
  ge_cached tca;

  int pos;
  for(pos = 0; pos < 4; pos++)
    ge_p3_0(&sum[pos]);

  ge_p3 P;
  ge_p3_copy(&P, H);


  for(pos = 0; pos < 255; ++pos){
    unsigned int b1 = 1 & (scalar1[pos / 8] >> (pos & 7));
    unsigned int b2 = 1 & (scalar2[pos / 8] >> (pos & 7));

    ge_dedicated_add(&sum[b1+2*b2], &sum[b1+2*b2], &P);

    ge_p3_dbl(&tp1, &P);
    ge_p1p1_to_p3(&P, &tp1);

  }

  ge_dedicated_add(&sum[1], &sum[1], &sum[3]);
  ge_dedicated_add(&sum[2], &sum[2], &sum[3]);

  ge_p3_merge_two_tobytes(out1, out2, &sum[1], &sum[2]);
}

static void montgomery_ladder(fe out_x2[32], fe out_z2[32], fe out_x3[32], fe out_z3[32],
                                       const uint8_t scalar[32], const uint8_t in_x[32])
{
    fe x1, x2, z2, x3, z3, tmp0, tmp1;
    uint8_t e[32];
    unsigned swap = 0;
    int pos;

    memcpy(e, scalar, 32);
    //e[0] &= 248;
    //e[31] &= 127;
    //e[31] |= 64;
    fe_frombytes(x1, in_x);
    fe_1(x2);
    fe_0(z2);
    fe_copy(x3, x1);
    fe_1(z3);

    for (pos = 254; pos >= 0; --pos) {
        unsigned b = 1 & (e[pos / 8] >> (pos & 7));
        //printf("%d", b);
        swap ^= b;
        fe_cswap(x2, x3, swap);
        fe_cswap(z2, z3, swap);
        swap = b;
        fe_sub(tmp0, x3, z3);
        fe_sub(tmp1, x2, z2);
        fe_add(x2, x2, z2);
        fe_add(z2, x3, z3);
        fe_mul(z3, tmp0, x2);
        fe_mul(z2, z2, tmp1);
        fe_sq(tmp0, tmp1);
        fe_sq(tmp1, x2);
        fe_add(x3, z3, z2);
        fe_sub(z2, z3, z2);
        fe_mul(x2, tmp1, tmp0);
        fe_sub(tmp1, tmp1, tmp0);
        fe_sq(z2, z2);
        fe_mul121666(z3, tmp1);
        fe_sq(x3, x3);
        fe_add(tmp0, tmp0, z3);
        fe_mul(z3, x1, z2);
        fe_mul(z2, tmp1, tmp0);
    }
    fe_cswap(x2, x3, swap);
    fe_cswap(z2, z3, swap);

    fe_copy(out_x2, x2);
    fe_copy(out_z2, z2);
    fe_copy(out_x3, x3);
    fe_copy(out_z3, z3);

    OPENSSL_cleanse(e, sizeof(e));
}

//can't find any way to combine square rooting
//inversion and square root can be combined by (x^8 y)^(2^252-3)*x^3 for sqr, and then finding y^(2^255-2^252-17)=y^(7*2^252-17)=(y^(2^252-3)^7)*y^4 for the inv
//can save one inversion combining two ladders
//with inversion optimization, total added cost is 3 exponentiations, with inv/sqr total is 2
//tobytes adds 2 more inversions to the total
static void montgomery_ladder_to_edwards(ge_p3 *out,
                                       const uint8_t scalar[32],
                                       const uint8_t in_x[32], const ge_p3 *in_P)
{

  fe x2, z2, x3, z3;
  montgomery_ladder(x2, z2, x3, z3, scalar, in_x);
  by x2_by, x3_by, x3n_by;
  fe x2_ed, x3_ed, x3n_ed;
  fe temp;
  fe t0, t1;
  fe_add(t0, x2, z2);
  fe_add(t1, x3, z3);
  fe_mul(temp, t0, t1);
  fe_invert(temp, temp);

  fe_mul(t1, temp, t1);
  fe_sub(x2_ed, x2, z2);
  fe_mul(x2_ed, t1, x2_ed);

  fe_mul(t0, temp, t0);
  fe_sub(x3_ed, x3, z3);
  fe_mul(x3_ed, t0, x3_ed);

  fe_tobytes(x2_by, x2_ed);

  ge_p3 x2_p3, x3n_p3;
  ge_p1p1 x3n_p1;
  ge_frombytes_vartime(&x2_p3, x2_by);

  ge_cached p_ca;
  ge_p3_to_cached(&p_ca, in_P);

  ge_add(&x3n_p1, &x2_p3, &p_ca);
  ge_p1p1_to_p3(&x3n_p3, &x3n_p1);
  ge_p3_tobytes(x3n_by, &x3n_p3);
  fe_frombytes(x3n_ed, x3n_by);

  fe f;
  fe_sub(f, x3_ed, x3n_ed);
  // eq = 0   when x3_ed - x3n_ed   = 0
  // eq = 1   when x3_ed - x3n_ed  != 0
  unsigned eq = fe_isnonzero(f);
  if(eq == 1){
    fe_neg(x2_p3.X, x2_p3.X);
    fe_neg(x2_p3.T, x2_p3.T);
  }
  ge_p3_copy(out, &x2_p3);

}

static int ge_p3_frombytes(ge_p3 *h, const uint8_t *s)
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

  uint8_t gamma[32], c[16], s[32];

  memcpy(gamma, pi, 32);
  memcpy(c, pi+32, 16);
  memcpy(s, pi+48, 32);

  static const uint8_t SUITE  = 0x04;
  static const uint8_t THREE  = 0x03;

  uint8_t gamma_cofactor[32];
  ge_p3 gamma_p3;
  ge_p2 gamma_p2;
  ge_p1p1 gamma_p1p1;
  if(ge_frombytes_vartime(&gamma_p3, gamma) == -1)
    return -1;


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
  return 1;
}

static void ECVRF_prove(double *t, uint8_t* pi, const uint8_t* SK,
                    const uint8_t* alpha, const size_t alpha_len)
{

  //1.  Use SK to derive the VRF secret scalar x and the VRF public key Y = x*B
  uint8_t hash[SHA512_DIGEST_LENGTH] = {0};
  SHA512_CTX hash_ctx;
  SHA512_Init(&hash_ctx);
  SHA512_Update(&hash_ctx, SK, 32);
  SHA512_Final(hash, &hash_ctx);

  by truncatedHash;
  memcpy(truncatedHash, hash, 32);

  truncatedHash[0]  &= 0xF8;
  truncatedHash[31] &= 0x7F;
  truncatedHash[31] |= 0x40;

  ge_p3 p3;
  by y;
  ge_scalarmult_base(&p3, truncatedHash);
  ge_p3_tobytes(y, &p3);


  //2.  H = ECVRF_hash_to_curve(suite_string, Y, alpha_string)
  //3.  h_string = point_to_string(H)
  ge_p3 H;
  by H_string;
  by mUb;
  elligator2_ed25519_fast(&H, H_string, alpha, alpha_len, y);


  //4.  gamma = x*H

  //5.  k = nonce = ECVRF_nonce_generation(SK, h_string)
  by kB, kH;
  uint8_t nonce[SHA512_DIGEST_LENGTH] = {0};
  SHA512_CTX nonce_ctx;
  SHA512_Init(&nonce_ctx);
  SHA512_Update(&nonce_ctx, hash + 32, 32);
  SHA512_Update(&nonce_ctx, H_string, 32);
  SHA512_Final(nonce, &nonce_ctx);


  x25519_sc_reduce(nonce);
  ge_scalarmult_base(&p3, nonce);
  ge_p3_tobytes(kB, &p3);

  by gamma;
  *t = (double)clock();
  double_scalar_fixed_point_mult(kH, gamma, &H, nonce, truncatedHash);
  *t = (double)clock()-(*t);
  //6.  c = ECVRF_hash_points(H, Gamma, k*B, k*H)
  by_print(kB);
  printf("\n");
  by_print(kH);
  printf("\n");
  static const uint8_t SUITE  = 0x04;
  static const uint8_t TWO    = 0x02;
  uint8_t c_string[SHA512_DIGEST_LENGTH] = {0};
  SHA512_CTX c_ctx;
  SHA512_Init(&c_ctx);
  SHA512_Update(&c_ctx, &SUITE, 1);
  SHA512_Update(&c_ctx, &TWO, 1);
  SHA512_Update(&c_ctx, H_string, 32);
  SHA512_Update(&c_ctx, gamma, 32);
  SHA512_Update(&c_ctx, kB, 32);
  SHA512_Update(&c_ctx, kH, 32);
  SHA512_Final(c_string, &c_ctx);

  uint8_t c[32];
  memset(c, 0, 32);
  memcpy(c, c_string, 16);

  //7.  s = (k + c*x) mod q
  uint8_t s[32];
  sc_muladd(s, truncatedHash, c, nonce);

  //8.  pi_string = point_to_string(Gamma) || int_to_string(c, n) || int_to_string(s, qLen)
  //9.  Output pi_string
  memcpy(pi, gamma, 32);
  memcpy(pi+32, c, 16);
  memcpy(pi+48, s, 32);

}

static int ECVRF_verify(const uint8_t *y, const uint8_t *pi,
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
  uint8_t c[32], s[32];
  ge_p3 G;
  int status = ECVRF_decode_proof(&G, c, s, pi);

  ge_p3 H;
  by H_string;
  elligator2_ed25519_fast(&H, H_string, alpha, alpha_len, y);
  ge_p2 U;
  ge_p3 Y;
  ge_frombytes_vartime(&Y, y);
  fe_neg(Y.X, Y.X);
  fe_neg(Y.T, Y.T);
  ge_double_scalarmult_vartime(&U, c, &Y, s);
  //fe_neg(U.X, U.X);

  ge_p3 s_H, c_G;
  fe mGx, mHx;
  fe t0, t1;
  fe_sub(mGx, G.Z, G.Y);
  fe_sub(mHx, H.Z, H.Y);
  fe_mul(t0, mHx, mGx);
  fe_invert(t0, t0);
  fe_add(t1, G.Z, G.Y);
  fe_mul(t1, t1, t0);
  fe_mul(t1, t1, mHx);  // u of G
  by Gby;
  fe_tobytes(Gby, t1);
  montgomery_ladder_to_edwards(&c_G, c, Gby, &G);
  fe_add(t1, H.Z, H.Y);
  fe_mul(t1, t1, t0);
  fe_mul(t1, t1, mGx);  // u of H
  by Hby;
  fe_tobytes(Hby, t1);
  montgomery_ladder_to_edwards(&s_H, s, Hby, &H);
  by s_Hby;
  ge_p3_tobytes(s_Hby, &s_H);
  //by_print(s);

  //by_print(s_Hby);
  //printf("\n");
  ge_cached c_G_cached;
  ge_p1p1 V_p1p1;
  ge_p2 V;
  ge_p3_to_cached(&c_G_cached, &c_G);
  ge_sub(&V_p1p1, &s_H, &c_G_cached);
  ge_p1p1_to_p2(&V, &V_p1p1);
  by Uby, Vby;
  ge_tobytes(Uby, &U);
  ge_tobytes(Vby, &V);


  //c' = ECVRF_hash_points(H, Gamma, U, V)
  //6.  c = ECVRF_hash_points(H, Gamma, k*B, k*H)
  static const uint8_t SUITE  = 0x04;
  static const uint8_t TWO    = 0x02;
  uint8_t cp_string[SHA512_DIGEST_LENGTH] = {0};
  SHA512_CTX cp_ctx;
  printf("\n");
  by_print(Uby);
  printf("\n");
  by_print(Vby);
  SHA512_Init(&cp_ctx);
  SHA512_Update(&cp_ctx, &SUITE, 1);
  SHA512_Update(&cp_ctx, &TWO, 1);
  SHA512_Update(&cp_ctx, H_string, 32);
  SHA512_Update(&cp_ctx, pi, 32);
  SHA512_Update(&cp_ctx, Uby, 32);
  SHA512_Update(&cp_ctx, Vby, 32);
  SHA512_Final(cp_string, &cp_ctx);

  uint8_t cp[32];
  memset(cp, 0, 32);
  memcpy(cp, cp_string, 16);

  int eq = 1;
  int i;
  for(i=0; i<32; i++){
    if(c[i] != cp[i]){
      eq = 0;
      break;
    }
  }

  return eq;
}
