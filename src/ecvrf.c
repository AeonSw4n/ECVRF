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
static void montgomery_ladder_to_edwards(uint8_t out[32],
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
  memcpy(out, x2_by, 32);
  out[31] ^= eq<<7;

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
  fe_copy(P.X, H->X);
  fe_copy(P.Y, H->Y);
  fe_copy(P.Z, H->Z);
  fe_copy(P.T, H->T);

  for(pos = 0; pos < 255; ++pos){
    unsigned int b1 = 1 & (scalar1[pos / 8] >> (pos & 7));
    unsigned int b2 = 1 & (scalar2[pos / 8] >> (pos & 7));

    /*
    ge_p3_to_cached(&tca, &P);
    ge_add(&tp1, &sum[b1 + 2*b2], &tca);
    ge_p1p1_to_p3(&sum[b1 + 2*b2], &tp1);
    */

    ge_dedicated_add(&sum[b1 + 2*b2], &sum[b1 + 2*b2], &P);

    ge_p3_dbl(&tp1, &P);
    ge_p1p1_to_p3(&P, &tp1);

  }

  /*
  ge_p3_to_cached(&tca, &sum[3]);
  ge_add(&tp1, &sum[1], &tca);
  ge_p1p1_to_p3(&sum[1], &tp1);
  ge_add(&tp1, &sum[2], &tca);
  ge_p1p1_to_p3(&sum[2], &tp1);
  */

  ge_dedicated_add(&sum[1], &sum[1], &sum[3]);
  ge_dedicated_add(&sum[2], &sum[2], &sum[3]);

  /*
    ge_p3_tobytes(out1, &sum[1]);
    ge_p3_tobytes(out2, &sum[2]);
  */
  ge_p3_merge_two_tobytes(out1, out2, &sum[1], &sum[2]);
}

static void ECVRF_prove(double *t, uint8_t* pi, const uint8_t* SK,
                    const uint8_t* alpha, const uint8_t alpha_len)
{
#ifdef DEBUG
  printf("----- SK -----\n");
  by_print(SK);
  printf("\n");
#endif

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

#ifdef DEBUG
  printf("----- x -----\n");
  by_print(truncatedHash);
  printf("\n");
#endif

  ge_p3 p3;
  by y;
  ge_scalarmult_base(&p3, truncatedHash);
  ge_p3_tobytes(y, &p3);

#ifdef DEBUG
  printf("----- PK -----\n");
  by_print(y);
  printf("\n");
#endif

  //2.  H = ECVRF_hash_to_curve(suite_string, Y, alpha_string)
  //3.  h_string = point_to_string(H)
  ge_p3 H;
  by H_string, mUb;
  elligator2_ed25519_fast(&H, H_string, mUb, alpha, alpha_len, y);
  //by H_string;
  //ge_p3_tobytes(H_string, &H);

#ifdef DEBUG
  printf("----- H -----\n");
  by_print(H_string);
  printf("\n");
#endif

  //4.  gamma = x*H
  //by gamma;
  //x25519_scalar_mult(gamma, truncatedHash, H);

  //5.  k = nonce = ECVRF_nonce_generation(SK, h_string)
  by kB, kH;
  uint8_t nonce[SHA512_DIGEST_LENGTH] = {0};
  SHA512_CTX nonce_ctx;
  SHA512_Init(&nonce_ctx);
  SHA512_Update(&nonce_ctx, hash + 32, 32);
  SHA512_Update(&nonce_ctx, H_string, 32);
  SHA512_Final(nonce, &nonce_ctx);

#ifdef DEBUG
  printf("----- k -----\n");
  by_print(nonce);
  printf("\n");
#endif

  x25519_sc_reduce(nonce);
  ge_scalarmult_base(&p3, nonce);
  ge_p3_tobytes(kB, &p3);

  *t = (double)clock();

  //montgomery_ladder_to_edwards(kH, nonce, mUb, &H);

  by gamma;
  //montgomery_ladder_to_edwards(gamma, truncatedHash, mUb, &H);
  double_scalar_fixed_point_mult(kH, gamma, &H, nonce, truncatedHash);

  *t = (double)clock()-(*t);

#ifdef DEBUG
  printf("----- U=k*B -----\n");
  by_print(kB);
  printf("\n");
  printf("----- V=k*H -----\n");
  by_print(kH);
  printf("\n");
#endif


  //6.  c = ECVRF_hash_points(H, Gamma, k*B, k*H)
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

#ifdef DEBUG
  printf("----- pi -----\n");
  for(int i=0; i<32; i++)
    printf("%x ", pi[i]);
  printf("|| ");
  for(int i=32; i<48; i++)
    printf("%x ", pi[i]);
  printf("|| ");
  for(int i=48; i<80; i++)
    printf("%x ", pi[i]);
  printf("\n");
#endif
}

static void ECVRF_proof2hash(uint8_t* beta, const uint8_t* pi)
{
  /*
    1.  D = ECVRF_decode_proof(pi)

    2.  If D is "INVALID", output "INVALID" and stop

    3.  (gamma, c, s) = D

    4.  three = 0x03 = I2OSP(3, 1), a single octet with value 3

    5.  preBeta = Hash(suite || three || EC2OSP(gamma^cofactor))

    6.  beta = first 2n octets of preBeta

    7.  Output beta
  */
}

static uint8_t ECVRF_verify(const by y, const uint8_t* pi,
                  const uint8_t* alpha)
{
  /*
    1.  D = ECVRF_decode_proof(pi)

    2.  If D is "INVALID", output "INVALID" and stop

    3.  (gamma, c, s) = D

    4.  u = g^s / y^c (where / denotes EC point subtraction, i.e. the
       group operation applied to g^s and the inverse of y^c)

    5.  h = ECVRF_hash_to_curve(suite, y, alpha)

    6.  v = h^s / gamma^c (where / again denotes EC point subtraction)

    7.  c' = ECVRF_hash_points(h, gamma, u, v)

    8.  If c and c' are equal, output "VALID"; else output "INVALID"
  */
  return 0;
}
