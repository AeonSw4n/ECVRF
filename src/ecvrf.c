#include <stdio.h>
#include "../lib/ecvrf.h"


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

static void montgomery_ladder(uint8_t out_x2[32], uint8_t out_x3[32], unsigned *out_b,
                                       const uint8_t scalar[32],
                                       const uint8_t in_x[32]) {
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
    //printf("\n");
    //printf("swap: %d \n", swap);
    fe_invert(z2, z2);
    fe_mul(x2, x2, z2);
    fe_invert(z3, z3);
    fe_mul(x3, x3, z3);
    fe_tobytes(out_x2, x2);
    fe_tobytes(out_x3, x3);
    *out_b = swap;

    OPENSSL_cleanse(e, sizeof(e));
}


static void double_add_scalar_mult(uint8_t out[32],
                                        const uint8_t scalar[32],
                                        const uint8_t point[32])
{
    ge_cached pc;
    ge_p3 pp3;
    ge_p1p1 pp1p1;
    ge_frombytes_vartime(&pp3, point);
    ge_p3_to_cached(&pc, &pp3);
    ge_p3 pout;
    for(int i=31;i>=0;i--){
      for(int j=7;j>=0;j--){
        if(i==31 && j==7)
          continue;
        ge_p3_dbl(&pp1p1, &pp3);
        ge_p1p1_to_p3(&pp3, &pp1p1);
        if((scalar[i]>>j) & 1 == 1){
          printf("1");
          fe_copy(pout.X, pp3.X);
          fe_copy(pout.Y, pp3.Y);
          fe_copy(pout.Z, pp3.Z);
          fe_copy(pout.T, pp3.T);
          ge_add(&pp1p1, &pp3, &pc);
          ge_p1p1_to_p3(&pp3, &pp1p1);
        }
        printf("0");
      }
    }
    printf("\n");
    ge_p3_tobytes(out, &pout);
}


static void fast_add(ge_p3 *out, ge_p3 *in, ge_p2 *h)
{
  fe A,B,C,D,E,F,G,H;
  fe_sub(A, in->Y, in->X);
  fe_add(B, h->Y, h->X);
  fe_mul(A, A, B);
  fe_add(B, in->Y, in->X);
  fe_sub(C, h->Y, h->X);
  fe_mul(B, B, C);
  fe_add(C, in->Z, in->Z);
  //fe_mul(C, C, h->T)
}


static void double_scalar_fixed_point_mult(ge_p3 out1, ge_p3 out2, ge_p2 H,
                                              const uint8_t scalar1[32], const uint8_t scalar2[32])
{
  ge_p3 sum[4];
  int pos;
  for(pos = 0; pos < 4; ++pos)
    ge_p3_0(&sum[pos]);


  for(pos = 0; pos < 255; ++pos){
    unsigned b1 = 1 & (scalar1[pos / 8] >> (pos & 7));
    unsigned b2 = 1 & (scalar2[pos / 8] >> (pos & 7));
    ge_p1p1 temp;
    ge_p2_dbl(&temp, &H);
    ge_p1p1_to_p2(&H, &temp);


  }
}

static void ECVRF_prove(uint8_t* pi, const uint8_t* SK,
                    const uint8_t* alpha, const uint8_t alpha_len)
{
  printf("----- SK -----\n");
  by_print(SK);
  printf("\n");
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

  printf("----- x -----\n");
  by_print(truncatedHash);
  printf("\n");

  ge_p3 p3;
  by y;
  ge_scalarmult_base(&p3, truncatedHash);
  ge_p3_tobytes(y, &p3);
  printf("----- PK -----\n");
  by_print(y);
  printf("\n");

  //2.  H = ECVRF_hash_to_curve(suite_string, Y, alpha_string)
  //3.  h_string = point_to_string(H)
  ge_p3 H;
  by H_string, mUb;
  elligator2_ed25519(&H, H_string, mUb, alpha, alpha_len, y);
  //by H_string;
  //ge_p3_tobytes(H_string, &H);
  printf("----- H -----\n");
  by_print(H_string);
  printf("\n");

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


  printf("----- k -----\n");
  by_print(nonce);
  x25519_sc_reduce(nonce);
  ge_scalarmult_base(&p3, nonce);
  ge_p3_tobytes(kB, &p3);

  fe t0, t1;
  fe_1(t1);
  /*fe Hy;
  fe_frombytes(Hy, H_string);
  fe t0, t1;
  fe U, Z;
  fe_1(t1);
  fe_sub(Z, t1, Hy);    // Z = 1 - Y
  fe_add(U, t1, Hy);    // U = 1 + Y
  fe U2, Z2;
  fe f4;
  fe_1(f4);
  fe_add(f4, f4, f4);
  fe_add(f4, f4, f4);   // f4 = 4
  fe s8;
  fe_add(s8, f4, f4);
  by s8b;
  fe_tobytes(s8b, s8);
  fe mU;
  by mUb;
  fe_invert(mU, Z);
  fe_mul(mU, mU, U);
  fe_tobytes(mUb, mU);*/
  by out_x2, out_x3;  //montgomery |nP|, |(n+1)P|
  unsigned out_b = 0;
  montgomery_ladder(out_x2, out_x3, &out_b, nonce, mUb);

  fe x2_fe, x3_fe;
  fe_frombytes(x2_fe, out_x2);
  fe_add(t0, x2_fe, t1);
  fe_invert(t0, t0);
  fe_sub(x2_fe, x2_fe, t1);
  fe_mul(x2_fe, x2_fe, t0);
  by x2_by, x3_by;
  fe_tobytes(x2_by, x2_fe);
  ge_p3 x2_p3, x3_p3;
  ge_p1p1 x3_p1;
  ge_frombytes_vartime(&x2_p3, x2_by);
  print_fe(x2_p3.X);
  //fe_neg(x2_p3.X, x2_p3.X);
  //fe_neg(x2_p3.T, x2_p3.T);
  print_fe(x2_p3.X);
  ge_cached p_ca;
  ge_p3 H_p3;
  ge_frombytes_vartime(&H_p3, H_string);
  ge_p3_to_cached(&p_ca, &H);

  ge_p3 x2_p3_test;
  ge_p1p1 x2_p1_test;
  by x2_by_test;
  fe_neg(x2_p3_test.X, x2_p3.X);
  fe_neg(x2_p3_test.T, x2_p3.T);
  fe_copy(x2_p3_test.Y, x2_p3.Y);
  fe_copy(x2_p3_test.Z, x2_p3.Z);
  ge_add(&x2_p1_test, &x2_p3_test, &p_ca);
  ge_p1p1_to_p3(&x2_p3_test, &x2_p1_test);
  ge_p3_tobytes(x2_by_test, &x2_p3_test);
  /*ge_p1p1 ap1, bp1;
  ge_p3 ap3, bp3;
  by aby, bby;
  by wtf;
  ge_p3_tobytes(wtf, &H);
  by_print(wtf);
  ge_add(&ap1, &H, &p_ca);
  ge_p3_dbl(&bp1, &H);
  ge_p1p1_to_p3(&ap3, &ap1);
  ge_p1p1_to_p3(&bp3, &bp1);
  ge_p3_tobytes(aby, &ap3);
  ge_p3_tobytes(bby, &bp3);
  printf("\n");
  by_print(aby);
  by_print(bby);*/
  ge_add(&x3_p1, &x2_p3, &p_ca);
  ge_p1p1_to_p3(&x3_p3, &x3_p1);
  fe_frombytes(x3_fe, out_x3);
  fe_add(t0, x3_fe, t1);
  fe_invert(t0, t0);
  fe_sub(x3_fe, x3_fe, t1);
  fe_mul(x3_fe, x3_fe, t0);
  fe_tobytes(x3_by, x3_fe);

  ge_p3 test_p3, test2_p3;
  ge_p1p1 test_p1p1, test2_p1p1;
  ge_cached test_ca, test2_ca;
  by test_by, test2_by;
  ge_frombytes_vartime(&test_p3, x3_by);
  //ge_p3_tobytes(test_by, &test_p3);
  //by_print(test_by);
  fe_copy(test2_p3.X, test_p3.X);
  fe_copy(test2_p3.Y, test_p3.Y);
  fe_copy(test2_p3.Z, test_p3.Z);
  fe_copy(test2_p3.T, test_p3.T);
  fe_neg(x2_p3.X, x2_p3.X);
  fe_neg(x2_p3.T, x2_p3.T);
  ge_p3_to_cached(&test_ca, &x2_p3);
  fe_neg(x2_p3.X, x2_p3.X);
  fe_neg(x2_p3.T, x2_p3.T);
  ge_p3_to_cached(&test2_ca, &x2_p3);
  ge_add(&test_p1p1, &test_p3, &test_ca);
  ge_add(&test2_p1p1, &test2_p3, &test2_ca);
  ge_p1p1_to_p3(&test_p3, &test_p1p1);
  ge_p1p1_to_p3(&test2_p3, &test2_p1p1);
  ge_p3_tobytes(test_by, &test_p3);
  ge_p3_tobytes(test2_by, &test2_p3);
  printf("\n HERE \n");
  by_print(test_by);
  by_print(test2_by);

  by x3n_by;
  ge_p3_tobytes(x3n_by, &x3_p3);
  x3n_by[31] &= 0x7f;
  unsigned eq = 1;
  printf("\n");
  by_print(x3n_by);
  by_print(x3_by);
  by_print(x2_by_test);
  for(int i=0; i<32; i++){
    if(x3_by[i] != x3n_by[i]){
      printf("\nUNEQ ON I: %d\n",i);
      eq = 0;
      break;
    }
  }

  if(eq == 1){
    if(out_b == 0){
      printf("case 1\n");
      ge_p3_tobytes(kH, &x2_p3);
    }
    else{
      printf("case 2\n");
      //ge_p3_tobytes(kH, &x3_p3);
      ge_p3_tobytes(kH, &x2_p3);
    }
  }
  else {
    if(out_b == 0){
      printf("case 3\n");
      ge_p3_tobytes(kH, &x2_p3);
      kH[31] ^= 0x80;
    }
    else{
      printf("0/1 case\n");
      /*by_print(x2_by);
        by_print(x3_by);
        fe_neg(x2_p3.X, x2_p3.X);
      /*ge_add(&x3_p1, &x2_p3, &p_ca);
      ge_p1p1_to_p3(&x3_p3, &x3_p1);
      ge_p3_tobytes(kH, &x3_p3);*/
      ge_p3_tobytes(kH, &x2_p3);
      kH[31] ^= 0x80;
    }
  }

  by gamma;
  by out_g2, out_g3;
  montgomery_ladder(out_g2, out_g3, &out_b, truncatedHash, mUb);
  fe_1(t1);
  fe g2_fe, g3_fe;
  fe_frombytes(g2_fe, out_g2);
  fe_add(t0, g2_fe, t1);
  fe_invert(t0, t0);
  fe_sub(g2_fe, g2_fe, t1);
  fe_mul(g2_fe, g2_fe, t0);
  by g2_by, g3_by;
  fe_tobytes(g2_by, g2_fe);
  ge_p3 g2_p3, g3_p3;
  ge_p1p1 g3_p1;
  ge_frombytes_vartime(&g2_p3, g2_by);
  ge_add(&g3_p1, &g2_p3, &p_ca);
  ge_p1p1_to_p3(&g3_p3, &g3_p1);
  fe_frombytes(g3_fe, out_g3);
  fe_add(t0, g3_fe, t1);
  fe_invert(t0, t0);
  fe_sub(g3_fe, g3_fe, t1);
  fe_mul(g3_fe, g3_fe, t0);
  fe_tobytes(g3_by, g3_fe);
  by g3n_by;
  ge_p3_tobytes(g3n_by, &g3_p3);
  g3n_by[31] &= 0x7f;
  eq = 1;
  for(int i=0; i<32; i++){
    if(g3_by[i] != g3n_by[i]){
      printf("\nUNEQ Gamma ON I: %d\n",i);
      eq = 0;
      break;
    }
  }

  if(eq == 1){
    if(out_b == 0){
      printf("case 1\n");
      ge_p3_tobytes(gamma, &g2_p3);
    }
    else{
      printf("case 2\n");
      ge_p3_tobytes(gamma, &g3_p3);
    }
  }
  else {
    if(out_b == 0){
      printf("case 3\n");
      ge_p3_tobytes(gamma, &g2_p3);
      gamma[31] ^= 0x80;
    }
    else{
      printf("0/1 case\n");
      fe_neg(g2_p3.X, g2_p3.X);
      ge_add(&g3_p1, &g2_p3, &p_ca);
      ge_p1p1_to_p3(&g3_p3, &g3_p1);
      ge_p3_tobytes(gamma, &g3_p3);
    }
  }
  //printf("8H\n");
  //by_print(YH);
  //by_print(Y2);

  //by_print(nonce);
  printf("----- U=k*B -----\n");
  by_print(kB);
  printf("\n");
  printf("----- V=k*H -----\n");
  by_print(kH);
  printf("\n");


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
