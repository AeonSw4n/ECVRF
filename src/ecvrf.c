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
                                       const uint8_t point[32]) {
    fe x1, x2, z2, x3, z3, tmp0, tmp1;
    uint8_t e[32];
    unsigned swap = 0;
    int pos;

    memcpy(e, scalar, 32);
    //e[0] &= 248;
    //e[31] &= 127;
    //e[31] |= 64;
    fe_frombytes(x1, point);
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


static void ECVRF_prove(uint8_t* pi, const uint8_t* SK,
                    const uint8_t* alpha)
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
  by H;
  elligator2_ed25519(H, alpha, 2, y);
  printf("----- H -----\n");
  by_print(H);
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
  SHA512_Update(&nonce_ctx, H, 32);
  SHA512_Final(nonce, &nonce_ctx);
  printf("----- k -----\n");
  by_print(nonce);
  x25519_sc_reduce(nonce);
  ge_scalarmult_base(&p3, nonce);
  ge_p3_tobytes(kB, &p3);
  ge_p3 Hp3;
  ge_p1p1 Hp1;
  ge_frombytes_vartime(&Hp3, H);
  ge_p3_dbl(&Hp1, &Hp3);
  ge_p1p1_to_p3(&Hp3, &Hp1);
  ge_p3_dbl(&Hp1, &Hp3);
  ge_p1p1_to_p3(&Hp3, &Hp1);
  ge_p3_dbl(&Hp1, &Hp3);
  ge_p1p1_to_p3(&Hp3, &Hp1);
  by YH;
  ge_p3_tobytes(YH, &Hp3);
/*  fe_copy(Hy2, Hp3.Z);
  fe_invert(Hy2, Hy2);
  fe_mul(Hy2, Hy2, Hp3.Y);
  printf("Y\n");
  //print_fe(Hy2);
  by testH;
  fe_tobytes(testH, Hy2);
  by_print(testH);
  fe_tobytes(testH, Hy);
  by_print(testH);
  printf("\n");*/
  fe Hy;
  fe_frombytes(Hy, H);
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
  fe_tobytes(mUb, mU);
  by out_x2, out_x3;  //montgomery |nP|, |(n+1)P|
  unsigned out_b = 0;
  montgomery_ladder(out_x2, out_x3, &out_b, nonce, mUb);


  for(int i=0;i<3;i++){
    fe_sq(U2, U);         // U2 = U^2
    fe_copy(t0, U2);      // t0 = U^2
    fe_sq(Z2, Z);         // Z2 = Z^2
    fe_sub(U2, U2, Z2);   // U2 = U^2 - Z^2
    fe_sq(U2, U2);        // U2 = (U^2 - Z^2)^2
    fe_add(Z2, Z2, t0);   // Z2 = Z^2 + U^2
    elligator_fe_A(t0);   // t0 = A
    fe_mul(t1, U, Z);     // t1 = U*Z
    fe_mul(t0, t0, t1);   // t0 = A*U*Z
    fe_add(Z2, Z2, t0);   // Z2 = Z^2 + A*U*Z + U^2
    fe_mul(Z2, Z2, t1);   // Z2 = U*Z*(Z^2 + A*U*Z + U^2)
    fe_mul(Z2, Z2, f4);   // Z2 = 4*U*Z*(Z^2 + A*U*Z + U^2)
    fe_copy(U, U2);
    fe_copy(Z, Z2);
  }

  fe_invert(Z, Z);
  fe_mul(U, U, Z);
  fe_1(t1);
  fe_add(Z, U, t1);
  fe_invert(Z, Z);
  fe_sub(U, U, t1);
  fe_mul(U, U, Z);
  fe Y;
  fe_copy(Y, U);
  by Y2;
  fe_tobytes(Y2, Y);

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
  ge_p3 p_p3;
  ge_cached p_ca;
  ge_frombytes_vartime(&p_p3, H);
  ge_p3_to_cached(&p_ca, &p_p3);
  ge_add(&x3_p1, &x2_p3, &p_ca);
  ge_p1p1_to_p3(&x3_p3, &x3_p1);
  fe_frombytes(x3_fe, out_x3);
  fe_add(t0, x3_fe, t1);
  fe_invert(t0, t0);
  fe_sub(x3_fe, x3_fe, t1);
  fe_mul(x3_fe, x3_fe, t0);
  fe_tobytes(x3_by, x3_fe);
  by x3n_by;
  ge_p3_tobytes(x3n_by, &x3_p3);
  x3n_by[31] &= 0x7f;
  unsigned eq = 1;
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
  SHA512_Update(&c_ctx, H, 32);
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
