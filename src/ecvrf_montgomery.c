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

  #ifdef DEBUG
    printf("----- r -----\n");
    by rby;
    fe_tobytes(rby, r);
    by_print(rby);
    printf("\n");
  #endif

    // Montgomery constant A
    fe A;
    elligator_fe_A(A);

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
    fe_u252m2(u252m2);

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

static void ECVRF_prove_montgomery(double *t, uint8_t* pi, const uint8_t* SK,
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
  elligator2_ed25519(&H, H_string, mUb, alpha, alpha_len, y);
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

  montgomery_ladder_to_edwards(kH, nonce, mUb, &H);

  by gamma;
  montgomery_ladder_to_edwards(gamma, truncatedHash, mUb, &H);
  //double_scalar_fixed_point_mult(kH, gamma, &H, nonce, truncatedHash);

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
