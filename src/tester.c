#include <stdio.h>
#include "../lib/tester.h"

static void by_print(by s){
  for(uint32_t i=0; i<32; i++){
    printf("%.2x ", s[i]);
  }
  printf("\n");
}

  static void print_fe(fe in)
  {
    for(uint32_t i=0;i<10;i++)
      printf("%d ", in[i]);
    printf("\n");
  }

static void test_fe_legendre()
{

 //TEST FROM BY
  uint8_t* test_str[] = {
                          "c7b5d6239e52a473a2b57a92825e0e5de4656e349bb198de5afd6a76e5a07066",
                          "7ff6d8b773bfbae57b2ab9d49f9d3cb7d9af40a03d3ed3c6beaaf2d486b1fe6e",
                          "2ceaa2c2ff3028c34f9fbe076ff99520b925f18d652285b4daad5ccc467e523b",
                          /*"6670a0e5766afd5ade98b19b346e65e45d0e5e82927ab5a273a4529e23d6b5c7",
                          "6efeb186d4f2aabec6d33e3da040afd9b73c9d9fd4b92a7be5babf73b7d8f67f",
                          "3b527e46cc5caddab48522658df125b92095f96f07be9f4fc32830ffc2a2ea2c",*/
                        };
  for(uint32_t i=0; i<3; i++){
    by a;
    by_fromstr(a, test_str[i]);
    printf("legendre from by test[%d]: %d\n", i+1, fe_legendre_ifromby(a));
  }

 //TEST FROM UINT32_T
  uint32_t count = 10;

  for(uint32_t x=1; x<=2; x++){
    uint32_t mul    = 4;
    uint32_t n      = x;
    for(uint32_t i=1; i<=count; i++){
      by byn;
      by_fromint(byn, n);
      printf("legendre from uint32_t num= %d, outcome: %d\n", n, fe_legendre_ifromby(byn));
      n *= mul;
    }
  }

 //TEST FOR 2ceaa2c2ff3028c34f9fbe076ff99520b925f18d652285b4daad5ccc467e523b
  uint8_t* test_str2[] = {"3", "7", "b", "1f", "7535049c595", "381d6e16ade488bccadf45f453bedf158dc09378d97fdb16c7"};
  uint8_t lens[] = {1, 1, 1, 2, 11, 50};

  fe mul;
  by a;
  by_fromstrc(a, test_str2[0], 1);
  fe_frombytes(mul, a);
  printf("legendre for %s, outcome: %d\n", test_str2[0], fe_legendre_ifromby(a));
  for(int i=1;i<6;i++){
    fe n;
    by b;
    by_fromstrc(b, test_str2[i], lens[i]);
    fe_frombytes(n, b);
    fe_mul(mul, mul, n);
    printf("legendre for %s, outome: %d\n", test_str2[i], fe_legendre_ifromby(b));
  }
  by mulby;
  fe_tobytes(mulby, mul);
  by_print(mulby);
}

static void fe_swap(fe f, fe g, unsigned int b)
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

//TODO CHECK IF WORKING!!!
static unsigned int fe_parity(const fe h){
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

static void test_elligator_2()
{
  fe r;
  uint8_t *input = "ff27b906c1dd409b64b455ee87529f967514bcbbbdbdd45f9c22eabc5de06455";
  int negX = 0;
  int sqrM = 0;
  by in;
  by_fromstr(in, input);
  fe_frombytes(r, in);

  // Montgomery constant A
  fe A;
  elligator_fe_A(A);

  // Legendre symbols
  fe e1, e2;

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
  fe_sqAp2(sqA);

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

  // 1. u = -A / (1 + 2 * (r ** 2))
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

  fe_tobytes(in, w);
  by_print(in);

  // 3. (sqr, e) = (v, legendre(w))
  fe_sqr_legendre(sqr, e1, w);  // (sqr, e1) = (v, lengendre(w))
  fe_copy(e2, e1);              // e2        = e1

  // 4. v   = u             e == 1
  //    w2  = w             e == 1
  //    sqr = sqr           e == 1
  //    v   = -A - u        e == -1
  //    w2  = ur^2w         e == -1
  //    sqr = sqr*r2^(1/2)  e == -1
  fe_sub(u1, u1, A);            // u1        = -A
  fe_sub(u1, u1, u0);           // u1        = -A - u

  e1[0] = (e1[0] + 1) & 2;
  unsigned int b = e1[0] >> 1;  // e         = {0,1} for cmove
  fe_cmov(u1, u0, b);
  fe_copy(uf, u1);              // v         = uf

  fe_sq2(w2, r);                // w2        = 2 * r ** 2
  fe_mul(w2, w2, w);            // w2        = 2 * r ** 2 * (u ** 3 + A * u ** 2 + u)
  fe_cmov(w2, w, b);            // w2


  fe_mul(u252m2, u252m2, r);
  fe_cmov(u252m2, o1, b);
  fe_mul(sqr, sqr, u252m2);     // sqr

  // 5. 
  fe_sq(sqr2, sqr);
  b = !(fe_parity(sqr2) ^ fe_parity(w2)); // 0 if they differ (good). 1 if they are the same (bad)
  fe_cmov(o1, sqrtm1, b);
  fe_mul(sqr, sqr, o1);


  e2[0] = -e2[0];
  fe_mul(sqr, sqr, e2);

  fe_neg(sqA, sqA);
  fe_mul(sqA, sqA, uf);

  ge_p1p1 piop1;
  ge_p2 piop2;

  fe_sub(t0, uf, o0); // t0 = uf - 1
  fe_add(t1, uf, o0); // t1 = uf + 1

  fe_copy(piop1.Y, t0);
  fe_copy(piop1.T, t1);
  fe_copy(piop1.X, sqA);
  fe_copy(piop1.Z, sqr);

  ge_p1p1_to_p2(&piop2, &piop1);
  //fe m2z, m2x, m2y;
  //if(negX)
  //  fe_neg(&piop2.X, &piop2.X);
  //fe_invert(m2z, piop2.Z);
  //fe_mul(m2x, piop2.X, m2z);
  //fe_mul(m2y, piop2.Y, m2z);

  //print_fe(m2x);
  //print_fe(m2y);
  //printf("\n");
  ge_p2_dbl(&piop1, &piop2);
  ge_p1p1_to_p2(&piop2, &piop1);
  ge_p2_dbl(&piop1, &piop2);
  ge_p1p1_to_p2(&piop2, &piop1);
  ge_p2_dbl(&piop1, &piop2);
  ge_p1p1_to_p2(&piop2, &piop1);

  by piopoint;
  ge_tobytes(piopoint, &piop2);
  by_print(piopoint);




  fe_invert(t1, t1);   // t1 = 1 / (uf + 1)
  fe_mul(y, t0, t1);   // y = (uf - 1) / (uf + 1)
  ge_p3 he;
  ge_p1p1 hc;
  ge_p2 hp;

  by point;
  fe_tobytes(point, y);
  ge_frombytes_vartime(&he, point);
  ge_p3_to_p2(&hp, &he);

  //fe m1z, m1x, m1y;
  //fe_invert(m1z, hp.Z);
  //fe_mul(m1x, hp.X, m1z);
  //fe_mul(m1y, hp.Y, m1z);
  //printf("\n");
  //print_fe(m1x);
  //print_fe(m1y);
  //printf("\n");
  ge_p2_dbl(&hc, &hp);
  ge_p1p1_to_p2(&hp, &hc);
  ge_p2_dbl(&hc, &hp);
  ge_p1p1_to_p2(&hp, &hc);
  ge_p2_dbl(&hc, &hp);
  ge_p1p1_to_p2(&hp, &hc);

  ge_tobytes(point, &hp);
  by_print(point);


}
