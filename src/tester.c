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
  uint8_t *input = "244df3f602b12172a6956cae6c484b3d7cb70c2688a4c68412481225c2848e76";
  int negX = 0;
  int sqrM = 0;
  by in;
  by_fromstr(in, input);
  fe_frombytes(r, in);
  by test;
  fe_tobytes(test, r);
  by_print(test);

  fe A = {0};
  fe e1 = {0};
  fe e2 = {0};

  fe t0 = {0};
  fe t1 = {0};

  fe o0 = {0};

  fe u0 = {0};
  fe u1 = {0};
  fe w = {0};
  fe sqr = {0};
  fe_1(o0);

  elligator_fe_A(A);

  fe_sq(t0, r);         // r ** 2
  fe_add(t0, t0, t0);   // 2 * (r ** 2)
  fe_add(t0, o0, t0);  // 1 + 2 * (r ** 2)
  fe_invert(t0, t0);    // 1 / (1 + 2 * (r ** 2))


  fe_0(u0);
  fe_sub(u0, u0, A);      // -A
  fe_mul(u0, u0, t0);     // -A / (1 + 2 * (r ** 2))

  // 11. w = u * (u^2 + A*u + 1) mod p

  fe_sq(t0, u0);       // u ** 2
  fe_mul(t1, A, u0);    // A * u
  fe_add(t0, t0, t1);  // u**2 + A*u
  fe_add(t0, t0, o0); // u**2 + A*u + 1

  fe_mul(w, u0, t0);    // u * (u**2 + A*u + 1)

  fe_tobytes(in, w);
  by_print(in);
  // 12. e = Legendre symbol of w and p

  fe_sqr_legendre(sqr, e1, w);

  fe_copy(e2, e1);

  fe_0(u1);
  fe_sub(u1, u1, A);   // -A
  fe_sub(u1, u1, u0);   // -A - u
  //print_fe(u1);
  //print_fe(u0);
  print_fe(e1);
  e1[0] = (e1[0] + 1) & 2;
  unsigned int b = e1[0] >> 1;
  fe_cmov(u1, u0, b);  // swaps if b == 1 (e=1), making uf = u; otherwise uf = -A-u (e=-1)
  fe uf = {0};
  fe_copy(uf, u1);

  fe w2;
  fe_sq2(w2, r);
  fe_mul(w2, w2, w);
  fe_cmov(w2, w, b);

  //e1[0] = (e1[0] + 1) & 2;
//  unsigned int b = e1[0] >> 1;
//  fe_cmov(u0, u1, b);  // swaps if b == 1 (e=1), making uf = u; otherwise uf = -A-u (e=-1)
  //printf("\n");
  //print_fe(uf);
  fe f2254;
  fe f1o2;
  fe_1(f1o2);
  fe_add(f1o2, f1o2, f1o2);
  fe_sq(f2254, f1o2);
  for(int i=2; i<254;i++)
    fe_mul(f2254, f2254, f1o2);
  //fe_neg(f2254, f2254);
  print_fe(f2254);
  fe_invert(f1o2, f1o2);
  print_fe(f1o2);
  fe m1f;
  fe_1(m1f);
  fe_sub(f1o2, f1o2, m1f);
  print_fe(f1o2);
  fe o2;
  fe u252m2;
  by fb4;
  by fb5;
  fe_tobytes(fb4, sqr);
  fe_tobytes(fb5, w);
  printf("\n");
  by_print(fb4);
  by_print(fb5);
  fe_1(o2);
  fe_u252m2(u252m2);
  fe_neg(u252m2, u252m2);
  printf("u252m2: %d\n", fe_isnegative(u252m2));
  fe_mul(u252m2, u252m2, r);
  fe_cmov(u252m2, o2, b);
  fe_mul(sqr, sqr, u252m2);
  fe sqr2;
  fe_sq(sqr2, sqr);
  b = !(fe_parity(sqr2) ^ fe_parity(w2)); // 0 if they differ (good). 1 if they are the same (bad)
  fe_cmov(o2, sqrtm1, b);
  fe_mul(sqr, sqr, o2);
  printf("Parity of sqrt: %d, differ: %d\n", fe_parity(sqr), b);

  printf("is negative: %d\n", fe_isnegative(sqr));
//  if(fe_isnegative(sqr))
//    fe_neg(sqr, sqr);
  //if(fe_isnegative(sqr))
  //  fe_neg(sqr, sqr);
  by fb3;
  fe_tobytes(fb3, sqr);
  fe testsqr2;
  fe_mul(testsqr2, sqr, sqr);
  //printf("\n w test\n");
  by fb1, fb2;
  //print_fe(testsqr2);
  fe_tobytes(fb1, testsqr2);
  //print_fe(w2);
  fe_tobytes(fb2, w2);
  printf("--------\n");
  //by_print(fb1);
  //by_print(fb2);
  //by_print(fb3);
  //printf("* differ: %d\n", !b);
  //by_print(fb4);
  //printf("parity sqr^2: %d, parity w: %d\n", fe_parity(testsqr2), fe_parity(w2));


  //print_fe(e2
  e2[0] = -e2[0];
  fe_mul(sqr, sqr, e2);
  fe sqA;
  fe_sqAp2(sqA);
  fe_neg(sqA, sqA);
  printf("sqA neg: %d\n", fe_isnegative(sqA));
  fe_mul(sqA, sqA, uf);
  //if(fe_isnegative(sqr) ^ fe_isnegative(sqA))
  //  fe_neg(sqr, sqr);

  ge_p1p1 piop1;

  //print_fe(o0);
  //print_fe(uf);
  fe_sub(t0, uf, o0); // t0 = uf - 1
  fe_copy(piop1.Y, t0);
  fe_add(t1, uf, o0); // t1 = uf + 1
  fe_copy(piop1.T, t1);
  printf("y neg: %d %d ", fe_parity(t0), fe_parity(t1));

  fe_invert(t1, t1);   // t1 = 1 / (uf + 1)
  printf("%d\n", fe_parity(t1));
  fe y = {0};
  fe_mul(y, t0, t1);   // y = (uf - 1) / (uf + 1)
  printf("y2 neg: %d\n", fe_parity(y));
  by testy;
  fe_tobytes(testy, y);
  printf("Edwards Y MSB: %d", testy[31]>>6);


  ge_p3 he;
  ge_p1p1 hc;
  ge_p2 hp;

  by point;
  fe_tobytes(point, y);
  ge_frombytes_vartime(&he, point);
  ge_p3_to_p2(&hp, &he);
  printf("\n");
  //print_fe(hp.X);
  //print_fe(hp.Y);
  //print_fe(hp.Z);
  fe m1z, m1x, m1y;
  fe_invert(m1z, hp.Z);
  fe_mul(m1x, hp.X, m1z);
  fe_mul(m1y, hp.Y, m1z);
  printf("\n");
  print_fe(m1x);
  print_fe(m1y);
  printf("\n");
  ge_p2_dbl(&hc, &hp);
  ge_p1p1_to_p2(&hp, &hc);
  ge_p2_dbl(&hc, &hp);
  ge_p1p1_to_p2(&hp, &hc);
  ge_p2_dbl(&hc, &hp);
  ge_p1p1_to_p2(&hp, &hc);

  ge_tobytes(point, &hp);
  by_print(point);

  /*fe_copy(hp.X, uf);
  fe_copy(hp.Y, y);
  fe_1(hp.Z);
*/

  // 16. out_point = (uf, y) ^ 8




  fe_copy(piop1.X, sqA);
  fe_copy(piop1.Z, sqr);
  ge_p2 piop2;
  ge_p1p1_to_p2(&piop2, &piop1);
  printf("\n");
  //print_fe(piop2.X);
  //print_fe(piop2.Y);
  //print_fe(piop2.Z);
  fe m2z, m2x, m2y;
  //if(negX)
  //  fe_neg(&piop2.X, &piop2.X);
  fe_invert(m2z, piop2.Z);
  fe_mul(m2x, piop2.X, m2z);
  fe_mul(m2y, piop2.Y, m2z);
  printf("%d %d %d | %d %d %d\n", fe_parity(piop2.Z), fe_parity(piop2.Y), fe_parity(piop2.X), fe_parity(m2x), fe_parity(m2y), fe_parity(m2z));
  print_fe(m2x);
  print_fe(m2y);
  printf("\n");
  ge_p2_dbl(&piop1, &piop2);
  ge_p1p1_to_p2(&piop2, &piop1);
  ge_p2_dbl(&piop1, &piop2);
  ge_p1p1_to_p2(&piop2, &piop1);
  ge_p2_dbl(&piop1, &piop2);
  ge_p1p1_to_p2(&piop2, &piop1);
  by piopoint;
  ge_tobytes(piopoint, &piop2);
  fe_invert(&piop2.Z, &piop2.Z);
  fe_mul(piop2.X, piop2.X, piop2.Z);
  fe_mul(piop2.Y, piop2.Y, piop2.Z);
  fe_tobytes(piopoint, piop2.Y);
  printf("(%d %d)\n", fe_parity(piop2.X), fe_parity(piop2.Y));
  //if(piopoint[31] >> 7 != fe_parity(piop2.X))
  //  piopoint[31] ^= 1<<7;
  //printf("\n");
  by_print(piopoint);
}
