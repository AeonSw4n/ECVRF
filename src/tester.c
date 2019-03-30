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
                          /*"c7b5d6239e52a473a2b57a92825e0e5de4656e349bb198de5afd6a76e5a07066",
                          "7ff6d8b773bfbae57b2ab9d49f9d3cb7d9af40a03d3ed3c6beaaf2d486b1fe6e",
                          "2ceaa2c2ff3028c34f9fbe076ff99520b925f18d652285b4daad5ccc467e523b",*/
                          "6670a0e5766afd5ade98b19b346e65e45d0e5e82927ab5a273a4529e23d6b5c7",
                          "6efeb186d4f2aabec6d33e3da040afd9b73c9d9fd4b92a7be5babf73b7d8f67f",
                          "3b527e46cc5caddab48522658df125b92095f96f07be9f4fc32830ffc2a2ea2c",
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
