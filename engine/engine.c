#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <openssl/engine.h>
#include "../lib/engine.h"
#include "../lib/env.h"

//#define DEBUG
#define TESTER

#define ELLIGATOR2_TESTS ROOT "/tests/tests_elligator2.txt"
#define PROVE_TESTS ROOT "/tests/tests_prove.txt"

static const char *engine_id = "ecvrf";
static const char *engine_name = "OpenSSL engine implementing ECVRF!";

int ecvrf_init(ENGINE *e)
{
  double t, t_total, t_average;
  t_average = 0.0;

 #ifndef TESTER
  uint8_t pi[80];
  const uint8_t x_raw[64] = "2b7794e737d0dffc7221360462617b256ddbec98018a2781152786f1560c4310";
  by x;
  by_fromstr(x, x_raw);
  uint32_t alpha_len = 26;
  const uint8_t alpha_raw[52] = "919687cfd226964cc3dc54aa62ea5793f17b567bd28c6a7f8c52";
  uint8_t alpha[26];
  for(uint8_t i=0; i<alpha_len;i++){
    const char b[2] = {alpha_raw[2*i], alpha_raw[2*i+1]};
    uint32_t xc32;
    sscanf(b, "%2x", &xc32);
    alpha[i] = (uint8_t)(xc32&255);
  }
  ECVRF_prove(pi, x, alpha, alpha_len);
  for(int i=0; i<80; i++){
    if(pi[i] < 16)
      printf("0%1x", pi[i]);
    else
      printf("%2x", pi[i]);
  }
  printf("\n\n");
  uint8_t beta[64];
  ECVRF_proof_to_hash(beta, pi);
  for(int i=0;i<64;i++){
    if(beta[i] < 16)
      printf("0%1x", beta[i]);
    else
      printf("%2x", beta[i]);
  }
  printf("\n");

  uint8_t hash[SHA512_DIGEST_LENGTH] = {0};
  SHA512_CTX hash_ctx;
  SHA512_Init(&hash_ctx);
  SHA512_Update(&hash_ctx, x, 32);
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

  //const uint8_t y_raw[64] = "d2c358b9b7937168ab115c6281272eed43d8a40cbaf77f3b09e72a815eae29a2";
  //by y;
  //by_fromstr(y, y_raw);
  int eq = ECVRF_verify(y, pi, alpha, alpha_len);
  printf("ECVRF_verify returns: %d\n", eq);
 #endif


  for(int rep=0; rep<1; rep++){
    t_total = 0.0;
    #ifdef TESTER
      FILE *fp;
      fp = fopen(PROVE_TESTS, "r");
      if (fp == NULL)
      {
         perror("Error while opening the file.\n");
         exit(EXIT_FAILURE);
      }
      int counter = 0;
      int passed = 0;
      int passed_verfity = 0;
      int passed_verfity_false = 0;
      uint8_t line[256];
      while(fgets(line, 256, fp) != NULL && counter < 10000){
        by SK;
        by_fromstr(SK, line);

        fgets(line, 256, fp);
        uint32_t len;
        sscanf(line, "%d %s", &len, line);

        uint8_t i;
        for(i=0; i<len;i++){
          const char b[2] = {line[2*i], line[2*i+1]};
          uint32_t xc32;
          sscanf(b, "%2x", &xc32);
          line[i] = (uint8_t)(xc32&255);
        }
        uint8_t pi[80];

        ECVRF_prove(pi, SK, line, len);
        t_total += t;

        uint8_t alpha[256];
        memcpy(alpha, line, len);

        fgets(line, 256, fp);
        for(i=0; i<80;i++){
          const char b[2] = {line[2*i], line[2*i+1]};
          uint32_t xc32;
          sscanf(b, "%2x", &xc32);
          line[i] = (uint8_t)(xc32&255);
        }
        uint8_t A = 1;
        for(i=0; i<80; i++){
          if(line[i] != pi[i]){
            A = 0;
            break;
          }
        }
        if(A)
          passed++;

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

        A = ECVRF_verify(y, pi, alpha, len);
        if(A)
          passed_verfity++;

        y[0] += 1;
        A = ECVRF_verify(y, pi, alpha, len);
        if(A)
          passed_verfity_false++;

      #ifdef DEBUG
        if(A == 0){
          printf("--------------------\n");
          printf("ERROR ON TEST: %d\n", counter);
          printf("RECEIVED: ");
          for(i=0; i<80; i++)
            printf("%2x", pi[i]);
          printf("\n");
          printf("EXPECTED: ");
          for(i=0; i<80; i++)
            printf("%2x", line[i]);
          printf("\n");
          printf("--------------------\n");
        }
        else{
          printf("---TEST %d PASSED---\n", counter);
        }
      #endif

        fgets(line, 256, fp);
        counter++;
      }
      printf("\n\n Passed Prove: %d/%d  ||  Passed Verify: %d/%d\n", passed, counter, passed_verfity-passed_verfity_false, counter);
    #endif
    t_average = ((double)rep*t_average + t_total)/((double)rep + 1.0);
  }
  double cpu_time_used = t_average/(double)CLOCKS_PER_SEC;
  printf("time %f seconds\n", cpu_time_used);
  return 1;

 end:
  return 0;
}


static int bind(ENGINE *e, const char *id)
{
  int ret = 0;
  if (!ENGINE_set_id(e, engine_id)) {
    fprintf(stderr, "ENGINE_set_id failed\n");
    goto end;
  }
  if (!ENGINE_set_name(e, engine_name)) {
    printf("ENGINE_set_name failed\n");
    goto end;
  }
  if (!ENGINE_set_init_function(e, ecvrf_init)){
    printf("ENGINE_set_init_function failed\n");
    goto end;
  }
  ecvrf_init(e);

  ret = 1;
 end:
  return ret;
}

IMPLEMENT_DYNAMIC_BIND_FN(bind)
IMPLEMENT_DYNAMIC_CHECK_FN()
