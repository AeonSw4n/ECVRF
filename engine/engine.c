#include <stdio.h>
#include <stdlib.h>
#include <openssl/engine.h>
#include "../lib/engine.h"
#include "../lib/env.h"
#define ELLIGATOR2_TESTS ROOT "/engine/tests.txt"

static const char *engine_id = "ecvrf";
static const char *engine_name = "OpenSSL engine implementing ECVRF!";

int ecvrf_init(ENGINE *e){
  uint8_t pi[80];
  const uint8_t x_raw[64] = "5735f28b16c6b0684e6787866ef296a0787f1b4760d8101e090a533e04df1015";
  by x;
  by_fromstr(x, x_raw);
  uint32_t alpha_len = 30;
  const uint8_t alpha_raw[60] = "f4599e95b49f18923bc5a7939e35980a062e2573445a79601964ed2753ba";
  uint8_t alpha[30];
  for(uint8_t i=0; i<alpha_len;i++){
    const char b[2] = {alpha_raw[2*i], alpha_raw[2*i+1]};
    uint32_t xc32;
    sscanf(b, "%2x", &xc32);
    alpha[i] = (uint8_t)(xc32&255);
  }
  ECVRF_prove(pi, x, alpha, alpha_len);
  /*
  FILE *fp;
  fp = fopen(ELLIGATOR2_TESTS, "r");
  if (fp == NULL)
  {
     perror("Error while opening the file.\n");
     exit(EXIT_FAILURE);
  }
  int counter = 0;
  int passed = 0;
  uint8_t line[256];
  while(fgets(line, 256, fp) != NULL && counter < 2000){
    by PK;
    by_fromstr(PK, line);

    fgets(line, 256, fp);
    uint32_t len;
    sscanf(line, "%d %s", &len, line);

    for(uint8_t i=0; i<len;i++){
      const char b[2] = {line[2*i], line[2*i+1]};
      uint32_t xc32;
      sscanf(b, "%2x", &xc32);
      line[i] = (uint8_t)(xc32&255);
    }
    by out;
    elligator2_ed25519(out, line, len, PK);

    by K;
    fgets(line, 256, fp);
    by_fromstr(K, line);

    uint8_t A = by_cmp(out, K);
    if(A == 0){
      printf("--------------------\n");
      printf("ERROR ON TEST: %d\n", counter);
      printf("RECEIVED: ");
      by_print(out);
      printf("EXPECTED: ");
      by_print(K);
      printf("--------------------\n");
    }
    else{
      printf("---TEST %d PASSED---\n", counter);
      passed++;
    }
    fgets(line, 256, fp);
    counter++;
  }
  printf("\n\n Passed/Counter: %d/%d\n", passed, counter);
  */
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
