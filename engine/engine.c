#include <stdio.h>
#include <openssl/engine.h>
#include "../lib/engine.h"

static const char *engine_id = "ecvrf";
static const char *engine_name = "OpenSSL engine implementing ECVRF!";

int ecvrf_init(ENGINE *e){
  //test_fe_legendre();
  //test_elligator_2();
  uint8_t data[2] = {175, 130};
  uint8_t *PK  = "fc51cd8e6218a1a38da47ed00230f0580816ed13ba3303ac5deb911548908025";
  by in;
  by_fromstr(in, PK);
  //by_print(in);
  by out;
  elligator2_ed25519(data, 2, in, out);
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
