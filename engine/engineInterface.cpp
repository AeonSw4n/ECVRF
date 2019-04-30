#include <openssl/engine.h>
#include <stdio.h>
#include <string.h>
#include "../lib/env.h"
#define ENGINE_SO ROOT "/engine/engine.so"

int main(int argc, const char* argv[] ) {
    OpenSSL_add_all_algorithms();

    ERR_load_crypto_strings();

    ENGINE_load_dynamic();
    ENGINE_load_builtin_engines();

    ENGINE *ecvrf = ENGINE_by_id("dynamic");
    ENGINE_ctrl_cmd_string(ecvrf, "SO_PATH", ENGINE_SO, 0);
    ENGINE_ctrl_cmd_string(ecvrf, "ID", "ecvrf", 0);
    ENGINE_ctrl_cmd_string(ecvrf, "LOAD", NULL, 0);

    if( ecvrf == NULL )
    {
        printf("Could not Load Engine!\n");
        exit(1);
    }
    printf("Engine successfully loaded\n");
    int init_res = ENGINE_init(ecvrf);
    printf("Engine name: %s init result : %d \n", ENGINE_get_id(ecvrf), init_res);
    return 0;
}
