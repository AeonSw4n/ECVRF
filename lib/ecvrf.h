#ifndef ECVRF_H_

#define ECVRF_H_
#include "ec.h"
#include "by.h"

static void ECVRF_prove(double *t, uint8_t* pi, const uint8_t* SK,
                    const uint8_t* alpha, const size_t alpha_len);
static void ECVRF_proof2hash(uint8_t* beta, const uint8_t* pi);
static int ECVRF_verify(const uint8_t *y,
                        const uint8_t *pi, const uint8_t *alpha, const size_t alpha_len);
static void double_add_scalar_mult(uint8_t out[32],
                                        const uint8_t scalar[32],
                                        const uint8_t point[32]);
#endif /* ECVRF_H_ */
