#ifndef ECVRF_MERGED_H_

#define ECVRF_MERGED_H_
#include "ec.h"
#include "by.h"
#include "tester.h"

static void fe_A(fe out);

static void fe_A3(fe out);

static void fe_sqrt_A_plus_2(fe out);

static void fe_sqrt_2(fe out);

static unsigned int fe_ispositive(const fe h);

static void fe_fraction_sqrt_and_legendre(fe sqrt, fe e, const fe w);

static void ECVRF_hash_to_curve_elligator2_25519(ge_p3 *out_point, uint8_t out_point_bytes[32],
                                                    const fe r);

static void ECVRF_alpha_to_curve(ge_p3 *out_point, uint8_t out_point_bytes[32],
                  const uint8_t public_key[32], const uint8_t *alpha, size_t alpha_len);
#endif
