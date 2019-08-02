#ifndef ECVRF_MERGED_H_

#define ECVRF_MERGED_H_
#include "ec.h"
#include "by.h"
#include "tester.h"

#include "../lib/ecvrf_merged.h"
#include <stdio.h>
#include <time.h>


static void fe_A(fe out);

static void fe_2A(fe out);

static void fe_A_cubed(fe out);

static void fe_sqrt_neg_A_plus_2(fe out);

static void fe_fourth_root_of_neg_4(fe out);

static unsigned int fe_ispositive(const fe h);

#if defined (BASE_2_51_IMPLEMENTED)
  static void fe_mul121666(fe h, const fe f);
  static void fe_cswap(fe f, fe g, unsigned int b);
#endif

static void fe_fraction_sqrt_and_legendre(fe sqrt_helper, fe e, const fe w);

static void cswap(uint8_t *t0, uint8_t *t1, unsigned int b);

static void ge_p2_to_p3(ge_p3 *r, const ge_p2 *p);

static void ge_p3_copy(ge_p3 *r, const ge_p3 *p);

static void ge_p3_cmov(ge_p3 *p, const ge_p3 *q, unsigned int b);

static void ge_p3_cswap(ge_p3 *p, ge_p3 *q, unsigned int b);

static void ge_p3_neg(ge_p3 *r, const ge_p3 *p);

static int ge_p3_frombytes(ge_p3 *h, const uint8_t s[32]);

static void ge_p3_two_tobytes(uint8_t out1[32], uint8_t out2[32],
                                const ge_p3 *p, const ge_p3 *q);

static void ge_p3_three_tobytes(uint8_t out1[32], uint8_t out2[32], uint8_t out3[32],
                                const ge_p3 *p, const ge_p3 *q, const ge_p3 *r);

static void ge_p3_two_to_montgomery(uint8_t u1[32], uint8_t v1[32], uint8_t u2[32], uint8_t v2[32],
                                              const ge_p3 *p, const ge_p3 *q);

static void ge_p1p1_0(ge_p1p1 *h);

static void ge_dedicated_add(ge_p3 *r, const ge_p3 *p, const ge_p3 *q);

static void double_scalar_fixed_point_mult(ge_p3 *p, ge_p3 *q, const ge_p3 *H,
                                  const uint8_t scalar1[32], const uint8_t scalar2[32]);

static void montgomery_p2_to_ge_p3(ge_p3 *r, const fe u, const fe v, const fe z);

static void montgomery_recover_point(fe U, fe V, fe Z, const uint8_t ub[32], const uint8_t vb[32],
                                const fe U1, const fe Z1, const fe U2, const fe Z2);

static void montgomery_ladder(fe out_u2, fe out_z2, fe out_u3, fe out_z3,
                                       const uint8_t * scalar, const size_t scalar_len, const uint8_t in_u[32]);

static void montgomery_ladder_scalar_mult(ge_p3 *r, const uint8_t *scalar, const size_t scalar_len,
                                                  const uint8_t u[32], const uint8_t v[32]);

static void montgomery_double_scalar_mult_difference(ge_p3 *r, const ge_p3 *p, const ge_p3 *q,
                                                  const uint8_t *scalar1, const size_t scalar1_len,
                                                  const uint8_t *scalar2, const size_t scalar2_len);

static void ECVRF_hash_to_curve_elligator2_25519(ge_p3 *out_point, uint8_t out_point_bytes[32], const fe r);

static void ECVRF_alpha_to_curve(ge_p3 *out_point, uint8_t out_point_bytes[32],
                  const uint8_t public_key[32], const uint8_t *alpha, size_t alpha_len);

static int ECVRF_decode_proof(ge_p3 *Gamma, uint8_t c[16], uint8_t s[32], const uint8_t pi[80]);

static int ECVRF_decode_proof_vartime(ge_p3 *Gamma, uint8_t c[16], uint8_t s[32], const uint8_t pi[80]);

static int ECVRF_proof_to_hash(uint8_t *beta, const uint8_t *pi);

static int ECVRF_proof_to_hash_vartime(uint8_t *beta, const uint8_t *pi);

static void ECVRF_prove(uint8_t pi[80], const uint8_t SK_bytes[32],
                    const uint8_t* alpha, const size_t alpha_len);

static int ECVRF_verify(const uint8_t *Y_bytes, const uint8_t *pi,
                  const uint8_t *alpha, const size_t alpha_len);

#endif
