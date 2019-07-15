#ifndef ECVRF_H_

#define ECVRF_H_
#include "ec.h"
#include "by.h"


static void cswap(uint8_t *t0, uint8_t *t1, unsigned int b)

static void fe_mul121666(fe h, fe f)

static void fe_cswap(fe f, fe g, unsigned int b)

static void ge_p3_copy(ge_p3 *r, ge_p3 *p)

static void ge_p3_cmov(ge_p3 *p, ge_p3 *q, unsigned int b)

static void ge_p3_cswap(ge_p3 *p, ge_p3 *q, unsigned int b)

static void ge_dedicated_add(ge_p3 *r, ge_p3 *p, ge_p3 *q)

static int ge_p3_frombytes(ge_p3 *h, const uint8_t *s)

static void ge_p3_merge_two_tobytes(by out1, by out2, ge_p3 *p, ge_p3 *q)

static void double_scalar_fixed_point_mult_xor(uint8_t out1[32], uint8_t out2[32], ge_p3 *H,
                                              const uint8_t scalar1[32], const uint8_t scalar2[32])

static void double_scalar_fixed_point_mult(uint8_t out1[32], uint8_t out2[32], ge_p3 *H,
                                              const uint8_t scalar1[32], const uint8_t scalar2[32])

static void double_scalar_fixed_point_mult_plain(uint8_t out1[32], uint8_t out2[32], ge_p3 *H,
                                              const uint8_t scalar1[32], const uint8_t scalar2[32])

static void montgomery_ladder(fe out_x2[32], fe out_z2[32], fe out_x3[32], fe out_z3[32],
                                       const uint8_t scalar[32], const uint8_t in_x[32])

static void montgomery_ladder_to_edwards(ge_p3 *out, const uint8_t scalar[32],
                                       const uint8_t in_x[32], const ge_p3 *in_P)

static int ECVRF_decode_proof(ge_p3 *Gamma, uint8_t *c, uint8_t *s,
                            const uint8_t* pi)

static int ECVRF_proof_to_hash(uint8_t *beta, const uint8_t *pi)

static int ECVRF_proof_to_hash_vartime(uint8_t *beta, const uint8_t *pi)

static void ECVRF_prove(double *t, uint8_t* pi, const uint8_t* SK,
                    const uint8_t* alpha, const size_t alpha_len)

static int ECVRF_verify(const uint8_t *y, const uint8_t *pi,
                  const uint8_t *alpha, const size_t alpha_len)


#endif /* ECVRF_H_ */
