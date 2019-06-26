#ifndef ECVRFELLIGATOR2_H_

#define ECVRFELLIGATOR2_H_
#include "ec.h"
#include "by.h"
#include "tester.h"

static void elligator_fe_A(fe out);
static void fe_legendre(fe out, const fe z);
static void fe_sqr_legendre(fe sqrt, fe e, const fe z);
static unsigned int fe_parity(const fe h);
static void fe_sqAp2(fe out);
static void fe_u252m2(fe out);
static void fe_legendre_fromby(fe a, by b);
static uint32_t fe_legendre_ifromby(by b);
static void elligator2_ed25519(ge_p3 *out_point, by out_point_bytes, by out_point_mont, const uint8_t *data,
                  size_t size, by public_key);
#endif /* ECVRFELLIGATOR2_H_ */
