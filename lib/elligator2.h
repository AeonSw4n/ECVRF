#ifndef ECVRFELLIGATOR2_H_

#define ECVRFELLIGATOR2_H_
#include "ec.h"
#include "by.h"

static void elligator_fe_A(fe out);
static void fe_legendre(fe out, const fe z);
static void fe_legendre_fromby(fe a, by b);
static uint32_t fe_legendre_ifromby(by b);
static void elligator2_ed25519(const uint8_t *data, size_t size,
                       const uint8_t public_key[32],
                       uint8_t out_point[32]);

#endif /* ECVRFELLIGATOR2_H_ */
