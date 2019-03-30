#include <stdio.h>
#include "../lib/ec.h"
#include "../lib/by.h"


/*
      F - finite field

      2n = 32 - length, in octets, of a field element in F; must be even

      E - elliptic curve (EC) defined over F

      m = 32 - length, in octets, of an EC point encoded as an octet string

      G - subgroup of E of large prime order

      q - prime order of group G; must be less than 2^{2n}

      cofactor = 8 - number of points on E divided by q

      g - generator of group G

      Hash - cryptographic hash function

      hLen - output length in octets of Hash; must be at least 2n

      suite = 0x03 - a single nonzero octet specifying the ECVRF ciphersuite,
      which determines the above options

*/

static void ECVRF_prove(uint8_t* pi, const uint8_t* y,
                const uint8_t* x, const uint8_t* alpha)
{
  /*
   1.  h = ECVRF_hash_to_curve(suite, y, alpha)

   2.  h1 = EC2OSP(h)

   3.  gamma = h^x

   4.  k = ECVRF_nonce_generation(SK, h1)

   5.  c = ECVRF_hash_points(h, gamma, g^k, h^k)

   6.  s = (k + c*x) mod q (where * denotes integer multiplication)

   7.  pi = EC2OSP(gamma) || I2OSP(c, n) || I2OSP(s, 2n)

   8.  Output pi
 */
}

static void ECVRF_proof2hash(uint8_t* beta, const uint8_t* pi)
{
  /*
    1.  D = ECVRF_decode_proof(pi)

    2.  If D is "INVALID", output "INVALID" and stop

    3.  (gamma, c, s) = D

    4.  three = 0x03 = I2OSP(3, 1), a single octet with value 3

    5.  preBeta = Hash(suite || three || EC2OSP(gamma^cofactor))

    6.  beta = first 2n octets of preBeta

    7.  Output beta
  */
}

static uint8_t ECVRF_verify(const by y, const uint8_t* pi,
                  const uint8_t* alpha)
{
  /*
    1.  D = ECVRF_decode_proof(pi)

    2.  If D is "INVALID", output "INVALID" and stop

    3.  (gamma, c, s) = D

    4.  u = g^s / y^c (where / denotes EC point subtraction, i.e. the
       group operation applied to g^s and the inverse of y^c)

    5.  h = ECVRF_hash_to_curve(suite, y, alpha)

    6.  v = h^s / gamma^c (where / again denotes EC point subtraction)

    7.  c' = ECVRF_hash_points(h, gamma, u, v)

    8.  If c and c' are equal, output "VALID"; else output "INVALID"
  */
  return 0;
}
