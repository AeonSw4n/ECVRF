#ifndef ECVRF_H_

#define ECVRF_H_

static void ECVRF_prove(uint8_t* pi, const uint8_t* y,
                const uint8_t* x, const uint8_t* alpha);
static void ECVRF_proof2hash(uint8_t* beta, const uint8_t* pi);
static uint8_t ECVRF_verify(const by y, const uint8_t* pi,
                  const uint8_t* alpha);


#endif /* ECVRF_H_
