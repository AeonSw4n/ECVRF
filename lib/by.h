#ifndef ECVRFBY_H_

#define ECVRFBY_H_
#include "ec.h"

typedef uint8_t by[32];

static void by_fromstr(by s, const uint8_t* str);
static void by_fromstrbe(by s, const uint8_t* str);
static void by_fromstrc(by s, const uint8_t* str,
                      const uint8_t count);
static void by_fromint(by s, int32_t x);

#endif /* ECVRFBY_H_ */
