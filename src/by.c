#include <stdio.h>
#include "../lib/by.h"


static void by_fromstr(by s, const uint8_t* str)
{
  memset(s, 0, 32);
  for(uint32_t i = 0; i<32; i++){
    const char b[2] = {str[2*i], str[2*i+1]};
    uint32_t xc32;
    sscanf(b, "%2x", &xc32);
    s[i] = (uint8_t)(xc32&255);
  }
}

static void by_fromstrbe(by s, const uint8_t* str)
{
  memset(s, 0, 32);
  for(uint32_t i = 0; i<32; i++){
    const char b[2] = {str[2*i], str[2*i+1]};
    uint32_t xc32;
    sscanf(b, "%2x", &xc32);
    s[31-i] = (uint8_t)(xc32&255);
  }
}

static void by_fromstrc(by s, const uint8_t* str,
                      const uint8_t count)
{
  uint8_t str2[64];
  memset(str2, '0', 64);

  for(uint8_t i=0; i<count;i++){
    str2[i+64-count] = str[i];
  }
  str = str2;

  by_fromstrbe(s, str);

}

static void by_fromint(by s, int32_t x)
{
  memset(s, 0, 32);
  for(uint32_t i=0;i<4;i++){
    s[i] = x&255;
    x >>= 8;
  }

}

static uint8_t by_cmp(by x, by y)
{
  for(uint32_t i=0; i<32; i++)
    if(x[i] != y[i])
      return 0;
  return 1;
}
