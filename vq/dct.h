/*
 * 1D DCT of the input data
 */
#define _USE_MATH_DEFINES
#include <cmath>

class DiscreteCosineTransform
{
public:
  static float* t (const float *x, const unsigned int N, const unsigned char D)
  {
    float  *c = new float[N*D];
    const float PI_N = M_PI/(float)N;

    for (unsigned char d=0; d<D; d++)
    {
      for (unsigned int k=d; k<N*D; k+=D)
      {
        c[k] = 0;

        for (unsigned int n=d; n<N*D; n+=D)
        {
          c[k] += x[n] * cos (PI_N*(n+1.0/2.0)*k);
        }
      }
    }

    return c;
  }
};
