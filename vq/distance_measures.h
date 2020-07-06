/*
 * bunch of distance measures for clustering in different
 * spaces.
 * Parameters are:
 *   T : primitive type
 *   D : arity of datum
 */
#ifndef DISTANCE_MEASURE_H
#define DISTANCE_MEASURE_H
#include <math.h>

template<class T, unsigned char D>
class DistanceMeasure
{
public:
  virtual float operator()(const T* a, const T *b) const = 0;
};

template<class T, unsigned char D>
class EuclideanDistance : public DistanceMeasure<T, D>
{
public:
  float operator ()(const T* a, const T* b) const
  {
    float  sum = 0;

    for (unsigned int i=0;i<D; i++)
    {
      sum += (a[i]-b[i]) * (a[i]-b[i]);
    }

    return static_cast<float>(sqrt (sum));
  }
};

template<class T, unsigned char D>
class MinkowskiDistance : public DistanceMeasure<T, D>
{
public:
  // FIXME: parameterise lambda
  float operator ()(const T* a, const T* b)
  {
    T  sum = 0;
    const float  ex = 1.0f;

    for (unsigned int i=0;i<D; i++)
    {
      sum += pow ((float)abs (a[i]-b[i]), 1.0f);
    }

    return static_cast<float>(pow (sum, 1.0f/1.0f));
  }
};

template<class T, unsigned char D>
class InnerProduct : public DistanceMeasure<T, D>
{
public:
  float operator ()(const T* a, const T* b) const
  {
    T  sum = 0;

    for (unsigned int i=0; i<D; i++)
    {
      sum += a[i]*b[i];
    }

    return sum;
  }
};
template<class T, unsigned char D>
class ManhattanDistance : public DistanceMeasure<T, D>
{
public:
  float operator ()(const T* a, const T* b) const
  {
    T  sum = 0;

    for (unsigned int i=0; i<D; i++)
    {
      sum += abs (a[i]-b[i]);
    }

    return sum;
  }
};

template<class T, unsigned char D>
class ChebyshevDistance : public DistanceMeasure<T, D>
{
public:
  float operator ()(const T* a, const T* b) const
  {
    T max = abs (a[0] - b[0]);

    for (unsigned int i=1; i<D; i++)
    {
      if (abs (a[i] - b[i]) > max)
        max = abs (a[i] - b[i]);
    }

    return max;
  }
};
template<class T, unsigned char D>
class SphericalDistance : public DistanceMeasure<T, D>
{
public:
  float operator ()(const T* a, const T* b) const
  {
    float  l = 0.0f;
    T    c[D];
    int    i;

    for (i=0; i<D; i++)
    {
      c[i] = acos (a[i]*b[i]);
    }

    for (i=0; i<D; i++)
    {
      l += c[i];
    }

    return l;
  }
};
#endif
