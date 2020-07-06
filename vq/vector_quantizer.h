/*
 * base class for computing cluster centroids
 */
#ifndef MEANS_QUANTIZER_H
#define MEANS_QUANTIZER_H

#include "distance_measures.h"

template<class T, unsigned char D>
class VectorQuantizer
{
public:
  const unsigned int    K;
  const float        delta;  // termination parameter (epsilon?)
  const unsigned int    maxIterations;

protected:
  DistanceMeasure<T, D>&  distance;

  T    *mean;

  // stats
  unsigned int  stat_iteration_changes;
  unsigned int  stat_ks_adjusted;
  unsigned int  stat_dist_calc;
  unsigned int  stat_iterations;

  void initialiseMeans ()
  {
    if (mean != NULL)
      delete[] mean;

    mean = new T[K*D];
  }

  void resetIterations ()
  {
    stat_iterations = 0;
  }

  void zeroStats ()
  {
    stat_iteration_changes = 0;
    stat_ks_adjusted = 0;
    stat_dist_calc = 0;
  }

  bool maxItersReached () const
  {
    return ((maxIterations > 0) && (stat_iterations == maxIterations ));
  }

  // select k random codewords from the sampledata
  void selectRandomMeans (const T *samples, const unsigned int N)
  {
    unsigned int  i, d, m;

    for (i=0; i<K; i++)
    {
      srand (i);
      m = static_cast<unsigned int>((static_cast<float>(rand ())/static_cast<float>(RAND_MAX))*N);

      for (d=0; d<D; d++)
        mean[(i*D)+d] = samples[(m*D)+d] + ((float)rand()/(float)RAND_MAX);
    }
  }
  virtual void adjustMeans (const T *samples, const unsigned int N, T* means) = 0;

  const float adjustDistance (const T* oldMeans, const T* newMeans)
  {
    unsigned int  i, j;
    float      l, d, e;

    l = 0.0f;

    for (i=0; i<K; i++)
    {
      d = 0.0f;

      for (j=0; j<D; j++)
      {
        e = (oldMeans[(i*D)+j] - newMeans[(i*D)+j]);
        d += e*e;
      }

      l += sqrt(d);
    }

    return l;
  }
public:
  VectorQuantizer (const unsigned int clusterCount, DistanceMeasure<T, D>& distanceMeasure, const unsigned int maximumIterations = 0, const float delta = 0.01f)
    : K(clusterCount), distance(distanceMeasure), mean(NULL), maxIterations(maximumIterations), delta(delta)
  {
    zeroStats ();
  }

  virtual ~VectorQuantizer ()
  {
    if (mean != NULL)
      delete[] mean;
  }

  // clustering
  virtual void computeMeans (const T* samples, const unsigned int N, const unsigned int maxIters = 0) = 0;
  // assignment info
  virtual unsigned int getCluster (const unsigned int n) const = 0;      // get cluster for sample index
  virtual float getClusterProbability (const unsigned int n) const = 0;    // get probability for sample index

  T* getClusterMean (const unsigned int n) const
  {
    return &mean[(n*D)];
  }
  // data output
  //void writeMeans (std::ostream& os)
  //{
  //  unsigned int  tK, tD;
  //  tK = K;
  //  tD = D;

  //  os.write ((char*)&tK, sizeof (unsigned int));
  //  os.write ((char*)&tD, sizeof (unsigned int));
  //  os.write ((char*)mean, sizeof (float) * K * D);
  //}
};
#endif
