/*
 * K-Means gives a hard cluster membership classification
 * Optmized method using triangle inequality is implemented below
 */
#ifndef KMEANSQUANTIZER_H
#define KMEANSQUANTIZER_H
#include <limits>
#include <time.h>

#include "vector_quantizer.h"

template<class T, unsigned char D>
class KMeansQuantizer : public VectorQuantizer<T, D>
{
protected:
  unsigned short  *map;    // lookup from samples id to cluster id
  unsigned int  *cnt;    // count of samples in the cluster (TEMP)

  void setTempTables (const unsigned int N)
  {
    deleteTempTables ();
    map = new unsigned short[N];  // map of samples -> cluster
    cnt = new unsigned int[K];    // count of samples in cluster
  }

  void deleteTempTables ()
  {
    if (map != NULL)
    {
      delete[] map;
      map = NULL;
    }

    if (cnt != NULL)
    {
      delete[] cnt;
      cnt = NULL;
    }
  }

  // calculate new centroids after clustering the data
  void adjustMeans (const T *samples, const unsigned int N, T* means)
  {
    unsigned int  n;
    unsigned short  k;
    unsigned char  d;

    // reset the means
    // FIXME: remove the memset for T datatypes (could be a class with a vtable)
    memset (means, 0, K*D*sizeof(T));
    memset (cnt, 0, K*sizeof(unsigned int));

    // get the averages
    for (n=0; n<N; n++)
    {
      for (d=0; d<D; d++)
        means[(map[n]*D)+d] += samples[(n*D)+d];

      ++cnt[map[n]];
    }

    // update the means
    for (k=0; k<K; k++)
    {
      if (cnt[k] == 0)
        continue;

      for (d=0; d<D; d++)
        means[(k*D)+d] /= static_cast<float>(cnt[k]);
    }
  }

public:
  KMeansQuantizer (unsigned int K, DistanceMeasure<T, D>& distanceMeasure, const unsigned int maxIters=0)
    : VectorQuantizer (K, distanceMeasure, maxIters), map(NULL), cnt(NULL)
  {}

  virtual ~KMeansQuantizer ()
  {
    deleteTempTables ();
  }

  virtual void computeMeans (const T* samples, const unsigned int N, const unsigned int maxIters = 0)
  {
    float      dist, min, uk = 10.0f;
    const float    *s;
    unsigned int  i, k, c;
    float      *tmpMeans = new float[K*D];

    initialiseMeans ();        // set up the mean array
    selectRandomMeans (samples, N);  // select the intial means randomly
    setTempTables (N);        // setup the map and count tables

    resetIterations ();

    // converge
    while (uk > delta)
    {
      if (maxItersReached ())
      {
        std::cout << "Max iters reached." << std::endl;
        break;
      }

      zeroStats ();

      std::cout << "ITER: " << stat_iterations++ << " {";

      std::cout << "a";
      for (i=0; i<N; i++)
      {
        // get the nearest cluster to the sample f
        s = &samples[i*D];

        // set the distance to zero to be the minimum
        c = 0;
        min = distance (s, &mean[0]);

        for (k=1; k<K; k++)
        {
          dist = distance (s, &mean[k*D]);

          if (dist < min)
          {
            min = dist;
            c = k;
          }
        }

        // update the assignment map
        if (map[i] != c)
        {
          map[i] = c;
          ++stat_iteration_changes;
        }
      }

      std::cout << "m";

      // adjust the means now that the assignments have changed
      memset (tmpMeans, 0, sizeof (float) * K * D);
      adjustMeans (samples, N, tmpMeans);
      uk = adjustDistance (mean, tmpMeans);
      memcpy (mean, tmpMeans, sizeof (float) * K * D);

      std::cout << "} " << stat_iteration_changes << " changes";
      std::cout << "\t" << uk << " > " << delta << std::endl;
    }

    std::cout << "\tclustered in " << stat_iterations << " iterations." << std::endl;

    delete[] tmpMeans;
  }

  // get the nearest cluster to a data value
  unsigned int getCluster (const unsigned int n) const
  {
#if 0
    float      d;
    unsigned int  mindex = 0;
    float      min = distance (f, &mean[0]);

    for (unsigned int i=1; i<K; i++)
    {
      d = distance (f, &mean[i*D]);

      if (d < min)
      {
        min = d;
        mindex = i;
      }
    }

    return mindex;
#else
    return map[n];
#endif
  }

  // assignment probabilities are binary in k-means
  float getClusterProbability (const unsigned int n) const
  {
    return 1.0f;
  }
};

#if 1
// KMeans accelerated with the triangle inequality measure
// Using the Triangle Inequality to Accelerate K-Means - Charles Elkan
// http://www-cse.ucsd.edu/~elkan/kmeansicml03.pdf
template<class T, unsigned char D>
class KMeansQuantizer_TriIneq : public KMeansQuantizer<T,D>
{
public:
  KMeansQuantizer_TriIneq (unsigned int K, DistanceMeasure<T, D>& distance, const unsigned int maxIters=0)
    : KMeansQuantizer<T,D>(K, distance, maxIters)
  {}

  void computeMeans (const T *samples, const unsigned int N, const unsigned int maxIters = 0)
  {
    float      *u = new float[N];      // upper bounds
    float      *sc = new float[K];      // half distance to nearest neighbouring cluster
    unsigned short  *s = new unsigned short[K];  // index of the nearest neighbour
    float      *nm = new float[K*D];    // newly calculated means
    // FIXME: change the two bool arrays to char arrays using bit flags
    bool      *r = new bool[N];      // indicates if a distance needs to be calcualted
    bool      *m = new bool[K];      // indicates if a cluster moved
    unsigned int  n;
    unsigned short  k, ki;
    unsigned char  d;
    unsigned short  min_k;
    float      min_d, dist;
    float      uk = 10.0f;

    // INITIALISATION
    initialiseMeans ();        // set up the mean array
    selectRandomMeans (samples, N);  // select the intial means randomly
    setTempTables (N);        // setup the map and count tables

    // assign each point to its closest cluster
    for (n=0; n<N; n++)
    {
      map[n] = 0; //getCluster (&samples[n*D]);    // clustering
      u[n] = std::numeric_limits<float>::max ();    // force calculation
      r[n] = true;
    }

    for (k=0; k<K; k++)
    {
      m[k] = true;
    }

    // CONVERGENCE
    while (uk > delta)
    {
      if (maxItersReached ())
        break;

      zeroStats ();

      std::cout << "ITER: " << stat_iterations++ << " {";

      // 1. For all centres c and c' compute d(c,c').
      // For all centers c, compute s(c) = 1/2 min_{c'!=c}d(c,c')
      std::cout << "h";
      for (k=0; k<K; k++)
      {
        min_d = std::numeric_limits<float>::max ();
        min_k = 0;

        for (ki=0; ki<K; ki++)
        {
          if (ki == k)
            continue;

          dist = distance (&mean[ki*D], &mean[k*D]);
          ++stat_dist_calc;

          if (dist < min_d)
          {
            min_d = dist;
            min_k = ki;
          }
        }

        // note the nearest cluster and half the distance
        s[k] = min_k;
        sc[k] = 0.5f * min_d;
      }

      // cluster the data
      std::cout << "c";
      for (n=0; n<N; n++)
      {
        // 2. Identify all points x such that u(x) <= s(c(x)).
        // i.e. if the upper bound of the sample is greater than the
        // half distance to the nearest cluster, keep going
        if (u[n] < sc[map[n]])
          continue;

        // 3.a if r(x) then compute d(x, c(x)) and assign r(x) = false.
        // otherwise d(x,c(x)) = u(x)
        if ((r[n]) || (m[map[n]]))
        {
          u[n] = distance (&samples[n*D], &mean[map[n]*D]);
          r[n] = false;
          ++stat_dist_calc;
        }

        // 3.b if d(x,c(x)) > l(x,c) or d(x,c(x)) > 1/2d(c(x),c) then
        //      compute d(x,c)
        //      if d(x,c) < d(x,c(x)) then assign c(x) = c
        // compute the distance from n to ALL other clusters
        // and take the minimum, not just the nearest cluster to c(n)
        min_k = map[n];
        min_d = u[n];

        for (k=0; k<K; k++)
        {
          if (!m[k])      // only calculate if the centroid has moved
            continue;    // this is a crucial speed up

                    dist = distance (&samples[n*D], &mean[k*D]);
          ++stat_dist_calc;

          if (dist < min_d)
          {
            min_k = k;
            min_d = dist;
          }
        }

        if (map[n] != min_k)
        {
          map[n] = min_k;
          u[n] = min_d;
          r[n] = false;
        }
      }
      // adjust the means

      // 4. for each center c, let m(c) be the mean of the points assigned to c.
      // calculate the new means and store them in the array nm
      std::cout << "m";
      // adjust the means
      memset (nm, 0, sizeof (float) * K * D);
      adjustMeans (samples, N, nm);
      uk = adjustDistance (mean, nm);

      // 7. replace each center c by m(c)
      for (k=0; k<K; k++)
      {
        m[k] = false;

        for (d=0; d<D; d++)
        {
          if (mean[(k*D)+d] != nm[(k*D)+d])
          {
            m[k] = true;
            break;
          }
        }

        if (m[k])
          ++stat_ks_adjusted;
      }

      // copy the new means into place
      memcpy (mean, nm, sizeof (float) * K * D);

      // output stats
      std::cout << "} " /*<< stat_iteration_changes << " changes"*/;
      std::cout << "\t" << uk << " > " << delta << std::endl;
    }

    // TERMINATION
    delete[] u;
    delete[] sc;
    delete[] s;
    delete[] nm;
    delete[] r;
    delete[] m;

    std::cout << "\tclustered in " << stat_iterations << " iterations" << std::endl;
  }
};
#endif
#endif
