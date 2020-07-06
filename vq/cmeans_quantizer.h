/*
 * C-Means is a centroid detection algorithm which
 * gives a membership vector for each datum to each
 * centroid.
 */
#include <math.h>
#include "vector_quantizer.h"

#define MULTITHREADED_CLUSTERING
//#undef MULTITHREADED_CLUSTERING

template<class T, unsigned char D>
class CMeansQuantizer : public VectorQuantizer<T, D>
{
protected:
  const float    m;    // fuzzy parameter

  float      *u;    // u[i,j] probability of sample i belonging to cluster j
              // order is i major (u[(i*K)+k])

#ifdef MULTITHREADED_CLUSTERING
  static const unsigned int  threadCount = 16;

  class ThreadDescription
  {
  public:
    // pointer to class data (shared)
    const CMeansQuantizer<T, D>& quant;
    const float      *samples;

    // volatile data (contained in the quantizer, but pointed to here
    // so we have a non-const pointer)
    float      *u;

    // start and end of the subregions of *samples this thread operates on
    unsigned int  start;
    unsigned int  end;

    // return values
    float      dtest;
    unsigned int  stat_dist_calcs;
    unsigned int  stat_ks_adjusted;

    ThreadDescription (const CMeansQuantizer<T, D>& quantizer, const float *samples)
      : dtest (0.0f), stat_dist_calcs(0), stat_ks_adjusted(0),
      quant(quantizer), samples(samples),
      start(0), end(0)
    {}

    void resetStats ()
    {
      dtest = 0.0f;
      stat_dist_calcs = 0;
      stat_ks_adjusted = 0;
    }
  };

  ThreadDescription  *threadDescription[threadCount];

  static DWORD WINAPI updateThread  (LPVOID lpParam)
  //  (const T* samples, float *u, float *mean, const unsigned int start, const unsigned int end)
  {
    float  kdist, dsum;
    ThreadDescription *td = static_cast<ThreadDescription*>(lpParam);

    td->dtest = 0.0f;
    td->stat_dist_calcs = 0;
    td->stat_ks_adjusted = 0;

    unsigned int  ki;

    for (unsigned int i=td->start; i<td->end; i++)
    {
      ki = i*td->quant.K;

      // sum the distances to all the centroids
      for (unsigned int k=0; k<td->quant.K; k++)
      {
        kdist = td->quant.distance (&td->quant.mean[k*D], &td->samples[i*D]);
        ++td->stat_dist_calcs;
        dsum = 0.0f;

        assert (kdist != std::numeric_limits<float>::quiet_NaN ());

        for (unsigned int l=0; l<td->quant.K; l++)
        {
          dsum += kdist / td->quant.distance (&td->quant.mean[l*D], &td->samples[i*D]);
          ++td->stat_dist_calcs;
        }

        assert (dsum != std::numeric_limits<float>::quiet_NaN ());

        dsum = pow (dsum, 2.0f/(td->quant.m-1.0f));

        assert (dsum != std::numeric_limits<float>::quiet_NaN ());

        dsum = 1.0f/dsum;

        if (dsum - td->u[ki+k] != 0)
          ++td->stat_ks_adjusted;

        td->dtest += abs (dsum - td->u[ki+k]);    // keep track of the difference

        td->u[ki+k] = dsum;        // assign new prob

        assert (td->u[ki+k] > -0.1f);
        assert (td->u[ki+k] < 1.1f);
      }
    }

    return 0;
  }
#endif

  float updateMembershipMatrix (const T* samples, const unsigned int N)
  {
    float  dsum = 0.0f;
    float  kdist = 0.0f;
    float  dtest = 0.0f;

#ifdef MULTITHREADED_CLUSTERING
    HANDLE        tHandles[threadCount];
    const unsigned int  partitionSize = N/threadCount;

    for (unsigned int t=0; t<threadCount; t++)
    {
      threadDescription[t]->resetStats ();
      threadDescription[t]->start = t*partitionSize;
      threadDescription[t]->end = (t+1)*partitionSize;
      threadDescription[t]->u = u;
    }

    for (unsigned int t=0; t<threadCount; t++)
      tHandles[t] = CreateThread (NULL, 0, &CMeansQuantizer<float, D>::updateThread, threadDescription[t], 0, NULL);

    // execute the threads
    WaitForMultipleObjects (threadCount, tHandles, TRUE, INFINITE);

    // close the threads
    for (unsigned int t=0; t<threadCount; t++)
      CloseHandle (tHandles[t]);

    // get the stats
    for (unsigned int t=0; t<threadCount; t++)
    {
      dtest += threadDescription[t]->dtest;
    //  stat_ks_adjusted += threadDescription[t]->stat_ks_adjusted;
    //  stat_dist_calc += threadDescription[t]->stat_dist_calcs;
    }
#else
    for (unsigned int i=0; i<N; i++)
    {
      // sum the distances to all the centroids
      for (unsigned int k=0; k<K; k++)
      {
        kdist = distance (&mean[k*D], &samples[i*D]);
      //  ++stat_dist_calc;
        dsum = 0.0f;

        assert (kdist != std::numeric_limits<float>::quiet_NaN ());

        for (unsigned int l=0; l<K; l++)
        {
          dsum += kdist / distance (&mean[l*D], &samples[i*D]);
      //    ++stat_dist_calc;
        }

        assert (dsum != std::numeric_limits<float>::quiet_NaN ());

        dsum = pow (dsum, 2.0f/(m-1.0f));

        assert (dsum != std::numeric_limits<float>::quiet_NaN ());

        dsum = 1.0f/dsum;

        dtest += dsum - u[(i*K)+k];    // keep track of the difference

        u[(i*K)+k] = dsum;        // assign new prob

        assert (u[(i*K)+k] > -0.1f);
        assert (u[(i*K)+k] < 1.1f);
      }
    }
#endif

    return dtest;
  }
  //void updatePrototypeMatrix (const T* samples, const unsigned int N)
  void adjustMeans (const T* samples, const unsigned int N, T* newMeans)
  {
    unsigned int  i, k, d;

    float  pmean[D], tmean;

    for (k=0; k<K; k++)
    {
      tmean = 0.0f;

      for (d=0; d<D; d++)
        pmean[d] = 0.0f;

      for (i=0; i<N; i++)
      {
        for (d=0; d<D; d++)
          pmean[d] += pow (u[(i*K)+k], m) * samples[(i*D)+d];

        tmean += pow (u[(i*K)+k], m);
      }

      for (d=0; d<D; d++)
        newMeans[(k*D)+d] = pmean[d]/tmean;
    }
  }

public:
  CMeansQuantizer (unsigned int clusterCount, DistanceMeasure<T, D>& distanceMeasure, const unsigned int maxIters=0)
    : VectorQuantizer (clusterCount, distanceMeasure, maxIters), m(2.0f), u(NULL)
  {}

  ~CMeansQuantizer ()
  {
    if (u != NULL)
      delete[] u;

#ifdef MULTITHREADED_CLUSTERING
    for (unsigned int t=0; t<threadCount; t++)
    {
      if (threadDescription[t] != NULL)
        delete threadDescription[t];
    }
#endif
  }

  // computes the means of the clusters and probability matrix
  void computeMeans (const T* samples, const unsigned int N, const unsigned int maxIters = 0)
  {
#ifdef MULTITHREADED_CLUSTERING
    for (unsigned int t=0; t<threadCount; t++)
      threadDescription[t] = new ThreadDescription (*this, samples);
#endif
    float    *tmpMeans = new float[K*D];
    float    uk = 10.0f;

    if (u != NULL)
      delete[] u;

    u = new float[K*N];

    for (unsigned int i=0; i<(K*N); i++)
      u[i] = 0.0f;//((float)rand ())/((float)RAND_MAX);

    // initialise the mean matrix
    initialiseMeans ();
    selectRandomMeans (samples, N);

    resetIterations ();

    while (uk > delta)
    {
      if (maxItersReached ())
      {
        std::cout << "Max iters reached." << std::endl;
        break;
      }

      zeroStats ();

      std::cout << "ITER: " << stat_iterations++ << " {u";

    //  uk = abs (updateMembershipMatrix (samples, N));
      updateMembershipMatrix (samples, N);

      std::cout << "m";
    //  adjustMeans (samples, N, mean);

      memset (tmpMeans, 0, sizeof (float) * K * D);
      adjustMeans (samples, N, tmpMeans);
      uk = adjustDistance (mean, tmpMeans);
      memcpy (mean, tmpMeans, sizeof (float) * K * D);

      std::cout << "} ";// << stat_ks_adjusted << " us adj";
      std::cout << "\t" << uk << " > " << delta << std::endl;
    }

#ifdef MULTITHREADED_CLUSTERING
    // delete the threads
    for (unsigned int t=0; t<threadCount; t++)
    {
      delete threadDescription[t];
      threadDescription[t] = NULL;
    }
#endif

    delete[] tmpMeans;
  }

  unsigned int getCluster (const unsigned int n) const
  {
    float  max = u[(n*K)];
    unsigned int  argmax = 0;

    for (unsigned int k=1; k<K; k++)
    {
      if (max < u[(n*K)+k])
      {
        max = u[(n*K)+k];
        argmax = k;
      }
    }

    return argmax;
  }

  float getClusterProbability (const unsigned int n) const
  {
    float  max = u[(n*K)];

    for (unsigned int k=1; k<K; k++)
    {
      if (max < u[(n*K)+k])
        max = u[(n*K)+k];
    }

    return max;
  }
  const float *getMVector (const unsigned int n) const
  {
    return &u[(n*K)];
  }
};
