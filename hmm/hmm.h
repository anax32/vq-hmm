// University of new south wales model http://www.cse.unsw.edu.au/~waleed/phd/tr9806/node12.html
#define  TEST_UNSW_MODEL    1

#ifdef RABINER_TEST_MODEL
#define  RAB_HEADS  0
#define RAB_TAILS  1
#endif

#define RANDOM_PROBS  1

#include <cmath>
#include <assert.h>

class HMM
{
protected:
  double      *A;    // probability of transiting from state i to j.
                     // i is the row, j is the col

  double      *B;    // probability of observing symbol k in state j.
                     // j is the row, k is the col

  double      *P;    // probability of intially being in state i

public:
  const unsigned int  N;  // state count
  const unsigned int  V;  // symbol count

  enum Type
  {
    RANDOM,
    BAKIS
  };

  const Type  type;

protected:
  void setInitialTransitionProbabilities ()
  {
    // set the initial probabilities to random value
    for (unsigned int i=0; i<N; i++)
    {
      for (unsigned int j=0; j<N; j++)
      {
        if (type == BAKIS)
        {
          if (j<i)
            setTransitionProbability (i, j, 0);
        }
        else
          setTransitionProbability (i, j, (double)rand () / (double)RAND_MAX);
      }
    }
  }

  void setInitialEmissionProbabilities ()
  {
    // set the initial probabilities to random value
    for (unsigned int j=0; j<N; j++)
      for (unsigned int k=0; k<V; k++)
        setObservationProbability (k, j, (double)rand () / (double)RAND_MAX);
  }
  void setInitialStateProbabilities ()
  {
    if (type == BAKIS)
    {
      P[0] = 1;

      for (unsigned int i=1; i<N; i++)
        P[i] = 0;
    }
    else
    {
      for (unsigned int i=0; i<N; i++)
        P[i] = (double)rand () / (double)RAND_MAX;
    }
  }

public:
  HMM (unsigned int stateCount, unsigned int symbolCount, Type modelType)
    : N(stateCount), V(symbolCount), type(modelType)
  {
    // create the probability matrices
    A = new double[N*N];
    B = new double[N*V];
    P = new double[N];

    setInitialTransitionProbabilities ();
    setInitialEmissionProbabilities ();
    setInitialStateProbabilities ();

    // always normalise
    normalise ();
  }

  ~HMM ()
  {
    delete[] A;
    delete[] B;
    delete[] P;
  }


  // get probabilities
  double transitionProbability (unsigned int i, unsigned int j)
  {
    return A[(i*N)+j];
  }

  double observationProbability (unsigned int k, unsigned int j)
  {
  //  return B[(j*V)+k];
    return B[(k*N)+j];
  }

  double initialProbability (unsigned int i)
  {
    return P[i];
  }

  // set probabilities
  void setTransitionProbability (unsigned int i, unsigned int j, double p)
  {
    A[(i*N)+j] = p;
  }

  void setObservationProbability (unsigned int k, unsigned int j, double p)
  {
    B[(k*N)+j] = p;
  }

  void setInitialProbability (unsigned int i, double p)
  {
    P[i] = p;
  }

  // normalise the hmm (after setting probabilities)
  void normalise ()
  {
    double      sum = 0;
    unsigned int  i, j, k;

    // normalise all the transition probabilties to total 1

    // Pi
    if (type == BAKIS)
    {
      P[0] = 1;

      for (i=1; i<N; i++)
        P[i] = 0;
    }
    else
    {
      for (i=0; i<N; i++)
        sum += P[i];

      for (i=0; i<N; i++)
        P[i] /= sum;
    }

    // A[i,j]
    for (i=0; i<N; i++)
    {
      sum = 0.0;

      if ((type == BAKIS) && (i==N-1))    // last state
      {
        setTransitionProbability (i, i, 1);

        for (j=0; j<N; j++)
          setTransitionProbability (i, j, 0);

        continue;
      }

      // first sum transitions
      for (j=0; j<N; j++)
        sum += transitionProbability (i, j);

      // normalise the values for this state
      // so the sum of probabilities is 1
      for (j=0; j<N; j++)
      {
        if ((type == BAKIS) && (j<i))
          continue;
        else
          setTransitionProbability (i, j, transitionProbability (i,j) / sum);
      }
    }
    // B[k,j]
    for (i=0; i<N; i++)
    {
      sum = 0.0;

      // first set transitions to a random value
      for (k=0; k<V; k++)
        sum += observationProbability (k, i);

      // normalise the values for this state
      for (k=0; k<V; k++)
        setObservationProbability (k, i, observationProbability (k, i) / sum);
    }
  }
  // get stats
  unsigned int getStateCount ()
  {
    return N;
  }

  unsigned int getSymbolCount ()
  {
    return V;
  }

  const double* getA ()
  {
    return A;
  }
  const double* getB ()
  {
    return B;
  }
  // FIXME: complete these by making the actual algorithms private
  // and only calling them when model variables/observations have changed
  // path alternatives
  template<class T>
  void createSequenceFromPath (ObservationSequence<T>& obs, unsigned int *statePath, T *string)
  {
    double  max = -1;
    double  prob;
    unsigned int maxarg;

    for (unsigned int i=0; i<obs.length; i++)
    {
      max = -1;

      for (unsigned int j=0; j<V; j++)
      {
        prob = observationProbability (j, path[i]);

        if (prob > max)
        {
          max = prob;
          maxarg = j;//obs.sequence[j]);
        }
      }

      string[i] = obs.getSymbol (maxarg);
    }
  }

  enum PATH_TYPE  {MIN, MAX};

  // FIXME: This function is wrong
  template<class T>
  void createStateSequence (unsigned int length, T *path, PATH_TYPE type)
  {
    double        best;
    unsigned int  bestarg;
    double        pr;
    unsigned int  state;

    // choose the initial state
    best = initialProbability (0);
    bestarg = 0;

    for (unsigned int i=1; i<N; i++)
    {
      pr = initialProbability (i);

      switch (type)
      {
        case MIN:
        {
          if (pr < best)
          {
            best = pr;
            bestarg = i;
          }
          break;
        }
        case MAX:
        {
          if (pr > best)
          {
            best = pr;
            bestarg = i;
          }
          break;
        }
      }
    }

    state = bestarg;

    // create the sequence
    for (unsigned int t=0; t<length; t++)
    {
      // get the most likely observation
      best = observationProbability (0, state);
      bestarg = 0;

      for (unsigned int k=1; k<V; k++)
      {
        pr = observationProbability (k, state);

        switch (type)
        {
          case MIN:
          {
            if (pr < best)
            {
              best = pr;
              bestarg = k;
            }
            break;
          }
          case MAX:
          {
            if (pr > best)
            {
              best = pr;
              bestarg = k;
            }
            break;
          }
        }
      }

      path[t] = bestarg;  // store the symbol

      // move to the most probable next state
      best = transitionProbability (state, 0);
      bestarg = 0;

      for (unsigned int i=1; i<N; i++)
      {
        pr = transitionProbability (state, i);

        switch (type)
        {
          case MIN:
          {
            if (pr < best)
            {
              best = pr;
              bestarg = i;
            }
            break;
          }
          case MAX:
          {
            if (pr > best)
            {
              best = pr;
              bestarg = i;
            }
            break;
          }
        }
      }

      state = bestarg;
    }
  }

  template<class T>
  void createMaxSequence(unsigned int length, T *path)
  {
    return createStateSequence(length, path, MAX);
  }
};
