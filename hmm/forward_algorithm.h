template<class F>
class Forward : public HMMAlgorithm<F>
{
protected:
  double    *alpha;
  double    probability;

  // initialise the first row
  void initialisation ()
  {
    for (unsigned int i=0; i<N; i++)
      alpha[i] = hmm.initialProbability (i) * hmm.observationProbability (obs.observation (0), i);
  }

  // compute the probabilities of paths through the graph
  void induction ()
  {
    unsigned int  t, j, i;
    double      sum;

    for (t=0; t<(T-1); t++)    // symbols (rows)
    {
      // loop through the states (cols)
      for (j=0; j<N; j++)
      {
        sum = 0.0;

        // sum the probabilites of the previous observation leading to this state
        for (i=0; i<N; i++)
          sum += alpha[(t*N)+i]*hmm.transitionProbability (i, j);

        // the probability of being in this state and observing the symbol at time t+1
        alpha[((t+1)*N)+j] = sum * hmm.observationProbability (obs.observation (t+1), j);
      }
    }
  }

  // sum all the final matrix row values to get the probability of observing the sequence
  void termination ()
  {
    probability = 0.0;

    for (unsigned int i=0; i<N; i++)
      probability += alpha[((T-1)*N)+i];
  }

public:
  Forward (HMM& hmm, ObservationSequence<F>& obs)
    : HMMAlgorithm (hmm, obs)
  {
    alpha = new double[N*T];
  }

  ~Forward ()
  {
    delete[] alpha;
  }

  double *getAlpha ()
  {
    return alpha;
  }
  double getProbability ()
  {
    return probability;
  }
};

template<class F>
class ForwardScaled : public Forward<F>
{
protected:
  double    *scale;

  // initialise the first row
  void initialisation ()
  {
    scale[0] = 0.0;

    // probability of observing o at t = 0 in state i
    for (unsigned int i=0; i<N; i++)
      alpha[i] = hmm.initialProbability (i) * hmm.observationProbability (obs.observation (0), i);

    // collect the first row values for scaling
    for (unsigned int i=0; i<N; i++)
      scale[0] += alpha[i];

    // reset the alpha variables according to the scale value
    for (unsigned int i=0; i<N; i++)
      alpha[i] /= scale[0];
  }

  // compute the probabilities of paths through the graph
  void induction ()
  {
    double  sum = 0.0;

    for (unsigned int t=0; t<(T-1); t++)    // symbols (rows)
    {
      // collect the alpha values from the previous row
      for (unsigned int j=0; j<N; j++)
      {
        sum = 0.0;

        // sum the probabilites of the previous observation leading to this state
        for (unsigned int i=0; i<N; i++)
          sum += alpha[(t*N)+i]*hmm.transitionProbability (i, j);

        // the probability of being in this state and observing the symbol at time t+1
        alpha[((t+1)*N)+j] = sum * hmm.observationProbability (obs.observation (t+1), j);
      }

      sum = 0.0;
      scale[t+1] = 0.0;

      // sum the alpha variable for this observation
      for (unsigned int i=0; i<N; i++)
        scale[t+1] += alpha[((t+1)*N)+i];

      // scale the alpha variables
      for (unsigned int i=0; i<N; i++)
        alpha[((t+1)*N)+i] /= scale[t+1]; //*= scale[t+1];
    }
  }

  void termination ()
  {
  //  return scale[T-1];
    probability = scale[T-1];
  }

public:
  ForwardScaled (HMM& hmm, ObservationSequence<F>& obs)
    : Forward (hmm, obs)
  {
    scale = new double[T];
  }

  ~ForwardScaled ()
  {
    delete[] scale;
  }

  double* getScale ()
  {
    return scale;
  }

  unsigned int getScaleLength ()
  {
    return T;
  }
};