template<class F>
class Backward : public HMMAlgorithm<F>
{
protected:
  double    *beta;

  // initialise the last row
  void initialisation ()
  {
    for (unsigned int i=0; i<N; i++)
      beta[((T-1)*N)+i] = 1;
  }

  void induction ()
  {
    double sum;

    for (unsigned int t=T-1; t>0; t--)
    {
      for (unsigned int i=0; i<N; i++)
      {
        sum = 0.0;

        for (unsigned int j=0; j<N; j++)
          sum += hmm.transitionProbability (i,j)*hmm.observationProbability (obs.observation (t), j)*beta[(t*N)+j];

        beta[((t-1)*N)+i] = sum;
      }
    }
  }

public:
  Backward (HMM& hmm, ObservationSequence<F>& obs)
    : HMMAlgorithm (hmm, obs)
  {
    beta = new double[N*T];
  }

  ~Backward ()
  {
    delete[] beta;
  }

  double* getBeta ()
  {
    return beta;
  }
};

template<class F>
class BackwardScaled : public Backward<F>
{
protected:
  double    *scale;

  // initialise the last row
  void initialisation ()
  {
    for (unsigned int i=0; i<N; i++)
    {
      beta[((T-1)*N)+i] = 1/scale[T-1];
    //  beta[((T-1)*N)+i] = 1*scale[T-1];
    //  beta[((T-1)*N)+i] = 1;
    }
  }

  // populate the data
  void induction ()
  {
    double  sum;

    for (unsigned int t=T-1; t>0; t--)
    {
      for (unsigned int i=0; i<N; i++)
      {
        sum = 0.0;

        for (unsigned int j=0; j<N; j++)
          sum += hmm.transitionProbability (i,j)*hmm.observationProbability (obs.observation (t), j)*beta[(t*N)+j];

      //  beta[((t-1)*N)+i] = sum;
        beta[((t-1)*N)+i] = sum/scale[t-1];
      }
    }
  }

public:
  BackwardScaled (HMM& hmm, ObservationSequence<F>& obs, double* forwardScaleVariables)
    : Backward (hmm, obs), scale(forwardScaleVariables)
  {}
};

#if 0
class BackwardLocalScaled : public Backward
{
  // Backward variable with a locally defined scale variable
  // Essentially the same as the scaled forward variables
protected:
  double    *scale;

  // initialise the last row
  void initialisation ()
  {
    scale[T-1] = 0.0;

    for (unsigned int i=0; i<N; i++)
      beta[((T-1)*N)+i] = 1;

    // collect the first row values for scaling
    for (unsigned int i=0; i<N; i++)
      scale[T-1] += beta[((T-1)*N)+i];

    // reset the alpha variables according to the scale value
    for (unsigned int i=0; i<N; i++)
      beta[((T-1)*N)+i] /= scale[T-1];
  }

  void induction ()
  {
    double sum;

    for (unsigned int t=T-1; t>0; t--)
    {
      for (unsigned int i=0; i<N; i++)
      {
        sum = 0.0;

        for (unsigned int j=0; j<N; j++)
          sum += hmm.transitionProbability (i,j)*hmm.observationProbability (obs.sequence[t],j)*beta[(t*N)+j];

        beta[((t-1)*N)+i] = sum;
      }

      // scale the variables
      sum = 0.0;
      scale[t-1] = 0.0;

      // sum the alpha variable for this observation
      for (unsigned int i=0; i<N; i++)
        scale[t-1] += beta[((t-1)*N)+i];

      // scale the alpha variables
      for (unsigned int i=0; i<N; i++)
        beta[((t-1)*N)+i] /= scale[t-1]; //*= scale[t+1];
    }
  }
public:
  BackwardLocalScaled (HMM& hmm, ObservationSequence& obs)
    : Backward (hmm, obs)
  {
    scale = new double[T];
  }

  ~BackwardLocalScaled ()
  {
    delete[] scale;
  }
};
#endif
