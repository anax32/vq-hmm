template<class F>
class Viterbi : public HMMAlgorithm<F>
{
protected:
  double      *delta;
  unsigned int  *psi;
  double      probability;
  unsigned int  finalState;

  // initial probabilities
  void initialisation ()
  {
    for (unsigned int i=0; i<N; i++)
    {
      delta[i] = hmm.initialProbability (i) * hmm.observationProbability (obs.sequence[0], i);
      psi[i] = 0;
    }
  }

  // fill the matrix
  void induction ()
  {
    unsigned int  maxarg;
    double      max, val;

    for (unsigned int t=1; t<T; t++)
    {
      for (unsigned int j=0; j<N; j++)
      {
        max = -1;
        maxarg = 1;

        for (unsigned int i=0; i<N; i++)
        {
          val = delta[((t-1)*N)+i]*hmm.transitionProbability (i, j);

          if (val > max)
          {
            max = val;
            maxarg = i;
          }
        }

        // add the observation probabilty to the delta value
        delta[(t*N)+j] = max*hmm.observationProbability (obs.sequence[t], j);
        psi[(t*N)+j] = maxarg;
      }
    }
  }

  // final probability and state
  void termination ()
  {
    probability = -1.0;

    // get the most probable final state
    for (unsigned int i=0; i<N; i++)
    {
      if (delta[((T-1)*N)+i] > probability)
      {
        probability = delta[((T-1)*N)+i];
        finalState = i;
      }
    }
  }

public:
  Viterbi (HMM& hmm, ObservationSequence<F>& obs)
    : HMMAlgorithm<F>(hmm, obs)
  {
    delta = new double[N*T];
    psi = new unsigned int[N*T];
  }
  ~Viterbi ()
  {
    delete[] delta;
    delete[] psi;
  }

  // viterbi algorithm
  double getProbability ()
  {
    return probability;
  }

  void getStateSequence (unsigned int *sequence)
  {
    unsigned int t = T-1;

    // set the last state in the sequence for backtracking
    sequence[t] = finalState;

    // build the path (state sequence) backtracking
    do
    {
      --t;  // decrement first so we can do t==0 with an unsigned int
      sequence[t] = psi[((t+1)*N)+sequence[t+1]];
    } while (t>0);
  }
};
