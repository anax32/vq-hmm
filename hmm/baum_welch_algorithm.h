#ifndef BAUMWELCH_H
#define BAUMWELCH_H

#undef TEST_NORMALISATION

template<class OT>
class BaumWelch : public HMMAlgorithm<OT>
{
protected:
  ForwardScaled<OT>  *forward;
  BackwardScaled<OT>  *backward;
  double        *gamma;
  double        prevProb, prob;

  const double    DELTA;
  const double    MIN_PROBABILITY;    // minimum probability to prevent null transitions

  // reestimate parameters for the HMM
  void reestimate ()
  {
    double  *alpha = forward->getAlpha ();
    double  *beta = backward->getBeta ();

    // reestimate frequency of state i in time t=0
    for (unsigned int i=0; i<N; i++)
    {
      hmm.setInitialProbability (i, std::max<double>(gamma[i], MIN_PROBABILITY));
    }

    // transition probabilities
    for (unsigned int i=0; i<N; i++)
    {
      double  num, denom=0.0;
#ifdef TEST_NORMALISATION
      double  row_sum = 0.0;
#endif
      // sum all the transitions from i
      for (unsigned int t=0; t<T-1; t++)
      {
        for (unsigned int js=0; js<N; js++)
        {
          denom += alpha[(t*N)+i] * hmm.transitionProbability (i, js) *
                hmm.observationProbability (obs.observation (t+1), js) *
                beta[((t+1)*N)+js];
        }
      }

      // sum the transitions between i and j
      for (unsigned int j=0; j<N; j++)
      {
        num = 0.0;

        for (unsigned int t=0; t<T-1; t++)
        {
          num += alpha[(t*N)+i] * hmm.transitionProbability (i, j) *
              hmm.observationProbability (obs.observation (t+1), j) *
              beta[((t+1)*N)+j];
        }

        double v = num/denom;

        // set the transition
        hmm.setTransitionProbability (i, j, v); //std::max<double>(v, MIN_PROBABILITY));
#ifdef TEST_NORMALISATION
        row_sum += v;
#endif
      }

#ifdef TEST_NORMALISATION
      if (row_sum != 1.0)
        row_sum -= 1.0;
#endif
    }

    // observation probabilities
    for (unsigned int i=0; i<N; i++)
    {
      double  denom=0.0;
      double  num, p;
#ifdef TEST_NORMALISATION
      double sym_sum=0.0;
#endif

      // sum all the observations in i
      for (unsigned int t=0; t<T; t++)
        denom += gamma[(t*N)+i];

      // sum the observations of a particular symbol
      for (unsigned int v=0; v<V; v++)
      {
        num = 0.0;

        for (unsigned int t=0; t<T; t++)
        {
          if (obs.observation (t) == v)
          {
            num += gamma[(t*N)+i];
          }
        }

        p = num/denom;
        hmm.setObservationProbability (v, i, std::max<double>(p, MIN_PROBABILITY));
#ifdef TEST_NORMALISATION
        sym_sum += p;
#endif
      }

#ifdef TEST_NORMALISATION
      if (sym_sum != 1.0)
        sym_sum -= 1.0;
#endif
    }

    // renormalise the hmm incase we used min-probabilities
    hmm.normalise ();
  }

  double evalreestimation ()
  {
    double  delta;

    // get the new probabilities
    // run the algorithms
    forward->execute ();
    backward->execute ();

    // get the new probability
    prob = forward->getProbability ();

    // recalculate the gamma table
    calculateGamma (forward->getAlpha (), backward->getBeta (), gamma);

    // compute difference between log probability of two iterations
    delta = abs (log (prob) - log (prevProb));
    prevProb = prob;

    return delta;
  }

  void calculateGamma (double *alpha, double *beta, double *gamma)
  {
    double  sum;

    // note we cannot calculate gamma from xi as gamma ranges
    // over all T xi only ranges over 0..T-1
    for (unsigned int t=0; t<T; t++)
    {
      sum = 0.0;

      for (unsigned int i=0; i<N; i++)
      {
        gamma[(t*N)+i] = alpha[(t*N)+i]*beta[(t*N)+i];
        sum += gamma[(t*N)+i];
      }

      for (unsigned int i=0; i<N; i++)
      {
        gamma[(t*N)+i] /= sum;
      }
#ifdef TEST_NORMALISATION
      sum = 0;

      for (unsigned int i=0; i<N; i++)
        sum += gamma[(t*N)+i];

      if (!((sum > 1-DELTA) && (sum < 1+DELTA)))
        assert (((sum > 1-DELTA) && (sum < 1+DELTA)));
#endif
    }
  }

  void initialisation ()
  {
    init ();
  }

  void induction ()
  {
    do
    {
      reestimate ();
    } while (!convergence ());
  }

  void termination ()      {}

public:
  BaumWelch (HMM& hmm, ObservationSequence<OT>& obs)
    : HMMAlgorithm (hmm, obs), DELTA(0.001), MIN_PROBABILITY(std::numeric_limits<double>::min ())
  {
    // link the algorithms
    forward = new ForwardScaled<OT> (hmm, obs);
    backward = new BackwardScaled<OT> (hmm, obs, forward->getScale ());

    gamma = new double[N*T];
  }

  ~BaumWelch ()
  {
    delete forward;
    delete backward;
    delete[] gamma;
  }


  // manual calling of the algorithm (to track iterative steps)
  void init ()
  {
    // run the algorithms
    forward->execute ();
    backward->execute ();

    // get the variables from the algorithms
  //  alpha = forward->getAlpha ();
  //  beta = backward->getBeta ();

    // get the probability
    prevProb = prob = forward->getProbability ();

    // calculate the gamma table
    calculateGamma (forward->getAlpha (), backward->getBeta (), gamma);
  }

  void step ()
  {
    reestimate ();
  }

  bool convergence ()
  {
    return (evalreestimation () < DELTA);
  }
};
#endif
