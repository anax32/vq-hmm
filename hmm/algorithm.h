template<class F>
class HMMAlgorithm
{
protected:
  HMM& hmm;
  ObservationSequence<F>& obs;
  unsigned int N, T, V;

  double probability;

  virtual void initialisation () {}
  virtual void induction () {}
  virtual void termination () {}

public:
  HMMAlgorithm (HMM& hmm, ObservationSequence<F>& obs)
    : hmm(hmm),
      obs(obs),
      N(hmm.getStateCount ()),
      T(obs.length),
      V(hmm.getSymbolCount ())
  {}

  void execute ()
  {
    initialisation ();
    induction ();
    termination ();
  }
};
