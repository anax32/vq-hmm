//#define _USE_MATH_DEFINES

#include <iostream>
#include <iterator>
#include <iomanip>

using std::cout;
using std::endl;
using std::setprecision;
using std::fixed;

#include "observation_sequence.h"
#include "hmm.h"
#include "algorithm.h"
#include "forward_algorithm.h"
#include "backward_algorithm.h"
#include "viterbi_algorithm.h"
#include "baum_welch_algorithm.h"

#ifdef TEST_UNSW_MODEL
unsigned int  test[]={0,0,1,2};
unsigned int  train[]={0,0,0};
unsigned int  set[]={'j', 'k', 'l'};
#endif

int main (int argc, char *argv[])
{
  cout << "Hidden Markov Model\0" << endl;

#ifdef TEST_UNSW_MODEL
  cout << "University of New South Wales example" << endl;

//  ObservationSequence obs (test, 4, 3, set);
//  ObservationSequence train (train, 3, 3, set);
  ObservationSequence<unsigned int> obs (train, 3, set, 3);

  HMM  hmm (obs.symbolCount, obs.symbolCount, HMM::RANDOM);

  unsigned int  iter = 0;

  do
  {
    cout << "Iteration " << iter << endl;

    ForwardScaled<unsigned int>  forward (hmm, obs);
    forward.execute ();
    cout << "forward: " << fixed << setprecision(8) << forward.getProbability () << endl;

    // backward probability
    BackwardScaled<unsigned int> backward (hmm, obs, forward.getScale ());
    backward.execute ();

    // viterbi path
    unsigned int *p = new unsigned int[obs.length];
    Viterbi<unsigned int>  viterbi (hmm, obs);
    viterbi.execute ();
    cout << "viterbi: " << fixed << setprecision(8) << viterbi.getProbability () << endl;

    viterbi.getStateSequence (p);

    cout << "s: ";
    std::copy(p, p+obs.length, std::ostream_iterator<int> (cout, ", "));
    cout << endl;

    hmm.createMaxSequence (obs.length, p);

    cout << "s: ";
    std::copy(p, p+obs.length, std::ostream_iterator<int> (cout, ", "));
    cout << endl;

    delete[] p;
  //  delete[] str;

    // baum-welch
  //  BaumWelch  baumWelch (hmm, obs);
  //  baumWelch.execute ();


  } while (++iter < 10);
#endif
  return 0;
}
