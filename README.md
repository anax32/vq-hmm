# Vector Quanitization and Hidden Markov Model Implementation

This is old code from ~2008 which needs some love, but might
be a useful guide for someone implementing something like this
again.

Few people will need to reimplement these methods in practice.

## Process

For a given dataset:
1) perform vector quantisation to reduce the data to symbols
2) transform the data to a symbol-sequence
3) train a hidden markov model on the symbol-sequence

## Notes

The HMM implementation was completed by following the paper:
+ [*A Tutorial on Hidden Markov Models and Selected Applications in Speech Recognition* by Lawrence Rabiner](https://ieeexplore.ieee.org/document/18626)
also [here](https://www.cs.ubc.ca/~murphyk/Bayes/rabiner.pdf)

The Triangle-Inequality optimization for K-Means was based on this paper:
+ [*Using the Triangle Inequality to Accelerate K-Means* by Charles Elkan](https://www.aaai.org/Papers/ICML/2003/ICML03-022.pdf)
