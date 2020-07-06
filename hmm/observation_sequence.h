template<class T>
class ObservationSequence
{
protected:
  unsigned int    stateChanges;  // counts the number of times symbol n-1 is different from n

  T          *sequence;          // indices into the symbol set
  T          *symbol;            // list of unique symbols

public:
  const unsigned int  length;       // length of the sequence
  const unsigned int  symbolCount;  // number of unique symbols

  ObservationSequence (const T *seq, const unsigned int len, const T *symbolSet, const unsigned int symbolSetSize)
    : length(len), symbolCount(symbolSetSize), stateChanges(0)
  {
    T  lastState = seq[0];

    sequence = new T[length];

    for (unsigned int i=0; i<length; i++)
    {
      sequence[i] = seq[i];

      // count how many times the symbol changes
      if (seq[i] != lastState)
      {
        ++stateChanges;
        lastState = seq[i];
      }
    }

    symbol = new T[symbolCount];

    for (unsigned int i=0; i<symbolCount; i++)
      symbol[i] = symbolSet[i];
  }

  ~ObservationSequence ()
  {
    if (sequence != NULL)
      delete[] sequence;

    if (symbol != NULL)
      delete[] symbol;
  }

  T getSymbol (unsigned int v) const
  {
    return symbol[v];
  }

  T* getSymbolSet () const
  {
    return symbol;
  }

  unsigned int getSymbolCount () const
  {
    return symbolCount;
  }

  unsigned int getNumberOfChanges () const
  {
    return stateChanges;
  }

  T observation (unsigned int i) const
  {
    return sequence[i];
  }
};
