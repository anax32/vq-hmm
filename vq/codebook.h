/*
 * codebook lookup for a set of symbols to a set of values
 *
 * type parameters are:
 *   T : data type
 *   D : arity of value
 */
template<class T, unsigned char D>
class Codebook
{
private:
  T          *data;
  const unsigned int  size;

public:
  Codebook (T* values, const unsigned int valueCount)
    : size(valueCount)
  {
    data = new float[size*D];

    for (unsigned int i=0; i<size*D; i++)
      for (unsigned int j=0; j<D; j++)
        data[i][j] = values[i][j];
  }

  ~Codebook ()
  {
    delete[] data;
  }

  float getValueForSymbol (const unsigned char symbol) const
  {
    return data[symbol-firstSymbol];
  }

  unsigned char getSymbolForValue (const float f) const
  {
    for (unsigned int i=0; i<size; i++)
    {
      if (data[i] == f)
        return symbols[i];
    }

    return NULL;
  }

  unsigned int getIndexForSymbol (unsigned char sym) const
  {
        for (unsigned int i=0; i<size; i++)
    {
      if (symbols[i] == sym)
        return i;
    }

    return 0xFFFFFFFF;
  }

  unsigned int getIndexForValue (const float f) const
  {
    return getIndexForSymbol (getSymbolForValue (f));
  }

  float getData (const unsigned int index) const
  {
    return data[index];
  }

  unsigned char getSymbol (const unsigned int index) const
  {
    return symbols[index];
  }

  unsigned int getSize () const
  {
    return size;
  }
  const unsigned char* getSymbolSet () const
  {
    return symbols;
  }
};
