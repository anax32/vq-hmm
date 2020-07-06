#ifndef CLUSTER_H
#define CLUSTER_H
#include <iostream>

/*
 * Cluster describes a cluster discovered in data
 * by one of the clustering algorithms.
 */
template<class T>
class Cluster
{
public:
  const static unsigned int  LABEL_LENGTH = 64;
  unsigned char        dimensionality;

  T        *centroid; // centroid of the cluster
  T        *min;      // minimum values
  T        *max;      // maximum values
  unsigned int  N;    // number of sample contributions to
                      // the cluster (for stat calculations)

  char      label[LABEL_LENGTH];    // label for this feature

  Cluster ()
    : centroid(NULL), min(NULL), max(NULL), dimensionality (0)
  {
    memset (&label[0], 0, sizeof (char) * LABEL_LENGTH);
    strcpy_s (&label[0], LABEL_LENGTH, "Unnamed feature\0");
  }

  ~Cluster ()
  {
    if (centroid != NULL)
      delete[] centroid;

    if (min != NULL)
      delete[] min;

    if (max != NULL)
      delete[] max;
  }

  void init (unsigned char dims)
  {
    dimensionality = dims;

    centroid = new T[dimensionality];
    min = new T[dimensionality];
    max = new T[dimensionality];

    for (unsigned int d=0; d<dimensionality; d++)
    {
      centroid[d] = 0.0f; //std::numeric_limits<T>::infinity ();
      min[d] = std::numeric_limits<T>::max ();
      max[d] = std::numeric_limits<T>::min ();
    }

    N=0;
  }

  void write (std::ostream& os)
  {
    os.write ((const char *)centroid, sizeof (T) * dimensionality);
    os.write ((const char *)min, sizeof (T) * dimensionality);
    os.write ((const char *)max, sizeof (T) * dimensionality);
    os.write ((const char *)&N, sizeof (unsigned int));
    os.write ((const char *)&label[0], sizeof (char) * LABEL_LENGTH);
  }

  void read (std::istream& is)
  {
    is.read ((char *)centroid, sizeof (T) * dimensionality);
    is.read ((char *)min, sizeof (T) * dimensionality);
    is.read ((char *)max, sizeof (T) * dimensionality);
    is.read ((char *)&N, sizeof (unsigned int));
    is.read ((char *)&label[0], sizeof (char) * LABEL_LENGTH);
  }
};

template<class T>
class ClusterSet
{
public:
  Cluster<T>    *clusters;
  unsigned int  K;
  const unsigned char  dimensionality;

  // hack for communication between the scatterplot and feature explorer
  unsigned int    currentFeature;

  ClusterSet (unsigned char dims)
    : clusters(NULL), K(0), dimensionality(dims), currentFeature(0)
  {}

  void read (std::istream& is)
  {
    unsigned int  tD;

    is.read ((char *)&K, sizeof (unsigned int));
    is.read ((char *)&tD, sizeof (unsigned int));

    if (tD != dimensionality)
      return;

    clusters = new Cluster<T>[K];

    for (unsigned int k=0; k<K; k++)
    {
      clusters[k].init (tD);
      clusters[k].read (is);
    }
  }

  unsigned char getDimensionality () const
  {
    return dimensionality;
  }
};
#endif
