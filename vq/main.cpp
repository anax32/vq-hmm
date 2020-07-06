#define _USE_MATH_DEFINES

#include <time.h>

#include <cmath>
#include <iostream>
#include <fstream>
#include <iterator>
#include <iomanip>
#include <map>

using std::cout;
using std::endl;
using std::setprecision;
using std::setw;
using std::setfill;
using std::ostream_iterator;
using std::ofstream;

#include "kmeans_quantizer.h"
#include "dct.h"
#include "codebook.h"
//#include "../Biodata/BioData.h"

// argv[1] == codebook size
int main (int argc, char *argv[])
{
  const unsigned char D = 3;
  const unsigned int  k = 30;

#ifdef BIODATA_H
  // read the biodata
  BioData::SampleSet *bds = BioData::SampleSet::memoryMapFile ("../Data/041207_Roy.bds\0");

  float      *accl_data = new float[bds->sampleCount*3];
  unsigned int  i;

  cout << "Reading " << bds->sampleCount << " samples" << endl;

  for (i=0; i < bds->sampleCount; i++)
  {
    accl_data[(i*3)+0] = bds->data[i].accl[0];
    accl_data[(i*3)+1] = bds->data[i].accl[1];
    accl_data[(i*3)+2] = bds->data[i].accl[2];

  //  cout << "[" << accl_data[(i*3)+0] << ", " << accl_data[(i*3)+1] << ", " << accl_data[(i*3)+2] << "]\t";
  //  cout << "[" << bds->data[i].accl[0] << ", " << bds->data[i].accl[1] << ", " << bds->data[i].accl[2] << "]";
  //  cout << endl;
  }

  // quantise the data
  cout << "Vector Quantisation\0" << endl;
  unsigned int s, e;
  s = clock ();
  EuclideanDistance<float, D>  euclid;
  KMeansQuantizer<float, D>  kmeans (k,  accl_data, bds->sampleCount);
  e = clock () - s;
  cout << "\tquantised in " << e/1000 << "s" << endl << endl;

  // output the codebook
  for (unsigned int i=0; i<kmeans.clusterCount; i++)
  {
    cout << setprecision (2);

    cout << i << "\t";// << kmeans.getCluster (i)->getCode () << "\t";

    float  *m = kmeans.getCluster(i)->getMean ();
    std::copy (m, m+D, ostream_iterator<float> (cout, ","));

    cout << endl;
  }

  return 0;

  // output the samples
  for (unsigned int i=0; i<bds->sampleCount; i++)
  {
    cout << setprecision (2);
    std::copy (&accl_data[(i*D)], &accl_data[(i*D)+D], ostream_iterator<float> (cout, ","));
    cout << "\t,";

    // print the nearest mean
    float  *m = kmeans.getNearestMean (&accl_data[i*D]);

  //  cout << "[";
    std::copy (m, m+D, ostream_iterator<float> (cout, ","));
  //  cout << "]\t";

    cout << kmeans.getNearestClusterIndex (&accl_data[i*D]) << endl;
  }

  return 0;
#else
  char      b[256];
  DistanceMeasure<float, D>*  dist = NULL;// = new EuclideanDistance<float, D> ();

  // create the sample data
  SampleData<float, D>    original_data (200, SampleData<float, D>::SIN, 8);
  float *dct = DiscreteCosineTransform::t (original_data.samples, original_data.sampleCount, D);

  for (int i=0; i<4; i++)
  {
    if (dist != NULL)
    {
      delete dist;
      dist = NULL;
    }

    switch (i)
    {
    case 0:
      dist = new EuclideanDistance<float, D> ();
      strncpy_s (b, 256, "Euclidean.txt\0", 256);
      break;
    case 1:
      dist = new MinkowskiDistance<float, D> ();
      strncpy_s (b, 256, "Minkowski.txt\0", 256);
      break;
    case 2:
      dist = new ManhattanDistance<float, D> ();
      strncpy_s (b, 256, "Manhattan.txt\0", 256);
      break;
    case 3:
      dist = new ChebyshevDistance<float, D> ();
      strncpy_s (b, 256, "Chebyshev.txt\0", 256);
      break;
  //  case 4:
  //    dist = new InnerProduct<float, D> ();
  //    strncpy_s (b, 256, "InnerProduct.txt\0", 256);
  //    break;
    }

    ofstream    of (&b[0]);

    // quantise the data
    KMeansQuantizer<float, D>  kmeans (k,  original_data.samples, original_data.sampleCount, *dist);

    // output the samples
    for (unsigned int i=0; i<original_data.sampleCount; i++)
    {
      of << setprecision (2);
      std::copy (&original_data.samples[(i*D)], &original_data.samples[(i*D)+D], ostream_iterator<float> (of, "\t"));
      of << "\t";

      of << setprecision (2);
      std::copy (&dct[(i*D)], &dct[(i*D)+D], ostream_iterator<float> (of, "\t"));
      of << "\t";

      // print the nearest mean
      float  *m = kmeans.getNearestMean (&original_data.samples[i*D]);
      std::copy (m, m+D, ostream_iterator<float> (of, "\t"));

      of << endl;
    }

    of.close ();
  }

  delete dct;

  return 0;
#endif
}
