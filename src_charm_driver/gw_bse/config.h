#ifndef CONFIG_H
#define CONFIG_H

#include "pup.h"

struct GWConfig {
  unsigned K, L, M;           // Number of k points, and psis
  unsigned n_elems;           // Number of elements in psi and f
  unsigned pipeline_stages;   // Number of stages in L pipeline
  unsigned rows_per_chare;    // Rows per PMatrix chare
};
PUPbytes(GWConfig);

#endif
