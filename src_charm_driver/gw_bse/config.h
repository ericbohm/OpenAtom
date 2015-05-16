#ifndef CONFIG_H
#define CONFIG_H

#include "pup.h"

struct GWConfig {
  unsigned K, Q, L, M;                      // Number of k points, and psis
  unsigned occupied_size, unoccupied_size;  // Sizes of psi arrays
  unsigned pipeline_stages;                 // Number of stages in L pipeline
  unsigned rows_of_p, rows_per_chare;       // Rows in PMatrix
};
PUPbytes(GWConfig);

#endif
