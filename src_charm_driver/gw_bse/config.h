#ifndef CONFIG_H
#define CONFIG_H

#include "pup.h"

struct GWConfig {
  unsigned K, Q, L, M;
  unsigned occupied_size, unoccupied_size;
  unsigned pipeline_stages;
};
PUPbytes(GWConfig);

#endif
