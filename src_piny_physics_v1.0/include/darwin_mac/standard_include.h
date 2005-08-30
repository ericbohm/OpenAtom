#define DARWIN_MAC
#define NO_PRAGMA
#define SIMP_NINT
#define CHARM_ON
#define PUP_ON
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <iostream>

using namespace std;
#define MACOSX_BLAS
#ifdef PUP_ON
#include <pup.h>
#endif
#ifdef CHARM_ON
#include "charm++.h"
#endif
#include "../class_defs/piny_constants.h"
#include "../../friend_lib/proto_friend_lib_entry.h"
