#define IBM_ESSL 
#define NO_PRAGMA
#define SIMP_NINT
#define CHARM_ON
#define PUP_ON
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <iostream>

using namespace std;

#ifdef PUP_ON
#include <pup.h>
#endif
#ifdef CHARM_ON
#include "charm++.h"
#endif
#include "../class_defs/piny_constants.h"
#include "../../friend_lib/proto_friend_lib_entry.h"
