#ifndef _STANDARD_INCLUDE_GWBSE_
#define _STANDARD_INCLUDE_GWBSE_


#define CHARM_ON
#define PUP_ON
#include <cstring>
#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <ctime>
#include <iostream>

#ifdef PUP_ON //you can turn off 
#include <pup.h>
#endif
#ifdef CHARM_ON  //you can turn off
#include "charm++.h"
#endif
#include "gwbse_constants.h"
#include "proto_friend_lib_entry.h"

#endif
