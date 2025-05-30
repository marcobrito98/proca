#ifndef __ns_analysis_h
#define __ns_analysis_h

#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "util_Table.h"

void ns_analysis_point_calc(CCTK_ARGUMENTS);

void ns_analysis_reduce(CCTK_ARGUMENTS);

#endif

