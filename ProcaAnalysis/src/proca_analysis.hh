#ifndef __proca_analysis_h
#define __proca_analysis_h

#include <math.h>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"
#include "util_Table.h"

void proca_analysis_point_calc(CCTK_ARGUMENTS);

void proca_analysis_reduce(CCTK_ARGUMENTS);

#endif

