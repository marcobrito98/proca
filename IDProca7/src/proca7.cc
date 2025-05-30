/*  File produced by Kranc */

#define KRANC_C

#include <algorithm>
#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "Kranc.hh"
#include "Differencing.h"
#include "loopcontrol.h"

namespace IDProca7 {

extern "C" void proca7_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % proca7_calc_every != proca7_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ADMBase::curv","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ADMBase::curv.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_curv","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN::ML_curv.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ProcaBase::E1i","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ProcaBase::E1i.");
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ProcaBase::E2i","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ProcaBase::E2i.");
  return;
}

static void proca7_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  /* Include user-supplied include files */
  /* Initialise finite differencing variables */
  const ptrdiff_t di CCTK_ATTRIBUTE_UNUSED = 1;
  const ptrdiff_t dj CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,1,0) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t dk CCTK_ATTRIBUTE_UNUSED = 
    CCTK_GFINDEX3D(cctkGH,0,0,1) - CCTK_GFINDEX3D(cctkGH,0,0,0);
  const ptrdiff_t cdi CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * di;
  const ptrdiff_t cdj CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dj;
  const ptrdiff_t cdk CCTK_ATTRIBUTE_UNUSED = sizeof(CCTK_REAL) * dk;
  const ptrdiff_t cctkLbnd1 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[0];
  const ptrdiff_t cctkLbnd2 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[1];
  const ptrdiff_t cctkLbnd3 CCTK_ATTRIBUTE_UNUSED = cctk_lbnd[2];
  const CCTK_REAL t CCTK_ATTRIBUTE_UNUSED = cctk_time;
  const CCTK_REAL cctkOriginSpace1 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(0);
  const CCTK_REAL cctkOriginSpace2 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(1);
  const CCTK_REAL cctkOriginSpace3 CCTK_ATTRIBUTE_UNUSED = 
    CCTK_ORIGIN_SPACE(2);
  const CCTK_REAL dt CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_TIME;
  const CCTK_REAL dx CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(0);
  const CCTK_REAL dy CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(1);
  const CCTK_REAL dz CCTK_ATTRIBUTE_UNUSED = CCTK_DELTA_SPACE(2);
  const CCTK_REAL dxi CCTK_ATTRIBUTE_UNUSED = pow(dx,-1);
  const CCTK_REAL dyi CCTK_ATTRIBUTE_UNUSED = pow(dy,-1);
  const CCTK_REAL dzi CCTK_ATTRIBUTE_UNUSED = pow(dz,-1);
  const CCTK_REAL khalf CCTK_ATTRIBUTE_UNUSED = 0.5;
  const CCTK_REAL kthird CCTK_ATTRIBUTE_UNUSED = 
    0.333333333333333333333333333333;
  const CCTK_REAL ktwothird CCTK_ATTRIBUTE_UNUSED = 
    0.666666666666666666666666666667;
  const CCTK_REAL kfourthird CCTK_ATTRIBUTE_UNUSED = 
    1.33333333333333333333333333333;
  const CCTK_REAL hdxi CCTK_ATTRIBUTE_UNUSED = 0.5*dxi;
  const CCTK_REAL hdyi CCTK_ATTRIBUTE_UNUSED = 0.5*dyi;
  const CCTK_REAL hdzi CCTK_ATTRIBUTE_UNUSED = 0.5*dzi;
  /* Initialize predefined quantities */
  const CCTK_REAL p1o1024dx CCTK_ATTRIBUTE_UNUSED = 0.0009765625*pow(dx,-1);
  const CCTK_REAL p1o1024dy CCTK_ATTRIBUTE_UNUSED = 0.0009765625*pow(dy,-1);
  const CCTK_REAL p1o1024dz CCTK_ATTRIBUTE_UNUSED = 0.0009765625*pow(dz,-1);
  const CCTK_REAL p1o120dx CCTK_ATTRIBUTE_UNUSED = 0.00833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL p1o120dy CCTK_ATTRIBUTE_UNUSED = 0.00833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL p1o120dz CCTK_ATTRIBUTE_UNUSED = 0.00833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL p1o12dx CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL p1o12dy CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL p1o12dz CCTK_ATTRIBUTE_UNUSED = 0.0833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL p1o144dxdy CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o144dxdz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o144dydz CCTK_ATTRIBUTE_UNUSED = 0.00694444444444444444444444444444*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o1680dx CCTK_ATTRIBUTE_UNUSED = 0.000595238095238095238095238095238*pow(dx,-1);
  const CCTK_REAL p1o1680dy CCTK_ATTRIBUTE_UNUSED = 0.000595238095238095238095238095238*pow(dy,-1);
  const CCTK_REAL p1o1680dz CCTK_ATTRIBUTE_UNUSED = 0.000595238095238095238095238095238*pow(dz,-1);
  const CCTK_REAL p1o180dx2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dx,-2);
  const CCTK_REAL p1o180dy2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dy,-2);
  const CCTK_REAL p1o180dz2 CCTK_ATTRIBUTE_UNUSED = 0.00555555555555555555555555555556*pow(dz,-2);
  const CCTK_REAL p1o24dx CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL p1o24dy CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL p1o24dz CCTK_ATTRIBUTE_UNUSED = 0.0416666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL p1o2dx CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dx,-1);
  const CCTK_REAL p1o2dy CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dy,-1);
  const CCTK_REAL p1o2dz CCTK_ATTRIBUTE_UNUSED = 0.5*pow(dz,-1);
  const CCTK_REAL p1o3600dxdy CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o3600dxdz CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o3600dydz CCTK_ATTRIBUTE_UNUSED = 0.000277777777777777777777777777778*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dx CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1);
  const CCTK_REAL p1o4dxdy CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o4dxdz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dy CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dy,-1);
  const CCTK_REAL p1o4dydz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o4dz CCTK_ATTRIBUTE_UNUSED = 0.25*pow(dz,-1);
  const CCTK_REAL p1o5040dx2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dx,-2);
  const CCTK_REAL p1o5040dy2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dy,-2);
  const CCTK_REAL p1o5040dz2 CCTK_ATTRIBUTE_UNUSED = 0.000198412698412698412698412698413*pow(dz,-2);
  const CCTK_REAL p1o560dx CCTK_ATTRIBUTE_UNUSED = 0.00178571428571428571428571428571*pow(dx,-1);
  const CCTK_REAL p1o560dy CCTK_ATTRIBUTE_UNUSED = 0.00178571428571428571428571428571*pow(dy,-1);
  const CCTK_REAL p1o560dz CCTK_ATTRIBUTE_UNUSED = 0.00178571428571428571428571428571*pow(dz,-1);
  const CCTK_REAL p1o60dx CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL p1o60dy CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL p1o60dz CCTK_ATTRIBUTE_UNUSED = 0.0166666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL p1o64dx CCTK_ATTRIBUTE_UNUSED = 0.015625*pow(dx,-1);
  const CCTK_REAL p1o64dy CCTK_ATTRIBUTE_UNUSED = 0.015625*pow(dy,-1);
  const CCTK_REAL p1o64dz CCTK_ATTRIBUTE_UNUSED = 0.015625*pow(dz,-1);
  const CCTK_REAL p1o705600dxdy CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dx,-1)*pow(dy,-1);
  const CCTK_REAL p1o705600dxdz CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dx,-1)*pow(dz,-1);
  const CCTK_REAL p1o705600dydz CCTK_ATTRIBUTE_UNUSED = 1.41723356009070294784580498866e-6*pow(dy,-1)*pow(dz,-1);
  const CCTK_REAL p1o840dx CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dx,-1);
  const CCTK_REAL p1o840dy CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dy,-1);
  const CCTK_REAL p1o840dz CCTK_ATTRIBUTE_UNUSED = 0.00119047619047619047619047619048*pow(dz,-1);
  const CCTK_REAL p1odx2 CCTK_ATTRIBUTE_UNUSED = pow(dx,-2);
  const CCTK_REAL p1ody2 CCTK_ATTRIBUTE_UNUSED = pow(dy,-2);
  const CCTK_REAL p1odz2 CCTK_ATTRIBUTE_UNUSED = pow(dz,-2);
  const CCTK_REAL pm1o120dx CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL pm1o120dy CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL pm1o120dz CCTK_ATTRIBUTE_UNUSED = -0.00833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL pm1o12dx CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dx,-1);
  const CCTK_REAL pm1o12dx2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dx,-2);
  const CCTK_REAL pm1o12dy CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dy,-1);
  const CCTK_REAL pm1o12dy2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dy,-2);
  const CCTK_REAL pm1o12dz CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dz,-1);
  const CCTK_REAL pm1o12dz2 CCTK_ATTRIBUTE_UNUSED = -0.0833333333333333333333333333333*pow(dz,-2);
  const CCTK_REAL pm1o16dx CCTK_ATTRIBUTE_UNUSED = -0.0625*pow(dx,-1);
  const CCTK_REAL pm1o16dy CCTK_ATTRIBUTE_UNUSED = -0.0625*pow(dy,-1);
  const CCTK_REAL pm1o16dz CCTK_ATTRIBUTE_UNUSED = -0.0625*pow(dz,-1);
  const CCTK_REAL pm1o256dx CCTK_ATTRIBUTE_UNUSED = -0.00390625*pow(dx,-1);
  const CCTK_REAL pm1o256dy CCTK_ATTRIBUTE_UNUSED = -0.00390625*pow(dy,-1);
  const CCTK_REAL pm1o256dz CCTK_ATTRIBUTE_UNUSED = -0.00390625*pow(dz,-1);
  const CCTK_REAL pm1o2dx CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dx,-1);
  const CCTK_REAL pm1o2dy CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dy,-1);
  const CCTK_REAL pm1o2dz CCTK_ATTRIBUTE_UNUSED = -0.5*pow(dz,-1);
  const CCTK_REAL pm1o4dx CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dx,-1);
  const CCTK_REAL pm1o4dy CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dy,-1);
  const CCTK_REAL pm1o4dz CCTK_ATTRIBUTE_UNUSED = -0.25*pow(dz,-1);
  const CCTK_REAL pm1o60dx CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL pm1o60dy CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL pm1o60dz CCTK_ATTRIBUTE_UNUSED = -0.0166666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL pm1o6dx CCTK_ATTRIBUTE_UNUSED = -0.166666666666666666666666666667*pow(dx,-1);
  const CCTK_REAL pm1o6dy CCTK_ATTRIBUTE_UNUSED = -0.166666666666666666666666666667*pow(dy,-1);
  const CCTK_REAL pm1o6dz CCTK_ATTRIBUTE_UNUSED = -0.166666666666666666666666666667*pow(dz,-1);
  const CCTK_REAL pm1o840dx CCTK_ATTRIBUTE_UNUSED = -0.00119047619047619047619047619048*pow(dx,-1);
  const CCTK_REAL pm1o840dy CCTK_ATTRIBUTE_UNUSED = -0.00119047619047619047619047619048*pow(dy,-1);
  const CCTK_REAL pm1o840dz CCTK_ATTRIBUTE_UNUSED = -0.00119047619047619047619047619048*pow(dz,-1);
  /* Assign local copies of arrays functions */
  
  
  /* Calculate temporaries and arrays functions */
  /* Copy local copies back to grid functions */
  /* Loop over the grid points */
  const int imin0=imin[0];
  const int imin1=imin[1];
  const int imin2=imin[2];
  const int imax0=imax[0];
  const int imax1=imax[1];
  const int imax2=imax[2];
  #pragma omp parallel
  CCTK_LOOP3(proca7,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL A1xL CCTK_ATTRIBUTE_UNUSED = A1x[index];
    CCTK_REAL A1yL CCTK_ATTRIBUTE_UNUSED = A1y[index];
    CCTK_REAL A1zL CCTK_ATTRIBUTE_UNUSED = A1z[index];
    CCTK_REAL A2xL CCTK_ATTRIBUTE_UNUSED = A2x[index];
    CCTK_REAL A2yL CCTK_ATTRIBUTE_UNUSED = A2y[index];
    CCTK_REAL A2zL CCTK_ATTRIBUTE_UNUSED = A2z[index];
    CCTK_REAL alpL CCTK_ATTRIBUTE_UNUSED = alp[index];
    CCTK_REAL Aphi1L CCTK_ATTRIBUTE_UNUSED = Aphi1[index];
    CCTK_REAL Aphi2L CCTK_ATTRIBUTE_UNUSED = Aphi2[index];
    CCTK_REAL At11L CCTK_ATTRIBUTE_UNUSED = At11[index];
    CCTK_REAL At12L CCTK_ATTRIBUTE_UNUSED = At12[index];
    CCTK_REAL At13L CCTK_ATTRIBUTE_UNUSED = At13[index];
    CCTK_REAL At22L CCTK_ATTRIBUTE_UNUSED = At22[index];
    CCTK_REAL At23L CCTK_ATTRIBUTE_UNUSED = At23[index];
    CCTK_REAL At33L CCTK_ATTRIBUTE_UNUSED = At33[index];
    CCTK_REAL betaxL CCTK_ATTRIBUTE_UNUSED = betax[index];
    CCTK_REAL betayL CCTK_ATTRIBUTE_UNUSED = betay[index];
    CCTK_REAL betazL CCTK_ATTRIBUTE_UNUSED = betaz[index];
    CCTK_REAL E1xL CCTK_ATTRIBUTE_UNUSED = E1x[index];
    CCTK_REAL E1yL CCTK_ATTRIBUTE_UNUSED = E1y[index];
    CCTK_REAL E1zL CCTK_ATTRIBUTE_UNUSED = E1z[index];
    CCTK_REAL E2xL CCTK_ATTRIBUTE_UNUSED = E2x[index];
    CCTK_REAL E2yL CCTK_ATTRIBUTE_UNUSED = E2y[index];
    CCTK_REAL E2zL CCTK_ATTRIBUTE_UNUSED = E2z[index];
    CCTK_REAL gxxL CCTK_ATTRIBUTE_UNUSED = gxx[index];
    CCTK_REAL gxyL CCTK_ATTRIBUTE_UNUSED = gxy[index];
    CCTK_REAL gxzL CCTK_ATTRIBUTE_UNUSED = gxz[index];
    CCTK_REAL gyyL CCTK_ATTRIBUTE_UNUSED = gyy[index];
    CCTK_REAL gyzL CCTK_ATTRIBUTE_UNUSED = gyz[index];
    CCTK_REAL gzzL CCTK_ATTRIBUTE_UNUSED = gzz[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    CCTK_REAL PDstandardNth1A1x CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2A1x CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3A1x CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1A1y CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2A1y CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3A1y CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1A1z CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2A1z CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3A1z CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1A2x CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2A2x CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3A2x CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1A2y CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2A2y CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3A2y CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1A2z CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2A2z CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3A2z CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3alp CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1Aphi1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2Aphi1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3Aphi1 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1Aphi2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2Aphi2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3Aphi2 CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3betaz CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDstandardNth1A1x = PDstandardNthfdOrder21(&A1x[index]);
        PDstandardNth2A1x = PDstandardNthfdOrder22(&A1x[index]);
        PDstandardNth3A1x = PDstandardNthfdOrder23(&A1x[index]);
        PDstandardNth1A1y = PDstandardNthfdOrder21(&A1y[index]);
        PDstandardNth2A1y = PDstandardNthfdOrder22(&A1y[index]);
        PDstandardNth3A1y = PDstandardNthfdOrder23(&A1y[index]);
        PDstandardNth1A1z = PDstandardNthfdOrder21(&A1z[index]);
        PDstandardNth2A1z = PDstandardNthfdOrder22(&A1z[index]);
        PDstandardNth3A1z = PDstandardNthfdOrder23(&A1z[index]);
        PDstandardNth1A2x = PDstandardNthfdOrder21(&A2x[index]);
        PDstandardNth2A2x = PDstandardNthfdOrder22(&A2x[index]);
        PDstandardNth3A2x = PDstandardNthfdOrder23(&A2x[index]);
        PDstandardNth1A2y = PDstandardNthfdOrder21(&A2y[index]);
        PDstandardNth2A2y = PDstandardNthfdOrder22(&A2y[index]);
        PDstandardNth3A2y = PDstandardNthfdOrder23(&A2y[index]);
        PDstandardNth1A2z = PDstandardNthfdOrder21(&A2z[index]);
        PDstandardNth2A2z = PDstandardNthfdOrder22(&A2z[index]);
        PDstandardNth3A2z = PDstandardNthfdOrder23(&A2z[index]);
        PDstandardNth1alp = PDstandardNthfdOrder21(&alp[index]);
        PDstandardNth2alp = PDstandardNthfdOrder22(&alp[index]);
        PDstandardNth3alp = PDstandardNthfdOrder23(&alp[index]);
        PDstandardNth1Aphi1 = PDstandardNthfdOrder21(&Aphi1[index]);
        PDstandardNth2Aphi1 = PDstandardNthfdOrder22(&Aphi1[index]);
        PDstandardNth3Aphi1 = PDstandardNthfdOrder23(&Aphi1[index]);
        PDstandardNth1Aphi2 = PDstandardNthfdOrder21(&Aphi2[index]);
        PDstandardNth2Aphi2 = PDstandardNthfdOrder22(&Aphi2[index]);
        PDstandardNth3Aphi2 = PDstandardNthfdOrder23(&Aphi2[index]);
        PDstandardNth1betax = PDstandardNthfdOrder21(&betax[index]);
        PDstandardNth2betax = PDstandardNthfdOrder22(&betax[index]);
        PDstandardNth3betax = PDstandardNthfdOrder23(&betax[index]);
        PDstandardNth1betay = PDstandardNthfdOrder21(&betay[index]);
        PDstandardNth2betay = PDstandardNthfdOrder22(&betay[index]);
        PDstandardNth3betay = PDstandardNthfdOrder23(&betay[index]);
        PDstandardNth1betaz = PDstandardNthfdOrder21(&betaz[index]);
        PDstandardNth2betaz = PDstandardNthfdOrder22(&betaz[index]);
        PDstandardNth3betaz = PDstandardNthfdOrder23(&betaz[index]);
        break;
      }
      
      case 4:
      {
        PDstandardNth1A1x = PDstandardNthfdOrder41(&A1x[index]);
        PDstandardNth2A1x = PDstandardNthfdOrder42(&A1x[index]);
        PDstandardNth3A1x = PDstandardNthfdOrder43(&A1x[index]);
        PDstandardNth1A1y = PDstandardNthfdOrder41(&A1y[index]);
        PDstandardNth2A1y = PDstandardNthfdOrder42(&A1y[index]);
        PDstandardNth3A1y = PDstandardNthfdOrder43(&A1y[index]);
        PDstandardNth1A1z = PDstandardNthfdOrder41(&A1z[index]);
        PDstandardNth2A1z = PDstandardNthfdOrder42(&A1z[index]);
        PDstandardNth3A1z = PDstandardNthfdOrder43(&A1z[index]);
        PDstandardNth1A2x = PDstandardNthfdOrder41(&A2x[index]);
        PDstandardNth2A2x = PDstandardNthfdOrder42(&A2x[index]);
        PDstandardNth3A2x = PDstandardNthfdOrder43(&A2x[index]);
        PDstandardNth1A2y = PDstandardNthfdOrder41(&A2y[index]);
        PDstandardNth2A2y = PDstandardNthfdOrder42(&A2y[index]);
        PDstandardNth3A2y = PDstandardNthfdOrder43(&A2y[index]);
        PDstandardNth1A2z = PDstandardNthfdOrder41(&A2z[index]);
        PDstandardNth2A2z = PDstandardNthfdOrder42(&A2z[index]);
        PDstandardNth3A2z = PDstandardNthfdOrder43(&A2z[index]);
        PDstandardNth1alp = PDstandardNthfdOrder41(&alp[index]);
        PDstandardNth2alp = PDstandardNthfdOrder42(&alp[index]);
        PDstandardNth3alp = PDstandardNthfdOrder43(&alp[index]);
        PDstandardNth1Aphi1 = PDstandardNthfdOrder41(&Aphi1[index]);
        PDstandardNth2Aphi1 = PDstandardNthfdOrder42(&Aphi1[index]);
        PDstandardNth3Aphi1 = PDstandardNthfdOrder43(&Aphi1[index]);
        PDstandardNth1Aphi2 = PDstandardNthfdOrder41(&Aphi2[index]);
        PDstandardNth2Aphi2 = PDstandardNthfdOrder42(&Aphi2[index]);
        PDstandardNth3Aphi2 = PDstandardNthfdOrder43(&Aphi2[index]);
        PDstandardNth1betax = PDstandardNthfdOrder41(&betax[index]);
        PDstandardNth2betax = PDstandardNthfdOrder42(&betax[index]);
        PDstandardNth3betax = PDstandardNthfdOrder43(&betax[index]);
        PDstandardNth1betay = PDstandardNthfdOrder41(&betay[index]);
        PDstandardNth2betay = PDstandardNthfdOrder42(&betay[index]);
        PDstandardNth3betay = PDstandardNthfdOrder43(&betay[index]);
        PDstandardNth1betaz = PDstandardNthfdOrder41(&betaz[index]);
        PDstandardNth2betaz = PDstandardNthfdOrder42(&betaz[index]);
        PDstandardNth3betaz = PDstandardNthfdOrder43(&betaz[index]);
        break;
      }
      
      case 6:
      {
        PDstandardNth1A1x = PDstandardNthfdOrder61(&A1x[index]);
        PDstandardNth2A1x = PDstandardNthfdOrder62(&A1x[index]);
        PDstandardNth3A1x = PDstandardNthfdOrder63(&A1x[index]);
        PDstandardNth1A1y = PDstandardNthfdOrder61(&A1y[index]);
        PDstandardNth2A1y = PDstandardNthfdOrder62(&A1y[index]);
        PDstandardNth3A1y = PDstandardNthfdOrder63(&A1y[index]);
        PDstandardNth1A1z = PDstandardNthfdOrder61(&A1z[index]);
        PDstandardNth2A1z = PDstandardNthfdOrder62(&A1z[index]);
        PDstandardNth3A1z = PDstandardNthfdOrder63(&A1z[index]);
        PDstandardNth1A2x = PDstandardNthfdOrder61(&A2x[index]);
        PDstandardNth2A2x = PDstandardNthfdOrder62(&A2x[index]);
        PDstandardNth3A2x = PDstandardNthfdOrder63(&A2x[index]);
        PDstandardNth1A2y = PDstandardNthfdOrder61(&A2y[index]);
        PDstandardNth2A2y = PDstandardNthfdOrder62(&A2y[index]);
        PDstandardNth3A2y = PDstandardNthfdOrder63(&A2y[index]);
        PDstandardNth1A2z = PDstandardNthfdOrder61(&A2z[index]);
        PDstandardNth2A2z = PDstandardNthfdOrder62(&A2z[index]);
        PDstandardNth3A2z = PDstandardNthfdOrder63(&A2z[index]);
        PDstandardNth1alp = PDstandardNthfdOrder61(&alp[index]);
        PDstandardNth2alp = PDstandardNthfdOrder62(&alp[index]);
        PDstandardNth3alp = PDstandardNthfdOrder63(&alp[index]);
        PDstandardNth1Aphi1 = PDstandardNthfdOrder61(&Aphi1[index]);
        PDstandardNth2Aphi1 = PDstandardNthfdOrder62(&Aphi1[index]);
        PDstandardNth3Aphi1 = PDstandardNthfdOrder63(&Aphi1[index]);
        PDstandardNth1Aphi2 = PDstandardNthfdOrder61(&Aphi2[index]);
        PDstandardNth2Aphi2 = PDstandardNthfdOrder62(&Aphi2[index]);
        PDstandardNth3Aphi2 = PDstandardNthfdOrder63(&Aphi2[index]);
        PDstandardNth1betax = PDstandardNthfdOrder61(&betax[index]);
        PDstandardNth2betax = PDstandardNthfdOrder62(&betax[index]);
        PDstandardNth3betax = PDstandardNthfdOrder63(&betax[index]);
        PDstandardNth1betay = PDstandardNthfdOrder61(&betay[index]);
        PDstandardNth2betay = PDstandardNthfdOrder62(&betay[index]);
        PDstandardNth3betay = PDstandardNthfdOrder63(&betay[index]);
        PDstandardNth1betaz = PDstandardNthfdOrder61(&betaz[index]);
        PDstandardNth2betaz = PDstandardNthfdOrder62(&betaz[index]);
        PDstandardNth3betaz = PDstandardNthfdOrder63(&betaz[index]);
        break;
      }
      
      case 8:
      {
        PDstandardNth1A1x = PDstandardNthfdOrder81(&A1x[index]);
        PDstandardNth2A1x = PDstandardNthfdOrder82(&A1x[index]);
        PDstandardNth3A1x = PDstandardNthfdOrder83(&A1x[index]);
        PDstandardNth1A1y = PDstandardNthfdOrder81(&A1y[index]);
        PDstandardNth2A1y = PDstandardNthfdOrder82(&A1y[index]);
        PDstandardNth3A1y = PDstandardNthfdOrder83(&A1y[index]);
        PDstandardNth1A1z = PDstandardNthfdOrder81(&A1z[index]);
        PDstandardNth2A1z = PDstandardNthfdOrder82(&A1z[index]);
        PDstandardNth3A1z = PDstandardNthfdOrder83(&A1z[index]);
        PDstandardNth1A2x = PDstandardNthfdOrder81(&A2x[index]);
        PDstandardNth2A2x = PDstandardNthfdOrder82(&A2x[index]);
        PDstandardNth3A2x = PDstandardNthfdOrder83(&A2x[index]);
        PDstandardNth1A2y = PDstandardNthfdOrder81(&A2y[index]);
        PDstandardNth2A2y = PDstandardNthfdOrder82(&A2y[index]);
        PDstandardNth3A2y = PDstandardNthfdOrder83(&A2y[index]);
        PDstandardNth1A2z = PDstandardNthfdOrder81(&A2z[index]);
        PDstandardNth2A2z = PDstandardNthfdOrder82(&A2z[index]);
        PDstandardNth3A2z = PDstandardNthfdOrder83(&A2z[index]);
        PDstandardNth1alp = PDstandardNthfdOrder81(&alp[index]);
        PDstandardNth2alp = PDstandardNthfdOrder82(&alp[index]);
        PDstandardNth3alp = PDstandardNthfdOrder83(&alp[index]);
        PDstandardNth1Aphi1 = PDstandardNthfdOrder81(&Aphi1[index]);
        PDstandardNth2Aphi1 = PDstandardNthfdOrder82(&Aphi1[index]);
        PDstandardNth3Aphi1 = PDstandardNthfdOrder83(&Aphi1[index]);
        PDstandardNth1Aphi2 = PDstandardNthfdOrder81(&Aphi2[index]);
        PDstandardNth2Aphi2 = PDstandardNthfdOrder82(&Aphi2[index]);
        PDstandardNth3Aphi2 = PDstandardNthfdOrder83(&Aphi2[index]);
        PDstandardNth1betax = PDstandardNthfdOrder81(&betax[index]);
        PDstandardNth2betax = PDstandardNthfdOrder82(&betax[index]);
        PDstandardNth3betax = PDstandardNthfdOrder83(&betax[index]);
        PDstandardNth1betay = PDstandardNthfdOrder81(&betay[index]);
        PDstandardNth2betay = PDstandardNthfdOrder82(&betay[index]);
        PDstandardNth3betay = PDstandardNthfdOrder83(&betay[index]);
        PDstandardNth1betaz = PDstandardNthfdOrder81(&betaz[index]);
        PDstandardNth2betaz = PDstandardNthfdOrder82(&betaz[index]);
        PDstandardNth3betaz = PDstandardNthfdOrder83(&betaz[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL detg CCTK_ATTRIBUTE_UNUSED = 2*gxyL*gxzL*gyzL - 
      gzzL*pow(gxyL,2) + gyyL*(gxxL*gzzL - pow(gxzL,2)) - gxxL*pow(gyzL,2);
    
    CCTK_REAL gu11 CCTK_ATTRIBUTE_UNUSED = (gyyL*gzzL - 
      pow(gyzL,2))*pow(detg,-1);
    
    CCTK_REAL gu12 CCTK_ATTRIBUTE_UNUSED = (gxzL*gyzL - 
      gxyL*gzzL)*pow(detg,-1);
    
    CCTK_REAL gu13 CCTK_ATTRIBUTE_UNUSED = (-(gxzL*gyyL) + 
      gxyL*gyzL)*pow(detg,-1);
    
    CCTK_REAL gu22 CCTK_ATTRIBUTE_UNUSED = (gxxL*gzzL - 
      pow(gxzL,2))*pow(detg,-1);
    
    CCTK_REAL gu23 CCTK_ATTRIBUTE_UNUSED = (gxyL*gxzL - 
      gxxL*gyzL)*pow(detg,-1);
    
    CCTK_REAL gu33 CCTK_ATTRIBUTE_UNUSED = (gxxL*gyyL - 
      pow(gxyL,2))*pow(detg,-1);
    
    CCTK_REAL Efield1u1 CCTK_ATTRIBUTE_UNUSED = (E1xL*omega + 
      gu11*(betaxL*PDstandardNth1A1x - Aphi1L*PDstandardNth1alp - 
      alpL*PDstandardNth1Aphi1 + A1xL*PDstandardNth1betax + 
      A1yL*PDstandardNth1betay + A1zL*PDstandardNth1betaz + 
      betayL*PDstandardNth2A1x + betazL*PDstandardNth3A1x) + 
      gu12*(betaxL*PDstandardNth1A1y + betayL*PDstandardNth2A1y - 
      Aphi1L*PDstandardNth2alp - alpL*PDstandardNth2Aphi1 + 
      A1xL*PDstandardNth2betax + A1yL*PDstandardNth2betay + 
      A1zL*PDstandardNth2betaz + betazL*PDstandardNth3A1y) + 
      gu13*(betaxL*PDstandardNth1A1z + betayL*PDstandardNth2A1z + 
      betazL*PDstandardNth3A1z - Aphi1L*PDstandardNth3alp - 
      alpL*PDstandardNth3Aphi1 + A1xL*PDstandardNth3betax + 
      A1yL*PDstandardNth3betay + A1zL*PDstandardNth3betaz))*pow(alpL,-1);
    
    CCTK_REAL Efield1u2 CCTK_ATTRIBUTE_UNUSED = (E1yL*omega + 
      gu12*(betaxL*PDstandardNth1A1x - Aphi1L*PDstandardNth1alp - 
      alpL*PDstandardNth1Aphi1 + A1xL*PDstandardNth1betax + 
      A1yL*PDstandardNth1betay + A1zL*PDstandardNth1betaz + 
      betayL*PDstandardNth2A1x + betazL*PDstandardNth3A1x) + 
      gu22*(betaxL*PDstandardNth1A1y + betayL*PDstandardNth2A1y - 
      Aphi1L*PDstandardNth2alp - alpL*PDstandardNth2Aphi1 + 
      A1xL*PDstandardNth2betax + A1yL*PDstandardNth2betay + 
      A1zL*PDstandardNth2betaz + betazL*PDstandardNth3A1y) + 
      gu23*(betaxL*PDstandardNth1A1z + betayL*PDstandardNth2A1z + 
      betazL*PDstandardNth3A1z - Aphi1L*PDstandardNth3alp - 
      alpL*PDstandardNth3Aphi1 + A1xL*PDstandardNth3betax + 
      A1yL*PDstandardNth3betay + A1zL*PDstandardNth3betaz))*pow(alpL,-1);
    
    CCTK_REAL Efield1u3 CCTK_ATTRIBUTE_UNUSED = (E1zL*omega + 
      gu13*(betaxL*PDstandardNth1A1x - Aphi1L*PDstandardNth1alp - 
      alpL*PDstandardNth1Aphi1 + A1xL*PDstandardNth1betax + 
      A1yL*PDstandardNth1betay + A1zL*PDstandardNth1betaz + 
      betayL*PDstandardNth2A1x + betazL*PDstandardNth3A1x) + 
      gu23*(betaxL*PDstandardNth1A1y + betayL*PDstandardNth2A1y - 
      Aphi1L*PDstandardNth2alp - alpL*PDstandardNth2Aphi1 + 
      A1xL*PDstandardNth2betax + A1yL*PDstandardNth2betay + 
      A1zL*PDstandardNth2betaz + betazL*PDstandardNth3A1y) + 
      gu33*(betaxL*PDstandardNth1A1z + betayL*PDstandardNth2A1z + 
      betazL*PDstandardNth3A1z - Aphi1L*PDstandardNth3alp - 
      alpL*PDstandardNth3Aphi1 + A1xL*PDstandardNth3betax + 
      A1yL*PDstandardNth3betay + A1zL*PDstandardNth3betaz))*pow(alpL,-1);
    
    CCTK_REAL Efield2u1 CCTK_ATTRIBUTE_UNUSED = (E2xL*omega + 
      gu11*(betaxL*PDstandardNth1A2x - Aphi2L*PDstandardNth1alp - 
      alpL*PDstandardNth1Aphi2 + A2xL*PDstandardNth1betax + 
      A2yL*PDstandardNth1betay + A2zL*PDstandardNth1betaz + 
      betayL*PDstandardNth2A2x + betazL*PDstandardNth3A2x) + 
      gu12*(betaxL*PDstandardNth1A2y + betayL*PDstandardNth2A2y - 
      Aphi2L*PDstandardNth2alp - alpL*PDstandardNth2Aphi2 + 
      A2xL*PDstandardNth2betax + A2yL*PDstandardNth2betay + 
      A2zL*PDstandardNth2betaz + betazL*PDstandardNth3A2y) + 
      gu13*(betaxL*PDstandardNth1A2z + betayL*PDstandardNth2A2z + 
      betazL*PDstandardNth3A2z - Aphi2L*PDstandardNth3alp - 
      alpL*PDstandardNth3Aphi2 + A2xL*PDstandardNth3betax + 
      A2yL*PDstandardNth3betay + A2zL*PDstandardNth3betaz))*pow(alpL,-1);
    
    CCTK_REAL Efield2u2 CCTK_ATTRIBUTE_UNUSED = (E2yL*omega + 
      gu12*(betaxL*PDstandardNth1A2x - Aphi2L*PDstandardNth1alp - 
      alpL*PDstandardNth1Aphi2 + A2xL*PDstandardNth1betax + 
      A2yL*PDstandardNth1betay + A2zL*PDstandardNth1betaz + 
      betayL*PDstandardNth2A2x + betazL*PDstandardNth3A2x) + 
      gu22*(betaxL*PDstandardNth1A2y + betayL*PDstandardNth2A2y - 
      Aphi2L*PDstandardNth2alp - alpL*PDstandardNth2Aphi2 + 
      A2xL*PDstandardNth2betax + A2yL*PDstandardNth2betay + 
      A2zL*PDstandardNth2betaz + betazL*PDstandardNth3A2y) + 
      gu23*(betaxL*PDstandardNth1A2z + betayL*PDstandardNth2A2z + 
      betazL*PDstandardNth3A2z - Aphi2L*PDstandardNth3alp - 
      alpL*PDstandardNth3Aphi2 + A2xL*PDstandardNth3betax + 
      A2yL*PDstandardNth3betay + A2zL*PDstandardNth3betaz))*pow(alpL,-1);
    
    CCTK_REAL Efield2u3 CCTK_ATTRIBUTE_UNUSED = (E2zL*omega + 
      gu13*(betaxL*PDstandardNth1A2x - Aphi2L*PDstandardNth1alp - 
      alpL*PDstandardNth1Aphi2 + A2xL*PDstandardNth1betax + 
      A2yL*PDstandardNth1betay + A2zL*PDstandardNth1betaz + 
      betayL*PDstandardNth2A2x + betazL*PDstandardNth3A2x) + 
      gu23*(betaxL*PDstandardNth1A2y + betayL*PDstandardNth2A2y - 
      Aphi2L*PDstandardNth2alp - alpL*PDstandardNth2Aphi2 + 
      A2xL*PDstandardNth2betax + A2yL*PDstandardNth2betay + 
      A2zL*PDstandardNth2betaz + betazL*PDstandardNth3A2y) + 
      gu33*(betaxL*PDstandardNth1A2z + betayL*PDstandardNth2A2z + 
      betazL*PDstandardNth3A2z - Aphi2L*PDstandardNth3alp - 
      alpL*PDstandardNth3Aphi2 + A2xL*PDstandardNth3betax + 
      A2yL*PDstandardNth3betay + A2zL*PDstandardNth3betaz))*pow(alpL,-1);
    
    E1xL = Efield1u1;
    
    E1yL = Efield1u2;
    
    E1zL = Efield1u3;
    
    E2xL = Efield2u1;
    
    E2yL = Efield2u2;
    
    E2zL = Efield2u3;
    
    CCTK_REAL kxxL CCTK_ATTRIBUTE_UNUSED = At11L;
    
    CCTK_REAL kxyL CCTK_ATTRIBUTE_UNUSED = At12L;
    
    CCTK_REAL kxzL CCTK_ATTRIBUTE_UNUSED = At13L;
    
    CCTK_REAL kyyL CCTK_ATTRIBUTE_UNUSED = At22L;
    
    CCTK_REAL kyzL CCTK_ATTRIBUTE_UNUSED = At23L;
    
    CCTK_REAL kzzL CCTK_ATTRIBUTE_UNUSED = At33L;
    
    At11L = 0;
    
    At12L = 0;
    
    At13L = 0;
    
    At22L = 0;
    
    At23L = 0;
    
    At33L = 0;
    /* Copy local copies back to grid functions */
    At11[index] = At11L;
    At12[index] = At12L;
    At13[index] = At13L;
    At22[index] = At22L;
    At23[index] = At23L;
    At33[index] = At33L;
    E1x[index] = E1xL;
    E1y[index] = E1yL;
    E1z[index] = E1zL;
    E2x[index] = E2xL;
    E2y[index] = E2yL;
    E2z[index] = E2zL;
    kxx[index] = kxxL;
    kxy[index] = kxyL;
    kxz[index] = kxzL;
    kyy[index] = kyyL;
    kyz[index] = kyzL;
    kzz[index] = kzzL;
  }
  CCTK_ENDLOOP3(proca7);
}
extern "C" void proca7(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering proca7_Body");
  }
  if (cctk_iteration % proca7_calc_every != proca7_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::curv",
    "ADMBase::lapse",
    "ADMBase::metric",
    "ADMBase::shift",
    "ML_BSSN::ML_curv",
    "ProcaBase::A1i",
    "ProcaBase::A2i",
    "ProcaBase::Aphi1",
    "ProcaBase::Aphi2",
    "ProcaBase::E1i",
    "ProcaBase::E2i"};
  AssertGroupStorage(cctkGH, "proca7", 11, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "proca7", 1, 1, 1);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "proca7", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "proca7", 3, 3, 3);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "proca7", 4, 4, 4);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, proca7_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving proca7_Body");
  }
}

} // namespace IDProca7
