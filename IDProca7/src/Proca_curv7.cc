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

extern "C" void Proca_curv7_SelectBCs(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (cctk_iteration % Proca_curv7_calc_every != Proca_curv7_calc_offset)
    return;
  CCTK_INT ierr CCTK_ATTRIBUTE_UNUSED = 0;
  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, GetBoundaryWidth(cctkGH), -1 /* no table */, "ML_BSSN::ML_curv","flat");
  if (ierr < 0)
    CCTK_WARN(CCTK_WARN_ALERT, "Failed to register flat BC for ML_BSSN::ML_curv.");
  return;
}

static void Proca_curv7_Body(const cGH* restrict const cctkGH, const int dir, const int face, const CCTK_REAL normal[3], const CCTK_REAL tangentA[3], const CCTK_REAL tangentB[3], const int imin[3], const int imax[3], const int n_subblock_gfs, CCTK_REAL* restrict const subblock_gfs[])
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
  CCTK_LOOP3(Proca_curv7,
    i,j,k, imin0,imin1,imin2, imax0,imax1,imax2,
    cctk_ash[0],cctk_ash[1],cctk_ash[2])
  {
    const ptrdiff_t index CCTK_ATTRIBUTE_UNUSED = di*i + dj*j + dk*k;
    /* Assign local copies of grid functions */
    
    CCTK_REAL alpL CCTK_ATTRIBUTE_UNUSED = alp[index];
    CCTK_REAL betaxL CCTK_ATTRIBUTE_UNUSED = betax[index];
    CCTK_REAL betayL CCTK_ATTRIBUTE_UNUSED = betay[index];
    CCTK_REAL betazL CCTK_ATTRIBUTE_UNUSED = betaz[index];
    CCTK_REAL gxxL CCTK_ATTRIBUTE_UNUSED = gxx[index];
    CCTK_REAL gxyL CCTK_ATTRIBUTE_UNUSED = gxy[index];
    CCTK_REAL gxzL CCTK_ATTRIBUTE_UNUSED = gxz[index];
    CCTK_REAL gyyL CCTK_ATTRIBUTE_UNUSED = gyy[index];
    CCTK_REAL gyzL CCTK_ATTRIBUTE_UNUSED = gyz[index];
    CCTK_REAL gzzL CCTK_ATTRIBUTE_UNUSED = gzz[index];
    CCTK_REAL kxxL CCTK_ATTRIBUTE_UNUSED = kxx[index];
    CCTK_REAL kxyL CCTK_ATTRIBUTE_UNUSED = kxy[index];
    CCTK_REAL kxzL CCTK_ATTRIBUTE_UNUSED = kxz[index];
    CCTK_REAL kyyL CCTK_ATTRIBUTE_UNUSED = kyy[index];
    CCTK_REAL kyzL CCTK_ATTRIBUTE_UNUSED = kyz[index];
    CCTK_REAL kzzL CCTK_ATTRIBUTE_UNUSED = kzz[index];
    
    /* Include user supplied include files */
    /* Precompute derivatives */
    CCTK_REAL PDstandardNth1betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3betax CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3betay CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3betaz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth1gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth3gzz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2kxx CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2kxy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2kxz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2kyy CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2kyz CCTK_ATTRIBUTE_UNUSED;
    CCTK_REAL PDstandardNth2kzz CCTK_ATTRIBUTE_UNUSED;
    
    switch (fdOrder)
    {
      case 2:
      {
        PDstandardNth1betax = PDstandardNthfdOrder21(&betax[index]);
        PDstandardNth2betax = PDstandardNthfdOrder22(&betax[index]);
        PDstandardNth3betax = PDstandardNthfdOrder23(&betax[index]);
        PDstandardNth1betay = PDstandardNthfdOrder21(&betay[index]);
        PDstandardNth2betay = PDstandardNthfdOrder22(&betay[index]);
        PDstandardNth3betay = PDstandardNthfdOrder23(&betay[index]);
        PDstandardNth1betaz = PDstandardNthfdOrder21(&betaz[index]);
        PDstandardNth2betaz = PDstandardNthfdOrder22(&betaz[index]);
        PDstandardNth3betaz = PDstandardNthfdOrder23(&betaz[index]);
        PDstandardNth1gxx = PDstandardNthfdOrder21(&gxx[index]);
        PDstandardNth2gxx = PDstandardNthfdOrder22(&gxx[index]);
        PDstandardNth3gxx = PDstandardNthfdOrder23(&gxx[index]);
        PDstandardNth1gxy = PDstandardNthfdOrder21(&gxy[index]);
        PDstandardNth2gxy = PDstandardNthfdOrder22(&gxy[index]);
        PDstandardNth3gxy = PDstandardNthfdOrder23(&gxy[index]);
        PDstandardNth1gxz = PDstandardNthfdOrder21(&gxz[index]);
        PDstandardNth2gxz = PDstandardNthfdOrder22(&gxz[index]);
        PDstandardNth3gxz = PDstandardNthfdOrder23(&gxz[index]);
        PDstandardNth1gyy = PDstandardNthfdOrder21(&gyy[index]);
        PDstandardNth2gyy = PDstandardNthfdOrder22(&gyy[index]);
        PDstandardNth3gyy = PDstandardNthfdOrder23(&gyy[index]);
        PDstandardNth1gyz = PDstandardNthfdOrder21(&gyz[index]);
        PDstandardNth2gyz = PDstandardNthfdOrder22(&gyz[index]);
        PDstandardNth3gyz = PDstandardNthfdOrder23(&gyz[index]);
        PDstandardNth1gzz = PDstandardNthfdOrder21(&gzz[index]);
        PDstandardNth2gzz = PDstandardNthfdOrder22(&gzz[index]);
        PDstandardNth3gzz = PDstandardNthfdOrder23(&gzz[index]);
        PDstandardNth2kxx = PDstandardNthfdOrder22(&kxx[index]);
        PDstandardNth2kxy = PDstandardNthfdOrder22(&kxy[index]);
        PDstandardNth2kxz = PDstandardNthfdOrder22(&kxz[index]);
        PDstandardNth2kyy = PDstandardNthfdOrder22(&kyy[index]);
        PDstandardNth2kyz = PDstandardNthfdOrder22(&kyz[index]);
        PDstandardNth2kzz = PDstandardNthfdOrder22(&kzz[index]);
        break;
      }
      
      case 4:
      {
        PDstandardNth1betax = PDstandardNthfdOrder41(&betax[index]);
        PDstandardNth2betax = PDstandardNthfdOrder42(&betax[index]);
        PDstandardNth3betax = PDstandardNthfdOrder43(&betax[index]);
        PDstandardNth1betay = PDstandardNthfdOrder41(&betay[index]);
        PDstandardNth2betay = PDstandardNthfdOrder42(&betay[index]);
        PDstandardNth3betay = PDstandardNthfdOrder43(&betay[index]);
        PDstandardNth1betaz = PDstandardNthfdOrder41(&betaz[index]);
        PDstandardNth2betaz = PDstandardNthfdOrder42(&betaz[index]);
        PDstandardNth3betaz = PDstandardNthfdOrder43(&betaz[index]);
        PDstandardNth1gxx = PDstandardNthfdOrder41(&gxx[index]);
        PDstandardNth2gxx = PDstandardNthfdOrder42(&gxx[index]);
        PDstandardNth3gxx = PDstandardNthfdOrder43(&gxx[index]);
        PDstandardNth1gxy = PDstandardNthfdOrder41(&gxy[index]);
        PDstandardNth2gxy = PDstandardNthfdOrder42(&gxy[index]);
        PDstandardNth3gxy = PDstandardNthfdOrder43(&gxy[index]);
        PDstandardNth1gxz = PDstandardNthfdOrder41(&gxz[index]);
        PDstandardNth2gxz = PDstandardNthfdOrder42(&gxz[index]);
        PDstandardNth3gxz = PDstandardNthfdOrder43(&gxz[index]);
        PDstandardNth1gyy = PDstandardNthfdOrder41(&gyy[index]);
        PDstandardNth2gyy = PDstandardNthfdOrder42(&gyy[index]);
        PDstandardNth3gyy = PDstandardNthfdOrder43(&gyy[index]);
        PDstandardNth1gyz = PDstandardNthfdOrder41(&gyz[index]);
        PDstandardNth2gyz = PDstandardNthfdOrder42(&gyz[index]);
        PDstandardNth3gyz = PDstandardNthfdOrder43(&gyz[index]);
        PDstandardNth1gzz = PDstandardNthfdOrder41(&gzz[index]);
        PDstandardNth2gzz = PDstandardNthfdOrder42(&gzz[index]);
        PDstandardNth3gzz = PDstandardNthfdOrder43(&gzz[index]);
        PDstandardNth2kxx = PDstandardNthfdOrder42(&kxx[index]);
        PDstandardNth2kxy = PDstandardNthfdOrder42(&kxy[index]);
        PDstandardNth2kxz = PDstandardNthfdOrder42(&kxz[index]);
        PDstandardNth2kyy = PDstandardNthfdOrder42(&kyy[index]);
        PDstandardNth2kyz = PDstandardNthfdOrder42(&kyz[index]);
        PDstandardNth2kzz = PDstandardNthfdOrder42(&kzz[index]);
        break;
      }
      
      case 6:
      {
        PDstandardNth1betax = PDstandardNthfdOrder61(&betax[index]);
        PDstandardNth2betax = PDstandardNthfdOrder62(&betax[index]);
        PDstandardNth3betax = PDstandardNthfdOrder63(&betax[index]);
        PDstandardNth1betay = PDstandardNthfdOrder61(&betay[index]);
        PDstandardNth2betay = PDstandardNthfdOrder62(&betay[index]);
        PDstandardNth3betay = PDstandardNthfdOrder63(&betay[index]);
        PDstandardNth1betaz = PDstandardNthfdOrder61(&betaz[index]);
        PDstandardNth2betaz = PDstandardNthfdOrder62(&betaz[index]);
        PDstandardNth3betaz = PDstandardNthfdOrder63(&betaz[index]);
        PDstandardNth1gxx = PDstandardNthfdOrder61(&gxx[index]);
        PDstandardNth2gxx = PDstandardNthfdOrder62(&gxx[index]);
        PDstandardNth3gxx = PDstandardNthfdOrder63(&gxx[index]);
        PDstandardNth1gxy = PDstandardNthfdOrder61(&gxy[index]);
        PDstandardNth2gxy = PDstandardNthfdOrder62(&gxy[index]);
        PDstandardNth3gxy = PDstandardNthfdOrder63(&gxy[index]);
        PDstandardNth1gxz = PDstandardNthfdOrder61(&gxz[index]);
        PDstandardNth2gxz = PDstandardNthfdOrder62(&gxz[index]);
        PDstandardNth3gxz = PDstandardNthfdOrder63(&gxz[index]);
        PDstandardNth1gyy = PDstandardNthfdOrder61(&gyy[index]);
        PDstandardNth2gyy = PDstandardNthfdOrder62(&gyy[index]);
        PDstandardNth3gyy = PDstandardNthfdOrder63(&gyy[index]);
        PDstandardNth1gyz = PDstandardNthfdOrder61(&gyz[index]);
        PDstandardNth2gyz = PDstandardNthfdOrder62(&gyz[index]);
        PDstandardNth3gyz = PDstandardNthfdOrder63(&gyz[index]);
        PDstandardNth1gzz = PDstandardNthfdOrder61(&gzz[index]);
        PDstandardNth2gzz = PDstandardNthfdOrder62(&gzz[index]);
        PDstandardNth3gzz = PDstandardNthfdOrder63(&gzz[index]);
        PDstandardNth2kxx = PDstandardNthfdOrder62(&kxx[index]);
        PDstandardNth2kxy = PDstandardNthfdOrder62(&kxy[index]);
        PDstandardNth2kxz = PDstandardNthfdOrder62(&kxz[index]);
        PDstandardNth2kyy = PDstandardNthfdOrder62(&kyy[index]);
        PDstandardNth2kyz = PDstandardNthfdOrder62(&kyz[index]);
        PDstandardNth2kzz = PDstandardNthfdOrder62(&kzz[index]);
        break;
      }
      
      case 8:
      {
        PDstandardNth1betax = PDstandardNthfdOrder81(&betax[index]);
        PDstandardNth2betax = PDstandardNthfdOrder82(&betax[index]);
        PDstandardNth3betax = PDstandardNthfdOrder83(&betax[index]);
        PDstandardNth1betay = PDstandardNthfdOrder81(&betay[index]);
        PDstandardNth2betay = PDstandardNthfdOrder82(&betay[index]);
        PDstandardNth3betay = PDstandardNthfdOrder83(&betay[index]);
        PDstandardNth1betaz = PDstandardNthfdOrder81(&betaz[index]);
        PDstandardNth2betaz = PDstandardNthfdOrder82(&betaz[index]);
        PDstandardNth3betaz = PDstandardNthfdOrder83(&betaz[index]);
        PDstandardNth1gxx = PDstandardNthfdOrder81(&gxx[index]);
        PDstandardNth2gxx = PDstandardNthfdOrder82(&gxx[index]);
        PDstandardNth3gxx = PDstandardNthfdOrder83(&gxx[index]);
        PDstandardNth1gxy = PDstandardNthfdOrder81(&gxy[index]);
        PDstandardNth2gxy = PDstandardNthfdOrder82(&gxy[index]);
        PDstandardNth3gxy = PDstandardNthfdOrder83(&gxy[index]);
        PDstandardNth1gxz = PDstandardNthfdOrder81(&gxz[index]);
        PDstandardNth2gxz = PDstandardNthfdOrder82(&gxz[index]);
        PDstandardNth3gxz = PDstandardNthfdOrder83(&gxz[index]);
        PDstandardNth1gyy = PDstandardNthfdOrder81(&gyy[index]);
        PDstandardNth2gyy = PDstandardNthfdOrder82(&gyy[index]);
        PDstandardNth3gyy = PDstandardNthfdOrder83(&gyy[index]);
        PDstandardNth1gyz = PDstandardNthfdOrder81(&gyz[index]);
        PDstandardNth2gyz = PDstandardNthfdOrder82(&gyz[index]);
        PDstandardNth3gyz = PDstandardNthfdOrder83(&gyz[index]);
        PDstandardNth1gzz = PDstandardNthfdOrder81(&gzz[index]);
        PDstandardNth2gzz = PDstandardNthfdOrder82(&gzz[index]);
        PDstandardNth3gzz = PDstandardNthfdOrder83(&gzz[index]);
        PDstandardNth2kxx = PDstandardNthfdOrder82(&kxx[index]);
        PDstandardNth2kxy = PDstandardNthfdOrder82(&kxy[index]);
        PDstandardNth2kxz = PDstandardNthfdOrder82(&kxz[index]);
        PDstandardNth2kyy = PDstandardNthfdOrder82(&kyy[index]);
        PDstandardNth2kyz = PDstandardNthfdOrder82(&kyz[index]);
        PDstandardNth2kzz = PDstandardNthfdOrder82(&kzz[index]);
        break;
      }
      default:
        CCTK_BUILTIN_UNREACHABLE();
    }
    /* Calculate temporaries and grid functions */
    CCTK_REAL K11 CCTK_ATTRIBUTE_UNUSED = 0.5*(2*(gxxL*PDstandardNth1betax 
      + gxyL*PDstandardNth1betay + gxzL*PDstandardNth1betaz) + 
      betaxL*PDstandardNth1gxx + betayL*PDstandardNth2gxx - PDstandardNth2kxx 
      + betazL*PDstandardNth3gxx)*pow(alpL,-1);
    
    CCTK_REAL K12 CCTK_ATTRIBUTE_UNUSED = 0.5*(gyyL*PDstandardNth1betay + 
      gyzL*PDstandardNth1betaz + betaxL*PDstandardNth1gxy + 
      gxxL*PDstandardNth2betax + gxyL*(PDstandardNth1betax + 
      PDstandardNth2betay) + gxzL*PDstandardNth2betaz + 
      betayL*PDstandardNth2gxy - PDstandardNth2kxy + 
      betazL*PDstandardNth3gxy)*pow(alpL,-1);
    
    CCTK_REAL K13 CCTK_ATTRIBUTE_UNUSED = 0.5*(gyzL*PDstandardNth1betay + 
      gzzL*PDstandardNth1betaz + betaxL*PDstandardNth1gxz + 
      betayL*PDstandardNth2gxz - PDstandardNth2kxz + gxxL*PDstandardNth3betax 
      + gxyL*PDstandardNth3betay + gxzL*(PDstandardNth1betax + 
      PDstandardNth3betaz) + betazL*PDstandardNth3gxz)*pow(alpL,-1);
    
    CCTK_REAL K22 CCTK_ATTRIBUTE_UNUSED = 0.5*(betaxL*PDstandardNth1gyy + 
      2*(gxyL*PDstandardNth2betax + gyyL*PDstandardNth2betay + 
      gyzL*PDstandardNth2betaz) + betayL*PDstandardNth2gyy - 
      PDstandardNth2kyy + betazL*PDstandardNth3gyy)*pow(alpL,-1);
    
    CCTK_REAL K23 CCTK_ATTRIBUTE_UNUSED = 0.5*(betaxL*PDstandardNth1gyz + 
      gxzL*PDstandardNth2betax + gzzL*PDstandardNth2betaz + 
      betayL*PDstandardNth2gyz - PDstandardNth2kyz + gxyL*PDstandardNth3betax 
      + gyyL*PDstandardNth3betay + gyzL*(PDstandardNth2betay + 
      PDstandardNth3betaz) + betazL*PDstandardNth3gyz)*pow(alpL,-1);
    
    CCTK_REAL K33 CCTK_ATTRIBUTE_UNUSED = 0.5*(betaxL*PDstandardNth1gzz + 
      betayL*PDstandardNth2gzz - PDstandardNth2kzz + 
      2*(gxzL*PDstandardNth3betax + gyzL*PDstandardNth3betay + 
      gzzL*PDstandardNth3betaz) + betazL*PDstandardNth3gzz)*pow(alpL,-1);
    
    CCTK_REAL At11L CCTK_ATTRIBUTE_UNUSED = K11;
    
    CCTK_REAL At12L CCTK_ATTRIBUTE_UNUSED = K12;
    
    CCTK_REAL At13L CCTK_ATTRIBUTE_UNUSED = K13;
    
    CCTK_REAL At22L CCTK_ATTRIBUTE_UNUSED = K22;
    
    CCTK_REAL At23L CCTK_ATTRIBUTE_UNUSED = K23;
    
    CCTK_REAL At33L CCTK_ATTRIBUTE_UNUSED = K33;
    /* Copy local copies back to grid functions */
    At11[index] = At11L;
    At12[index] = At12L;
    At13[index] = At13L;
    At22[index] = At22L;
    At23[index] = At23L;
    At33[index] = At33L;
  }
  CCTK_ENDLOOP3(Proca_curv7);
}
extern "C" void Proca_curv7(CCTK_ARGUMENTS)
{
  DECLARE_CCTK_ARGUMENTS;
  DECLARE_CCTK_PARAMETERS;
  
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Entering Proca_curv7_Body");
  }
  if (cctk_iteration % Proca_curv7_calc_every != Proca_curv7_calc_offset)
  {
    return;
  }
  
  const char* const groups[] = {
    "ADMBase::curv",
    "ADMBase::lapse",
    "ADMBase::metric",
    "ADMBase::shift",
    "ML_BSSN::ML_curv"};
  AssertGroupStorage(cctkGH, "Proca_curv7", 5, groups);
  
  switch (fdOrder)
  {
    case 2:
    {
      EnsureStencilFits(cctkGH, "Proca_curv7", 1, 1, 1);
      break;
    }
    
    case 4:
    {
      EnsureStencilFits(cctkGH, "Proca_curv7", 2, 2, 2);
      break;
    }
    
    case 6:
    {
      EnsureStencilFits(cctkGH, "Proca_curv7", 3, 3, 3);
      break;
    }
    
    case 8:
    {
      EnsureStencilFits(cctkGH, "Proca_curv7", 4, 4, 4);
      break;
    }
    default:
      CCTK_BUILTIN_UNREACHABLE();
  }
  
  LoopOverInterior(cctkGH, Proca_curv7_Body);
  if (verbose > 1)
  {
    CCTK_VInfo(CCTK_THORNSTRING,"Leaving Proca_curv7_Body");
  }
}

} // namespace IDProca7
