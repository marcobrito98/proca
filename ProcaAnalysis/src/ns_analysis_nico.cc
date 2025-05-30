#include <string.h>
#include <assert.h>
#include <string>
#include <iomanip>
#include <math.h>
#include <vector>
#include <iostream>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

#include "loopcontrol.h"
#include "cctk_Loop.h"
#include "CoordBase.h"

#include <CarpetReduce_bits.h>


#include "ns_analysis.hh"

#include "carpet.hh"


#ifdef HAVE_CARPET
using namespace Carpet;
#endif


using namespace std;

void ns_analysis_point_calc(CCTK_ARGUMENTS)
{
     DECLARE_CCTK_ARGUMENTS
     DECLARE_CCTK_PARAMETERS

     int const ni = cctkGH->cctk_lsh[0];
     int const nj = cctkGH->cctk_lsh[1];
     int const nk = cctkGH->cctk_lsh[2];
       
     #pragma omp parallel
     {
     	LC_LOOP3(global_quantities, i ,j, k, 0, 0 ,0 ,ni, nj ,nk ,ni, nj, nk)
     	{
	
 		CCTK_INT i3D = CCTK_GFINDEX3(cctkGH, i, j, k);
		CCTK_INT idx1 = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,0);
                CCTK_INT idx2 = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,1);
                CCTK_INT idx3 = CCTK_VECTGFINDEX3D(cctkGH, i,j,k,2);
                
		CCTK_REAL h;
		CCTK_REAL u0, ul0;
		CCTK_REAL T0l0, Tmulmu;
		
		
		
		if(rho[i3D] <= (rho_atmo*10.0))
		{
		
			M0_point[i3D] = 0.0;
			M_point[i3D]  = 0.0;

		}
		else
		{
                        h = 1.0 + eps[i3D] + press[i3D] / rho[i3D];

                        u0 = w_lorentz[i3D] / alp[i3D];

			ul0 = w_lorentz[i3D] * ( -alp[i3D]
   				 + gxx[i3D] * vel[idx1] * betax[i3D] + gxy[i3D] * vel[idx1] * betay[i3D] + gxz[i3D] * vel[idx1] * betaz[i3D]
			         + gxy[i3D] * vel[idx2] * betax[i3D] + gyy[i3D] * vel[idx2] * betay[i3D] + gyz[i3D] * vel[idx2] * betaz[i3D]
			         + gxz[i3D] * vel[idx3] * betax[i3D] + gyz[i3D] * vel[idx3] * betay[i3D] + gzz[i3D] * vel[idx3] * betaz[i3D]
			    );

			T0l0 = rho[i3D] * h * u0 * ul0 + press[i3D];
			Tmulmu = - rho[i3D] * (1.0 + eps[i3D]) + 3.0 * press[i3D];			
			
			//rest mass (D \sqrt{\det{\gamma_ij}}
			M0_point[i3D] = sdetg[i3D] * rho[i3D] * w_lorentz[i3D]; 

			//gravitational mass
			M_point[i3D]  = (-2.0 * T0l0 + Tmulmu) * alp[i3D] * sdetg[i3D];
		}		
	
     	}LC_ENDLOOP3(global_quantities);
     
     }
}

void ns_analysis_reduce(CCTK_ARGUMENTS)
{
	DECLARE_CCTK_ARGUMENTS
	DECLARE_CCTK_PARAMETERS

	CCTK_INT reduction_handle;

        CCTK_REAL coarse_grid_coord_volume = cctk_delta_space[0] * cctk_delta_space[1] * cctk_delta_space[2];

        char *x_reflection, *y_reflection, *z_reflection;


        if(CCTK_IsThornActive("ReflectionSymmetry")  !=0 ) 
        {
                        x_reflection=CCTK_ParameterValString("reflection_x","ReflectionSymmetry");
                        y_reflection=CCTK_ParameterValString("reflection_y","ReflectionSymmetry");
                        z_reflection=CCTK_ParameterValString("reflection_z","ReflectionSymmetry");
                        if (CCTK_Equals(x_reflection,"yes") == 1 ) 
                        {
                                coarse_grid_coord_volume *= 2.0;
                        }
                        if (CCTK_Equals(y_reflection,"yes") == 1 ) 
                        {
                                coarse_grid_coord_volume *= 2.0;
                        }
                        if (CCTK_Equals(z_reflection,"yes") == 1 ) 
                        {
                                coarse_grid_coord_volume *= 2.0;
                        }
         }

         if(CCTK_IsThornActive("RotatingSymmetry180") !=0 ) 
         {
                        coarse_grid_coord_volume *= 2.0;
         }

         if(CCTK_IsThornActive("RotatingSymmetry90") !=0 ) 
         {
                        coarse_grid_coord_volume *= 4.0;
         }

         reduction_handle = CCTK_ReductionHandle("sum");
         if (reduction_handle < 0)
         	CCTK_WARN(0, "Unable to get reduction handle.");

	 if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    NS_rest_mass, 1,
                    CCTK_VarIndex("NS_Analysis::M0_point")))
                CCTK_WARN(0, "Error while reducing NS rest mass");

         if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    NS_gravitational_mass, 1,
                    CCTK_VarIndex("NS_Analysis::M_point")))
                CCTK_WARN(0, "Error while reducing NS gravitational mass");

	 *NS_rest_mass *= coarse_grid_coord_volume;
	 *NS_gravitational_mass *= coarse_grid_coord_volume;

}
