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


#include "proca_analysis.hh"

#include "carpet.hh"


#ifdef HAVE_CARPET
using namespace Carpet;
#endif


using namespace std;

static inline CCTK_REAL SQUARE(CCTK_REAL x) {return x*x;}

void proca_analysis_point_calc(CCTK_ARGUMENTS)
{
	DECLARE_CCTK_ARGUMENTS
        DECLARE_CCTK_PARAMETERS

	int const ni = cctkGH->cctk_lsh[0];
	int const nj = cctkGH->cctk_lsh[1];
	int const nk = cctkGH->cctk_lsh[2];

/*	#pragma omp parallel
	{*/
		LC_LOOP3(global_quantities, i, j, k, 0, 0, 0, ni, nj, nk, ni, nj, nk) 	     
		{
		
			CCTK_INT i3D = CCTK_GFINDEX3D(cctkGH, i, j, k);

			CCTK_REAL det_metric, sdetg;
			CCTK_REAL guxx, guxy, guxz, guyy, guyz, guzz;
            CCTK_REAL xx, yy, zz, R, rho, costheta, sintheta, rho_z, costheta_y, sintheta_y;//the radial coordinate introduced by Liu et al 2010
            CCTK_REAL cosphi, sinphi, jac21, jac22, jac23, jac31, jac32, phi, phi_y, cosphi_y, sinphi_y, jac31_y, jac32_y;
            
            det_metric = gxx[i3D]*gyy[i3D]*gzz[i3D]-gxx[i3D]*gyz[i3D]*gyz[i3D]-gxy[i3D]*gxy[i3D]*gzz[i3D]+gxy[i3D]*gyz[i3D]*gxz[i3D]+gxz[i3D]*gxy[i3D]*gyz[i3D]-gxz[i3D]*gyy[i3D]*gxz[i3D];
            
            guxx = (-gyz[i3D]*gyz[i3D] + gyy[i3D]*gzz[i3D])/det_metric;
            guxy = ( gxz[i3D]*gyz[i3D] - gxy[i3D]*gzz[i3D])/det_metric;
            guxz = (-gxz[i3D]*gyy[i3D] + gxy[i3D]*gyz[i3D])/det_metric;
            guyy = (-gxz[i3D]*gxz[i3D] + gxx[i3D]*gzz[i3D])/det_metric;
            guyz = ( gxy[i3D]*gxz[i3D] - gxx[i3D]*gyz[i3D])/det_metric;
            guzz = (-gxy[i3D]*gxy[i3D] + gxx[i3D]*gyy[i3D])/det_metric;
            
            sdetg = sqrt(det_metric);
            
            xx = x[i3D];
            yy = y[i3D];
            zz = z[i3D];
            

				//rest mass (D \sqrt{\det{\gamma_ij}})
				Proca_energy_point[i3D] = - alp[i3D]*sdetg*(-1.0/(alp[i3D]*alp[i3D])*eTtt[i3D] - ((guxx - betax[i3D]*betax[i3D]/(alp[i3D]*alp[i3D]))*eTxx[i3D] + 2.0*(guxy - betax[i3D]*betay[i3D]/(alp[i3D]*alp[i3D]))*eTxy[i3D] + 2.0*(guxz - betax[i3D]*betaz[i3D]/(alp[i3D]*alp[i3D]))*eTxz[i3D] + (guyy - betay[i3D]*betay[i3D]/(alp[i3D]*alp[i3D]))*eTyy[i3D] + 2.0*(guyz - betay[i3D]*betaz[i3D]/(alp[i3D]*alp[i3D]))*eTyz[i3D] + (guzz - betaz[i3D]*betaz[i3D]/(alp[i3D]*alp[i3D]))*eTzz[i3D]));
            
       /*     //rest mass (D \sqrt{\det{\gamma_ij}})
            Proca_energy_point[i3D] = - alp[i3D]*sdetg*(-1.0/(alp[i3D]*alp[i3D])*eTtt[i3D] - ((guxx)*eTxx[i3D] + 2.0*(guxy)*eTxy[i3D] + 2.0*(guxz)*eTxz[i3D] + (guyy)*eTyy[i3D] + 2.0*(guyz)*eTyz[i3D] + (guzz)*eTzz[i3D]));
        
        rho = n^mu n^nu T_{mu nu}
        !     = (T_{00} - 2 beta^i T_{i0} + beta^i beta^j T_{ij})/alph^2
*/
				//grav mass
               R = sqrt(SQUARE(xx) + SQUARE(yy) + SQUARE(zz));
               rho = sqrt(SQUARE(xx) + SQUARE(yy));
               rho_z = sqrt(SQUARE(xx) + SQUARE(zz));
            
            if(rho<1e-10)
            
            {
                xx = 1e-10;
                rho = sqrt(SQUARE(xx) + SQUARE(yy));
                rho_z = sqrt(SQUARE(xx) + SQUARE(zz));
                R = sqrt(SQUARE(xx) + SQUARE(yy) + SQUARE(zz));
            }
            
               sintheta = rho / R;
               costheta = zz / R;

               sintheta_y = rho_z / R;
               costheta_y = yy / R;
            
            phi = atan(yy/xx);

            phi_y = atan(zz/xx);
            
            if((x[i3D]<0.0)){
                phi = atan(yy/xx) + 3.14159265358979323846;
            }

            if((x[i3D]<0.0)){
                phi_y = atan(zz/xx) + 3.14159265358979323846;
            }
            
            cosphi = cos(phi);
            sinphi = sin(phi);

            cosphi_y = cos(phi_y);
            sinphi_y = sin(phi_y);

            jac21 = + R*costheta_y*cosphi_y;
            jac22 = + R*costheta_y*sinphi_y;
            jac23 = - R*sintheta_y;
            
            jac31 = - R*sintheta*sinphi;
            jac32 = + R*sintheta*cosphi;

            jac31_y = - R*sintheta_y*sinphi_y;
            jac32_y = + R*sintheta_y*cosphi_y;
            
            Proca_angmomentum_point[i3D] = sdetg*((-1/alp[i3D]*(eTtx[i3D] - betax[i3D]*eTxx[i3D] - betay[i3D]*eTxy[i3D] - betaz[i3D]*eTxz[i3D])*jac31) 
						+ (-1/alp[i3D]*(eTty[i3D] - betax[i3D]*eTxy[i3D] - betay[i3D]*eTyy[i3D] - betaz[i3D]*eTyz[i3D])*jac32));

            Proca_angmomentum_point_y[i3D] = sdetg*((-1/alp[i3D]*(eTtx[i3D] - betax[i3D]*eTxx[i3D] - betay[i3D]*eTxy[i3D] - betaz[i3D]*eTxz[i3D])*jac31_y) 
						  + (-1/alp[i3D]*(eTtz[i3D] - betax[i3D]*eTxz[i3D] - betay[i3D]*eTyz[i3D] - betaz[i3D]*eTzz[i3D])*jac32_y));
            
				Proca_rho_point[i3D] =  1.0/(alp[i3D]*alp[i3D])*sdetg*(eTtt[i3D]
                                                                 - 2.0*betax[i3D]*eTtx[i3D]
                                                                 - 2.0*betay[i3D]*eTty[i3D]
                                                                 - 2.0*betaz[i3D]*eTtz[i3D]
                                                                 + betax[i3D]*betax[i3D]*eTxx[i3D]
                                                                 + betay[i3D]*betay[i3D]*eTyy[i3D]
                                                                 + betaz[i3D]*betaz[i3D]*eTzz[i3D]
                                                                 + 2.0*betax[i3D]*betay[i3D]*eTxy[i3D]
                                                                 + 2.0*betax[i3D]*betaz[i3D]*eTxz[i3D]
                                                                 + 2.0*betay[i3D]*betaz[i3D]*eTyz[i3D]);


        Proca_rho_point_r30[i3D] = 0.0;
        Proca_energy_point_r30[i3D] = 0.0;
        Proca_angmomentum_point_r30[i3D] = 0.0;

        Proca_rho_point_r60[i3D] = 0.0;
        Proca_energy_point_r60[i3D] = 0.0;
        Proca_angmomentum_point_r60[i3D] = 0.0;
		
	if (R<150.0)
	{
                if (R<30.0)
                {
                Proca_rho_point_r30[i3D] = Proca_rho_point[i3D];

                Proca_energy_point_r30[i3D] = Proca_energy_point[i3D];

                Proca_angmomentum_point_r30[i3D] = Proca_angmomentum_point[i3D];
                }
                if (R<60.0)
                {
                Proca_rho_point_r60[i3D] = Proca_rho_point[i3D];

                Proca_energy_point_r60[i3D] = Proca_energy_point[i3D];

                Proca_angmomentum_point_r60[i3D] = Proca_angmomentum_point[i3D];
                }
		if (R<100.0)
		{
	
		Proca_rho_point_r100[i3D] = Proca_rho_point[i3D];
                Proca_rho_point_r150[i3D] = Proca_rho_point[i3D];

                Proca_energy_point_r100[i3D] = Proca_energy_point[i3D];
                Proca_energy_point_r150[i3D] = Proca_energy_point[i3D];

                Proca_angmomentum_point_r100[i3D] = Proca_angmomentum_point[i3D];
                Proca_angmomentum_point_r150[i3D] = Proca_angmomentum_point[i3D];


		}else {

                Proca_rho_point_r150[i3D] = Proca_rho_point[i3D];
		Proca_rho_point_r100[i3D] = 0.0;

                Proca_energy_point_r150[i3D] = Proca_energy_point[i3D];
                Proca_energy_point_r100[i3D] = 0.0;

                Proca_angmomentum_point_r150[i3D] = Proca_angmomentum_point[i3D];
                Proca_angmomentum_point_r100[i3D] = 0.0;

		}
	
	}else {

                Proca_rho_point_r100[i3D] = 0.0;
                Proca_rho_point_r150[i3D] = 0.0;

                Proca_energy_point_r100[i3D] = 0.0;
                Proca_energy_point_r150[i3D] = 0.0;

                Proca_angmomentum_point_r100[i3D] = 0.0;
                Proca_angmomentum_point_r150[i3D] = 0.0;

	}



		}LC_ENDLOOP3(global_quantities);

/*	}*/
	
}


void proca_analysis_reduce(CCTK_ARGUMENTS)
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
                    Proca_energy, 1,
                    CCTK_VarIndex("ProcaAnalysis::Proca_energy_point")))
        	CCTK_WARN(0, "Error while reducing Proca energy mass");

        if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    Proca_rho, 1,
                    CCTK_VarIndex("ProcaAnalysis::Proca_rho_point")))
        	CCTK_WARN(0, "Error while reducing Proca rho energy");
    
    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                   CCTK_VARIABLE_REAL,
                   Proca_angularmomentum, 1,
                   CCTK_VarIndex("ProcaAnalysis::Proca_angmomentum_point")))
    CCTK_WARN(0, "Error while reducing Proca angular momentum");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                   CCTK_VARIABLE_REAL,
                   Proca_angularmomentum_y, 1,
                   CCTK_VarIndex("ProcaAnalysis::Proca_angmomentum_point_y")))
    CCTK_WARN(0, "Error while reducing Proca angular momentum y");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    Proca_energy_r30, 1,
                    CCTK_VarIndex("ProcaAnalysis::Proca_energy_point_r30")))
                CCTK_WARN(0, "Error while reducing Proca energy mass");

        if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    Proca_rho_r30, 1,
                    CCTK_VarIndex("ProcaAnalysis::Proca_rho_point_r30")))
                CCTK_WARN(0, "Error while reducing Proca rho energy r30");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                   CCTK_VARIABLE_REAL,
                   Proca_angularmomentum_r30, 1,
                   CCTK_VarIndex("ProcaAnalysis::Proca_angmomentum_point_r30")))
    CCTK_WARN(0, "Error while reducing Proca angular momentum r30");


    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    Proca_energy_r60, 1,
                    CCTK_VarIndex("ProcaAnalysis::Proca_energy_point_r60")))
                CCTK_WARN(0, "Error while reducing Proca energy mass");

        if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    Proca_rho_r60, 1,
                    CCTK_VarIndex("ProcaAnalysis::Proca_rho_point_r60")))
                CCTK_WARN(0, "Error while reducing Proca rho energy r60");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                   CCTK_VARIABLE_REAL,
                   Proca_angularmomentum_r60, 1,
                   CCTK_VarIndex("ProcaAnalysis::Proca_angmomentum_point_r60")))
    CCTK_WARN(0, "Error while reducing Proca angular momentum r60");

        if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    Proca_energy_r100, 1,
                    CCTK_VarIndex("ProcaAnalysis::Proca_energy_point_r100")))
                CCTK_WARN(0, "Error while reducing Proca energy mass");

        if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    Proca_rho_r100, 1,
                    CCTK_VarIndex("ProcaAnalysis::Proca_rho_point_r100")))
                CCTK_WARN(0, "Error while reducing Proca rho energy r100");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                   CCTK_VARIABLE_REAL,
                   Proca_angularmomentum_r100, 1,
                   CCTK_VarIndex("ProcaAnalysis::Proca_angmomentum_point_r100")))
    CCTK_WARN(0, "Error while reducing Proca angular momentum r100");


        if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    Proca_energy_r150, 1,
                    CCTK_VarIndex("ProcaAnalysis::Proca_energy_point_r150")))
                CCTK_WARN(0, "Error while reducing Proca energy mass r150");

        if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                    CCTK_VARIABLE_REAL,
                    Proca_rho_r150, 1,
                    CCTK_VarIndex("ProcaAnalysis::Proca_rho_point_r150")))
                CCTK_WARN(0, "Error while reducing Proca rho energy r150");

    if (CCTK_Reduce(cctkGH, -1, reduction_handle, 1,
                   CCTK_VARIABLE_REAL,
                   Proca_angularmomentum_r150, 1,
                   CCTK_VarIndex("ProcaAnalysis::Proca_angmomentum_point_r150")))
    CCTK_WARN(0, "Error while reducing Proca angular momentum r150");


	*Proca_energy *= coarse_grid_coord_volume;
    *Proca_rho *= coarse_grid_coord_volume;
	*Proca_angularmomentum *= coarse_grid_coord_volume;
        *Proca_angularmomentum_y *= coarse_grid_coord_volume;

    *Proca_energy_r30 *= coarse_grid_coord_volume;
    *Proca_rho_r30 *= coarse_grid_coord_volume;
        *Proca_angularmomentum_r30 *= coarse_grid_coord_volume;

    *Proca_energy_r60 *= coarse_grid_coord_volume;
    *Proca_rho_r60 *= coarse_grid_coord_volume;
        *Proca_angularmomentum_r60 *= coarse_grid_coord_volume;

        *Proca_energy_r100 *= coarse_grid_coord_volume;
    *Proca_rho_r100 *= coarse_grid_coord_volume;
        *Proca_angularmomentum_r100 *= coarse_grid_coord_volume;

        *Proca_energy_r150 *= coarse_grid_coord_volume;
    *Proca_rho_r150 *= coarse_grid_coord_volume;
        *Proca_angularmomentum_r150 *= coarse_grid_coord_volume;



}


