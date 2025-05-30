#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <fstream>

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"

#include "loopcontrol.h"

#include "proca.hh"

using namespace std;


static inline CCTK_REAL SQUARE(CCTK_REAL x) {return x*x;}


void IDProca6_setup(CCTK_ARGUMENTS)
{
    DECLARE_CCTK_ARGUMENTS;
    DECLARE_CCTK_PARAMETERS;
    
    
    int const ni = cctkGH->cctk_lsh[0];
    int const nj = cctkGH->cctk_lsh[1];
    int const nk = cctkGH->cctk_lsh[2];
    
    /*CCTK_REAL const R22 = cos(rotation_angle);
     CCTK_REAL const R23 = -sin(rotation_angle);
     CCTK_REAL const R32 = sin(rotation_angle);
     CCTK_REAL const R33 = cos(rotation_angle);*/
    
    
    double a;
    double data[num_var][num_r][num_th];
    CCTK_REAL sh1[num_r];
    CCTK_REAL th1[num_th];


  CCTK_REAL *F1_out, *F2_out, *F0_out, *H1_out, *H2_out, *Vproca_out;
    
  F1_out    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  F2_out    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  F0_out    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  H1_out    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  H2_out    = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
  Vproca_out     = (CCTK_REAL *) malloc(maxNF * sizeof(CCTK_REAL));
    
    ifstream infile(filename);
    
    CCTK_VInfo(CCTK_THORNSTRING,"Reading file");
    
    for(int k = 0; k <num_th; k++){
        for(int j = 0; j <num_r; j++){
            for(int i = 0; i <num_var; i++){
                infile >> a;
                data[i][j][k]=a;
            }
        }
    }
    
    
    for(int k = 0; k <num_th; k++){
        for(int j = 0; j <num_r; j++){
            sh1[j] = data[0][j][0];
            th1[k] = data[1][0][k];
            F1_out[num_r*k+j] = data[2][j][k];
            F2_out[num_r*k+j] = data[3][j][k];
            F0_out[num_r*k+j] = data[4][j][k];
            H1_out[num_r*k+j] = data[5][j][k];
            H2_out[num_r*k+j] = data[6][j][k];
            Vproca_out[num_r*k+j] = data[7][j][k];
            
            
        }
    }
    
    // the spacing in each coordinate is
    CCTK_REAL dX     = (sh1[num_r-1] - sh1[0])/(num_r-1);
    CCTK_REAL dtheta = (th1[num_th-1] - th1[0])/(num_th-1);
    
    // make sure spacing is uniform in the provided grid
    for (int i = 1; i < num_r; i++) {
        if (fabs(sh1[i] - sh1[i - 1] - dX) > epsilon) {
            printf("i = %d\n", i);
            CCTK_WARN(0, "X grid is not uniformly spaced. Aborting.");
        }
    }
    for (int j = 1; j < num_th; j++) {
        if (fabs(th1[j] - th1[j - 1] - dtheta) > epsilon)
        CCTK_WARN(0, "theta grid is not uniformly spaced. Aborting.");
    }
    
    CCTK_INT N_interp_points = cctk_lsh[0]*cctk_lsh[1]*cctk_lsh[2]; // total points
    
    CCTK_REAL *X_g, *theta_g;
    X_g     = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
    theta_g = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
    
    for (int k = 0; k < cctk_lsh[2]; ++k)
    for (int j = 0; j < cctk_lsh[1]; ++j)
    for (int i = 0; i < cctk_lsh[0]; ++i) {
        
        const CCTK_INT i3D  = CCTK_GFINDEX3D (cctkGH, i, j, k);
        
        const CCTK_REAL xx  = x[i3D];
        const CCTK_REAL yy  = y[i3D];
        const CCTK_REAL zz  = z[i3D];
        
        const CCTK_REAL R_sq = xx*xx + yy*yy + zz*zz;
        const CCTK_REAL R  = sqrt(R_sq);
        
        // X compactified radial coordinate (used in input files) as function of isotropic
        // radial coordinate R
        CCTK_REAL lX = R / (1.0 + R);
        
        CCTK_REAL ltheta = acos( zz/(R + epsilon));
        
       /* if((z[i3D]<0.0)){
            /* sintheta = sin(3.14159265358979323846 + atan(rho/zz));
            ltheta = 3.14159265358979323846 + atan((SQUARE(xx)+SQUARE(yy))/zz);
        }*/
        
        X_g[i3D]     = lX;
        theta_g[i3D] = ltheta;
    }
    
    /* now for the interpolation */
    
    const CCTK_INT N_dims  = 2;   // 2-D interpolation
    
    const CCTK_INT N_input_arrays  = 6;
    const CCTK_INT N_output_arrays = 6;
    
    /* origin and stride of the input coordinates. with this Cactus reconstructs
     the whole X and theta array. */
    CCTK_REAL origin[N_dims];
    CCTK_REAL delta [N_dims];
    origin[0] = sh1[0];  origin[1] = th1[0];
    delta[0]  = dX;    delta[1]  = dtheta;
    
    /* points onto which we want to interpolate, ie, the grid points themselves in
     (X, theta) coordinates (computed above) */
    const void *interp_coords[N_dims];
    interp_coords[0] = (const void *) X_g;
    interp_coords[1] = (const void *) theta_g;
    
    
    /* input arrays */
    const void *input_arrays[N_input_arrays];
    CCTK_INT input_array_type_codes[N_input_arrays];
    CCTK_INT input_array_dims[N_dims];
    input_array_dims[0] = num_r;
    input_array_dims[1] = num_th;
    
    input_array_type_codes[0] = CCTK_VARIABLE_REAL;
    input_array_type_codes[1] = CCTK_VARIABLE_REAL;
    input_array_type_codes[2] = CCTK_VARIABLE_REAL;
    input_array_type_codes[3] = CCTK_VARIABLE_REAL;
    input_array_type_codes[4] = CCTK_VARIABLE_REAL;
    input_array_type_codes[5] = CCTK_VARIABLE_REAL;
    
    /* Cactus stores and expects arrays in Fortran order, that is, faster in the
     first index. this is compatible with our input file, where the X coordinate
     is faster. */
    input_arrays[0] = (const void *) F1_out;
    input_arrays[1] = (const void *) F2_out;
    input_arrays[2] = (const void *) F0_out;
    input_arrays[3] = (const void *) H1_out;
    input_arrays[4] = (const void *) H2_out;
    input_arrays[5] = (const void *) Vproca_out;
    
    
    /* output arrays */
    void *output_arrays[N_output_arrays];
    CCTK_INT output_array_type_codes[N_output_arrays];
    
    CCTK_REAL *Vproca_aux, *F2_aux, *F0_aux, *F1_aux, *H1_aux, *H2_aux;
    
    F2_aux    = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
    F0_aux    = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
    F1_aux    = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
    H1_aux    = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
    H2_aux    = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
    Vproca_aux    = (CCTK_REAL *) malloc(N_interp_points * sizeof(CCTK_REAL));
    
    output_array_type_codes[0] = CCTK_VARIABLE_REAL;
    output_array_type_codes[1] = CCTK_VARIABLE_REAL;
    output_array_type_codes[2] = CCTK_VARIABLE_REAL;
    output_array_type_codes[3] = CCTK_VARIABLE_REAL;
    output_array_type_codes[4] = CCTK_VARIABLE_REAL;
    output_array_type_codes[5] = CCTK_VARIABLE_REAL;
    
    output_arrays[0] = (void *) F1_aux;
    output_arrays[1] = (void *) F2_aux;
    output_arrays[2] = (void *) F0_aux;
    output_arrays[3] = (void *) H1_aux;
    output_arrays[4] = (void *) H2_aux;
    output_arrays[5] = (void *) Vproca_aux;
    
    
    /* handle and settings for the interpolation routine */
    int operator_handle, param_table_handle;
    operator_handle    = CCTK_InterpHandle("Lagrange polynomial interpolation");
    param_table_handle = Util_TableCreateFromString("order=4 boundary_extrapolation_tolerance={0.1 1.0 0.05 0.05}");
    CCTK_INFO("Interpolating result...");
    
    /* do the actual interpolation, and check for error returns */
    int status = CCTK_InterpLocalUniform(N_dims, operator_handle,
                                         param_table_handle,
                                         origin, delta,
                                         N_interp_points,
                                         CCTK_VARIABLE_REAL,
                                         interp_coords,
                                         N_input_arrays, input_array_dims,
                                         input_array_type_codes,
                                         input_arrays,
                                         N_output_arrays, output_array_type_codes,
                                         output_arrays);
    if (status < 0) {
        CCTK_VWarn(0, __LINE__, __FILE__, CCTK_THORNSTRING,
                   "interpolation screwed up!");
    }
    
    free(X_g); free(theta_g);
    
    free(F1_out); free(F2_out); free(F0_out); free(H1_out); free(H2_out); free(Vproca_out);
    
    
    
    CCTK_VInfo(CCTK_THORNSTRING,"Setting up Proca ID");
    
    
    for (int k = 0; k < cctk_lsh[2]; ++k)
    for (int j = 0; j < cctk_lsh[1]; ++j)
    for (int i = 0; i < cctk_lsh[0]; ++i) {
            
			CCTK_INT i3D = CCTK_GFINDEX3D(cctkGH, i, j, k);
			CCTK_REAL xx, yy, zz;
			
			CCTK_REAL R, costheta, sintheta;//the radial coordinate introduced by Liu et al 2010
            CCTK_REAL cosphi, sinphi;
			CCTK_REAL rho, rho_sq;
			CCTK_REAL r_plus, r_minus, r_BL, Sigma, A, Delta;
			CCTK_REAL gii, gjj, gkk;//spatial metric, i=r, j=theta, k=phi
			CCTK_REAL Kik, Kjk;//extrinsic curvature
			CCTK_REAL betak;//shift
			CCTK_REAL a_sq, costh_sq, sinth_sq;
			CCTK_REAL r_BL_sq, R_sq;
			CCTK_REAL detg;
                        CCTK_REAL det_metric, sdetg;
                        CCTK_REAL guxx, guxy, guxz, guyy, guyz, guzz;
            
            CCTK_REAL theta, phi;
            
            /*CCTK_REAL F1_aux, F2_aux, F0_aux, Wproca_aux, H1_aux, H2_aux, H3_aux, Vproca_aux;*/
            CCTK_REAL sv1, sv2, av1li, av1lj, av1lk, av2li, av2lj, av2lk; //auxProca Varibles.
            CCTK_REAL av1ui, av1uj, av1uk, av2ui, av2uj, av2uk; //auxProca Varibles.
            CCTK_REAL Efield1li, Efield1lj, Efield1lk, Efield2li, Efield2lj, Efield2lk;
            CCTK_REAL Efield1lx, Efield1ly, Efield1lz, Efield2lx, Efield2ly, Efield2lz;
            CCTK_REAL alp_aux;

			CCTK_REAL jac11, jac12, jac13, jac21, jac22, jac23, jac31, jac32;//elements
			//of the coordinate transformation matrix from spherical polar to Cartesian coordinates
			//jac33 is zero for this transformation!

			CCTK_REAL g11, g12, g13, g22, g23, g33;////////////////////////////////////////////////
			CCTK_REAL k11, k12, k13, k22, k23, k33;//variables for the rotation of the spacetime!//
			CCTK_REAL beta1, beta2, beta3;/////////////////////////////////////////////////////////
	
        
				xx = x[i3D];
				yy = y[i3D];
				zz = z[i3D];
            
                    /*   r_plus = mass + sqrt(SQUARE(mass) - a_sq);
            r_minus = mass - sqrt(SQUARE(mass) - a_sq);
            
            r_BL = R * SQUARE(1.0 + r_plus/(4.0 * R));//the choice in Liu et al 2010, eq 11
            r_BL_sq = SQUARE(r_BL);
            
            Sigma = r_BL_sq + a_sq * costh_sq;
            Delta = r_BL_sq - 2.0 * mass * r_BL + a_sq;
            A = SQUARE(r_BL_sq + a_sq) - Delta * a_sq * sinth_sq; */
            
            rho_sq = SQUARE(xx)+SQUARE(yy);
            
            rho = sqrt(rho_sq);
            
            
            
            if(rho<epsilon)
                
            {
                xx = epsilon;
                rho_sq = SQUARE(xx)+SQUARE(yy);
                rho = sqrt(rho_sq);
            }
            
            
            R_sq = rho_sq+SQUARE(zz);
            
            R = sqrt(R_sq);
            
            if(R == 0.0)
                
            {
                CCTK_VInfo(CCTK_THORNSTRING,"R is still zero, something went wrong!");
            }
            
            
            
            /*costheta = zz / R;
            sintheta = rho / R;*/
            
            theta = atan(rho/zz);
            /*theta = abs(acos(costheta));*/

            if((z[i3D]<0.0)){
                 /* sintheta = sin(3.14159265358979323846 + atan(rho/zz));*/
		   theta = 3.14159265358979323846 + atan(rho/zz);
            }

            costheta = zz / R;
            sintheta = rho / R;

           /* printf("%f%f%f%f%f\n",xx,yy,zz,theta,abs(atan(rho/zz)));*/
            
            
            phi = atan(yy/xx);
            
            if((xx<0.0)){
                phi = atan(yy/xx) + 3.14159265358979323846;
            }
            
            cosphi = cos(mode*phi);
            sinphi = sin(mode*phi);
            
            /*  a_sq = SQUARE(a);*/
            costh_sq = SQUARE(costheta);
            sinth_sq = SQUARE(sintheta);
            
            
            if(CCTK_Equals(initial_lapse,"proca-ID"))	{
                alp[i3D] = exp(F0_aux[i3D]);
                alp_aux = alp[i3D];
            }
            
            if(CCTK_Equals(initial_lapse,"proca-one"))	{
                alp[i3D] = 1.0;	
            }
            
            betak = 0.0;
            
            gii = exp(2.0*F1_aux[i3D]);
            gjj = exp(2.0*F1_aux[i3D])*R_sq;
            gkk = exp(2.0*F2_aux[i3D])*R_sq*sinth_sq;

            if(CCTK_Equals(perturbation,"cv"))    {
            betak = 0.0;
            gii = 1.0;
            gjj = R_sq;
            gkk = R_sq*sinth_sq;
            alp_aux = 1.0;
            alp[i3D] = alp_aux;
            }
  
            sv1 = 0.0; /*+ factor*(Vproca_aux)*sinphi/alp_aux;*/
            sv2 = - (Vproca_aux[i3D])/alp_aux; /**cosphi/alp_aux;*/
            
            av1li = H1_aux[i3D]/R; /* *cosphi;*/
            av2li = 0.0; /*factor*H1_aux/R*sinphi;*/
            
            av1lj = H2_aux[i3D]; /* *cosphi;*/
            av2lj = 0.0; /*factor*H2_aux*sinphi;*/
            
            av1lk = 0.0;
            av2lk = 0.0;
            
         /*   av1li = gii*av1ui;
            av1lj = gjj*av1uj;
            av1lk = gkk*av1uk;
            
            av2li = gii*av2ui;
            av2lj = gjj*av2uj;
            av2lk = gkk*av2uk;*/
            
            Efield1li = 0.0; /*- factor*H1_aux/R*sinphi;*/
            Efield1lj = 0.0; /*- factor*H2_aux*sinphi;*/
            Efield1lk = 0.0;
            
            Efield2li = + H1_aux[i3D]/R; /**cosphi;*/
            Efield2lj = + H2_aux[i3D]; /**cosphi;*/
            Efield2lk = 0.0;

            
          /*  El1ux[i3D] = - Efield1ui*sinphi;
            El1uy[i3D] = - Efield1uj*sinphi;
            El1uz[i3D] = - Efield1uk*cosphi;
            
            El2ux[i3D] = + Efield2ui*cosphi;
            El2uy[i3D] = + Efield2uj*cosphi;
            El2uz[i3D] = - Efield2uk*sinphi;*/
            
            

 			

		/*	Kik = mass*a*sinth_sq/(Sigma*sqrt(A*Sigma)) *
						(3.0*SQUARE(r_BL_sq) + 2.0*a_sq*r_BL_sq -
						SQUARE(a_sq) -a_sq*(r_BL_sq - a_sq)*sinth_sq) *
						(1.0+r_plus/(4.0*R)) * 
						(1.0/sqrt(R*(r_BL - r_minus)));
			Kjk = -2.0*a*a_sq*mass*r_BL*costheta*sintheta*sinth_sq / 
						(Sigma*sqrt(A*Sigma)) * (R-(r_plus/4.0)) * 
						(sqrt((r_BL-r_minus)/R));*/

		/*	betak = - 2.0 * mass * a * r_BL / A; */

			//coordinate transformation to cactus cartesian grid
			jac11 = xx/R;//////////////////////////////////////////////////////////////
			jac12 = xx*zz/(R_sq*rho);//////////////////////////////////////////////////
			jac13 = -yy/rho_sq;////////////////////////////////////////////////////////
			jac21 = yy/R;//////////////////////////////////////////////////////////////
			jac22 = yy*zz/(R_sq*rho);/////////////Jacobian from Spherical//////////////  
			jac23 = xx/rho_sq;/////////////Polar Coordinates to Cartesian Coordinates//
			jac31 = zz/R;//for COVARIANT tensors!!!////////////////////////////////////
			jac32 = -rho/R_sq;/////////////////////////////////////////////////////////
			//jac33 = 0.0;/////////////////////////////////////////////////////////////

			gxx[i3D] = gii*SQUARE(jac11) + gjj*SQUARE(jac12) + gkk*SQUARE(jac13);//Jac*g_ij*Tranpose[Jac]
			gxy[i3D] = gii*jac11*jac21 + gjj*jac12*jac22 + gkk*jac13*jac23;
			gxz[i3D] = gii*jac11*jac31 + gjj*jac12*jac32;
			gyy[i3D] = gii*SQUARE(jac21) + gjj*SQUARE(jac22) + gkk*SQUARE(jac23);
			gyz[i3D] = gii*jac21*jac31 + gjj*jac22*jac32;
			gzz[i3D] = gii*SQUARE(jac31) + gjj*SQUARE(jac32);
                          
            det_metric = gxx[i3D]*gyy[i3D]*gzz[i3D]-gxx[i3D]*gyz[i3D]*gyz[i3D]-gxy[i3D]*gxy[i3D]*gzz[i3D]+gxy[i3D]*gyz[i3D]*gxz[i3D]+gxz[i3D]*gxy[i3D]*gyz[i3D]-gxz[i3D]*gyy[i3D]*gxz[i3D];
            
            guxx = (-gyz[i3D]*gyz[i3D] + gyy[i3D]*gzz[i3D])/det_metric;
            guxy = ( gxz[i3D]*gyz[i3D] - gxy[i3D]*gzz[i3D])/det_metric;
            guxz = (-gxz[i3D]*gyy[i3D] + gxy[i3D]*gyz[i3D])/det_metric;
            guyy = (-gxz[i3D]*gxz[i3D] + gxx[i3D]*gzz[i3D])/det_metric;
            guyz = ( gxy[i3D]*gxz[i3D] - gxx[i3D]*gyz[i3D])/det_metric;
            guzz = (-gxy[i3D]*gxy[i3D] + gxx[i3D]*gyy[i3D])/det_metric;

			kxx[i3D] = 0.0;//Jac*K_ij*Tranpose[Jac]
			kxy[i3D] = 0.0;
			kxz[i3D] = 0.0;
			kyy[i3D] = 0.0;
            kyz[i3D] = 0.0;
			kzz[i3D] = 0.0;
            
            /*av1lx[i3D] = av1li;
            av1ly[i3D] = av1lj;
            av1lz[i3D] = av1lk;
            
            av2lx[i3D] = av2li;
            av2ly[i3D] = av2lj;
            av2lz[i3D] = av2lk;*/

                if(CCTK_Equals(perturbation,"new"))    {
                sv1 = + (Vproca_aux[i3D])/alp_aux*(sinphi + fac*sin(2.0*phi));
                sv2 = - (Vproca_aux[i3D])/alp_aux*(cosphi + fac*cos(2.0*phi));
                av1li   = + H1_aux[i3D]/R*(cosphi + fac*cos(modepert*phi));
                av2li   = + H1_aux[i3D]/R*(sinphi + fac*sin(modepert*phi));
                av1lj   = + H2_aux[i3D]*(cosphi + fac*cos(modepert*phi));
                av2lj   = + H2_aux[i3D]*(sinphi + fac*sin(modepert*phi));
                av1lk   = 0.0;
                av2lk   = 0.0;
                Efield1li   = - H1_aux[i3D]/R*(sinphi + fac*sin(modepert*phi));
                Efield1lj   = - H2_aux[i3D]*(sinphi + fac*sin(modepert*phi));
                Efield1lk   = 0.0;
                Efield2li   = + H1_aux[i3D]/R*(cosphi + fac*cos(modepert*phi));
                Efield2lj   = + H2_aux[i3D]*(cosphi + fac*cos(modepert*phi));
                Efield2lk   = 0.0;
            }
            
            A1x[i3D] = av1li*jac11 + av1lj*jac12 + av1lk*jac13;
            A1y[i3D] = av1li*jac21 + av1lj*jac22 + av1lk*jac23;
            A1z[i3D] = av1li*jac31 + av1lj*jac32;
            
            A2x[i3D] = av2li*jac11 + av2lj*jac12 + av2lk*jac13;
            A2y[i3D] = av2li*jac21 + av2lj*jac22 + av2lk*jac23;
            A2z[i3D] = av2li*jac31 + av2lj*jac32;
            
            Efield1lx = Efield1li*jac11 + Efield1lj*jac12 + Efield1lk*jac13;
            Efield1ly = Efield1li*jac21 + Efield1lj*jac22 + Efield1lk*jac23;
            Efield1lz = Efield1li*jac31 + Efield1lj*jac32;
            
            Efield2lx = Efield2li*jac11 + Efield2lj*jac12 + Efield2lk*jac13;
            Efield2ly = Efield2li*jac21 + Efield2lj*jac22 + Efield2lk*jac23;
            Efield2lz = Efield2li*jac31 + Efield2lj*jac32;

            E1x[i3D] = Efield1lx*guxx + Efield1ly*guxy + Efield1lz*guxz;
            E1y[i3D] = Efield1lx*guxy + Efield1ly*guyy + Efield1lz*guyz;
            E1z[i3D] = Efield1lx*guxz + Efield1ly*guyz + Efield1lz*guzz;

            E2x[i3D] = Efield2lx*guxx + Efield2ly*guxy + Efield2lz*guxz;
            E2y[i3D] = Efield2lx*guxy + Efield2ly*guyy + Efield2lz*guyz;
            E2z[i3D] = Efield2lx*guxz + Efield2ly*guyz + Efield2lz*guzz;
            
            Aphi1[i3D] = sv1;
            Aphi2[i3D] = sv2;
            
            Zeta1[i3D] = 0.0;
            Zeta2[i3D] = 0.0;
            
            
            
			if(CCTK_Equals(initial_shift,"proca-ID"))	{
				betax[i3D] = -yy*betak;
				betay[i3D] = xx*betak;
/*                betax[i3D] = betak*gkk*jac13*guxx + betak*gkk*jac23*guxy;
                betay[i3D] = betak*gkk*jac13*guxy + betak*gkk*jac23*guyy;*/
				betaz[i3D] = 0.0;
			}
			
			if(CCTK_Equals(initial_shift,"zero"))	{
				betax[i3D] = 0.0;
				betay[i3D] = 0.0;
				betaz[i3D] = 0.0;
			}
            
            if(CCTK_Equals(perturbation,"m2"))    {
                Aphi1[i3D] = Aphi1[i3D]*(1.0 + fac*sinth_sq*cos(2.0*phi));
                Aphi2[i3D] = Aphi2[i3D]*(1.0 + fac*sinth_sq*sin(2.0*phi));
                E1x[i3D]   =   E1x[i3D]*(1.0 + fac*sinth_sq*cos(2.0*phi));
                E1y[i3D]   =   E1y[i3D]*(1.0 + fac*sinth_sq*cos(2.0*phi));
                E1z[i3D]   =   E1z[i3D]*(1.0 + fac*sinth_sq*cos(2.0*phi));
                E2x[i3D]   =   E2x[i3D]*(1.0 + fac*sinth_sq*sin(2.0*phi));
                E2y[i3D]   =   E2y[i3D]*(1.0 + fac*sinth_sq*sin(2.0*phi));
                E2z[i3D]   =   E2z[i3D]*(1.0 + fac*sinth_sq*sin(2.0*phi));
                A1x[i3D]   =   A1x[i3D]*(1.0 + fac*sinth_sq*cos(2.0*phi));
                A1y[i3D]   =   A1y[i3D]*(1.0 + fac*sinth_sq*cos(2.0*phi));
                A1z[i3D]   =   A1z[i3D]*(1.0 + fac*sinth_sq*cos(2.0*phi));
                A2x[i3D]   =   A2x[i3D]*(1.0 + fac*sinth_sq*sin(2.0*phi));
                A2y[i3D]   =   A2y[i3D]*(1.0 + fac*sinth_sq*sin(2.0*phi));
                A2z[i3D]   =   A2z[i3D]*(1.0 + fac*sinth_sq*sin(2.0*phi));
            }

            if(CCTK_Equals(perturbation,"rad"))    {
                Aphi1[i3D] = Aphi1[i3D]*(1.0 + fac);
                Aphi2[i3D] = Aphi2[i3D]*(1.0 + fac);
                E1x[i3D]   =   E1x[i3D]*(1.0 + fac);
                E1y[i3D]   =   E1y[i3D]*(1.0 + fac);
                E1z[i3D]   =   E1z[i3D]*(1.0 + fac);
                E2x[i3D]   =   E2x[i3D]*(1.0 + fac);
                E2y[i3D]   =   E2y[i3D]*(1.0 + fac);
                E2z[i3D]   =   E2z[i3D]*(1.0 + fac);
                A1x[i3D]   =   A1x[i3D]*(1.0 + fac);
                A1y[i3D]   =   A1y[i3D]*(1.0 + fac);
                A1z[i3D]   =   A1z[i3D]*(1.0 + fac);
                A2x[i3D]   =   A2x[i3D]*(1.0 + fac);
                A2y[i3D]   =   A2y[i3D]*(1.0 + fac);
                A2z[i3D]   =   A2z[i3D]*(1.0 + fac);
            }
            
            if(CCTK_Equals(perturbation,"cos"))    {
                Aphi1[i3D] = Aphi1[i3D]*(1.0 + fac*cos(2.0*phi));
                Aphi2[i3D] = Aphi2[i3D]*(1.0 + fac*cos(2.0*phi));
                E1x[i3D]   =   E1x[i3D]*(1.0 + fac*cos(2.0*phi));
                E1y[i3D]   =   E1y[i3D]*(1.0 + fac*cos(2.0*phi));
                E1z[i3D]   =   E1z[i3D]*(1.0 + fac*cos(2.0*phi));
                E2x[i3D]   =   E2x[i3D]*(1.0 + fac*cos(2.0*phi));
                E2y[i3D]   =   E2y[i3D]*(1.0 + fac*cos(2.0*phi));
                E2z[i3D]   =   E2z[i3D]*(1.0 + fac*cos(2.0*phi));
                A1x[i3D]   =   A1x[i3D]*(1.0 + fac*cos(2.0*phi));
                A1y[i3D]   =   A1y[i3D]*(1.0 + fac*cos(2.0*phi));
                A1z[i3D]   =   A1z[i3D]*(1.0 + fac*cos(2.0*phi));
                A2x[i3D]   =   A2x[i3D]*(1.0 + fac*cos(2.0*phi));
                A2y[i3D]   =   A2y[i3D]*(1.0 + fac*cos(2.0*phi));
                A2z[i3D]   =   A2z[i3D]*(1.0 + fac*cos(2.0*phi));
            }

            if(CCTK_Equals(perturbation,"sin"))    {
                Aphi1[i3D] = Aphi1[i3D]*(1.0 + fac*sin(2.0*phi));
                Aphi2[i3D] = Aphi2[i3D]*(1.0 + fac*sin(2.0*phi));
                E1x[i3D]   =   E1x[i3D]*(1.0 + fac*sin(2.0*phi));
                E1y[i3D]   =   E1y[i3D]*(1.0 + fac*sin(2.0*phi));
                E1z[i3D]   =   E1z[i3D]*(1.0 + fac*sin(2.0*phi));
                E2x[i3D]   =   E2x[i3D]*(1.0 + fac*sin(2.0*phi));
                E2y[i3D]   =   E2y[i3D]*(1.0 + fac*sin(2.0*phi));
                E2z[i3D]   =   E2z[i3D]*(1.0 + fac*sin(2.0*phi));
                A1x[i3D]   =   A1x[i3D]*(1.0 + fac*sin(2.0*phi));
                A1y[i3D]   =   A1y[i3D]*(1.0 + fac*sin(2.0*phi));
                A1z[i3D]   =   A1z[i3D]*(1.0 + fac*sin(2.0*phi));
                A2x[i3D]   =   A2x[i3D]*(1.0 + fac*sin(2.0*phi));
                A2y[i3D]   =   A2y[i3D]*(1.0 + fac*sin(2.0*phi));
                A2z[i3D]   =   A2z[i3D]*(1.0 + fac*sin(2.0*phi));
            }

            if(CCTK_Equals(perturbation,"cos1"))    {
                Aphi1[i3D] = Aphi1[i3D]*(1.0 + fac*cos(1.0*phi));
                Aphi2[i3D] = Aphi2[i3D]*(1.0 + fac*cos(1.0*phi));
                E1x[i3D]   =   E1x[i3D]*(1.0 + fac*cos(1.0*phi));
                E1y[i3D]   =   E1y[i3D]*(1.0 + fac*cos(1.0*phi));
                E1z[i3D]   =   E1z[i3D]*(1.0 + fac*cos(1.0*phi));
                E2x[i3D]   =   E2x[i3D]*(1.0 + fac*cos(1.0*phi));
                E2y[i3D]   =   E2y[i3D]*(1.0 + fac*cos(1.0*phi));
                E2z[i3D]   =   E2z[i3D]*(1.0 + fac*cos(1.0*phi));
                A1x[i3D]   =   A1x[i3D]*(1.0 + fac*cos(1.0*phi));
                A1y[i3D]   =   A1y[i3D]*(1.0 + fac*cos(1.0*phi));
                A1z[i3D]   =   A1z[i3D]*(1.0 + fac*cos(1.0*phi));
                A2x[i3D]   =   A2x[i3D]*(1.0 + fac*cos(1.0*phi));
                A2y[i3D]   =   A2y[i3D]*(1.0 + fac*cos(1.0*phi));
                A2z[i3D]   =   A2z[i3D]*(1.0 + fac*cos(1.0*phi));
            }
        
            if(CCTK_Equals(perturbation,"cos3"))    {
                Aphi1[i3D] = Aphi1[i3D]*(1.0 + fac*cos(3.0*phi));
                Aphi2[i3D] = Aphi2[i3D]*(1.0 + fac*cos(3.0*phi));
                E1x[i3D]   =   E1x[i3D]*(1.0 + fac*cos(3.0*phi));
                E1y[i3D]   =   E1y[i3D]*(1.0 + fac*cos(3.0*phi));
                E1z[i3D]   =   E1z[i3D]*(1.0 + fac*cos(3.0*phi));
                E2x[i3D]   =   E2x[i3D]*(1.0 + fac*cos(3.0*phi));
                E2y[i3D]   =   E2y[i3D]*(1.0 + fac*cos(3.0*phi));
                E2z[i3D]   =   E2z[i3D]*(1.0 + fac*cos(3.0*phi));
                A1x[i3D]   =   A1x[i3D]*(1.0 + fac*cos(3.0*phi));
                A1y[i3D]   =   A1y[i3D]*(1.0 + fac*cos(3.0*phi));
                A1z[i3D]   =   A1z[i3D]*(1.0 + fac*cos(3.0*phi));
                A2x[i3D]   =   A2x[i3D]*(1.0 + fac*cos(3.0*phi));
                A2y[i3D]   =   A2y[i3D]*(1.0 + fac*cos(3.0*phi));
                A2z[i3D]   =   A2z[i3D]*(1.0 + fac*cos(3.0*phi));
            }

     
		}

    CCTK_VInfo(CCTK_THORNSTRING,"Done IDProca6");

		
		
	/*if(fill_past_timelevels){
	CCTK_VInfo(CCTK_THORNSTRING,"Setting past timelevels by copying!");
	
			LC_LOOP3(fill_timelevels, i, j, k, 0, 0, 0, ni, nj, nk, ni, nj, nk)
			{
				CCTK_INT i3D = CCTK_GFINDEX3D(cctkGH, i, j, k);
				
				gxx_p[i3D]   = gxx[i3D];
				gxx_p_p[i3D] = gxx[i3D];
				
				gxy_p[i3D]   = gxy[i3D];
				gxy_p_p[i3D] = gxy[i3D];
				
				gxz_p[i3D]   = gxz[i3D];
				gxz_p_p[i3D] = gxz[i3D];
				
				gyy_p[i3D]   = gyy[i3D];
				gyy_p_p[i3D] = gyy[i3D];
				
				gyz_p[i3D]   = gyz[i3D];
				gyz_p_p[i3D] = gyz[i3D];
				
				gzz_p[i3D]   = gzz[i3D];
				gzz_p_p[i3D] = gzz[i3D];
				
				
				/*kxx_p[i3D]   = kxx[i3D];
				kxx_p_p[i3D] = kxx[i3D];
				
				kxy_p[i3D]   = kxy[i3D];
				kxy_p_p[i3D] = kxy[i3D];
				
				kxz_p[i3D]   = kxz[i3D];
				kxz_p_p[i3D] = kxz[i3D];
				
				kyy_p[i3D]   = kyy[i3D];
				kyy_p_p[i3D] = kyy[i3D];
				
				kyz_p[i3D]   = kyz[i3D];
				kyz_p_p[i3D] = kyz[i3D];
				
				kzz_p[i3D]   = kzz[i3D];
				kzz_p_p[i3D] = kzz[i3D];*//*
				
				
				betax_p[i3D]   = betax[i3D];
				betax_p_p[i3D] = betax[i3D];
				
				betay_p[i3D]   = betay[i3D];
				betay_p_p[i3D] = betay[i3D];
				
				betaz_p[i3D]   = betaz[i3D];
				betaz_p_p[i3D] = betaz[i3D];
				
				alp_p[i3D]   = alp[i3D]; 
				alp_p_p[i3D] = alp[i3D];
                
                Phis1_p[i3D] = Phis1[i3D];
                Phis1_p_p[i3D] = Phis1[i3D];
                
                Phis2_p[i3D] = Phis2[i3D];
                Phis2_p_p[i3D] = Phis2[i3D];
                
                av1lx_p[i3D] = av1lx[i3D];
                av1ly_p[i3D] = av1ly[i3D];
                av1lz_p[i3D] = av1lz[i3D];
            
                av2lx_p[i3D] = av2lx[i3D];
                av2ly_p[i3D] = av2ly[i3D];
                av2lz_p[i3D] = av2lz[i3D];
                
                av1lx_p_p[i3D] = av1lx[i3D];
                av1ly_p_p[i3D] = av1ly[i3D];
                av1lz_p_p[i3D] = av1lz[i3D];
                
                av2lx_p_p[i3D] = av2lx[i3D];
                av2ly_p_p[i3D] = av2ly[i3D];
                av2lz_p_p[i3D] = av2lz[i3D];
                
                
				
			}LC_ENDLOOP3(fill_timelevels);
		

	}*/
	return;
}
