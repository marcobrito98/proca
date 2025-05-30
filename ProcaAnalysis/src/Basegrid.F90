! Basegrid.F90 : Register symmetries of the grid functions
!
!=============================================================================

#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"


subroutine procaanalysis_symmetries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS

  CCTK_INT ierr

 !  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaAnalysis::Proca_mass_point" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaAnalysis::Proca_energy_point" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaAnalysis::Proca_angmomentum_point" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaAnalysis::Proca_rho_point" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/), "ProcaAnalysis::Proca_angmomentum_point_y" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/),"ProcaAnalysis::Proca_energy_point_r150" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/),"ProcaAnalysis::Proca_angmomentum_point_r150" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/),"ProcaAnalysis::Proca_rho_point_r150" )

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,1/),"ProcaAnalysis::Proca_energy_point_r100" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,1/),"ProcaAnalysis::Proca_angmomentum_point_r100" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/),"ProcaAnalysis::Proca_rho_point_r100")

  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,1/),"ProcaAnalysis::Proca_energy_point_r60" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,1/),"ProcaAnalysis::Proca_angmomentum_point_r60" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/),"ProcaAnalysis::Proca_rho_point_r60")


  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,1/),"ProcaAnalysis::Proca_energy_point_r30" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1,1/),"ProcaAnalysis::Proca_angmomentum_point_r30" )
  call SetCartSymVN( ierr, cctkGH, (/ 1, 1, 1/),"ProcaAnalysis::Proca_rho_point_r30")

end subroutine procaanalysis_symmetries
