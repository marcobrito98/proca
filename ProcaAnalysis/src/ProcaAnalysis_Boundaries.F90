#include "cctk.h"
#include "cctk_Arguments.h"
#include "cctk_Parameters.h"
#include "cctk_Functions.h"

subroutine ProcaAnalysis_Boundaries( CCTK_ARGUMENTS )

  implicit none
  DECLARE_CCTK_ARGUMENTS
  DECLARE_CCTK_PARAMETERS
  DECLARE_CCTK_FUNCTIONS

  CCTK_INT ierr

  CCTK_INT, parameter :: one = 1

  ! The outgoing (radiative) boundary conditions are being handled from calls to
  ! the NewRad infraestructure. Here we register all BCs as 'none', which
  ! enforces all the symmetry BCs.

  ierr = Boundary_SelectGroupForBC(cctkGH, CCTK_ALL_FACES, one, -one,        &
       "ProcaAnalysis::Proca_mass_point", "flat")
  if (ierr < 0)                                                            &
       call CCTK_WARN(0, "Failed to register BC for ProcaAnalysis::Proca_energy_point!")

end subroutine ProcaAnalysis_Boundaries
