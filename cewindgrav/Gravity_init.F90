!!****if* source/physics/Gravity/GravityMain/CeWindGrav/Gravity_init
!!
!! NAME
!!
!!  Gravity_init
!!  
!! SYNOPSIS
!!
!!  Gravity_init()
!!
!! DESCRIPTION
!!
!!  This routine initializes the gravitational physics unit for Pointmass.
!!
!! ARGUMENTS
!!
!!  
!!
!!***

subroutine Gravity_init()

  use Gravity_data
  use Driver_interface, ONLY : Driver_getMype, Driver_getNumProcs

  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use PhysicalConstants_interface, ONLY:  PhysicalConstants_get

  implicit none

#include "constants.h"

  ! Everybody should know these
  call Driver_getMype(MESH_COMM,grv_meshMe)
  call Driver_getNumProcs(MESH_COMM,grv_meshMe)

  call PhysicalConstants_get("newton", grv_newton)

  call RuntimeParameters_get("ptxpos", grv_ptxpos)
  call RuntimeParameters_get("ptypos", grv_ptypos)
  call RuntimeParameters_get("ptzpos", grv_ptzpos)
  call RuntimeParameters_get("ptmass", grv_ptmass)
  call RuntimeParameters_get("useGravity", useGravity)

  grv_factor = -grv_newton * grv_ptmass


  
  !print*, "CEWINDGRAV has been initialized"
  
  
!==============================================================================

!==============================================================================

  return
end subroutine Gravity_init
