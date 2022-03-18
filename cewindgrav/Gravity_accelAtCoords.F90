!!****if* source/physics/Gravity/GravityMain/CeWindGrav/Gravity_accelAtCoords
!!
!! NAME
!!
!!  Gravity_accelAtCoords 
!!
!! SYNOPSIS
!!
!!  Gravity_accelAtCoords(integer(IN) :: numPoints,
!!                      real(IN)      :: iCoords(:),
!!                      real(IN)      :: jCoords(:),
!!                      real(IN)      :: kCoords(:),
!!                      integer(IN)   :: accelDir,
!!                      real(OUT)     :: accel(numPoints),
!!                      integer(IN)   :: blockID,
!!                      integer(IN),optional :: potentialIndex)
!!
!! DESCRIPTION
!!
!!  This routine computes the gravitational acceleration in a
!!  specified direction for a vector of points given by their
!!  coordinates.
!!
!! ARGUMENTS
!!
!!  iCoords,jCoords,kCoords: coordinates of the points where the
!!                           gravitational accelation is requested.
!!                           Each of these arrays should either be
!!                           of lenght numPoints (or more), in which
!!                           case its nth value is used for the nth
!!                           point; or else of dimension 1, in which
!!                           case the value is used for all points.
!!  accelDir :    The acceleration direction:  allowed values are 
!!              IAXIS, JAXIS and IAXIS. These values are defined
!!              in constants.h.
!!  numPoints :  Number of cells to update in accel()
!!  accel     :   Array to receive results
!!  blockID  :  The local identifier of the block to work on,
!!                not applicable in pointmass gravity.
!!  potentialIndex :  optional, not applicable in pointmass gravity
!! 
!!***

subroutine Gravity_accelAtCoords (numPoints, iCoords,jCoords,kCoords, accelDir,&
     accel, blockID, &
     potentialIndex)

!=======================================================================

  use Gravity_data, ONLY: grv_ptxpos, grv_ptypos, grv_ptzpos, grv_factor, &
       useGravity, grv_newton

  use Simulation_data, ONLY : sim_mprimary, sim_yPrim , sim_radius

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: accelDir, numPoints
  real, dimension(:),INTENT(in) :: iCoords,jCoords,kCoords
  real, dimension(numPoints),INTENT(OUT) :: accel
  integer, intent(IN),optional :: blockID
  integer, intent(IN),optional :: potentialIndex

!==========================================================================

#ifdef FIXEDBLOCKSIZE
  real,dimension(numPoints) ::xCenter,yCenter,zCenter
#else
  real,allocatable,dimension(:) ::xCenter,yCenter,zCenter
#endif
  real :: dr32, tmpdr32

  integer :: ii

  real :: accel_primary, radius

!==============================================================================
  
  if (.NOT.useGravity) then
     accel(1:numPoints) = 0.0
     return
  end if

#ifndef FIXEDBLOCKSIZE
  allocate(xCenter(numPoints))
  allocate(yCenter(numPoints))
  allocate(zCenter(numPoints))
#endif
  zCenter = 0.
  yCenter = 0.
  if (NDIM == 3) then 
     if (size(kCoords) .GE. numPoints) then
        zCenter(1:numPoints) = kCoords(1:numPoints) - grv_ptzpos
     else
        zCenter(1:numPoints) = kCoords(1) - grv_ptzpos
     end if

  endif
  if (NDIM >= 2) then
     if (size(jCoords) .GE. numPoints) then
        yCenter(1:numPoints) = jCoords(1:numPoints) - grv_ptypos
     else
        yCenter(1:numPoints) = jCoords(1) - grv_ptypos
     end if

  endif
  if (size(iCoords) .GE. numPoints) then
     xCenter = iCoords(1:numPoints) - grv_ptxpos
  else
     xCenter = iCoords(1) - grv_ptxpos
  end if

  if (accelDir .eq. IAXIS) then                       ! x-component
     do ii = 1, numPoints
        tmpdr32 = xCenter(ii)*xCenter(ii) + yCenter(ii)*yCenter(ii)  + zCenter(ii)*zCenter(ii) 
        dr32 = sqrt(tmpdr32)*tmpdr32

        accel(ii) = grv_factor*xCenter(ii)/dr32
     end do

  else if (accelDir .eq. JAXIS) then          ! y-component

     do ii = 1, numPoints
        tmpdr32 = xCenter(ii)*xCenter(ii) + yCenter(ii)*yCenter(ii)  + zCenter(ii)*zCenter(ii) 
        dr32 = sqrt(tmpdr32)*tmpdr32
        
        accel_primary = -grv_newton*sim_mprimary*(yCenter(ii) - sim_yPrim)**(-2)

        accel(ii) = grv_factor*yCenter(ii)/dr32 + accel_primary
     end do

  else if (accelDir .eq. KAXIS) then          ! z-component

     do ii = 1, numPoints
        tmpdr32 = xCenter(ii)*xCenter(ii) + yCenter(ii)*yCenter(ii)  + zCenter(ii)*zCenter(ii) 
        dr32 = sqrt(tmpdr32)*tmpdr32

        accel(ii) = grv_factor*zCenter(ii)/dr32
     end do

  end if

!==============================================================================

  ! Turn off gravity within the sink
  DO ii=1,numPoints
     radius = sqrt(xCenter(ii)**2.0+yCenter(ii)**2.0+zCenter(ii)**2.)
     if (radius .le. sim_radius) then
        accel(ii) = 0.0
     endif
  ENDDO




#ifndef FIXEDBLOCKSIZE
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
#endif

  return

end subroutine Gravity_accelAtCoords
