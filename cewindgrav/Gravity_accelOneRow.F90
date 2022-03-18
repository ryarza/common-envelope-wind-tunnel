!!****if* source/physics/Gravity/GravityMain/CeWindGrav/Gravity_accelOneRow
!!
!! NAME
!!
!!  Gravity_accelOneRow 
!!
!! SYNOPSIS
!!
!!  call Gravity_accelOneRow(integer(IN)  :: pos(2),
!!                           integer(IN)  :: sweepDir,
!!                           integer(IN)  :: blockID,
!!                           integer(IN)  :: numCells,
!!                           real(INOUT)  :: grav(numCells),
!!                           integer(IN),optional :: potentialIndex,
!!                           integer(IN),optional :: extraAccelVars(MDIM))
!!
!! DESCRIPTION
!!
!!  This routine computes the gravitational acceleration for a row
!!  of cells in a specified direction in a given block.
!!
!! ARGUMENTS
!!
!!  pos      :  Row indices transverse to the sweep direction
!!  sweepDir :    The sweep direction:  allowed values are 
!!              SWEEP_X, SWEEP_Y, and SWEEP_Z. These values are defined
!!              in constants.h.
!!  blockID  :  The local identifier of the block to work on
!!  numCells :  Number of cells to update in grav()
!!  grav()   :   Array to receive result
!!  potentialIndex :  optional, not applicable in pointmass gravity
!!  extraAccelVars :  optional, ignored in this implementation
!! 
!!***

subroutine Gravity_accelOneRow (pos, sweepDir, blockID, numCells, grav, &
                                potentialIndex, extraAccelVars)

!=======================================================================

  use Gravity_data, ONLY: grv_ptxpos, grv_ptypos, grv_ptzpos, grv_factor, &
       useGravity, grv_newton
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, &
    Grid_getCellCoords

  use Simulation_data, ONLY : sim_mprimary, sim_yPrim, sim_radius

  implicit none

#include "Flash.h"
#include "constants.h"

  integer, intent(IN) :: sweepDir,blockID,numCells
  integer, dimension(2),INTENT(in) ::pos
  real, dimension(numCells),INTENT(inout) :: grav
  integer,intent(IN),optional :: potentialIndex
  integer,intent(IN),OPTIONAL :: extraAccelVars(MDIM)

!==========================================================================


#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_KHI_GC) :: zCenter
  real,dimension(GRID_JHI_GC) :: yCenter
  real,dimension(GRID_IHI_GC) :: xCenter
#else
  real,allocatable,dimension(:) ::xCenter,yCenter,zCenter
  integer, dimension(LOW:HIGH,MDIM):: blkLimits, blkLimitsGC
#endif
  real :: dr32, tmpdr32

  integer :: sizeX,sizeY,sizez

  integer :: ii,j,k
  logical :: gcell = .true.

  real :: accel_primary, radius

  
!==============================================================================

  if (.NOT.useGravity) return

  j=pos(1)
  k=pos(2)
#ifndef FIXEDBLOCKSIZE
  call Grid_getBlkIndexLimits(blockID,blkLimits,blkLimitsGC)
  sizeX=blkLimitsGC(HIGH,IAXIS)
  sizeY=blkLimitsGC(HIGH,JAXIS)
  sizeZ=blkLimitsGC(HIGH,KAXIS)
  allocate(xCenter(sizeX))
  allocate(yCenter(sizeY))
  allocate(zCenter(sizeZ))
#else
  sizeX=GRID_IHI_GC
  sizeY=GRID_JHI_GC
  sizeZ=GRID_KHI_GC
#endif
  zCenter = 0.
  yCenter = 0.
  if (NDIM == 3) then 
     call Grid_getCellCoords(KAXIS, blockID, CENTER, gcell, zCenter, sizeZ)
     zCenter = zCenter - grv_ptzpos
  endif
  if (NDIM >= 2) then
     call Grid_getCellCoords(JAXIS, blockID, CENTER, gcell, yCenter, sizeY)
     yCenter = yCenter - grv_ptypos
  endif
  call Grid_getCellCoords(IAXIS, blockID, CENTER, gcell, xCenter, sizeX)
  xCenter = xCenter - grv_ptxpos
  

  if (sweepDir .eq. SWEEP_X) then                       ! x-component

     tmpdr32 = yCenter(j)*yCenter(j) + zCenter(k)*zCenter(k) 

     do ii = 1, numCells

        dr32 = sqrt(xCenter(ii)*xCenter(ii) + tmpdr32)
        dr32 = dr32*dr32*dr32

        grav(ii) = grv_factor*xCenter(ii)/dr32
     enddo


  else if (sweepDir .eq. SWEEP_Y) then          ! y-component

     tmpdr32 = xCenter(j)*xCenter(j) + zCenter(k)*zCenter(k) 

     do ii = 1, numCells
        
        dr32 = sqrt(yCenter(ii)*yCenter(ii) + tmpdr32)
        dr32 = dr32*dr32*dr32

        accel_primary = -grv_newton*sim_mprimary*(yCenter(ii) - sim_yPrim)**(-2)

        grav(ii) = grv_factor*yCenter(ii)/dr32 + accel_primary
        
     enddo

  else if (sweepDir .eq. SWEEP_Z) then          ! z-component

     tmpdr32 = xCenter(j)*xCenter(j) + yCenter(k)*yCenter(k) 

     do ii = 1, numCells
        
        dr32 = sqrt(zCenter(ii)*zCenter(ii) + tmpdr32)           
        dr32 = dr32*dr32*dr32
        
        grav(ii) = grv_factor*zCenter(ii)/dr32
     enddo

  endif

!==============================================================================

  ! Turn off gravity within the sink
  DO ii=1,numCells
     radius = sqrt(xCenter(ii)**2.0+yCenter(ii)**2.0+zCenter(ii)**2.)
     if (radius .le. sim_radius) then
        grav(ii) = 0.0
     endif
  ENDDO


#ifndef FIXEDBLOCKSIZE
  deallocate(xCenter)
  deallocate(yCenter)
  deallocate(zCenter)
#endif

  return

end subroutine Gravity_accelOneRow
