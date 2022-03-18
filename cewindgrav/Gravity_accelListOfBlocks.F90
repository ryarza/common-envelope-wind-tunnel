!!****if* source/physics/Gravity/GravityMain/CeWindGrav/Gravity_accelListOfBlocks
!!
!! NAME
!!
!!  Gravity_accelListOfBlocks  
!!
!!
!! SYNOPSIS
!!
!!  Gravity_accelListOfBlocks(integer(IN) :: blockCount,
!!                         integer(IN)    :: blockList(blockCount),
!!                         integer(IN)    :: component,
!!                         integer(IN)    :: accelIndex,
!!                         integer(IN),optional :: potentialIndex)
!!
!! DESCRIPTION
!!
!!  Compute components of the zone-averaged gravitational
!!  acceleration on all mesh blocks.  Either a single component
!!  of the acceleration or all three can be computed.
!!
!! ARGUMENTS
!!
!!   blockCount   - The number of blocks in the list
!!   blockList    - The list of blocks on which to calculate acceleration.
!!   component         - The component of the acceleration to compute.
!!                  Permitted values are IAXIS, JAXIS, KAXIS.  In
!!                     theory, ALLDIR will be implemented.  At the moment
!!                     this routine aborts in a messy way if called with ALLDIR.
!!   accelIndex - variable # to store the acceleration
!!   potentialIndex :   Variable # to take as potential if present, 
!!                     not applicable to pointmass
!!
!!***

subroutine Gravity_accelListOfBlocks (blockCount,blockList,component, &
     accelIndex, potentialIndex)

!==============================================================================

  use Gravity_data, ONLY: grv_ptxpos, grv_ptypos, grv_ptzpos, grv_factor, &
       useGravity, grv_newton
  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, &
    Grid_getCellCoords, Grid_releaseBlkPtr
  use Driver_interface, ONLY: Driver_abortFlash
  use Simulation_data, ONLY : sim_mprimary, sim_yPrim, sim_radius

  implicit none

#include "Flash.h"
#include "constants.h"

  integer,intent(IN)                      :: blockCount
  integer,dimension(MAXBLOCKS), intent(IN)     :: blockList
  integer, INTENT(in) ::  component
  integer, INTENT(in) ::  accelIndex
  integer,intent(IN),optional :: potentialIndex

!==============================================================================



  real          :: dr32, y2, z2
  integer       :: i, j, k, lb
  integer       :: csize
  logical       :: gcell = .true.
#ifdef FIXEDBLOCKSIZE
  real,dimension(GRID_IHI_GC) :: xCenter
  real,dimension(GRID_JHI_GC) :: yCenter
  real,dimension(GRID_KHI_GC) :: zCenter
#else
  real,allocatable,dimension(:) ::xCenter,yCenter,zCenter
#endif
  real, pointer :: solnVec(:,:,:,:)
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC

  real :: accel_primary, radius
!==============================================================================

  if (.NOT.useGravity) return

  do lb = 1, blockCount
     call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
     
#ifndef FIXEDBLOCKSIZE
     allocate(xCenter(blkLimitsGC(HIGH,IAXIS)))
     allocate(yCenter(blkLimitsGC(HIGH,JAXIS)))
     allocate(zCenter(blkLimitsGC(HIGH,KAXIS)))
#endif
     call Grid_getBlkPtr(blockList(lb), solnVec)
!------------------------------------------------------------------------------

     csize = blkLimitsGC(HIGH, IAXIS)
     call Grid_getCellCoords(IAXIS,blockList(lb),CENTER, gcell,xCenter,csize)
     xCenter = xCenter - grv_ptxpos
     yCenter = 0.
     zCenter = 0.
     if (NDIM >= 2) then
        csize = blkLimitsGC(HIGH, JAXIS)
        call Grid_getCellCoords(JAXIS,blockList(lb),CENTER,gcell,yCenter,csize)
        yCenter = yCenter - grv_ptypos
     endif
     if (NDIM == 3) then 
        csize = blkLimitsGC(HIGH, KAXIS)
        call Grid_getCellCoords(KAXIS, blockList(lb), CENTER, gcell, zCenter, csize)
        zCenter = zCenter - grv_ptzpos
     endif

!------------------------------------------------------------------------------
  
     if (component == IAXIS) then                    ! x-axis

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           z2 = zCenter(k)**2
           do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              y2 = yCenter(j)**2
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
                 dr32 = sqrt(xCenter(i)**2 + y2 + z2)
                 radius = dr32
                 dr32 = dr32**3
                 if (radius .gt. sim_radius) then
                    solnVec(accelIndex, i, j, k) = grv_factor*xCenter(i)/dr32
                 else
                    solnVec(accelIndex, i, j, k) = 0.0
                 endif
              enddo
           enddo
        enddo

!------------------------------------------------------------------------------
  
     Elseif (component == JAXIS) then                ! y-axis

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           z2 = zCenter(k)**2
           do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              y2 = yCenter(j)**2
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
                 dr32 = sqrt(xCenter(i)**2 + y2 + z2)
                 radius = dr32
                 dr32 = dr32**3
                 accel_primary = -grv_newton*sim_mprimary*(yCenter(j) - sim_yPrim)**(-2)
                 if (radius .gt. sim_radius) then
                    solnVec(accelIndex, i, j, k) = grv_factor*yCenter(j)/dr32 + accel_primary
                 else
                    solnVec(accelIndex, i, j, k) = 0.0
                 endif
              enddo
           enddo
        enddo

!------------------------------------------------------------------------------
  
     elseif (component == KAXIS) then                ! z-axis

        do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
           z2 = zCenter(k)**2
           do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
              y2 = yCenter(j)**2
              do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)
                 dr32 = sqrt(xCenter(i)**2 + y2 + z2)
                 radius = dr32
                 dr32 = dr32**3
                 if(radius .gt. sim_radius) then
                    solnVec(accelIndex, i, j, k) = grv_factor*zCenter(k)/dr32
                 else
                    solnVec(accelIndex, i, j, k) = 0.0
                 endif
              enddo
           enddo
        enddo

!------------------------------------------------------------------------------
  
     else                                        ! ALLAXIS
        call Driver_abortFlash &
             ('[Gravity_accelListOfBlocks] ALLAXIS not supported!')
     endif

!------------------------------------------------------------------------------
     
#ifndef FIXEDBLOCKSIZE
    deallocate(xCenter)
    deallocate(yCenter)
    deallocate(zCenter)
#endif
    call Grid_releaseBlkPtr(blockList(lb),solnVec)
    
 enddo

!==============================================================================

 return

end subroutine Gravity_accelListOfBlocks
