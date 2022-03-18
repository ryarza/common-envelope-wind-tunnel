!!****if* source/IO/IOMain/IO_writeIntegralQuantities
!!
!!
!!  NAME
!!    IO_writeIntegralQuantities
!!
!!  SYNOPSIS
!!    call IO_writeIntegralQuantities(integer(in) :: isFirst,
!!                                    real(in)    :: simTime)
!!
!!  DESCRIPTION
!!
!!   Compute the values of integral quantities (eg. total energy)
!!   and write them to an ASCII file.  If this is the initial step,
!!   create the file and write a header to it before writing the data.
!!
!!   Presently, this supports 1, 2, and 3-d Cartesian geometry and 2-d
!!   cylindrical geometry (r,z).  More geometries can be added by
!!   modifying the volume of each zone (dvol).
!!
!!   Users should modify this routine if they want to store any
!!   quantities other than default values in the flash.dat.  Make sure
!!   to modify the nGlobalSum parameter to match the number of
!!   quantities written.  Also make sure to modify the header to match
!!   the names of quantities with those calculated in the lsum and
!!   gsum arrays.
!!  
!!  ARGUMENTS
!!    
!!   isFirst - if 1 then write header info plus data, otherwise just write data
!!   simTime - simulation time
!!
!!
!!***

!!REORDER(4):solnData

subroutine IO_writeIntegralQuantities ( isFirst, simTime)

  use IO_data, ONLY : io_restart, io_statsFileName, io_globalComm
  use Grid_interface, ONLY : Grid_getListOfBlocks, &
    Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_getSingleCellVol, &
    Grid_releaseBlkPtr

   use IO_data, ONLY : io_globalMe
   use Simulation_data, ONLY : mdot_pres, sim_gcell, jdot_pres, eint_dot, ke_dot, &
       sim_radius, sim_xCtr, sim_yCtr, sim_zCtr, & !AA
       df_radii, df_radii_labels, NDFPOINTS, xmomdot_pres  !AA
   
   use Gravity_data, ONLY : grv_newton, grv_factor !AA

  implicit none

#include "Flash_mpi.h"
#include "constants.h"
#include "Flash.h"
  
  
  real, intent(in) :: simTime

  integer, intent(in) :: isFirst

  integer :: lb, count
  
  integer :: funit = 99
  integer :: error
  
  character (len=MAX_STRING_LENGTH), save :: fname 
  
  integer :: blockList(MAXBLOCKS)

  integer :: blkLimits(HIGH, MDIM), blkLimitsGC(HIGH, MDIM)

  ! MM increasing nGlobalSum from 8
  !AA: increasing nGlobalSum to accomodate df_x_i
  !RY: increase nGlobalSum by 3 to accomodate for each component of the total pressure force on the object
  integer, parameter ::  nGlobalSum = 12 + NDFPOINTS + 3 ! Number of globally-summed quantities
  real :: gsum(nGlobalSum) !Global summed quantities
  real :: lsum(nGlobalSum) !Global summed quantities

  !RY: PFXIDX = Index (IDX) within the lsum and gsum arrays of the pressure (P) force (F) in the x (X) direction. The first 12+NDFPOINTS elements are occupied so we start with that+1.
  integer, parameter :: PFXIDX = 12 + NDFPOINTS + 1
  !RY: Just do +1 to make sure they're contiguous
  integer, parameter :: PFYIDX = PFXIDX + 1
  integer, parameter :: PFZIDX = PFYIDX + 1

  integer :: i, j, k
  real :: dvol             , del(MDIM)
  real, DIMENSION(:,:,:,:), POINTER :: solnData

  integer :: point(MDIM)
  integer :: ioStat

  ! MM defining a few variables
  real :: mdotface
  real,dimension(3) :: rad,vel,radxvel
  integer :: sizeX, sizeY, sizeZ, istat, blockID
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord

  real :: radius

  !AA
  real :: xx, yy, zz
  real :: df_drag_temp = 0.0
  integer :: df_i

  ! Sum quantities over all locally held leaf-node blocks.
  gsum  = 0.
  lsum = 0.
  
  call Grid_getListOfBlocks(LEAF, blockList, count)
  
  do lb = 1, count
     !get the index limits of the block
     blockID = blocklist(lb)
     call Grid_getBlkIndexLimits(blockID, blkLimits, blkLimitsGC)

     !MM:  get the coordinates
     sizeX = blkLimitsGC(HIGH,IAXIS)-blkLimitsGC(LOW,IAXIS)+1
     sizeY = blkLimitsGC(HIGH,JAXIS)-blkLimitsGC(LOW,JAXIS)+1
     sizeZ = blkLimitsGC(HIGH,KAXIS)-blkLimitsGC(LOW,KAXIS)+1
     
     allocate(xCoord(sizeX),stat=istat)
     allocate(yCoord(sizeY),stat=istat)
     allocate(zCoord(sizeZ),stat=istat)

     xCoord = 0.0
     yCoord = 0.0
     zCoord = 0.0
     

     if (NDIM == 3) call Grid_getCellCoords(KAXIS,blockID,CENTER,sim_gCell,zCoord,sizeZ)
     if (NDIM >= 2) call Grid_getCellCoords(JAXIS,blockID,CENTER,sim_gCell,yCoord,sizeY)
     call Grid_getCellCoords(IAXIS,blockID,CENTER,sim_gCell,xCoord,sizeX)


     ! get a pointer to the current block of data
     call Grid_getBlkPtr(blockID, solnData)

     ! Sum contributions from the indicated blkLimits of cells.
     do k = blkLimits(LOW,KAXIS), blkLimits(HIGH,KAXIS)
      ! AA: if not 3d, zz = 0
     if (NDIM < 3) zCoord(k) = sim_zCtr
     zz = zCoord(k) - sim_zCtr
        do j = blkLimits(LOW,JAXIS), blkLimits(HIGH,JAXIS)
        ! AA: if 1D, yy = 0
        if (NDIM < 2) yCoord(j) = sim_yCtr
        yy = yCoord(j) - sim_yCtr
           do i = blkLimits(LOW,IAXIS), blkLimits(HIGH,IAXIS)
           xx = xCoord(i) - sim_xCtr
              
              point(IAXIS) = i
              point(JAXIS) = j
              point(KAXIS) = k

!! Get the cell volume for a single cell
              call Grid_getSingleCellVol(blockID, EXTERIOR, point, dvol)
!! MM: also get the deltas
              call Grid_getDeltas(blockID,del)

              ! mass   
#ifdef DENS_VAR
              lsum(1) = lsum(1) + solnData(DENS_VAR,i,j,k)*dvol 
#endif           


#ifdef DENS_VAR
#ifdef VELX_VAR      
              ! momentum
              lsum(2) = lsum(2) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELX_VAR,i,j,k)*dvol 
           
#endif
#ifdef VELY_VAR      

              lsum(3) = lsum(3) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELY_VAR,i,j,k)*dvol
           
#endif
#ifdef VELZ_VAR      
              lsum(4) = lsum(4) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(VELZ_VAR,i,j,k)*dvol
#endif

              ! total energy
#ifdef ENER_VAR
              lsum(5) = lsum(5) + solnData(ENER_VAR,i,j,k) * & 
                   &                                solnData(DENS_VAR,i,j,k)*dvol
#ifdef MAGP_VAR
              ! total plasma energy
!!$              lsum(5) = lsum(5) + (solnData(ENER_VAR,i,j,k) * & 
!!$                   &    solnData(DENS_VAR,i,j,k) + solnData(MAGP_VAR,i,j,k))*dvol

              lsum(5) = lsum(5) + solnData(MAGP_VAR,i,j,k)*dvol
#endif
#endif

           
#ifdef VELX_VAR      
#ifdef VELY_VAR      
#ifdef VELZ_VAR      
              ! kinetic energy
              lsum(6) = lsum(6) + 0.5*solnData(DENS_VAR,i,j,k) * & 
                   &                             (solnData(VELX_VAR,i,j,k)**2+ & 
                   &                              solnData(VELY_VAR,i,j,k)**2+ & 
                   &                              solnData(VELZ_VAR,i,j,k)**2)*dvol           

#endif
#endif
#endif


#ifdef EINT_VAR
              ! internal energy
              lsum(7) = lsum(7) + solnData(DENS_VAR,i,j,k) * & 
                   &                                solnData(EINT_VAR,i,j,k)*dvol
#endif
#endif ! ifdef DENS_VAR

#ifdef MAGP_VAR
              ! magnetic energy
              lsum(8) = lsum(8) + solnData(MAGP_VAR,i,j,k)*dvol
#endif


! MM mass accreted into the inner BC
#ifdef BDRY_VAR
              !if (io_globalMe  == MASTER_PE) print*, radxvel
              mdotface = 0.d0
              rad = 0.d0
              vel = 0.d0
              radxvel = 0.d0

              ! CASES of face flux
              !xleft
              if(solnData(BDRY_VAR,i,j,k) > 0 .and. solnData(BDRY_VAR,i-1,j,k) < 0 ) then

                 !RY: Add pressure times area (i.e. force). Pressure on this face will give a force to the right. Add regardless of sign(velocity)
                 lsum(PFXIDX) = lsum(PFXIDX) + solnData(PRES_VAR, i - 1, j, k) * del(JAXIS) * del(KAXIS)

                 ! if flowing towards the origin
                 if(solnData(VELX_VAR,i-1,j,k) > 0) then
                    ! add mdot through that face
                    mdotface = solnData(DENS_VAR,i-1,j,k)*del(JAXIS)*del(KAXIS)*solnData(VELX_VAR,i-1,j,k)
                    
                    ! mdot
                    lsum(9) = lsum(9) + mdotface

                    ! set the radius and velocity and do the cross product here
                    rad = (/ xCoord(i-1),yCoord(j),zCoord(k) /)
                    vel = (/ solnData(VELX_VAR,i-1,j,k), solnData(VELY_VAR,i-1,j,k), solnData(VELZ_VAR,i-1,j,k) /)
                    call cross(rad,vel,radxvel)

                    ! jx,jy,jz
                    lsum(10) = lsum(10) + mdotface*radxvel(1)
                    lsum(11) = lsum(11) + mdotface*radxvel(2)
                    lsum(12) = lsum(12) + mdotface*radxvel(3)


                 end if
              end if
              
              !xright
              if(solnData(BDRY_VAR,i,j,k) > 0 .and. solnData(BDRY_VAR,i+1,j,k) < 0 ) then

                 !RY: Add pressure times area (i.e. force). Pressure on this face will give a force to the left
                 lsum(PFXIDX) = lsum(PFXIDX) - solnData(PRES_VAR, i + 1, j, k) * del(JAXIS) * del(KAXIS)

                 if(solnData(VELX_VAR,i+1,j,k) < 0) then
                    ! minus sign is because inward radial flow corresponds to VEL < 0
                    mdotface = - solnData(DENS_VAR,i+1,j,k)*del(JAXIS)*del(KAXIS)*solnData(VELX_VAR,i+1,j,k)
                    
                    ! mdot
                    lsum(9) = lsum(9) + mdotface

                    ! set the radius and velocity and do the cross product here
                    rad = (/ xCoord(i+1),yCoord(j),zCoord(k) /)
                    vel = (/ solnData(VELX_VAR,i+1,j,k), solnData(VELY_VAR,i+1,j,k), solnData(VELZ_VAR,i+1,j,k) /)
                    call cross(rad,vel,radxvel)

                    ! jx,jy,jz
                    lsum(10) = lsum(10) + mdotface*radxvel(1)
                    lsum(11) = lsum(11) + mdotface*radxvel(2)
                    lsum(12) = lsum(12) + mdotface*radxvel(3)


                 end if
              end if

              !yleft
              if(solnData(BDRY_VAR,i,j,k) > 0 .and. solnData(BDRY_VAR,i,j-1,k) < 0 ) then

                !RY: Add pressure times area (i.e. force). Pressure on this face will give a force in the +y direction
                lsum(PFYIDX) = lsum(PFYIDX) + solnData(PRES_VAR, i, j - 1, k) * del(IAXIS) * del(KAXIS)

                if(solnData(VELY_VAR,i,j-1,k) > 0) then
                   mdotface =  solnData(DENS_VAR,i,j-1,k)*del(IAXIS)*del(KAXIS)*solnData(VELY_VAR,i,j-1,k)
                   
                   ! mdot
                   lsum(9) = lsum(9) + mdotface

                   ! set the radius and velocity and do the cross product here
                   rad = (/ xCoord(i),yCoord(j-1),zCoord(k) /)
                   vel = (/ solnData(VELX_VAR,i,j-1,k), solnData(VELY_VAR,i,j-1,k), solnData(VELZ_VAR,i,j-1,k) /)
                   call cross(rad,vel,radxvel)

                   ! jx,jy,jz
                   lsum(10) = lsum(10) + mdotface*radxvel(1)
                   lsum(11) = lsum(11) + mdotface*radxvel(2)
                   lsum(12) = lsum(12) + mdotface*radxvel(3)
 

                 end if
              end if

              !yright
              if(solnData(BDRY_VAR,i,j,k) > 0 .and. solnData(BDRY_VAR,i,j+1,k) < 0 ) then

                 !RY: Add pressure times area (i.e. force). Pressure on this face will give a force in the +y direction
                 lsum(PFYIDX) = lsum(PFYIDX) - solnData(PRES_VAR, i, j + 1, k) * del(IAXIS) * del(KAXIS)

                 if(solnData(VELY_VAR,i,j+1,k) < 0) then
                    ! minus for rhs
                    mdotface =  - solnData(DENS_VAR,i,j+1,k)*del(IAXIS)*del(KAXIS)*solnData(VELY_VAR,i,j+1,k)
                    
                    !mdot
                    lsum(9) = lsum(9) + mdotface

                    ! set the radius and velocity and do the cross product here
                    rad = (/ xCoord(i),yCoord(j+1),zCoord(k) /)
                    vel = (/ solnData(VELX_VAR,i,j+1,k), solnData(VELY_VAR,i,j+1,k), solnData(VELZ_VAR,i,j+1,k) /)
                    call cross(rad,vel,radxvel)

                    ! jx,jy,jz
                    lsum(10) = lsum(10) + mdotface*radxvel(1)
                    lsum(11) = lsum(11) + mdotface*radxvel(2)
                    lsum(12) = lsum(12) + mdotface*radxvel(3)


                 end if
              end if
              
              !zleft
              if(solnData(BDRY_VAR,i,j,k) > 0 .and. solnData(BDRY_VAR,i,j,k-1) < 0 ) then

                 !RY: Add pressure times area (i.e. force). Pressure on this face will give a force in the +y direction
                 lsum(PFZIDX) = lsum(PFZIDX) + solnData(PRES_VAR, i, j, k - 1) * del(IAXIS) * del(JAXIS)

                 if(solnData(VELZ_VAR,i,j,k-1) > 0) then
                    mdotface =  solnData(DENS_VAR,i,j,k-1)*del(IAXIS)*del(JAXIS)*solnData(VELZ_VAR,i,j,k-1)
                    
                    ! mdot
                    lsum(9) = lsum(9) + mdotface

                    ! set the radius and velocity and do the cross product here
                    rad = (/ xCoord(i),yCoord(j),zCoord(k-1) /)
                    vel = (/ solnData(VELX_VAR,i,j,k-1), solnData(VELY_VAR,i+1,j,k-1), solnData(VELZ_VAR,i+1,j,k-1) /)
                    call cross(rad,vel,radxvel)

                    ! jx,jy,jz
                    lsum(10) = lsum(10) + mdotface*radxvel(1)
                    lsum(11) = lsum(11) + mdotface*radxvel(2)
                    lsum(12) = lsum(12) + mdotface*radxvel(3)


                 end if
              end if

              !zright
              if(solnData(BDRY_VAR,i,j,k) > 0 .and. solnData(BDRY_VAR,i,j,k+1) < 0 ) then

                 !RY: Add pressure times area (i.e. force). Pressure on this face will give a force in the +y direction
                 lsum(PFZIDX) = lsum(PFZIDX) - solnData(PRES_VAR, i, j, k + 1) * del(IAXIS) * del(JAXIS)

                 if(solnData(VELZ_VAR,i,j,k+1) < 0) then
                    ! minus for rhs
                    mdotface =  - solnData(DENS_VAR,i,j,k+1)*del(IAXIS)*del(JAXIS)*solnData(VELZ_VAR,i,j,k+1)
                    
                    ! mdot
                    lsum(9) = lsum(9) + mdotface

                    ! set the radius and velocity and do the cross product here
                    rad = (/ xCoord(i),yCoord(j),zCoord(k+1) /)
                    vel = (/ solnData(VELX_VAR,i,j,k+1), solnData(VELY_VAR,i,j,k+1), solnData(VELZ_VAR,i,j,k+1) /)
                    call cross(rad,vel,radxvel)

                    ! jx, jy, jz
                    lsum(10) = lsum(10) + mdotface*radxvel(1)
                    lsum(11) = lsum(11) + mdotface*radxvel(2)
                    lsum(12) = lsum(12) + mdotface*radxvel(3)


                 end if
              end if

#endif
              
              ! AA: dynamical friction force - This is the gravitational drag force 
              ! on the particle (by the gas) in x-hat direction
              
              radius  = sqrt(xx**2 + yy**2 + zz**2)
              
              if( (radius .gt. sim_radius) .and. (radius .le. df_radii(NDFPOINTS)) ) then

                 df_drag_temp = -grv_factor*solnData(DENS_VAR,i,j,k)*dvol * xx / radius**3	

                 do df_i = 1, NDFPOINTS
                    if( radius .le. df_radii(df_i) ) then
                       lsum(df_i + 12) = lsum(df_i + 12) + df_drag_temp
                    endif
                 enddo
               endif



           enddo
        enddo
     enddo
     call Grid_releaseBlkPtr(blockList(lb), solnData)

  enddo
  

  
  ! Now the MASTER_PE sums the local contributions from all of
  ! the processors and writes the total to a file.
  
  call MPI_Reduce (lsum, gsum, nGlobalSum, FLASH_REAL, MPI_SUM, & 
       &                MASTER_PE, io_globalComm, error)
  

  if (io_globalMe  == MASTER_PE) then
     
     ! create the file from scratch if it is a not a restart simulation, 
     ! otherwise append to the end of the file
     
     !No mater what, we are opening the file. Check to see if already there
     ioStat = 0
     open(funit, file=trim(io_statsFileName), position='APPEND', status='OLD', iostat=ioStat)
     if(ioStat .NE. 0) then
        !print *, 'FILE FOUND'
        open(funit, file=trim(io_statsFileName), position='APPEND')
     endif
     
     if (isFirst .EQ. 1 .AND. (.NOT. io_restart .or. ioStat .NE. 0)) then
        
#ifndef MAGP_VAR
        write (funit, 10)               &
             '#time                     ', &
             'mass                      ', &
             'total-mom-x               ', &
             'total-mom-y               ', & 
             'total-mom-z               ', &
             'E_total                   ', &
             'E_kinetic                 ', &
             'E_internal                ', &
             'nothing                   ', &
             'mdot_accr_flux            ', &
             'jdot_accr_x_flux          ', &
             'jdot_accr_y_flux          ', &
             'jdot_accr_z_flux          ', &
             df_radii_labels(:)          , &  !AA
             'pres_force_bdry_x         ', &  !RY
             'pres_force_bdry_y         ', &  !RY
             'pres_force_bdry_z         ', &  !RY
             'mdot_accr_integ           ', &
             'jdot_accr_x_integ         ', &
             'jdot_accr_y_integ         ', &
             'jdot_accr_z_integ         ', &
             'eintdot_accr              ', &
             'kedot_accr                ', &
             'momdrag_accr_x            '    !AA
#endif
        
#ifdef MAGP_VAR
        write (funit, 10)               &
             '#time                     ', &
             'mass                      ', &
             'x-momentum                ', &
             'y-momentum                ', & 
             'z-momentum                ', &
             'E_total                   ', &
             'E_kinetic                 ', &
             'E_internal                ', &
             'MagEnergy                 '
#endif
        
10         format (2x,50(a25, :, 1X))

     else if(isFirst .EQ. 1) then
        write (funit, 11) 
11      format('# simulation restarted')
     endif
     
     ! Write the global sums to the file
     write (funit, 12) simtime, gsum, mdot_pres, &
          jdot_pres(1),jdot_pres(2),jdot_pres(3), eint_dot, ke_dot, &
          xmomdot_pres !AA
12   format (1x, 50(es25.18, :, 1x))
 
     close (funit)          ! Close the file.
     
  endif
  
  call MPI_Barrier (io_globalComm, error)
  
  !=============================================================================
  
  return
end subroutine IO_writeIntegralQuantities


! MM: A little cross product routine
subroutine cross(a,b,axb)
 
  implicit none
  real,dimension(3),intent(in) :: a
  real,dimension(3),intent(in) :: b 
  real,dimension(3),intent(out) :: axb
  
  axb(1) = a(2)*b(3) - a(3)*b(2)
  axb(2) = a(3)*b(1) - a(1)*b(3)
  axb(3) = a(1)*b(2) - a(2)*b(1)
 
end subroutine cross
