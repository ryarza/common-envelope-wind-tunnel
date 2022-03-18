!!****if* source/Simulation/SimulationMain/FlatPlate/Simulation_initBlock
!!
!! NAME
!!
!!  Simulation_initBlock
!!
!!
!! SYNOPSIS
!!
!!  Simulation_initBlock(integer(IN) :: blockID) 
!!                       
!!
!!
!!
!! DESCRIPTION
!!
!!  Initializes fluid (initialized with -1.0) and 
!!  stationary rigid body data (initialized with 1.0)
!!
!!
!! ARGUMENTS
!!
!!  blockID -           the number of the block to update
!!***

!!REORDER(4): solnData

subroutine Simulation_initBlock(blockID)

  use Simulation_data
  use Grid_interface,  ONLY : Grid_getBlkIndexLimits, &
                              Grid_getBlkPtr, &
                              Grid_releaseBlkPtr, &
                              Grid_getCellCoords

  implicit none

#include "constants.h"
#include "Flash.h"

  integer,intent(IN) :: blockID
  real,pointer :: solnData(:,:,:,:)

  real :: rho_zone, velx_zone, vely_zone, velz_zone, pres_zone, &
       ener_zone, ekin_zone, eint_zone
  integer :: i, j, k, n, istat, sizeX, sizeY, sizeZ
  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  real :: radius
  integer, dimension(2,MDIM) :: blkLimits, blkLimitsGC
  real, dimension(MDIM) :: del
  real :: xu,xl,yu,yl
  real :: rhotemp
  integer :: jlo
!===============================================================================

  call Grid_getBlkIndexLimits(blockId,blkLimits,blkLimitsGC)

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

  call Grid_getDeltas(blockID,del)

  ! wind in the x direction only
  velx_zone = sim_windVelx
  vely_zone = 0.d0
  velz_zone = 0.d0

  ! Compute the gas energy and set the gamma-values needed for
  ! the equation of state.
  ekin_zone = 0.5 * (velx_zone**2 + vely_zone**2 + velz_zone**2)

! SET A FEW THINGS FOR THE WHOLE BLOCK
  call Grid_getBlkPtr(blockID, solnData, CENTER)
#if NSPECIES > 0
  solnData(SPECIES_BEGIN,:,:,:) =  1.0-(NSPECIES-1)*sim_smallX
  solnData(SPECIES_BEGIN+1:SPECIES_END,:,:,:) =     sim_smallX
#endif

  solnData(VELX_VAR,:,:,:) = velx_zone
  solnData(VELY_VAR,:,:,:) = vely_zone
  solnData(VELZ_VAR,:,:,:) = velz_zone

  solnData(BDRY_VAR,:,:,:) = -1. !! -1 for fluid cells
  
  solnData(GAMC_VAR,:,:,:) = sim_gamma
  solnData(GAME_VAR,:,:,:) = sim_gamma

  ! Now loop and set some position-dependant things
  do k = blkLimitsGC(LOW,KAXIS),blkLimitsGC(HIGH,KAXIS)
     do j = blkLimitsGC(LOW,JAXIS),blkLimitsGC(HIGH,JAXIS)
        do i = blkLimitsGC(LOW,IAXIS),blkLimitsGC(HIGH,IAXIS)

           ! rigid sphere
           radius = sqrt((xCoord(i)-sim_xCtr)**2 + (yCoord(j)-sim_yCtr)**2 + (zCoord(k)-sim_zCtr)**2)
           if (radius .le. sim_radius .and. sim_bctype .eq. 1) then
              rho_zone = sim_rhoBulk
              pres_zone = sim_pBulk
              eint_zone = pres_zone / (sim_gamma-1.)
              eint_zone = eint_zone / rho_zone
              ener_zone = eint_zone ! no velocity
              ener_zone = max(ener_zone, sim_smallP)


              ! we need to overwrite the previous values of vel and bdry
              solnData(VELX_VAR,i,j,k) = 0.d0
              solnData(VELY_VAR,i,j,k) = 0.d0
              solnData(VELZ_VAR,i,j,k) = 0.d0

              solnData(BDRY_VAR,i,j,k) = 1. ! 1 for solid cells
           elseif (radius .le. sim_radius .and. sim_bctype .eq. 2) then
              rho_zone = sim_rhoBulk
              pres_zone = sim_pBulk*sim_presfactor
              eint_zone = pres_zone / (sim_gamma-1.)
              eint_zone = eint_zone / rho_zone
              ener_zone = eint_zone ! no velocity
              ener_zone = max(ener_zone, sim_smallP)


              ! we need to overwrite the previous values of vel and bdry
              solnData(VELX_VAR,i,j,k) = 0.d0
              solnData(VELY_VAR,i,j,k) = 0.d0
              solnData(VELZ_VAR,i,j,k) = 0.d0
              
              solnData(BDRY_VAR,i,j,k) = 1.
              
           else ! exterior
        
              ! introduce a gradient
              !rho_zone = sim_rhoAmbient *exp( -  sim_epsilon_grad * yCoord(j) )
              !pres_zone = sim_pAmbient
              
              ! read from the profiles
              ! get the lower index jlo
              call ut_hunt(sim_yprof,sim_nprof, yCoord(j), jlo)
              ! interpolate to find rho_zone
              call ut_quadraticInterpol(sim_yprof(jlo-1:jlo+1),sim_rhoprof(jlo-1:jlo+1),yCoord(j), rho_zone)
              ! interpolate to find pres_zone
              call ut_quadraticInterpol(sim_yprof(jlo-1:jlo+1),sim_pprof(jlo-1:jlo+1),yCoord(j), pres_zone)
              
              eint_zone = pres_zone / (sim_gamma-1.)
              eint_zone = eint_zone / rho_zone
              ener_zone = eint_zone + ekin_zone
              ener_zone = max(ener_zone, sim_smallP)
           endif
             
           ! store the variables in the block's unk data
           solnData(DENS_VAR,i,j,k) = rho_zone
           solnData(PRES_VAR,i,j,k) = pres_zone
           solnData(ENER_VAR,i,j,k) = ener_zone
#ifdef EINT_VAR
           solnData(EINT_VAR,i,j,k) = eint_zone
#endif
           
        

        end do
     end do
  end do

  call Grid_releaseBlkPtr(blockID, solnData, CENTER)
 

  return
end subroutine Simulation_initBlock



