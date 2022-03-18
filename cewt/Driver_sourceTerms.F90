!!****if* source/Driver/DriverMain/Driver_sourceTerms
!!
!! NAME
!!
!!  Driver_sourceTerms
!!
!! SYNOPSIS
!!
!!  Driver_sourceTerms(integer(IN)::blockCount,
!!                     integer(IN)::blockList(blockCount),
!!                     real(IN) :: dt)
!!
!! DESCRIPTION
!!
!!  Driver for source terms. Instead of calling all these routines 
!!  from Driver_evolveFlash we call Driver_sourceTerms which then
!!  makes the calls to Cool, Burn, Heat and Stir.  If a unit is not
!!  included in the simulation, the routine will be a stub and return
!!  without doing anything.
!! 
!!
!! ARGUMENTS
!!  blockCount   : The number of blocks in the list
!!  blockList    : The list of blocks on which to apply the stirring operator
!!  dt           : the current timestep
!!
!!***



subroutine Driver_sourceTerms(blockCount, blockList, dt, pass)

  use Polytrope_interface, ONLY : Polytrope
  use Driver_data, ONLY: dr_simTime, dr_meshComm
  use Flame_interface, ONLY : Flame_step
  use Stir_interface, ONLY : Stir
  use Heat_interface, ONLY : Heat
  use Heatexchange_interface, ONLY : Heatexchange
  use Burn_interface, ONLY : Burn
  use Cool_interface, ONLY : Cool
  use Ionize_interface, ONLY : Ionize
  use EnergyDeposition_interface, ONLY : EnergyDeposition
  use Deleptonize_interface, ONLY : Deleptonize


  use Grid_interface, ONLY : Grid_getBlkIndexLimits, Grid_getBlkPtr, Grid_releaseBlkPtr,&
       Grid_getCellCoords, Grid_putPointData, Grid_getMinCellSize
  use Eos_interface, ONLY : Eos_wrapped, Eos
  use Eos_data, ONLY: eos_smalle
  use PhysicalConstants_interface, ONLY : PhysicalConstants_get
  use RuntimeParameters_interface, ONLY : RuntimeParameters_mapStrToInt, &
       RuntimeParameters_get
  use Driver_interface, ONLY: Driver_getMype
  use Gravity_data, ONLY: grv_ptmass, grv_factor, grv_newton
  use Simulation_data

  implicit none

#include "Eos.h"
#include "constants.h"
#include "Flash.h"
#include "Flash_mpi.h"


  real, intent(IN)    :: dt
  integer, intent(IN) :: blockCount
  integer, dimension(blockCount), intent(IN):: blockList
  integer, OPTIONAL, intent(IN):: pass

  integer  ::  i, j, k, l, myPE, lb, ierr, cnt, tot_cnt
  real     ::  xx, yy, zz, c2
  real     ::  G, falling_time, relax_rate, delta, weight, local_delta
  real     ::  tau, dist, v, newrho, rho_avg, tot_rho_avg, ke
  real     ::  eint_avg, tot_eint_avg, t_avg, tot_t_avg, p_avg, tot_p_avg
  real     ::  v_avg, tot_v_avg, dist_avg, tot_dist_avg
  real     ::  P, rho, b, bzero, ei, xvel, yvel, zvel, rhotemp,vel
  real     ::  mnew_loc, mold_loc
  real     ::  eintnew_loc,eintold_loc, kenew_loc
  real, dimension(3) :: jnew_loc,radvec,velvec,radxvelvec
  integer  ::  istat
  real :: dvol
  real :: xmomnew_loc  !AA

  real,allocatable,dimension(:) :: xCoord,yCoord,zCoord
  integer,dimension(2,MDIM) :: blkLimits,blkLimitsGC
  integer :: sizeX,sizeY,sizeZ
  real, dimension(:,:,:,:),pointer :: solnData
  real, dimension(EOS_NUM) :: eosData
  integer,dimension(MDIM) :: axis

  logical :: gcell = .true.


  call Driver_getMype(MESH_COMM,myPE)

  ! VERY FIRST THING TO DO: set pt mass as a function of time
  if(dr_simTime .lt. sim_tRamp) then
     !grv_factor = - 0.5d0 * dr_simTime / sim_tRamp
     grv_factor =  -grv_newton * grv_ptmass * ( dr_simTime / sim_tRamp )
  else
     grv_factor =  -grv_newton * grv_ptmass
  endif




  ! get delta 
  call Grid_getMinCellSize(delta)

  
  ! USE AVG VALUES TO UPDATE BULK VALUES?
  if(sim_avgs) then
     !initialize avgs to zero
     rho_avg = 0
     p_avg = 0
     cnt = 0
  
     ! LOOP over ring surrounding sink, compute avg values
     do lb = 1, blockCount
        call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
        sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
        allocate(xCoord(sizeX),stat=istat)
        sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
        allocate(yCoord(sizeY),stat=istat)
        sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
        allocate(zCoord(sizeZ),stat=istat)
        
        if (NDIM == 3) call Grid_getCellCoords&
             (KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
        if (NDIM >= 2) call Grid_getCellCoords&
             (JAXIS, blockList(lb), CENTER,gcell, yCoord, sizeY)
        call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
        call Grid_getBlkPtr(blockList(lb),solnData)
        
     ! Loop z
        do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
           ! if not 3D then set coord to zcenter
           if (NDIM < 3) zCoord(k) = sim_zCtr 
           zz = zCoord(k) - sim_zCtr
           ! Loop y
           do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
              ! if 1D, set yCoord = yCenter
              if (NDIM < 2) yCoord(j) = sim_yCtr
              yy = yCoord(j) - sim_yCtr
              ! Loop x
              do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
                 xx = xCoord(i) - sim_xCtr
                 
                 ! Compute the average values in an annulus surrounding the sink
                 dist = sqrt(xx**2. + yy**2. + zz**2.)
                 if (dist .gt. sim_radius .and. dist .lt. sim_avgs_radius * sim_radius &
                      .and. xx .lt. 0) then
                    rho_avg = rho_avg + solnData(DENS_VAR,i,j,k)
                    p_avg = p_avg + solnData(PRES_VAR,i,j,k)                    
                    cnt = cnt + 1
                 endif !ring
                 
              enddo !xx
           enddo !yy
        enddo !zz
        
        call Grid_releaseBlkPtr(blockList(lb), solnData)
        deallocate(xCoord)
        deallocate(yCoord)
        deallocate(zCoord)
     enddo !blocks
     
     ! Sum all processors, then divide by number of calls to get average
     call MPI_ALLREDUCE(rho_avg, tot_rho_avg, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(p_avg, tot_p_avg, 1, FLASH_REAL, MPI_SUM, MPI_COMM_WORLD, ierr)
     call MPI_ALLREDUCE(cnt, tot_cnt, 1, FLASH_INTEGER, MPI_SUM, MPI_COMM_WORLD, ierr)
     tot_rho_avg = tot_rho_avg / tot_cnt
     tot_p_avg = tot_p_avg / tot_cnt
     

     !print*, myPE, pass, blockCount, tot_cnt

     ! SET THE BULK VALUES BASED ON THE AVERAGE OF THEIR SURROUNDINGS
     sim_rhoBulk = tot_rho_avg * sim_avgs_mult
     sim_pBulk   = tot_p_avg * sim_avgs_mult
  endif

  
  ! NOW LOOP THROUGH AND APPLY BC
  mnew_loc = 0.d0
  mold_loc = 0.d0
  jnew_loc(:) = 0.d0
  eintnew_loc = 0.d0
  eintold_loc = 0.d0
  kenew_loc = 0.d0
  xmomnew_loc = 0.d0  !AA
  
  do lb = 1, blockCount
     call Grid_getBlkIndexLimits(blockList(lb),blkLimits,blkLimitsGC)
     sizeX = blkLimitsGC(HIGH,IAXIS) - blkLimitsGC(LOW,IAXIS) + 1
     allocate(xCoord(sizeX),stat=istat)
     sizeY = blkLimitsGC(HIGH,JAXIS) - blkLimitsGC(LOW,JAXIS) + 1
     allocate(yCoord(sizeY),stat=istat)
     sizeZ = blkLimitsGC(HIGH,KAXIS) - blkLimitsGC(LOW,KAXIS) + 1
     allocate(zCoord(sizeZ),stat=istat)

     if (NDIM == 3) call Grid_getCellCoords&
          (KAXIS, blockList(lb), CENTER, gcell, zCoord, sizeZ)
     if (NDIM >= 2) call Grid_getCellCoords&
          (JAXIS, blockList(lb), CENTER,gcell, yCoord, sizeY)
     call Grid_getCellCoords(IAXIS, blockList(lb), CENTER, gcell, xCoord, sizeX)
     call Grid_getBlkPtr(blockList(lb),solnData)

     ! Loop z                                                         
     do k = blkLimits(LOW, KAXIS), blkLimits(HIGH, KAXIS)
        ! if not 3D then set coord to zcenter                         
        if (NDIM < 3) zCoord(k) = sim_zCtr
        zz = zCoord(k) - sim_zCtr
        ! Loop y                                                     
        do j = blkLimits(LOW, JAXIS), blkLimits(HIGH, JAXIS)
           ! if 1D, set yCoord = yCenter       
           if (NDIM < 2) yCoord(j) = sim_yCtr
           yy = yCoord(j) - sim_yCtr
           ! Loop x                       
           do i = blkLimits(LOW, IAXIS), blkLimits(HIGH, IAXIS)
              xx = xCoord(i) - sim_xCtr
              ! Compute the average values in an annulus surrounding the sink 
              dist = sqrt(xx**2. + yy**2. + zz**2.)
              dvol = delta**3.

              ! bdry var sink
              if (dist .le. sim_radius .and. sim_bctype .eq. 1) then
                 solnData(BDRY_VAR,i,j,k) = 1.0
                 solnData(DENS_VAR,i,j,k) = sim_rhoBulk
                 solnData(PRES_VAR,i,j,k) = sim_pBulk
                 solnData(VELX_VAR,i,j,k) = 0.d0
                 solnData(VELY_VAR,i,j,k) = 0.d0
                 solnData(VELZ_VAR,i,j,k) = 0.d0

              
              ! pressure sink
              elseif (dist .le. sim_radius .and. sim_bctype .eq. 2) then
                 ! integrate the mass inflow
                 mold_loc = mold_loc + sim_rhosink_old*dvol
                 mnew_loc = mnew_loc + solnData(DENS_VAR,i,j,k)*dvol
                  
                 ! angular momentum (old = 0)
                 radvec = (/ xx, yy, zz /)
                 velvec = (/ solnData(VELX_VAR,i,j,k), &
                      solnData(VELY_VAR,i,j,k), &
                      solnData(VELZ_VAR,i,j,k) /)
                 call crossprod(radvec,velvec,radxvelvec)
                 
                 jnew_loc(1) = jnew_loc(1) + solnData(DENS_VAR,i,j,k)*radxvelvec(1)*dvol
                 jnew_loc(2) = jnew_loc(2) + solnData(DENS_VAR,i,j,k)*radxvelvec(2)*dvol
                 jnew_loc(3) = jnew_loc(3) + solnData(DENS_VAR,i,j,k)*radxvelvec(3)*dvol
                 
                 ! internal energy 
                 eintnew_loc = eintnew_loc + &
                      solnData(DENS_VAR,i,j,k)*solnData(EINT_VAR,i,j,k)*dvol
                 eintold_loc = eintold_loc + sim_eintsink_old*sim_rhosink_old*dvol

                 !kinetic energy (old = 0)
                 kenew_loc = kenew_loc + solnData(DENS_VAR,i,j,k)* &
                      (velvec(1)**2 + velvec(2)**2 +velvec(3)**2 )*dvol
                      
                 ! x-momentum (old = 0)
                 xmomnew_loc = xmomnew_loc + &
                      solnData(DENS_VAR,i,j,k)*solnData(VELX_VAR,i,j,k)*dvol
                      
                 ! AA: Changing this to -1.0 since these are flow cells
                 ! when running the sink with the uhd solver.  If we decide
                 ! to use Morgans 'homemade' boundary, must change to 1.0
                 solnData(BDRY_VAR,i,j,k) = -1.0
                 ! floor the pres and density
                 solnData(DENS_VAR,i,j,k) = sim_rhoBulk*sim_presfactor
                 solnData(PRES_VAR,i,j,k) = sim_pBulk*sim_presfactor
                 solnData(VELX_VAR,i,j,k) = 0.d0
                 solnData(VELY_VAR,i,j,k) = 0.d0
                 solnData(VELZ_VAR,i,j,k) = 0.d0
              
              ! not in the sink
              else
                 solnData(BDRY_VAR,i,j,k) = -1.0

                 ! apply gravitational acceleration in the y direction from primary object
                 ! v = a * dt
                 !solnData(VELY_VAR,i,j,k) = solnData(VELY_VAR,i,j,k) + &
                 !     -grv_newton*sim_mprimary*(yCoord(j) - sim_yPrim)**(-2) * dt


              endif
             
              ! update ener var everywhere
              !solnData(ENER_VAR,i,j,k) = solnData(EINT_VAR,i,j,k) + 
              !     0.5*(solnData(VELX_VAR,i,j,k)**2 + solnData(VELY_VAR,i,j,k)**2 + solnData(VELZ_VAR,i,j,k)**2 )


           enddo ! xx
        enddo ! yy
     enddo ! zz
          
     call Grid_releaseBlkPtr(blockList(lb), solnData)        
     deallocate(xCoord)
     deallocate(yCoord)
     deallocate(zCoord)
     ! Calculate total energy, gamma_e, gamma_c ... 
     call Eos_wrapped(MODE_DENS_PRES,blkLimitsGC,blockList(lb))
  enddo
       

  ! REDUCE THE MASS INTEGRATION
  ! mass
  call MPI_ALLREDUCE(mnew_loc, sink_mass_new, 1, FLASH_REAL, MPI_SUM,  dr_meshComm, ierr)
  call MPI_ALLREDUCE(mold_loc, sink_mass_old, 1, FLASH_REAL, MPI_SUM,  dr_meshComm, ierr)
  ! internal energy
  call MPI_ALLREDUCE(eintnew_loc, sink_eint_new, 1, FLASH_REAL, MPI_SUM,  dr_meshComm, ierr)
  call MPI_ALLREDUCE(eintold_loc, sink_eint_old, 1, FLASH_REAL, MPI_SUM,  dr_meshComm, ierr)
  ! kinetic energy
  call MPI_ALLREDUCE(kenew_loc, sink_ke_new, 1, FLASH_REAL, MPI_SUM,  dr_meshComm, ierr)
  ! angular momentum
  call MPI_ALLREDUCE(jnew_loc, sink_j_new, 3, FLASH_REAL, MPI_SUM,  dr_meshComm, ierr)
  !AA: x-momentum
  call MPI_ALLREDUCE(xmomnew_loc, sink_xmom_new, 1, FLASH_REAL, MPI_SUM,  dr_meshComm, ierr)

  mdot_pres = (sink_mass_new - sink_mass_old)/dt
  jdot_pres = sink_j_new/dt
  eint_dot = (sink_eint_new - sink_eint_old)/dt
  ke_dot = sink_ke_new/dt
  xmomdot_pres = sink_xmom_new/dt  !AA
  
  !if(myPE .eq. MASTER_PE) then
  !print*, myPE, pass,  tot_cnt, mnew
  !endif

  ! set the old values... only do this on the second pass (?)
  ! AA:  ! unsplit solver only makes one pass
  if(present(pass)) then
     if(pass .eq. 2) then
        sim_rhosink_old = sim_rhoBulk*sim_presfactor
        sim_psink_old = sim_pBulk*sim_presfactor
        sim_eintsink_old = sim_psink_old/((sim_gamma -1.d0)*sim_rhosink_old) 
     endif
  else
     sim_rhosink_old = sim_rhoBulk*sim_presfactor
     sim_psink_old = sim_pBulk*sim_presfactor
     sim_eintsink_old = sim_psink_old/((sim_gamma -1.d0)*sim_rhosink_old) 
  endif
 
  call Polytrope(blockCount, blockList, dt)
  call Stir(blockCount, blockList, dt) 
  call Flame_step(blockCount, blockList, dt)
  call Burn(blockCount, blockList, dt) 
  call Heat(blockCount, blockList, dt, dr_simTime) 
  call Heatexchange(blockCount, blockList, dt)
  call Cool(blockCount, blockList, dt, dr_simTime)
  call Ionize(blockCount, blockList, dt, dr_simTime)
  call EnergyDeposition(blockCount, blockList, dt, dr_simTime, pass)
  call Deleptonize(blockCount, blockList, dt, dr_simTime)

  return
end subroutine Driver_sourceTerms


! MM: A little cross product routine
subroutine crossprod(a,b,axb)
 
  implicit none
  real,dimension(3),intent(in) :: a
  real,dimension(3),intent(in) :: b 
  real,dimension(3),intent(out) :: axb
  
  axb(1) = a(2)*b(3) - a(3)*b(2)
  axb(2) = a(3)*b(1) - a(1)*b(3)
  axb(3) = a(1)*b(2) - a(2)*b(1)
 
end subroutine crossprod
