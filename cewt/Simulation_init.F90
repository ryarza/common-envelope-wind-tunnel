!!****if* source/Simulation/SimulationMain/FlatPlate/Simulation_init
!!
!! NAME
!!
!!  Simulation_init
!!
!!
!! SYNOPSIS
!!
!!  Simulation_init()
!!
!!
!! DESCRIPTION
!!
!!  Initializes all the parameters needed for the wind tunnel with a step problem
!!
!! ARGUMENTS
!!
!!  
!!
!! PARAMETERS
!!
!!  sim_pAmbient    Initial ambient pressure
!!  sim_rhoAmbient  Initial ambient density
!!  sim_windVel     Inflow velocity (parallel to x-axis)
!!  gamma           the Gamma EOS thing
!!  smallp          minimum for pressure
!!  smallx          minimum for abundance
!!
!!***

subroutine Simulation_init()

  use Simulation_data
  use RuntimeParameters_interface, ONLY : RuntimeParameters_get
  use Logfile_interface, ONLY : Logfile_stamp
  use grid_data, ONLY : gr_globalMe
  use Gravity_data, ONLY: grv_ptmass,grv_newton
  use Grid_data, ONLY: gr_imin,gr_imax,gr_jmin,gr_jmax,gr_kmin,gr_kmax  !AA

  implicit none

#include "constants.h"
#include "Flash.h"

  ! a few local variables
  real :: soundspeed,vwind, racc, Ggrav
  real :: presBase, dp, dr
  integer :: i, mid
  integer :: ioStat !AA
  real :: df_max_radius, df_min_radius    !AA
  character(10) :: label, string, prefix = "df_x_"  !AA

  ! integration variables for rk4 HSE profile
  real, dimension(2) :: rpvec, rpvec0, k1,k2,k3,k4, deriv_loc
  
  call RuntimeParameters_get('sim_rhoAmbient', sim_rhoAmbient)
  call RuntimeParameters_get('sim_rhoBulk',sim_rhoBulk)
  call RuntimeParameters_get('sim_pBulk',sim_pBulk)
  call RuntimeParameters_get('gamma', sim_gamma)
  call RuntimeParameters_get('smallx', sim_smallX)
  call RuntimeParameters_get('smallp', sim_smallP)
  call RuntimeParameters_get('sim_radius', sim_radius)
  call RuntimeParameters_get('sim_xCtr', sim_xCtr)
  call RuntimeParameters_get('sim_yCtr', sim_yCtr)
  call RuntimeParameters_get('sim_zCtr', sim_zCtr)
  !call RuntimeParameters_get('sim_Mach', sim_Mach)
  call RuntimeParameters_get("sim_bctype", sim_bctype)
  call RuntimeParameters_get("sim_presfactor", sim_presfactor)
  call RuntimeParameters_get("sim_epsilon_grad",sim_epsilon_grad)
  call RuntimeParameters_get("sim_max_ref_rad",sim_max_ref_rad)
  call RuntimeParameters_get("sim_grad_ymax",sim_grad_ymax)
  call RuntimeParameters_get("sim_avgs",sim_avgs)
  call RuntimeParameters_get("sim_avgs_radius",sim_avgs_radius)
  call RuntimeParameters_get("sim_avgs_mult",sim_avgs_mult)
  call RuntimeParameters_get("sim_tRamp",sim_tRamp)
  call RuntimeParameters_get("sim_q",sim_q)
  call RuntimeParameters_get("sim_gamma_structure",sim_gamma_structure)
  call RuntimeParameters_get("sim_rhoBackground",sim_rhoBackground)
  call RuntimeParameters_get("sim_fkep",sim_fkep)

  call RuntimeParameters_get("ptmass", grv_ptmass)
  !print*,"grv_ptmass =", grv_ptmass
  
  Ggrav = 6.67428e-8

  sim_gCell = .true.


  ! set the mach number  based on the HSE condition. 
  sim_mach = sqrt( (1+sim_q)**2/(2.*sim_q)  * sim_epsilon_grad * (sim_gamma_structure/sim_gamma) ) * sim_fkep**2
  if (gr_globalMe  == MASTER_PE) then
    print *, 'sim_q', sim_q
    print *, 'sim_q', sim_epsilon_grad
    print *, 'sim_q', sim_gamma_structure
    print *, 'sim_q', sim_gamma
    print *, 'sim_q', sim_fkep
  endif 

  ! set pAmbient such that cs = 1/sim_mach
  ! this specifies that the time unit be Ra=1
  sim_pAmbient = sim_rhoAmbient/(sim_gamma*sim_mach*sim_mach)
  

  ! set the Background values
  sim_pBackground = sim_pAmbient*(sim_rhoBackground/sim_rhoAmbient)**(sim_gamma_structure)

  !set wind only in the x-direction 
  soundspeed = sqrt(sim_gamma * sim_pAmbient / sim_rhoAmbient)
  vwind = sim_Mach*soundspeed
  sim_windVelx = vwind
  sim_windVely = 0.d0
  sim_windVelz = 0.d0

  !binary parameters
  sim_mprimary = grv_ptmass * sim_q**(-1)
  sim_a = 0.5 * sim_fkep**2 * (1.d0 + sim_q**(-1) )
  sim_yPrim = sim_yCtr - sim_a

 
  ! establish the y-direction profile
  ! first the coordinate & density
  do i=1,sim_nprof
     ! make the initial profile slightly wider than the domain
     sim_yprof(i) = (i-1.0)/(sim_nprof-1.0) * 1.25*(gr_jmax-gr_jmin) + 1.25*gr_jmin
  end do

  
  ! central value
  mid = (sim_nprof-1)/2+1
  sim_pprof(mid)    = sim_pAmbient
  sim_rhoprof(mid)  = sim_rhoAmbient
  sim_erhoprof(mid) = sim_epsilon_grad


  ! integrate the HSE pressure profile
  dr = sim_yprof(2) - sim_yprof(1)
  ! first from midpoint to upper boundary
  rpvec0 = (/ sim_rhoprof(mid), sim_pprof(mid) /)
  do i=mid+1,sim_nprof
     k1 = dr*derivs( sim_yprof(i-1), rpvec0 )
     k2 = dr*derivs( sim_yprof(i-1) + dr/2., rpvec0 + k1/2. )
     k3 = dr*derivs( sim_yprof(i-1) + dr/2., rpvec0 + k2/2. )
     k4 = dr*derivs( sim_yprof(i-1) + dr, rpvec0 + k3 )
     
     rpvec = rpvec0 + (k1/6. + k2/3. + k3/3. + k4/6.)
     
     sim_rhoprof(i) = rpvec(1)
     sim_pprof(i)   = rpvec(2)
     deriv_loc = derivs(sim_yprof(i) , rpvec)
     sim_erhoprof(i) = -1.0/sim_rhoprof(i)  * deriv_loc(1)

     rpvec0 = rpvec
  end do
  
  ! then backwards to the lower boundary
  dr = -dr 
  rpvec0 = (/ sim_rhoprof(mid), sim_pprof(mid) /)
  do i=mid-1,1,-1
     k1 = dr*derivs( sim_yprof(i+1), rpvec0 )
     k2 = dr*derivs( sim_yprof(i+1) + dr/2., rpvec0 + k1/2. )
     k3 = dr*derivs( sim_yprof(i+1) + dr/2., rpvec0 + k2/2. )
     k4 = dr*derivs( sim_yprof(i+1) + dr, rpvec0 + k3 )
     
     rpvec = rpvec0 + (k1/6. + k2/3. + k3/3. + k4/6.)
     
     sim_rhoprof(i) = rpvec(1)
     sim_pprof(i)   = rpvec(2)
     deriv_loc = derivs(sim_yprof(i) , rpvec)
     sim_erhoprof(i) = -1.0/sim_rhoprof(i)  * deriv_loc(1)

     rpvec0 = rpvec
  end do

  ! clean up the profile by applying background values
  do i=1,sim_nprof
     if ( (sim_pprof(i) < sim_pBackground) .or. (sim_rhoprof(i) < sim_rhoBackground) ) then
        sim_rhoprof(i) = sim_rhoBackground
        sim_pprof(i) = sim_pBackground
        sim_erhoprof(i) = 0.0
     end if
  end do
  

  ! init xmomentum sink variables
  sink_xmom_new = 0.d0 !AA
  xmomdot_pres = 0.d0 !AA
  
  ! AA: fill df_radii with NDFPOINTS log spaced points
  ! from 2.0*sim_radius and to the closest domain boundary
  
  df_max_radius = MIN(gr_imax,gr_jmax,gr_kmax,&
                  ABS(gr_imin),ABS(gr_jmin),ABS(gr_kmin))
  df_min_radius = 2.0 * sim_radius
  
  call logspaced(df_min_radius,df_max_radius,NDFPOINTS,df_radii)
  
  ! make the array of df column names
  do i = 1, NDFPOINTS 
     write (string, "(I10)") i
     label = trim(prefix)//trim(ADJUSTL(string))
     df_radii_labels(i) = label	
  enddo



  ! PRINT some information to stdout and file
  if (gr_globalMe  == MASTER_PE) then

     ! write the profile to file
     open(unit=73,file="profile.dat")
     write(73,"(7a20)") "i","y","rho","pres","cs","mach", "erho"

     do i=1,sim_nprof
        write(73,"(i20,6f20.15)") i, sim_yprof(i), &
             sim_rhoprof(i), sim_pprof(i), sqrt(sim_gamma * sim_pprof(i)/sim_rhoprof(i)), &
             vwind / sqrt(sim_gamma * sim_pprof(i)/sim_rhoprof(i)), &
             sim_erhoprof(i)
        
     enddo

     close(73)
     
     
     ! write df_radii to file
     open (40, file="df_sum_radii.dat", iostat=ioStat,  & 
        status='unknown', position='append')
     write(40,*)  "i     df_radius"
     do i = 1, NDFPOINTS
        write(40,*) df_radii_labels(i),"    ", df_radii(i)
     enddo
     close (40)


     racc = 2 * Ggrav * grv_ptmass / vwind**2

     print*,"##############################################"
     print*,"sim_rhoAmbient = ", sim_rhoAmbient
     print*,"sim_pAmbient = ",sim_pAmbient
     print*,"sim_gamma = ",sim_gamma
     print*,"sim_gamma_structure = ",sim_gamma_structure
     print*,"sim_Mach = ",sim_Mach
     print*,"sim_epsilon_grad =",sim_epsilon_grad
     print*,"Sound Speed = ", soundspeed
     print*,"Wind Velocity =",vwind
     print*,"Point mass = ", grv_ptmass
     print*,"Accretion radius Racc = 2GM/v_w^2 = ", racc
     print*,"time unit (Racc/vw) = ", racc/vwind
     print*,"Inner BC type = ", sim_bctype
     print*,"sim_pBackground =",sim_pBackground
     print*,"sim_rhoBackground =",sim_rhoBackground
     print*,"##############################################"
     print*,"Binary params:"
     print*,"Primary Mass = ",sim_mprimary
     print*,"mass ratio (param, calc) = ",grv_ptmass/sim_mprimary, sim_q
     print*,"separation = ", sim_a
     print*,"primary y-coord = ",sim_yPrim
     print*,"Hills Radius (accretor) = ",sim_q**(1./3.) * sim_a
     print*,"##############################################"
 endif
 
 ! example of using interpolation
 !call ut_hunt(sim_yprof,sim_nprof,0.25,i)
 !call ut_quadraticInterpol(sim_yprof(i-1:i+1),sim_rhoprof(i-1:i+1),0.25,dp)
 !print*, "interpolation test = ", dp

 


contains
  ! Define a derivs function for the HSE envelope
  function derivs(ycoord, rpvec)
    use Simulation_data, ONLY : sim_mprimary, sim_yPrim, sim_gamma_structure

    real,dimension(2) :: derivs
    real,dimension(2) :: rpvec
    real :: ycoord
    
    ! local vars
    real :: p,rho, dpdr,drhodr, g, Ggrav
    Ggrav = 6.67428e-8
    
    rho = rpvec(1)
    p   = rpvec(2)

    g = Ggrav * sim_mprimary* (ycoord - sim_yPrim)**(-2)
    drhodr = - g*((rho**2) /(sim_gamma_structure * p) ) 
    dpdr   = - g*rho

    derivs = (/ drhodr,dpdr /)
  end function derivs




end subroutine Simulation_init


! Returns a log spaced array in linear space

subroutine logspaced(min,max,npoints,linarray)
  
  implicit none
  
  real, intent (in) :: min
  real, intent (in) :: max
  integer, intent (in) :: npoints
  real, dimension(npoints), intent(out) :: linarray
  
  ! local vars
  real :: step
  integer :: i
  
  linarray(1) = min
  linarray(npoints) = max
  
  step = (LOG(max) - LOG(min)) / (npoints - 1)
  
  do i = 2, npoints - 1
    linarray(i) = LOG(min) + step * (real(i-1))
    linarray(i) = EXP( linarray(i) )
  enddo
  
end subroutine logspaced




