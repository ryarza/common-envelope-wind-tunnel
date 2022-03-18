!!****if* source/Simulation/SimulationMain/FlatPlate/Simulation_data
!!
!! NAME
!!  Simulation_data
!!
!! SYNOPSIS
!!
!!  use Simulation_data
!!
!! DESCRIPTION
!!
!!  Store the simulation data for the wind tunnel problem with a step
!!
!! ARGUMENTS
!!
!!
!! PARAMETERS
!!
!!
!!***

module Simulation_data

  implicit none

  !! *** Runtime Parameters *** !!

  ! wind parameters
  real, save :: sim_pAmbient, sim_rhoAmbient
  real, save :: sim_windVelx, sim_windVely, sim_windVelz, sim_Mach
  real, save :: sim_gamma, sim_smallP, sim_smallX
  real, save :: sim_epsilon_grad, sim_grad_ymax
  integer, parameter :: sim_nprof = 3001  !should be odd!! 
  real, dimension(sim_nprof) :: sim_yprof, sim_rhoprof,sim_pprof, sim_erhoprof
  real, save :: sim_pBackground, sim_rhoBackground 


  ! sink parameters
  real, save :: sim_rhoBulk, sim_pBulk  
  real, save :: sim_radius, sim_xCtr, sim_yCtr, sim_zCtr
  real, save :: sim_presfactor
  integer, save :: sim_bctype  
  logical, save :: sim_avgs
  real, save :: sim_avgs_radius, sim_avgs_mult
  real, save :: sim_tRamp
  real, save :: sim_max_ref_rad

  ! binary parameter
  real, save :: sim_q ! mass ratio
  real, save :: sim_a, sim_mprimary, sim_yPrim
  real, save :: sim_gamma_structure
  real, save :: sim_fkep

  ! integral variables for diagnositcs
  real, save :: sink_mass_old, sink_mass_new, mdot_pres
  real, save,dimension(3) :: sink_j_new, jdot_pres
  real, save :: sink_eint_new,sink_eint_old,eint_dot
  real, save :: sink_ke_new, ke_dot
  real, save :: sim_rhosink_old, sim_psink_old, sim_eintsink_old
  real, save :: sink_xmom_new, xmomdot_pres !AA
  
  ! AA: dynamical friction radii variable
  integer, parameter :: NDFPOINTS = 10
  real, dimension(NDFPOINTS) :: df_radii
  character (len=10), dimension(NDFPOINTS) :: df_radii_labels
  
  ! general
  logical, save :: sim_gCell
  integer, save :: sim_meshMe
  



  
end module Simulation_data


