!!****if* source/physics/Gravity/GravityMain/CeWindGrav/Gravity_data
!!
!! NAME
!!
!!  Gravity_data
!!
!! SYNOPSIS
!!
!!  use Gravity_data
!!
!! DESCRIPTION
!!
!!  Stores the local data for point-mass gravity.
!!
!! PARAMETERS
!!   grv_ptxpos, grv_ptypos, grv_ptzpos, grv_ptmass, grv_factor
!!
!!***

module Gravity_data

  implicit none

  !! *** Runtime Parameters *** !!

  real, save :: grv_ptxpos, grv_ptypos, grv_ptzpos, grv_ptmass, grv_factor
  

  !! *** Physical Constants *** !!

  real, save :: grv_newton


  integer, save :: grv_meshMe, grv_meshNumProcs
  logical, save :: useGravity
end module Gravity_data
