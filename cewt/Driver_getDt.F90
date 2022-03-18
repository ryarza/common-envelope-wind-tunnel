!!****if* source/Driver/DriverMain/Driver_getDt
!!
!! NAME
!!  Driver_getDt
!!
!! SYNOPSIS
!!  
!!
!!  Driver_getDt(real(out) :: dt)
!!  
!! DESCRIPTION 
!!
!! The current simulation timestep "dt" is known to and owned by the
!! Driver unit. This accessor function makes it available to the
!! other units.
!!
!! ARGUMENTS
!!  dt - returned value, current run simulation dt
!!
!! NOTES
!!
!! NOTES
!!
!! The Driver unit uses a few unit scope variables that are
!! accessible to all routines within the unit, but not to the
!! routines outside the unit. These variables begin with "dr_"
!! like, dr_globalMe or dr_dt, dr_beginStep, and are stored in fortran
!! module Driver_data (in file Driver_data.F90. The other variables
!! are local to the specific routine and do not have the prefix "dr_"
!!
!!
!!***

subroutine Driver_getDt(dt)

  use Driver_data, ONLY : dr_dt

implicit none
  real, intent(out) :: dt
  
  dt = dr_dt
  

end subroutine Driver_getDt

