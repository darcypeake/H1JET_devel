!========================================================
!--------------------------------------------------------
! Module containing commonly used public variables  
!--------------------------------------------------------
!========================================================

module common_vars 

  use hoppet_v1

  implicit none

  ! Output level 
  integer, public, parameter :: stdout = 6
  integer, public, parameter :: stderr = 0

  ! Maximum number of warnings
  integer :: max_warns = 3

  ! Process id 
  integer, public, parameter :: id_H = 1
  integer, public, parameter :: id_Z = 2
  integer, public, parameter :: id_bbH = 3
  integer, public, parameter :: id_user = 4
  integer, public :: iproc

  ! User interface related id's 
  integer, public, parameter :: id_noImplementation = 0
  integer, public, parameter :: id_missingTotXsec = 1 
  integer, public, parameter :: id_fullImplementation = 2

  ! Hack: Resummation coefficients control flag
  integer, public, parameter :: resum_none = 0
  integer, public, parameter :: resum_h12 = 1
  integer, public, parameter :: resum_h11 = 2
  integer, public, parameter :: resum_h24 = 3
  integer, public, parameter :: resum_h23 = 4
  integer, public :: resum_flag

  ! Hack: Observable flags
  integer, public, parameter :: obs_onejet   = 1
  integer, public, parameter :: obs_Tminor = 2
  integer, public :: obs_flag 



  ! Kinematics related variables 
  real(dp), public :: M
  real(dp), public :: roots, pt 

  ! Strong coupling alpha_s 
  real(dp), public :: alphas 

end module common_vars 



