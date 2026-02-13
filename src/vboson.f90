!========================================================
!--------------------------------------------------------
! Module for calculating the cross-section for vector
! boson production with a jet 
!--------------------------------------------------------
!========================================================

module vboson

  use hoppet_v1
  use ew_parameters

  implicit none

contains

!=======================================================================================
! Cross-section 

  function vboson_cross_section(lumiqq) result(res)

    real(dp), intent(in) :: lumiqq  
    real(dp) :: res

    res = pi / three * sqrt(two) * GF_GeVm2_in * lumiqq * invGev2_to_nb

  end function vboson_cross_section

!=======================================================================================
! Matrix element squared 

  subroutine vboson_Msquared(s, t, u, wtgg, wtqg, wtgq, wtqq)

    real(dp), intent(in) :: s, t, u    
    real(dp), intent(out):: wtqq, wtqg, wtgq, wtgg

    wtqq = cf * ((s + t)**2 + (s + u)**2) / (u * t)
    wtgq = tr * ((s + t)**2 + (u + t)**2) / (-s * u)
    wtqg = tr * ((s + u)**2 + (t + u)**2) / (-s * t)
    wtgg = zero

  end subroutine vboson_Msquared



!=======================================================================================
! Calculate the expansion of the resummation coefficients h11, h12, h23, h24
!=======================================================================================
  subroutine vboson_resum_coeffs(s, t, u, &
     & wtgg_coeff, wtqg_coeff, wtgq_coeff, wtqq_coeff)
    
    use common_vars
    implicit none
    
    ! Input variables
    real(dp), intent(in)  :: s, t, u
    
    ! Output variables - weight coefficients
    real(dp), intent(out) :: wtgg_coeff, wtqg_coeff, wtgq_coeff, wtqq_coeff
    real(dp)  :: h11_fac
    real(dp)  :: h12qg_coeff, h12gq_coeff, h12qq_coeff, h12gg_coeff
    real(dp)  :: g23qg, g23gq, g23qq, g23gg

    ! Internal variables
    real(dp) :: aleg(3), bleg(3), dleg(3), Eleg(3)
    real(dp) :: Q2
    real(dp) :: Cqg(3), Cgq(3), Cqq(3)
    real(dp) :: Bqg(3), Bgq(3), Bqq(3)
    real(dp), parameter :: catalan = 0.9159655941772_dp
    real(dp), parameter :: cTT = -(four * catalan + pi * log(two)) / pi

    ! Initialize scale and energy variables
    Q2 = M ** 2
   

    Eleg(1) = sqrt(s / Q2)
    Eleg(2) = sqrt(s / Q2)
    Eleg(3) = sqrt((u + t)**2 / (s * Q2))

    ! Set observable-dependent parameters
    select case (obs_flag)

    case (obs_onejet)
      aleg(:) = one
      bleg(:) = one
      dleg(1) = sqrt(Q2 * s / t**2)
      dleg(2) = sqrt(Q2 * s / u**2)
      dleg(3) = sqrt(Q2 * s * (t + u)**2 / (four * (u * t)**2))
      h11_fac = one

    case (obs_Tminor)
      aleg(:) = one
      bleg(:) = zero
      dleg(:) = sqrt(Q2 * s / (u * t))
      h11_fac = two

    case (obs_Tthrust)
      aleg(:) = one
      bleg(1) = zero
      bleg(2) = zero
      bleg(3) = one
      dleg(1) = exp(cTT) * sqrt(Q2 * s / (u * t))
      dleg(2) = exp(cTT) * sqrt(Q2 * s / (u * t))
      dleg(3) = sqrt(Q2 * (u + t)**2/(s * four)) * s/(four * u * t)
      h11_fac = two

    case default
      call wae_error('vboson_resum_coeffs', 'Unknown observable flag')

    end select

    ! Initialize color structure arrays
    Cqg = (/ cf, ca, cf /)
    Cgq = (/ ca, cf, cf /)
    Cqq = (/ cf, cf, ca /)


    Bqg = (/ Bq, Bg, Bq /)
    Bgq = (/ Bg, Bq, Bq /)
    Bqq = (/ Bq, Bq, Bg /)

    ! Calculate h12 coefficients
    h12qg_coeff = -two * sum(Cqg / (aleg + bleg))
    h12gq_coeff = -two * sum(Cgq / (aleg + bleg))
    h12qq_coeff = -two * sum(Cqq / (aleg + bleg))
    h12gg_coeff = zero

    ! Calculate g23 coefficients
    g23qg = -twopi_beta0 * (four / three) * sum(Cqg * (two * aleg + bleg) / (aleg + bleg)**2)
    g23gq = -twopi_beta0 * (four / three) * sum(Cgq * (two * aleg + bleg) / (aleg + bleg)**2)
    g23qq = -twopi_beta0 * (four / three) * sum(Cqq * (two * aleg + bleg) / (aleg + bleg)**2)
    g23gg = zero

    ! Initialize weight coefficients
    wtgg_coeff = zero
    wtqg_coeff = zero
    wtgq_coeff = zero
    wtqq_coeff = zero

    ! Calculate weight coefficients based on resummation type
    select case (resum_flag)

    case (resum_h12)
      wtqg_coeff = -sum(Cqg * two / (aleg + bleg))
      wtgq_coeff = -sum(Cgq * two / (aleg + bleg))
      wtqq_coeff = -sum(Cqq * two / (aleg + bleg))

    case (resum_h11, resum_h23)

      ! qqbar channel
      wtqq_coeff = -sum(Cqq * (four * Bqq / (aleg + bleg)))
      wtqq_coeff = wtqq_coeff - sum(Cqq * four / (aleg + bleg) * (log(dleg) - bleg * log(Eleg)))
      wtqq_coeff = wtqq_coeff - two * ca * log(t * u / (s * Q2)) - four * cf * log(s / Q2)

      ! qg channel
      wtqg_coeff = -sum(Cqg * (four * Bqg / (aleg + bleg)))
      wtqg_coeff = wtqg_coeff - sum(Cqg * four / (aleg + bleg) * (log(dleg) - bleg * log(Eleg)))
      wtqg_coeff = wtqg_coeff - two * ca * log(s * u / (t * Q2)) - four * cf * log(-t / Q2)

      ! gq channel
      wtgq_coeff = -sum(Cgq * (four * Bgq / (aleg + bleg)))
      wtgq_coeff = wtgq_coeff - sum(Cgq * four / (aleg + bleg) * (log(dleg) - bleg * log(Eleg)))
      wtgq_coeff = wtgq_coeff - two * ca * log(t * s / (u * Q2)) - four * cf * log(-u / Q2)


    case (resum_h24)
      wtqg_coeff = half * (-sum(Cqg * two / (aleg + bleg)))**2
      wtgq_coeff = half * (-sum(Cgq * two / (aleg + bleg)))**2
      wtqq_coeff = half * (-sum(Cqq * two / (aleg + bleg)))**2
      wtgg_coeff = zero


    case (resum_none)
      ! No modification to coefficients

    end select

  end subroutine vboson_resum_coeffs
        

end module vboson


