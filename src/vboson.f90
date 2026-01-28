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
    wtqg = tr * ((s + t)**2 + (u + t)**2) / (-s * u)
    wtgq = tr * ((s + u)**2 + (t + u)**2) / (-s * t)
    wtgg = zero

  end subroutine vboson_Msquared



!=======================================================================================
! Calculate the expansion of the resummation coefficients h11, h12, h23, h24
!=======================================================================================
  subroutine vboson_resum_coeffs( &
      s, t, u, &
      wtgg_coeff, wtqg_coeff, wtgq_coeff, wtqq_coeff, &
      h11_fac, &
      h12qg_coeff, h12gq_coeff, h12qq_coeff, h12gg_coeff, &
      g23qg, g23gq, g23qq, g23gg )
    
    use common_vars
    implicit none
    
    ! Input variables
    real(dp), intent(in)  :: s, t, u
    
    ! Output variables - weight coefficients
    real(dp), intent(out) :: wtgg_coeff, wtqg_coeff, wtgq_coeff, wtqq_coeff
    real(dp), intent(out) :: h11_fac
    real(dp), intent(out) :: h12qg_coeff, h12gq_coeff, h12qq_coeff, h12gg_coeff
    real(dp), intent(out) :: g23qg, g23gq, g23qq, g23gg

    ! Internal variables
    real(dp) :: a(3), b(3), d(3), E(3)
    real(dp) :: Bq, Bg, Q2
    real(dp) :: Cqg(3), Cgq(3), Cqq(3)
    real(dp) :: Bqg(3), Bgq(3), Bqq(3)
    real(dp), parameter :: catalan = 0.9159655941772_dp
    real(dp), parameter :: cTT = -(four * catalan + pi * log(two)) / pi

    ! Initialize scale and energy variables
    Q2 = M ** 2

    E(1) = sqrt(s / Q2)
    E(2) = sqrt(s / Q2)
    E(3) = sqrt((u + t)**2 / (s * Q2))

    ! Initialize beta function coefficients
    Bq = -three / four
    Bg = -(11.0_dp * ca - four * tr * nf) / (12.0_dp * ca)

    ! Set observable-dependent parameters
    select case (obs_flag)

    case (obs_onejet)
      a(:) = one
      b(:) = one
      d(1) = sqrt(Q2 * s / t**2)
      d(2) = sqrt(Q2 * s / u**2)
      d(3) = sqrt(Q2 * s * (t + u)**2 / (four * (u * t)**2))
      h11_fac = one

    case (obs_Tminor)
      a(:) = one
      b(:) = zero
      d(:) = sqrt(Q2 * s / (u * t))
      h11_fac = two

    case (obs_Tthrust)
      a(:) = one
      b(1) = zero
      b(2) = zero
      b(3) = one
      d(1) = exp(cTT) * sqrt(Q2 * s / (u * t))
      d(2) = exp(cTT) * sqrt(Q2 * s / (u * t))
      d(3) = (sqrt(Q2) * s * E(3)) / (4 * t * u)
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
    h12qg_coeff = -two * sum(Cqg / (a + b))
    h12gq_coeff = -two * sum(Cgq / (a + b))
    h12qq_coeff = -two * sum(Cqq / (a + b))
    h12gg_coeff = zero

    ! Calculate g23 coefficients
    g23qg = -twopi_beta0 * (four / three) * sum(Cqg * (two * a + b) / (a + b)**2)
    g23gq = -twopi_beta0 * (four / three) * sum(Cgq * (two * a + b) / (a + b)**2)
    g23qq = -twopi_beta0 * (four / three) * sum(Cqq * (two * a + b) / (a + b)**2)
    g23gg = zero

    ! Initialize weight coefficients
    wtgg_coeff = one
    wtqg_coeff = one
    wtgq_coeff = one
    wtqq_coeff = one

    ! Calculate weight coefficients based on resummation type
    select case (resum_flag)

    case (resum_h12)
      wtqg_coeff = -sum(Cqg * two / (a + b))
      wtgq_coeff = -sum(Cgq * two / (a + b))
      wtqq_coeff = -sum(Cqq * two / (a + b))
      wtgg_coeff = zero

    case (resum_h11, resum_h23)
      wtgg_coeff = zero

      ! qqbar channel
      wtqq_coeff = -sum(Cqq * (four * Bqq / (a + b)))
      wtqq_coeff = wtqq_coeff - sum(Cqq * four / (a + b) * (log(d) - b * log(E)))
      wtqq_coeff = wtqq_coeff - two * ca * log(t * u / (s * Q2)) - four * cf * log(s / Q2)

      ! qg channel
      wtqg_coeff = -sum(Cqg * (four * Bqg / (a + b)))
      wtqg_coeff = wtqg_coeff - sum(Cqg * four / (a + b) * (log(d) - b * log(E)))
      wtqg_coeff = wtqg_coeff - two * ca * log(s * u / (t * Q2)) - four * cf * log(-t / Q2)

      ! gq channel
      wtgq_coeff = -sum(Cgq * (four * Bgq / (a + b)))
      wtgq_coeff = wtgq_coeff - sum(Cgq * four / (a + b) * (log(d) - b * log(E)))
      wtgq_coeff = wtgq_coeff - two * ca * log(t * s / (u * Q2)) - four * cf * log(-u / Q2)

    case (resum_h24)
      wtqg_coeff = half * (-sum(Cqg * two / (a + b)))**2
      wtgq_coeff = half * (-sum(Cgq * two / (a + b)))**2
      wtqq_coeff = half * (-sum(Cqq * two / (a + b)))**2
      wtgg_coeff = zero

    case (resum_none)
      ! No modification to coefficients

    end select

  end subroutine vboson_resum_coeffs
        

end module vboson


