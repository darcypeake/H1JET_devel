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
!Calculate the expansion of the resummation h11, h12, h23, h24
 subroutine vboson_resum_coeffs( s, t, u, wtgg_coeff, wtqg_coeff, wtgq_coeff, wtqq_coeff, &
      h11_fac, h12qg_coeff, h12gq_coeff, h12qq_coeff, &
      g23qg, g23gq, g23qq )
    use common_vars

    implicit none
    
    real(dp), intent(in)  :: s, t, u
    real(dp), intent(out) :: wtgg_coeff, wtqg_coeff, wtgq_coeff, wtqq_coeff
    real(dp), intent(out) :: h11_fac
    real(dp), intent(out) :: h12qg_coeff, h12gq_coeff, h12qq_coeff
    real(dp), intent(out) :: g23qg, g23gq, g23qq

    real(dp) :: a(3), b(3), d(3), E(3)
    real(dp) :: Cleg(3), Bleg(3)
    real(dp) :: Bq, Bg, Q2, aobs
    real(dp) :: Cqg(3), Cgq(3), Cqq(3)
    real(dp), parameter :: Catalan

    Q2 = M**2

    E(1) = sqrt(s/Q2)
    E(2) = sqrt(s/Q2)
    E(3) = sqrt((u + t)**2/(s*Q2))

    Bq = -three/four
    Bg = -(11.0_dp*ca - four*tr*nf)/(12.0_dp*ca)
    Catalan = 0.915966_dp

    select case (obs_flag)

    case (obs_onejet)
      a(:) = one
      b(:) = one
      d(1) = sqrt(Q2*s/t**2)
      d(2) = sqrt(Q2*s/u**2)
      d(3) = sqrt(Q2*s*(t+u)**2/(four*(u*t)**2))

    case (obs_Tminor)
      a(:) = one
      b(:) = zero
      d(:) = sqrt(Q2*s/(u*t))


    case default
      call wae_error('vboson_resum_coeffs', 'Unknown observable flag')

    end select

    aobs = a(1)

    ! Hadronic collisions
    h11_fac = one/(aobs + b(1)) + one/(aobs + b(2))

    Cqg = (/ cf, ca, cf /)
    Cgq = (/ ca, cf, cf /)
    Cqq = (/ cf, cf, ca /)

    h12qg_coeff = -(two/aobs) * sum( Cqg/(aobs + b) )
    h12gq_coeff = -(two/aobs) * sum( Cgq/(aobs + b) )
    h12qq_coeff = -(two/aobs) * sum( Cqq/(aobs + b) )

    g23qg = -twopi_beta0 * (4.0_dp/(3.0_dp*aobs*aobs)) * &
            sum( Cqg * (two*aobs + b)/(aobs + b)**2 )
    g23gq = -twopi_beta0 * (4.0_dp/(3.0_dp*aobs*aobs)) * &
            sum( Cgq * (two*aobs + b)/(aobs + b)**2 )
    g23qq = -twopi_beta0 * (4.0_dp/(3.0_dp*aobs*aobs)) * &
            sum( Cqq * (two*aobs + b)/(aobs + b)**2 )


    wtgg_coeff = one
    wtqg_coeff = one
    wtgq_coeff = one
    wtqq_coeff = one

    select case (resum_flag)
            
    case (resum_h12)
      wtqg_coeff = -sum( (/cf, ca, cf/) * two/(a+b) )
      wtgq_coeff = -sum( (/ca, cf, cf/) * two/(a+b) )
      wtqq_coeff = -sum( (/cf, cf, ca/) * two/(a+b) )
      wtgg_coeff = zero

    case (resum_h11, resum_h23)
     
      ! qg channel
      Cleg = (/ cf, ca, cf /)
      Bleg = (/ Bq, Bg, Bq /)
      wtqg_coeff = -sum( Cleg * (four * Bleg/(a + b) + four/(a + b) * (log(d)-b * log(E)) ) )
      wtqg_coeff = wtqg_coeff - two*ca*log(s*u/(t*Q2)) - four*cf*log(-t/Q2)

      ! gq channel
      Cleg = (/ ca, cf, cf /)
      Bleg = (/ Bg, Bq, Bq /)
      wtgq_coeff = -sum( Cleg * (four * Bleg/(a + b) + four/(a + b)*(log(d)-b * log(E)) ) )
      wtgq_coeff = wtgq_coeff - two*ca*log(t*s/(u*Q2)) - four*cf*log(-u/Q2)

      ! qqbar channel
      Cleg = (/ cf, cf, ca /)
      Bleg = (/ Bq, Bq, Bg /)
      wtqq_coeff = -sum( Cleg * (four*Bleg/(a+b) + four/(a+b)*(log(d)-b*log(E)) ) )
      wtqq_coeff = wtqq_coeff - two*cf*log(t*u/(s*Q2)) - four*cf*log(s/Q2)

    case (resum_h24)
      wtqg_coeff = half * (-sum( (/cf, ca, cf/) * two/(a+b))) ** 2
      wtgq_coeff = half * (-sum( (/ca, cf, cf/) * two/(a+b))) ** 2
      wtqq_coeff = half * (-sum( (/cf, cf, ca/) * two/(a+b))) ** 2


    case (resum_none)

    end select
    
    
       end subroutine vboson_resum_coeffs
        

end module vboson


