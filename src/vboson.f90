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
  subroutine vboson_resum_coeffs(s, t, u, dlumi_coeff,&
     & wtgg_coeff, wtqg_coeff, wtgq_coeff, wtqq_coeff)
    
    use common_vars
    use input
    implicit none
    
    ! Input variables
    real(dp), intent(in)  :: s, t, u
    real(dp), intent(out)  :: dlumi_coeff(2) 
    
    ! Output variables - weight coefficients
    real(dp), intent(out) :: wtgg_coeff, wtqg_coeff, wtgq_coeff, wtqq_coeff

    real(dp)  :: h11_fac, g_avg
    real(dp)  :: h12qg_coeff, h12gq_coeff, h12qq_coeff, h12gg_coeff
    real(dp)  :: g23qg, g23gq, g23qq, g23gg

    ! Internal variables
    real(dp) :: aleg(3), bleg(3), dleg(3), Eleg(3), Qleg(3)
    real(dp) :: Q2
    real(dp) :: Cqg(3), Cgq(3), Cqq(3)
    real(dp) :: Bqg(3), Bgq(3), Bqq(3)
    real(dp), parameter :: catalan = 0.9159655941772_dp
    real(dp), parameter :: cTT = -(four * catalan + pi * log(two)) / pi

    !--------------------------------------
    ! Set hard scale
    Q2 = M ** 2 
    dlumi_coeff(1) = zero
    dlumi_coeff(2) = zero
    !--------------------------------------
    ! Set observable-dependent parameters
    select case (obs_flag)

    case (obs_onejet)
      aleg(:) = one
      bleg(:) = one
      dleg(1) = sqrt(Q2 * s / t**2)
      dleg(2) = sqrt(Q2 * s / u**2)
      dleg(3) = sqrt(Q2 * s * (t + u)**2 / (four * (u * t)**2))
      h11_fac = one
      g_avg = one

    case (obs_Tminor)
      aleg(:) = one
      bleg(:) = zero
      dleg(:) = sqrt(Q2 * s / (u * t))
      h11_fac = two
      g_avg = half

    case (obs_Tthrust)
      aleg(:) = one
      bleg(1) = zero
      bleg(2) = zero
      bleg(3) = one
      dleg(1) = exp(cTT) * sqrt(Q2 * s / (u * t))
      dleg(2) = exp(cTT) * sqrt(Q2 * s / (u * t))
      dleg(3) = sqrt(Q2 * (u + t)**2/(s * four)) * s/(four * u * t)
      h11_fac = two
      g_avg = 2 / pi

    case default
      call wae_error('vboson_resum_coeffs', 'Unknown observable flag')

    end select

    !--------------------------------------
    Eleg(1) = sqrt(s / Q2)
    Eleg(2) = sqrt(s / Q2)
    Eleg(3) = sqrt((u + t)**2 / (s * Q2))

    Qleg(1) = (Eleg(1) * Q2 * sqrt(s))/ (dleg(1)* sqrt(t * u))
    Qleg(2) = (Eleg(2) * Q2 * sqrt(s))/ (dleg(2)* sqrt(t * u))
    Qleg(3) = (Eleg(3) * Q2 * sqrt(s))/ (dleg(3)* sqrt(t * u))

    !--------------------------------------
    ! Initialize color structure arrays
    Cqg = (/ cf, ca, cf /)
    Cgq = (/ ca, cf, cf /)
    Cqq = (/ cf, cf, ca /)

    Bqg = (/ Bq, Bg, Bq /)
    Bgq = (/ Bg, Bq, Bq /)
    Bqq = (/ Bq, Bq, Bg /)

    !--------------------------------------
    ! For h23
    h12qg_coeff = -two * sum(Cqg / (aleg + bleg))
    h12gq_coeff = -two * sum(Cgq / (aleg + bleg))
    h12qq_coeff = -two * sum(Cqq / (aleg + bleg))
    h12gg_coeff = zero

    g23qg = -twopi_beta0 * (four / three) * sum(Cqg * (two * aleg + bleg) / (aleg + bleg)**2)
    g23gq = -twopi_beta0 * (four / three) * sum(Cgq * (two * aleg + bleg) / (aleg + bleg)**2)
    g23qq = -twopi_beta0 * (four / three) * sum(Cqq * (two * aleg + bleg) / (aleg + bleg)**2)
    g23gg = zero

    !--------------------------------------
    ! Initialize weight coefficients
    wtgg_coeff = zero
    wtqg_coeff = zero
    wtgq_coeff = zero
    wtqq_coeff = zero

    !--------------------------------------
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

    case(resum_h10)
      ! --- qqbar channel ---
      ! Rs
      wtqq_coeff = -half * (two * cf - ca) * (log(dleg(1) * (sqrt(s)/sqrt(Q2)) * (sqrt(s)/(Eleg(1) * sqrt(Q2)))**bleg(1))) ** 2 &
                    - half * (two * cf - ca) * (log(dleg(2) * (sqrt(s)/sqrt(Q2)) * (sqrt(s)/(Eleg(2) * sqrt(Q2)))**bleg(2))) ** 2 & 
                    - half * ca *  (log(dleg(1) * (sqrt(-t)/sqrt(Q2)) * (sqrt(-t)/(Eleg(1) * sqrt(Q2)))**bleg(1))) ** 2 &
                    - half * ca *  (log(dleg(3) * (sqrt(-t)/sqrt(Q2)) * (sqrt(-t)/(Eleg(3) * sqrt(Q2)))**bleg(3))) ** 2 &
                    - half * ca *  (log(dleg(2) * (sqrt(-u)/sqrt(Q2)) * (sqrt(-u)/(Eleg(2) * sqrt(Q2)))**bleg(2))) ** 2 &
                    - half * ca *  (log(dleg(3) * (sqrt(-u)/sqrt(Q2)) * (sqrt(-u)/(Eleg(3) * sqrt(Q2)))**bleg(3))) ** 2
        
  
      ! Chc 
      wtqq_coeff = wtqq_coeff + (67.0_dp/18.0_dp*ca - 13.0_dp/9.0_dp*tr*nf)*bleg(3)/(aleg(3)+bleg(3))
      wtqq_coeff = wtqq_coeff + (11.0_dp/3.0_dp*ca - 4.0_dp/3.0_dp*tr*nf)/(aleg(3)+bleg(3)) * (log(dleg(3)) - bleg(3)*log(Eleg(3)))
      wtqq_coeff = wtqq_coeff + (one/3.0_dp) * tr * nf

      ! dFrec: only for one-jettiness
      wtqq_coeff = wtqq_coeff + ca*(67.0_dp/36.0_dp - pi**2/3.0_dp) - tr*nf*13.0_dp/18.0_dp ! one-jettiness 

      ! DFrec: only for thrust minor
      !wtqq_coeff = wtqq_coeff  + (two/three) * (half * ca - tr * nf) ! thrust minor

      ! ISR
      wtqq_coeff = wtqq_coeff + (three/(aleg(1) + bleg(1))) * cf * (log(dleg(1)) - bleg(1)*log(Eleg(1))) +&
                    + (three/(aleg(2) + bleg(2))) * cf * (log(dleg(2)) - bleg(2)*log(Eleg(2)))


      ! dFwa: only for one-jettiness
      wtqq_coeff = wtqq_coeff + half * (two * cf - ca) * ((log(Qleg(2)/Qleg(1)))**2 + 2.682881831153185_dp)
      wtqq_coeff = wtqq_coeff + half * ca * ((log(Qleg(3)/Qleg(1)))**2 + 0.9935003483842395_dp)
      wtqq_coeff = wtqq_coeff + half * ca * ((log(Qleg(3)/Qleg(2)))**2 + 0.9935003483842395_dp)


      ! --- qg channel ---
      ! Rs
      wtqg_coeff = - half * ca * (log(dleg(1) * (sqrt(s)/sqrt(Q2)) * (sqrt(s)/(Eleg(1) * sqrt(Q2)))**bleg(1))) ** 2 &
                    - half * ca * (log(dleg(2) * (sqrt(s)/sqrt(Q2)) * (sqrt(s)/(Eleg(2) * sqrt(Q2)))**bleg(2))) ** 2 & 
              -half * (two * cf - ca) *  (log(dleg(1) * (sqrt(-t)/sqrt(Q2)) * (sqrt(-t)/(Eleg(1) * sqrt(Q2)))**bleg(1))) ** 2 &
              -half * (two * cf - ca) *  (log(dleg(3) * (sqrt(-t)/sqrt(Q2)) * (sqrt(-t)/(Eleg(3) * sqrt(Q2)))**bleg(3))) ** 2 &
                    - half * ca *  (log(dleg(2) * (sqrt(-u)/sqrt(Q2)) * (sqrt(-u)/(Eleg(2) * sqrt(Q2)))**bleg(2))) ** 2 &
                    - half * ca *  (log(dleg(3) * (sqrt(-u)/sqrt(Q2)) * (sqrt(-u)/(Eleg(3) * sqrt(Q2)))**bleg(3))) ** 2

      ! Chc
      wtqg_coeff = wtqg_coeff + cf * (7.0_dp/2.0_dp * bleg(3)/(aleg(3)+bleg(3)))
      wtqg_coeff = wtqg_coeff + cf * (three/(aleg(3)+bleg(3))*(log(dleg(3)) - bleg(3)*log(Eleg(3))) + half)
   
      ! dFrec: only for one-jettiness
      wtqg_coeff = wtqg_coeff + cf*(5.0_dp/4.0_dp - pi**2/3.0_dp) ! one-jettiness

      ! dFwa: only for one-jettiness
      wtqg_coeff = wtqg_coeff + half * ca * ((log(Qleg(2)/Qleg(1)))**2 + 2.682881831153185_dp)
      wtqg_coeff = wtqg_coeff + half * (two * cf - ca) * ((log(Qleg(3)/Qleg(1)))**2 + 0.9935003483842395_dp)
      wtqg_coeff = wtqg_coeff + half * ca * ((log(Qleg(3)/Qleg(2)))**2 + 0.9935003483842395_dp)

      ! ISR 
      wtqg_coeff = wtqg_coeff + (three/(aleg(1) + bleg(1))) * cf * (log(dleg(1)) - bleg(1)*log(Eleg(1))) + &
                   + (two/(aleg(2) + bleg(2))) * twopi_beta0 * (log(dleg(2)) - bleg(2)*log(Eleg(2))) 
       
      ! --- gq channel ---
      ! Rs
      wtgq_coeff = - half * ca * (log(dleg(1) * (sqrt(s)/sqrt(Q2)) * (sqrt(s)/(Eleg(1) * sqrt(Q2)))**bleg(1))) ** 2 &
                    - half * ca * (log(dleg(2) * (sqrt(s)/sqrt(Q2)) * (sqrt(s)/(Eleg(2) * sqrt(Q2)))**bleg(2))) ** 2 & 
                    - half * ca *  (log(dleg(1) * (sqrt(-t)/sqrt(Q2)) * (sqrt(-t)/(Eleg(1) * sqrt(Q2)))**bleg(1))) ** 2 &
                    - half * ca *  (log(dleg(3) * (sqrt(-t)/sqrt(Q2)) * (sqrt(-t)/(Eleg(3) * sqrt(Q2)))**bleg(3))) ** 2 &
                  - half * (two * cf - ca) *  (log(dleg(2) * (sqrt(-u)/sqrt(Q2)) * (sqrt(-u)/(Eleg(2) * sqrt(Q2)))**bleg(2))) ** 2 &
                  - half * (two * cf - ca) *  (log(dleg(3) * (sqrt(-u)/sqrt(Q2)) * (sqrt(-u)/(Eleg(3) * sqrt(Q2)))**bleg(3))) ** 2      
    

      ! Chc
      wtgq_coeff = wtgq_coeff + cf * (7.0_dp/2.0_dp * bleg(3)/(aleg(3)+bleg(3)))
      wtgq_coeff = wtgq_coeff + cf * (three/(aleg(3)+bleg(3))*(log(dleg(3)) - bleg(3)*log(Eleg(3))) + half)

      ! dFrec: only for one-jettiness
      wtgq_coeff = wtgq_coeff + cf*(5.0_dp/4.0_dp - pi**2/3.0_dp) ! one-jettiness

      ! dFwa: only for one-jettiness
      wtgq_coeff = wtgq_coeff + half * ca * ((log(Qleg(2)/Qleg(1)))**2 + 2.682881831153185_dp)
      wtgq_coeff = wtgq_coeff + half * ca * ((log(Qleg(3)/Qleg(1)))**2 + 0.9935003483842395_dp)
      wtgq_coeff = wtgq_coeff + half * (two * cf - ca)* ((log(Qleg(3)/Qleg(2)))**2 + 0.9935003483842395_dp)

      ! ISR 
      wtgq_coeff = wtgq_coeff + (two/(aleg(1) + bleg(1))) * twopi_beta0 * (log(dleg(1)) - bleg(1)*log(Eleg(1))) + &
                   + (three/(aleg(2) + bleg(2))) * cf * (log(dleg(2)) - bleg(2)*log(Eleg(2)))  
      

      dlumi_coeff(1) = -(log(dleg(1)) - bleg(1)*log(Eleg(1)))
      dlumi_coeff(2) = -(log(dleg(2)) - bleg(2)*log(Eleg(2)))

      

    case (resum_none)
      ! No modification to coefficients

    end select

  end subroutine vboson_resum_coeffs
        

end module vboson


