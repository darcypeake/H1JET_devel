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
  subroutine vboson_resum_coeffs(s, t, u, wtgg_coeff, wtqg_coeff, wtgq_coeff, wtqq_coeff)

         use common_vars

         implicit none
    
         real(dp), intent(in) :: s, t, u
         real(dp), intent(inout) :: wtgg_coeff, wtqg_coeff, wtgq_coeff, wtqq_coeff
         real(dp) :: h12, h11, h24
    
         real(dp) :: h11gg, h11qg, h11gq, h11qqbar
         real(dp) :: y
         real(dp) :: E, Eb, Ebeam
         real(dp) :: p(4,4)
         real(dp) :: dl, d1, d2, d3
         real(dp) :: Q2, E1, E2, E3


         ! Define variables
          Q2 = M ** 2
          d1 = sqrt(Q2 * s/t**2)
          d2 = sqrt(Q2 * s/u**2)
          d3 = sqrt((Q2 * s * (t + u)**2)/(four * (u * t)**2))
          E1 = sqrt(s/Q2)
          E2 = sqrt(s/Q2)
          E3 = sqrt((u + t)**2/s/Q2)


         ! Initialise coefficients
          wtgg_coeff = one
          wtqg_coeff = one
          wtgq_coeff = one
          wtqq_coeff = one


          select case (resum_flag)
            
            case (resum_h12)
              wtqq_coeff = -(two * cf + ca)
              wtqg_coeff = -(two * cf + ca)
              wtgq_coeff = -(two * cf + ca)

            case (resum_h11, resum_h23)
              wtgg_coeff = zero
              wtqg_coeff = zero
              wtgq_coeff = zero
              wtqq_coeff = zero

              wtqg_coeff = wtqg_coeff + three * cf + (11.0_dp * ca - four * tr * nf)/6.0_dp
              wtqg_coeff = wtqg_coeff - two * cf * (log(d1) - log(E1))
              wtqg_coeff = wtqg_coeff - two * ca * (log(d2) - log(E2))
              wtqg_coeff = wtqg_coeff - two * cf * (log(d3) - log(E3))
              wtqg_coeff = wtqg_coeff - two * ca * log(s * u/(t * Q2)) - four * cf * log(-t/Q2)

              wtgq_coeff = wtgq_coeff + three * cf + (11.0_dp * ca - four * tr * nf)/6.0_dp
              wtgq_coeff = wtgq_coeff - two * ca * (log(d1) - log(E1))
              wtgq_coeff = wtgq_coeff - two * cf * (log(d2) - log(E2))
              wtgq_coeff = wtgq_coeff - two * cf * (log(d3) - log(E3))
              wtgq_coeff = wtgq_coeff - two * ca * log(t * s/(u * Q2)) - four * cf * log(-u/Q2)

              wtqq_coeff = wtqq_coeff + three * cf + (11.0_dp * ca - four * tr * nf)/6.0_dp 
              wtqq_coeff = wtqq_coeff - two * cf * (log(d1) - log(E1))
              wtqq_coeff = wtqq_coeff - two * cf * (log(d2) - log(E2))
              wtqq_coeff = wtqq_coeff - two * ca * (log(d3) - log(E3))
              wtqq_coeff = wtqq_coeff - two * cf * log(t * u/(s * Q2)) - four * cf * log(s/Q2)

            case (resum_h24)
              wtqq_coeff = half * (-(two * cf + ca))**2
              wtqg_coeff = half * (-(two * cf + ca))**2
              wtgq_coeff = half * (-(two * cf + ca))**2


            case (resum_none)

          end select
    
    ! -----------------------------------------------------------
        ! ! h12  
        ! wtqq_coeff = -(two * cf + ca)
        ! wtqg_coeff = -(two * cf + ca)
        ! wtgq_coeff = -(two * cf + ca)

    
    ! -----------------------------------------------------------    
       ! h24
      !  wtqq_coeff = half * (-(two * cf + ca)) ** 2
      !  wtqg_coeff = half * (-(two * cf + ca)) ** 2
      !  wtgq_coeff = half * (-(two * cf + ca)) ** 2
    
    ! -----------------------------------------------------------    
        ! h11 (and h23)
        ! ! Initialise h11
        !  h11 = zero
    
        ! !Initialise the individual channels
        !  h11gg = zero
        !  h11qg = zero
        !  h11gq = zero
        !  h11qqbar = zero
  
    
        !Hard scales
        ! Q2 = M ** 2
        !  Q2 = M ** 2 + (u * t)/s
        !  Q2 = (sqrt(M ** 2 + (u * t)/s) + sqrt(u * t/s))**2
        !  Q2 = u * t/s
    
         !Thrust minor d_ell
          ! dl = two * sqrt(Q2 * s/(u * t))

         ! One-jettiness d_ell
          ! d1 = sqrt(Q2 * s/t**2)
          ! d2 = sqrt(Q2 * s/u**2)
          ! d3 = sqrt((Q2 * s * (t + u)**2)/(four * (u * t)**2))

          ! E1 = sqrt(s/Q2)
          ! E2 = sqrt(s/Q2)
          ! E3 = sqrt((u + t)**2/s/Q2)
        
    
        ! qg - channel 

        ! h11qg = h11qg + three * cf + (11.0_dp * ca - four * tr * nf)/6.0_dp
        ! h11qg = h11qg - two * cf * (log(d1) - log(E1))
        ! h11qg = h11qg - two * ca * (log(d2) - log(E2))
        ! h11qg = h11qg - two * cf * (log(d3) - log(E3))
        ! h11qg = h11qg - two * ca * log(s * u/(t * Q2)) - four * cf * log(-t/Q2)
        ! write(*,*) -two * cf * (log(d1) - log(E1)), -two * ca * (log(d2) - log(E2)), -two * cf * (log(d3) - log(E3))
        
        ! gq - channel

        ! h11gq = h11gq + three * cf + (11.0_dp * ca - four * tr * nf)/6.0_dp
        ! h11gq = h11gq - two * ca * (log(d1) - log(E1))
        ! h11gq = h11gq - two * cf * (log(d2) - log(E2))
        ! h11gq = h11gq - two * cf * (log(d3) - log(E3))
        ! h11gq = h11gq - two * ca * log(t * s/(u * Q2)) - four * cf * log(-u/Q2)
        ! write(*,*) - two * ca * (log(d1) - log(E1)), - two * cf * (log(d2) - log(E1)), - two * cf * (log(d3) - log(E3))
    
        ! qqbar - channel

        ! h11qqbar = h11qqbar + three * cf + (11.0_dp * ca - four * tr * nf)/6.0_dp 
        ! h11qqbar = h11qqbar - two * cf * (log(d1) - log(E1))
        ! h11qqbar = h11qqbar - two * cf * (log(d2) - log(E2))
        ! h11qqbar = h11qqbar - two * ca * (log(d3) - log(E3))
        ! h11qqbar = h11qqbar - two * cf * log(t * u/(s * Q2)) - four * cf * log(s/Q2)
        ! write(*,*) -two * cf * (log(d1) - log(E1)), -two * cf * (log(d2) - log(E2)), -two * ca * (log(d3) - log(E3))
        
        ! h11qqbar = h11qqbar - two * ca * log(u * t/(s * Q2)) - four * cf * log(s/Q2)
        
        ! Thrust minor
        ! h11qg = h11qg - (two * cf + ca) *  four * (log(dl)-log(two))
        !  h11gq = h11gq - (two * cf + ca) *  four * (log(dl)-log(two))
        ! h11qqbar = h11qqbar - (two * cf + ca) *  four * (log(dl)-log(two))

      

      !  wtqg = wtqg * h11qg
      !  wtgq = wtgq * h11gq
      !  wtqq = wtqq * h11qqbar
    
    
       end subroutine vboson_resum_coeffs
        

end module vboson


