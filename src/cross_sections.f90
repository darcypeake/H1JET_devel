!========================================================
!--------------------------------------------------------
! Module containing general functions to handle the 
! calculation of the cross-section  
!--------------------------------------------------------
!========================================================

module cross_sections

  use hoppet_v1 
  use pdfs_tools
  use common_vars
  use hboson; use vboson; use user_interface 
  
  implicit none

  ! Only used here and in h1jet.f90 
  real(dp), pointer, public :: lumi_gg(:), lumi_qg(:), lumi_gq(:), lumi_qqbar(:)
  ! Hack
  real(dp), pointer, public :: dlumi_gg(:), dlumi_qg(:), dlumi_gq(:), dlumi_qqbar(:)

  public :: h1jet_prefactor, cross_section, dsigma_dptdy  

contains

!=======================================================================================
! Dot product between two four-vectors 

  function dot(p, q) result(res)
    real(dp), intent(in) :: p(4), q(4)
    real(dp) :: res

    res = p(4) * q(4) - dot_product(p(:3), q(:3))

  end function dot

!=======================================================================================
! Process-dependent Electroweak Prefactor 

  subroutine h1jet_prefactor(factor)
    real(dp), intent(out) :: factor

    select case(iproc)
      case (id_H)
        ! alpha_W = g_W^2 / (4 * pi) 
        ! and M_W = (1/2) * g_W * vev 
        factor = mw_in**2 / (pi * higgs_vev_in**2)
      case (id_Z)
        factor = 8.0_dp * pi / three * sqrt(two) * GF_GeVm2_in * mz_in**2  
      case (id_bbH) 
        factor = four * pi / higgs_vev_in**2 / three
      case (id_user)
        ! EW Prefactor already included in user code 
        factor = 1.0_dp 
      case default
        call wae_error('h1jet_prefactor', 'Unrecognised process')
    end select
      
    factor = factor * invGev2_to_fb 

  end subroutine h1jet_prefactor

!=======================================================================================
! Calculate the born-level total cross-section 

  function cross_section(lumi_gg, lumi_qg, lumi_gq, lumi_qqbar, tau) result(res)
    real(dp), intent(in) :: lumi_gg(0:), lumi_qg(0:), lumi_gq(0:), lumi_qqbar(0:)
    type(gdval), intent(in) :: tau

    real(dp) :: res

    select case(iproc)
      case (id_H)
        mh2 = M**2
        res = hboson_cross_section(lumi_gg .atx. tau)
      case (id_Z)
        res = vboson_cross_section(lumi_qqbar .atx. tau)
      case (id_bbH)
        mh2 = M**2
        res = bbH_cross_section(lumi_qqbar .atx. tau)
      case (id_user) 
        res = user_cross_section(lumi_gg .atx. tau, lumi_qg .atx. tau, lumi_gq .atx. tau, lumi_qqbar .atx. tau) 
      case default
        call wae_error('cross_section', 'Unrecognised process')
    end select

  end function cross_section

!=======================================================================================
! Generate particle momenta in the centre-of-mass frame 

  subroutine gen_momenta(y, p)
    real(dp), intent(in) :: y
    real(dp), intent(out) :: p(4,4)

    real(dp) :: E, pz, Eb, Ebeam

    E = pt * cosh(y) 
    pz = pt * sinh(y)
    p(3,:) = (/zero, pt, pz, E/)
    Eb = sqrt(M**2 + E**2)

    p(4,:) = (/zero, -pt, -pz, Eb/)

    Ebeam = E + Eb

    p(1,:) = Ebeam / two * (/zero, zero, one, one/)
    p(2,:) = Ebeam / two * (/zero, zero, -one, one/)
    ! write(*,*) 'Energies', half * Ebeam, half * Ebeam, E
  
  end subroutine gen_momenta

  
!=======================================================================================
! Calculate the differential cross-section dsigma / dpT dy 
! as a function of rapidity y for the specified process 

  function dsigma_dptdy(y) result(res)
    real(dp), intent(in) :: y

    real(dp) :: res
    real(dp) :: p(4,4), rootshat, Eb
    type(gdval) :: tauhat
    real(dp) :: s, t, u
    real(dp) :: lumigg, lumiqg, lumigq, lumiqqbar
    real(dp)  :: dlumigg, dlumiqg, dlumigq, dlumiqqbar
    real(dp) :: wtqq, wtqg, wtgq, wtgg
    real(dp) :: wtqq_coeff, wtqg_coeff, wtgq_coeff, wtgg_coeff
    real(dp) :: jakob
    real(dp) :: h12, h11

    call gen_momenta(y, p)

    rootshat = two * p(1,4)
    Eb = p(4,4)
    jakob = pt / (8.0_dp * pi * Eb * rootshat**3)
    
    tauhat    = (rootshat**2 / roots**2) .with. grid
    lumigg    = lumi_gg .atx. tauhat
    lumiqg    = lumi_qg .atx. tauhat
    lumigq    = lumi_gq .atx. tauhat
    lumiqqbar = lumi_qqbar .atx. tauhat

    ! Hack
    dlumigg    = dlumi_gg .atx. tauhat
    dlumiqg    = dlumi_qg .atx. tauhat
    dlumigq    = dlumi_gq .atx. tauhat
    dlumiqqbar = dlumi_qqbar .atx. tauhat

    s =  two * dot(p(1,:), p(2,:))
    t = -two * dot(p(1,:), p(3,:))
    u = -two * dot(p(2,:), p(3,:))

    !Hack from CAESAR
    ! s = 1385041.5765690478 
    ! t = -710546.52241731714
    ! u = -666179.87575797061

    select case(iproc)
      case (id_H)
        if (cpodd) then
           call hboson_cpodd_Msquared(s, t, u, wtqq, wtqg, wtgq, wtgg)
        else
           call hboson_Msquared(s, t, u, wtqq, wtqg, wtgq, wtgg)
        end if
      case (id_Z)
        call vboson_Msquared(s, t, u, wtgg, wtqg, wtgq, wtqq)
        !Hack
        ! wtgg_coeff = wtgg; wtqg_coeff = wtqg
        ! wtgq_coeff = wtgq; wtqq_coeff = wtqq
       call vboson_resum_coeffs(s, t, u, wtgg_coeff, wtqg_coeff, wtgq_coeff, &
              & wtqq_coeff)
              ! wtgg = wtgg * wtgg_coeff; wtqg = wtqg * wtqg_coeff
              ! wtgq = wtgq * wtgq_coeff; wtqq = wtqq * wtqq_coeff
      case (id_bbH)
        call bbH_Msquared(s, t, u, wtqq, wtqg, wtgq, wtgg)
      case (id_user)
        call user_Msquared(s, t, u, wtqq, wtqg, wtgq, wtgg) 
      case default
        call wae_error('dsigma_dptdy', 'Unrecognised process')
    end select

   ! -----------------------------------------------------------    
  ! h12 and h24 res calculation
  ! res = (wtgg * lumigg + lumiqg * wtqg + lumigq * wtgq + wtqq * lumiqqbar) * jakob 

    select case (resum_flag)
    case (resum_none, resum_h12, resum_h24)
      wtgg = wtgg * wtgg_coeff
      wtqg = wtqg * wtqg_coeff
      wtgq = wtgq * wtgq_coeff
      wtqq = wtqq * wtqq_coeff

      res = (wtgg * lumigg + lumiqg * wtqg + &
            lumigq * wtgq + wtqq * lumiqqbar) * jakob

    case (resum_h11)

      res = wtgg * wtgg_coeff * lumigg + &
            wtqg * wtqg_coeff * lumiqg + &
            wtgq * wtgq_coeff * lumigq + &
            wtqq * wtqq_coeff * lumiqqbar

      res = res - wtgg * dlumigg - &
                  wtqg * dlumiqg - &
                  wtgq * dlumigq - &
                  wtqq * dlumiqqbar

      res = res * jakob

    case (resum_h23)

      res = wtgg * wtgg_coeff * lumigg + &
            wtqg * wtqg_coeff * lumiqg + &
            wtgq * wtgq_coeff * lumigq + &
            wtqq * wtqq_coeff * lumiqqbar

      res = res - wtgg * dlumigg - &
                  wtqg * dlumiqg - &
                  wtgq * dlumigq - &
                  wtqq * dlumiqqbar

      res = -(two * cf + ca) * res

      res = res - twopi_beta0 * (two * cf + ca) * (wtgg * lumigg + &
            & lumiqg * wtqg + lumigq * wtgq + &
            & wtqq * lumiqqbar)

      res = res * jakob



    end select

! -----------------------------------------------------------   
  ! h11 res calculation
    ! res = wtgg_coeff * lumigg + lumiqg * wtqg_coeff + lumigq * wtgq_coeff +&
    !        & wtqq_coeff * lumiqqbar
    ! res = res-wtgg * dlumigg - dlumiqg * wtqg - dlumigq * wtgq - wtqq * dlumiqqbar
    ! res = res * jakob

  ! Code to test the individual channels
    !  res = wtqq_coeff * lumiqqbar * jakob
    !  res = wtqg_coeff * lumiqg * jakob
!     res = wtgq_coeff * lumigq * jakob

  ! Code to test the luminosity of the individual channels
    ! res = - wtqq * dlumiqqbar * jakob
    ! res = - dlumiqg * wtqg * jakob
    ! res = - dlumigq * wtgq * jakob
   
    
! -----------------------------------------------------------    
  ! h23 res calculation
    ! res = wtgg_coeff * lumigg + lumiqg * wtqg_coeff + lumigq * wtgq_coeff +&
    !       & wtqq_coeff * lumiqqbar
    ! res = res-wtgg * dlumigg - dlumiqg * wtqg - dlumigq * wtgq - wtqq * dlumiqqbar
    ! res = -(two * cf + ca) * res
    ! res = res - twopi_beta0 * (two * cf + ca) * (wtgg * lumigg + lumiqg * wtqg + lumigq * wtgq +&
    !       & wtqq * lumiqqbar)
    ! res = res * jakob

  end function dsigma_dptdy


!=======================================================================================

end module cross_sections 


