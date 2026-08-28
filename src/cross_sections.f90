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
  real(dp), pointer, public :: dlumi_gg(:,:), dlumi_qg(:,:), dlumi_gq(:,:), dlumi_qqbar(:,:)
  ! Hack
  real(dp), pointer, public :: zlumi_gg(:), zlumi_qg(:), zlumi_gq(:), zlumi_qqbar(:)
  real(dp), pointer, public :: DFrec_lumi_gg(:), DFrec_lumi_qg(:), DFrec_lumi_gq(:), DFrec_lumi_qqbar(:)

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
  use input
    real(dp), intent(in) :: lumi_gg(0:), lumi_qg(0:), lumi_gq(0:), lumi_qqbar(0:)
    type(gdval), intent(in) :: tau
    real(dp) :: jetveto_coeff(2)
    real(dp) :: res

    jetveto_coeff = log((M/2)**2/muF**2)
    ! jetveto_coeff = one
    
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
    real(dp)  :: dlumigg(2), dlumiqg(2), dlumigq(2), dlumiqqbar(2)
    real(dp)  :: zlumigg, zlumiqg, zlumigq, zlumiqqbar
    real(dp)  :: DFreclumigg, DFreclumiqg, DFreclumigq, DFreclumiqqbar
    real(dp) :: wtqq, wtqg, wtgq, wtgg
    real(dp) :: wtqq_coeff, wtqg_coeff, wtgq_coeff, wtgg_coeff
    real(dp) :: jakob
    real(dp) :: h12, h11
    real(dp) :: h11_fac
    real(dp) :: h12qg_coeff, h12gq_coeff, h12qq_coeff, h12gg_coeff
    real(dp) :: g23qg, g23gq, g23qq, g23gg
    real(dp) :: muF, muR, dlumi_coeff(2)
    real(dp) :: PL_gq, PL_qg

    call gen_momenta(y, p)

    rootshat = two * p(1,4)
    Eb = p(4,4)
    jakob = pt / (8.0_dp * pi * Eb * rootshat**3)

    tauhat = (rootshat**2 / roots**2) .with. grid

    lumigg    = lumi_gg .atx. tauhat
    lumiqg    = lumi_qg .atx. tauhat
    lumigq    = lumi_gq .atx. tauhat
    lumiqqbar = lumi_qqbar .atx. tauhat

    dlumigg(1)    = dlumi_gg(1, :) .atx. tauhat
    dlumigg(2)    = dlumi_gg(2, :) .atx. tauhat
    dlumiqg(1)    = dlumi_qg(1, :) .atx. tauhat
    dlumiqg(2)    = dlumi_qg(2, :) .atx. tauhat
    dlumigq(1)    = dlumi_gq(1, :) .atx. tauhat
    dlumigq(2)    = dlumi_gq(2, :) .atx. tauhat
    dlumiqqbar(1) = dlumi_qqbar(1, :) .atx. tauhat
    dlumiqqbar(2) = dlumi_qqbar(2, :) .atx. tauhat

    ! Hack
    zlumigg    = zlumi_gg .atx. tauhat
    zlumiqg    = zlumi_qg .atx. tauhat
    zlumigq    = zlumi_gq .atx. tauhat
    zlumiqqbar = zlumi_qqbar .atx. tauhat
  

    ! Hack
    DFreclumigg    = DFrec_lumi_gg .atx. tauhat
    DFreclumiqg    = DFrec_lumi_qg .atx. tauhat
    DFreclumigq    = DFrec_lumi_gq .atx. tauhat
    DFreclumiqqbar = DFrec_lumi_qqbar .atx. tauhat

    s =  two * dot(p(1,:), p(2,:))
    t = -two * dot(p(1,:), p(3,:))
    u = -two * dot(p(2,:), p(3,:))   

    ! DFrec pre-factors
    PL_gq = (two * s * u) / (s**2 + u**2 + two * M**2 * t) 
    PL_qg = (two * s * t) / (s**2 + t**2 + two * M**2 * u)

    select case(iproc)
      case (id_H)
        if (cpodd) then
           call hboson_cpodd_Msquared(s, t, u, wtqq, wtqg, wtgq, wtgg)
        else
           call hboson_Msquared(s, t, u, wtqq, wtqg, wtgq, wtgg)
        end if
      case (id_Z)
        call vboson_Msquared(s, t, u, wtgg, wtqg, wtgq, wtqq)
        ! Hack
        call vboson_resum_coeffs( s, t, u, dlumi_coeff, &
          wtgg_coeff, wtqg_coeff, wtgq_coeff, wtqq_coeff)! &
!          h11_fac, h12qg_coeff, h12gq_coeff, h12qq_coeff, h12gg_coeff, &
!          g23qg, g23gq, g23qq, g23gg )
      case (id_bbH)
        call bbH_Msquared(s, t, u, wtqq, wtqg, wtgq, wtgg)
      case (id_user)
        call user_Msquared(s, t, u, wtqq, wtqg, wtgq, wtgg) 
      case default
        call wae_error('dsigma_dptdy', 'Unrecognised process')
    end select

      ! res = (wtgg * lumigg + lumiqg * wtqg + &
      !       lumigq * wtgq + wtqq * lumiqqbar) * jakob

    !=======================================================================================
    ! Resummation coefficients for V boson production
    if (resum_flag > 0) then
       res =  wtgg * wtgg_coeff * lumigg + &
            & wtqg * wtqg_coeff * lumiqg + &
            & wtgq * wtgq_coeff * lumigq + &
            & wtqq * wtqq_coeff * lumiqqbar
    else
       res =  wtgg * lumigg +  lumiqg * wtqg +  lumigq * wtgq +  wtqq * lumiqqbar
      
    end if


    select case (resum_flag)
    case (resum_h11,resum_h23)
      ! one-jettiness  
      res = res - (wtgg * sum(dlumigg) + wtqg * sum(dlumiqg) + wtgq * sum(dlumigq) + wtqq * sum(dlumiqqbar))
      
      ! thrust minor
      !res = res - two * (wtgg * sum(dlumigg) + wtqg * sum(dlumiqg) + wtgq * sum(dlumigq) + wtqq * sum(dlumiqqbar))
    
    case (resum_h10)

      ! one-jettiness  
      res = res  + (wtgg * zlumigg + wtqg * zlumiqg + wtgq * zlumigq + wtqq * zlumiqqbar)
      res = res + (wtgg * dot_product(dlumigg,dlumi_coeff) + &
        & wtqg * dot_product(dlumiqg,dlumi_coeff) + &
        & wtgq * dot_product(dlumigq,dlumi_coeff) + &
        & wtqq * dot_product(dlumiqqbar,dlumi_coeff))
     
      ! thrust minor
      !res = res  + (wtgg * zlumigg + wtqg * zlumiqg + wtgq * zlumigq + wtqq * zlumiqqbar)
      !res = res + two * (wtgg * dot_product(dlumigg,dlumi_coeff) + &
      !  & wtqg * dot_product(dlumiqg,dlumi_coeff) + &
      !  & wtgq * dot_product(dlumigq,dlumi_coeff) + &
      !  & wtqq * dot_product(dlumiqqbar,dlumi_coeff))
      !res = res + (PL_qg * wtqg * DFreclumiqg + PL_gq * wtgq * DFreclumigq)

    
    end select

    res = res * jakob

  end function dsigma_dptdy


!=======================================================================================

end module cross_sections 


