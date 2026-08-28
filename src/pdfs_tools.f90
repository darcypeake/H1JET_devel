!========================================================
!--------------------------------------------------------
! Includes a module providing the parton luminosities
! and their convolutions with splitting and coefficient
! functions @ LO, NLO, NNLO
!--------------------------------------------------------
!========================================================

module pdfs_tools
  ! If using LHAPDF, rename a couple of hoppet functions which
  ! would otherwise conflict with LHAPDF 
  use hoppet_v1, EvolvePDF_hoppet => EvolvePDF, InitPDF_hoppet => InitPDF

  implicit none

  type(dglap_holder), save, public :: dglap_h ! Splitting function holder
  type(pdf_table) :: PDFs ! Parton densities
  type(grid_def), public :: grid
  type(split_mat), save, public :: ISR_matrix !Hack
  type(split_mat), save, public :: DFrec_matrix !Hack

  private 

  public :: init_pdfs_from_LHAPDF
  public :: init_pdfs_hoppet_evolution
  public :: get_pdfs, lumi_at_x, lumi_all_x

  public :: luminosities, luminosities_expanded
  public :: lumi_Z_coeffs !Hack 
  public :: lumi_DFrec_coeffs !Hack 

contains

!=======================================================================================
! Set up Hoppet's x grid, splitting functions and PDF tabulation object 
! (just the memory allocation -- it gets filled below) 

  subroutine init_grid_and_dglap(Qmax) 
    real(dp), intent(in), optional :: Qmax
    type(grid_def) ::  gdarray(4) ! grid 
    integer :: pdf_fit_order, nf_lcl
    real(dp) :: dy, ymax, dlnlnQ

    ! Build the PDF grid 
    dy = 0.10_dp ! Grid spacing, 0.1 is a sensible value 
    dlnlnQ = dy / 4.0_dp ! Table spacing in lnlnQ 
    ymax = 15.0_dp ! Highest value of ln(1/x) we want to access 
    pdf_fit_order = -5 ! Order of numerical interpolation 
    
    ! Create nested grid 
    !call InitGridDef(gdarray(5), dy / 81.0_dp, 0.1_dp, order = pdf_fit_order)
    call InitGridDef(gdarray(4), dy / 27.0_dp, 0.2_dp, order = pdf_fit_order)
    call InitGridDef(gdarray(3), dy / 9.0_dp,  0.5_dp, order = pdf_fit_order)
    call InitGridDef(gdarray(2), dy / 3.0_dp,  2.0_dp, order = pdf_fit_order)
    call InitGridDef(gdarray(1), dy,           ymax  , order = pdf_fit_order)
    call InitGridDef(grid, gdarray(1:4), locked = .true.)

    nf_lcl = 5
    call qcd_SetNf(nf_lcl)
    call InitDglapHolder(grid, dglap_h,  factscheme = factscheme_MSbar, nloop = 3)

    if (present(Qmax)) then
      call AllocPdfTable(grid, PDFs, 1d0, Qmax, dlnlnQ = dlnlnQ, freeze_at_Qmin = .true.)
    else
      call AllocPdfTable(grid, PDFs, 1d0, 2d4, dlnlnQ = dlnlnQ, freeze_at_Qmin = .true.)
    end if

  end subroutine init_grid_and_dglap

!=======================================================================================
! Initialise Hoppet's PDF table from LHAPDF
   
  subroutine init_pdfs_from_LHAPDF(pdf_name, pdf_set)
    use ew_parameters ! For mz and mz_in 

    character(len=*), intent(in) :: pdf_name
    integer,          intent(in) :: pdf_set

    real(dp) :: alphasMZ, alphasPDF
    !real(dp) :: mz = 91.1876_dp

    interface
      subroutine evolvePDF(x,Q,res)
        use types; implicit none
        real(dp), intent(in)  :: x,Q
        real(dp), intent(out) :: res(*)
      end subroutine evolvePDF
    end interface

    call init_grid_and_dglap()

    ! Set up LHAPDF
    call InitPDFsetByName(trim(pdf_name))
    call InitPDF(pdf_set)

    ! Sort out the coupling: get it from LHAPDF, which will have been initialised in init_pdfs
    alphasMZ = alphasPDF(mz_in)
    call InitRunningCoupling(alfas = alphasMZ, Q = mz_in, nloop = 3, fixnf = nf_int)

    call FillPdfTable_LHAPDF(PDFs, evolvePDF)

  end subroutine init_pdfs_from_LHAPDF

!=======================================================================================
! Initialise Hoppet's PDF table by evolution from an LHAPDF input
! at some reference scale.

  subroutine init_pdfs_hoppet_evolution(pdf_name, pdf_set, ref_Q, alphas_Q, rts)
    character(len=*), intent(in) :: pdf_name
    integer,          intent(in) :: pdf_set
    real(dp),         intent(in) :: ref_Q, alphas_Q, rts

    real(dp), pointer :: pdf_Q(:,:)
    real(dp) :: Qmax
    interface
       subroutine evolvePDF(x,Q,res)
         use types; implicit none
         real(dp), intent(in)  :: x,Q
         real(dp), intent(out) :: res(*)
       end subroutine evolvePDF
    end interface

    Qmax = max(2e4_dp, rts)
    call init_grid_and_dglap(rts)
    call AllocPDF(grid, pdf_Q)

    ! Set up LHAPDF
    if (pdf_name == "dummy") then
      pdf_Q = unpolarized_dummy_pdf(xValues(grid))
    else
      call InitPDFsetByName(trim(pdf_name))
      call InitPDF(pdf_set)
      ! Extract the PDF at our reference Q value
      call InitPDF_LHAPDF(grid, pdf_Q, evolvePDF, ref_Q)
    end if

    ! Set up a running coupling, set up the PDF tabulation and perform the evolution
    call InitRunningCoupling(alfas = alphas_Q, Q = ref_Q, nloop = 3, fixnf = nf_int)
    call EvolvePdfTable(PDFs, ref_Q, pdf_Q, dglap_h, ash_global)

  end subroutine init_pdfs_hoppet_evolution

!=======================================================================================

  function lumi_all_x(pdf1, pdf2) result(res)
    real(dp), intent(in) :: pdf1(:,-6:), pdf2(:,-6:)
    real(dp) :: res(size(pdf1, 1))

    integer :: i

    res = PartonLuminosity(grid, pdf1(:, iflv_g), pdf2(:, iflv_g))

  end function lumi_all_x

!=======================================================================================
! This works out the partonic luminosity for the LO of process proc, 
! returns at a fixed value x = M^2/rts^2; gives a luminosity
! based on the two input PDFs

  function lumi_at_x(pdf1, pdf2, rtshat, rts) result(res)
    real(dp), intent(in) :: pdf1(:,-6:), pdf2(:,-6:)
    real(dp), intent(in) :: rtshat, rts
    real(dp) :: res, shat_rts2

    shat_rts2 = rtshat**2 / rts**2
    res = lumi_all_x(pdf1, pdf2) .atx. (shat_rts2 .with. grid)

  end function lumi_at_x

!=======================================================================================
!Hack - calculate the function of the expanded luminosity

  function lumi_at_x_expanded(pdf1, pdf2) result(res)
    real(dp), intent(in) :: pdf1(0:grid%ny,-6:7), pdf2(0:grid%ny,-6:7)
    real(dp) :: P_x_pdf1(0:grid%ny,-6:7), P_x_pdf2(0:grid%ny,-6:7)
    real(dp) :: lumi_gg(0:grid%ny), lumi_qg(0:grid%ny), lumi_qqbar(0:grid%ny)
    real(dp) :: res
    type(gdval) :: x

    P_x_pdf1 = dglap_h%P_LO * pdf1
    P_x_pdf2 = dglap_h%P_LO * pdf2

    res = res + lumi_all_x(P_x_pdf1, pdf2).atx.x
    res = res + lumi_all_x(pdf1, P_x_pdf2).atx.x

  end function lumi_at_x_expanded

!=======================================================================================
! Return the two incoming PDFs evaluated at scale muF, taking into
! account whether the collider is pp or ppbar

  subroutine get_pdfs(muF, collider, pdf1, pdf2)
    real(dp),         intent(in)  :: muF
    character(len=*), intent(in)  :: collider
    real(dp),         intent(out) :: pdf1(:,-6:), pdf2(:,-6:)

    call EvalPdfTable_Q(PDFs, muF, pdf1)
    select case(trim(collider))
      case("pp")
        pdf2 = pdf1
      case("ppbar")
        pdf2(:,-6:6) = pdf1(:,6:-6:-1)
        pdf2(:,7)    = pdf1(:,7)          ! index 7 contains info on representation in flavour space
      case default
        call wae_error("get_pdfs: unrecognized collider "//collider)
    end select

  end subroutine get_pdfs

!=======================================================================================
! The dummy PDF suggested by Vogt as the initial condition for the 
! unpolarized evolution (as used in hep-ph/0511119).

  function unpolarized_dummy_pdf(xvals) result(pdf)
    real(dp), intent(in) :: xvals(:)
    real(dp)             :: pdf(size(xvals), ncompmin:ncompmax)
    real(dp) :: uv(size(xvals)), dv(size(xvals))
    real(dp) :: ubar(size(xvals)), dbar(size(xvals))

    real(dp), parameter :: N_g = 1.7_dp, N_ls = 0.387975_dp
    real(dp), parameter :: N_uv = 5.107200_dp, N_dv = 3.064320_dp
    real(dp), parameter :: N_db = half * N_ls

    pdf = zero
    ! Clean method for labelling as PDF as being in the human representation
    ! (not actually needed after setting pdf = 0) 
    call LabelPdfAsHuman(pdf)

    ! Remember that these are all xvals*q(xvals)
    uv = N_uv * xvals**0.8_dp * (1 - xvals)**3
    dv = N_dv * xvals**0.8_dp * (1 - xvals)**4
    dbar = N_db * xvals**(-0.1_dp) * (1 - xvals)**6
    ubar = dbar * (1 - xvals)

    ! The labels iflv_g, etc., comes from the hoppet_v1 module, 
    ! inherited from the main program
    pdf(:,  iflv_g) = N_g * xvals**(-0.1_dp) * (1 - xvals)**5
    pdf(:, -iflv_s) = 0.2_dp * (dbar + ubar)
    pdf(:,  iflv_s) = pdf(:, -iflv_s)
    pdf(:,  iflv_u) = uv + ubar
    pdf(:, -iflv_u) = ubar
    pdf(:,  iflv_d) = dv + dbar
    pdf(:, -iflv_d) = dbar

  end function unpolarized_dummy_pdf

!=======================================================================================
! Hack for pdf conv splitting functions
! Evaluated at shat/s
! Used for h11, h23, h10

  subroutine luminosities_expanded(muF, collider, lumi_gg, lumi_qg, lumi_gq, lumi_qqbar)
    use ew_parameters, only : gv2_ga2
    !use hboson
    use user_interface, only : user_luminosities
    use common_vars

    real(dp), intent(in) :: muF
    character(len=*), intent(in) :: collider
    !integer, intent(in) :: iproc

    real(dp), intent(out) :: lumi_gg(2,0:grid%ny), lumi_qg(2, 0:grid%ny), lumi_gq(2, 0:grid%ny), lumi_qqbar(2, 0:grid%ny)

    real(dp) :: pdf1(0:grid%ny,-6:7), pdf2(0:grid%ny,-6:7)
    integer  :: i

    !Hack
    real(dp) :: P_x_pdf1(0:grid%ny,-6:7), P_x_pdf2(0:grid%ny,-6:7)
    real(dp) :: res
    type(gdval) :: x


    call get_pdfs(muF, collider, pdf1, pdf2)

    P_x_pdf1 = dglap_h%P_LO .conv. pdf1
    P_x_pdf2 = dglap_h%P_LO .conv. pdf2

    select case(iproc)
      case(id_H)
        lumi_gg(1,:) = PartonLuminosity(grid, P_x_pdf1(:,iflv_g), pdf2(:,iflv_g)) 
        lumi_gg(2,:) = PartonLuminosity(grid, pdf1(:,iflv_g), P_x_pdf2(:,iflv_g))
        lumi_gq(1,:) = PartonLuminosity(grid, P_x_pdf1(:,iflv_g), &
              & sum(pdf2(:,-6:-1), dim = 2) + sum(pdf2(:, 1:6 ), dim = 2) )
        lumi_gq(2,:) = PartonLuminosity(grid, pdf1(:,iflv_g), &
              & sum(P_x_pdf2(:,-6:-1), dim = 2) + sum(P_x_pdf2(:, 1:6 ), dim = 2) )
        lumi_qg(1,:) = PartonLuminosity(grid, &
              & sum(P_x_pdf1(:,-6:-1),dim = 2) + sum(P_x_pdf1(:, 1:6 ), dim = 2),pdf2(:,iflv_g)) 
        lumi_qg(2,:) = PartonLuminosity(grid, &
              & sum(pdf1(:,-6:-1),dim = 2) + sum(pdf1(:, 1:6 ), dim = 2),P_x_pdf2(:,iflv_g))
        lumi_qqbar(1,:) = 0 
        lumi_qqbar(2,:) = 0
        do i = -6, 6
          if (i == 0) cycle
          lumi_qqbar(1,:) = lumi_qqbar(1,:) + PartonLuminosity(grid, P_x_pdf1(:,i), pdf2(: ,-i)) 
          lumi_qqbar(2,:) = lumi_qqbar(2,:) + PartonLuminosity(grid, pdf1(:,i), P_x_pdf2(: ,-i))
        end do
      case(id_Z)
        lumi_qqbar(1,:) = 0
        lumi_qqbar(2,:) = 0
        do i = 1, 6
          lumi_qqbar(1,:) = lumi_qqbar(1,:) + gv2_ga2(i) * PartonLuminosity(grid, P_x_pdf1(:, i), pdf2(:,-i))
          lumi_qqbar(2, :) = lumi_qqbar(2,:) + gv2_ga2(i) * PartonLuminosity(grid, pdf1(:, i), P_x_pdf2(:,-i))
          lumi_qqbar(1,:) = lumi_qqbar(1,:) + gv2_ga2(i) * PartonLuminosity(grid, P_x_pdf1(:,-i), pdf2(:, i))
          lumi_qqbar(2,:) = lumi_qqbar(2,:) + gv2_ga2(i) * PartonLuminosity(grid, pdf1(:,-i), P_x_pdf2(:, i))
        end do
        lumi_gq(1,:) = 0
        lumi_gq(2,:) = 0
        lumi_qg(1,:) = 0
        lumi_qg(2,:) = 0
        do i = 1, 6
          lumi_gq(1,:) = lumi_gq(1,:) + gv2_ga2(i) * (PartonLuminosity(grid, P_x_pdf1(:,iflv_g), &
               & pdf2(:,-i)+pdf2(:,i)))
          lumi_gq(2,:) = lumi_gq(2,:) + gv2_ga2(i) * (PartonLuminosity(grid, &
               & pdf1(:,iflv_g), P_x_pdf2(:,-i) + P_x_pdf2(:,i)))
          lumi_qg(1,:) = lumi_qg(1,:) + gv2_ga2(i) * (PartonLuminosity(grid, &
               & P_x_pdf1(:,-i) + P_x_pdf1(:, i), pdf2(:,iflv_g)))
          lumi_qg(2,:) = lumi_qg(2,:) + gv2_ga2(i) * (PartonLuminosity(grid, &
               & pdf1(:,-i) + pdf1(:, i), P_x_pdf2(:,iflv_g)))
        end do
        lumi_gg(1,:) = 0
        lumi_gg(2,:) = 0
      case(id_bbH)
        lumi_qqbar(1,:) = PartonLuminosity(grid, P_x_pdf1(:,5), pdf2(:,-5)) + &
              & PartonLuminosity(grid, P_x_pdf1(:,-5), pdf2(:,5))
        lumi_qqbar(2,:) = PartonLuminosity(grid, pdf1(:,5), P_x_pdf2(:,-5)) + &
              & PartonLuminosity(grid, pdf1(:,-5), P_x_pdf2(:,5))
        lumi_gq(1,:) = PartonLuminosity(grid, P_x_pdf1(:,iflv_g), &
              & pdf2(:,-5) + pdf2(:,5)) 
        lumi_gq(2,:) = PartonLuminosity(grid, &
              & pdf1(:,iflv_g), P_x_pdf2(:,-5) +P_x_pdf2(:,5))
        lumi_qg(1,:) = PartonLuminosity(grid, &
              & P_x_pdf1(:,-5) + P_x_pdf1(:, 5), pdf2(:,iflv_g)) 
        lumi_qg(2,:) = PartonLuminosity(grid, pdf1(:,-5) + pdf1(:, 5), &
              & P_x_pdf2(:,iflv_g))
        lumi_gg(1,:) = 0
        lumi_gg(2,:) = 0
      case (id_user)
        ! Call user interface 
        ! call user_luminosities(grid, pdf1, pdf2, lumi_gg, lumi_qg, lumi_gq, lumi_qqbar)
        call wae_error('Process not implemented')
    end select

  end subroutine luminosities_expanded

!=======================================================================================
! Returns NLO luminosities for Higgs production (or call user defined function) 
! Evaluated at shat/s

  subroutine luminosities(muF, collider, lumi_gg, lumi_qg, lumi_gq, lumi_qqbar)
    use ew_parameters, only : gv2_ga2 
    !use hboson
    use user_interface, only : user_luminosities 
    use common_vars 

    real(dp), intent(in) :: muF
    character(len=*), intent(in) :: collider
    !integer, intent(in) :: iproc

    real(dp), intent(out) :: lumi_gg(0:), lumi_qg(0:), lumi_gq(0:), lumi_qqbar(0:)

    real(dp) :: pdf1(0:grid%ny,-6:7), pdf2(0:grid%ny,-6:7)
    integer  :: i

    call get_pdfs(muF, collider, pdf1, pdf2)

    select case(iproc)
      case(id_H)
        lumi_gg = PartonLuminosity(grid, pdf1(:,iflv_g), pdf2(:,iflv_g))
        lumi_gq = PartonLuminosity(grid, pdf1(:,iflv_g), &
              & sum(pdf2(:,-6:-1), dim = 2) + sum(pdf2(:, 1:6 ), dim = 2) )
        lumi_qg = PartonLuminosity(grid, &
              & sum(pdf1(:,-6:-1), dim = 2) + sum(pdf1(:, 1:6 ), dim = 2), pdf2(:,iflv_g))
        lumi_qqbar = 0
        do i = -6, 6
          if (i == 0) cycle
          lumi_qqbar = lumi_qqbar + PartonLuminosity(grid, pdf1(:,i), pdf2(:,-i))
        end do
      case(id_Z)
        lumi_qqbar = 0
        do i = 1, 6
          lumi_qqbar = lumi_qqbar + gv2_ga2(i) * PartonLuminosity(grid, pdf1(:, i), pdf2(:,-i))
          lumi_qqbar = lumi_qqbar + gv2_ga2(i) * PartonLuminosity(grid, pdf1(:,-i), pdf2(:, i))
        end do
        lumi_gq = 0; lumi_qg = 0;
        do i = 1, 6
          lumi_gq = lumi_gq+gv2_ga2(i) * PartonLuminosity(grid, pdf1(:,iflv_g), &
                & pdf2(:,-i)+pdf2(:,i))
          lumi_qg = lumi_qg + gv2_ga2(i) * PartonLuminosity(grid, &
                & pdf1(:,-i) + pdf1(:, i), pdf2(:,iflv_g))
        end do
        lumi_gg = 0
      case(id_bbH)
        lumi_qqbar = PartonLuminosity(grid, pdf1(:,5), pdf2(:,-5)) + &
              & PartonLuminosity(grid, pdf1(:,-5), pdf2(:,5))
        lumi_gq = PartonLuminosity(grid, pdf1(:,iflv_g), & 
              & pdf2(:,-5) + pdf2(:,5))
        lumi_qg = PartonLuminosity(grid, &
              & pdf1(:,-5) + pdf1(:, 5), pdf2(:,iflv_g))
        lumi_gg = 0
      case (id_user)
        ! Call user interface 
        call user_luminosities(grid, pdf1, pdf2, lumi_gg, lumi_qg, lumi_gq, lumi_qqbar) 
    end select

  end subroutine luminosities

! =======================================================================================
! Hack: for ISR terms

  subroutine lumi_Z_coeffs(muF,collider, lumi_gg, lumi_qg, lumi_gq, lumi_qqbar)
    use ew_parameters, only : gv2_ga2
    use common_vars
    use user_interface, only : user_luminosities
    use initial_state_rad
    
    implicit none

    real(dp), intent(in)  :: muF
    character(len=*), intent(in) :: collider
    real(dp), intent(out) :: lumi_gg(0:grid%ny), lumi_qg(0:grid%ny), &
                           lumi_gq(0:grid%ny), lumi_qqbar(0:grid%ny)

    real(dp) :: pdf1(0:grid%ny,-6:7), pdf2(0:grid%ny,-6:7)
    real(dp) :: Z_x_pdf1(0:grid%ny,-6:7), Z_x_pdf2(0:grid%ny,-6:7)
    type(split_mat) :: ISR_matrix, ISR
    integer :: i
  
    call get_pdfs(muF, collider, pdf1, pdf2)
    call InitISRMatrix(grid, ISR_matrix)
  
    Z_x_pdf1 = ISR_matrix .conv. pdf1
    Z_x_pdf2 = ISR_matrix .conv. pdf2 



    select case(iproc)
      case(id_H)
        lumi_gg = PartonLuminosity(grid, Z_x_pdf1(:,iflv_g), pdf2(:,iflv_g)) + &
          & PartonLuminosity(grid, pdf1(:,iflv_g), Z_x_pdf2(:,iflv_g))
        lumi_gq = PartonLuminosity(grid, Z_x_pdf1(:,iflv_g), &
          & sum(pdf2(:,-6:-1), dim = 2) + sum(pdf2(:, 1:6 ), dim = 2) ) + &
          & PartonLuminosity(grid, pdf1(:,iflv_g), &
          & sum(Z_x_pdf2(:,-6:-1), dim = 2) + sum(Z_x_pdf2(:, 1:6 ), dim = 2) )
        lumi_qg = PartonLuminosity(grid, &
          & sum(Z_x_pdf1(:,-6:-1),dim = 2) + sum(Z_x_pdf1(:, 1:6 ), dim = 2),pdf2(:,iflv_g)) + &
          & PartonLuminosity(grid, &
          & sum(pdf1(:,-6:-1),dim = 2) + sum(pdf1(:, 1:6 ), dim = 2),Z_x_pdf2(:,iflv_g))
        lumi_qqbar = 0
        do i = -6, 6
          if (i == 0) cycle
          lumi_qqbar = lumi_qqbar + PartonLuminosity(grid, Z_x_pdf1(:,i), pdf2(: ,-i)) + &
          & PartonLuminosity(grid, pdf1(:,i), Z_x_pdf2(: ,-i))
        end do
    
      case(id_Z)
        lumi_qqbar = 0
        do i = 1, 5
          lumi_qqbar = lumi_qqbar + gv2_ga2(i) * (PartonLuminosity(grid, Z_x_pdf1(:, i), pdf2(:,-i)) + &
             & PartonLuminosity(grid, pdf1(:, i), Z_x_pdf2(:,-i)))
           lumi_qqbar = lumi_qqbar + gv2_ga2(i) * (PartonLuminosity(grid, Z_x_pdf1(:,-i), pdf2(:, i)) + &
             & PartonLuminosity(grid, pdf1(:,-i),Z_x_pdf2(:, i)))
        end do
      
        lumi_gq = 0; lumi_qg = 0;
        
        do i = 1, 5
          lumi_gq = lumi_gq + gv2_ga2(i) * (PartonLuminosity(grid, Z_x_pdf1(:,iflv_g), &
            & pdf2(:,-i)+pdf2(:,i)) + PartonLuminosity(grid, &
            & pdf1(:,iflv_g), Z_x_pdf2(:,-i) + Z_x_pdf2(:,i)))
          lumi_qg = lumi_qg + gv2_ga2(i) * (PartonLuminosity(grid, &
            & Z_x_pdf1(:,-i) + Z_x_pdf1(:, i), pdf2(:,iflv_g)) + &
            & PartonLuminosity(grid, &
            & pdf1(:,-i) + pdf1(:, i),Z_x_pdf2(:,iflv_g)))
        end do
          lumi_gg = 0 
    end select

  end subroutine lumi_Z_coeffs

! =======================================================================================
! Hack: for DFrec terms (thrust minor only)

  subroutine lumi_DFrec_coeffs(muF,collider, lumi_gg, lumi_qg, lumi_gq, lumi_qqbar)
    use ew_parameters, only : gv2_ga2
    use common_vars
    use user_interface, only : user_luminosities
    use initial_state_rad
    
    implicit none

    real(dp), intent(in)  :: muF
    character(len=*), intent(in) :: collider
    real(dp), intent(out) :: lumi_gg(0:grid%ny), lumi_qg(0:grid%ny), &
                           lumi_gq(0:grid%ny), lumi_qqbar(0:grid%ny)

    real(dp) :: pdf1(0:grid%ny,-6:7), pdf2(0:grid%ny,-6:7)
    real(dp) :: DFrec_x_pdf1(0:grid%ny,-6:7), DFrec_x_pdf2(0:grid%ny,-6:7)
    type(split_mat) :: DFrec_matrix, DFrec
    integer :: i
  
    call get_pdfs(muF, collider, pdf1, pdf2)
    call Init_DFrec_Matrix(grid, DFrec_matrix)
  
    DFrec_x_pdf1 = DFrec_matrix .conv. pdf1
    DFrec_x_pdf2 = DFrec_matrix .conv. pdf2 



    select case(iproc)
      case(id_H)
        lumi_gg = PartonLuminosity(grid, DFrec_x_pdf1(:,iflv_g), pdf2(:,iflv_g)) + &
          & PartonLuminosity(grid, pdf1(:,iflv_g), DFrec_x_pdf2(:,iflv_g))
        lumi_gq = PartonLuminosity(grid, DFrec_x_pdf1(:,iflv_g), &
          & sum(pdf2(:,-6:-1), dim = 2) + sum(pdf2(:, 1:6 ), dim = 2) ) + &
          & PartonLuminosity(grid, pdf1(:,iflv_g), &
          & sum(DFrec_x_pdf2(:,-6:-1), dim = 2) + sum(DFrec_x_pdf2(:, 1:6 ), dim = 2) )
        lumi_qg = PartonLuminosity(grid, &
          & sum(DFrec_x_pdf1(:,-6:-1),dim = 2) + sum(DFrec_x_pdf1(:, 1:6 ), dim = 2),pdf2(:,iflv_g)) + &
          & PartonLuminosity(grid, &
          & sum(pdf1(:,-6:-1),dim = 2) + sum(pdf1(:, 1:6 ), dim = 2),DFrec_x_pdf2(:,iflv_g))
        lumi_qqbar = 0
        do i = -6, 6
          if (i == 0) cycle
          lumi_qqbar = lumi_qqbar + PartonLuminosity(grid, DFrec_x_pdf1(:,i), pdf2(: ,-i)) + &
          & PartonLuminosity(grid, pdf1(:,i), DFrec_x_pdf2(: ,-i))
        end do
    
      case(id_Z)
        lumi_qqbar = 0
        do i = 1, 5
          lumi_qqbar = lumi_qqbar + gv2_ga2(i) * (PartonLuminosity(grid, DFrec_x_pdf1(:, i), pdf2(:,-i)) + &
            & PartonLuminosity(grid, pdf1(:, i), DFrec_x_pdf2(:,-i)))
          lumi_qqbar = lumi_qqbar + gv2_ga2(i) * (PartonLuminosity(grid, DFrec_x_pdf1(:,-i), pdf2(:, i)) + &
            & PartonLuminosity(grid, pdf1(:,-i),DFrec_x_pdf2(:, i)))
        end do
        lumi_gq = 0; lumi_qg = 0;
        do i = 1, 5
          lumi_gq = lumi_gq + gv2_ga2(i) * (PartonLuminosity(grid, DFrec_x_pdf1(:,iflv_g), &
            & pdf2(:,-i)+pdf2(:,i)) + PartonLuminosity(grid, &
            & pdf1(:,iflv_g), DFrec_x_pdf2(:,-i) + DFrec_x_pdf2(:,i)))
          lumi_qg = lumi_qg + gv2_ga2(i) * (PartonLuminosity(grid, &
            & DFrec_x_pdf1(:,-i) + DFrec_x_pdf1(:, i), pdf2(:,iflv_g)) + &
            & PartonLuminosity(grid, &
            & pdf1(:,-i) + pdf1(:, i),DFrec_x_pdf2(:,iflv_g)))
        end do
          lumi_gg = 0 
    end select

  end subroutine lumi_DFrec_coeffs

end module pdfs_tools

