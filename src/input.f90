!========================================================
!--------------------------------------------------------
! Module handling various user inputs  
!--------------------------------------------------------
!========================================================

module input

  use hoppet_v1 
  use sub_defs_io
  use common_vars 

  implicit none

  private

  ! The process as a character string -- only used here 
  ! iproc (defined in common_vars.f90) is used in the code instead 
  character(len=8) :: proc 

  ! Strategy for the setting of the scales mu_R and mu_F 
  character(len=2) :: scale_strategy

  ! Collider type (pp or ppbar) (only needed here and in h1jet.f90) 
  character(len=5), public :: collider 
  ! pt_min and pt_max of the histogram (only needed here and in h1jet.f90) 
  real(dp), public :: ptmin, ptmax
  ! Number of bins (only needed here and h1jet in.f90)
  integer, public :: nbins 
  ! Strong coupling exponent power (only needed here and in h1jet.f90) 
  integer, public :: as_pow = 0
  ! Total born-level cross-section (only needed here and in h1jet.f90) 
  real(dp), public :: sigma0 
  ! Factorisation and renormalisation factors (only needed here and in h1jet.f90) 
  real(dp), public :: muF, muR
  ! Scale factors (only needed here and in h1jet.f90) 
  real(dp), public :: xmur, xmuf
  ! PDF related variables (only needed here and in h1jet.f90) 
  character(len=30), public :: pdf_name
  integer, public :: pdf_mem

  ! Model dependent parameters 
  real(dp) :: tbeta, sth2
 
  ! Identifiers for specific composite Higgs models  
  integer, parameter :: M1_5 = 1, M1_14 = 2, M4_5 = 3, M4_14 = 4

  public :: input_handler, set_scale, print_settings 

contains

!======================================================================================= 
! Handle all user input 

  subroutine input_handler(idev)
    use hboson; use ew_parameters 
    integer, intent(in) :: idev

    ! Process 
    proc = trim(string_val_opt('--proc', 'H'))
    iproc = proc_to_id(proc)

    ! Set the alpha_s power depending on the process 
    if (iproc == id_H) then 
      ! Result multiplied by factor alpha_s^(as_pow + 1) = alpha_s^3 
      as_pow = 2 
    else if (iproc == id_user) then 
      ! Result multiplied by factor alpha_s^0 = 1, since the 
      ! appropriate alpha_s factor is already included 
      as_pow = -1 
    end if 

    ! Centre-of-mass collision energy 
    roots = dble_val_opt('--roots', 13000._dp) 
    
    ! Set up collider type 
    collider = trim(string_val_opt('--collider', 'pp'))
    if (collider /= 'pp' .and. collider /= 'ppbar') then
      call wae_error('handle_input', 'Unrecognised collider type &
                &-- Please specify collider as pp or ppbar')
    end if 

    ! Set SM and EW parameters if changed (from ew_parameters.f90) 
    ! G_mu scheme 
    mz_in = dble_val_opt('--mZ', mz)
    mw_in = dble_val_opt('--mW', mw)
    GF_GeVm2_in = dble_val_opt('--GF', GF_GeVm2)
    sinwsq_in = one - (mw_in**2) / (mz_in**2) 
    higgs_vev_in = one / dsqrt(dsqrt(two) * GF_GeVm2_in) 
    call set_gv2 

    ! Set heavy quark masses if changed (from mass_helper.f90) 
    mt_in = dble_val_opt('--mt', mt)
    mb_in = dble_val_opt('--mb', mb)
    
    ! Set PDF set 
    pdf_name = trim(string_val_opt('--pdf_name', 'MSTW2008nlo68cl')) !//'.LHgrid'
    pdf_mem  = int_val_opt('--pdf_mem', 0)

    call set_masses_and_couplings

    ! Set up factorisation and renormalisation scales 
    scale_strategy = trim(string_val_opt('--scale_strategy', 'HT'))

    muF = set_scale(zero)
    muR = muF

    xmur = dble_val_opt('--xmur', 0.5_dp)
    xmuf = dble_val_opt('--xmuf', 0.5_dp)

    muF = muF * xmuf
    muR = muR * xmur

    ! Histogram options 
    nbins = int_val_opt('--nbins', 400)
    ptmin = dble_val_opt('--ptmin', zero)
    ptmax = dble_val_opt('--ptmax', 4e3_dp)

    if (log_val_opt('--log') .and. (ptmin == zero)) then 
      call wae_error('handle_input', 'Set minimum pT value (with --ptmin) &
                    &when using --log')
    end if

    ! Hack: Resummation coefficients control flag
    select case (trim(string_val_opt('--resum', 'none')))
    case ('none')
      resum_flag = resum_none
    case ('h12')
      resum_flag = resum_h12
    case ('h11')
      resum_flag = resum_h11
    case ('h24')
      resum_flag = resum_h24
    case ('h23')
      resum_flag = resum_h23
    case default
      call wae_error('input_handler', 'Unrecognised --resum option (use none/h12/h11/h24)')
    end select

  end subroutine input_handler
  
!======================================================================================= 
! Convert --proc command option argument to internal id 

  function proc_to_id(proc) result(res)
    character(len=*) :: proc
    integer :: res
    
    if (proc == 'H') then
      res = id_H
    else if (proc == 'Z') then
      res = id_Z
    else if (proc == 'bbH') then
      res = id_bbH
    else  if (proc == 'user') then
      res = id_user
    else
      call wae_error('proc_to_id', 'Unrecognised process ', proc)
    end if

  end function proc_to_id
  
!======================================================================================= 
! Convert --model command option argument to internal id 

  function model_to_id(model) result(res)
    character(len=*) :: model
    integer :: res

    if (model == 'M1_5') then
      res = M1_5
    else if (model == 'M1_14') then
      res = M1_14
    else if (model == 'M4_5') then
      res = M4_5
    else  if (model == 'M4_14') then
      res = M4_14
    else
      call wae_error('model_to_id', 'Unrecognised model ', model)
    end if
    
  end function model_to_id
  
!=======================================================================================
! Function to set the factorisation and renormalisation scales according to the 
! specified scale strategy and given pT 

  function set_scale(ptval) result(res)
    real(dp), intent(in) :: ptval
    real(dp) :: res

    if (scale_strategy == 'M') then
      res = M
    else if (scale_strategy == 'HT') then
      res = ptval + sqrt(ptval**2 + M**2)
    else if (scale_strategy == 'MT') then
      res = sqrt(ptval**2 + M**2)
    else
      call wae_error('set_scale', 'Unrecognised scale strategy: ', scale_strategy)
    end if

  end function set_scale

!======================================================================================= 
! Function to set the relevant masses and parameters depending on the selected process 

  subroutine set_masses_and_couplings
    use hboson; use ew_parameters
    real(dp) :: invfscale, imc1
    integer  :: model
    integer  :: indev, ios, nq
    character(len=7) :: nq_eq

    ! Model-specific parameters 
    real(dp) :: yt, yb, ytp
    real(dp) :: mtp, cth2, s2th2
    real(dp) :: yst1, yst2
    real(dp) :: delta, alpha1, alpha2, c2beta
    real(dp) :: mst1, mst2 

    select case(iproc)

      ! pp -> H + jet 
      case(id_H)

        ! Check for CP odd Higgs 
        cpodd = log_val_opt('--cpodd') 
        
        ! Standard parameters 
        M = dble_val_opt('--mH', mh) ! Higgs mass 
        yt = dble_val_opt('--yt', one) ! Top Yukawa factor 
        if (cpodd) then
          yb = dble_val_opt('--yb', zero) ! Bottom CP-odd Yukawa factor
        else
          yb = dble_val_opt('--yb', one) ! Bottom Yukawa factor
        end if
        sth2 = dble_val_opt('--sth2', zero) ! Stop/Top-partner mixing angle 

        ! Top partner mass and Yukawa coupling 
        mtp = dble_val_opt('--mtp', zero)
        ytp = dble_val_opt('--ytp', one)
        
        ! Stop squark mass and SUSY parameters 
        mst1 = dble_val_opt('--mst', zero) ! Stop mass 
        delta = dble_val_opt('--delta', zero) ! SUSY fine-tuning parameter 
        mst2 = sqrt(mst1**2 + delta**2)
        tbeta = dble_val_opt('--tbeta', zero) ! Ratio of VEV's in MHDMs, tan(beta) 
        c2beta = (one - tbeta**2) / (one + tbeta**2) ! cos(2 * beta) 

        if ((mst1 * mtp) /= zero) then 
          call wae_error('set_masses_and_couplings', 'Simoultaneous top partners &
                    &and stops not allowed') 
        end if 

        if (log_val_opt('--in') .or. log_val_opt('-i')) then
          ! Include multiple top-partners in quark loops 
          if (mst1 /= zero) then 
            call wae_error('set_masses_and_couplings', 'Simoultaneous top partners &
                    &and stops not allowed') 
          end if

          if (log_val_opt('--in')) then
            indev = idev_open_opt('--in', status = "old")
          else if (log_val_opt('-i')) then
            indev = idev_open_opt('-i', status = "old")
          end if

          read(indev, *, iostat = ios) nq_eq, nq
          if (ios /= 0) then 
            call wae_error('set_masses_and_couplings', 'Incorrect file format, see SM.dat for an example')
          end if 
           
          call read_top_partners(indev, nq)

        else if (mtp /= zero) then
          ! Include top partner in quark loops 
          allocate(mass_array(3), yukawa_array(3), iloop_array(3))

          mass_array = (/mt_in, mb_in, mtp/)

          ! The inverse of fscale
          if (log_val_opt('--fscale')) then
            invfscale = one / dble_val_opt('--fscale', zero)
          else if (log_val_opt('-f')) then
            invfscale = one / dble_val_opt('-f', zero)
          else
            invfscale = zero
            yt = yt * (one - sth2)
            ytp = ytp * sth2
          end if

          ! In case of top partner for composite Higgs models 
          if (invfscale /= zero) then
            model = model_to_id(string_val_opt('--model', 'M1_5'))
            imc1 = dble_val_opt('--imc1', zero)
            call set_yukawas(model, invfscale, imc1, mtp, yt, yb, ytp)
          end if

          yukawa_array = (/yt, yb, ytp/)

          iloop_array = iloop_fm_fermion

        else if (mst1 /= zero) then
          ! Include SUSY stop squark in loops 
          cth2 = one - sth2
          s2th2 = four * sth2 * cth2
          alpha1 = mz_in**2 / mt_in**2 * c2beta * (one - four / three * sinwsq_in)
          alpha2 = mz_in**2 / mt_in**2 * c2beta * four / three * sinwsq_in
          yst1 = (mt_in / mst1)**2 * (alpha1 * cth2 + alpha2 * sth2 + two - delta**2 / two / mt_in**2 * s2th2)
          yst2 = (mt_in / mst2)**2 * (alpha1 * sth2 + alpha2 * cth2 + two + delta**2 / two / mt_in**2 * s2th2)
          allocate(mass_array(4), yukawa_array(4), iloop_array(4))
          mass_array = (/mt_in, mb_in, mst1, mst2/)
          yukawa_array = (/yt, yb, yst1, yst2/)
          iloop_array(1:2) = iloop_fm_fermion
          iloop_array(3:4) = iloop_fm_scalar
        else
          ! SM, only bottom and top in quark loops 
          allocate(mass_array(2), yukawa_array(2), iloop_array(2))
          mass_array = (/mt_in, mb_in/)
          yukawa_array = (/yt, yb/)
          iloop_array = iloop_fm_fermion
        end if

        ! Reset approximation used to compute loops of particles, the
        ! user needs to know the meaning of the various approximations
        call reset_iloop_array
        
      ! pp -> Z + jet   
      case(id_Z)

        M = mz_in
      
      ! bbbar -> H + jet 
      case(id_bbH)

        M = dble_val_opt('--mH', mh) ! Higgs mass 
        allocate(mass_array(1), yukawa_array(1))
        mass_array(1) = dble_val_opt('--mbmb', mbmb)

      ! User specified process 
      case (id_user)

        ! Relevant mass for the user specified amplitude 
        ! Important to set for correct scales muR and muF 
        M = zero
        if (log_val_opt('--mass')) then
          M = dble_val_opt('--mass', zero) 
        else if (log_val_opt('-M')) then
          M = dble_val_opt('-M', zero)
        else 
          call wae_warn(max_warns, 'set_masses_and_couplings', 'It is &
          &recommended to set relevant mass with option -M or --mass for user-&
          &specified processes')
        end if

      case default

        call wae_error('set_masses_and_couplings', 'Unrecognised process ', proc)

    end select

  end subroutine set_masses_and_couplings

!======================================================================================= 
! Set Yukawa couplings for specific composite Higgs models

  subroutine set_yukawas(model, invfscale, imc1, mtp, yt, yb, ytp)
    use hboson
    use ew_parameters
    integer, intent(in) :: model
    real(dp), intent(in) :: invfscale, imc1, mtp
    real(dp), intent(inout) :: yt, yb, ytp
    !-----------------------------------------  
    real(dp) :: seps, ceps
    real(dp) :: sthRsq, cthRsq, tanthRsq
    real(dp) :: sthLsq, cthLsq, tanthLsq

    seps = higgs_vev_in * invfscale
    ceps = sqrt(one - seps**2)

    select case(model)

      case(M1_5)
        sthLsq = sth2            
        cthLsq = one - sth2 
        if (cpodd) then
          yt = zero
          ytp = zero
        else
          yt = yt * cthLsq * ceps
          ytp = ytp * sthLsq * ceps
          yb = yb * ceps
        end if

      case(M1_14)
        sthLsq = sth2            
        cthLsq = one - sth2 
        if (cpodd) then
          yt = zero
          ytp = zero
          yb = zero
        else
          yt = yt * cthLsq * (two * ceps**2 - one) / ceps
          ytp = ytp * sthLsq * (two * ceps**2 - one) / ceps
          yb = yb * (two * ceps**2 - one) / ceps
        end if
        
      case(M4_5)
        sthRsq = sth2            
        cthRsq = one - sth2 
        tanthLsq = (mtp / mt)**2 * sthRsq / cthRsq
        cthLsq = one / (one + tanthLsq)
        sthLsq = one - cthLsq 
        if (cpodd) then
          yt = four * ceps * seps / sqrt(two * (one + ceps**2)) * &
                & imc1 * sqrt(sthRsq * cthRsq)
          ytp = -yt
          yb = zero
        else
          yt = yt * ceps * (cthRsq - seps**2 / (one + ceps**2) * (cthLsq - cthRsq))
          ytp = ytp * ceps * (sthRsq - seps**2 / (one + ceps**2) * (sthLsq - sthRsq))
          yb = yb * ceps
        end if

      case(M4_14)
        sthRsq = sth2            
        cthRsq = one - sth2 
        tanthLsq = (mtp / mt)**2 * sthRsq / cthRsq
        cthLsq = one / (one + tanthLsq)
        sthLsq = one - cthLsq 
        if (cpodd) then
          yt = four * seps * (one - two * seps**2) / &
                & sqrt(two * (four * ceps**4 - three * ceps**2 + one)) * &
                & imc1 * sqrt(sthRsq * cthRsq)
          ytp = -yt
          yb = zero
        else
          yt = yt * (cthRsq * (two * ceps**2 - one) / ceps - &
                & seps**2 * ceps / (four * ceps**4 - three * ceps**2 + one) * &
                & (8._dp * ceps**2 - three) * (cthLsq - cthRsq))
          ytp = ytp * (sthRsq * (two * ceps**2 - one) / ceps - &
                & seps**2 * ceps / (four * ceps**4 - three * ceps**2 + one) * &
                & (8._dp * ceps**2 - three) * (sthLsq - sthRsq))
          yb = yb * (two * ceps**2 - one) / ceps
        end if

      case default
        call wae_error('set_yukawas', 'Unrecognised model: ', &
              & intval = model)

    end select

  end subroutine set_yukawas

!======================================================================================= 
! Read top partners from input file 

  subroutine read_top_partners(indev, nq)
    use hboson
    integer, intent(in) :: indev, nq
    !---------------------------
    integer :: i, ios
    real(dp) :: kappa, kappa_tilde
    
    allocate(mass_array(nq), yukawa_array(nq), iloop_array(nq))

    do i = 1, nq
      read(indev, *, iostat = ios) mass_array(i), kappa, kappa_tilde, iloop_array(i)

      if (ios /= 0) then 
        call wae_error('read_top_partners', 'Incorrect file format, see SM.dat for an example')
      end if 

      if (cpodd) then
        yukawa_array(i) = kappa_tilde
      else
        yukawa_array(i) = kappa
      end if
    end do

    close(indev)

  end subroutine read_top_partners

!======================================================================================
! Change approximations used to compute loops in Higgs production
! The user needs to change those at compile time, because knowledge
! of the meaning of the approximation is strictly required
  
  subroutine reset_iloop_array
    use hboson
    
    !! Example: effective theory with SM top and bottom, and one
    !! infinitely heavy top partner, equivalent to adding a dimension-6
    !! contact interaction
    !    if (size(iloop_array) == 3) then
    !       iloop_array = (/ iloop_fm_fermion, iloop_fm_fermion, iloop_lm_fermion /)
    !    else
    !       call wae_error('reset_iloop_array', 'Expected size of iloop_array is 3, whereas actual &
    !            &one is', intval = size(iloop_array))
    !    end if

  end subroutine reset_iloop_array
  
!======================================================================================= 
! Print the settings from user input 

  subroutine print_settings(idev)
    use hboson
    integer, intent(in) :: idev
    
    write(idev,*) ! Blank line for nicer output 
    write(idev,*) 'collider       = ', collider
    write(idev,*) 'roots(GeV)     =', roots
    write(idev,*) 'process        = ', proc
    select case(iproc)
      case(id_user)
        write(idev,*) 'model          = User-defined from included custom code' 
      case(id_H)
        if (cpodd) then
          write(idev,*) 'model          = CP-odd Higgs'
        else if (size(mass_array)==2) then
          write(idev,*) 'model          = SM' 
        else if (iloop_array(3) <= iloop_lm_fermion) then
          write(idev,*) 'model          = SM with top partners' 
        else
          write(idev,*) 'model          = SUSY' 
        end if
      end select
    write(idev,*) 'mass           =', M
    write(idev,*) 'muR            =', muR
    write(idev,*) 'muF            =', muF
    write(idev,*) 'alphas(muR)    =', alphas
    write(idev,*) 'scale strategy = ', scale_strategy 
    if (iproc == id_H) then 
      write(idev,*) 'masses         =', mass_array
      write(idev,*) 'tan(beta)      =', tbeta
      write(idev,*) 'sin^2(theta)   =', sth2
    end if
    write(idev,*) 'sigma0 [nb]    =', sigma0 
    write(idev,*) 'log            =', log_val_opt('--log')
    write(idev,*) ! Blank line for nicer output 

  end subroutine print_settings 

!=======================================================================================

end module input 


