!========================================================
!--------------------------------------------------------
!                         H1jet 
! Tool for the fast and accurate calculation of the pT 
! distribution in Higgs + jet production at NLO 
!--------------------------------------------------------
!========================================================

program h1jet

  use pdfs_tools 
  use hoppet_v1
  use sub_defs_io
  use ew_parameters
  use hboson
  use input 
  use cross_sections
  use banner
  use common_vars 
  use warnings_and_errors
  use gauss_integrator 
  
  implicit none

  ! Loop index 
  integer :: i
  ! I/O unit identifier
  integer :: idev
  ! Local kinematical variables
  real(dp) :: lnpt, xmom, ymin, ymax 
  ! Variables related to bbH 
  real(dp), pointer :: mb0
  real(dp) :: running_mb
  real(dp) :: as0
  ! Histogram related variables 
  real(dp), allocatable :: binmin(:), binmed(:), binmax(:)
  real(dp) :: bin_width
  ! Cross-sections 
  real(dp), allocatable :: dsigma_dpt(:), sigma(:)
  real(dp) :: dsigmadpt, ew_prefactor
  ! Integration accuracy 
  real(dp) :: accuracy
  ! tau = M^2 / s 
  type(gdval) :: tau

  ! Print help message or version number if invoked 
  if (log_val_opt('-h') .or. log_val_opt('--help')) then
    call print_help_message
  else if (log_val_opt('-v') .or. log_val_opt('--version')) then 
    call print_version_number 
  end if

  ! Set output level 
  if (log_val_opt('--out')) then
    idev = idev_open_opt('--out')
  else if (log_val_opt('-o')) then
    idev = idev_open_opt('-o')
  else
    idev = stdout 
  end if

  ! Print welcome banner
  call print_welcome_banner(idev)

  ! Print commandline options
  write(idev,'(a)') 'Provided command options: '
  write(idev,'(a)') '# '//trim(command_line()) 
  write(idev,*) ! Blank line for nicer output 

  ! Handle user input 
  call input_handler(idev) 
  ! Set the desired integration accuracy 
  accuracy = dble_val_opt('--accuracy', 1e-3_dp) 

  ! Set up PDFs 
  call init_pdfs_from_LHAPDF(pdf_name, pdf_mem)
  allocate(lumi_gg(0:grid%ny), lumi_qg(0:grid%ny), lumi_gq(0:grid%ny),&
       & lumi_qqbar(0:grid%ny))

  ! Hack - define the dlumi's       
       allocate(dlumi_gg(0:grid%ny), dlumi_qg(0:grid%ny), dlumi_gq(0:grid%ny),&
       & dlumi_qqbar(0:grid%ny))

  ! If proc == user then check if it is set up 
  if (iproc == id_user) then 
    select case(user_included())
      case (id_noImplementation)
        call wae_error('h1jet', 'User-specified process not set -- please recompile &
             &H1jet with USERPATH set')
      case (id_missingTotXsec)
        call wae_warn(max_warns, 'h1jet', 'Total born-level cross-section not &
             &included in user interface')
    end select
  end if

  ! Evaluate the running alpha_s for the born-level cross-section 
  alphas = RunningCoupling(muR)
  
  ! Set tau needed for the luminosities 
  tau = (M**2 / roots**2) .with. grid 

  ! Calculate born-level cross-section 
  call luminosities(muF, collider, lumi_gg, lumi_qg, lumi_gq, lumi_qqbar)
  sigma0 = cross_section(lumi_gg, lumi_qg, lumi_gq, lumi_qqbar, tau)
  if (iproc /= id_user) then 
    sigma0 = sigma0 * alphas**as_pow
  end if 

  ! bbH specific set up 
  if (iproc == id_bbH) then

    ! Get mb(mb) 
    mb0 => mass_array(1)

    ! Get alpha_s(mb(mb))
    as0 = RunningCoupling(mb0)
    
! Calculate running bottom mass 
    running_mb = RunningMass(muR, mb0, alphas, as0)

    write(idev,*) ! Blank line for nicer output 
    write(idev,*) "m_b(mu_R)      =", running_mb
    write(idev,*) "as(MuRbbh)     =", as0

    ! Scale born-level cross-section with running bottom mass 
    sigma0 = sigma0 * running_mb**2

  end if

  ! Allocate histogram and cross-section related arrays 
  allocate(binmin(nbins), binmed(nbins), binmax(nbins))
  allocate(dsigma_dpt(nbins), sigma(nbins))
 
  ! Output settings used 
  call print_settings(idev) 
  
  write(idev,*) ! Blank line for nicer output
  write(idev,*) '# cols are:   binmin   binmed   binmax   dsigma/dpt [fb/GeV]   sigma(pt) [fb]'

  ! Evaluate the process-dependent EW prefactor 
  call h1jet_prefactor(ew_prefactor)

  ! Loop over pT values (bins) 
  do i = 1, nbins

    if (log_val_opt('--log')) then
      bin_width = (log(ptmax) - log(ptmin)) / nbins
      lnpt = log(ptmin) + bin_width * i
      binmax(i) = lnpt
      binmed(i) = binmax(i) - half * bin_width
      binmin(i) = binmed(i) - half * bin_width
      pt = exp(binmed(i))
    else
      bin_width = (ptmax - ptmin) / nbins
      pt = ptmin + bin_width * i
      binmax(i) = pt
      binmed(i) = binmax(i) - half * bin_width
      binmin(i) = binmed(i) - half * bin_width
      pt = binmed(i)
    end if

    ! Update factorisation and renormalisation scales with new pT 
    muF = set_scale(pt)
    muR = muF

    muF = muF * xmuf
    muR = muR * xmur
    
    call luminosities(muF, collider, lumi_gg, lumi_qg, lumi_gq, lumi_qqbar)
    ! Hack
    call luminosities_expanded(muF, collider, dlumi_gg, dlumi_qg, dlumi_gq, dlumi_qqbar)

    ! Update running alpha_s 
    alphas = RunningCoupling(muR)

    ! Gaussian integration
    xmom = (roots**2 - M**2) / roots / pt
    ymax = log(xmom * half + half * sqrt(xmom**2 - four))
    ymin = -ymax
    dsigmadpt = gauss_integrate(dsigma_dptdy, ymin, ymax, accuracy)

    ! Include couplings and prefactors
    dsigmadpt = dsigmadpt * ew_prefactor * alphas**(as_pow + 1)
    ! Hack
    dsigmadpt = dsigmadpt * alphas/twopi
    if (iproc == id_bbH) then
      running_mb = RunningMass(muR, mb0, alphas, as0)
      dsigmadpt = dsigmadpt * running_mb**2
    end if

    dsigma_dpt(i) = dsigmadpt

  enddo

  ! Compute the integrated distribution as a function of ptmin 
  ! and print the result on the screen
  !Hack but needs to be git pushed as necessary for the programme. 
  do i = 1, nbins
    if (log_val_opt('--log')) then
       sigma(i) = sum(exp(binmed(i:nbins))*dsigma_dpt(i:nbins)) * bin_width
       write(idev,*) exp(binmin(i)), exp(binmed(i)), exp(binmax(i)), dsigma_dpt(i), sigma(i)
    else
       sigma(i) = sum(dsigma_dpt(i:nbins)) * bin_width
       write(idev,*) binmin(i), binmed(i), binmax(i), dsigma_dpt(i), sigma(i)
    end if
 end do
 write(*,*) ! Blank line for nicer output    

  ! Close output file (if output file has been specified) 
  if (idev /= stdout) then 
    close(idev)
  end if 

end program h1jet

