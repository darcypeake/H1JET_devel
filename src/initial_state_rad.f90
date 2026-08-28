!! Code to calculate the initial-state radiation terms using hoppet

module initial_state_rad
  use types; use consts_dp; use convolution_communicator; use convolution
  use dglap_objects
  use qcd
  implicit none
  private

  public :: InitISRMatrix, Init_DFrec_Matrix

contains

  function Z1qq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
      ! ------- one - jettiness  ----------
      res = cf * (1+x**2)/(1-x) * log(1-x) - cf * (1+x**2)/(1-x)*log(x)
      res = res + cf * (1-x)
      

      ! ---- thrust minor ----
      !res = cf * (1-x) 

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
    ! -------- one - jettiness ---------
    res = res - two * cf/(1-x) * log(1-x)
    
    case(cc_DELTA)
    res = zero
    end select    
    
    if (cc_piece /= cc_DELTA) res = res * x
  end function Z1qq

  ! --------------------------------------------------------------------!

  function Z1qg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
      ! -------- one - jettiness ---------
      res = tr * (x**2 + (1 - x)**2)*log((1-x)/x) 
      res = res + two * tr * x * (1-x)

      ! ---- thrust minor ----
      !res = two  * tr * x * (1-x)

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
    res = zero
    end select 

    if (cc_piece /= cc_DELTA) res = res * x
  end function Z1qg
 ! --------------------------------------------------------------------!
  function Z1gq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)                 :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
      ! -------- one - jettiness ---------
      res = cf * (1 + (1 - x)**2)/x * log((1-x)/x) 
      res = res + cf * x

      ! --- thrust minor ----
      !res = cf * x

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
    res = zero
    end select 

    if (cc_piece /= cc_DELTA) res = res * x

  end function Z1gq

   ! --------------------------------------------------------------------!
  function Z1gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)                 :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
    ! --- thrust minor ----
    !res = zero 

    ! -------- one - jettiness ---------
    res = two * ca * (1-x+x**2)**2/x * (log(1-x)/(1-x)-log(x)/(1-x))

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
      ! -------- one - jettiness ---------
      res = res - two * ca * log(1-x)/(1-x)
      
    case(cc_DELTA)
     res = zero
    end select 

    if (cc_piece /= cc_DELTA) res = res * x


  end function Z1gg
  
  
  ! --------------------------------------------------------------------!
  subroutine InitISRMatrix(grid, ISR)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: ISR

    ISR%nf_int = nf_int

    call cobj_InitSplitLinks(ISR)

    call InitGridConv(grid, ISR%qq, Z1qq)
    call InitGridConv(grid, ISR%qg, Z1qg)
    call InitGridConv(grid, ISR%gq, Z1gq)
    call InitGridConv(grid, ISR%gg, Z1gg)

    call Multiply(ISR%qg, 2*nf)

    call InitGridConv(ISR%NS_plus,  ISR%qq)
    call InitGridConv(ISR%NS_minus, ISR%qq)
    
    call InitGridConv(ISR%NS_V, ISR%NS_minus)

  end subroutine InitISRMatrix

  ! --------------------------------------------------------------------!

  function DFrec_qq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
      res = zero
    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
    ! no virtual piece
    case(cc_DELTA)
    res = zero
    end select    
    
    if (cc_piece /= cc_DELTA) res = res * x
  end function DFrec_qq

  ! --------------------------------------------------------------------!

  function DFrec_qg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)             :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
      res = zero

    end select
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
    res = zero
    end select 

    if (cc_piece /= cc_DELTA) res = res * x
  end function DFrec_qg
 ! --------------------------------------------------------------------!
  function DFrec_gq(y) result(res)
    real(dp), intent(in) :: y
    real(dp)                 :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
      res = res + cf * four * (1 - x)/x ! DFrec
    end select
    
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
       ! no virtual piece
    case(cc_DELTA)
    res = zero
    end select 

    if (cc_piece /= cc_DELTA) res = res * x

  end function DFrec_gq

   ! --------------------------------------------------------------------!
  function DFrec_gg(y) result(res)
    real(dp), intent(in) :: y
    real(dp)                 :: res
    real(dp)             :: x
    x = exp(-y)
    res = zero

    select case(cc_piece)
    case(cc_REAL,cc_REALVIRT)
      res = res - ca * four * (1 - x)/x
    end select
    
    select case(cc_piece)
    case(cc_VIRT,cc_REALVIRT)
    ! no virtual piece
    case(cc_DELTA)
     res = zero
    end select 

    if (cc_piece /= cc_DELTA) res = res * x


  end function DFrec_gg

  ! --------------------------------------------------------------------!

  subroutine Init_DFrec_Matrix(grid, DFrec)
    type(grid_def),  intent(in)    :: grid
    type(split_mat), intent(inout) :: DFrec

    DFrec%nf_int = nf_int

    call cobj_InitSplitLinks(DFrec)

    call InitGridConv(grid, DFrec%qq, DFrec_qq)
    call InitGridConv(grid, DFrec%qg, DFrec_qg)
    call InitGridConv(grid, DFrec%gq, DFrec_gq)
    call InitGridConv(grid, DFrec%gg, DFrec_gg)

    call Multiply(DFrec%qg, 2*nf)

    call InitGridConv(DFrec%NS_plus,  DFrec%qq)
    call InitGridConv(DFrec%NS_minus, DFrec%qq)
    
    call InitGridConv(DFrec%NS_V, DFrec%NS_minus)

  end subroutine Init_DFrec_Matrix


end module initial_state_rad
