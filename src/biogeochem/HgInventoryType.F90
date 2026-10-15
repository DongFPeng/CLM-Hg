module HgInventoryType

    !-----------------------------------------------------------------------
    ! !DESCRIPTION:
    ! Bring Hg land/water inventory to CLM 
    !
    use shr_kind_mod      , only : r8 => shr_kind_r8
    use shr_log_mod       , only : errMsg => shr_log_errMsg
    use clm_varcon        , only : spval
    use decompMod         , only : bounds_type
    use abortutils        , only : endrun
    use ColumnType        , only : col

    !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
   ! !PUBLIC TYPES:
  type, public :: HgInvent_type

  ! Hg emission flux
  real(r8), pointer ::Hg_emission_grc(:)         ! grid Hg emission (g/km2/yr)
  real(r8), pointer ::Hg_emission_col(:)         ! col Hg emission (g/km2/yr)
  real(r8), pointer ::Hg_emission_rof(:)         ! col Hg emission to rof (ng/m2/s)


  
contains

  procedure, public  :: Init
  procedure, private :: InitAllocate
  procedure, private :: InitHistory
  procedure, private :: InitCold
  procedure, public  :: Restart

end type HgInvent_type
  !------------------------------------------------------------------------
contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)
    
    class(HgInvent_type)            :: this
    type(bounds_type) , intent(in) :: bounds

    call this%InitAllocate ( bounds )
    call this%InitHistory ( bounds )
    call this%InitCold ( bounds )

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(HgInvent_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer :: begc, endc, begg, endg
    !------------------------------------------------------------------------
    begg = bounds%begg; endg= bounds%endg
    begc = bounds%begc; endc= bounds%endc

    allocate( this%Hg_emission_grc       (begg:endg))      ; this%Hg_emission_grc         (:) = nan
    allocate( this%Hg_emission_col       (begc:endc))      ; this%Hg_emission_col         (:) = nan
    allocate( this%Hg_emission_rof       (begc:endc))      ; this%Hg_emission_rof         (:) = nan

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use histFileMod    , only : hist_addfld1d
    use shr_infnan_mod , only : nan => shr_infnan_nan, assignment(=)
    !
    ! !ARGUMENTS:
    class(HgInvent_type) :: this
    type(bounds_type), intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer           :: begc,endc,begg,endg
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc

    this%Hg_emission_grc(begg:endg) = spval
    call hist_addfld1d (fname='Hgemission_grc', units='gHg/km^2/yr', &
         avgflag='A', long_name='grid intial Hg emission to CLM', &
         ptr_gcell=this%Hg_emission_grc)

    this%Hg_emission_col(begc:endc) = spval
    call hist_addfld1d (fname='Hgemission_col', units='gHg/km^2/yr', &
         avgflag='A', long_name='column intial Hg emission to CLM', &
         ptr_col=this%Hg_emission_col)
    
    this%Hg_emission_rof(begc:endc) = spval
    call hist_addfld1d (fname='Hgemission_rof', units='ngHg/m^2/s', &
        avgflag='A', long_name='column Hg emission to rof', &
        ptr_col=this%Hg_emission_rof)

  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize cold start conditions for module variables
    !
    ! !USES:
    use shr_kind_mod   , only : r8 => shr_kind_r8
    use clm_varcon     , only : grlnd
    use fileutils      , only : getfil
    !
    ! !ARGUMENTS:
    class(HgInvent_type)         :: this
    type(bounds_type) , intent(in) :: bounds
    !
    ! !LOCAL VARIABLES:
    integer            :: c, g
    logical            :: readvar
    character(len=256) :: locfn

    !-----------------------------------------------------------------------

  end subroutine InitCold

  !------------------------------------------------------------------------
  subroutine Restart(this, bounds, ncid, flag)
   !
    ! !DESCRIPTION: 
    ! Read/write CN restart data for nitrogen state
    !
    ! !USES:
    use shr_infnan_mod       , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
    use clm_time_manager     , only : is_restart, get_nstep
    use restUtilMod
    use ncdio_pio
    use shr_const_mod        , only : SHR_CONST_PDB
    use clm_varcon       , only : denice, denh2o, pondmx, watmin, spval, nameg
    !
    ! !ARGUMENTS:
    class (HgInvent_type) :: this
    type(bounds_type)          , intent(in)    :: bounds 
    type(file_desc_t)          , intent(inout) :: ncid   
    character(len=*)           , intent(in)    :: flag   !'read' or 'write' or 'define'
    !real(r8)                   , intent(in)    :: totvegc_col(bounds%begc:bounds%endc) ! (gC/m2) total vegetation carbon
    !
    ! !LOCAL VARIABLES:
    logical            :: readvar
    !------------------------------------------------------------------------

    
    call restartvar(ncid=ncid, flag=flag, varname='HgEmission', xtype=ncd_double,  &
    dim1name=nameg, long_name='grid intial Hg emission to CLM', units='g/km2/yr', &
    interpinic_flag='interp', readvar=readvar, data=this%Hg_emission_grc) 
    

  end subroutine Restart

end module HgInventoryType
