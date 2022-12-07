module hg0depStreamMod

  !----------------------------------------------------------------------- 
  ! !DESCRIPTION: 
  ! Contains methods for reading in nitrogen deposition data file
  ! Also includes functions for dynamic ndep file handling and 
  ! interpolation.
  !
  ! !USES
  use shr_kind_mod, only: r8 => shr_kind_r8, CL => shr_kind_cl 
  use shr_strdata_mod, only: shr_strdata_type, shr_strdata_create
  use shr_strdata_mod, only: shr_strdata_print, shr_strdata_advance
  use mct_mod     , only: mct_ggrid
  use spmdMod     , only: mpicom, masterproc, comp_id, iam
  use clm_varctl  , only: iulog
  use abortutils  , only: endrun
  use fileutils   , only: getavu, relavu
  use decompMod   , only: bounds_type, ldecomp
  use domainMod   , only: ldomain

  ! !PUBLIC TYPES:
  implicit none
  private
  save

  ! !PUBLIC MEMBER FUNCTIONS:
  public :: hg0dep_init      ! position datasets for dynamic ndep
  public :: hg0dep_interp    ! interpolates between two years of ndep file data
  public :: clm_domain_mct ! Sets up MCT domain for this resolution

  ! !PRIVATE MEMBER FUNCTIONS:
  private :: check_units   ! Check the units and make sure they can be used
  ! ! PRIVATE TYPES
  type(shr_strdata_type)  :: sdat           ! input data stream
  integer :: stream_year_first_hg0dep         ! first year in stream to use
  integer :: stream_year_last_hg0dep          ! last year in stream to use
  integer :: model_year_align_hg0dep          ! align stream_year_firsthg0dep with 
  logical :: divide_by_secs_per_yr = .true. ! divide by the number of seconds per year

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !==============================================================================

contains

  !==============================================================================

  subroutine hg0dep_init(bounds, NLFilename)
   !    
   ! Initialize data stream information.  
   !
   ! Uses:
   use shr_kind_mod     , only : CS => shr_kind_cs
   use clm_varctl       , only : inst_name
   use clm_time_manager , only : get_calendar
   use ncdio_pio        , only : pio_subsystem
   use shr_pio_mod      , only : shr_pio_getiotype
   use shr_nl_mod       , only : shr_nl_find_group_name
   use shr_log_mod      , only : errMsg => shr_log_errMsg
   use shr_mpi_mod      , only : shr_mpi_bcast
   use decompMod        , only : gsmap_lnd_gdc2glo 
   !
   ! arguments
   implicit none
   type(bounds_type), intent(in) :: bounds  
   character(len=*),  intent(in) :: NLFilename   ! Namelist filename
   !
   ! local variables
   integer            :: nu_nml    ! unit for namelist file
   integer            :: nml_error ! namelist i/o error flag
   type(mct_ggrid)    :: dom_clm   ! domain information 
   character(len=CL)  :: stream_fldFileName_hg0dep
   character(len=CL)  :: hg0depmapalgo = 'bilinear'
   character(len=CL)  :: hg0dep_tintalgo = 'linear'
   character(len=CS)  :: hg0dep_taxmode = 'cycle'
   character(len=CL)  :: hg0dep_varlist = 'Hg0dry'
   character(*), parameter :: shr_strdata_unset = 'NOT_SET'
   character(*), parameter :: subName = "('hg0depdyn_init')"
   character(*), parameter :: F00 = "('(hg0depdyn_init) ',4a)"
   !-----------------------------------------------------------------------

   namelist /hg0depdyn_nml/          &
        stream_year_first_hg0dep,    &
        stream_year_last_hg0dep,     &
        model_year_align_hg0dep,     &
        hg0depmapalgo, hg0dep_taxmode, &
        hg0dep_varlist,              &
        stream_fldFileName_hg0dep,   &
        hg0dep_tintalgo

   ! Default values for namelist
    stream_year_first_hg0dep  = 1                ! first year in stream to use
    stream_year_last_hg0dep   = 1                ! last  year in stream to use
    model_year_align_hg0dep   = 1                ! align stream_year_first_hg0dep with this model year
    stream_fldFileName_hg0dep = ' '

   ! Read hg0depdyn_nml namelist
   if (masterproc) then
      nu_nml = getavu()
      open( nu_nml, file=trim(NLFilename), status='old', iostat=nml_error )
      call shr_nl_find_group_name(nu_nml, 'hg0depdyn_nml', status=nml_error)
      if (nml_error == 0) then
         read(nu_nml, nml=hg0depdyn_nml,iostat=nml_error)
         if (nml_error /= 0) then
            call endrun(msg=' ERROR reading hg0depdyn_nml namelist'//errMsg(sourcefile, __LINE__))
         end if
      else
         call endrun(msg=' ERROR finding hg0depdyn_nml namelist'//errMsg(sourcefile, __LINE__))
      end if
      close(nu_nml)
      call relavu( nu_nml )
   endif

   call shr_mpi_bcast(stream_year_first_hg0dep , mpicom)
   call shr_mpi_bcast(stream_year_last_hg0dep  , mpicom)
   call shr_mpi_bcast(model_year_align_hg0dep  , mpicom)
   call shr_mpi_bcast(stream_fldFileName_hg0dep, mpicom)
   call shr_mpi_bcast(hg0dep_varlist           , mpicom)
   call shr_mpi_bcast(hg0dep_taxmode           , mpicom)
   call shr_mpi_bcast(hg0dep_tintalgo          , mpicom)

   if (masterproc) then
      write(iulog,*) ' '
      write(iulog,*) 'hg0depdyn stream settings:'
      write(iulog,*) '  stream_year_first_hg0dep  = ',stream_year_first_hg0dep  
      write(iulog,*) '  stream_year_last_hg0dep   = ',stream_year_last_hg0dep   
      write(iulog,*) '  model_year_align_hg0dep   = ',model_year_align_hg0dep   
      write(iulog,*) '  stream_fldFileName_hg0dep = ',stream_fldFileName_hg0dep
      write(iulog,*) '  hg0dep_varList            = ',hg0dep_varList
      write(iulog,*) '  hg0dep_taxmode            = ',hg0dep_taxmode
      write(iulog,*) '  hg0dep_tintalgo           = ',hg0dep_tintalgo
      write(iulog,*) ' '
   endif
   ! Read in units
   call check_units( stream_fldFileName_hg0dep, hg0dep_varList )

   ! Set domain and create streams
   call clm_domain_mct (bounds, dom_clm)

   call shr_strdata_create(sdat,name="clmhg0dep",    &
        pio_subsystem=pio_subsystem,               & 
        pio_iotype=shr_pio_getiotype(inst_name),   &
        mpicom=mpicom, compid=comp_id,             &
        gsmap=gsmap_lnd_gdc2glo, ggrid=dom_clm,    &
        nxg=ldomain%ni, nyg=ldomain%nj,            &
        yearFirst=stream_year_first_hg0dep,          &
        yearLast=stream_year_last_hg0dep,            &
        yearAlign=model_year_align_hg0dep,           &
        offset=0,                                  &
        domFilePath='',                            &
        domFileName=trim(stream_fldFileName_hg0dep), &
        domTvarName='time',                        &
        domXvarName='lon' ,                        &
        domYvarName='lat' ,                        &  
        domAreaName='area',                        &
        domMaskName='mask',                        &
        filePath='',                               &
        filename=(/trim(stream_fldFileName_hg0dep)/),&
        fldListFile=hg0dep_varlist,                  &
        fldListModel=hg0dep_varlist,                 &
        fillalgo='none',                           &
        mapalgo=hg0depmapalgo,                       &
        tintalgo=hg0dep_tintalgo,                    &
        calendar=get_calendar(),                   &
        taxmode=hg0dep_taxmode                       )


   if (masterproc) then
      call shr_strdata_print(sdat,'CLMhg0dep data')
   endif

 end subroutine hg0dep_init
 !================================================================
  
 subroutine check_units( stream_fldFileName_hg0dep, hg0dep_varList )
   !-------------------------------------------------------------------
   ! Check that units are correct on the file and if need any conversion
   use ncdio_pio     , only : ncd_pio_openfile, ncd_inqvid, ncd_getatt, ncd_pio_closefile, ncd_nowrite
   use ncdio_pio     , only : file_desc_t, var_desc_t
   use shr_kind_mod  , only : CS => shr_kind_cs
   use shr_log_mod   , only : errMsg => shr_log_errMsg
   use shr_string_mod, only : shr_string_listGetName
   implicit none

   !-----------------------------------------------------------------------
   !
   ! Arguments
   character(len=*), intent(IN)  :: stream_fldFileName_hg0dep  ! hg0dep filename
   character(len=*), intent(IN)  :: hg0dep_varList             ! hg0dep variable list to examine
   !
   ! Local variables
   type(file_desc_t) :: ncid     ! NetCDF filehandle for hg0dep file
   type(var_desc_t)  :: vardesc  ! variable descriptor
   integer           :: varid    ! variable index
   logical           :: readvar  ! If variable was read
   character(len=CS) :: hg0depunits! hg0dep units
   character(len=CS) :: fname    ! hg0dep field name
   !-----------------------------------------------------------------------
   call ncd_pio_openfile( ncid, trim(stream_fldFileName_hg0dep), ncd_nowrite )
   call shr_string_listGetName( hg0dep_varList, 1, fname )
   call ncd_inqvid(       ncid, fname, varid, vardesc, readvar=readvar )
   if ( readvar ) then
      call ncd_getatt(    ncid, varid, "units", hg0depunits )
   else
      call endrun(msg=' ERROR finding variable: '//trim(fname)//" in file: "// &
                      trim(stream_fldFileName_hg0dep)//errMsg(sourcefile, __LINE__))
   end if
   call ncd_pio_closefile( ncid )

   ! Now check to make sure they are correct
   if (      trim(hg0depunits) == "g(Hg)/m2/s"  )then
      divide_by_secs_per_yr = .false.
   else if ( trim(hg0depunits) == "g(Hg)/m2/yr" )then
      divide_by_secs_per_yr = .true.
   else
      call endrun(msg=' ERROR in units for hg deposition equal to: '//trim(hg0depunits)//" not units expected"// &
                      errMsg(sourcefile, __LINE__))
   end if

 end subroutine check_units

 !================================================================
 subroutine hg0dep_interp(bounds, atm2lnd_inst)

   !-----------------------------------------------------------------------
   use clm_time_manager, only : get_curr_date, get_days_per_year
   use clm_varcon      , only : secspday
   use atm2lndType     , only : atm2lnd_type
   !
   ! Arguments
   type(bounds_type) , intent(in)    :: bounds  
   type(atm2lnd_type), intent(inout) :: atm2lnd_inst
   !
   ! Local variables
   integer :: g, ig 
   integer :: year    ! year (0, ...) for nstep+1
   integer :: mon     ! month (1, ..., 12) for nstep+1
   integer :: day     ! day of month (1, ..., 31) for nstep+1
   integer :: sec     ! seconds into current date for nstep+1
   integer :: mcdate  ! Current model date (yyyymmdd)
   integer :: dayspyr ! days per year
   !-----------------------------------------------------------------------

   call get_curr_date(year, mon, day, sec)
   mcdate = year*10000 + mon*100 + day

   call shr_strdata_advance(sdat, mcdate, sec, mpicom, 'hg0depdyn')

   if ( divide_by_secs_per_yr )then
      ig = 0
      dayspyr = get_days_per_year( )
      do g = bounds%begg,bounds%endg
         ig = ig+1
         atm2lnd_inst%forc_hg0dep_grc(g) = sdat%avs(1)%rAttr(1,ig) / (secspday * dayspyr)
      end do
   else
      ig = 0
      do g = bounds%begg,bounds%endg
         ig = ig+1
         atm2lnd_inst%forc_hg0dep_grc(g) = sdat%avs(1)%rAttr(1,ig)
         !atm2lnd_inst%forc_hg0dep_grc(g) = sdat%avs(1)%rAttr(1,ig)*2.23596_r8
      end do
   end if
   
 end subroutine hg0dep_interp

 !==============================================================================
  subroutine clm_domain_mct(bounds, dom_clm, nlevels)

    !-------------------------------------------------------------------
    ! Set domain data type for internal clm grid
    use clm_varcon  , only : re
    use domainMod   , only : ldomain
    use seq_flds_mod
    use mct_mod     , only : mct_ggrid, mct_gsMap_lsize, mct_gGrid_init
    use mct_mod     , only : mct_gsMap_orderedPoints, mct_gGrid_importIAttr
    use mct_mod     , only : mct_gGrid_importRAttr
    use mct_mod     , only : mct_gsMap
    use decompMod   , only : gsmap_lnd_gdc2glo, gsMap_lnd2Dsoi_gdc2glo
    implicit none
    ! 
    ! arguments
    type(bounds_type), intent(in) :: bounds  
    type(mct_ggrid), intent(out)   :: dom_clm     ! Output domain information for land model
    integer, intent(in), optional :: nlevels      ! Number of levels if this is a 3D field
    !
    ! local variables
    integer :: g,i,j,k            ! index
    integer :: lsize              ! land model domain data size
    real(r8), pointer :: data(:)  ! temporary
    integer , pointer :: idata(:) ! temporary
    integer :: nlevs              ! Number of vertical levels
    type(mct_gsMap), pointer :: gsmap => null() ! MCT GS map
    !-------------------------------------------------------------------
    ! SEt number of levels, and get the GS map for either the 2D or 3D grid
    nlevs = 1
    if ( present(nlevels) ) nlevs = nlevels
    if ( nlevs == 1 ) then
       gsmap => gsmap_lnd_gdc2glo
    else
       gsmap => gsMap_lnd2Dsoi_gdc2glo
    end if
    !
    ! Initialize mct domain type
    ! lat/lon in degrees,  area in radians^2, mask is 1 (land), 0 (non-land)
    ! Note that in addition land carries around landfrac for the purposes of domain checking
    ! 
    lsize = mct_gsMap_lsize(gsmap, mpicom)
    call mct_gGrid_init( GGrid=dom_clm, CoordChars=trim(seq_flds_dom_coord), &
                         OtherChars=trim(seq_flds_dom_other), lsize=lsize )
    !
    ! Allocate memory
    !
    allocate(data(lsize))
    !
    ! Determine global gridpoint number attribute, GlobGridNum, which is set automatically by MCT
    !
    call mct_gsMap_orderedPoints(gsmap, iam, idata)
    gsmap => null()
    call mct_gGrid_importIAttr(dom_clm,'GlobGridNum',idata,lsize)
    !
    ! Determine domain (numbering scheme is: West to East and South to North to South pole)
    ! Initialize attribute vector with special value
    !
    data(:) = -9999.0_R8 
    call mct_gGrid_importRAttr(dom_clm,"lat"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_clm,"lon"  ,data,lsize) 
    call mct_gGrid_importRAttr(dom_clm,"area" ,data,lsize) 
    call mct_gGrid_importRAttr(dom_clm,"aream",data,lsize) 
    data(:) = 0.0_R8     
    call mct_gGrid_importRAttr(dom_clm,"mask" ,data,lsize) 
    !
    ! Determine bounds
    !
    ! Fill in correct values for domain components
    ! Note aream will be filled in in the atm-lnd mapper
    !
    do k = 1, nlevs
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%lonc(g)
    end do
    end do
    call mct_gGrid_importRattr(dom_clm,"lon",data,lsize) 

    do k = 1, nlevs
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%latc(g)
    end do
    end do
    call mct_gGrid_importRattr(dom_clm,"lat",data,lsize) 

    do k = 1, nlevs
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = ldomain%area(g)/(re*re)
    end do
    end do
    call mct_gGrid_importRattr(dom_clm,"area",data,lsize) 

    do k = 1, nlevs
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = real(ldomain%mask(g), r8)
    end do
    end do
    call mct_gGrid_importRattr(dom_clm,"mask",data,lsize) 

    do k = 1, nlevs
    do g = bounds%begg,bounds%endg
       i = 1 + (g - bounds%begg)
       data(i) = real(ldomain%frac(g), r8)
    end do
    end do
    call mct_gGrid_importRattr(dom_clm,"frac",data,lsize) 

    deallocate(data)
    deallocate(idata)

  end subroutine clm_domain_mct
    
end module hg0depStreamMod

