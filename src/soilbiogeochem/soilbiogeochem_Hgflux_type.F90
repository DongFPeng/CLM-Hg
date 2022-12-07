module SoilBiogeochemHgFluxType

  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_infnan_mod                     , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use clm_varpar                         , only : ndecomp_cascade_transitions, ndecomp_pools
  use clm_varpar                         , only : nlevdecomp_full, nlevdecomp
  use clm_varcon                         , only : spval, ispval, dzsoi_decomp
  use decompMod                          , only : bounds_type
  use clm_varctl                         , only : use_nitrif_denitrif, use_vertsoilc, use_crop
  use CNSharedParamsMod                  , only : use_fun
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use abortutils                         , only : endrun
  use LandunitType                       , only : lun                
  use ColumnType                         , only : col  
  use clm_varpar                         , only : nlevdecomp_full, nlevdecomp, nlevsoi
  use subgridAveMod                      , only : p2c, c2g
  use PhotosynthesisMod                  , only : photosyns_type
  use spmdMod        , only : masterproc
  use clm_varctl     , only : iulog, use_vertsoilc, use_fates
  use decompMod      , only : bounds_type
  use PatchType                          , only : patch
! !PUBLIC TYPES:
  implicit none
  private
  !
  type, public :: SoilBiogeochem_Hgflux_type

     ! deposition fluxes
     real(r8), pointer :: decomp_hgpools_sourcesink_col              (:,:,:) ! change in decomposing c pools. Used to update concentrations concurrently with vertical transport (gC/m3/timestep)  
     
     real(r8), pointer :: decomp_cascade_hr_hg_vr_col               (:,:,:) ! added hg: vertically-resolved het. resp. from decomposing hg pools (gC/m3/s)
     real(r8), pointer :: decomp_cascade_hr_hg_col                     (:,:) ! added hg: vertically-integrated (diagnostic) het. resp. from decomposing hg pools (gC/m2/s)
     real(r8), pointer :: decomp_cascade_hgtransfer_vr_col           (:,:,:) ! added hg:vertically-resolved hg transferred along deomposition cascade (gC/m3/s)
     real(r8), pointer :: decomp_cascade_hgtransfer_col              (:,:)   ! vertically-integrated (diagnostic) hg transferred along decomposition cascade (gC/m2/s)



     real(r8), pointer :: Hg0dep_to_sminn_col                         (:)     ! col atmospheric Hg deposition to soil mineral Hg0 (gN/m2/s)
     real(r8), pointer :: hg2dep_to_sminn_col                        (:)     ! col atmospheric Hg deposition to soil mineral hg2 (gN/m2/s)
     real(r8), pointer :: hg2_dry_dep_to_sminn_col                   (:)
     !leaching fluxes
     real(r8), pointer :: soil_hg2_leached_vr_col                   (:,:)   ! col vertically-resolved soil mineral hg loss to leaching (gN/m3/s)
     real(r8), pointer :: soil_hg2_leached_col                      (:)     ! col soil mineral hg pool loss to leaching (gN/m2/s)
     !fire fluxes
     !real(r8), pointer :: m_decomp_hgpools_to_fire_col               (:,:)   ! vertically-integrated (diagnostic) decomposing hg fire loss (gC/m2/s)
     !real(r8), pointer :: m_decomp_hgpools_to_fire_vr_col            (:,:,:) ! vertically-resolved decomposing hg fire loss (ghg/m3/s)
     real(r8), pointer :: somhg_fire_col                             (:)     ! (gC/m2/s) hg emissions due to peat burning
     

     real(r8), pointer :: rrsun_p2c_col     (:)   !(ghg/m2/s) patch 2 col averaged column-level (p2c)
     real(r8), pointer :: rrsun_tot_col     (:)   ! (ghg/m2/s) total patch-level 

     real(r8), pointer :: hg0_rrsun                   (:)
     real(r8), pointer :: hg_mt_litter                (:)
     real(r8), pointer :: hg0_uptake_veg              (:) 
     real(r8), pointer :: hg0_uptake_out_veg          (:)
     real(r8), pointer :: hg0_uptake_fixed            (:) 
     real(r8), pointer :: hg0_reemission_veg          (:)  
 contains

     procedure , public  :: Init   
     
     procedure , public  :: SetValues
     procedure , public  :: Summary
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory
     procedure , private :: InitCold

 end type SoilBiogeochem_Hgflux_type 
 character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------
 contains

  !------------------------------------------------------------------------
  subroutine Init(this, bounds)

    class(soilbiogeochem_Hgflux_type) :: this
    type(bounds_type), intent(in) :: bounds  

    call this%InitAllocate (bounds)
    call this%InitHistory (bounds)
    call this%InitCold (bounds)

  end subroutine Init    
  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize Hg flux
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_Hgflux_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer           :: begc,endc
    integer           :: begp,endp
    !------------------------------------------------------------------------

    begc = bounds%begc; endc = bounds%endc
    begp = bounds%begp; endp = bounds%endp
    allocate(this%decomp_hgpools_sourcesink_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                  
    this%decomp_hgpools_sourcesink_col(:,:,:)= nan
  
    

     allocate(this%decomp_cascade_hr_hg_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions))        
     this%decomp_cascade_hr_hg_vr_col(:,:,:)= spval !added hg


     allocate(this%decomp_cascade_hr_hg_col(begc:endc,1:ndecomp_cascade_transitions))                             
     this%decomp_cascade_hr_hg_col(:,:)= spval

     allocate(this%decomp_cascade_hgtransfer_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_cascade_transitions)) 
     this%decomp_cascade_hgtransfer_vr_col(:,:,:)= nan

     allocate(this%decomp_cascade_hgtransfer_col(begc:endc,1:ndecomp_cascade_transitions))                      
     this%decomp_cascade_hgtransfer_col(:,:)= nan
    
     !allocate(this%m_decomp_hgpools_to_fire_col(begc:endc,1:ndecomp_pools))                                     
   ! this%m_decomp_hgpools_to_fire_col(:,:)= nan

     !allocate(this%m_decomp_hgpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools))                
     !this%m_decomp_hgpools_to_fire_vr_col(:,:,:)= nan
     allocate(this%somhg_fire_col     (begc:endc))                  ; this%somhg_fire_col     (:)   =nan
     
    allocate(this%Hg0dep_to_sminn_col                (begc:endc))                   ; this%Hg0dep_to_sminn_col         (:)   = nan
    allocate(this%hg2dep_to_sminn_col               (begc:endc))                   ; this%hg2dep_to_sminn_col        (:)   = nan
    allocate(this%hg2_dry_dep_to_sminn_col           (begc:endc))                   ; this%hg2_dry_dep_to_sminn_col        (:)   = nan
    allocate(this%soil_hg2_leached_vr_col           (begc:endc,1:nlevdecomp_full)) ; this%soil_hg2_leached_vr_col    (:,:) = nan
    allocate(this%soil_hg2_leached_col              (begc:endc))                   ; this%soil_hg2_leached_col       (:)   = nan
    allocate(this%rrsun_tot_col                      (begc:endc))                  ; this%rrsun_tot_col              (:) = nan
    allocate(this%rrsun_p2c_col                      (begc:endc))                  ; this%rrsun_p2c_col              (:) = nan

    allocate(this%hg0_rrsun                   (begp:endp)) ; this%hg0_rrsun                   (:) = nan
    allocate(this%hg0_uptake_veg                   (begp:endp)) ; this%hg0_uptake_veg                   (:) = nan
    allocate(this%hg0_uptake_out_veg                   (begp:endp)) ; this%hg0_uptake_out_veg                   (:) = nan
    allocate(this%hg0_uptake_fixed                   (begp:endp)) ; this%hg0_uptake_fixed                   (:) = nan
    allocate(this%hg0_reemission_veg                   (begp:endp)) ; this%hg0_reemission_veg          (:) = nan
    allocate(this%hg_mt_litter                   (begp:endp)) ; this%hg_mt_litter                   (:) = nan
  end subroutine InitAllocate  
  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initialize module data structure
    !
    ! !USES:
    use histFileMod    , only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_Hgflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer        :: k,l,ii,jj,c
    integer        :: begc, endc
    integer        :: begp,endp
    character(24)  :: fieldname
    character(100) :: longname
    character(8)   :: vr_suffix,default
    real(r8), pointer :: data2dptr(:,:), data1dptr(:) ! temp. pointers for slicing larger arrays
    !------------------------------------------------------------------------

    begc = bounds%begc; endc= bounds%endc
    begp = bounds%begp; endp = bounds%endp
    ! add suffix if number of soil decomposition depths is greater than 1
    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif

    !-------------------------------
    ! Hg flux variables - native to column
    !-------------------------------
    !added Hg
    this%Hg0dep_to_sminn_col(begc:endc) = spval
    call hist_addfld1d (fname='Hg0DEP_TO_SMINN', units='ghg/m^2/s', &
         avgflag='A', long_name='atmospheric Hg0 deposition to soil  Hg', &
         ptr_col=this%Hg0dep_to_sminn_col, default='active')



    this%hg2dep_to_sminn_col(begc:endc) = spval
    call hist_addfld1d (fname='hg2DEP_TO_SMINN', units='ghg/m^2/s', &
         avgflag='A', long_name='atmospheric hg2 deposition to soil  Hg', &
         ptr_col=this%hg2dep_to_sminn_col, default='active')
    
    this%hg2_dry_dep_to_sminn_col(begc:endc) = spval
    call hist_addfld1d (fname='hg2_dry_DEP_TO_SMINN', units='ghg/m^2/s', &
         avgflag='A', long_name='atmospheric hg2 dry deposition to soil  Hg', &
         ptr_col=this%hg2_dry_dep_to_sminn_col, default='active')


    if (use_nitrif_denitrif) then
       this%soil_hg2_leached_col(begc:endc) = spval
       call hist_addfld1d (fname='soil_hg2_LEACHED', units='ghg/m^2/s', &
           avgflag='A', long_name='soil hg2 pool loss to leaching', &
           ptr_col=this%soil_hg2_leached_col)
    end if

    if (use_nitrif_denitrif .and.  nlevdecomp_full > 1 ) then 
       this%soil_hg2_leached_vr_col(begc:endc,:) = spval
       call hist_addfld_decomp (fname='soil_hg2_LEACHED'//trim(vr_suffix), units='ghg/m^3/s', type2d='levdcmp', &
            avgflag='A', long_name='soil hg2 pool loss to leaching', &
            ptr_col=this%soil_hg2_leached_vr_col, default='inactive')
    end if


     !--------------added Hg-------------------------------------------------------
     this%decomp_hgpools_sourcesink_col(begc:endc,:,:)= spval 
     do l  = 1, ndecomp_pools
      if ( nlevdecomp_full > 1 ) then
          data2dptr => this%decomp_hgpools_sourcesink_col(:,1:nlevsoi,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'SShg_vr'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' SS_hg (vertically resolved)'
          call hist_addfld2d (fname=fieldname, units='ghg/m^3',  type2d='levsoi', &
              avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       endif
     end do

      this%hg0_rrsun(begp:endp) = spval
       call hist_addfld1d (fname='M_LEAFhg_TO_LITTER', units='g/m^2/s', &
            avgflag='A', long_name='leaf hg mortality', &
            ptr_patch=this%hg0_rrsun, default='inactive')

    this%hg_mt_litter(begp:endp) = spval
       call hist_addfld1d (fname='leafhg_TO_LITTER', units='g/m^2/s', &
            avgflag='A', long_name='leafhg mortality to litter', &
            ptr_patch=this%hg_mt_litter, default='inactive')

    this%hg0_uptake_veg(begp:endp) = spval
       call hist_addfld1d (fname='hg0_uptake_veg', units='g/m^2/s', &
            avgflag='A', long_name='hg uptake via vegetation', &
            ptr_patch=this%hg0_uptake_veg, default='inactive')

    this%hg0_uptake_out_veg(begp:endp) = spval
       call hist_addfld1d (fname='hg0_uptake_out_veg', units='g/m^2/s', &
            avgflag='A', long_name='hg uptake into veg and re_emission ', &
            ptr_patch=this%hg0_uptake_out_veg, default='inactive')
   this%hg0_uptake_fixed(begp:endp) = spval
       call hist_addfld1d (fname='hg0_uptake_fixed', units='g/m^2/s', &
            avgflag='A', long_name='hg uptake and fixed into veg  ', &
            ptr_patch=this%hg0_uptake_fixed, default='inactive')

    this%hg0_reemission_veg(begp:endp) = spval
       call hist_addfld1d (fname='hg0_reemission_veg', units='g/m^2/s', &
            avgflag='A', long_name='hg reemission via vegetation', &
            ptr_patch=this%hg0_reemission_veg, default='inactive')
      !--------------added Hg end ---------------------------------------------  

        this%decomp_cascade_hr_hg_col(begc:endc,:)             = spval
       
        this%decomp_cascade_hr_hg_vr_col(begc:endc,:,:)     = spval !added hg
        this%decomp_cascade_hgtransfer_col(begc:endc,:)      = spval
        this%decomp_cascade_hgtransfer_vr_col(begc:endc,:,:) = spval
        do l = 1, ndecomp_cascade_transitions

           ! output the vertically integrated fluxes only as  default
           !-- HR fluxes (none from CWD)
           if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
              data1dptr => this%decomp_cascade_hr_hg_col(:,l)
              ! check to see if there are multiple pathways that include respiration, and if so, note that in the history file
              ii = 0
              do jj = 1, ndecomp_cascade_transitions
                 if ( decomp_cascade_con%cascade_donor_pool(jj) == decomp_cascade_con%cascade_donor_pool(l) ) ii = ii+1
              end do
              if ( ii == 1 ) then
                 fieldname = &
                      trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'hg_HR'
              else
                 fieldname = &
                      trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'hg_HR_'//&
                      trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_receiver_pool(l)))
              endif
              longname =  'Hg Het. Resp. from '//&
                   trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))
              call hist_addfld1d (fname=fieldname, units='ghg/m^2/s',  &
                   avgflag='A', long_name=longname, &
                   ptr_col=data1dptr, default='inactive')
           endif

           !-- transfer fluxes (none from terminal pool, if present)
           if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
              data1dptr => this%decomp_cascade_hgtransfer_col(:,l)
              fieldname = &
                   trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'hg_TO_'//&
                   trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))//'hg'
              longname =  'decomp. of '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                   ' hg to '//trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' hg'
              call hist_addfld1d (fname=fieldname, units='ghg/m^2/s', &
                   avgflag='A', long_name=longname, &
                   ptr_col=data1dptr, default='inactive')
           endif

          
             
      !----------------------------------------added hg --------------------------------------------------
       ! output the vertically resolved fluxes 
       if ( nlevdecomp_full > 1 ) then  
          if ( .not. decomp_cascade_con%is_cwd(decomp_cascade_con%cascade_donor_pool(l)) ) then
                 data2dptr => this%decomp_cascade_hr_hg_vr_col(:,:,l)
                 ! check to see if there are multiple pathways that include respiration, and if so, note that in the history file
                 ii = 0
                 do jj = 1, ndecomp_cascade_transitions
                    if ( decomp_cascade_con%cascade_donor_pool(jj) == decomp_cascade_con%cascade_donor_pool(l) ) ii = ii+1
                 end do
                 if ( ii == 1 ) then
                    fieldname = &
                         trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))&
                         //'_HR_hg'//trim(vr_suffix)
                 else
                    fieldname = &
                         trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'_HR_hg'//&
                         trim(decomp_cascade_con%decomp_pool_name_short(decomp_cascade_con%cascade_receiver_pool(l)))&
                         //trim(vr_suffix)
                 endif
                 longname =  'Hg Het. Resp. from '//&
                      trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))
                 call hist_addfld_decomp (fname=fieldname, units='gC/m^3/s',  type2d='levdcmp', &
                      avgflag='A', long_name=longname, &
                      ptr_col=data2dptr, default='inactive')
         endif
      !------------------------------------------added end------------------------------------------------
              !-- transfer fluxes (none from terminal pool, if present)
              if ( decomp_cascade_con%cascade_receiver_pool(l) /= 0 ) then
                 data2dptr => this%decomp_cascade_hgtransfer_vr_col(:,:,l)
                 fieldname = &
                      trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_donor_pool(l)))//'Hg_TO_'//&
                      trim(decomp_cascade_con%decomp_pool_name_history(decomp_cascade_con%cascade_receiver_pool(l)))&
                      //'Hg'//trim(vr_suffix)
                 longname =  'decomp. of '//&
                      trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_donor_pool(l)))//&
                      ' Hg to '//&
                      trim(decomp_cascade_con%decomp_pool_name_long(decomp_cascade_con%cascade_receiver_pool(l)))//' Hg'
                 call hist_addfld_decomp (fname=fieldname, units='ghg/m^3/s',  type2d='levdcmp', &
                      avgflag='A', long_name=longname, &
                      ptr_col=data2dptr, default='inactive')
              endif
      end if

        end do

        

     !end if



!------------------------------------------added hg--------------------------------------------------
   !this%m_decomp_hgpools_to_fire_col(begc:endc,:)      = spval
       !this%m_decomp_hgpools_to_fire_vr_col(begc:endc,:,:) = spval
      ! do k = 1, ndecomp_pools
         ! if ( decomp_cascade_con%is_litter(k) .or. decomp_cascade_con%is_cwd(k) ) then
            ! data1dptr => this%m_decomp_hgpools_to_fire_col(:,k)
            ! fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'hg_TO_FIRE'
             !longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' hg fire loss'
            ! call hist_addfld1d (fname=fieldname, units='ghg/m^2/s',  &
                 ! avgflag='A', long_name=longname, &
                 ! ptr_col=data1dptr, default='inactive')

            ! if ( nlevdecomp_full > 1 ) then
               ! data2dptr => this%m_decomp_hgpools_to_fire_vr_col(:,:,k)
               ! fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'hg_TO_FIRE'//trim(vr_suffix)
                !longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' hg fire loss'
               ! call hist_addfld_decomp (fname=fieldname, units='ghg/m^3/s', type2d='levdcmp', &
                     !vgflag='A', long_name=longname, &
                     !ptr_col=data2dptr, default='inactive')
             !endif
          !endif
      ! end do

     this%somhg_fire_col(begc:endc) = spval
        call hist_addfld1d (fname='SOMHG_FIRE', units='ghg/m^2/s', &
             avgflag='A', long_name='hg loss due to peat burning', &
             ptr_col=this%somhg_fire_col)

    
 !------------------------------------------------added end -------------------------------------------------
  end subroutine InitHistory  

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use landunit_varcon , only : istsoil, istcrop
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_Hgflux_type) :: this
    type(bounds_type), intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer :: c,l,p
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: num_special_patch
    integer :: special_col(bounds%endc-bounds%begc+1)    ! special landunit filter - columns
    integer :: special_patch(bounds%endp-bounds%begp+1)  ! special landunit filter - patches
    !---------------------------------------------------------------------

    ! Set column filters

    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

! Set patch filters

    num_special_patch = 0
    do p = bounds%begp,bounds%endp
       l = patch%landunit(p)

       if (lun%ifspecial(l)) then
          num_special_patch = num_special_patch + 1
          special_patch(num_special_patch) = p
       end if
    end do
    !-----------------------------------------------
    ! initialize Hg flux variables
    !-----------------------------------------------

    call this%SetValues (&
         num_patch=num_special_patch, filter_patch=special_patch, value_patch=0._r8, &
         num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, &
       num_patch,filter_patch, value_patch, &
       num_column, filter_column, value_column)
    !
    ! !DESCRIPTION:
    ! Set nitrogen flux variables
    !
    ! !ARGUMENTS:
    ! !ARGUMENTS:
    class(soilbiogeochem_Hgflux_type) :: this
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    integer , intent(in) :: num_patch
    integer , intent(in) :: filter_patch(:)
    real(r8), intent(in) :: value_patch
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i,j,k,l     ! loop index
    !------------------------------------------------------------------------
    
    


    do l = 1, ndecomp_cascade_transitions
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_cascade_hr_hg_col(i,l)             = value_column
            
             this%decomp_cascade_hr_hg_vr_col(i,j,l)     = value_column !added hg
             this%decomp_cascade_hgtransfer_col(i,l)      = value_column
             this%decomp_cascade_hgtransfer_vr_col(i,j,l) = value_column
             !this%decomp_k_col(i,j,l)                    = value_column
          end do
       end do
    end do

    do k = 1, ndecomp_pools
       !do fi = 1,num_column
         ! i = filter_column(fi)
          !this%decomp_cpools_leached_col(i,k) = value_column
       !end do
       do j = 1, nlevdecomp_full
          do fi = 1,num_column
             i = filter_column(fi)
            ! this%decomp_cpools_transport_tendency_col(i,j,k) = value_column
             this%decomp_hgpools_sourcesink_col(i,j,k)         = value_column  
          end do
       end do
    end do


   



    do j=1, nlevdecomp_full

       do fi = 1,num_column
          i = filter_column(fi)
 
       this%Hg0dep_to_sminn_col(i)             = value_column
       this%hg2dep_to_sminn_col(i)            = value_column
       this%hg2_dry_dep_to_sminn_col(i)           = value_column
       this%soil_hg2_leached_vr_col(i,j)       = value_column
       this%soil_hg2_leached_col(i)          = value_column
    end do    
  end do


   !do k = 1, ndecomp_pools
       !do j = 1, nlevdecomp_full
          !do fi = 1,num_column
            ! i = filter_column(fi)
            ! this%m_decomp_hgpools_to_fire_vr_col(i,j,k) = value_column
          !end do
       !end do
    !end do

  !do k = 1, ndecomp_pools
       !do fi = 1,num_column
         ! i = filter_column(fi)
          !this%m_decomp_hgpools_to_fire_col(i,k) = value_column
       !end do
    !end do



    do fi = 1,num_column
       i = filter_column(fi)
       !this%hr_col(i)            = value_column
       this%somhg_fire_col(i)     = value_column
       this%rrsun_p2c_col(i)      = value_column
       this%rrsun_tot_col(i)      = value_column
       
       !this%som_c_leached_col(i) = value_column
       !this%somhr_col(i)         = value_column
       !this%lithr_col(i)         = value_column
       !this%soilc_change_col(i)  = value_column
    end do

     do fi = 1,num_patch
       i = filter_patch(fi)

       this%hg0_rrsun(i)  = value_patch
       this%hg0_uptake_veg(i) = value_patch
       this%hg0_uptake_out_veg(i) = value_patch
       this%hg0_uptake_fixed(i) = value_patch
       this%hg0_reemission_veg(i) = value_patch
       this%hg_mt_litter(i)  = value_patch
    end do 


  end subroutine SetValues
!-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_soilc, filter_soilc,photosyns_inst)
    !
    ! !USES:
    use clm_varpar , only: nlevdecomp, ndecomp_cascade_transitions,ndecomp_pools
    use clm_varctl , only: use_nitrif_denitrif
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_hgflux_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    type(photosyns_type), intent(in):: photosyns_inst
    integer           , intent(in) :: num_soilc       ! number of soil columns in filter
    integer           , intent(in) :: filter_soilc(:) ! filter for soil columns
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,k,l   ! indices
    integer  :: fc        ! filter indices
    !-----------------------------------------------------------------------
     !associate( &
            !fstom                =>photosyns_inst%rssun_patch  )

    

    !call p2c(bounds, num_soilc, filter_soilc, &
         !fstom (bounds%begp:bounds%endp), &
         !this%rrsun_p2c_col(bounds%begc:bounds%endc))

  !do fc = 1,num_soilc
   !   !c = filter_soilc(fc)
    !this%rrsun_tot_col(c) = this%rrsun_p2c_col(c)
 ! end do 
 

    !do fc = 1,num_soilc
       !c = filter_soilc(fc)

    !this%tot_lithg_fire(c)=0._r8
  ! end do 




     
     !do l = 1, ndecomp_pools
         !if ( decomp_cascade_con%is_litter(l) ) then
           ! do fc = 1,num_soilc
               !c = filter_soilc(fc)
              ! this%tot_lithg_fire(c) = &
                   ! this%tot_lithg_fire(c) + &
                    !this%m_decomp_hgpools_to_fire_col(c,l)
          !  end do
       !  end if
     !end do
  

  do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
             this%decomp_cascade_hr_hg_col(c,k) = &
                  this%decomp_cascade_hr_hg_col(c,k) + &
                  this%decomp_cascade_hr_hg_vr_col(c,j,k) * dzsoi_decomp(j) 

             this%decomp_cascade_hgtransfer_col(c,k) = &
                  this%decomp_cascade_hgtransfer_col(c,k) + &
                  this%decomp_cascade_hgtransfer_vr_col(c,j,k) * dzsoi_decomp(j) 
          end do
       end do
    end do




  do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
          this%soil_hg2_leached_col(c) = &
                  this%soil_hg2_leached_col(c) + &
                  this%soil_hg2_leached_vr_col(c,j) * dzsoi_decomp(j)
          end do
  end do  
  

   ! vertically integrate column-level Hg fire losses
   

       ! do l = 1, ndecomp_pools
        
          !do fc = 1,num_soilc
            ! c = filter_soilc(fc)
   
         !this%m_decomp_hgpools_to_fire_col(c,l) = 0._r8
          !end do
     !end do   

    !do l = 1, ndecomp_pools
     ! if ( decomp_cascade_con%is_litter(l) ) then
       !do j = 1,nlevdecomp
         ! do fc = 1,num_soilc
           !  c = filter_soilc(fc)
          	
            ! this%m_decomp_hgpools_to_fire_col(c,l) = &
                !  this%m_decomp_hgpools_to_fire_col(c,l) + &
                !  this%m_decomp_hgpools_to_fire_vr_col(c,j,l)*dzsoi_decomp(j)
         ! end do
       !end do
   !  end if   
   ! end do

 ! end associate

  end subroutine Summary      

end module SoilBiogeochemHgFluxType

