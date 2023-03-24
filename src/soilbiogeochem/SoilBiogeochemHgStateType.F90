module SoilBiogeochemHgStateType

#include "shr_assert.h"

  use shr_kind_mod                       , only : r8 => shr_kind_r8
  use shr_infnan_mod                     , only : nan => shr_infnan_nan, assignment(=)
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  use decompMod                          , only : bounds_type
  use abortutils                         , only : endrun
  use spmdMod                            , only : masterproc 
  use clm_varpar                         , only : ndecomp_cascade_transitions, ndecomp_pools, nlevcan
  use clm_varpar                         , only : nlevdecomp_full, nlevdecomp, nlevsoi
  use clm_varcon                         , only : spval, dzsoi_decomp, zisoi, ispval, zsoi, c3_r2
  use clm_varctl                         , only : use_nitrif_denitrif, use_vertsoilc, use_century_decomp
  use clm_varctl                         , only : iulog, override_bgc_restart_mismatch_dump, spinup_state
  use landunit_varcon                    , only : istcrop, istsoil 
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use LandunitType                       , only : lun                
  use ColumnType                         , only : col                
  use GridcellType                       , only : grc
  use SoilBiogeochemStateType            , only : get_spinup_latitude_term
  ! 
  ! !PUBLIC TYPES:
  implicit none
  private

  type, public :: soilbiogeochem_Hgstate_type
    !Hg leaf out 

     real(r8), pointer :: hg_leafout_pools_col        (:)   ! col (gHg/m2)  leaf out Hg pools
     
     real(r8), pointer :: decomp_hg_rest_vr_col         (:,:) ! col (gN/m3) vertically-resolved decomposing (litter, cwd, soil) N pools
     real(r8), pointer :: decomp_hg_rest_col          (:)   ! col (gN/m3) vertically-resolved decomposing total soil N pool  
     
     real(r8), pointer :: decomp_hgpools_vr_col         (:,:,:) ! col (gN/m3) vertically-resolved decomposing (litter, cwd, soil) Hg pools
     real(r8), pointer :: decomp_soilHg_vr_col           (:,:)   ! col (gN/m3) vertically-resolved decomposing total soil Hg pool

     !wild fire
     real(r8), pointer :: m_decomp_hgpools_to_fire_col               (:,:)   ! vertically-integrated (diagnostic) decomposing hg fire loss (gC/m2/s)
     real(r8), pointer :: m_decomp_hgpools_to_fire_vr_col            (:,:,:) ! vertically-resolved decomposing hg fire loss (ghg/m3/s)
     real(r8), pointer :: tot_lithg_fire                               (:)     ! total litter Hg fire emissions

     ! nitrif_denitrif
     
     real(r8), pointer :: soil_hg2_vr_col              (:,:)  ! col (gHg/m3)vertically-resolved soil hg2
     real(r8), pointer :: soil_hg2_col                 (:)    ! col (gHg/m3)soil hg2 pool

     real(r8), pointer :: soil_hg2_dry_vr_col              (:,:)  ! col (gHg/m3)vertically-resolved soil hg2_dry
     real(r8), pointer :: soil_hg2_dry_col                 (:)    ! col (gHg/m3)soil hg2_dry pool

     real(r8), pointer :: soil_hg_vr_col              (:,:)   ! col (gHg/m3)vertically-resolved soil hg
     real(r8), pointer :: soil_hg_col                 (:)     ! col (gHg/m3)soil hg pool
     real(r8), pointer :: soilHg_to_col                 (:)   !  soil Hg0
     real(r8), pointer :: soilHg_grc                       (:)     
     real(r8), pointer :: initial_soilHg_col            (:,:)   ! col (gHg/m3)vertically-resolved soil hg
     ! summary (diagnostic) state variables, not   involved in mass balance
     !real(r8), pointer :: decomp_npools_col            (:,:)   ! col (gN/m2)  decomposing (litter, cwd, soil) N pools
    !real(r8), pointer :: decomp_npools_1m_col         (:,:)   ! col (gN/m2)  diagnostic: decomposing (litter, cwd, soil) N pools to 1 meter
     !--------------added Hg--------------------------
     
     real(r8), pointer :: decomp_hgpools_col            (:,:)   ! col (gHg/m2)  decomposing (litter, cwd, soil) Hg pools
     real(r8), pointer :: decomp_hgpools_30cm_col            (:,:)
     !---------------added Hg end----------------------  

    
     !----------------------------added Hg----------------------------------
     real(r8), pointer :: totlithg_col                  (:)     ! col (gHg/m2) total litter Hg
     real(r8), pointer :: totsomhg_col                  (:)     ! col (gHg/m2) total soil organic matter Hg
     real(r8), pointer :: totsomhg_30cm_col             (:)
     !----------------------------added Hg end------------------------------
     real(r8), pointer :: totlitn_1m_col               (:)     ! col (gN/m2) total litter nitrogen to 1 meter
     real(r8), pointer :: totsomn_1m_col               (:)     ! col (gN/m2) total soil organic matter nitrogen to 1 meter
     
     ! Track adjustments to no3 and nh4 pools separately, since those aren't included in
     ! the N balance check
     real(r8), pointer :: hg2dry_leaf__patch        (:) 
     real(r8), pointer :: hg2dry_leafwash_col       (:)
      real(r8), pointer :: hg2dry_photored_col     (:)   ! photored part of hg2 dry
      real(r8), pointer :: tot_hg2photored_col     (:)   ! photored budget

      real(r8), pointer :: hg2wet_surf_col(:)            ! finally into the soil surface from hg2wet
      real(r8), pointer :: hg2wet_photored_col(:)        ! photoreducted part form hg2wet
      real(r8), pointer::  hg2wet_heatre_col(:)          ! reemitted of        from hg2wet

      real(r8), pointer :: hg2dry_surf_col(:)            ! finally into the soil surface from hg2dry

      real(r8), pointer :: soilhg_surf_col(:)   
      real(r8), pointer :: soilhg_surf_vr_col(:,:)


     real(r8), pointer :: dyn_hgbal_adjustments_col (:) ! (gHg/m2) hg adjustments to each column made in this timestep via dynamic column adjustments (only makes sense at the column-level: meaningless if averaged to the gridcell-level)

     real(r8)          :: totvegcthresh                  ! threshold for total vegetation carbon to zero out decomposition pools
     integer  :: restart_file_spinup_state             ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
   contains

     procedure , public  :: Init   
     procedure , public  :: Restart
     procedure , public  :: SetValues
     procedure , public  :: Summary
     procedure , public  :: DynamicColumnAdjustments  ! adjust state variables when column areas change
     procedure , public  :: SetTotVgCThresh           ! Set value for totvegcthresh needed in Restart
     procedure , private :: InitAllocate 
     procedure , private :: InitHistory  
     procedure , private :: InitCold     

  end type soilbiogeochem_Hgstate_type

  character(len=*), parameter, private :: sourcefile = &
       __FILE__
  !------------------------------------------------------------------------

contains
  
  !------------------------------------------------------------------------
  subroutine Init(this, bounds,                          &
       soilbiogeochem_hgstate_inst)

    class(soilbiogeochem_hgstate_type)                           :: this
    type(bounds_type)                     , intent(in)           :: bounds  
    !character(len=3)                      , intent(in)           :: carbon_type
    !real(r8)                              , intent(in)           :: ratio
    type(soilbiogeochem_hgstate_type) , intent(in)               :: soilbiogeochem_hgstate_inst
    !real(r8)                              , intent(in)           :: soilHg(bounds%begg)

    this%totvegcthresh = nan
    call this%InitAllocate (bounds )

    call this%InitHistory (bounds)

    call this%InitCold ( bounds, & 
         soilbiogeochem_hgstate_inst)

  end subroutine Init

  !------------------------------------------------------------------------
  subroutine InitAllocate(this, bounds)
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_Hgstate_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    !
    ! !LOCAL VARIABLES:
    integer  :: begc,endc
    integer  :: begg, endg, begp,endp
    real(r8) :: ival  = 0.0_r8
    !------------------------------------------------------------------------
    begg = bounds%begg; endg= bounds%endg
    begc = bounds%begc; endc = bounds%endc
    begp = bounds%begp; endp = bounds%endp
   
    allocate(this%soil_hg2_vr_col     (begc:endc,1:nlevdecomp_full)) ; this%soil_hg2_vr_col     (:,:) = nan
    allocate(this%soil_hg_vr_col      (begc:endc,1:nlevdecomp_full)) ; this%soil_hg_vr_col      (:,:) = nan
    allocate(this%soil_hg2_dry_vr_col (begc:endc,1:nlevdecomp_full)) ; this%soil_hg2_dry_vr_col  (:,:) = nan

    
    allocate(this%soil_hg2_col         (begc:endc))                  ; this%soil_hg2_col        (:)   = nan
    allocate(this%soil_hg_col         (begc:endc))                   ; this%soil_hg_col       (:)   = nan
    allocate(this%soil_hg2_dry_col    (begc:endc))                   ; this%soil_hg2_dry_col   (:)   = nan
    allocate(this%soilHg_grc             (begg:endg))                  ; this%soilHg_grc              (:)   = nan
    !--------------------added Hg------------------------

    allocate(this%totlithg_col          (begc:endc))                   ; this%totlithg_col          (:)   = nan
    allocate(this%tot_lithg_fire         (begc:endc))                   ; this%tot_lithg_fire         (:)   =nan
    allocate(this%totsomhg_col          (begc:endc))                   ; this%totsomhg_col          (:)   = nan
    allocate(this%totsomhg_30cm_col          (begc:endc))              ; this%totsomhg_30cm_col     (:)   = nan
    !----------------------added Hg end------------------

    
    !--------------added Hg-----------------------
    allocate(this%dyn_hgbal_adjustments_col (begc:endc)) ; this%dyn_hgbal_adjustments_col (:) = nan
    !---------------added Hg end------------------

   ! allocate(this%decomp_npools_col    (begc:endc,1:ndecomp_pools))   ; this%decomp_npools_col    (:,:) = nan
    !allocate(this%decomp_npools_1m_col (begc:endc,1:ndecomp_pools))   ; this%decomp_npools_1m_col (:,:) = nan


    !--------------added Hg-----------------------
     allocate(this%decomp_hgpools_col    (begc:endc,1:ndecomp_pools))   ; this%decomp_hgpools_col    (:,:) = nan
     allocate(this%decomp_hgpools_30cm_col    (begc:endc,1:ndecomp_pools))   ; this%decomp_hgpools_30cm_col    (:,:) = nan
     allocate(this%hg_leafout_pools_col    (begc:endc))   ; this%hg_leafout_pools_col    (:) = nan
    !---------------added Hg end------------------

    !allocate(this%decomp_npools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools));
    !this%decomp_npools_vr_col(:,:,:)= nan
   
    allocate(this%decomp_hg_rest_vr_col(begc:endc,1:nlevdecomp_full)) ;this%decomp_hg_rest_vr_col(:,:)= nan
    allocate(this%decomp_hg_rest_col(begc:endc))                      ;this%decomp_hg_rest_col(:)= nan



!------------------------added Hg-------------------------
    allocate(this%decomp_hgpools_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools));
    this%decomp_hgpools_vr_col(:,:,:)= nan
    allocate(this%decomp_soilHg_vr_col(begc:endc,1:nlevdecomp_full));
    this%decomp_soilHg_vr_col(:,:)= nan


    allocate(this%m_decomp_hgpools_to_fire_col(begc:endc,1:ndecomp_pools)); 
    this%m_decomp_hgpools_to_fire_col(:,:)= nan
   
    allocate(this%m_decomp_hgpools_to_fire_vr_col(begc:endc,1:nlevdecomp_full,1:ndecomp_pools)); 
    this%m_decomp_hgpools_to_fire_vr_col(:,:,:)= nan


    allocate(this%hg2dry_leafwash_col (begc:endc))                      ;this%hg2dry_leafwash_col(:)= nan
    
    allocate(this%soilHg_to_col                 (begc:endc))                   ; this%soilHg_to_col  (:)   = nan
    allocate(this%hg2dry_leaf__patch  (begp:endp))                      ;this%hg2dry_leaf__patch(:) = nan
    allocate(this%hg2dry_photored_col (begc:endc))                      ;this%hg2dry_photored_col(:)= nan
    allocate(this%tot_hg2photored_col (begc:endc))                      ;this%tot_hg2photored_col(:)= nan
    allocate(this%hg2wet_surf_col (begc:endc))                          ;this%hg2wet_surf_col(:)= nan
    allocate(this%hg2dry_surf_col (begc:endc))                          ;this%hg2dry_surf_col(:)= nan
    allocate(this%soilhg_surf_col (begc:endc))                          ;this%soilhg_surf_col(:)= nan 
    allocate(this%hg2wet_photored_col (begc:endc))                      ;this%hg2wet_photored_col(:)= nan
    allocate(this%hg2wet_heatre_col (begc:endc))                        ;this%hg2wet_heatre_col(:)= nan
    allocate(this%soilhg_surf_vr_col     (begc:endc,1:nlevdecomp_full)) ; this%soilhg_surf_vr_col     (:,:) = nan
    allocate(this%initial_soilHg_col    (begc:endc,1:nlevdecomp_full)) ; this%initial_soilHg_col      (:,:) = nan
!------------------------added Hg end-------------------------
   this%restart_file_spinup_state = huge(1)

  end subroutine InitAllocate

  !------------------------------------------------------------------------
  subroutine InitHistory(this, bounds)
    !
    ! !DESCRIPTION:
    ! add history fields for all CN variables, always set as default='inactive'
    !
    ! !USES:
    use clm_varpar , only : ndecomp_cascade_transitions, ndecomp_pools
    use clm_varpar , only : nlevdecomp, nlevdecomp_full, nlevgrnd
    use histFileMod, only : hist_addfld1d, hist_addfld2d, hist_addfld_decomp 
    use decompMod  , only : bounds_type
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_Hgstate_type) :: this
    type(bounds_type)         , intent(in) :: bounds 
    !
    ! !LOCAL VARIABLES:
    integer           :: k,l,ii,jj 
    character(10)     :: active
    character(8)      :: vr_suffix
    integer           :: begc,endc,begg,endg,begp,endp
    character(24)     :: fieldname
    character(100)    :: longname
    real(r8), pointer :: data1dptr(:)   ! temp. pointer for slicing larger arrays
    real(r8), pointer :: data2dptr(:,:) ! temp. pointer for slicing larger arrays
    !---------------------------------------------------------------------
    begp = bounds%begp; endp = bounds%endp
    begc = bounds%begc; endc = bounds%endc
    begg = bounds%begg; endg= bounds%endg
    !if ( nlevdecomp_full > 1 ) then
       !this%decomp_soiln_vr_col(begc:endc,:) = spval
       !call hist_addfld2d (fname='SOILN_vr', units='gN/m^3',  type2d='levdcmp', &
            !avgflag='A', long_name='SOIL N (vertically resolved)', &
            !ptr_col=this%decomp_soiln_vr_col)
        !---------------added Hg------------------------
      !this%decomp_soilHg_vr_col(begc:endc,:) = spval
       !call hist_addfld2d (fname='SOILHg_vr', units='gHg/m^3',  type2d='levdcmp', &
            !avgflag='A', long_name='SOIL Hg (vertically resolved)', &
            !ptr_col=this%decomp_soilHg_vr_col)    
        !---------------added Hg end------------------------
    !end if

    
    

    !--------------added Hg-----------------------
    
        do l  = 1, ndecomp_pools
       if ( nlevdecomp_full > 1 ) then
          data2dptr => this%decomp_hgpools_vr_col(:,1:nlevsoi,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'hg_vr'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' hg (vertically resolved)'
          call hist_addfld2d (fname=fieldname, units='ghg/m^3',  type2d='levsoi', &
               avgflag='A', long_name=longname, &
               ptr_col=data2dptr)
       endif
     


    this%decomp_hgpools_col(begc:endc,:) = spval
    data1dptr => this%decomp_hgpools_col(:,l)
          fieldname = trim(decomp_cascade_con%decomp_pool_name_history(l))//'hg'
          longname =  trim(decomp_cascade_con%decomp_pool_name_history(l))//' hg'
          call hist_addfld1d (fname=fieldname, units='ghg/m^2', &
               avgflag='A', long_name=longname, &
               ptr_col=data1dptr)
    !---------------added Hg end------------------

     end do

    

    ! add suffix if number of soil decomposition depths is greater than 1
    if (nlevdecomp > 1) then
       vr_suffix = "_vr"
    else 
       vr_suffix = ""
    endif

    if (use_nitrif_denitrif) then
       if ( nlevdecomp_full > 1 ) then
          

!----------------------------------added Hg----------------------------
          data2dptr => this%soil_hg2_vr_col(begc:endc,1:nlevsoi)
          call hist_addfld_decomp (fname='Soil_hg2'//trim(vr_suffix), units='gHg/m^3',  type2d='levsoi', &
               avgflag='A', long_name='vertically-resolved soil hg2 (vert. res.)', &
               ptr_col=data2dptr)
         
         data2dptr => this%soil_hg2_dry_vr_col(begc:endc,1:nlevsoi)
          call hist_addfld_decomp (fname='Soil_hg2_dry'//trim(vr_suffix), units='gHg/m^3',  type2d='levsoi', &
               avgflag='A', long_name='vertically-resolved soil hg2_dry (vert. res.)', &
               ptr_col=data2dptr)

          data2dptr => this%soil_hg_vr_col(begc:endc,1:nlevsoi)
          call hist_addfld_decomp (fname='Soil_hg'//trim(vr_suffix), units='gHg/m^3',  type2d='levsoi', &
               avgflag='A', long_name='vertically-resolved soil hg (vert. res.)', &
               ptr_col=data2dptr)  
         
       this%m_decomp_hgpools_to_fire_col(begc:endc,:)      = spval
       this%m_decomp_hgpools_to_fire_vr_col(begc:endc,:,:) = spval
       do k = 1, ndecomp_pools
          if ( decomp_cascade_con%is_litter(k) .or. decomp_cascade_con%is_cwd(k) ) then
             data1dptr => this%m_decomp_hgpools_to_fire_col(:,k)
             fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'hg_TO_FIRE'
             longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' hg fire loss'
             call hist_addfld1d (fname=fieldname, units='ghg/m^2/s',  &
                  avgflag='A', long_name=longname, &
                  ptr_col=data1dptr, default='inactive')

             if ( nlevdecomp_full > 1 ) then
                data2dptr => this%m_decomp_hgpools_to_fire_vr_col(:,:,k)
                fieldname = 'M_'//trim(decomp_cascade_con%decomp_pool_name_history(k))//'hg_TO_FIRE'//trim(vr_suffix)
                longname =  trim(decomp_cascade_con%decomp_pool_name_long(k))//' hg fire loss'
                call hist_addfld_decomp (fname=fieldname, units='ghg/m^3/s', type2d='levdcmp', &
                     avgflag='A', long_name=longname, &
                     ptr_col=data2dptr, default='inactive')
             endif
          endif
       end do    
!----------------------------------added Hg end----------------------------
         data2dptr => this%decomp_hg_rest_vr_col(begc:endc,1:nlevsoi)
          call hist_addfld_decomp (fname='decomp_hg_rest'//trim(vr_suffix), units='gHg/m^3',  type2d='levsoi', &
               avgflag='A', long_name='vertically-resolved remained hg (vert. res.)', &
               ptr_col=data2dptr)

         
         this%decomp_hg_rest_col(begc:endc) = spval
          call hist_addfld1d (fname='decomp_hg_pool', units='gHg/m^2', &
               avgflag='A', long_name='total remain decomposed hg', &
               ptr_col=this%decomp_hg_rest_col)


          

          !----------------------------------added Hg----------------------------

          this%soil_hg2_col(begc:endc) = spval
          call hist_addfld1d (fname='Soil_hg2_pool', units='gHg/m^2', &
               avgflag='A', long_name='soil hg2 pool', &
               ptr_col=this%soil_hg2_col)

          this%soil_hg2_dry_col(begc:endc) = spval
          call hist_addfld1d (fname='Soil_hg2_dry_pool', units='gHg/m^2', &
               avgflag='A', long_name='soil hg2 dry pool', &
               ptr_col=this%soil_hg2_col)

          this%soil_hg_col(begc:endc) = spval
          call hist_addfld1d (fname='Soil_hg_pool', units='gHg/m^2', &
               avgflag='A', long_name='soil hg pool', &
               ptr_col=this%soil_hg_col)
       
          this%hg_leafout_pools_col(begc:endc) = spval
          call hist_addfld1d (fname='hg_leafout_pool', units='gHg/m^2', &
               avgflag='A', long_name='hg from leafout pool', &
               ptr_col=this%hg_leafout_pools_col)

          this%soilHg_grc(begg:endg) = spval
          call hist_addfld1d (fname='soilHg', units='gHg/m^2', &
               avgflag='A', long_name='initial soilHg pool', &
               ptr_gcell=this%soilHg_grc)
          !----------------------------------added Hg end----------------------------

       end if
    
       

    end if

    

    !-----------added Hg ---------------------------------------

     this%totlithg_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTLIThg', units='gHg/m^2', &
         avgflag='A', long_name='total litter Hg', &
         ptr_col=this%totlithg_col)

    this%totsomhg_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTSOMhg', units='gHg/m^2', &
         avgflag='A', long_name='total soil organic matter Hg', &
         ptr_col=this%totsomhg_col)

    this%totsomhg_30cm_col(begc:endc) = spval
    call hist_addfld1d (fname='TOTSOMhg_30cm', units='gHg/m^2', &
         avgflag='A', long_name='total soil organic matter Hg in depth of 30cm', &
         ptr_col=this%totsomhg_30cm_col)

    this%tot_lithg_fire(begc:endc) = spval
        call hist_addfld1d (fname='tot_lithg_fire', units='ghg/m^2/s', &
             avgflag='A', long_name='hg loss due to litter burning', &
             ptr_col=this%tot_lithg_fire)



   this%soilHg_to_col(begc:endc) = spval
    call hist_addfld1d (fname='soilHg to column', units='ghg/m2', &
         avgflag='A', long_name='initail soil Hg concentration in column', &
         ptr_col=this%soilHg_to_col)
    
   this%hg2dry_leaf__patch(begp:endp) = spval
       call hist_addfld1d (fname='hg2dry_leaf_p', units='ghg/m^2', &
            avgflag='A', long_name='hg2 dry into leaf', &
            ptr_patch=this%hg2dry_leaf__patch)

   this%hg2dry_leafwash_col(begc:endc) = spval
       call hist_addfld1d (fname='hg2dry_leafwash_c', units='ghg/m^2', &
            avgflag='A', long_name='hg2 dry washed off from leaf', &
            ptr_col=this%hg2dry_leafwash_col)

   this%hg2dry_photored_col(begc:endc) = spval
       call hist_addfld1d (fname='hg2dry_photored_c', units='ghg/m^2', &
            avgflag='A', long_name='hg2 dry photoreducted from leaf', &
            ptr_col=this%hg2dry_photored_col)
   
   this%hg2wet_photored_col(begc:endc) = spval
       call hist_addfld1d (fname='hg2wet_photored_c', units='ghg/m^2', &
            avgflag='A', long_name='hg2 wet photoreducted from leaf', &
            ptr_col=this%hg2wet_photored_col)

   this%hg2wet_heatre_col(begc:endc) = spval
       call hist_addfld1d (fname='hg2wet_heatre_c', units='ghg/m^2', &
            avgflag='A', long_name='hg2 wet heated from leaf', &
            ptr_col=this%hg2wet_heatre_col)

   this%hg2wet_surf_col(begc:endc) = spval
       call hist_addfld1d (fname='hg2wet_surf_c', units='ghg/m^2', &
            avgflag='A', long_name='hg2wet into soil surface', &
            ptr_col=this%hg2wet_surf_col)

   this%tot_hg2photored_col(begc:endc) = spval
       call hist_addfld1d (fname='tot_hg2photored_c', units='ghg/m^2', &
            avgflag='A', long_name='total hg2  photoreducted', &
            ptr_col=this%tot_hg2photored_col)

   this%hg2dry_surf_col(begc:endc) = spval
       call hist_addfld1d (fname='hg2dry_surf_c', units='ghg/m^2', &
            avgflag='A', long_name='hg2dry into soil surface', &
            ptr_col=this%hg2dry_surf_col)

   this%soilhg_surf_col(begc:endc) = spval
       call hist_addfld1d (fname='soilhg_surf_c', units='ghg/m^2', &
            avgflag='A', long_name='all hg into soil surface', &
            ptr_col=this%soilhg_surf_col)

   data2dptr => this%soilhg_surf_vr_col(begc:endc,1:nlevsoi)
          call hist_addfld_decomp (fname='soilhg_fromsurtodeep'//trim(vr_suffix), units='gHg/m^3',  type2d='levsoi', &
               avgflag='A', long_name='vertically-resolved soil hg (vert. res.)', &
               ptr_col=data2dptr)

    !added Hg end---------------------------------------------



      if (use_nitrif_denitrif) then

      call hist_addfld1d (fname='DYN_COL_SOIL_ADJUSTMENTS_hg', units='gHg/m^2', &
            avgflag='SUM', &
            long_name='Adjustments in soil Hg due to dynamic column areas; &
            &only makes sense at the column level: should not be averaged to gridcell', &
            ptr_col=this%dyn_hgbal_adjustments_col, default='inactive')


    end if
  end subroutine InitHistory

  !-----------------------------------------------------------------------
  subroutine InitCold(this, bounds, &
      soilbiogeochem_hgstate_inst)
    !
    ! !DESCRIPTION:
    ! Initializes time varying variables used only in coupled carbon-nitrogen mode (CN):
    !
    ! !USES:
    use decompMod      , only : bounds_type
    !
    ! !ARGUMENTS:
    class(soilbiogeochem_Hgstate_type)            :: this
    type(bounds_type) , intent(in)                :: bounds  
    type(soilbiogeochem_hgstate_type), intent(in) :: soilbiogeochem_hgstate_inst
 
    !
    ! !LOCAL VARIABLES:
    integer :: fc,g,l,c,j,k,n                              ! indices
    integer :: num_special_col                           ! number of good values in special_col filter
    integer :: special_col   (bounds%endc-bounds%begc+1) ! special landunit filter - columns
    !------------------------------------------------------------------------
  
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
      if (lun%itype(l) == istsoil .or. lun%itype(l) == istcrop) then
         do j = 1, nlevdecomp
                do k = 1, ndecomp_pools
                      this%decomp_hgpools_vr_col(c,j,k) = 0._r8
                      this%m_decomp_hgpools_to_fire_vr_col(c,j,k) = 0._r8
                end do
         end do      
           if ( nlevdecomp > 1 ) then
                do j = nlevdecomp+1, nlevdecomp_full
                   do k = 1, ndecomp_pools
                      this%decomp_hgpools_vr_col(c,j,k) = 0._r8
                      this%m_decomp_hgpools_to_fire_vr_col(c,j,k) = 0._r8
                   end do
              end do
             end if

            this%decomp_hgpools_col(c,1:ndecomp_pools)    = 0._r8
            this%decomp_hgpools_30cm_col(c,1:ndecomp_pools)    = 0._r8
            this%m_decomp_hgpools_to_fire_col(c,1:ndecomp_pools) = 0._r8 
      end if 

          if (use_nitrif_denitrif) then
             do j = 1, nlevdecomp_full
               this%soil_hg_vr_col(c,j) = 0._r8
                this%soil_hg2_vr_col(c,j) = 0._r8
                this%soil_hg2_dry_vr_col(c,j) = 0._r8
                this%soilhg_surf_vr_col(c,j) = 0._r8
                this%decomp_hg_rest_vr_col= 0._r8
                this%initial_soilHg_col(c,j) = 0._r8
             end do
               this%soil_hg_col(c) = 0._r8
                this%soil_hg2_col(c) = 0._r8
                this%soil_hg2_dry_col(c) = 0._r8
          end if
               this%totlithg_col(c)    = 0._r8
               this%totsomhg_col(c)    = 0._r8
               this%totsomhg_30cm_col(c)    = 0._r8
               this%hg_leafout_pools_col(c) = 0._r8
               this%tot_lithg_fire(c) = 0._r8
     end do


    ! initialize fields for special filters
    num_special_col = 0
    do c = bounds%begc, bounds%endc
       l = col%landunit(c)
       if (lun%ifspecial(l)) then
          num_special_col = num_special_col + 1
          special_col(num_special_col) = c
       end if
    end do

    call this%SetValues (num_column=num_special_col, filter_column=special_col, value_column=0._r8)

  end subroutine InitCold

  !-----------------------------------------------------------------------
  subroutine Restart ( this,  bounds, ncid, flag )
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
    class (soilbiogeochem_Hgstate_type) :: this
    type(bounds_type)          , intent(in)    :: bounds 
    type(file_desc_t)          , intent(inout) :: ncid   
    character(len=*)           , intent(in)    :: flag   !'read' or 'write' or 'define'
    !real(r8)                   , intent(in)    :: totvegc_col(bounds%begc:bounds%endc) ! (gC/m2) total vegetation carbon
    !type(soilbiogeochem_Hgstate_type), intent(in) :: soilbiogeochem_hgrestart_inst
    !
    ! !LOCAL VARIABLES:
    integer            :: i,j,k,l,c
    logical            :: readvar
    integer            :: idata
    logical            :: exit_spinup = .false.
    logical            :: enter_spinup = .false.
    real(r8)           :: m          ! multiplier for the exit_spinup code
    real(r8), pointer  :: ptr2d(:,:) ! temp. pointers for slicing larger arrays
    real(r8), pointer  :: ptr1d(:)   ! temp. pointers for slicing larger arrays
    character(len=128) :: varname    ! temporary
    integer            :: itemp      ! temporary 
    integer , pointer  :: iptemp(:)  ! pointer to memory to be allocated
    ! spinup state as read from restart file, for determining whether to enter or exit spinup mode.
     !integer            :: restart_file_spinup_state 
    ! flags for comparing the model and restart decomposition cascades
    integer            :: decomp_cascade_state, restart_file_decomp_cascade_state 
    !------------------------------------------------------------------------

    

  !-------------------------------------------added hg-----------------------------------------------
   call restartvar(ncid=ncid, flag=flag, varname='soilHg', xtype=ncd_double,  &
         dim1name=nameg, long_name='', units='', &
         interpinic_flag='interp', readvar=readvar, data=this%soilHg_grc) 



!    do k = 1, ndecomp_pools
!           varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'hg'
!           if (use_vertsoilc) then
!              ptr2d => this%decomp_hgpools_vr_col(:,:,k)
!              call restartvar(ncid=ncid, flag=flag, varname=trim(varname)//"_vr", xtype=ncd_double,  &
!                   dim1name='column', dim2name='levgrnd', switchdim=.true., &
!                   long_name='',  units='', fill_value=spval, &
!                   interpinic_flag='interp', readvar=readvar, data=ptr2d)
!           else
!              ptr1d => this%decomp_hgpools_vr_col(:,1,k) ! nlevdecomp = 1; so treat as 1D variable
!              call restartvar(ncid=ncid, flag=flag, varname=varname, xtype=ncd_double,  &
!                   dim1name='column', long_name='',  units='', fill_value=spval, &
!                   interpinic_flag='interp' , readvar=readvar, data=ptr1d)
!           end if
!     end do 

! ! !------------------------------------------added end --------------------------------------------------------


! !     ! decomposing N pools
    

    
!     !--------------------------------
!     ! Spinup state
!     !--------------------------------

!     ! Do nothing for write
!     ! Note that the call to write spinup_state out was done in soilbiogeochem_carbonstate_inst and
!     ! cannot be called again because it will try to define the variable twice
!     ! when the flag below is set to define
!       if (flag == 'write') idata = spinup_state
!            call restartvar(ncid=ncid, flag=flag, varname='spinup_state', xtype=ncd_int,  &
!              long_name='Spinup state of the model that wrote this restart file: ' &
!              // ' 0 = normal model mode, 1 = AD spinup', units='', &
!              interpinic_flag='copy', readvar=readvar,  data=idata)
!            if (flag == 'read') then
!               if (readvar) then
!                  this%restart_file_spinup_state = idata
!               else
!                  call endrun(msg=' CN-hgRest: spinup_state was not on the restart file and is required' // &
!                    errMsg(sourcefile, __LINE__))
!               end if
!            end if
!         else
!           this%restart_file_spinup_state = soilbiogeochem_hgstate_inst%restart_file_spinup_state
!      end if

!         ! now compare the model and restart file spinup states, and either take the 
!         ! model into spinup mode or out of it if they are not identical
!         ! taking model out of spinup mode requires multiplying each decomposing pool 
!         ! by the associated AD factor.
!         ! putting model into spinup mode requires dividing each decomposing pool 
!         ! by the associated AD factor.
!         ! only allow this to occur on first timestep of model run.
        
!         if (flag == 'read' .and. spinup_state /= this%restart_file_spinup_state ) then
!            if (spinup_state == 0 .and. this%restart_file_spinup_state >= 1 ) then
!               if ( masterproc ) write(iulog,*) ' CNRest: taking ','hg spin',' SOM pools out of AD spinup mode'
!               exit_spinup = .true.
!            else if (spinup_state >= 1 .and. this%restart_file_spinup_state == 0 ) then
!               if ( masterproc ) write(iulog,*) ' CNRest: taking ','hg spin',' SOM pools into AD spinup mode'
!               enter_spinup = .true.
!            else
!               call endrun(msg=' CNRest: error in entering/exiting spinup.  spinup_state ' &
!                    // ' != restart_file_spinup_state, but do not know what to do'//&
!                    errMsg(sourcefile, __LINE__))
!            end if
!            if (get_nstep() >= 2) then
!               call endrun(msg=' CNRest: error in entering/exiting spinup - should occur only when nstep = 1'//&
!                    errMsg(sourcefile, __LINE__))
!            endif
!            if ( exit_spinup .and. isnan(this%totvegcthresh) )then
!               call endrun(msg=' CNRest: error in exit spinup - totvegcthresh was not set with SetTotVgCThresh'//&
!                    errMsg(sourcefile, __LINE__))
!            end if
!            do k = 1, ndecomp_pools
!               if ( exit_spinup ) then
!                  m = decomp_cascade_con%spinup_factor(k)
!               else if ( enter_spinup ) then
!                  m = 1. / decomp_cascade_con%spinup_factor(k)
!               end if
!               do c = bounds%begc, bounds%endc
!                  l = col%landunit(c)
!                  do j = 1, nlevdecomp_full
!                     if ( abs(m - 1._r8) .gt. 0.000001_r8 .and. exit_spinup) then
!                        this%decomp_hgpools_vr_col(c,j,k) = this%decomp_hgpools_vr_col(c,j,k) * m * &
!                             get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
!                        ! If there is no vegetation carbon, implying that all vegetation has died, then 
!                        ! reset decomp pools to near zero during exit_spinup to avoid very 
!                        ! large and inert soil carbon stocks; note that only pools with spinup factor > 1 
!                        ! will be affected, which means that total SOMC and LITC pools will not be set to 0.
!                        !if (totvegc_col(c) <= this%totvegcthresh .and. lun%itype(l) /= istcrop) then 
!                          !this%decomp_hgpools_vr_col(c,j,k) = 0.0_r8
!                        !endif
!                     elseif ( abs(m - 1._r8) .gt. 0.000001_r8 .and. enter_spinup) then
!                        this%decomp_hgpools_vr_col(c,j,k) = this%decomp_hgpools_vr_col(c,j,k) * m / &
!                             get_spinup_latitude_term(grc%latdeg(col%gridcell(c)))
!                     else
!                        this%decomp_hgpools_vr_col(c,j,k) = this%decomp_hgpools_vr_col(c,j,k) * m 
!                     endif
!                  end do
!               end do
!            end do
!         end if

  end subroutine Restart

  !-----------------------------------------------------------------------
  subroutine SetValues ( this, num_column, filter_column, value_column )
    !
    ! !DESCRIPTION:
    ! Set hg state variables
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_Hgstate_type) :: this
    integer , intent(in) :: num_column
    integer , intent(in) :: filter_column(:)
    real(r8), intent(in) :: value_column
    !
    ! !LOCAL VARIABLES:
    integer :: fi,i     ! loop index
    integer :: j,k      ! indices
    !------------------------------------------------------------------------

    do fi = 1,num_column
       i = filter_column(fi)

       
          !--------------added Hg-----------------------
            this%soil_hg2_col(i) = value_column
            this%soil_hg2_dry_col(i) = value_column
            this%soil_hg_col(i) = value_column
            !this%soilHg(i) = value_column
          !---------------added Hg end------------------
   
      ! end if
             
       !--------------added Hg-----------------------
         this%totlithg_col(i)     = value_column
         this%totsomhg_col(i)     = value_column
         this%totsomhg_30cm_col(i)     = value_column
         this%hg_leafout_pools_col(i) = value_column
         this%tot_lithg_fire(i)     = value_column
       !---------------added Hg end------------------

      
    end do

    do j = 1,nlevdecomp_full
       do fi = 1,num_column
          i = filter_column(fi)
          
          if (use_nitrif_denitrif) then
             
             !--------------added Hg-----------------------
            this%soil_hg2_vr_col(i,j) = value_column
            this%soilhg_surf_vr_col(i,j) = value_column
            this%soil_hg2_dry_vr_col(i,j) = value_column
            this%soil_hg_vr_col(i,j) = value_column
            this%decomp_Hg_rest_vr_col(i,j) = value_column
            this%initial_soilHg_col(i,j) = value_column
            
           end if
         end do
        end do
   
      do k = 1, ndecomp_pools
       do fi = 1,num_column
          i = filter_column(fi)
          this%decomp_hgpools_col(i,k) = value_column
          this%decomp_hgpools_30cm_col(i,k) = value_column
          this%m_decomp_hgpools_to_fire_col(i,k) = value_column
       end do
    end do


  do k = 1, ndecomp_pools

    do j = 1,nlevdecomp_full
      
          do fi = 1,num_column
             i = filter_column(fi)
             this%decomp_hgpools_vr_col(i,j,k) = value_column
            
          end do
       end do
    end do


 do k = 1, ndecomp_pools
     if ( decomp_cascade_con%is_litter(k) ) then
    do j = 1,nlevdecomp_full
      
          do fi = 1,num_column
             i = filter_column(fi)
             
             this%m_decomp_hgpools_to_fire_vr_col(i,j,k) = value_column 
          end do
       end do
     end if
    end do




          !---------------added Hg end------------------
             

          
       

    ! column and decomp_pools
    

    


  end subroutine SetValues

  !-----------------------------------------------------------------------
  subroutine Summary(this, bounds, num_allc, filter_allc)
    !
    ! !ARGUMENTS:
    class (soilbiogeochem_Hgstate_type) :: this
    type(bounds_type) , intent(in) :: bounds  
    integer           , intent(in) :: num_allc       ! number of columns in allc filter
    integer           , intent(in) :: filter_allc(:) ! filter for all active columns
    !
    ! !LOCAL VARIABLES:
    integer  :: c,j,k,l     ! indices
    integer  :: fc          ! lake filter indices
    real(r8) :: maxdepth    ! depth to integrate soil variables
    !-----------------------------------------------------------------------

   ! vertically integrate NO3 NH4 N2O pools
   

  
     !--------------added Hg-----------------------
    if (use_nitrif_denitrif) then
      do fc = 1,num_allc
         c = filter_allc(fc)
                  this%soil_hg_col(c) = 0._r8
                  this%soil_hg2_col(c) = 0._r8
                  this%soil_hg2_dry_col(c) = 0._r8
                  this%decomp_hg_rest_col(c) = 0._r8
      end do
       do j = 1, nlevdecomp
       
         do fc = 1,num_allc
            c = filter_allc(fc)
            
                 this%soil_hg_col(c) =&
                 this%soil_hg_col(c) + &
                 this%soil_hg_vr_col(c,j) * dzsoi_decomp(j) !土壤厚度


                 this%soil_hg2_col(c) =&
                 this%soil_hg2_col(c) + &
                 this%soil_hg2_vr_col(c,j) * dzsoi_decomp(j) !土壤厚度
                 
                 this%soil_hg2_dry_col(c) =&
                 this%soil_hg2_dry_col(c) + &
                 this%soil_hg2_dry_vr_col(c,j) * dzsoi_decomp(j) !土壤厚度

                 this%decomp_hg_rest_col(c) =&
                 this%decomp_hg_rest_col(c) + &
                 this%decomp_hg_rest_vr_col(c,j) * dzsoi_decomp(j) !土壤厚度

                 

               

          end do 
       end do

      end if

 ! vertically integrate each of the decomposing hg pools
    do l = 1, ndecomp_pools
       do fc = 1,num_allc
          c = filter_allc(fc)
          this%decomp_hgpools_col(c,l) = 0._r8
          this%decomp_hgpools_30cm_col(c,l) = 0._r8
       end do
    end do

    do l = 1, ndecomp_pools
       do j = 1, nlevdecomp
          do fc = 1,num_allc
             c = filter_allc(fc)
             this%decomp_hgpools_col(c,l) = &
                  this%decomp_hgpools_col(c,l) + &
                  this%decomp_hgpools_vr_col(c,j,l) * dzsoi_decomp(j)

             
          end do
       end do
    end do


 do l = 1, ndecomp_pools
       do j = 1, 5
          do fc = 1,num_allc
             c = filter_allc(fc)
             this%decomp_hgpools_30cm_col(c,l) = &
                  this%decomp_hgpools_30cm_col(c,l) + &
                  this%decomp_hgpools_vr_col(c,j,l) * dzsoi_decomp(j)

             
          end do
       end do
    end do



do l = 1, ndecomp_pools
       do fc = 1,num_allc
          c = filter_allc(fc)
          this%m_decomp_hgpools_to_fire_col(c,l) = 0._r8
       end do
    end do

      do l = 1, ndecomp_pools
        ! if ( decomp_cascade_con%is_litter(l) ) then
       do j = 1, nlevdecomp
          do fc = 1,num_allc
             c = filter_allc(fc)
 

               this%m_decomp_hgpools_to_fire_col(c,l) = &
               this%m_decomp_hgpools_to_fire_col(c,l) + &
               this%m_decomp_hgpools_to_fire_vr_col(c,j,l) * dzsoi_decomp(j)
    
            
        end do
       end do
     !end if
    end do
 


      !do l = 1, ndecomp_pools
       !do j = 1, nlevdecomp
          !do fc = 1,num_allc
           !  c = filter_allc(fc)
            ! this%decomp_hgpools_col(c,l) = &
                 ! this%decomp_hgpools_col(c,l) + &
                 ! this%decomp_hgpools_vr_col(c,j,l) * dzsoi_decomp(j)
         ! end do
      ! end do
    !end do
     !---------------added Hg end------------------


   ! vertically integrate each of the decomposing N pools
  


   !--------------added Hg-----------------------
   ! vertically integrate each of the decomposing Hg pools
   !do l = 1, ndecomp_pools
     ! do fc = 1,num_allc
        ! c = filter_allc(fc)
         !this%decomp_Hgpools_col(c,l) = 0._r8
      !end do
     ! do j = 1, nlevdecomp
         !do fc = 1,num_allc
           ! c = filter_allc(fc)
            !this%decomp_Hgpools_col(c,l) = &
                ! this%decomp_Hgpools_col(c,l) + &
                 !this%decomp_hgpools_vr_col(c,j,l) * dzsoi_decomp(j)
         !end do
      !end do
  ! end do
   !---------------added Hg end------------------


   ! for vertically-resolved soil biogeochemistry, calculate some diagnostics of carbon pools to a given depth
  
  !----------------------summary the total decompsed remian soil hg----------------------

   



      !--------------added Hg-----------------------
       if ( nlevdecomp_full > 1 ) then
         do j = 1, nlevdecomp
            do fc = 1,num_allc
               c = filter_allc(fc)
               this%decomp_soilHg_vr_col(c,j) = 0._r8
            end do
         end do
         do l = 1, ndecomp_pools
            if ( decomp_cascade_con%is_soil(l) ) then
               do j = 1, nlevdecomp
                  do fc = 1,num_allc
                     c = filter_allc(fc)
                     this%decomp_soilHg_vr_col(c,j) = this%decomp_soilHg_vr_col(c,j) + &
                          this%decomp_hgpools_vr_col(c,j,l)
                  end do
               end do
            end if
         end do
      end if
      !---------------added Hg end------------------


     
      ! total soil organic matter nitrogen to 1 meter (TOTSOMN_1m)
      
   
    

  

  !--------------added Hg-----------------------
   ! total litter Hg (TOTLITHg)
   do fc = 1,num_allc
      c = filter_allc(fc)
      this%totlithg_col(c)    = 0._r8
   end do
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_litter(l) ) then
       do fc = 1,num_allc
            c = filter_allc(fc)
            this%totlithg_col(c) = &
                 this%totlithg_col(c) + &
                 this%decomp_hgpools_col(c,l)
         end do
     end if
   end do
   
  ! total litter fire hg 
   do fc = 1,num_allc
      c = filter_allc(fc)
      this%tot_lithg_fire(c)    = 0._r8
   end do

   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_litter(l) ) then
       do fc = 1,num_allc
            c = filter_allc(fc)
            this%tot_lithg_fire(c) = &
                 this%tot_lithg_fire(c) + &
                 this%m_decomp_hgpools_to_fire_col(c,l)
        end do
     end if
   end do


   ! total soil organic matter Hg (TOTSOMHg)
   do fc = 1,num_allc
      c = filter_allc(fc)
      this%totsomhg_col(c)    = 0._r8
      this%totsomhg_30cm_col(c)    = 0._r8
   end do

   
   do l = 1, ndecomp_pools
      if ( decomp_cascade_con%is_soil(l) ) then
         do fc = 1,num_allc
            c = filter_allc(fc)
            this%totsomhg_col(c) = this%totsomhg_col(c) + &
               this%decomp_hgpools_col(c,l)

            this%totsomhg_30cm_col(c) = this%totsomhg_30cm_col(c) + &
               this%decomp_hgpools_30cm_col(c,l)


         end do
      end if
   end do
   
   
  !---------------added Hg end------------------ 


 end subroutine Summary

 !-----------------------------------------------------------------------
 subroutine DynamicColumnAdjustments(this, bounds, clump_index, column_state_updater)
   !
   ! !DESCRIPTION:
   ! Adjust state variables when column areas change due to dynamic landuse
   !
   ! !USES:
   use dynColumnStateUpdaterMod, only : column_state_updater_type
   !
   ! !ARGUMENTS:
   class(soilbiogeochem_Hgstate_type) , intent(inout) :: this
   type(bounds_type)                        , intent(in)    :: bounds

   ! Index of clump on which we're currently operating. Note that this implies that this
   ! routine must be called from within a clump loop.
   integer                                  , intent(in)    :: clump_index

   type(column_state_updater_type)          , intent(in)    :: column_state_updater
   !
   ! !LOCAL VARIABLES:
   integer :: j  ! level
   integer :: l  ! decomp pool
   real(r8) :: adjustment_one_level(bounds%begc:bounds%endc)
   integer :: begc, endc

   character(len=*), parameter :: subname = 'DynamicColumnAdjustments'
   !-----------------------------------------------------------------------

   begc = bounds%begc
   endc = bounds%endc

   


!--------------------------------added hg ---------------------------------------------------------------
  do l = 1, ndecomp_pools
      do j = 1, nlevdecomp
         call column_state_updater%update_column_state_no_special_handling( &
              bounds = bounds, &
              clump_index = clump_index, &
              var    = this%decomp_hgpools_vr_col(begc:endc, j, l), &
              adjustment = adjustment_one_level(begc:endc))
         
      end do
   end do

  !---------------------------------end add--------------------------------------------------------------------------



   
   if (use_nitrif_denitrif) then

      do j = 1, nlevdecomp
         ! Separately track adjustments made to no3 and nh4 pools, since those aren't included
         ! in the overall N balance (totn)
         
         this%dyn_hgbal_adjustments_col(begc:endc) = 0._r8


         

        
        !----------------------added Hg------------------------------
         !call column_state_updater%update_column_state_no_special_handling( &
              !bounds = bounds, &
              !clump_index = clump_index, &
             ! var    = this%soil_hg_vr_col(begc:endc, j), &
              !adjustment = adjustment_one_level(begc:endc))
        ! this%dyn_hgbal_adjustments_col(begc:endc) = &
              !this%dyn_hgbal_adjustments_col(begc:endc) + &
              !adjustment_one_level(begc:endc) * dzsoi_decomp(j)
        !----------------------added Hg end------------------------------    
      end do

   end if

 end subroutine DynamicColumnAdjustments

  !------------------------------------------------------------------------
  subroutine SetTotVgCThresh(this, totvegcthresh)

    class(soilbiogeochem_Hgstate_type)           :: this
    real(r8)                              , intent(in) :: totvegcthresh

    if ( totvegcthresh <= 0.0_r8 )then
        call endrun(msg=' Error totvegcthresh is zero or negative and should be > 0'//&
               errMsg(sourcefile, __LINE__))
    end if
    this%totvegcthresh = totvegcthresh

  end subroutine SetTotVgCThresh

end module SoilBiogeochemHgStateType
