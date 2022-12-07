module SoilBiogeochemHgStateUpdate1Mod
#include "shr_assert.h"
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable updates, non-mortality fluxes.
  !
  ! !USES:
  use atm2lndType                        , only : atm2lnd_type
  use shr_kind_mod                       , only: r8 => shr_kind_r8
  use clm_time_manager                   , only : get_step_size, get_days_per_year
  use clm_varpar                         , only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions,nlevsoi
  use clm_varpar                         , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varctl                         , only : iulog, use_nitrif_denitrif, use_crop,use_cndv
  use clm_varcon                         , only : nitrif_n2o_loss_frac, dzsoi_decomp, zisoi,secspday
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemNitrogenStateType    , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenfluxType     , only : soilbiogeochem_nitrogenflux_type
  use pftconMod               , only : pftcon
  
  use SoilBiogeochemHgStateType          , only : soilbiogeochem_Hgstate_type
  use SoilBiogeochemHgfluxType           , only : soilbiogeochem_Hgflux_type

  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use CNSharedParamsMod                  , only : use_fun
  use ColumnType                         , only : col 
  use CNCStateUpdate1Mod                 , only : f_surf_met, f_surf_cel, f_surf_lig
  use SoilBiogeochemCarbonStateType      , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemCarbonFluxType       , only : soilbiogeochem_carbonflux_type
  use spmdMod        , only : masterproc
  use clm_varctl     , only : iulog, use_vertsoilc, use_fates
  use decompMod      , only : bounds_type

  use SoilBiogeochemCarbonFluxType , only : soilbiogeochem_carbonflux_type
  use CNVegCarbonFluxType , only : cnveg_carbonflux_type

  use PhotosynthesisMod, only : photosyns_type
  use SoilBiogeochemHgFluxType, only : SoilBiogeochem_Hgflux_type
  use shr_log_mod             , only : errMsg => shr_log_errMsg
  use subgridAveMod           , only :   p2c, c2g, p2g
  use SoilBiogeochemHgStateType, only : soilbiogeochem_Hgstate_type
  use WaterStateType           , only : waterstate_type
  use WaterFluxType            , only : waterflux_type
  use PatchType     , only : patch
  use CanopyStateType     , only : canopystate_type
  use CNDVType                , only : dgvs_type
  use clm_varcon     , only : spval
  !use Hg_DryDepVelocity      , only :  hg_drydepvel_type
  !use CNCStateUpdate1Mod, only : soilbiogeochem_carbonflux_inst
  implicit none
  private

  public :: readParams
  public:: SoilBiogeochemHgStateUpdate1

  type, private :: params_type
     real(r8):: sf        ! soluble fraction of mineral N (unitless)
     real(r8):: sf_no3    ! soluble fraction of NO3 (unitless)
     real(r8):: am     ! mortality rate based on annual rate, fractional mortality (1/yr)
     real(r8):: k_mort
  end type params_type
  
  type(params_type), private ::  params_inst

  !
  ! !PUBLIC MEMBER FUNCTIONS:
  
  !-----------------------------------------------------------------------
  character(len=*), parameter, private :: sourcefile = &
       __FILE__
contains


 !-------------------------------读入参数------------------------------
   subroutine readParams ( ncid )
    !
    ! !DESCRIPTION:
    ! Read in parameters
    !
    ! !USES:
    use ncdio_pio   , only : file_desc_t,ncd_io
    use abortutils  , only : endrun
    use shr_log_mod , only : errMsg => shr_log_errMsg
    !
    ! !ARGUMENTS:
    type(file_desc_t),intent(inout) :: ncid   ! pio netCDF file id
    !
    ! !LOCAL VARIABLES:
    character(len=32)  :: subname = 'hg_leachandlitter_type'
    character(len=100) :: errCode = '-Error reading in parameters file:'
    logical            :: readv ! has variable been read in or not
    real(r8)           :: tempr ! temporary to read in constant
    character(len=100) :: tString ! temp. var for reading
    real(r8):: am     ! mortality rate based on annual rate, fractional mortality (1/yr)
    real(r8):: k_mort
    !-----------------------------------------------------------------------
    
    tString='sf_minn'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%sf=tempr

    tString='sf_no3'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%sf_no3=tempr
    
    tString='r_mort'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%am=tempr
    
    tString='k_mort'
    call ncd_io(varname=trim(tString),data=tempr, flag='read', ncid=ncid, readvar=readv)
    if ( .not. readv ) call endrun(msg=trim(errCode)//trim(tString)//errMsg(sourcefile, __LINE__))
    params_inst%k_mort=tempr 


  end subroutine readParams

!-----------------------------------------------------------------------------------



  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemHgStateUpdate1(bounds,num_soilc, filter_soilc,&
       soilbiogeochem_state_inst, soilbiogeochem_Hgflux_inst, soilbiogeochem_Hgstate_inst, atm2lnd_inst, soilbiogeochem_nitrogenflux_inst,&
        soilbiogeochem_nitrogenstate_inst,soilbiogeochem_carbonstate_inst,soilbiogeochem_carbonflux_inst,cnveg_carbonflux_inst,photosyns_inst,&
        waterstate_inst, waterflux_inst,canopystate_inst,dgvs_inst,num_soilp,filter_soilp)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic nitrogen state
    ! variables (except for gap-phase mortality and fire fluxes)
    !
    ! !ARGUMENTS:
    integer                                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilbiogeochem_state_type)         , intent(in)    :: soilbiogeochem_state_inst
    type(soilbiogeochem_nitrogenflux_type)  , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(in) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(in) :: soilbiogeochem_carbonflux_inst
    !--------------------------added Hg----------------------
     type(soilbiogeochem_Hgflux_type)  , intent(inout) :: soilbiogeochem_Hgflux_inst
     type(bounds_type) , intent(in)                :: bounds
     type(cnveg_carbonflux_type), intent(in) :: cnveg_carbonflux_inst
     type(photosyns_type), intent(in):: photosyns_inst
     type(waterstate_type)                   , intent(in)    :: waterstate_inst
     type(waterflux_type)                    , intent(in)    :: waterflux_inst

    !---------------------------added Hg end ----------------

    

    !--------------------------------added Hg-------------------------------

    type(soilbiogeochem_Hgstate_type) , intent(inout) :: soilbiogeochem_Hgstate_inst
    type(atm2lnd_type)                , intent(in)    :: atm2lnd_inst
    type(canopystate_type)            , intent(in)    :: canopystate_inst
    type(dgvs_type)                   , intent(inout) :: dgvs_inst
    integer                           , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                         , intent(in)    :: filter_soilp(:) ! patch filter for soil points
    !type(hg_drydepvel_type)         ,intent(inout)      :: hg_drydepvel_inst
     !---------------------------------added Hg end--------------------------------


    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k,g ! indices
    integer :: fp,fc     ! lake filter indices
   
    real(r8):: dt        ! radiation time step (seconds)
    integer :: begc,endc
    !real(r8)  , intent(in)    :: rrsun_p2c_col(bounds%begc:)
    !real(r8):: hg0_rrsun(bounds%begp:bounds%endp) 
    real(r8):: hg2_rrsun(bounds%begp:bounds%endp)
    real(r8):: hg2_dry_rrsun(bounds%begp:bounds%endp)
    real(r8):: hg_p2c(bounds%begc:bounds%endc)
    real(r8):: hg2_p2c(bounds%begc:bounds%endc)
    real(r8):: hg2_dry_p2c(bounds%begc:bounds%endc)
    real(r8):: leaf_out_1(bounds%begc:bounds%endc)
    real(r8):: leaf_out_2(bounds%begc:bounds%endc)
    real(r8):: leaf_out_3(bounds%begc:bounds%endc)
    real(r8):: hg2_rrsun_col(bounds%begc:bounds%endc)
    real(r8):: hg2dry_leaf_out_patch(bounds%begp:bounds%endp)
    real(r8):: hg2dry_leaf_out_col(bounds%begc:bounds%endc)
    real(r8):: hg2dry_2_leaf_patch(bounds%begp:bounds%endp)
    real(r8):: hg2dry_2_hg2wet_patch(bounds%begp:bounds%endp)
    real(r8):: hg2dry_2_leaf_col(bounds%begc:bounds%endc)
    real(r8):: hg2dry_2_hg2wet_col(bounds%begc:bounds%endc)
    real(r8):: hg2_dry_dep_to_wetHg2pool(bounds%begc:bounds%endc)
    real(r8):: f_soil_5
    real(r8):: f_soil_6
    real(r8):: f_soil_7
    real(r8) :: sf                                     ! soluble fraction of mineral N (unitless)
    real(r8) :: sf_no3                                 ! soluble fraction of NO3 (unitless)
    real(r8) :: disn_conc                              ! dissolved mineral N concentration (gN/kg water)
    real(r8) :: tot_water(bounds%begc:bounds%endc)     ! total column liquid water (kg water/m2)
    real(r8) :: surface_water(bounds%begc:bounds%endc) ! liquid water to shallow surface depth (kg water/m2)
    real(r8) :: drain_tot(bounds%begc:bounds%endc)     ! total drainage flux (mm H2O /s)
    real(r8) :: f_fra_lit(bounds%begc:bounds%endc)
    real(r8) :: f_fra_soil(bounds%begc:bounds%endc)
    real(r8), parameter :: depth_runoff_Hgloss = 0.05
    !real(r8) :: hg_fstom(bounds%begp:bounds%endp)
    real(r8):: k_mort = 0.3
    real(r8):: am            ! rate for fractional mortality (1/yr)
    real(r8):: m             ! rate for fractional mortality (1/s)
    real(r8):: mort_max
    real(r8) :: rs    ! stomatal resistance associated with dry deposition velocity (s/m)
    !real(r8):: leaf_out_7(bounds%begc:bounds%endc,1:nlevdecomp)
    !-----------------------------------------------------------------------
     !SHR_ASSERT_ALL((ubound(rrsun_p2c_col)          == (/bounds%endc/))                               , errMsg(sourcefile, __LINE__))

    

     associate(                                                                  & 
         cascade_donor_pool    => decomp_cascade_con%cascade_donor_pool        , & ! Input:  [integer  (:)     ]  which pool is C taken from for a given decomposition step
         cascade_receiver_pool => decomp_cascade_con%cascade_receiver_pool     , & ! Input:  [integer  (:)     ]  which pool is C added to for a given decomposition step
         ivt                   =>    patch%itype                               , & ! Input:  [integer  (:) ]  patch vegetation type 
         greffic               =>    dgvs_inst%greffic_patch                   , & 
         ndep_prof             => soilbiogeochem_state_inst%ndep_prof_col      , & ! Input:  [real(r8) (:,:)   ]  profile over which N deposition is distributed through column (1/m)
         nfixation_prof        => soilbiogeochem_state_inst%nfixation_prof_col , &
         Hgdep_prof            => soilbiogeochem_state_inst%Hgdep_prof_col     , & ! Input:  [real(r8) (:,:)   ]  profile over which Hg deposition is distributed through column (1/m)
         heatstress            =>    dgvs_inst%heatstress_patch                , &
         woody                 =>    pftcon%woody                              , & ! Input:  binary flag for woody lifeform
        
         rssun_hg                =>photosyns_inst%rssun_patch                   , &
         rssha_hg                =>photosyns_inst%rssha_patch                   , &
         hg_rs                =>photosyns_inst%rs_hg_patch                   , &
         hg_fstom             =>photosyns_inst%fstom_hg_patch                , &
         decomp_cpools_to_fire_vr_col  =>cnveg_carbonflux_inst%m_decomp_cpools_to_fire_vr_col, & 
         nf                    => soilbiogeochem_nitrogenflux_inst             , & ! Output:
         Hgf                   => soilbiogeochem_Hgflux_inst                   , & ! Output:
         Hgs                   => soilbiogeochem_Hgstate_inst                  , &
         cs_soil               => soilbiogeochem_carbonstate_inst              , &
         cf_soil               => soilbiogeochem_carbonflux_inst               , &
         hg_leafout            => soilbiogeochem_Hgstate_inst%Hg_leafout_pools_col, & 
         h2osoi_liq          => waterstate_inst%h2osoi_liq_col                 , & ! Input:  [real(r8) (:,:) ]  liquid water (kg/m2) (new) (-nlevsno+1:nlevgrnd)

         qflx_drain          => waterflux_inst%qflx_drain_col                  , & ! Input:  [real(r8) (:)   ]  sub-surface runoff (mm H2O /s)                    
         qflx_surf           => waterflux_inst%qflx_surf_col                   , & ! Input:  [real(r8) (:)   ]  surface runoff (mm H2O /s)                        
         
         !sminn_vr            => soilbiogeochem_nitrogenstate_inst%sminn_vr_col           , & ! Input:  [real(r8) (:,:) ]  (gN/m3) soil mineral N                          
         !soil_hg2_vr         => soilbiogeochem_Hgstate_inst%soil_hg2_vr_col  ,       & ! Input:  [real(r8) (:,:) ]                                                  

         !sminn_leached_vr    => soilbiogeochem_nitrogenflux_inst%sminn_leached_vr_col    , & ! Output: [real(r8) (:,:) ]  rate of mineral N leaching (gN/m3/s)            
          soil_hg2_leached_vr => soilbiogeochem_Hgflux_inst%soil_hg2_leached_vr_col, &
          forc_Hg0dep      =>  atm2lnd_inst%forc_Hg0dep_grc ,                       & ! Input:  [real(r8) (:)]  nitrogen deposition rate (gN/m2/s)
          forc_hg2dep      =>  atm2lnd_inst%forc_hg2dep_grc ,                       & 
          forc_hg2_dry_dep      =>  atm2lnd_inst%forc_hg2_dry_dep_grc,              &
          laisun_hg0                 => canopystate_inst%laisun_patch             , & ! Input:  [real(r8) (:)   ]  sunlit leaf area                                                      
          laisha_hg0                 => canopystate_inst%laisha_patch             , &
          !velocity_hg           => hg_drydepvel_inst%hg_velocity_patch             , &
          fsun       =>    canopystate_inst%fsun_patch                             , &
          elai       =>    canopystate_inst%elai_patch                             ,&
          ns                    => soilbiogeochem_nitrogenstate_inst                & ! Output:
         )

      ! set time steps
   !-----------------------------------------------计算凋落速率--------------------------------------------------------------
       am = params_inst%am
      ! set coeff of growth efficiency in mortality equation 
       k_mort = params_inst%k_mort
       
    !-----------------------------------------植物叶片气孔吸收参数化---------------------------
 !     do p = bounds%begp,bounds%endp
    
!     end do
   !-----------------------------------------植物叶片气孔吸收参数化---------------------------





      ! patch loop
      

!       if (masterproc) then
!      write(iulog,*)' : The  Output of resistance velocity = ',velocity_hg
!      end if
!     if (masterproc) then
!       write(iulog,*)' : The  Output of resistance sha = ',c_rs2
!       end if

   !---------------------------------------------------------------------------------------------------------------------------------
    

     !do p = bounds%begp, bounds%endp
     do fp = 1,num_soilp
         p = filter_soilp(fp)
         
         hg2dry_leaf_out_patch(p) = 0._r8
         hg2dry_2_leaf_patch(p) = 0._r8
         hg2dry_2_hg2wet_patch(p) = 0._r8
         hg_fstom(p) = 0._r8
         hg2_dry_rrsun(p) = 0._r8
         Hgf%hg0_rrsun(p) = 0._r8
         Hgf%hg0_uptake_veg(p) = 0._r8
         Hgf%hg0_uptake_out_veg(p) = 0._r8
         Hgf%hg0_uptake_fixed(p) = 0._r8
         Hgf%hg0_reemission_veg(p) = 0._r8
    end do
    

      !do p = bounds%begp, bounds%endp  
    dt = real( get_step_size(), r8 ) 
      do fp = 1,num_soilp

         p = filter_soilp(fp)
         
           g = patch%gridcell(p)
           !hg_fstom(p) = ((laisun_hg0(p) + (laisha_hg0(p)/10._r8))/2._r8)/5._r8
            !hg_fstom(p) = (laisun_hg0(p) + laisha_hg0(p))/2._r8/5._r8
!         if (hg_fstom(p)>0.4_r8) then
!             hg_fstom(p) = 0.4_r8
!         end if

!          !------------------------------------计算气孔导度------------------------------
!            if (c_rs1(p)>1._r8 .and. c_rs1(p)<1E4_r8 .and. c_rs2(p)>1._r8 .and. c_rs2(p)<1E4_r8) then
!             hg_fstom(p) = ((1._r8 / c_rs1(p))*0.9_r8 + (1._r8/c_rs2(p))*0.1_r8)*600._r8 
!               if (hg_fstom(p)>1._r8) then
!                   hg_fstom(p) = 0.91_r8
!               end if    
!            else
!              hg_fstom(p) = 0._r8       
!           end if   
         !------------------------------------计算气孔导度------------------------------


         
            Hgf%hg0_rrsun(p) = forc_Hg0dep(g)
            hg2_rrsun(p) = forc_Hg2dep(g)
            hg2_dry_rrsun(p) = forc_hg2_dry_dep(g)
            
            Hgs%hg2dry_leaf__patch(p) = forc_hg2_dry_dep(g)*(3._r8/4._r8)
            
            !Hgs%hg2dry_leaf__patch(p) = Hgs%hg2dry_leaf__patch(p)
          end do 



 do c = bounds%begc, bounds%endc 
  hgs%hg2dry_leafwash_col(c) = 0._r8
               hgs%hg2dry_photored_col(c)= 0._r8
               hgs%tot_hg2photored_col(c)= 0._r8
               hgs%hg2wet_surf_col(c)= 0._r8
               hgs%soilhg_surf_col(c)= 0._r8
               hgs%hg2wet_photored_col(c)= 0._r8
               hgs%hg2wet_heatre_col(c)= 0._r8
               hgf%soil_hg2_leached_col(c) = 0._r8
 end do
  !--------------------------------hg2的相关过程------------------------------------------------------   
  dt = real( get_step_size(), r8 )   
        do fc = 1,num_soilc
         c = filter_soilc(fc)
         g =col%gridcell(c)
         hgs%hg2dry_leafwash_col(c) = forc_hg2_dry_dep(g)*0.25_r8*0.8_r8
         hgs%hg2dry_surf_col(c) =    hgs%hg2dry_leafwash_col(c)

         hgs%hg2dry_photored_col(c) = forc_hg2_dry_dep(g)*(0.25_r8)*0.2_r8

              forc_hg2dep(g) = -forc_hg2dep(g)!转化负号

              hgs%hg2wet_photored_col(c) = forc_hg2dep(g)*0.1_r8
              hgs%hg2wet_heatre_col(c) = forc_hg2dep(g)*0.15_r8
              hgs%hg2wet_surf_col(c)   = forc_hg2dep(g)*(1._r8-0.1_r8-0.15_r8)
 
              hgs%tot_hg2photored_col(c) = hgs%hg2wet_photored_col(c) + hgs%hg2dry_photored_col(c)


      end do


 !--------------------------------------------------------------------------------------------------



      do fp = 1,num_soilp
         p = filter_soilp(fp)
        !if (Hgf%hg0_rrsun(p) > 0._r8 .and. Hgf%hg0_rrsun(p)<1E1_r8) then
          !if (hg2_dry_rrsun(p) > 0._r8 .and. hg2_dry_rrsun(p) < 1._r8) then

!--------------------------------------凋落物周期--------------------------------------------------------
           if (use_cndv) then
            ! Stress mortality from lpj's subr Mortality.

            if (woody(ivt(p)) == 1._r8) then

               if (ivt(p) == 8) then
                  mort_max = 0.03_r8 ! BDT boreal
               else
                  mort_max = 0.01_r8 ! original value for all patches
               end if

               ! heatstress and greffic calculated in Establishment once/yr

               ! Mortality rate inversely related to growth efficiency
               ! (Prentice et al 1993)
               am = mort_max / (1._r8 + k_mort * greffic(p))

               ! Mortality rate inversely related to growth efficiency
               ! (Prentice et al 1993)
               am = mort_max / (1._r8 + k_mort * greffic(p))

               am = min(1._r8, am + heatstress(p))
            else ! lpj didn't set this for grasses; cn does
               ! set the mortality rate based on annual rate
               am = params_inst%am
            end if

         end if

         m  = am/(get_days_per_year() * secspday)
  !--------------------------------------凋落物周期end--------------------------------------------------------   
 !          if (masterproc) then
!       write(iulog,*)' : The  Output of m in hg = ',m
!       end if
                 if (rssun_hg(p)>0._r8 .and. rssun_hg(p)<1.e30 .and. rssha_hg(p)>0._r8 .and. rssha_hg(p)<1.e30.and.elai(p)>0._r8 ) then
                     !LKE: corrected rs to add rssun and rssha in parallel (11/30/2017)
                     rs=1._r8/(fsun(p)*elai(p)/rssun_hg(p) + (1.-fsun(p))*elai(p)/rssha_hg(p))
                  else
                     rs=spval
                  endif

        Hgf%hg0_uptake_veg(p) = Hgf%hg0_rrsun(p)* (1._r8/rs)*200._r8*10._r8*2._r8
        
        Hgf%hg0_uptake_out_veg(p) = Hgf%hg0_uptake_veg(p)*0.15_r8 !计算进到植物体内被还原释放出来的量
        Hgf%hg0_uptake_fixed(p) = Hgf%hg0_uptake_veg(p) - Hgf%hg0_uptake_out_veg(p) ! 减掉还原后存留在叶片内的量

        Hgf%hg0_reemission_veg(p) = (Hgf%hg0_rrsun(p)* (1._r8/rs)*200._r8*10._r8*2._r8)*0.1_r8!未被气孔吸收的部分，再挥发掉

        !Hgf%hg0_rrsun(p) = Hgf%hg0_rrsun(p)*hg_fstom(p)
      
        !叶片凋落后
        Hgf%hg_mt_litter(p) = (Hgf%hg0_uptake_veg(p)+ Hgs%hg2dry_leaf__patch(p))*6.34E-1_r8
        
       !hg2dry_leaf_out_patch(p) = hg2_dry_rrsun(p) *0.2_r8 
       
        !hg2dry_2_leaf_patch(p) = hg2_dry_rrsun(p) * (1._r8/rs)*200._r8*10._r8*0.5_r8
        

       !没有进入气孔的hg2干沉降除了被光解掉的其余进入hg2wet库  
        !hg2dry_2_hg2wet_patch(p) = hg2_dry_rrsun(p) - hg2dry_2_leaf_patch(p)- hg2dry_leaf_out_patch(p)

       !hg2凋落后
       ! hg2dry_2_leaf_patch(p) = hg2dry_2_leaf_patch(p)*6.34E-1_r8/2._r8


     !end if 
  !end if
   end do

   

       do c = bounds%begc, bounds%endc 

         hg_p2c(c) = 0._r8
         hg2_dry_p2c(c) = 0._r8
         hg2dry_2_leaf_col(c) = 0._r8
         hg2dry_2_hg2wet_col(c) = 0._r8
         hg2dry_leaf_out_col(c) = 0._r8

      end do

    call p2c(bounds, num_soilc, filter_soilc, &
         Hgf%hg_mt_litter(bounds%begp:bounds%endp), &
         hg_p2c(bounds%begc:bounds%endc))
         
!     call p2c(bounds, num_soilc, filter_soilc, &
!          hg2_rrsun(bounds%begp:bounds%endp), &
!          hg2_p2c(bounds%begc:bounds%endc))

    call p2c(bounds, num_soilc, filter_soilc, &
         hg2_dry_rrsun(bounds%begp:bounds%endp), &
         hg2_dry_p2c(bounds%begc:bounds%endc)) 
     

    call p2c(bounds, num_soilc, filter_soilc, &
          hg2dry_2_leaf_patch(bounds%begp:bounds%endp), &
          hg2dry_2_leaf_col(bounds%begc:bounds%endc))     

   call p2c(bounds, num_soilc, filter_soilc, &
          hg2dry_2_hg2wet_patch(bounds%begp:bounds%endp), &
          hg2dry_2_hg2wet_col(bounds%begc:bounds%endc)) 

   call p2c(bounds, num_soilc, filter_soilc, &
          hg2dry_leaf_out_patch(bounds%begp:bounds%endp), &
          hg2dry_leaf_out_col(bounds%begc:bounds%endc))

   
!     call p2g(bounds, n_drydep, &
!             drydepvel_inst%velocity_patch (bounds%begp:bounds%endp, :), &
!             lnd2atm_inst%ddvel_grc        (bounds%begg:bounds%endg, :), &
!             p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')
!     endif
   

    call p2g(bounds, &
             Hgf%hg0_uptake_veg (bounds%begp:bounds%endp), &
             atm2lnd_inst%forc_hg0dep_uptake_grc(bounds%begg:bounds%endg), &
            p2c_scale_type='unity', c2l_scale_type= 'unity', l2g_scale_type='unity')

   do c = bounds%begc, bounds%endc
         g = col%gridcell(c)
         !Hg0dep_to_sminn(c) = forc_Hg0dep(g)
         Hgs%soilHg_to_col(c) = Hgs%soilHg_grc(g)
         !hg2_dry_dep_to_sminn(c) = forc_hg2_dry_dep(g)
      end do
    

    !------------------------added Hg---------------------------------------
    dt = real( get_step_size(), r8 )
     do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
   !设soucesinke为0   
        
                
                 
    !到地面的汞
      !Hg2wet经历还原(光解和热蒸发)
      !Hgf%hg2dep_to_sminn_col(c)= Hgf%hg2dep_to_sminn_col(c)*(1._r8-0.2_r8)

      ! 叶片未被吸收而挥发的部分; 根据呼吸排放比，气孔吸收碳汞的阻力系数比为 1/0.16
      !if ( rrsun_p2c(c) > 0._r8 ) then 
      !leaf_out_hg2dry_photored(c) = hg2_rrsun(p) * 0.2
                   

                   !Hgs%soil_Hg0_vr_col(c,j) =  Hgs%soil_Hg0_vr_col(c,j) + hg_p2c(c) *dt*nfixation_prof(c,j)

                  
                  !进入叶子的hg2干沉降
                 ! Hgs%soil_hg2_dry_vr_col(c,j) = Hgs%soil_hg2_dry_vr_col(c,j) + hg2dry_2_leaf_col(c) *dt * nfixation_prof(c,j)
                  
                
                  !汇集到土壤表层的所有汞
                   hgs%soilhg_surf_col(c) = hgs%hg2dry_surf_col(c) + hgs%hg2wet_surf_col(c) + hg_p2c(c)
                   !分配到垂直土壤层
                   hgs%soilhg_surf_vr_col(c,j) = hgs%soilhg_surf_vr_col(c,j) + hgs%soilhg_surf_col(c)*dt*nfixation_prof(c,j)

                   !土壤初始汞浓度
                  Hgs%initial_soilHg_col(c,j) =  Hgs%soilHg_to_col(c)*nfixation_prof(c,j)




               !soil_hg2_vr_col代表所有汇聚到地表的hg2湿(hg2前加负号统一单位)
             !Hgs%soil_hg2_vr_col(c,j) = Hgs%soil_hg2_vr_col(c,j) + (-Hgf%hg2dep_to_sminn_col(c)+hg2dry_2_hg2wet_col(c)) *dt * nfixation_prof(c,j)


           !经过淋溶

            !Hgs%soil_hg2_vr_col(c,j) = (-Hgs%soil_hg2_vr_col(c,j)) - (hgf%soil_hg2_leached_vr_col(c,j) * dt)





!     Hgs%decomp_hg_rest_vr_col(c,j) = (Hgs%soil_Hg0_vr_col(c,j) - Hgs%soil_hg2_vr_col(c,j)) !保留到土壤中的汞- ((Hgs%soil_Hg0_vr_col(c,j) + Hgs%soil_hg2_vr_col(c,j))*nf%anaerobic_frac_col)
!     Hgs%decomp_hg_rest_vr_col(c,j) = Hgs%decomp_hg_rest_vr_col(c,j) - (Hgs%decomp_hg_rest_vr_col(c,j)*nf%anaerobic_frac_col(c,j)*1.2)!减去异养呼吸所消耗的汞
    
                
                  



      

                  
                 !-------------------------added Hg end ---------------------------------



                
        !end if   
          
         end do
         
      end do


   
   



   !---------------------Hg0干沉降到达地表前的过程-------------------- ！  
    !Hg0经历挥发(没进入叶片气孔的落到地表后挥发掉)
    !Hgf%Hg0dep_to_sminn_col = Hgf%Hg0dep_to_sminn_col*(1._r8-0.7_r8)

    !Hg0受气孔阻力后进入叶片
     !Hgf%Hg0dep_to_sminn_col = Hgf%Hg0dep_to_sminn_col*(1._r8/rrsun_p2c)
   !----------------------------------------------------------------------
   
   !---------------------------计算 淋溶过程--------------------------------------------

    ! set time steps
      dt = real( get_step_size(), r8 )

      if (.not. use_nitrif_denitrif) then
         ! set constant sf 
         sf = params_inst%sf
      else
         ! Assume that 100% of the soil NO3 is in a soluble form
         sf_no3 =  params_inst%sf_no3 
      end if

      ! calculate the total soil water
      tot_water(bounds%begc:bounds%endc) = 0._r8
      do j = 1,nlevsoi
         do fc = 1,num_soilc
            c = filter_soilc(fc)
            tot_water(c) = tot_water(c) + h2osoi_liq(c,j)
         end do
      end do

      ! for runoff calculation; calculate total water to a given depth
      surface_water(bounds%begc:bounds%endc) = 0._r8
      do j = 1,nlevsoi
         if ( zisoi(j) <= depth_runoff_Hgloss)  then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               surface_water(c) = surface_water(c) + h2osoi_liq(c,j)
            end do
         elseif ( zisoi(j-1) < depth_runoff_Hgloss)  then
            do fc = 1,num_soilc
               c = filter_soilc(fc)
               surface_water(c) = surface_water(c) + h2osoi_liq(c,j) * ( (depth_runoff_Hgloss - zisoi(j-1)) / col%dz(c,j))
            end do
         endif
      end do

      ! Loop through columns
      do fc = 1,num_soilc
         c = filter_soilc(fc)
         drain_tot(c) = qflx_drain(c)
      end do


      !if (masterproc) then
      !write(iulog,*)' : The  Output of drain_tot = ',drain_tot
      !end if

     !if (masterproc) then
      !write(iulog,*)' : The  Output of tot_water = ',tot_water
     ! end if



     !if (masterproc) then
      !write(iulog,*)' : The  Output of sf_no3 = ',sf_no3
      !end if

     !if (masterproc) then
      !write(iulog,*)' : The  output of soil_hg2_vr_col = ',Hgs%soil_hg2_vr_col
      !end if

 
         do j = 1,nlevdecomp
            ! Loop through columns
            do fc = 1,num_soilc
               c = filter_soilc(fc)

              
                  if (tot_water(c) > 0._r8) then
                     disn_conc = (sf_no3 * (hgs%soilhg_surf_vr_col(c,j)) )/tot_water(c)
                  end if
                
     
                  ! calculate the N leaching flux as a function of the dissolved
                  ! concentration and the sub-surface drainage flux
                  soil_hg2_leached_vr(c,j) = disn_conc * drain_tot(c)
              
                  !
                 
                  !
                  ! ensure that leaching rate isn't larger than soil N pool
                 soil_hg2_leached_vr(c,j) = min(soil_hg2_leached_vr(c,j), (hgs%soilhg_surf_vr_col(c,j)) / dt )
                  !
                  ! limit the leaching flux to a positive value
                 soil_hg2_leached_vr(c,j) = max(soil_hg2_leached_vr(c,j), 0._r8)
                  !
                  !
                  

               
               ! limit the flux based on current soil_hg2 state
               ! only let at most the assumed soluble fraction
               ! of soil_hg2 be leached on any given timestep
               soil_hg2_leached_vr(c,j) = min(soil_hg2_leached_vr(c,j), (sf_no3 * (hgs%soilhg_surf_vr_col(c,j))/dt))

               ! limit the flux to a positive value
               soil_hg2_leached_vr(c,j) = max(soil_hg2_leached_vr(c,j), 0._r8)

               


            end do
         end do

   !整合算一下所有土壤层的淋溶
     do j = 1, nlevdecomp
          do fc = 1,num_soilc
             c = filter_soilc(fc)
          hgf%soil_hg2_leached_col(c) = &
                  hgf%soil_hg2_leached_col(c) + &
                  hgf%soil_hg2_leached_vr_col(c,j) * dzsoi_decomp(j)
          end do
    end do  

   !--------------------------------淋溶过程计算结束------------------------ -----------------！
  

    do k = 1, 7
      do j = 1, 25
         do c = bounds%begc, bounds%endc 
             
        Hgf%decomp_hgpools_sourcesink_col(c,j,k) = 0._r8

          end do 
       end do
   end do    



     dt = real( get_step_size(), r8 )
      
       do j = 1, nlevdecomp
         do fc = 1,num_soilc
            c = filter_soilc(fc)
      
           
     
     !if ( rrsun_p2c(c) > 0._r8  ) then
        
      
     !! Hgs%decomp_hg0pools_vr_col(c,j,i_met_lit) = Hgs%decomp_hg0pools_vr_col(c,j,i_met_lit) + (Hgf%Hg0dep_to_sminn_col(c)*(1._r8/rrsun_p2c(c))*(1._r8/0.16_r8)*dt *nfixation_prof(c,j))*f_surf_met
     !! Hgs%decomp_hg0pools_vr_col(c,j,i_cel_lit) = Hgs%decomp_hg0pools_vr_col(c,j,i_cel_lit) + (Hgf%Hg0dep_to_sminn_col(c)*(1._r8/rrsun_p2c(c))*(1._r8/0.16_r8)*dt *nfixation_prof(c,j))*f_surf_cel
     !! Hgs%decomp_hg0pools_vr_col(c,j,i_lig_lit) = Hgs%decomp_hg0pools_vr_col(c,j,i_lig_lit) + ((Hgf%Hg0dep_to_sminn_col(c)+Hgf%hg2dep_to_sminn_col(c))*(1._r8/rrsun_p2c(c))*(1._r8/0.16_r8)*dt * nfixation_prof(c,j))*f_surf_lig
     
    !!  Hgs%decomp_hg0pools_vr_col(c,j,5) = cs_soil%decomp_cpools_vr_col(c,j,5)*3E-7_r8
     !! Hgs%decomp_hg0pools_vr_col(c,j,6) = cs_soil%decomp_cpools_vr_col(c,j,6)*3E-7_r8
    !!  Hgs%decomp_hg0pools_vr_col(c,j,7) = cs_soil%decomp_cpools_vr_col(c,j,7)*3E-7_r8
    !! Hgs%decomp_hg0pools_vr_col(c,j,i_cwd) = cs_soil%decomp_cpools_vr_col(c,j,i_cwd)*3E-7_r8

      !if (cs_soil%decomp_cpools_vr_col(c,j,5) > 0._r8 .and. cs_soil%decomp_cpools_vr_col(c,j,6) > 0._r8 .and. cs_soil%decomp_cpools_vr_col(c,j,7) > 0._r8) then
      ! f_soil_6(c,j,6) = (cs_soil%decomp_cpools_vr_col(c,j,6))/(cs_soil%decomp_cpools_vr_col(c,j,5) + cs_soil%decomp_cpools_vr_col(c,j,6) + cs_soil%decomp_cpools_vr_col(c,j,7))
     !f_soil_7(c,j,7) = (cs_soil%decomp_cpools_vr_col(c,j,7))/(cs_soil%decomp_cpools_vr_col(c,j,5) + cs_soil%decomp_cpools_vr_col(c,j,6) + cs_soil%decomp_cpools_vr_col(c,j,7))
      ! end if 
     
     !---------------hg2结合土壤有机质之前被淋溶的部分得到最终保留到土壤表层的汞-----------------------------
         hgs%soilhg_surf_vr_col(c,j) = hgs%soilhg_surf_vr_col(c,j) - &
                    hgf%soil_hg2_leached_vr_col(c,j) * dt
      !--------------------------------------------------------------------------


     !保留的格式
     !Hgs%decomp_hg0pools_vr_col(c,j,i_met_lit) = Hgs%decomp_hg0pools_vr_col(c,j,i_met_lit) + (Hgs%soil_Hg0_vr_col(c,j)+Hgs%soil_hg2_dry_vr_col(c,j))*f_surf_met
     Hgs%decomp_hgpools_vr_col(c,j,i_met_lit) = Hgs%decomp_hgpools_vr_col(c,j,i_met_lit) + hgs%soilhg_surf_vr_col(c,j)*f_surf_met
     Hgs%decomp_hgpools_vr_col(c,j,i_cel_lit) = Hgs%decomp_hgpools_vr_col(c,j,i_cel_lit) + hgs%soilhg_surf_vr_col(c,j)*f_surf_cel
     Hgs%decomp_hgpools_vr_col(c,j,i_lig_lit) = Hgs%decomp_hgpools_vr_col(c,j,i_lig_lit) + hgs%soilhg_surf_vr_col(c,j)*f_surf_lig
     Hgs%decomp_hgpools_vr_col(c,j,i_cwd) = cs_soil%decomp_cpools_vr_col(c,j,i_cwd)*3E-7_r8
     
     

     f_soil_5 = 0.33_r8
     f_soil_6 = 0.33_r8
     f_soil_7 = 0.33_r8


      Hgs%decomp_hgpools_vr_col(c,j,5) = Hgs%initial_soilHg_col(c,j)*f_soil_5
      Hgs%decomp_hgpools_vr_col(c,j,6) = Hgs%initial_soilHg_col(c,j)*f_soil_6
      Hgs%decomp_hgpools_vr_col(c,j,7) = Hgs%initial_soilHg_col(c,j)*f_soil_7

      
     
!         if (cs_soil%totsomc_col(c)>0 .and. cs_soil%totlitc_col(c)>0) then

!           f_fra_lit(c) = cs_soil%totlitc_col(c)/(cs_soil%totsomc_col(c)+cs_soil%totlitc_col(c))

!            f_fra_soil(c) = cs_soil%totsomc_col(c)/(cs_soil%totsomc_col(c)+cs_soil%totlitc_col(c))  

    !hg2进来后
          !Hgs%decomp_hg0pools_vr_col(c,j,i_met_lit) = Hgs%decomp_hg0pools_vr_col(c,j,i_met_lit) + (Hgs%soil_hg2_vr_col(c,j))*f_surf_met

          !Hgs%decomp_hg0pools_vr_col(c,j,i_cel_lit) = Hgs%decomp_hg0pools_vr_col(c,j,i_cel_lit) + (Hgs%soil_hg2_vr_col(c,j))*f_surf_cel
      
          !Hgs%decomp_hg0pools_vr_col(c,j,i_lig_lit) = Hgs%decomp_hg0pools_vr_col(c,j,i_lig_lit) + (Hgs%soil_hg2_vr_col(c,j))*f_surf_lig  

 
!       Hgs%decomp_hg0pools_vr_col(c,j,i_met_lit) = Hgs%decomp_hg0pools_vr_col(c,j,i_met_lit) 

!      Hgs%decomp_hg0pools_vr_col(c,j,i_cel_lit) = Hgs%decomp_hg0pools_vr_col(c,j,i_cel_lit) 
      
!      Hgs%decomp_hg0pools_vr_col(c,j,i_lig_lit) = Hgs%decomp_hg0pools_vr_col(c,j,i_lig_lit) 




!      Hgs%decomp_hg0pools_vr_col(c,j,5) = Hgs%decomp_hg0pools_vr_col(c,j,5) + (Hgs%soil_hg2_vr_col(c,j))*f_soil_5

!      Hgs%decomp_hg0pools_vr_col(c,j,6) = Hgs%decomp_hg0pools_vr_col(c,j,6) + (Hgs%soil_hg2_vr_col(c,j))*f_soil_6
      
!      Hgs%decomp_hg0pools_vr_col(c,j,7) = Hgs%decomp_hg0pools_vr_col(c,j,7) + (Hgs%soil_hg2_vr_col(c,j))*f_soil_7


     
     

    
      !Hgf%decomp_hgpools_sourcesink_col(c,j,i_met_lit) = Hgs%soil_Hg0_vr_col(c,j)*f_surf_met*dt 
      !Hgf%decomp_hgpools_sourcesink_col(c,j,i_cel_lit) = Hgs%soil_Hg0_vr_col(c,j)*f_surf_cel*dt 
      !Hgf%decomp_hgpools_sourcesink_col(c,j,i_lig_lit) = Hgs%soil_Hg0_vr_col(c,j)*f_surf_lig*dt 
     
      !Hgf%decomp_hgpools_sourcesink_col(c,j,5) = Hgs%decomp_hg0pools_vr_col(c,j,5)*dt 
      !Hgf%decomp_hgpools_sourcesink_col(c,j,6) = Hgs%decomp_hg0pools_vr_col(c,j,6)*dt 
      !Hgf%decomp_hgpools_sourcesink_col(c,j,7) = Hgs%decomp_hg0pools_vr_col(c,j,7)*dt 
      !Hgf%decomp_hgpools_sourcesink_col(c,j,i_cwd) = Hgs%decomp_hg0pools_vr_col(c,j,i_cwd)*dt 

     
     
       Hgf%decomp_hgpools_sourcesink_col(c,j,i_met_lit) = 0._r8
       Hgf%decomp_hgpools_sourcesink_col(c,j,i_cel_lit) = 0._r8 
       Hgf%decomp_hgpools_sourcesink_col(c,j,i_lig_lit) = 0._r8 
     
       Hgf%decomp_hgpools_sourcesink_col(c,j,5) = 0._r8 
       Hgf%decomp_hgpools_sourcesink_col(c,j,6) = 0._r8 
       Hgf%decomp_hgpools_sourcesink_col(c,j,7) = 0._r8 
       Hgf%decomp_hgpools_sourcesink_col(c,j,i_cwd) = 0._r8
 
    !end if

      !end if



         end do
    end do
  

 !---------------------------------Hg0没有进入叶片的部分------------------------------------------
 !    dt = real( get_step_size(), r8 )
 !        do j = 1._r8, nlevdecomp
 !          do fc = 1._r8,num_soilc
 !             c = filter_soilc(fc)
! !   if ( rrsun_p2c(c) > 0._r8 ) then
! !     leaf_out_1(c) = leaf_out_1(c) + Hgs%decomp_hg0pools_vr_col(c,j,i_met_lit) + &
! !                    (Hgf%Hg0dep_to_sminn_col(c)*(1._r8-(1._r8/rrsun_p2c(c)))*(1._r8/0.16_r8)*dt *nfixation_prof(c,j))*f_surf_met * dzsoi_decomp(j)

! !     leaf_out_2(c) = leaf_out_2(c) + Hgs%decomp_hg0pools_vr_col(c,j,i_cel_lit) + &
! !                    (Hgf%Hg0dep_to_sminn_col(c)*(1._r8-(1._r8/rrsun_p2c(c)))*(1._r8/0.16_r8)*dt *nfixation_prof(c,j))*f_surf_cel * dzsoi_decomp(j)
    
! !     leaf_out_3(c) = leaf_out_3(c) + Hgs%decomp_hg0pools_vr_col(c,j,i_lig_lit) + &
! !                    (Hgf%Hg0dep_to_sminn_col(c)*(1._r8-(1._r8/rrsun_p2c(c)))*(1._r8/0.16_r8)*dt *nfixation_prof(c,j))*f_surf_lig * dzsoi_decomp(j)
! !   end if
!     end do
!   end do 
! hg_leafout = leaf_out_1(c) + leaf_out_2(c) + leaf_out_3(c) +leaf_out_hg2dry_photored(c)
!---------------------------------Hg0没有进入叶片的部分-------------------------------------------------------
   
  


  !if (masterproc) then
    !write(iulog,*)' : The  result of hgpool0220 = ',Hgs%decomp_hg0pools_vr_col
  ! end if
   !野火过程
     
  !if (masterproc) then
   ! write(iulog,*)' : The  Output of cpools_to_fire_vr = ', decomp_cpools_to_fire_vr_col
  !end if
    !Hgs%decomp_hg0pools_vr_col = Hgs%decomp_hg0pools_vr_col-(decomp_cpools_to_fire_vr_col*0.16_r8)


   !do k = 1, 7
        ! do j = 1, nlevdecomp
            !do fc = 1,num_soilc
               !c = filter_soilc(fc)
              ! Hgs%decomp_hg0pools_vr_col(c,j,l) = Hgs%decomp_hg0pools_vr_col(c,j,l) - &
                    !(decomp_cpools_to_fire_vr_col(c,j,l)*3E-7_r8*0.16_r8) * dt
           ! end do
         !end do
    !  end do


  !淋溶过程






    !if (.not. use_fates) then  copy待处理
       

    !dt = real( get_step_size(), r8 )
     !do j = 1, nlevdecomp
        ! do fc = 1,num_soilc
           ! c = filter_soilc(fc)
      

       !Hgs%decomp_hgpools_vr_col(c,j,i_cel_lit) = Hgs%decomp_hgpools_vr_col(c,j,i_cel_lit) +Hgs%soil_Hg2_vr_col(c,j) * dt

      !end do
         
      !end do




  dt = real( get_step_size(), r8 )
  do k = 1, ndecomp_cascade_transitions
       do j = 1,nlevdecomp
            do fc = 1,num_soilc
               c = filter_soilc(fc)
              Hgf%decomp_hgpools_sourcesink_col(c,j,cascade_donor_pool(k)) = &
                    Hgf%decomp_hgpools_sourcesink_col(c,j,cascade_donor_pool(k)) &
                   - ( hgf%decomp_cascade_hr_hg_vr_col(c,j,k)*(0.16_r8) + Hgf%decomp_cascade_hgtransfer_vr_col(c,j,k)) *dt
   
          !if (Hgf%decomp_hgpools_sourcesink_col(c,j,cascade_donor_pool(k)) < 0._r8) then
             !Hgf%decomp_hgpools_sourcesink_col(c,j,cascade_donor_pool(k)) = 0._r8
         ! end if   
          if (abs(Hgf%decomp_hgpools_sourcesink_col(c,j,cascade_donor_pool(k))) > 1.e20_r8) then
            Hgf%decomp_hgpools_sourcesink_col(c,j,cascade_donor_pool(k)) = 0._r8
            end if
           end do
         end do
      end do

    do k = 1, ndecomp_cascade_transitions
         if ( cascade_receiver_pool(k) /= 0 ) then  ! skip terminal transitions
          do j = 1,nlevdecomp
           do fc = 1,num_soilc
               c = filter_soilc(fc)
               Hgf%decomp_hgpools_sourcesink_col(c,j,cascade_receiver_pool(k)) = &
                     Hgf%decomp_hgpools_sourcesink_col(c,j,cascade_receiver_pool(k)) &
                      + Hgf%decomp_cascade_hgtransfer_vr_col(c,j,k)*dt
!
              if (abs(Hgf%decomp_hgpools_sourcesink_col(c,j,cascade_receiver_pool(k))) > 1.e20_r8) then
             Hgf%decomp_hgpools_sourcesink_col(c,j,cascade_receiver_pool(k)) = 0._r8
             end if
              
             end do
            end do
        end if
    end do



 !if (masterproc) then
    ! write(iulog,*)' : The  output of single culumn hgpool-source = ',Hgf%decomp_hgpools_sourcesink_col(1,:,:)
  ! end if  
   
!if (masterproc) then
     !write(iulog,*)' : The  output of single layer hgpool-source = ',Hgf%decomp_hgpools_sourcesink_col(:,1,:)
  ! end if  

!if (masterproc) then
    ! write(iulog,*)' : The  output of seven hgpool-source = ',Hgf%decomp_hgpools_sourcesink_col(1:,1,:)
   !end if  


!if (masterproc) then
     !write(iulog,*)' : The  output of mixed hgpool-source = ',Hgf%decomp_hgpools_sourcesink_col(:,:,1)
  ! end if  

!if (masterproc) then
     !write(iulog,*)' : The  output of mixed hgpool-source = ',Hgf%decomp_hgpools_sourcesink_col(:,:,6)
   !end if  






!if (masterproc) then
     !write(iulog,*)' : The  output of all column layer hgpool-source = ',Hgf%decomp_hgpools_sourcesink_col(:,:,:)
  ! end if











  



  end associate

  end subroutine SoilBiogeochemHgStateUpdate1

end module SoilBiogeochemHgStateUpdate1Mod
