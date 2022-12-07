module SoilBiogeochemHgStateUpdate3Mod

#include "shr_assert.h"

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for nitrogen state variable updates, non-mortality fluxes.
  !
  ! !USES:
  use shr_kind_mod                       , only: r8 => shr_kind_r8
  use clm_time_manager                   , only : get_step_size
  use clm_varpar                         , only : nlevdecomp, ndecomp_pools, ndecomp_cascade_transitions
  use clm_varpar                         , only : i_met_lit, i_cel_lit, i_lig_lit, i_cwd
  use clm_varctl                         , only : iulog, use_nitrif_denitrif, use_crop
  use clm_varcon                         , only : nitrif_n2o_loss_frac, dzsoi_decomp
  use SoilBiogeochemStateType            , only : soilbiogeochem_state_type
  use SoilBiogeochemNitrogenStateType    , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemNitrogenfluxType     , only : soilbiogeochem_nitrogenflux_type
  
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
  use CNVegStateType  , only : cnveg_state_type
  use CNFireMethodMod                    , only : cnfire_method_type
  use GridcellType                       , only : grc
  use CNFireBaseMod                      , only : cnfire_const
  use shr_log_mod                        , only : errMsg => shr_log_errMsg
  !use CNCStateUpdate1Mod, only : soilbiogeochem_carbonflux_inst

  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: SoilBiogeochemHgStateUpdate3
  !-----------------------------------------------------------------------
   character(len=*), parameter, private :: sourcefile = &
       __FILE__
   contains

  !-----------------------------------------------------------------------
  subroutine SoilBiogeochemHgStateUpdate3(bounds,num_soilc, filter_soilc,&
       soilbiogeochem_state_inst, soilbiogeochem_Hgflux_inst, soilbiogeochem_Hgstate_inst,&
        soilbiogeochem_nitrogenflux_inst, &
        soilbiogeochem_carbonstate_inst,soilbiogeochem_carbonflux_inst,cnveg_state_inst,totsomhg_col)
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
    !type(soilbiogeochem_nitrogenstate_type) , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_carbonstate_type)   , intent(in) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_carbonflux_type)    , intent(in) :: soilbiogeochem_carbonflux_inst
    !--------------------------added Hg----------------------
     type(soilbiogeochem_Hgflux_type)  , intent(inout) :: soilbiogeochem_Hgflux_inst
     type(bounds_type) , intent(in)                :: bounds
     !type(cnveg_carbonflux_type), intent(in) :: cnveg_carbonflux_inst
    !---------------------------added Hg end ----------------

    

    !--------------------------------added Hg-------------------------------
    type(cnveg_state_type)         , intent(inout) :: cnveg_state_inst
    type(soilbiogeochem_Hgstate_type) , intent(inout) :: soilbiogeochem_Hgstate_inst

    !---------------------------------added Hg end--------------------------------


    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k,g ! indices
    integer :: fp,fc     ! lake filter indices
    real(r8):: dt        ! radiation time step (seconds)
    integer :: begc,endc
    real(r8):: f                    ! rate for fire effects (1/s)
    real(r8) , intent(in)    :: totsomhg_col(bounds%begc:) 
    !-----------------------------------------------------------------------

      
    SHR_ASSERT_ALL((ubound(totsomhg_col)          == (/bounds%endc/))                               , errMsg(sourcefile, __LINE__))
    
    associate(                                                                   & 
         cascade_donor_pool    => decomp_cascade_con%cascade_donor_pool        , & ! Input:  [integer  (:)     ]  which pool is C taken from for a given decomposition step
         cascade_receiver_pool => decomp_cascade_con%cascade_receiver_pool     , & ! Input:  [integer  (:)     ]  which pool is C added to for a given decomposition step
        
         ndep_prof             => soilbiogeochem_state_inst%ndep_prof_col      , & ! Input:  [real(r8) (:,:)   ]  profile over which N deposition is distributed through column (1/m)
         nfixation_prof        => soilbiogeochem_state_inst%nfixation_prof_col , &
         Hgdep_prof            => soilbiogeochem_state_inst%Hgdep_prof_col     , & ! Input:  [real(r8) (:,:)   ]  profile over which Hg deposition is distributed through column (1/m)

         
         !decomp_cpools_to_fire_vr_col  =>cnveg_carbonflux_inst%m_decomp_cpools_to_fire_vr_col, & 
         nf                    => soilbiogeochem_nitrogenflux_inst             , & ! Output:
         Hgf                   => soilbiogeochem_Hgflux_inst                   , & ! Output:
         somhg_fire            => soilbiogeochem_Hgflux_inst%somhg_fire_col    , & 

         Hgs                   => soilbiogeochem_Hgstate_inst                  , &
         totsomhg             => soilbiogeochem_Hgstate_inst%totsomhg_col    , &
         decomp_hgpools_vr    => soilbiogeochem_Hgstate_inst%decomp_hgpools_vr_col, &
         is_cwd                => decomp_cascade_con%is_cwd                    , &                                                         
         is_litter             => decomp_cascade_con%is_litter                 , &
         baf_crop              => cnveg_state_inst%baf_crop_col                , &
         baf_peatf             => cnveg_state_inst%baf_peatf_col               , &
         !cmb_cmplt_fact        => cnfire_const%cmb_cmplt_fact                  , &
         farea_burned          => cnveg_state_inst%farea_burned_col            , &
         m_decomp_hgpools_to_fire_vr => soilbiogeochem_Hgstate_inst%m_decomp_hgpools_to_fire_vr_col, &
         cs_soil               => soilbiogeochem_carbonstate_inst              , &
         cf_soil               => soilbiogeochem_carbonflux_inst                 &
         )
 
      
      ! set time steps
    
          !f (masterproc) then
      !write(iulog,*)' : The  Output of soil_hg_vr_col = ',Hgs%soil_hg_vr_col
      !end if


      !if (masterproc) then
    !write(iulog,*)' : The  Output of totsomhg_col = ', totsomhg
   ! end if

   !野火过程
    dt = real( get_step_size(), r8 )
     do fc = 1,num_soilc
        c = filter_soilc(fc)

        f = farea_burned(c) 
    
        do j = 1, nlevdecomp
           ! hg fluxes
           do l = 1, ndecomp_pools
              if ( is_litter(l) ) then
                 m_decomp_hgpools_to_fire_vr(c,j,l) = decomp_hgpools_vr(c,j,l) * f * &
                      0.5_r8
              end if
              if ( is_cwd(l) ) then
                 m_decomp_hgpools_to_fire_vr(c,j,l) = decomp_hgpools_vr(c,j,l) * &
                      (f-baf_crop(c)) * 0.25_r8

           decomp_hgpools_vr(c,j,l) = decomp_hgpools_vr(c,j,l) - &
                    (m_decomp_hgpools_to_fire_vr(c,j,l)*0.16_r8) * dt


              end if
           end do
        end do 
      end do   

!-----------------------------------汇总----------------------------------------------

   
!-----------------------------------汇总end----------------------------------------------





   !泥炭地燃烧（peat burning)
    

    do fc = 1,num_soilc
        c = filter_soilc(fc)
        g = col%gridcell(c)
        !if (masterproc) then
        !write(iulog,*)' : The  Output of baf in column= ', baf_peatf(c)
         !end if
        if( grc%latdeg(g)  <  40._r8)then
           somhg_fire(c)= totsomhg(c)*baf_peatf(c)*0.16_r8*(6.0_r8/33.9_r8)
        else
           somhg_fire(c)= baf_peatf(c)*2.2e3_r8*0.16_r8
        end if
     end do

  

  !if (masterproc) then
   ! write(iulog,*)' : The  Output of somhg_fire = ', somhg_fire
  !end if

   !if (masterproc) then
   ! write(iulog,*)' : The  Output of m_decomp_hgpools_to_fire_vr0220 = ', m_decomp_hgpools_to_fire_vr
   !end if
    !Hgs%decomp_hgpools_vr_col = Hgs%decomp_hgpools_vr_col-(decomp_cpools_to_fire_vr_col*0.16_r8)

 
     
    !do l = 1, ndecomp_pools
        !do j = 1, nlevdecomp
            !do fc = 1,num_soilc
              ! c = filter_soilc(fc)
               !decomp_hgpools_vr(c,j,l) = decomp_hgpools_vr(c,j,l) - &
                   ! (m_decomp_hgpools_to_fire_vr(c,j,l)*0.16_r8) * dt
            !end do
        !end do
     !end do

    
     !if (masterproc) then
    !write(iulog,*)' : The  Output of final decomp_hgpools_vr = ', decomp_hgpools_vr
  !end if

   




    



 
  



  end associate

  end subroutine SoilBiogeochemHgStateUpdate3

end module SoilBiogeochemHgStateUpdate3Mod
