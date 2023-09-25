module ErosionMod 

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate erosion induced soil particulate C, N, Hg fluxes 
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use clm_time_manager  , only : get_step_size
  use clm_varcon        , only : spval, dzsoi_decomp
  use clm_varpar        , only : ndecomp_pools, nlevdecomp
  use CNVegNitrogenFluxType           , only : cnveg_nitrogenflux_type
  use CNVegNitrogenStateType          , only : cnveg_nitrogenstate_type
  use CNVegCarbonFluxType             , only : cnveg_carbonflux_type
  use CNVegCarbonStateType            , only : cnveg_carbonstate_type
  use SoilBiogeochemNitrogenFluxType  , only : soilbiogeochem_nitrogenflux_type
  use SoilBiogeochemCarbonFluxType    , only : soilbiogeochem_carbonflux_type
  use SoilBiogeochemCarbonStateType   , only : soilbiogeochem_carbonstate_type
  use SoilBiogeochemNitrogenStateType  , only : soilbiogeochem_nitrogenstate_type
  use SoilBiogeochemDecompCascadeConType , only : decomp_cascade_con
  use SedFluxType       , only : sedflux_type
  use SoilStateType     , only : soilstate_type
  use SoilBiogeochemHgStateType       , only : soilbiogeochem_Hgstate_type
  use SoilBiogeochemHgFluxType        , only : soilbiogeochem_Hgflux_type
  use shr_infnan_mod       , only : isnan => shr_infnan_isnan, nan => shr_infnan_nan, assignment(=)
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ErosionFluxes           ! Calculate erosion introduced CN, Hg fluxes 
  ! !MODULE CONSTANTS
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine ErosionFluxes (bounds, num_soilc, filter_soilc, &
    soilstate_inst, sedflux_vars, cnveg_carbonstate_inst, cnveg_carbonflux_inst,cnveg_nitrogenstate_inst,   &
                                         soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst,                          &              
                                         cnveg_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst,                      &
                                         soilbiogeochem_nitrogenflux_inst, soilbiogeochem_Hgflux_inst, soilbiogeochem_Hgstate_inst)
    !
    ! !DESCRIPTION:
    ! Calculate erosion introduced soil C, N, P fluxes 
    !
    ! !USES:
    use clm_varctl      , only : iulog
    use spmdMod         , only : masterproc
    !
    ! !ARGUMENTS:
    type(bounds_type)        , intent(in)    :: bounds
    integer                  , intent(in)    :: num_soilc       ! number of column soil points in column filter
    integer                  , intent(in)    :: filter_soilc(:) ! filter for soil columns
    type(soilstate_type)     , intent(in)    :: soilstate_inst
    type(sedflux_type)       , intent(in)    :: sedflux_vars
    type(cnveg_carbonflux_type)            , intent(in)    :: cnveg_carbonflux_inst
    type(cnveg_carbonstate_type)           , intent(inout) :: cnveg_carbonstate_inst
    type(cnveg_nitrogenflux_type)            , intent(in)    :: cnveg_nitrogenflux_inst
    type(cnveg_nitrogenstate_type)           , intent(inout) :: cnveg_nitrogenstate_inst
    type(soilbiogeochem_carbonstate_type)    , intent(inout) :: soilbiogeochem_carbonstate_inst
    type(soilbiogeochem_nitrogenstate_type)  , intent(inout) :: soilbiogeochem_nitrogenstate_inst
    type(soilbiogeochem_carbonflux_type) , intent(inout) :: soilbiogeochem_carbonflux_inst
    type(soilbiogeochem_nitrogenflux_type) , intent(inout) :: soilbiogeochem_nitrogenflux_inst
    type(soilbiogeochem_Hgflux_type),   intent(inout) :: soilbiogeochem_Hgflux_inst
    type(soilbiogeochem_Hgstate_type),   intent(inout) :: soilbiogeochem_Hgstate_inst
    
    !
    ! !LOCAL VARIABLES:
    integer  :: c, fc, l, j                              ! indices
    integer  :: kc, kn, khg                              ! new layer indices
    integer  :: k1, k2                                   ! new layer indices     
    real(r8) :: dh, dhn                                  ! erosion height (m)
    real(r8) :: dm_ero, dm_yld, dm                       ! erosion mass (kg/m2)
    real(r8) :: soc_yld                                  ! SOC erosion (gC/m2)
    real(r8) :: sohg_yld                                 ! SOhg erosion (gC/m2)
    real(r8) :: zsoi_tot, zsoi_top, temp_zsoi            ! soil layer height (m)
    real(r8) :: soc, son, sohg                           ! soil OC and ON, Hg content (%)
    real(r8) :: ctot, ntot, hgtot                        ! total particulate C, N density (g/m3)
    real(r8) :: pn2poc                                   ! particulate particle N/C  ratios
    !real(r8) :: decomp_cpools_vr_col_new(1:ndecomp_pools)      ! new SOM C pools (gC/m3)
    !real(r8) :: decomp_npools_vr_col_new(1:ndecomp_pools)      ! new SOM N pools (gN/m3)
    !real(r8) :: decomp_hgpools_vr_col_new(1:ndecomp_pools)     ! new SOM hg pools (gN/m3)
    real(r8) :: dt                                       ! radiation time step (seconds)
    !real(r8), parameter :: TINYVALUE = 1.0e-14_r8  ! double precision variable has a significance of about 16 decimal digits
    character(len=32) :: subname = 'ErosionFluxes'       ! subroutine name
    !-----------------------------------------------------------------------

    associate(                                                        &
         decomp_cpools_vr    =>    soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col       , & ! Input: [real(r8) (:,:,:) ] SOM C pools [gC/m3]
         decomp_soilc_vr     =>    soilbiogeochem_carbonstate_inst%decomp_soilc_vr_col       , & ! Input: [real(r8) (:,:) ] (gC/m3) vertically-resolved decomposing total soil c pool

         decomp_npools_vr    =>    soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col      , & ! Input: [real(r8) (:,:,:) ] SOM N pools [gN/m3]
         decomp_soiln_vr     =>    soilbiogeochem_nitrogenstate_inst%decomp_soiln_vr_col       , & ! Input: [real(r8) (:,:) ] (gN/m3) vertically-resolved decomposing total soil n pool

         decomp_hgpools_vr    =>    soilbiogeochem_Hgstate_inst%decomp_hgpools_vr_col      , & ! Input: [real(r8) (:,:,:) ] SOM Hg pools [gN/m3]
         decomp_soilhg_vr     =>    soilbiogeochem_Hgstate_inst%decomp_soilHg_vr_col       , & ! Input: [real(r8) (:,:) ] (gHg/m3) vertically-resolved decomposing total soil Hg pool

         flx_sed_ero      =>    sedflux_vars%sed_ero_col      , & ! Input: [real(r8) (:) ] sed total detachment (kg/m2/s)
         flx_sed_yld      =>    sedflux_vars%sed_yld_col      , & ! Input: [real(r8) (:) ] sed total yield (kg/m2/s)
         
         bd               =>    soilstate_inst%bd_col         , & ! Input: [real(r8) (:,:) ] soil bulk density (kg/m3)

         cpools_erode     =>    soilbiogeochem_carbonflux_inst%decomp_cpools_erode    , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing C detachment (gC/m2/s)
         cpools_deposit   =>    soilbiogeochem_carbonflux_inst%decomp_cpools_deposit  , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing C redeposition (gC/m2/s)
         cpools_yield_vr  =>    soilbiogeochem_carbonflux_inst%decomp_cpools_yield_vr , & ! Output: [real(r8) (:,:,:) ] vertically-resolved decomposing C loss (gC/m3/s)

         somc_yield       =>    soilbiogeochem_carbonflux_inst%somc_yield ,             & ! Output: [real(r8) (:) ] total SOM C loss (gC/m^2/s)

         npools_erode     =>    soilbiogeochem_nitrogenflux_inst%decomp_npools_erode    , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing N detachment (gN/m2/s)
         npools_deposit   =>    soilbiogeochem_nitrogenflux_inst%decomp_npools_deposit  , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing N redeposition (gN/m2/s)
         npools_yield_vr  =>    soilbiogeochem_nitrogenflux_inst%decomp_npools_yield_vr  , & ! Output: [real(r8) (:,:,:) ] vertically-resolved decomposing N loss (gN/m3/s)
         
         hgpools_erode     =>    soilbiogeochem_Hgflux_inst%decomp_hgpools_erode    , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing hg detachment (gN/m2/s)
         hgpools_deposit   =>    soilbiogeochem_Hgflux_inst%decomp_hgpools_deposit  , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing hg redeposition (gN/m2/s)
         hgpools_yield_vr  =>    soilbiogeochem_Hgflux_inst%decomp_hgpools_yield_vr,  & ! Output: [real(r8) (:,:,:) ] vertically-resolved decomposing hg loss (gN/m3/s)
         hgpools_yield     =>    soilbiogeochem_Hgflux_inst%decomp_hgpools_yield      & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing hg loss (gC/m^2/s)
         )

         dt = real( get_step_size(), r8 )

         ! soil col filters
         do fc = 1, num_soilc
            c = filter_soilc(fc)

            ! initialization
            somc_yield(c)           = 0._r8
            cpools_erode(c,:)       = 0._r8
            cpools_deposit(c,:)     = 0._r8
            cpools_yield_vr(c,:,:)  = 0._r8
            
            npools_erode(c,:)       = 0._r8
            npools_deposit(c,:)     = 0._r8
            npools_yield_vr(c,:,:)  = 0._r8

            hgpools_erode(c,:)       = 0._r8
            hgpools_deposit(c,:)     = 0._r8
            hgpools_yield_vr(c,:,:)  = 0._r8
            hgpools_yield(c,:)       = 0._r8

            dm_ero = max(0._r8, flx_sed_ero(c)) * dt
            dm_yld = max(0._r8, flx_sed_yld(c)) * dt

            ! for Carbon
            ! assume top layer was the main source of Erosion
                  !somc_yield(c) = somc_yield(c) + max(0._r8, decomp_soilc_vr(c,1)/bd(c,1)*flx_sed_yld(c))

                  !do j = 1, nlevdecomp
                     do l = 1, ndecomp_pools
                        if ( decomp_cascade_con%is_soil(l)) then
                           cpools_yield_vr(c,1,l) = cpools_yield_vr(c,1,l) + max(0._r8, decomp_cpools_vr(c,1,l)/dzsoi_decomp(1)/bd(c,1)*flx_sed_yld(c)) !all from to 1st layer
                           !cpools_erode(c,l) = cpools_erode(c,l) + max(0._r8, cpools_yield_vr(c,j,l)* &
                              !dzsoi_decomp(j)*dm_ero/dm_yld)
                           cpools_erode(c,l) = cpools_erode(c,l) + cpools_yield_vr(c,1,l)* &
                              dzsoi_decomp(j)*dm_ero/dm_yld
                           !cpools_deposit(c,l) = cpools_deposit(c,l) + max(0._r8, cpools_erode(c,l) - &
                              !cpools_yield_vr(c,j,l)*dzsoi_decomp(j))
                           cpools_deposit(c,l) = cpools_deposit(c,l) + cpools_erode(c,l) - &
                              cpools_yield_vr(c,1,l)*dzsoi_decomp(1)
                        end if
                     end do
                  !end do

            ! for Nitrogen
            ! assume top layer was the main source of Erosion
                  !somn_yield(c) = somn_yield(c) + max(0._r8, decomp_soiln_vr(c,1)/bd(c,1)*flx_sed_yld(c))

                  !do j = 1, nlevdecomp
                     do l = 1, ndecomp_pools
                        if ( decomp_cascade_con%is_soil(l)) then
                           npools_yield_vr(c,1,l) = npools_yield_vr(c,1,l) + max(0._r8, decomp_npools_vr(c,1,l)/dzsoi_decomp(1)/bd(c,1)*flx_sed_yld(c)) !all from to 1st layer
                           !npools_erode(c,l) = npools_erode(c,l) + max(0._r8, npools_yield_vr(c,j,l)* &
                              !dzsoi_decomp(j)*dm_ero/dm_yld)
                           npools_erode(c,l) = npools_erode(c,l) + npools_yield_vr(c,1,l)* &
                              dzsoi_decomp(j)*dm_ero/dm_yld
                           !npools_deposit(c,l) = npools_deposit(c,l) + max(0._r8, npools_erode(c,l) - &
                              !npools_yield_vr(c,j,l)*dzsoi_decomp(j))
                           npools_deposit(c,l) = npools_deposit(c,l) + npools_erode(c,l) - &
                              npools_yield_vr(c,1,l)*dzsoi_decomp(1)
                        end if
                     end do
                  !end do

            ! for Hg
            ! assume top layer was the main source of Erosion
                  !somhg_yield(c) = somhg_yield(c) + max(0._r8, decomp_soilhg_vr(c,1)/bd(c,1)*flx_sed_yld(c))

                  !do j = 1, nlevdecomp
                     do l = 1, ndecomp_pools
                        if ( decomp_cascade_con%is_soil(l)) then
                           hgpools_yield_vr(c,1,l) = hgpools_yield_vr(c,1,l) + max(0._r8, decomp_hgpools_vr(c,1,l)/dzsoi_decomp(1)/bd(c,1)*flx_sed_yld(c)) !all from to 1st layer
                           !hgpools_erode(c,l) = hgpools_erode(c,l) + max(0._r8, hgpools_yield_vr(c,j,l)* &
                              !dzsoi_decomp(j)*dm_ero/dm_yld)
                           hgpools_erode(c,l) = hgpools_erode(c,l) + hgpools_yield_vr(c,1,l)* &
                              dzsoi_decomp(1)*dm_ero/dm_yld 
                           !hgpools_deposit(c,l) = hgpools_deposit(c,l) + max(0._r8, hgpools_erode(c,l) - &
                              !hgpools_yield_vr(c,j,l)*dzsoi_decomp(j))
                           hgpools_deposit(c,l) = hgpools_deposit(c,l) + hgpools_erode(c,l) - &
                              hgpools_yield_vr(c,1,l)*dzsoi_decomp(1)
                           hgpools_yield(c,l) = hgpools_yield(c,l) + hgpools_yield_vr(c,1,l)*dzsoi_decomp(1)
                        end if
                     end do
                  !end do

                  !open(unit=10, file='/data4/home/dpeng/CESM/cesm2.1.3/outputdata/clm_hg_test1/run/variables.log')
                  !write(10,*) '####'
                  !write(10,*)  'decomp_cpools_vr(c,j,l)=', decomp_cpools_vr(c,j,l)
                  !write(10,*)  'cpools_yield_vr(c,j,l)=', cpools_yield_vr(c,j,l)
                  !write(10,*)  'cpools_erode(c,l)=', cpools_erode(c,l)
                  !write(10,*)  'cpools_deposit(c,l)=', cpools_deposit(c,l)
                  !write(10,*)  'decomp_npools_vr(c,j,l)=', decomp_npools_vr(c,j,l)
                  !write(10,*)  'npools_yield_vr(c,j,l)=', npools_yield_vr(c,j,l)
                  !write(10,*)  'npools_erode(c,l)=', npools_erode(c,l)
                  !write(10,*)  'npools_deposit(c,l)=', npools_deposit(c,l)
                  !write(10,*)  'decomp_hgpools_vr(c,j,l)=', decomp_hgpools_vr(c,j,l)
                  !write(10,*)  'hgpools_yield_vr(c,j,l)=', hgpools_yield_vr(c,j,l)
                  !write(10,*)  'hgpools_erode(c,l)=', hgpools_erode(c,l)
                  !write(10,*)  'hgpools_deposit(c,l)=', hgpools_deposit(c,l)
                  !write(10,*)  'hgpools_yield(c,j)=', hgpools_yield(c,j)
                  !write(10,*)  'dzsoi_decomp(j)=', dzsoi_decomp(j)
                  !write(10,*) '####'
         end do
              
    end associate

  end subroutine ErosionFluxes

end module ErosionMod 
