module ErosionMod 

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Calculate erosion induced soil particulate C, N and P fluxes 
  !
  use shr_kind_mod      , only : r8 => shr_kind_r8
  use shr_log_mod       , only : errMsg => shr_log_errMsg
  use decompMod         , only : bounds_type
  use clm_time_manager  , only : get_step_size
  use clm_varcon        , only : dzsoi_decomp
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
  !
  ! !PUBLIC TYPES:
  implicit none
  save
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: ErosionFluxes           ! Calculate erosion introduced CNP fluxes 
  ! !MODULE CONSTANTS
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine ErosionFluxes (bounds, num_soilc, filter_soilc, &
    soilstate_inst, sedflux_vars, cnveg_carbonstate_inst, cnveg_carbonflux_inst,cnveg_nitrogenstate_inst,   &
                                         soilbiogeochem_carbonflux_inst, soilbiogeochem_carbonstate_inst,                          &              
                                         cnveg_nitrogenflux_inst, soilbiogeochem_nitrogenstate_inst,                      &
                                         soilbiogeochem_nitrogenflux_inst)
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
    
    !
    ! !LOCAL VARIABLES:
    integer  :: c, fc, l, j                              ! indices
    integer  :: kc, kn                                   ! new layer indices
    integer  :: k1, k2                                   ! new layer indices     
    real(r8) :: dh, dhn                                  ! erosion height (m)
    real(r8) :: dm_ero, dm_yld, dm                       ! erosion mass (kg/m2)
    real(r8) :: soc_yld                                  ! SOC erosion (gC/m2)
    real(r8) :: zsoi_tot, zsoi_top                       ! soil layer height (m)
    real(r8) :: soc, son                                 ! soil OC and ON content (%)
    real(r8) :: ctot, ntot                        ! total particulate C, N density (g/m3)
    real(r8) :: pn2poc                          ! particulate particle N/C  ratios
    real(r8) :: decomp_cpools_vr_col_new(1:ndecomp_pools)    ! new SOM C pools (gC/m3)
    real(r8) :: decomp_npools_vr_col_new(1:ndecomp_pools)    ! new SOM N pools (gN/m3)
    real(r8) :: dt                                       ! radiation time step (seconds)
    character(len=32) :: subname = 'ErosionFluxes'       ! subroutine name
    !-----------------------------------------------------------------------

    associate(                                                        &
         cpools_vr_col    =>    soilbiogeochem_carbonstate_inst%decomp_cpools_vr_col       , & ! Input: [real(r8) (:,:,:) ] SOM C pools [gC/m3]
         npools_vr_col    =>    soilbiogeochem_nitrogenstate_inst%decomp_npools_vr_col      , & ! Input: [real(r8) (:,:,:) ] SOM N pools [gN/m3]

         flx_sed_ero      =>    sedflux_vars%sed_ero_col      , & ! Input: [real(r8) (:) ] sed total detachment (kg/m2/s)
         flx_sed_yld      =>    sedflux_vars%sed_yld_col      , & ! Input: [real(r8) (:) ] sed total yield (kg/m2/s)
         
         bd               =>    soilstate_inst%bd_col         , & ! Input: [real(r8) (:,:) ] soil bulk density (kg/m3)

         cpools_erode     =>    soilbiogeochem_carbonflux_inst%decomp_cpools_erode    , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing C detachment (gC/m2/s)
         cpools_deposit   =>    soilbiogeochem_carbonflux_inst%decomp_cpools_deposit  , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing C redeposition (gC/m2/s)
         cpools_yield_vr  =>    soilbiogeochem_carbonflux_inst%decomp_cpools_yield_vr , & ! Output: [real(r8) (:,:,:) ] vertically-resolved decomposing C loss (gC/m3/s)
                  
         npools_erode     =>    soilbiogeochem_nitrogenflux_inst%decomp_npools_erode    , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing N detachment (gN/m2/s)
         npools_deposit   =>    soilbiogeochem_nitrogenflux_inst%decomp_npools_deposit  , & ! Output: [real(r8) (:,:) ] vertically-integrated decomposing N redeposition (gN/m2/s)
         npools_yield_vr  =>    soilbiogeochem_nitrogenflux_inst%decomp_npools_yield_vr  & ! Output: [real(r8) (:,:,:) ] vertically-resolved decomposing N loss (gN/m3/s)
         )

         dt = real( get_step_size(), r8 )

         ! soil col filters
         do fc = 1, num_soilc
            c = filter_soilc(fc)

            ! initialization

            cpools_erode(c,:)       = 0._r8
            cpools_deposit(c,:)     = 0._r8
            cpools_yield_vr(c,:,:)  = 0._r8
            
            npools_erode(c,:)       = 0._r8
            npools_deposit(c,:)     = 0._r8
            npools_yield_vr(c,:,:)  = 0._r8
              

            dm_ero = max(0._r8, flx_sed_ero(c)) * dt
            dm_yld = max(0._r8, flx_sed_yld(c)) * dt
            ! calculate erosion heights 
            dh = 0._r8
            dm = dm_yld
            soc_yld = 0._r8
            kc = 1
            do j = 1, nlevdecomp
               ctot = 0._r8
               do l = 1, ndecomp_pools
                  if ( decomp_cascade_con%is_soil(l) .and. all(cpools_vr_col(c,:,l)>=0) ) then
                     ctot = ctot + cpools_vr_col(c,j,l)
                  end if
               end do
               if (dm<=bd(c,j)*dzsoi_decomp(j)) then
                  dh = dh + dm/bd(c,j) 
                  soc_yld = soc_yld + ctot*dm/bd(c,j)
                  exit
               end if
               dh = dh + dzsoi_decomp(j) 
               dm = dm - bd(c,j)*dzsoi_decomp(j)
               soc_yld = soc_yld + ctot*dzsoi_decomp(j)
               kc = kc + 1
            end do

            ! SOC content in eroded soil (only mud in surface)
            if (dh>0._r8 .and. dh<=dzsoi_decomp(1)) then
               soc = 0.1_r8 * soc_yld / 1460._r8 / dh
            else if (dh>dzsoi_decomp(1)) then
               soc = 0.1_r8 * soc_yld / dm_yld
            else
               soc = 0._r8
            end if
            son = 0.116_r8 * soc - 0.019_r8     ! Beusen et al. (2005)
            if (soc<=0._r8 .and. son<=0._r8) then
               pn2poc = 1._r8 / 10.2_r8         ! Beusen et al. (2005)
            end if

            ! equivalent erosion heights for N
            dm = pn2poc * soc_yld
            dhn = 0._r8
            kn = 1
            do j = 1, nlevdecomp
               ntot = 0._r8
               do l = 1, ndecomp_pools
                  if ( decomp_cascade_con%is_soil(l) .and. all(npools_vr_col(c,:,l)>=0._r8) ) then
                     ntot = ntot + npools_vr_col(c,j,l)   
                  end if
               end do
               if (dm<=ntot*dzsoi_decomp(j) .and. ntot>0._r8) then
                  open(unit=10, file='final_v2.log' )
                  write(10,*) 'dhn = ', dhn, 'dm = ', dm, 'ntot =', ntot
                  dhn = dhn + dm/ntot
                  exit
               end if
               dhn = dhn + dzsoi_decomp(j)
               dm = dm - ntot*dzsoi_decomp(j)
               kn = kn + 1
            end do


            ! particulate C loss
            if ( soc_yld > 0._r8 .and. kc <= nlevdecomp) then
               k1 = kc
               do j = 1, nlevdecomp
                  ! update bottom index
                  if (k1<=nlevdecomp) then
                     zsoi_tot = sum(dzsoi_decomp(1:k1))
                     zsoi_top = sum(dzsoi_decomp(1:j-1)) + dh
                     if (zsoi_top+dzsoi_decomp(j)<zsoi_tot) then
                        k2 = k1
                     else
                        k2 = k1 + 1
                     end if
                  end if
                  if (k1==k2) then
                     do l = 1, ndecomp_pools
                        if ( all(cpools_vr_col(c,:,l)>=0._r8) ) then
                           if (k1<=nlevdecomp) then
                              decomp_cpools_vr_col_new(l) = cpools_vr_col(c,k1,l)
                           else
                              decomp_cpools_vr_col_new(l) = 0._r8
                           end if
                        else
                           decomp_cpools_vr_col_new(l) = cpools_vr_col(c,j,l)
                        end if
                     end do
                  else
                     do l = 1, ndecomp_pools
                        if ( all(cpools_vr_col(c,:,l)>=0._r8) ) then
                           if (k2<=nlevdecomp) then
                              decomp_cpools_vr_col_new(l) = cpools_vr_col(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j) + &
                                 cpools_vr_col(c,k2,l)*(1._r8-(zsoi_tot-zsoi_top)/dzsoi_decomp(j)) 
                           else
                              decomp_cpools_vr_col_new(l) = cpools_vr_col(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j)
                           end if
                        else
                           decomp_cpools_vr_col_new(l) = cpools_vr_col(c,j,l)
                        end if
                     end do
                  end if
                  do l = 1, ndecomp_pools
                     if ( decomp_cascade_con%is_soil(l) ) then
                        cpools_yield_vr(c,j,l) = (cpools_vr_col(c,j,l)-decomp_cpools_vr_col_new(l))/dt
                        cpools_erode(c,l) = cpools_erode(c,l) + cpools_yield_vr(c,j,l)* &
                           dzsoi_decomp(j)*dm_ero/dm_yld
                        cpools_deposit(c,l) = cpools_deposit(c,l) + cpools_erode(c,l) - &
                           cpools_yield_vr(c,j,l)*dzsoi_decomp(j)
                     end if
                  end do
                  k1 = k2     ! update top index
               end do
            end if

            ! particulate N loss
            if ( soc_yld > 0._r8 .and. kn <= nlevdecomp ) then
               k1 = kn
               do j = 1, nlevdecomp 
                  ! update bottom index
                  if (k1<=nlevdecomp) then
                     zsoi_tot = sum(dzsoi_decomp(1:k1))
                     zsoi_top = sum(dzsoi_decomp(1:j-1)) + dhn
                     if (zsoi_top+dzsoi_decomp(j)<zsoi_tot) then
                        k2 = k1
                     else
                        k2 = k1 + 1
                     end if
                  end if
                  if (k1==k2) then
                     do l = 1, ndecomp_pools
                        if ( all(npools_vr_col(c,:,l)>=0._r8) ) then
                           if (k1<=nlevdecomp) then
                              decomp_npools_vr_col_new(l) = npools_vr_col(c,k1,l)
                           else
                              decomp_npools_vr_col_new(l) = 0._r8
                           end if
                        else
                           decomp_npools_vr_col_new(l) = npools_vr_col(c,j,l)
                        end if
                     end do
                  else
                     do l = 1, ndecomp_pools
                        if ( all(npools_vr_col(c,:,l)>=0._r8) ) then
                           if (k2<=nlevdecomp) then
                              decomp_npools_vr_col_new(l) = npools_vr_col(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j) + & 
                                 npools_vr_col(c,k2,l)*(1._r8-(zsoi_tot-zsoi_top)/dzsoi_decomp(j))
                           else
                              decomp_npools_vr_col_new(l) = npools_vr_col(c,k1,l)*(zsoi_tot-zsoi_top)/dzsoi_decomp(j)
                           end if
                        else
                           decomp_npools_vr_col_new(l) = npools_vr_col(c,j,l)
                        end if
                     end do
                  end if
                  do l = 1, ndecomp_pools
                     if ( decomp_cascade_con%is_soil(l) ) then
                        npools_yield_vr(c,j,l) = (npools_vr_col(c,j,l)-decomp_npools_vr_col_new(l))/dt 
                        npools_erode(c,l) = npools_erode(c,l) + npools_yield_vr(c,j,l)* &
                           dzsoi_decomp(j)*dm_ero/dm_yld 
                        npools_deposit(c,l) = npools_deposit(c,l) + npools_erode(c,l) - &
                           npools_yield_vr(c,j,l)*dzsoi_decomp(j)
                     end if
                  end do
                  k1 = k2     ! update top index
               end do
            end if


         end do
    
    end associate

  end subroutine ErosionFluxes

end module ErosionMod 
