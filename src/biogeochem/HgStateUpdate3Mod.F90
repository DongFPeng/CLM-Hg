module HgStateUpdate3Mod

  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! Module for Hg state variable update, mortality fluxes.
  ! Also, sminn leaching flux.
  !
  ! !USES:
  use shr_kind_mod                    , only: r8 => shr_kind_r8
  use clm_varpar                      , only: nlevdecomp, ndecomp_pools
  use clm_time_manager                , only : get_step_size
  use clm_varctl                      , only : iulog, use_nitrif_denitrif
  use clm_varpar                      , only : i_cwd, i_met_lit, i_cel_lit, i_lig_lit
  use CNVegNitrogenStateType          , only : cnveg_nitrogenstate_type
  use CNVegNitrogenFluxType           , only : cnveg_nitrogenflux_type
  use SoilBiogeochemHgStateType       , only : soilbiogeochem_Hgstate_type
  use SoilBiogeochemHgFluxType        , only : soilbiogeochem_Hgflux_type
  !
  implicit none
  private
  !
  ! !PUBLIC MEMBER FUNCTIONS:
  public:: HgStateUpdate3
  !-----------------------------------------------------------------------

contains

  !-----------------------------------------------------------------------
  subroutine HgStateUpdate3(num_soilc, filter_soilc, num_soilp, filter_soilp, &
       
       soilbiogeochem_Hgflux_inst, soilbiogeochem_Hgstate_inst)
    !
    ! !DESCRIPTION:
    ! On the radiation time step, update all the prognostic nitrogen state
    ! variables affected by gap-phase mortality fluxes. Also the Sminn leaching flux.
    ! NOTE - associate statements have been removed where there are
    ! no science equations. This increases readability and maintainability.
    !
    ! !ARGUMENTS:
    integer                                 , intent(in)    :: num_soilc       ! number of soil columns in filter
    integer                                 , intent(in)    :: filter_soilc(:) ! filter for soil columns
    integer                                 , intent(in)    :: num_soilp       ! number of soil patches in filter
    integer                                 , intent(in)    :: filter_soilp(:) ! filter for soil patches
    
    type(soilbiogeochem_hgflux_type)  , intent(in)    :: soilbiogeochem_hgflux_inst
    type(soilbiogeochem_hgstate_type) , intent(inout) :: soilbiogeochem_hgstate_inst
    !
    ! !LOCAL VARIABLES:
    integer :: c,p,j,l,k        ! indices
    integer :: fp,fc      ! lake filter indices
    real(r8):: dt         ! radiation time step (seconds)
    !-----------------------------------------------------------------------

    associate(                                         & 
        
         hgf => soilbiogeochem_hgflux_inst , & ! Input
         hgs => soilbiogeochem_hgstate_inst  & ! Output
         )

      ! set time steps
      dt = real( get_step_size(), r8 )
      !column loop
      do j = 1, nlevdecomp
         
         do fc = 1,num_soilc
            c = filter_soilc(fc)

            !!if (.not. use_nitrif_denitrif) then
            

               !!hgs%soil_hg2_vr_col(c,j) = (-hgs%soil_hg2_vr_col(c,j)) - &
               !!hgf%soil_hg2_leached_vr_col(c,j) * dt

              

           
         end do ! end of column loop
      end do

     

     

    end associate 

  end subroutine HgStateUpdate3

end module HgStateUpdate3Mod
