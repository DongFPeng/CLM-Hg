module HgInventoryMod
  !-----------------------------------------------------------------------
  ! !DESCRIPTION:
  ! module to cal Hg inventory
  !
  ! !USES:
    use shr_kind_mod      , only : r8 => shr_kind_r8
    use shr_log_mod       , only : errMsg => shr_log_errMsg
    use clm_varcon        , only : spval
    use decompMod         , only : bounds_type
    use ColumnType        , only : col
    use HgInventoryType   , only : HgInvent_type

    !
    implicit none
    save
    private
    !  
  ! !PUBLIC MEMBER FUNCTIONS:
  public :: HgEmissionProcesses        ! Calculate HgEmissionProcesses
  ! !MODULE CONSTANTS
  !-----------------------------------------------------------------------

contains

    subroutine HgEmissionProcesses(bounds, HgInvent_inst)
         !
    ! !DESCRIPTION:
    ! cal Hg emission from grid to column format
    !
    ! !ARGUMENTS:
        type(bounds_type)   , intent(in)    :: bounds
        type(HgInvent_type) , intent(inout) :: HgInvent_inst

    !
    ! !LOCAL VARIABLES:
        integer :: c,g,f ! indices
        character(len=32) :: subname = 'HgEmissionProcesses'           ! subroutine name
    !---------
        associate(                                                                      & 
            Hg_emission_grc    => HgInvent_inst%Hg_emission_grc            , &   ! Output: [real(r8) (:) ] grid Hg emission (g/km2/yr)
            Hg_emission_rof    => HgInvent_inst%Hg_emission_rof            , &   ! Output: [real(r8) (:) ] col Hg emission to rof (g/km2/yr)   
            Hg_emission_col    => HgInvent_inst%Hg_emission_col            &     ! Output: [real(r8) (:) ] col Hg emission (ng/m2/s)
            )
        !------------------------

            do c = bounds%begc, bounds%endc
               g = col%gridcell(c)
                    !initialization
                    Hg_emission_rof(c) = 0.0_r8
                    !grid to column
                    Hg_emission_col(c) = Hg_emission_grc(g)
                    Hg_emission_rof(c) = Hg_emission_rof(c) + Hg_emission_col(c) ! unit g/km2/yr
           end do

           do c = bounds%begc, bounds%endc
                f = 31536.0_r8
                Hg_emission_rof(c) = Hg_emission_rof(c)/f    ! unit g/km2/yr -> ng/m2/s

                ! open(unit=10, file='/data4/home/dpeng/waterloo/cesm2.1.3/outputdata/log/log_1.log')
                ! write(10,*) '####'
                ! write(10,*) 'Hg_emission_rof(c)_before =', Hg_emission_rof(c)
                ! write(10,*) '(c) =', c

                if (Hg_emission_rof(c)>1e30_r8) then
                    Hg_emission_rof(c) = 0.0_r8
                end if

                ! !open(unit=10, file='/data4/home/dpeng/waterloo/cesm2.1.3/outputdata/log/log_1.log')
                ! write(10,*) '####'
                ! write(10,*) 'Hg_emission_rof(c)_after =', Hg_emission_rof(c)
                ! write(10,*) '(c) =', c

           end do



       end associate

    end subroutine HgEmissionProcesses
end module HgInventoryMod