subroutine helmholtz_calc_sub(rho_in, p_in, mu_in, helmholtz_out)
    ! The calculation of helmholtz energy (not just residual) at given density
    ! the subroutine contain rough vectorization at input density
    ! - prepare the necesarry components such as all compresibility factor and Helmholtz energy contributions
    !	* calling only what is necesarry
    ! - Helmholtz energy calculation
    !
    !	== LAST MODIFICATIONS ==
    !	18.5.2017 - D.C - implementation, descriptions and comments
    !	19.5.2017 - D.C - abandon density vectorisation at this level
    !	25.5.2017 - D.C - add pressure and chemical potential as input parameters

    use control_mod, only : dp, r_gas, n_a, m2angstrom, &
            & stdout, stderr, stddeb, stdlog, &
            & mol_rat, &
            & get_ares_fun
    !  	use contrib_mod ,only: m_mean
    implicit none

    ! variable section IN/OUT
    real(dp), intent(in) :: rho_in ! input density
    real(dp), dimension(:), intent(in) :: mu_in ! input chemical potential for sigle density
    real(dp), intent(in) :: p_in ! input pressure
    real(dp), intent(out) :: helmholtz_out ! output helmholtz energy
    ! local variables

    !------- PREPARING ------
    write (stddeb,*) "=> helmholtz_calc_sub" !DEBUG
    ! based on VV implementation + conversion of density
    ! formula is based on thermodynamical derivation A = U - TS = -PV + mu.x
    ! and fact that chemical potential contains the ideal part
    helmholtz_out = sum(mol_rat * mu_in) - p_in / rho_in * n_a / m2angstrom**3

    write (stddeb,*) "<= helmholtz_calc_sub" !DEBUG
end subroutine helmholtz_calc_sub

