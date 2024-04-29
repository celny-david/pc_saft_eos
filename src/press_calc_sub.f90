subroutine press_calc_sub(rho_in, p_out, p_drho_out, p_drho2_out)
    ! The calculation of the pressure and pressure derivative(s)
    ! the subroutine contain rough vectorization at input density
    ! - prepare the necesarry components such as all compresibility factor contributions
    !   * calling only what is necesarry
    ! - chemical potential and fugacity computation
    ! - pressure computation and first(second) density derivative of pressure
    !   * calling only what is necesarry
    !
    !   == LAST MODIFICATIONS ==
    !   15.5.2017 - D.C - descriptions and comments
    !   19.5.2017 - D.C - abandon density vectorisation at this level

    use control_mod, only : dp, n_a, k_b, tt, m2angstrom, &
            &  stdout, stderr, stddeb, stdlog, &
            &  sec_der, &
            &  get_zres_fun, get_zres_drho_fun, get_zres_drho2_fun

    implicit none
    ! variable section IN/OUT
    real(dp), intent(in) :: rho_in ! input density number
    real(dp), intent(out) :: p_out ! output pressure
    real(dp), intent(out) :: p_drho_out ! output pressure density derivative
    real(dp), intent(out) :: p_drho2_out ! output pressure second density derivative
    ! local variables
    real(dp) :: z_res, z_res_drho, z_res_drho2 ! the local residual and density derivative(s) of residual compresibility factor

    !------- PREPARING ------
    write (stddeb,*) "=> press_calc_sub" !DEBUG
    ! initialization of return variables to error output
    p_out = -1.0;
    p_drho_out = -1.0;
    p_drho2_out = -1.0;
    ! initialization of local variables
    z_res = get_zres_fun(rho_in)
    z_res_drho = get_zres_drho_fun(rho_in)

    ! computation according to paper Gross-Sadowski 2001 application of PC-SAFT eq (A.23, A.24)
    ! z_hc already contain the z_hs contribution
    p_out = (1.0 + z_res) * rho_in * k_b * tt * m2angstrom**3 ! z * rho * Kb * T * 1e30
    p_drho_out = ((1.0 + z_res) + (z_res_drho) * rho_in) * k_b * tt * m2angstrom**3 ! z_res_drho~ * Kb * T * 1e30
    p_drho_out = p_drho_out * n_a / m2angstrom**3! correction for the number density derivation [Pa*m**3/mol]
            ! write (*, "('P:',e18.6,'/','P_drho:',e18.6,' Correction term:',e18.6,'/')") &
            !         &   p_out, p_drho_out, n_a/m2angstrom**3
    
    if (sec_der .eqv. .true.) then
        z_res_drho2 = get_zres_drho2_fun(rho_in)

        p_drho2_out = (2.0 * (z_res_drho) + (z_res_drho2) * rho_in) * k_b * tt * m2angstrom**3 ! z_res_drho2~ * Kb * T * 1e30
        p_drho2_out = p_drho2_out * (n_a / m2angstrom**3)**2 ! correction for the number density derivation [Pa*m**6/mol**2]
                    ! write (*, "('P_drho2:',e18.6,'/')") &
                    !     &   p_drho2_out
    end if

    write (stddeb,*) "<= press_calc_sub" !DEBUG
end subroutine press_calc_sub
