subroutine enthalpy_entropy_calc_sub(rho_in, enthalpy_out, entropy_out, gibbs_out)
    ! The calculation of enthalpy and entropy at given density
    ! the subroutine contain rough vectorization at input density
    ! - prepare the necesarry components such as all compresibility factor and Helmholtz energy contributions
    !	* calling only what is necesarry
    ! - Enthalpy and Entropy computation
    !
    !	== LAST MODIFICATIONS ==
    !	15.5.2017 - D.C - descriptions and comments
    !	19.5.2017 - D.C - abandon density vectorisation at this level

    use control_mod, only : dp, r_gas, tt, &
            &    stdout, stderr, stddeb, stdlog, &
            &   get_zres_fun, get_ares_fun, get_ares_dtt_fun

    implicit none
    ! variable section IN/OUT
    real(dp), intent(in) :: rho_in ! input density in assumed type as 1d array
    real(dp), intent(out) :: enthalpy_out ! output residual enthalpy in assumed type as 1d array
    real(dp), intent(out) :: entropy_out ! output residual entropy in assumed type as 1d array
    real(dp), intent(out) :: gibbs_out ! output residual gibbs energy in assumed type as 1d array
    ! local variables
    real(dp) :: z_res ! the local residual compresibility factor
    real(dp) :: a_res ! the local reduced helmholtz energy
    real(dp) :: a_res_dtt ! the local reduced helmholtz energy temperature derivative

    !------- PREPARING ------
    write (stddeb,*) "=> enthalpy_entropy_calc_sub" !DEBUG
    ! initialization of local variables
    z_res = get_zres_fun(rho_in)
    a_res = get_ares_fun(rho_in)
    a_res_dtt = get_ares_dtt_fun(rho_in)

    ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.46)
    enthalpy_out = z_res - tt * (a_res_dtt)
    enthalpy_out = enthalpy_out * r_gas * tt
    ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.48)
    entropy_out = log(1.0_dp + z_res) - tt * (-a_res)
    entropy_out = entropy_out * r_gas

    ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.49)
    gibbs_out = enthalpy_out - entropy_out * tt
    ! 	write (*,"('! gibbs energy verification:',f)"), ( (a_res) + (z_res) - log( 1.0_dp + z_res ) )*r_gas*tt !DEBUG just for verification of alternative computation

    write (stddeb,*) "<= enthalpy_entropy_calc_sub" !DEBUG
end subroutine enthalpy_entropy_calc_sub
