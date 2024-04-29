subroutine chem_pot_calc_sub(rho_in, chem_pot_out, fugacity_out)
    ! The calculation of the chemical potential and fugacities of individual components from given density
    ! the subroutine contain rough vectorization at input density
    ! - prepare the necesarry components such as all compresibility factor and Helmholtz energy contributions
    !	* calling only what is necesarry
    ! - chemical potential and fugacity computation
    ! - INTERFACE REQUIRED WHEN CALLING
    !	== LAST MODIFICATIONS ==
    !	15.5.2017 - D.C - descriptions and comments
    !	2?.8.2017 - D.C - dynamic n_comp implementation via allocatable approach

    use control_mod, only : dp, n_comp, r_gas, n_a, m2angstrom, tt, mol_rat, &
            &   stdout, stderr, stddeb, stdlog, &
            &   get_zres_fun, get_ares_fun, get_ares_dx_fun

    implicit none
    ! variable section IN/OUT
    real(dp), intent(in) :: rho_in ! input density
    real(dp), dimension(:), intent(out) :: chem_pot_out ! output chemical potential
    real(dp), dimension(:), intent(out) :: fugacity_out ! output fugacity coeffitient
    ! local variables
    integer :: i, j ! counter
    real(dp) :: z_res ! the local residual compresibility factor
    real(dp) :: a_res ! the local reduced helmholtz energy
    real(dp), dimension(:), allocatable :: a_res_dx ! the local reduced helmholtz energy composition derivative
    real(dp), dimension(:), allocatable :: chem_pot_id ! ideal contribution to the chemical potential
    ! - can be fixed dimension of number of components

    !------- PREPARING ------
    write (stddeb,*) "=> chem_pot_calc_sub" !DEBUG
    ! allocation
    allocate(a_res_dx(n_comp))
    allocate(chem_pot_id(n_comp))
    ! initialization of local variables
    z_res = get_zres_fun(rho_in)
    a_res = get_ares_fun(rho_in)
    call get_ares_dx_fun(rho_in, a_res_dx)

    ! computation according to paper Gross-Sadowski 2001 application of PC-SAFT eq (A.23, A.24)
    ! z_hc already contain the z_hs contribution, a_hc already contain the a_hs contribution
    do i = 1, n_comp
        ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.33)
        chem_pot_out(i) = a_res + z_res + a_res_dx(i)

        ! beware following incosistency with VV implementation
        ! BUG beware that the sum from A.33 is not present in VV implementation
        do j = 1, n_comp
            chem_pot_out(i) = chem_pot_out(i) - mol_rat(j) * (a_res_dx(j)) ! DEBUG added 0*...
        end do

        ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.32)	+edit
        fugacity_out(i) = exp (chem_pot_out(i) - log(1 + z_res))

        ! ideal gas chemical potential contribution - according VV implementation
        chem_pot_id(i) = log (mol_rat(i) * rho_in / n_a * m2angstrom**3)

        chem_pot_out(i) = (chem_pot_out(i) + chem_pot_id(i)) * r_gas * tt
        ! 			write (*,"('mu_',i2,':',e,' | fug_',i2,':',e,/)"), &
        ! 				&   i, chem_pot_out(i), i, fugacity_out(i)
    end do

    ! ----- deallocation section -----
    deallocate(a_res_dx)
    deallocate(chem_pot_id)

    write (stddeb,*) "<= chem_pot_calc_sub" !DEBUG
end subroutine chem_pot_calc_sub
