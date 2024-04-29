! ======================================================= !
! This is implementation of PC-SAFT equation of state     !
!  according to the paper of Gross and Sadowski from 2001 !
!                                                         !
!      Code was written by:     David Celný               !
!      Contact:               celny@it.cas.cz             !
!      Affiliation:                                       !
!       Institute of Thermomechanics of the CAS, v. v. i. !
!       Dolejškova 1402/5, 18200 Praha 8, Czech Republic  !
! ------------------------------------------------------- !
!      Present version:            BETA  0.2.0            !
!      Latest edit:              18.8.2017                !
! ------------------------------------------------------- !
!      CHANGES:                                           !
!         - 23.6.16 - Command line parameters input       !
!         - 22.7.16 - entropy correction + gibbs energy   !
!         - 15.5.17 - description and input parameter help!
!         - 23.5.17 - interface implementation and parser !
!         - ??.7.17 - polarity implementation             !
!         - 2?.8.17 - dynamic n_comp implementation       !
!         - 02.1.20 - list input & gfort compiler         !
! ======================================================= !


program pcsaft_main_prog
    ! the frontend where all the calls of primary methods used
    ! program accept comand line parameters in following format:
    ! spaces in the input are denoted as SPACE=' '
    ! "number of component(s)"SPACE"component(s)"SPACE"molar ratio(s)"SPACE"temperature"SPACE"density"SPACE"k_ij"SPACE"computed properties flag(s)"
    ! in brackets are the default values used when empty "" is found
    ! - component denotes computed components
    !	* in form "first component, second component2,..." i.e. "'methane','carbon dioxide'"
    !	* has to be equal to n_comp in control_mod (this will be modified)
    ! - molar ratio is always required! it is i.e. "1.0" or "0.57,0.43" ...
    ! - temperature is single positive number representing temperature in kelvins i.e. "200"
    !	* parser is capable of handling both integer 200 and float 200.0 input
    ! - density is single positive number representing computation density in mol/m^3 i.e. "997"
    !	* parser is capable of handling both integer 997 and float 997.0 input
    ! - k_ij is binary interaction coeffitient - it is single number i.e. "0"
    !	* parser is capable of handling both integer 0 and float 0.0 input
    ! - request for the computed outputs - only the requested will be printed/outputed
    !		even when other useful parameters are also computed i.e. "p"
    !	*	p	 - pressure computation
    !	*	pp	 - first density derivative of pressure
    !	*	ppp	 - second density derivative of pressure
    !	*	mu	 - chemical potentials in same order as names
    !	*	fug	 - fugacity coeffitient in same order as names
    !	*	a	 - Helmholtz energy A
    !	*	h	 - Entalpy H - residual part
    !	*	s	 - Entropy S - residual part
    !	*	g	 - Gibbs energy G - residual part
    !	*	1/2/3- control over numerical difference (if enabled in code)

    use control_mod, only : dp, n_comp, rho, ch512_length, &
            & stdout, stderr, stddeb, stdlog, init_ouput_units, &
            &  n_a, m2angstrom, &
            &  z_hs, z_hc, z_disp, z_dd, z_qq, z_dq, z_asoc, & ! DEBUG for residuals print
            &  z_hs_drho, z_hc_drho, z_disp_drho, z_dd_drho, z_qq_drho, z_dq_drho, z_asoc_drho, & ! DEBUG for residuals print
            &  z_hs_drho2, z_hc_drho2, z_disp_drho2, z_dd_drho2, z_qq_drho2, z_dq_drho2, z_asoc_drho2, & ! DEBUG for residuals print
            &  a_hs, a_hc, a_disp, a_dd, a_qq, a_dq, a_asoc  ! DEBUG for residuals print
    ! 						 &  a_hs_dtt, a_hc_dtt, a_disp_dtt, a_dd_dtt, a_qq_dtt, a_dq_dtt, a_asoc_dtt, &
    ! 						 &  a_hs_dx, a_hc_dx, a_disp_dx, a_dd_dx, a_qq_dx, a_dq_dx, a_asoc_dx &
    use interface_mod, only : comp_opt, option, &
            & is_dense, &
            & parse_input_universal
implicit none

    interface
        subroutine press_calc_sub(rho_in, p_out, p_drho_out, p_drho2_out)
            use control_mod
            real(dp), intent(in) :: rho_in ! input density
            real(dp), intent(out) :: p_out ! output pressure
            real(dp), intent(out) :: p_drho_out ! output pressure density derivative
            real(dp), intent(out) :: p_drho2_out ! output pressure second density derivative
        end subroutine press_calc_sub

        subroutine chem_pot_calc_sub(rho_in, chem_pot_out, fugacity_out)
            use control_mod
            real(dp), intent(in) :: rho_in ! input density
            real(dp), dimension(:), intent(out) :: chem_pot_out ! output chemical potential
            real(dp), dimension(:), intent(out) :: fugacity_out ! output fugacity coeffitient
        end subroutine chem_pot_calc_sub

        subroutine enthalpy_entropy_calc_sub(rho_in, enthalpy_out, entropy_out, gibbs_out)
            use control_mod
            real(dp), intent(in) :: rho_in ! input density
            real(dp), intent(out) :: enthalpy_out ! output enthalpy
            real(dp), intent(out) :: entropy_out ! output entropy
            real(dp), intent(out) :: gibbs_out ! output entropy
        end subroutine enthalpy_entropy_calc_sub

        subroutine helmholtz_calc_sub(rho_in, p_in, mu_in, helmholtz_out)
            use control_mod
            real(dp), intent(in) :: rho_in ! input density
            real(dp), dimension(:), intent(in) :: mu_in ! input chemical potential for sigle density
            real(dp), intent(in) :: p_in ! input pressure
            real(dp), intent(out) :: helmholtz_out ! output helmholtz energy
        end subroutine helmholtz_calc_sub
    end interface
    ! internal variable section
    ! OLD interface moved to the separate module
    ! 	logical, dimension(9) :: comp_opt 	! compute options that holds what user requested to be computed
    ! 										! contains: | p | dp/drho | dp/drho2 | mu | fug | A | H | S | G |
    integer :: ii ! counter
    real(dp) :: p, p1, p2, a, h, s, g ! output values of computational procedures
    real(dp), dimension(:), allocatable :: mu, fug ! output vectors of computational procedures
    real(dp), dimension(14) :: out
    ! VV modif - added for density iteration
    !real(dp) :: p_it, p1_it, p2_it
    !real(dp) :: rho_it

    call init_ouput_units()

    call parse_input_universal()

    write (*, "('!=====            RESULTS             =====')") ! NOTE input and output section divisor

    ! 	write (*,*), '7++++++++++++++++++'
    if (is_dense) then
        write (*, "('density = ',e16.8,' [mol/m**3]')") rho / n_a * m2angstrom**3
    end if

    if (comp_opt(1) .or. comp_opt(2) .or. comp_opt(3) .or. comp_opt(6)) then
        call press_calc_sub(rho, p, p1, p2)
        if (comp_opt(1)) then
            write (*, "('pressure = ',e16.8,' [Pa]')") p
        end if
        if (comp_opt(2)) then
            write (*, "('dp/drho = ',e16.8,' [Pa*m**3/mol]')") p1
        end if
        if (comp_opt(3)) then
            write (*, "('dp/drho2 = ',e16.8,' [Pa*m**6/mol**2]')") p2
        end if
    end if

    ! 	write (*,*), '8++++++++++++++++++'
    allocate(mu(n_comp), fug(n_comp))

    if (comp_opt(4).or.comp_opt(5).or.comp_opt(6)) then
        call chem_pot_calc_sub(rho, mu, fug)
        if (comp_opt(4)) then
            do ii = 1, n_comp
                write (*, "('chemical potential (',i2,') = ',e16.8)") ii, mu(ii)
            end do
        end if
        if (comp_opt(5)) then
            do ii = 1, n_comp
                write (*, "('fugacity (',i2,') = ',e16.8)") ii, fug(ii)
            end do
        end if
    end if

    ! 	write (*,*), '9++++++++++++++++++'
    if (comp_opt(6)) then
        call helmholtz_calc_sub(rho, p, mu, a)
        write (*, "('helmholtz energy = ',e16.8)") a
    end if
    deallocate(mu, fug)

    ! 	write (*,*), '10++++++++++++++++++'
    if (comp_opt(7).or.comp_opt(8).or.comp_opt(9)) then
        call enthalpy_entropy_calc_sub(rho, h, s, g)
        if (comp_opt(7)) then
            write (*, "('enthalpy (res) = ',e16.8)") h
        end if
        if (comp_opt(8)) then
            write (*, "('entropy (res) = ',e16.8)") s
        end if
        if (comp_opt(9)) then
            write (*, "('gibbs energy (res) = ',e16.8)") g
        end if
    end if
    ! 	write (*,*), '11+++++++++++++++++'
    if (option/=0) then
        call num_diff_calc_sub(rho, option, 0.00001_dp, out)
    end if

    ! DEBUG the basic info output of residual properties
!    write (*, "('!=====            RESIDUALS             =====')") ! NOTE input and output section divisor
!    write(*, "('z_res = ',e16.8)") z_hs + z_hc + z_disp + z_dd + z_qq + z_dq + z_asoc ! DEBUG
!    write(*, "('z_res_drho = ',e16.8)") z_hs_drho + z_hc_drho + z_disp_drho + z_dd_drho + z_qq_drho + z_dq_drho + z_asoc_drho ! DEBUG
!    write(*, "('z_res_drho2 = ',e16.8)") z_hs_drho2 + z_hc_drho2 + z_disp_drho2 + z_dd_drho2 + z_qq_drho2 + z_dq_drho2 + z_asoc_drho2 ! DEBUG
!
!    write(*, "('a_res = ',e16.8)") a_hs + a_hc + a_disp + a_dd + a_qq + a_dq + a_asoc ! DEBUG

    !read(*,*)

    close(stddeb)
    close(stdlog)
    close(stderr)

end program pcsaft_main_prog
