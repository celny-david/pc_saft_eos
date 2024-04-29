
! minimal testing example of external library usage
program test_prog
    use control_mod, only : dp, n_comp, rho, ch512_length, &
                &  n_a, m2angstrom
    use interface_mod, only : comp_opt, option, &
                & is_dense, &
                & parse_input_universal
    implicit none
    integer :: success
    integer :: ii ! counter
    real(dp) :: p, p1, p2, a, h, s, g ! output values of computational procedures
    real(dp), dimension(:), allocatable :: mu, fug ! output vectors of computational procedures
    real(dp), dimension(14) :: out

    call parse_input_universal()

    write (*, "('!=====            RESULTS             =====')") ! NOTE input and output section divisor

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

    !   write (*,*), '9++++++++++++++++++'
    if (comp_opt(6)) then
        call helmholtz_calc_sub(rho, p, mu, a)
        write (*, "('helmholtz energy = ',e16.8)") a
    end if
    deallocate(mu, fug)

    !   write (*,*), '10++++++++++++++++++'
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
    !   write (*,*), '11+++++++++++++++++'
    if (option/=0) then
        call num_diff_calc_sub(rho, option, 0.00001_dp, out)
    end if

end program test_prog