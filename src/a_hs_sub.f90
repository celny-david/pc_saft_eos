subroutine a_hs_sub (rho_in)
    ! function for the reduced helmholtz energy contribution from hard sphere HS term calculation
    ! - preapre the necesarry components such as dzeta_fun
    !	* calling only when necesarry
    ! - the actual initialization
    ! - computation of the hard sphere Helmholtz energy contribution
    !	* in loop below compute composition or plainly the temperature derivative
    !	* calling only when necesarry
    ! - parameter update in control module
    !
    !	== LAST MODIFICATIONS ==
    !	15.5.2017 - D.C - descriptions and comments
    !	20.6.2017 - D.C - init_flag_g dependency omitment -> moved to parent control module -> replaced by single flag
    !	04.7.2017 - D.C - moving the z_hs contribution to get method in control mod

    use control_mod, only : dp, n_comp, composition_der, temperature_der, &
            &    stdout, stderr, stddeb, stdlog, &
            &   a_hs, a_hs_dx, a_hs_dtt, &
            &    init_flag, init_flag_ahs
    use contrib_mod, only : m_mean, dzeta, dzeta_dx, dzeta_dtt, &
            & initialize_dzeta_fun, check_dzeta_sub
    use param_list_mod, only : param_m
    implicit none
    ! variable section IN/OUT
    real(dp), intent(in) :: rho_in ! the input density
    ! local variables outputed into control_mod
    real(dp) :: a_hs_ ! the hard sphere compressibility contribution
    real(dp), dimension(:), allocatable :: a_hs_dx_ ! the first composition derivative of hard sphere compressibility contribution, size [n_comp]
    real(dp) :: a_hs_dtt_ ! the first temperature derivative of hard sphere compressibility contribution
    ! local variables
    real(dp), dimension(0:3) :: dzeta_r ! the dzeta vector dzeta_n for n=0,1,2,3 * rho !!
    real(dp), dimension(:, :), allocatable :: dzeta_dx_r ! the dzeta composition derivation vector dzeta_dx_r for n=0,1,2,3 * rho !!, size [0:3,1:n_comp]
    real(dp), dimension(0:3) :: dzeta_dtt_r ! the dzeta temperature derivation vector dzeta_dtt_r for n=0,1,2,3 * rho !!
    integer :: i ! counters

    !------- PREPARING ------
    write (stddeb,*) "==> a_hs_sub" !DEBUG
    ! call for the required helper procedures
    if (init_flag(2) .eqv. .false.) then
        ! 	write (*,*), 'Call the initialize_dzeta_fun from a_hs_sub' ! DEBUG control of computation
        if (initialize_dzeta_fun() .eqv. .false.) then
            write (*, *) 'Error in helper function initialize_dzeta_fun'
            stop
        end if
    end if

    write (stddeb,*) '-----  HS  -----'  ! DEBUG just for brevity of debug
    ! ----- allocation section -----
    allocate(a_hs_dx_(n_comp))
    allocate(dzeta_dx_r(0:3, 1:n_comp))
    ! the actual dzeta (multiplied by density)
    dzeta_r = dzeta * rho_in
    dzeta_dx_r = dzeta_dx * rho_in
    dzeta_dtt_r = dzeta_dtt * rho_in

    call check_dzeta_sub(dzeta_r)
    write (stddeb,*) "dzeta_r    : ", dzeta_r
    write (stddeb,*) "dzeta_dx_r : ", dzeta_dx_r
    write (stddeb,*) "dzeta_dtt_r: ", dzeta_dtt_r
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    ! 	write (*,*), '--------- A-HS ---------'  ! DEBUG just for brevity of debug
    ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.6)
    a_hs_ = (3.0_dp * dzeta_r(1) * dzeta_r(2) / (1.0_dp - dzeta_r(3)) &
            & + dzeta_r(2)**3 / (dzeta_r(3) * (1.0_dp - dzeta_r(3))**2) &
            & + (dzeta_r(2) * (dzeta_r(2) / dzeta_r(3))**2 - dzeta_r(0)) * log(1.0_dp - dzeta_r(3))) / dzeta_r(0)

    ! 	write (*,"('test : ',e,/)"), dzeta_dx(0,1) ! DEBUG just for brevity of debug

    ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.36)
    if (composition_der .eqv. .true.) then
        do i = 1, n_comp
            a_hs_dx_(i) = -a_hs_ * dzeta_dx_r(0, i) / dzeta_r(0) &
                    & + (3.0_dp * (dzeta_dx_r(1, i) * dzeta_r(2) &
                            & + dzeta_r(1) * dzeta_dx_r(2, i)) / (1.0_dp - dzeta_r(3)) &
                            & + 3.0_dp * dzeta_r(1) * dzeta_r(2) * dzeta_dx_r(3, i) / (1.0_dp - dzeta_r(3))**2 &
                            & + 3.0_dp * dzeta_r(2)**2 * dzeta_dx_r(2, i) / (dzeta_r(3) * (1.0_dp - dzeta_r(3))**2) &
                            & + dzeta_r(2)**3 * dzeta_dx_r(3, i) * (3.0_dp * dzeta_r(3) - 1.0_dp) / (dzeta_r(3)**2 * (1.0_dp - dzeta_r(3))**3) &
                            & + (3.0_dp * dzeta_r(2)**2 * dzeta_dx_r(2, i) / dzeta_r(3)**2 &
                                    & - 2.0_dp * dzeta_r(2)**3 * dzeta_dx_r(3, i) / dzeta_r(3)**3 &
                                    & - dzeta_dx_r(0, i)) * log(1.0_dp - dzeta_r(3)) &
                            & - (dzeta_r(2) * (dzeta_r(2) / dzeta_r(3))**2 - dzeta_r(0)) * dzeta_dx_r(3, i) / (1.0_dp - dzeta_r(3))) / dzeta_r(0)
        end do
    end if
    ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.55)
    ! the shape of derivative is identical with used fact that dzeta_dx(0)=0
    if (temperature_der .eqv. .true.) then
        a_hs_dtt_ = (3.0_dp * (dzeta_dtt_r(1) * dzeta_r(2) &
                & + dzeta_r(1) * dzeta_dtt_r(2)) / (1.0_dp - dzeta_r(3)) &
                & + 3.0_dp * dzeta_r(1) * dzeta_r(2) * dzeta_dtt_r(3) / (1.0_dp - dzeta_r(3))**2 &
                & + 3.0_dp * dzeta_r(2)**2 * dzeta_dtt_r(2) / (dzeta_r(3) * (1.0_dp - dzeta_r(3))**2) &
                & - dzeta_r(2)**3 * dzeta_dtt_r(3) / (dzeta_r(3) * (1.0_dp - dzeta_r(3)))**2 &
                & + 2.0_dp * dzeta_r(2)**3 * dzeta_dtt_r(3) / (dzeta_r(3) * (1.0_dp - dzeta_r(3))**3) &
                & + (3.0_dp * dzeta_r(2)**2 * dzeta_dtt_r(2) / dzeta_r(3)**2 &
                        & - 2.0_dp * dzeta_r(2)**3 * dzeta_dtt_r(3) / dzeta_r(3)**3) * log(1.0_dp - dzeta_r(3)) &
                & - (dzeta_r(2) * (dzeta_r(2) / dzeta_r(3))**2 - dzeta_r(0)) * dzeta_dtt_r(3) / (1.0_dp - dzeta_r(3))) / dzeta_r(0)
    end if

    write (stddeb,*) "a_hs_    : ",a_hs_ ! DC DEBUG
    if (composition_der .eqv. .true.) then
        write (stddeb,*) "a_hs_dx_ : ",a_hs_dx_ ! DC DEBUG
    end if
    if (temperature_der .eqv. .true.) then
        write (stddeb,*) "a_hs_dtt_: ",a_hs_dtt_ ! DC DEBUG
    end if
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    ! output the results into control_mod for easy access
    a_hs = m_mean * a_hs_ ! BUG beware potential inconsistency in output
    !	a_hs = a_hs_ ! BUG deprechated -it require cyclic inclusion of modules 12.3.2018
    if (composition_der .eqv. .true.) then
        a_hs_dx = param_m * a_hs_ + m_mean * a_hs_dx_ ! BUG beware potential inconsistency in output
        ! 		a_hs_dx = a_hs_dx_ ! BUG deprechated -it require cyclic inclusion of modules 12.3.2018
    end if
    if (temperature_der .eqv. .true.) then
        a_hs_dtt = m_mean * a_hs_dtt_ ! BUG beware potential inconsistency in output
        ! 		a_hs_dtt = a_hs_dtt_ ! BUG deprechated -it require cyclic inclusion of modules 12.3.2018
    end if

    init_flag_ahs = .true.

    ! ----- deallocation section -----
    deallocate(a_hs_dx_)
    deallocate(dzeta_dx_r)

    write (stddeb,*) "<== a_hs_sub" !DEBUG
end subroutine a_hs_sub
