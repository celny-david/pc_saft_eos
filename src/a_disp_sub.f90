subroutine a_disp_sub (rho_in)
    ! subroutine for the reduced helmholtz energy contribution from disperison DISP term calculation
    ! - preapre the necesarry components such as disperion parameters a0 b0, dzeta_fun, combination rules and ab dispersion computation
    !	* calling only what is necesarry
    ! - the actual initialization
    ! - the preparation of intermediate computation variables
    !	* also includes cc1, cc2 terms used further down in computation section
    ! - computation of the dispersion Helmholtz energy contribution
    !	* in loop below compute composition or plainly the temperature derivative
    !	* calling only when necesarry
    ! - parameter update in control module
    !
    !	== LAST MODIFICATIONS ==
    !	15.5.2017 - D.C - descriptions and comments
    !	20.6.2017 - D.C - init_flag_g dependency omitment(partial) -> moved to parent control module -> replaced by single flag

    use control_mod, only : dp, n_comp, composition_der, temperature_der, pi, &
            &    stdout, stderr, stddeb, stdlog, &
            &    a_disp, a_disp_dx, a_disp_dtt, tt, &
            &    init_flag, init_flag_adisp
    use contrib_mod, only : a__disp, b__disp, a__disp_dx, b__disp_dx, &
            & dzeta, dzeta_dx, dzeta_dtt, m_mean, &
            & m2es3, m2e2s3, m2es3_dx, m2e2s3_dx, &
            & initialize_a0b0_disp_fun, initialize_dzeta_fun, initialize_combrules_fun, compute_ab_disp_fun, &
            & a0_disp, b0_disp
    use param_list_mod, only : param_m
    implicit none
    ! variable section IN/OUT
    real(dp), intent(in) :: rho_in ! the input density
    ! local variables outputed into control_mod
    real(dp) :: a_disp_       ! the dispersive helmholtz contribution
    real(dp), dimension(:), allocatable :: a_disp_dx_  ! the first composition derivative of dispersive helmholtz contribution , size [n_comp]
    real(dp) :: a_disp_dtt_ ! the first temperature derivative of dispersive helmholtz contribution
    ! local variables
    real(dp), dimension(0:3) :: dzeta_r ! the dzeta vector dzeta_n for n=0,1,2,3 * rho !!
    real(dp), dimension(:, :), allocatable :: dzeta_dx_r ! the dzeta composition derivative vector dzeta_dx_r for n=0,1,2,3 * rho !!, size [0:3,1:n_comp]
    real(dp), dimension(0:3) :: dzeta_dtt_r ! the dzeta temperature derivative vector dzeta_dtt_r for n=0,1,2,3 * rho !!
    real(dp) :: ii1, ii2 ! the approximated integrals computed as power series
    real(dp), dimension(:), allocatable :: ii1_dx, ii2_dx ! the approximated integrals first composition derivatives, size [n_comp]
    real(dp) :: ii1_dtt, ii2_dtt ! the approximated integrals first temperature derivatives
    real(dp) :: cc1, cc2 ! the C1 C2 abbreiations used in a_disp computation
    real(dp), dimension(:), allocatable :: cc1_dx ! the C1 first composition derivative (cc2_dx is not required) , size [n_comp]
    integer :: i, j ! counters

    !------- PREPARING ------
    write (stddeb,*) "==> a_disp_sub" !DEBUG
    ! call for the required helper procedures
    if (init_flag(1) .eqv. .false.) then
        ! 		write (stddeb,*), 'Call the initialize_a0b0_disp_fun from a_disp_sub' ! DEBUG control of computation
        if (initialize_a0b0_disp_fun() .eqv. .false.) then
            write (stderr, *) 'Error in helper function initialize_a0b0_disp_fun'
            stop
        end if
    end if
    if (init_flag(2) .eqv. .false.) then
        ! 		write (stddeb,*), 'Call the initialize_dzeta_fun from a_disp_sub' ! DEBUG control of computation
        if (initialize_dzeta_fun() .eqv. .false.) then
            write (stderr, *) 'Error in helper function initialize_dzeta_fun'
            stop
        end if
    end if
    if (init_flag(3) .eqv. .false.) then
        ! 		write (stddeb,*), 'Call the initialize_combrules_fun from a_disp_sub' ! DEBUG control of computation
        if (initialize_combrules_fun() .eqv. .false.) then
            write (stderr, *) 'Error in helper function initialize_combrules_fun'
            stop
        end if
    end if
    if (init_flag(4) .eqv. .false.) then
        ! 		write (stddeb,*), 'Call the compute_ab_disp_fun from a_disp_sub' ! DEBUG control of computation
        if (compute_ab_disp_fun() .eqv. .false.) then
            write (stderr, *) 'Error in helper function compute_ab_disp_fun'
            stop
        end if
    end if

    ! ----- allocation section -----
    allocate(a_disp_dx_(n_comp))
    allocate(dzeta_dx_r(0:3, 1:n_comp))
    allocate(ii1_dx(n_comp))
    allocate(ii2_dx(n_comp))
    allocate(cc1_dx(n_comp))
    ! the actual dzeta (multiplied by density)
    dzeta_r = dzeta * rho_in
    dzeta_dx_r = dzeta_dx * rho_in
    dzeta_dtt_r = dzeta_dtt * rho_in

    ! initialization
    ii1 = 0.0_dp
    ii1_dx = 0.0_dp
    ii1_dtt = 0.0_dp

    ii2 = 0.0_dp
    ii2_dx = 0.0_dp
    ii2_dtt = 0.0_dp

    write (stddeb,*) '----- DISP -----'  ! DEBUG just for brevity of debug
    write (stddeb,*) "param_m : ", param_m
    write (stddeb,*) "a0_disp : ", a0_disp
    write (stddeb,*) "b0_disp : ", b0_disp
    write (stddeb,*) "a0__disp: ", a__disp
    write (stddeb,*) "b0__disp: ", b__disp
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug
    write (stddeb,*) "m2es3_dx  : ", m2es3_dx
    write (stddeb,*) "m2e2s3_dx : ", m2e2s3_dx
    write (stddeb,*) "a__disp_dx: ", a__disp_dx
    write (stddeb,*) "b__disp_dx: ", b__disp_dx
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    ! ----- PREPARING INTERMEDIATE VARIABLES -----
    ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.16, A.17)
    do i = 1, 7 ! TODO modify for speed if need arises
        ! derived by maple with the property d/deta = rho/dzeta_r(3) *d/drho
        ii1 = ii1 + a__disp(i) * dzeta_r(3)**(i - 1.0_dp) ! A.16
        ii2 = ii2 + b__disp(i) * dzeta_r(3)**(i - 1.0_dp) ! A.17

        ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.42, A.43)
        if (composition_der .eqv. .true.) then
            do j = 1, n_comp
                ii1_dx(j) = ii1_dx(j) + a__disp_dx(i, j) * dzeta_r(3)**(i - 1.0_dp) &
                        & + a__disp(i) * (i - 1.0_dp) * dzeta_dx_r(3, j) * dzeta_r(3)**(i - 2.0_dp)
                ii2_dx(j) = ii2_dx(j) + b__disp_dx(i, j) * dzeta_r(3)**(i - 1.0_dp) &
                        & + b__disp(i) * (i - 1.0_dp) * dzeta_dx_r(3, j) * dzeta_r(3)**(i - 2.0_dp)
            end do
        end if
        ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.59, A.60)
        if (temperature_der .eqv. .true.) then
            ii1_dtt = ii1_dtt + a__disp(i) * (i - 1.0_dp) * dzeta_dtt_r(3) * dzeta_r(3)**(i - 2.0_dp)
            ii2_dtt = ii2_dtt + b__disp(i) * (i - 1.0_dp) * dzeta_dtt_r(3) * dzeta_r(3)**(i - 2.0_dp)
        end if
    end do

    write (stddeb,*) "ii1_dx: ", ii1_dx ! DEBUG just for brevity of debug
    write (stddeb,*) "ii2_dx: ", ii2_dx ! DEBUG just for brevity of debug
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    ! ----- C1 C2 VV implementation -----
    ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.11) RHS
    ! used the VV implementation
    !	c1_con= 1/(1 + m_mean*(8*dense-2*dense^2)/zms^4 +
    !			+ (1-m_mean)*(20*dense-27*dense^2 + 12*dense^3-2*dense^4)/(zms*(2-dense))^2); ! VV code reference
    cc1 = 1.0_dp / (1.0_dp + m_mean * (8.0_dp * dzeta_r(3)&
            & - 2.0_dp * dzeta_r(3)**2) / (1.0_dp - dzeta_r(3))**4 &
            & + (1.0_dp - m_mean) * (20.0_dp * dzeta_r(3) &
                            & - 27.0_dp * dzeta_r(3)**2 &
                            & + 12.0_dp * dzeta_r(3)**3 &
                            & - 2.0_dp * dzeta_r(3)**4) / ((1.0_dp - dzeta_r(3)) * (2.0_dp - dzeta_r(3)))**2)

    !	c2_con= - c1_con*c1_con*(m_mean*(-4*dense^2+20*dense+8)/zms^5 +
    !			+(1-m_mean)*(2*dense^3+12*dense^2-48*dense+40)/(zms*(2-dense))^3); ! VV code reference
    cc2 = -(cc1**2) * (m_mean * (-4.0_dp * dzeta_r(3)**2 &
                                & + 20.0_dp * dzeta_r(3) &
                                & + 8.0_dp) / (1.0_dp - dzeta_r(3))**5 &
                            & + (1.0_dp - m_mean) * (2.0_dp * dzeta_r(3)**3 &
                                & + 12.0_dp * dzeta_r(3)**2 &
                                & - 48.0_dp * dzeta_r(3) &
                                & + 40.0_dp) / ((1.0_dp - dzeta_r(3)) * (2.0_dp - dzeta_r(3)))**3)


    ! ----- COMPUTATION OF A_disp VV implementation -----
    ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.10)
    a_disp_ = -2.0_dp * pi * rho_in * ii1 * m2es3 &
            & - pi * rho_in * m_mean * cc1 * ii2 * m2e2s3

    write (stddeb,*) "disp pi    : ",pi ! DC DEBUG
    write (stddeb,*) "disp Dens  : ",rho_in ! DC DEBUG
    write (stddeb,*) "disp mmean : ",m_mean ! DC DEBUG
    write (stddeb,*) "disp C1    : ",cc1 ! DC DEBUG
    write (stddeb,*) "disp I1    : ",ii1 ! DC DEBUG
    write (stddeb,*) "disp I2    : ",ii2 ! DC DEBUG
    write (stddeb,*) "disp m2es2 : ",m2es3 ! DC DEBUG
    write (stddeb,*) "disp m2e2s2: ",m2e2s3 ! DC DEBUG
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug
    write (stddeb,*) "a_disp_    : ",a_disp_ ! DC DEBUG

    ! NOTE this implementation saves use of one DO cycle
    if (composition_der .eqv. .true.) then
        do i = 1, n_comp
            ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.41)
            cc1_dx(i) = cc2 * dzeta_dx_r(3, i) &
                    & - cc1**2 * param_m(i) &
                      & * ((8.0_dp * dzeta_r(3) &
                         & - 2.0_dp * dzeta_r(3)**2) / (1.0_dp - dzeta_r(3))**4 &
                         & - (20.0_dp * dzeta_r(3) &
                           & - 27.0_dp * dzeta_r(3)**2 &
                           & + 12.0_dp * dzeta_r(3)**3 &
                           & - 2.0_dp * dzeta_r(3)**4) / ((1.0_dp - dzeta_r(3)) * (2 - dzeta_r(3)))**2)

            ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.38)
            a_disp_dx_(i) = -2.0_dp * pi * rho_in * (ii1_dx(i) * m2es3  &
                    & + ii1 * m2es3_dx(i)) &
                    & - pi * rho_in * (m2e2s3 * (param_m(i) * cc1 * ii2 &
                                    & + m_mean * cc1_dx(i) * ii2 &
                                    & + m_mean * cc1 * ii2_dx(i))&
                                    & + m2e2s3_dx(i) * m_mean * cc1 * ii2)
        end do
        write (stddeb,*) "cc1_dx     : ", cc1_dx ! DEBUG just for brevity of debug
        write (stddeb,*) "a_disp_dx_ : ", a_disp_dx_ ! DEBUG just for brevity of debug
    end if

    if (temperature_der .eqv. .true.) then
        a_disp_dtt_ = -2.0_dp * pi * rho_in * m2es3 * (ii1_dtt - ii1 / tt) &
                & - pi * rho_in * m_mean * m2e2s3 * (dzeta_dtt_r(3) * cc2 * ii2 &
                        & + cc1 * ii2_dtt &
                        & - 2.0_dp * cc1 * ii2 / tt)
        write (stddeb,*) "a_disp_dtt_: ", a_disp_dtt_ ! DEBUG just for brevity of debug
    end if
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    a_disp = a_disp_
    if (composition_der .eqv. .true.) then
        a_disp_dx = a_disp_dx_
    end if
    if (temperature_der .eqv. .true.) then
        a_disp_dtt = a_disp_dtt_
    end if

    init_flag_adisp = .true.

    ! ----- deallocation section -----
    deallocate(a_disp_dx_)
    deallocate(dzeta_dx_r)
    deallocate(ii1_dx)
    deallocate(ii2_dx)
    deallocate(cc1_dx)

    write (stddeb,*) "<== a_disp_sub" !DEBUG
end subroutine a_disp_sub
