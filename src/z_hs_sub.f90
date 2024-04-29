subroutine z_hs_sub (rho_in)
    ! subroutine for the compresibility factor contribution from hard sphere HS term calculation
    ! - preapre the necesarry components such as dzeta_fun
    !	* calling only when necesarry
    ! - the actual initialization
    ! - computation of the hard sphere compresibility contribution
    !	* also second derivative if required
    ! - parameter update in control module
    !
    !	== LAST MODIFICATIONS ==
    !	15.5.2017 - D.C - descriptions and comments
    !	20.6.2017 - D.C - init_flag_g dependency omitment -> moved to parent control module -> replaced by single flag module
    !	04.7.2017 - D.C - moving the z_hs contribution to get method in control mod

    use control_mod, only : dp, sec_der, &
            &   stdout, stderr, stddeb, stdlog, &
            &   write_array_1d, &
            &   a_hs, a_hs_dx, a_hs_dtt, &
            &   z_hs, z_hs_drho, z_hs_drho2, &
            &   init_flag, init_flag_zhs
    use contrib_mod, only : dzeta, m_mean, &
            & initialize_dzeta_fun
    implicit none
    ! variable section IN/OUT
    real(dp), intent(in) :: rho_in ! the input density
    ! local variables outputed into control_mod
    real(dp) :: z_hs_ ! the hard sphere compressibility contribution
    real(dp) :: z_hs_drho_ ! the hard sphere compressibility contribution
    real(dp) :: z_hs_drho2_ ! the hard sphere compressibility contribution
    ! local variables
    real(dp), dimension(0:3) :: dzeta_r ! the dzeta vector dzeta_n for n=0,1,2,3 * rho !!

    !------- PREPARING ------
    write (stddeb,*) "==> z_hs_sub" !DEBUG
    ! call for the required helper procedures
    if (init_flag(2) .eqv. .false.) then
        ! 	write (*,*), 'Call the initialize_dzeta_fun from z_hs_sub' ! DEBUG control of computation
        if (initialize_dzeta_fun() .eqv. .false.) then
            write (*, *) 'Error in helper function initialize_dzeta_fun'
            stop
        end if
    end if

    ! the actual dzeta (multiplied by density)
    dzeta_r = dzeta * rho_in

    ! DEBUG control computation
    write (stddeb,*) '-----  HS  -----'  ! DEBUG just for brevity of debug
    call write_array_1d (dzeta,"dzeta      : ",stddeb)
    call write_array_1d (dzeta_r,"dzeta_r    : ",stddeb)
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug
    ! DEBUG end of control computation

    ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.26)
    ! inconsistency within the paper eq 5 from 1999 and A.26
    z_hs_ = dzeta_r(3) / (1.0_dp - dzeta_r(3)) &
            & + (3.0_dp * dzeta_r(1) * dzeta_r(2)) / (dzeta_r(0) * (1.0_dp - dzeta_r(3))**2) &
            & + (3.0_dp * dzeta_r(2)**3 &
                    & - dzeta_r(3) * dzeta_r(2)**3) / (dzeta_r(0) * (1.0_dp - dzeta_r(3))**3)
    z_hs_drho_ = (dzeta(0) * dzeta(3) * dzeta_r(3)**2 &
            & - 3.0_dp * dzeta(1) * dzeta(2) * dzeta_r(3)**2 &
            & - 2.0_dp * dzeta_r(0) * dzeta(3)**2 &
            & + 6.0_dp * dzeta_r(2) * dzeta(2)**2 &
            & + 3.0_dp * dzeta(1) * dzeta(2)&
            & + dzeta(0) * dzeta(3)) / (dzeta(0) * (1.0_dp - dzeta_r(3))**4) ! computed and simplified by maple
    if (sec_der .eqv. .true.) then
        z_hs_drho2_ = 2.0_dp * (dzeta(0) * (dzeta_r(3) * dzeta(3))**2 &
                & - 3.0_dp * dzeta_r(1) * dzeta_r(2) * dzeta(3)**3 &
                & - 2.0_dp * dzeta_r(0) * dzeta(3)**3 &
                & - 3.0_dp * dzeta_r(1) * dzeta(2) * dzeta(3)**2 &
                & + 9.0_dp * dzeta_r(3) * dzeta(2)**3 &
                & + 6.0_dp * dzeta(1) * dzeta(2) * dzeta(3) &
                & + 3.0_dp * dzeta(2)**3 &
                & + dzeta(0) * dzeta(3)**2) / (dzeta(0) * (1.0_dp - dzeta_r(3))**5)! computed and simplified by maple
    end if

    write (stddeb,*) "z_hs_      : ",z_hs_ ! DC DEBUG
    write (stddeb,*) "z_hs_drho_ : ",z_hs_drho_ ! DC DEBUG
    if (sec_der .eqv. .true.) then
        write (stddeb,*) "z_hs_drho2_: ",z_hs_drho2_ ! DC DEBUG
    end if
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    ! output the results into control_mod for easy access
    z_hs = m_mean * z_hs_
    z_hs_drho = m_mean * z_hs_drho_
    if (sec_der .eqv. .true.) then
        z_hs_drho2 = m_mean * z_hs_drho2_
    end if

    init_flag_zhs = .true.

    write (stddeb,*) "<== z_hs_sub" !DEBUG
end subroutine z_hs_sub
