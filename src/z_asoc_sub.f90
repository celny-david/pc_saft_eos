subroutine z_asoc_sub (rho_in)
    ! subroutine for the compresibility factor contribution from dipolar-dipolar interaction
    ! - expect that necessary polar parametrs were computed beforehands
    !	* calling only what is necesarry
    !	* a,b,c and ee..., eee... parameters
    ! - initialisation and computation of the J2, J3 integral aproximation and respective derivatives
    !	* optimalized for symmetrical form of computation
    !	* symmetry correction of density derivatives with dzeta factor
    !		** this form saves multiplication
    ! - initialisation and computation of intermediate computation variables a2, a3 (derivatives)
    !	* optimalized for symmetrical form of computation
    ! - computation of compresibility factor helmholtz energy contribution and derivatives
    ! - parameter update in control module
    !	* due to the wasteful initialization of variable for temperature and composition derivative
    !	  these two are modified directly without substantial loss of debugging verbosity
    !
    !	== LAST MODIFICATIONS ==
    !	20.6.2017 - D.C - init_flag_g dependency omitment -> moved to parent control module -> replaced by single flag
    !	19.7.2017 - D.C - finished implementation also with all derivatives

    use control_mod, only : dp, n_comp, pi, tt, mol_rat, &
            &    stdout, stderr, stddeb, stdlog, &
            &    write_array_1d, write_array_2d, write_array_3d, &
            &    sec_der, temperature_der, composition_der, &
            &    init_flag, init_flag_aasoc, init_flag_zasoc, init_flag_gij, &
            &    g_ij, g_ij_drho, g_ij_drho2, g_ij_drho3, g_ij_dtt, g_ij_dx, &
            &    a_asoc, a_asoc_dtt, a_asoc_dx, &
            &    z_asoc, z_asoc_drho, z_asoc_drho2
    use contrib_mod, only : dzeta, dzeta_dtt, dzeta_dx, &
            &    sigma_ij, eps_asoc_ij, kapa_asoc_ij, &
            &    d_seg, d_seg_dtt, &
            &    initialize_dzeta_fun, initialize_combrules_fun, compute_qudratic_sub, compute_qudratic3_sub
    use param_list_mod, only : param_bond
    implicit none
    ! variable section IN/OUT
    real(dp), intent(in) :: rho_in ! the input density
    ! local variables outputed into control_mod
    real(dp), dimension(:), allocatable :: a_asoc_, z_asoc_, z_asoc_drho_, z_asoc_drho2_, a_asoc_dtt_ ! the association helmholtz contribution , size [n_comp]
    real(dp), dimension(:, :), allocatable :: a_asoc_dx_ ! the association contribution
    ! local variables
    real(dp), dimension(0:3) :: dzeta_r ! the dzeta vector dzeta_n for n=0,1,2,3 * rho !!
    real(dp), dimension(:, :), allocatable :: dzeta_dx_r ! the dzeta composition derivative vector dzeta_dx_r for n=0,1,2,3 * rho !!, size [0:3,1:n_comp]
    real(dp), dimension(0:3) :: dzeta_dtt_r ! the dzeta temperature derivative vector dzeta_dtt_r for n=0,1,2,3 * rho !!
    real(dp) :: delta_ij_tmp ! the temporarry value for delta_ij computation
    real(dp), dimension(:, :), allocatable :: delta_ij ! strength of association interaction, size [n_comp,n_comp]
    real(dp), dimension(:, :), allocatable :: delta_ij_drho ! strength of association interaction derivative, size [n_comp,n_comp]
    real(dp), dimension(:, :), allocatable :: delta_ij_drho2 ! strength of association interaction second derivative, size [n_comp,n_comp]
    real(dp), dimension(:, :), allocatable :: delta_ij_drho3 ! strength of association interaction third derivative, size [n_comp,n_comp]
    real(dp), dimension(:, :), allocatable :: delta_ij_dtt ! strength of association interaction temperature derivative, size [n_comp,n_comp]
    real(dp), dimension(:, :, :), allocatable :: delta_ij_dx ! strength of association interaction composition derivative, size [n_comp,n_comp,n_comp]
    real(dp), dimension(4) :: x_a, x_b, x_c ! quadratic equation cofitients and their first, second and third derivative
    ! the shape refers to x_a = (/a, da_drho, da_drho2, da_drho3/)
    real(dp), dimension(4) :: xx, xx_drho, xx_drho2, xx_drho3 ! the X^A parameter holding all
    real(dp), dimension(4) :: xx_dtt! the X^A temperature derivative
    real(dp), dimension(:, :), allocatable :: xx_dx! the X^A composition derivative, size [4,n_comp]
    integer, dimension(:), allocatable :: n_site ! the temporarry variable determining the number of sites, size [n_comp]
    integer :: tmp_bond ! the temporarry variable determining the type of bond from param_bond
    integer :: i, j, k ! counters

    !------- PREPARING ------
    write (stddeb,*) "==> z_assoc_sub" !DEBUG
    ! call for the required helper procedures
    if (init_flag(2) .eqv. .false.) then
        if (initialize_dzeta_fun() .eqv. .false.) then
            write (stderr, *) 'Error in helper function initialize_dzeta_fun'
            stop
        end if
    end if
    if (init_flag(3) .eqv. .false.) then
        if (initialize_combrules_fun() .eqv. .false.) then
            write (stderr, *) 'Error in helper function initialize_combrules_fun'
            stop
        end if
    end if
    ! call for the required level 1 subroutine rdf
    if (init_flag_gij .eqv. .false.) then
        call rdf_sub(rho_in)
    end if

    ! ----- allocation section -----
    allocate(dzeta_dx_r(0:3, 1:n_comp))
    allocate(a_asoc_(n_comp))
    allocate(a_asoc_dtt_(n_comp))
    allocate(a_asoc_dx_(n_comp, n_comp))
    allocate(z_asoc_(n_comp))
    allocate(z_asoc_drho_(n_comp))
    allocate(z_asoc_drho2_(n_comp))
    allocate(delta_ij(n_comp, n_comp))
    allocate(delta_ij_drho(n_comp, n_comp))
    allocate(delta_ij_drho2(n_comp, n_comp))
    allocate(delta_ij_drho3(n_comp, n_comp))
    allocate(delta_ij_dtt(n_comp, n_comp))
    allocate(delta_ij_dx(n_comp, n_comp, n_comp))
    allocate(xx_dx(4, n_comp))
    allocate(n_site(n_comp))

    ! the actual dzeta (multiplied by density)
    dzeta_r = dzeta * rho_in
    dzeta_dx_r = dzeta_dx * rho_in
    dzeta_dtt_r = dzeta_dtt * rho_in

    ! DEBUG control computation
    write (stddeb,*) '-----  QQ  -----'  ! DEBUG just for brevity of debug
    call write_array_1d (dzeta,  "dzeta      : ",stddeb)
    call write_array_1d (dzeta_r,"dzeta_r    : ",stddeb)

    call write_array_1d (dzeta_dtt,  "dzeta_dtt  : ",stddeb)
    call write_array_1d (dzeta_dtt_r,"dzeta_dtt_r: ",stddeb)

    call write_array_2d (dzeta_dx,  "dzeta_dx   : ",stddeb)
    call write_array_2d (dzeta_dx_r,"dzeta_dx_r : ",stddeb)
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug
    call write_array_2d (g_ij,    "g_ij     : ",stddeb)
    call write_array_2d (g_ij_dtt,"g_ij_dtt : ",stddeb)
    call write_array_3d (g_ij_dx, "g_ij_dx  : ",stddeb)
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug
    ! DEBUG end of control computation

    ! BUG major problem with inconsistency of number of sites expessed here via sigma_ij
    !     - i,j index sites A,B,C ...
    !	  and the RDF g_ij which express different components via g_ij(i,j)
    !	  - i,j index number of components
    !     this schema can work only when number of component match number of substances but results are still bad
    ! DEBUG - the problem with inconsistency in paper Huang-Radosz 1990 with sigma in delta_ij
    if (.true.) then ! NOTE the hadcoded lever for one or other delta_ij formulation
        ! NOTE paper Huang-Radosz 1990 version with sigma in delta_ij
        ! NOTE paper Huang-Radosz 1990 states version for pure hard-sphere liquid
        ! NOTE simple cycle version can utilise symmetry
        do i = 1, n_comp
            do j = 1, n_comp
                delta_ij_tmp = (exp(eps_asoc_ij(i, j) / tt) - 1.0_dp) * sigma_ij(i, j)**3 * kapa_asoc_ij(i, j)
                delta_ij(i, j) = g_ij(i, j) * delta_ij_tmp
                delta_ij_drho(i, j) = g_ij_drho(i, j) * delta_ij_tmp
                delta_ij_drho2(i, j) = g_ij_drho2(i, j) * delta_ij_tmp
                if (sec_der .eqv. .true.) then
                    delta_ij_drho3(i, j) = g_ij_drho3(i, j) * delta_ij_tmp
                end if
                if (temperature_der .eqv. .true.) then
                    delta_ij_dtt(i, j) = g_ij_dtt(i, j) * delta_ij_tmp  &
                            & - g_ij(i, j) * exp(eps_asoc_ij(i, j) / tt) * eps_asoc_ij(i, j) / tt**2 * sigma_ij(i, j)**3 * kapa_asoc_ij(i, j)
                end if
                if (composition_der .eqv. .true.) then
                    do k = 1, n_comp
                        delta_ij_dx(i, j, k) = g_ij_dx(i, j, k) * delta_ij_tmp ! OPTIM the third loop can be optimised
                    end do
                end if
            end do
        end do
    else
        ! NOTE paper Chapman-Gubbins-Jackson-Radosz 1990 version with d in delta_ij
        ! NOTE simple cycle version can utilise symmetry
        do i = 1, n_comp
            do j = 1, n_comp
                delta_ij_tmp = (exp(eps_asoc_ij(i, j) / tt) - 1.0_dp) * ((d_seg(i) + d_seg(j)) / 2.0_dp)**3 * kapa_asoc_ij(i, j)
                delta_ij(i, j) = g_ij(i, j) * delta_ij_tmp
                delta_ij_drho(i, j) = g_ij_drho(i, j) * delta_ij_tmp
                delta_ij_drho2(i, j) = g_ij_drho2(i, j) * delta_ij_tmp
                if (sec_der .eqv. .true.) then
                    delta_ij_drho3(i, j) = g_ij_drho3(i, j) * delta_ij_tmp
                end if
                if (temperature_der .eqv. .true.) then
                    delta_ij_dtt(i, j) = g_ij_dtt(i, j) * delta_ij_tmp  &
                            & - g_ij(i, j) * exp(eps_asoc_ij(i, j) / tt) * eps_asoc_ij(i, j) / tt**2 * ((d_seg(i) + d_seg(j)) / 2.0_dp)**3 * kapa_asoc_ij(i, j) &
                            & + 3.0_dp * g_ij(i, j) * delta_ij_tmp * (2.0_dp / (d_seg(i) + d_seg(j))) * ((d_seg_dtt(i) + d_seg_dtt(j)) / 2.0_dp)
                end if
                if (composition_der .eqv. .true.) then
                    do k = 1, n_comp
                        delta_ij_dx(i, j, k) = g_ij_dx(i, j, k) * delta_ij_tmp ! OPTIM the third loop can be optimised
                    end do
                end if
            end do
        end do
    end if

    call write_array_2d (delta_ij,      "delta_ij       :", stddeb) ! DEBUG control computation
    call write_array_2d (delta_ij_drho, "delta_ij_drho  :", stddeb) ! DEBUG control computation
    call write_array_2d (delta_ij_drho2,"delta_ij_drho2 :", stddeb) ! DEBUG control computation
    call write_array_2d (delta_ij_drho3,"delta_ij_drho3 :", stddeb) ! DEBUG control computation
    call write_array_2d (delta_ij_dtt,  "delta_ij_dtt   :", stddeb) ! DEBUG control computation
    call write_array_3d (delta_ij_dx,   "delta_ij_dx    :", stddeb) ! DEBUG control computation
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    ! the if statemets corresponds with paper Huang-Radosz 1990 table VII
    ! set up the quadratic equation parameters for xx(a) and compute it - the result is redistributed according to type of bond interaction
    do i = 1, n_comp
        write (stddeb,"('--CYCLE ',i2,' / ',i2,' --')") i, n_comp ! DEBUG just for brevity of debug

        if (kapa_asoc_ij(i, i)==0.0_dp) then ! NOTE when non-associating substance encountered then skip the computation with values set to zeros
            a_asoc_(i) = 0.0_dp
            z_asoc_(i) = 0.0_dp
            z_asoc_drho_(i) = 0.0_dp
            z_asoc_drho2_(i) = 0.0_dp
            a_asoc_dtt_(i) = 0.0_dp
            a_asoc_dx_(i, :) = 0.0_dp ! NOTE vector intialization
            cycle
        end if
        tmp_bond = int(param_bond(i))
        if (tmp_bond == 0) then
            n_site(i) = 1

            x_a = (/rho_in * delta_ij(i, 1), &
                    rho_in * delta_ij_drho(i, 1) + delta_ij(i, 1), &
                    rho_in * delta_ij_drho2(i, 1) + 2.0_dp * delta_ij_drho(i, 1), &
                    rho_in * delta_ij_drho3(i, 1) + 3.0_dp * delta_ij_drho2(i, 1)/)
            x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
            x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

            call compute_qudratic3_sub(x_a, x_b, x_c, xx(1), xx_drho(1), xx_drho2(1), xx_drho3(1)) ! compute xx, xx_drho...

            if (temperature_der .eqv. .true.) then
                x_a = (/rho_in * delta_ij(i, 1), &
                        rho_in * delta_ij_dtt(i, 1), &
                        0.0_dp, &
                        0.0_dp/)
                x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                call compute_qudratic_sub(x_a, x_b, x_c, xx_dtt(1))
            end if

            if (composition_der .eqv. .true.) then
                do j = 1, n_comp
                    x_a = (/rho_in * delta_ij(i, 1), &
                            rho_in * delta_ij_dx(i, 1, j), &
                            0.0_dp, &
                            0.0_dp/)
                    x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                    x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                    call compute_qudratic_sub(x_a, x_b, x_c, xx_dx(1, j))
                end do
            end if

        elseif (tmp_bond == 21) then
            n_site(i) = 2

            x_a = (/2.0_dp * rho_in * delta_ij(i, 1), &
                    2.0_dp * rho_in * delta_ij_drho(i, 1) + 2.0_dp * delta_ij(i, 1), &
                    2.0_dp * rho_in * delta_ij_drho2(i, 1) + 4.0_dp * delta_ij_drho(i, 1), &
                    2.0_dp * rho_in * delta_ij_drho3(i, 1) + 6.0_dp * delta_ij_drho2(i, 1)/)
            x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
            x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

            call compute_qudratic3_sub(x_a, x_b, x_c, xx(1), xx_drho(1), xx_drho2(1), xx_drho3(1)) ! compute xx, xx_drho...
            ! hold the approx equations between X^A X^B
            xx(2) = xx(1)
            xx_drho(2) = xx_drho(1)
            xx_drho2(2) = xx_drho2(1)
            xx_drho3(2) = xx_drho3(1)

            if (temperature_der .eqv. .true.) then
                x_a = (/2.0_dp * rho_in * delta_ij(i, 1), &
                        2.0_dp * rho_in * delta_ij_dtt(i, 1), &
                        0.0_dp, &
                        0.0_dp/)
                x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                call compute_qudratic_sub(x_a, x_b, x_c, xx_dtt(1))

                xx_dtt(2) = xx_dtt(1)
            end if

            if (composition_der .eqv. .true.) then
                do j = 1, n_comp
                    x_a = (/2.0_dp * rho_in * delta_ij(i, 1), &
                            2.0_dp * rho_in * delta_ij_dx(i, 1, j), &
                            0.0_dp, &
                            0.0_dp/)
                    x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                    x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                    call compute_qudratic_sub(x_a, x_b, x_c, xx_dx(1, j))

                    xx_dx(2, j) = xx_dx(1, j)
                end do
            end if

        elseif (tmp_bond == 22) then
            ! delta_11 = delta_22 = 0
            ! delta_12 <> 0
            n_site(i) = 2

            x_a = (/rho_in * delta_ij(1, 1), &
                    rho_in * delta_ij_drho(1, 1) + delta_ij(1, 1), &
                    rho_in * delta_ij_drho2(1, 1) + 2.0_dp * delta_ij_drho(1, 1), &
                    rho_in * delta_ij_drho3(1, 1) + 3.0_dp * delta_ij_drho2(1, 1)/)
            x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
            x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

            call compute_qudratic3_sub(x_a, x_b, x_c, xx(1), xx_drho(1), xx_drho2(1), xx_drho3(1)) ! compute xx, xx_drho...
            ! hold the approx equations between X^A X^B
            xx(2) = xx(1)
            xx_drho(2) = xx_drho(1)
            xx_drho2(2) = xx_drho2(1)
            xx_drho3(2) = xx_drho3(1)

            if (temperature_der .eqv. .true.) then
                x_a = (/rho_in * delta_ij(i, 1), &
                        rho_in * delta_ij_dtt(i, 1), &
                        0.0_dp, &
                        0.0_dp/)
                x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                call compute_qudratic_sub(x_a, x_b, x_c, xx_dtt(1))

                xx_dtt(2) = xx_dtt(1)
            end if

            if (composition_der .eqv. .true.) then
                do j = 1, n_comp
                    x_a = (/rho_in * delta_ij(i, 1), &
                            rho_in * delta_ij_dx(i, 1, j), &
                            0.0_dp, &
                            0.0_dp/)
                    x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                    x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                    call compute_qudratic_sub(x_a, x_b, x_c, xx_dx(1, j))

                    xx_dx(2, j) = xx_dx(1, j)
                end do
            end if

        elseif (tmp_bond == 31) then
            n_site(i) = 3

            x_a = (/3.0_dp * rho_in * delta_ij(i, 1), &
                    3.0_dp * rho_in * delta_ij_drho(i, 1) + 3.0_dp * delta_ij(i, 1), &
                    3.0_dp * rho_in * delta_ij_drho2(i, 1) + 6.0_dp * delta_ij_drho(i, 1), &
                    3.0_dp * rho_in * delta_ij_drho3(i, 1) + 9.0_dp * delta_ij_drho2(i, 1)/)
            x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
            x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

            call compute_qudratic3_sub(x_a, x_b, x_c, xx(1), xx_drho(1), xx_drho2(1), xx_drho3(1)) ! compute xx, xx_drho...
            ! hold the approx equations between X^A X^B X^C
            xx(2) = xx(1)
            xx_drho(2) = xx_drho(1)
            xx_drho2(2) = xx_drho2(1)
            xx_drho3(2) = xx_drho3(1)

            xx(3) = xx(1)
            xx_drho(3) = xx_drho(1)
            xx_drho2(3) = xx_drho2(1)
            xx_drho3(3) = xx_drho3(1)

            if (temperature_der .eqv. .true.) then
                x_a = (/3.0_dp * rho_in * delta_ij(i, 1), &
                        3.0_dp * rho_in * delta_ij_dtt(i, 1), &
                        0.0_dp, &
                        0.0_dp/)
                x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                call compute_qudratic_sub(x_a, x_b, x_c, xx_dtt(1))

                xx_dtt(2) = xx_dtt(1)
                xx_dtt(3) = xx_dtt(1)
            end if

            if (composition_der .eqv. .true.) then
                do j = 1, n_comp
                    x_a = (/3.0_dp * rho_in * delta_ij(i, 1), &
                            3.0_dp * rho_in * delta_ij_dx(i, 1, j), &
                            0.0_dp, &
                            0.0_dp/)
                    x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                    x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                    call compute_qudratic_sub(x_a, x_b, x_c, xx_dx(1, j))

                    xx_dx(2, j) = xx_dx(1, j)
                    xx_dx(3, j) = xx_dx(1, j)
                end do
            end if

        elseif (tmp_bond == 32) then
            n_site(i) = 3

            ! these are not same as VV coefitients
            x_a = (/2.0_dp * rho_in * delta_ij(i, 1), &
                    2.0_dp * rho_in * delta_ij_drho(i, 1) + 2.0_dp * delta_ij(i, 1), &
                    2.0_dp * rho_in * delta_ij_drho2(i, 1) + 4.0_dp * delta_ij_drho(i, 1), &
                    2.0_dp * rho_in * delta_ij_drho3(i, 1) + 6.0_dp * delta_ij_drho2(i, 1)/)
            x_b = (/1.0_dp - rho_in * delta_ij(i, 1), &
                    -rho_in * delta_ij_drho(i, 1) - delta_ij(i, 1), &
                    -rho_in * delta_ij_drho2(i, 1) - 2.0_dp * delta_ij_drho(i, 1), &
                    -rho_in * delta_ij_drho3(i, 1) - 3.0_dp * delta_ij_drho2(i, 1)/)
            x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

            call compute_qudratic3_sub(x_a, x_b, x_c, xx(1), xx_drho(1), xx_drho2(1), xx_drho3(1)) ! compute xx, xx_drho...
            ! hold the approx equations between X^A X^B X^C
            xx(2) = xx(1)
            xx_drho(2) = xx_drho(1)
            xx_drho2(2) = xx_drho2(1)
            xx_drho3(2) = xx_drho3(1)

            xx(3) = 2.0_dp * xx(1) - 1.0_dp
            xx_drho(3) = 2.0_dp * xx_drho(1)
            xx_drho2(3) = 2.0_dp * xx_drho2(1)
            xx_drho3(3) = 2.0_dp * xx_drho3(1)

            if (temperature_der .eqv. .true.) then
                x_a = (/2.0_dp * rho_in * delta_ij(i, 1), &
                        2.0_dp * rho_in * delta_ij_dtt(i, 1), &
                        0.0_dp, &
                        0.0_dp/)
                x_b = (/1.0_dp - rho_in * delta_ij(i, 1), &
                        -rho_in * delta_ij_dtt(i, 1), &
                        0.0_dp, &
                        0.0_dp/)
                x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                call compute_qudratic_sub(x_a, x_b, x_c, xx_dtt(1))

                xx_dtt(2) = xx_dtt(1)
                xx_dtt(3) = 2.0_dp * xx_dtt(1)
            end if

            if (composition_der .eqv. .true.) then
                do j = 1, n_comp
                    x_a = (/2.0_dp * rho_in * delta_ij(i, 1), &
                            2.0_dp * rho_in * delta_ij_dx(i, 1, j), &
                            0.0_dp, &
                            0.0_dp/)
                    x_b = (/1.0_dp - rho_in * delta_ij(i, 1), &
                            -rho_in * delta_ij_dx(i, 1, j), &
                            0.0_dp, &
                            0.0_dp/)
                    x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                    call compute_qudratic_sub(x_a, x_b, x_c, xx_dx(1, j))

                    xx_dx(2, j) = xx_dx(1, j)
                    xx_dx(3, j) = 2.0_dp * xx_dx(1, j)
                end do
            end if

        elseif (tmp_bond == 41) then
            n_site(i) = 4

            x_a = (/4.0_dp * rho_in * delta_ij(i, 1), &
                    4.0_dp * rho_in * delta_ij_drho(i, 1) + 4.0_dp * delta_ij(i, 1), &
                    4.0_dp * rho_in * delta_ij_drho2(i, 1) + 8.0_dp * delta_ij_drho(i, 1), &
                    4.0_dp * rho_in * delta_ij_drho3(i, 1) + 12.0_dp * delta_ij_drho2(i, 1)/)
            x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
            x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

            call compute_qudratic3_sub(x_a, x_b, x_c, xx(1), xx_drho(1), xx_drho2(1), xx_drho3(1)) ! compute xx, xx_drho...
            ! hold the approx equations between X^A X^B X^C X^D
            xx(2) = xx(1)
            xx_drho(2) = xx_drho(1)
            xx_drho2(2) = xx_drho2(1)
            xx_drho3(2) = xx_drho3(1)

            xx(3) = xx(1)
            xx_drho(3) = xx_drho(1)
            xx_drho2(3) = xx_drho2(1)
            xx_drho3(3) = xx_drho3(1)

            xx(4) = xx(1)
            xx_drho(4) = xx_drho(1)
            xx_drho2(4) = xx_drho2(1)
            xx_drho3(4) = xx_drho3(1)

            if (temperature_der .eqv. .true.) then
                x_a = (/4.0_dp * rho_in * delta_ij(i, 1), &
                        4.0_dp * rho_in * delta_ij_dtt(i, 1), &
                        0.0_dp, &
                        0.0_dp/)
                x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                call compute_qudratic_sub(x_a, x_b, x_c, xx_dtt(1))

                xx_dtt(2) = xx_dtt(1)
                xx_dtt(3) = xx_dtt(1)
                xx_dtt(4) = xx_dtt(1)
            end if

            if (composition_der .eqv. .true.) then
                do j = 1, n_comp
                    x_a = (/4.0_dp * rho_in * delta_ij(i, 1), &
                            4.0_dp * rho_in * delta_ij_dx(i, 1, j), &
                            0.0_dp, &
                            0.0_dp/)
                    x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                    x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                    call compute_qudratic_sub(x_a, x_b, x_c, xx_dx(1, j))

                    xx_dx(2, j) = xx_dx(1, j)
                    xx_dx(3, j) = xx_dx(1, j)
                    xx_dx(4, j) = xx_dx(1, j)
                end do
            end if

        elseif (tmp_bond == 42) then
            n_site(i) = 4

            ! these are not same as VV coefitients
            x_a = (/3.0_dp * rho_in * delta_ij(i, 1), &
                    3.0_dp * rho_in * delta_ij_drho(i, 1) + 3.0_dp * delta_ij(i, 1), &
                    3.0_dp * rho_in * delta_ij_drho2(i, 1) + 6.0_dp * delta_ij_drho(i, 1), &
                    3.0_dp * rho_in * delta_ij_drho3(i, 1) + 9.0_dp * delta_ij_drho2(i, 1)/)
            x_b = (/1.0_dp - 2.0_dp * rho_in * delta_ij(i, 1), &
                    -2.0_dp * rho_in * delta_ij_drho(i, 1) - 2.0_dp * delta_ij(i, 1), &
                    -2.0_dp * rho_in * delta_ij_drho2(i, 1) - 4.0_dp * delta_ij_drho(i, 1), &
                    -2.0_dp * rho_in * delta_ij_drho3(i, 1) - 6.0_dp * delta_ij_drho2(i, 1)/)
            x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

            call compute_qudratic3_sub(x_a, x_b, x_c, xx(1), xx_drho(1), xx_drho2(1), xx_drho3(1)) ! compute xx, xx_drho...
            ! hold the approx equations between X^A X^B X^C X^D
            xx(2) = xx(1)
            xx_drho(2) = xx_drho(1)
            xx_drho2(2) = xx_drho2(1)
            xx_drho3(2) = xx_drho3(1)

            xx(3) = xx(1)
            xx_drho(3) = xx_drho(1)
            xx_drho2(3) = xx_drho2(1)
            xx_drho3(3) = xx_drho3(1)

            xx(4) = 3.0_dp * xx(1) - 2.0_dp
            xx_drho(4) = 3.0_dp * xx_drho(1)
            xx_drho2(4) = 3.0_dp * xx_drho2(1)
            xx_drho3(4) = 3.0_dp * xx_drho3(1)

            if (temperature_der .eqv. .true.) then
                x_a = (/3.0_dp * rho_in * delta_ij(i, 1), &
                        3.0_dp * rho_in * delta_ij_dtt(i, 1), &
                        0.0_dp, &
                        0.0_dp/)
                x_b = (/1.0_dp - 2.0_dp * rho_in * delta_ij(i, 1), &
                        -2.0_dp * rho_in * delta_ij_dtt(i, 1), &
                        0.0_dp, &
                        0.0_dp/)
                x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                call compute_qudratic_sub(x_a, x_b, x_c, xx_dtt(1))

                xx_dtt(2) = xx_dtt(1)
                xx_dtt(3) = xx_dtt(1)
                xx_dtt(4) = 3.0_dp * xx_dtt(1)
            end if

            if (composition_der .eqv. .true.) then
                do j = 1, n_comp
                    x_a = (/3.0_dp * rho_in * delta_ij(i, 1), &
                            3.0_dp * rho_in * delta_ij_dx(i, 1, j), &
                            0.0_dp, &
                            0.0_dp/)
                    x_b = (/1.0_dp - 2.0_dp * rho_in * delta_ij(i, 1), &
                            -2.0_dp * rho_in * delta_ij_dx(i, 1, j), &
                            0.0_dp, &
                            0.0_dp/)
                    x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                    call compute_qudratic_sub(x_a, x_b, x_c, xx_dx(1, j))

                    xx_dx(2, j) = xx_dx(1, j)
                    xx_dx(3, j) = xx_dx(1, j)
                    xx_dx(4, j) = 3.0_dp * xx_dx(1, j)
                end do
            end if

        elseif (tmp_bond == 43) then
            n_site(i) = 4

            x_a = (/2.0_dp * rho_in * delta_ij(i, 1), &
                    2.0_dp * rho_in * delta_ij_drho(i, 1) + 2.0_dp * delta_ij(i, 1), &
                    2.0_dp * rho_in * delta_ij_drho2(i, 1) + 4.0_dp * delta_ij_drho(i, 1), &
                    2.0_dp * rho_in * delta_ij_drho3(i, 1) + 6.0_dp * delta_ij_drho2(i, 1)/)
            x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
            x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

            call compute_qudratic3_sub(x_a, x_b, x_c, xx(1), xx_drho(1), xx_drho2(1), xx_drho3(1)) ! compute xx, xx_drho...
            ! hold the approx equations between X^A X^B X^C X^D
            xx(2) = xx(1)
            xx_drho(2) = xx_drho(1)
            xx_drho2(2) = xx_drho2(1)
            xx_drho3(2) = xx_drho3(1)

            xx(3) = xx(1)
            xx_drho(3) = xx_drho(1)
            xx_drho2(3) = xx_drho2(1)
            xx_drho3(3) = xx_drho3(1)

            xx(4) = xx(1)
            xx_drho(4) = xx_drho(1)
            xx_drho2(4) = xx_drho2(1)
            xx_drho3(4) = xx_drho3(1)

            if (temperature_der .eqv. .true.) then
                x_a = (/2.0_dp * rho_in * delta_ij(i, 1), &
                        2.0_dp * rho_in * delta_ij_dtt(i, 1), &
                        0.0_dp, &
                        0.0_dp/)
                x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                call compute_qudratic_sub(x_a, x_b, x_c, xx_dtt(1))

                xx_dtt(2) = xx_dtt(1)
                xx_dtt(3) = xx_dtt(1)
                xx_dtt(4) = xx_dtt(1)
            end if

            if (composition_der .eqv. .true.) then
                do j = 1, n_comp
                    x_a = (/2.0_dp * rho_in * delta_ij(i, 1), &
                            2.0_dp * rho_in * delta_ij_dx(i, 1, j), &
                            0.0_dp, &
                            0.0_dp/)
                    x_b = (/1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)
                    x_c = (/-1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/)

                    call compute_qudratic_sub(x_a, x_b, x_c, xx_dx(1, j))

                    xx_dx(2, j) = xx_dx(1, j)
                    xx_dx(3, j) = xx_dx(1, j)
                    xx_dx(4, j) = xx_dx(1, j)
                end do
            end if

        else
            write (stderr, *) 'Error in z_asoc_sub function: unimplemented bond type'
            stop
        end if

        call write_array_1d (x_a,     "x_a      :", stddeb) ! DEBUG control computation
        call write_array_1d (x_b,     "x_b      :", stddeb) ! DEBUG control computation
        call write_array_1d (x_c,     "x_c      :", stddeb) ! DEBUG control computation
        write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug
        call write_array_1d (xx,      "xx       :", stddeb) ! DEBUG control computation
        call write_array_1d (xx_drho, "xx_drho  :", stddeb) ! DEBUG control computation
        call write_array_1d (xx_drho2,"xx_drho2 :", stddeb) ! DEBUG control computation
        call write_array_1d (xx_drho3,"xx_drho3 :", stddeb) ! DEBUG control computation
        call write_array_1d (xx_dtt,  "xx_dtt   :", stddeb) ! DEBUG control computation
        call write_array_2d (xx_dx,   "xx_dx    :", stddeb) ! DEBUG control computation
        write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

        ! ----- asociation term calculation section -----
        ! the pure components version evaluation
        a_asoc_(i) = n_site(i) / 2.0_dp ! the +n_site(i)/2.0_dp is for M_i/2
        z_asoc_(i) = 0.0_dp
        z_asoc_drho_(i) = 0.0_dp
        z_asoc_drho2_(i) = 0.0_dp

        ! j index iterate over sites
        do j = 1, n_site(i)
            a_asoc_(i) = a_asoc_(i) + log(xx(j)) - xx(j) / 2.0_dp

            z_asoc_(i) = z_asoc_(i) + (1.0_dp / xx(j) - 1.0_dp / 2.0_dp) * xx_drho(j)
            z_asoc_drho_(i) = z_asoc_drho_(i) - (xx_drho(j) / xx(j))**2 &
                    & + (1.0_dp / xx(j) - 1.0_dp / 2.0_dp) * xx_drho2(j)
            if (sec_der .eqv. .true.) then
                z_asoc_drho2_(i) = z_asoc_drho2_(i) + 2.0_dp * (xx_drho(j) / xx(j))**3 &
                        & - 3.0_dp * xx_drho(j) * xx_drho2(j) / xx(j)**2 &
                        & + (1.0_dp / xx(j) - 1.0_dp / 2.0_dp) * xx_drho3(j)
            end if
            if (temperature_der .eqv. .true.) then
                a_asoc_dtt_(i) = a_asoc_dtt_(i) + (1.0_dp / xx(j) - 1.0_dp / 2.0_dp) * xx_dtt(j)
            end if
            if (composition_der .eqv. .true.) then
                do k = 1, n_comp ! NOTE not efficient summation
                    a_asoc_dx_(i, k) = a_asoc_dx_(i, k) + (1.0_dp / xx(j) - 1.0_dp / 2.0_dp) * xx_dx(j, k)
                end do
            end if
        end do

        write(stddeb, "('a_asoc_(',i2,')      : ',e16.8)") i, a_asoc_(i) ! DEBUG
        write(stddeb, "('a_asoc_dtt_(',i2,')  : ',e16.8)") i, a_asoc_dtt_(i) ! DEBUG
        call write_array_1d(a_asoc_dx_(i,:), 'a_asoc_dx_(i)    : ') ! DEBUG

        write(stddeb, "('z_asoc_(',i2,')      : ',e16.8)") i, z_asoc_(i) ! DEBUG
        write(stddeb, "('z_asoc_drho_(',i2,') : ',e16.8)") i, z_asoc_drho_(i) ! DEBUG
        write(stddeb, "('z_asoc_drho2_(',i2,'): ',e16.8)") i, z_asoc_drho2_(i) ! DEBUG
    end do ! end of i cycle from line ~127

    call write_array_2d(a_asoc_dx_, '*** a_asoc_dx_ ***') ! DEBUG

    ! the finalization for mixtures
    a_asoc = 0.0_dp
    a_asoc_dtt = 0.0_dp
    a_asoc_dx = 0.0_dp
    z_asoc = 0.0_dp
    z_asoc_drho = 0.0_dp
    z_asoc_drho2 = 0.0_dp
    do i = 1, n_comp
        ! 		a_asoc = a_asoc + mol_rat(i)*a_asoc_(i)
        ! 		z_asoc = z_asoc + mol_rat(i)*z_asoc_(i)*rho_in
        ! 		z_asoc_drho = z_asoc_drho + mol_rat(i)*(z_asoc_drho_(i)*rho_in + z_asoc_(i))

        ! BUG this is only for single component (as it does not contain mol_rat)
        a_asoc = a_asoc + a_asoc_(i)
        z_asoc = z_asoc + z_asoc_(i) * rho_in
        z_asoc_drho = z_asoc_drho + (z_asoc_drho_(i) * rho_in + z_asoc_(i))

        if (sec_der .eqv. .true.) then
            z_asoc_drho2 = z_asoc_drho2 + mol_rat(i) * (z_asoc_drho2_(i) * rho_in + 2.0_dp * z_asoc_drho_(i))
        end if
        if (temperature_der .eqv. .true.) then
            a_asoc_dtt = a_asoc_dtt + mol_rat(i) * a_asoc_dtt_(i)
        end if
        if (composition_der .eqv. .true.) then
            a_asoc_dx(i) = a_asoc_(i)
            do j = 1, n_comp
                a_asoc_dx(i) = a_asoc_dx(i) + mol_rat(j) * a_asoc_dx_(j, i)
            end do
        end if
    end do

    ! 	a_asoc_dtt =
    ! 	a_asoc_dx =
    write (stddeb,*) "z_asoc_      : ", z_asoc_ ! DEBUG just for brevity of debug
    write (stddeb,*) "z_asoc_drho_ : ", z_asoc_drho_ ! DEBUG just for brevity of debug
    if (sec_der .eqv. .true.) then
        write (stddeb,*) "z_asoc_drho2_: ", z_asoc_drho2_ ! DEBUG just for brevity of debug
    end if

    write (stddeb,*) "a_asoc_      : ", a_asoc_ ! DEBUG just for brevity of debug
    if (temperature_der .eqv. .true.) then
        write (stddeb,*) "a_asoc_dtt   : ", a_asoc_dtt ! DEBUG just for brevity of debug
    end if
    if (composition_der .eqv. .true.) then
        write (stddeb,*) "a_asoc_dx    : ", a_asoc_dx ! DEBUG just for brevity of debug
    end if
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    init_flag_aasoc = .true.
    init_flag_zasoc = .true.

    ! ----- deallocation section -----
    deallocate(dzeta_dx_r)
    deallocate(a_asoc_)
    deallocate(a_asoc_dtt_)
    deallocate(a_asoc_dx_)
    deallocate(z_asoc_)
    deallocate(z_asoc_drho_)
    deallocate(z_asoc_drho2_)
    deallocate(delta_ij)
    deallocate(delta_ij_drho)
    deallocate(delta_ij_drho2)
    deallocate(delta_ij_drho3)
    deallocate(delta_ij_dtt)
    deallocate(delta_ij_dx)
    deallocate(xx_dx)
    deallocate(n_site)

    write (stddeb,*) "<== z_assoc_sub" !DEBUG
end subroutine z_asoc_sub

