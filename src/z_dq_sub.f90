subroutine z_dq_sub (rho_in)
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
    !	11.7.2017 - D.C - finished implementation also with all derivatives

    use control_mod, only : dp, n_comp, pi, tt, mol_rat, &
            &    stdout, stderr, stddeb, stdlog, &
            &    write_array_1d, write_array_2d, write_array_3d, &
            &    sec_der, temperature_der, composition_der, &
            &    init_flag, init_flag_adq, init_flag_zdq, &
            &    a_dq, a_dq_dtt, z_dq, z_dq_drho, z_dq_drho2, a_dq_dx
    use contrib_mod, only : a__dq, b__dq, c__dq, &
            &    dzeta, dzeta_dtt, dzeta_dx, &
            &    eps_ij, ees3s5_s5nq, eees4s4s4_s2s2s2nnqnqq, &
            &    initialize_dzeta_fun, initialize_polrules_fun
    use param_list_mod, only : param_m
    implicit none
    ! variable section IN/OUT
    real(dp), intent(in) :: rho_in ! the input density
    ! local variables outputed into control_mod
    real(dp) :: a_dq_, z_dq_, z_dq_drho_, z_dq_drho2_ ! the dispersive helmholtz contribution
    ! local variables
    real(dp) :: a2, a2_drho, a2_drho2, a2_drho3, a2_dtt ! the helper values for a_dq_ and derivatives computation
    real(dp) :: a3, a3_drho, a3_drho2, a3_drho3, a3_dtt ! the helper values for a_dq_ and derivatives computation
    real(dp), dimension(:), allocatable :: a2_dx, a3_dx ! helper values for a_dq_dx derivative, size [n_comp]
    real(dp), dimension(0:3) :: dzeta_r ! the dzeta vector dzeta_n for n=0,1,2,3 * rho !!
    real(dp), dimension(:, :), allocatable :: dzeta_dx_r ! the dzeta composition derivative vector dzeta_dx_r for n=0,1,2,3 * rho !!, size [0:3,1:n_comp]
    real(dp), dimension(0:3) :: dzeta_dtt_r ! the dzeta temperature derivative vector dzeta_dtt_r for n=0,1,2,3 * rho !!
    real(dp), dimension(:, :), allocatable :: jj2_dq, jj2_dq_drho, jj2_dq_drho2, jj2_dq_drho3, jj2_dq_dtt  ! the intermediate J2 terms and corresponding derivatives, size [n_comp,n_comp]
    real(dp), dimension(:, :, :), allocatable :: jj2_dq_dx ! the intermediate J2 composition derivative, size [n_comp,n_comp,n_comp]
    real(dp), dimension(:, :, :), allocatable :: jj3_dq, jj3_dq_drho, jj3_dq_drho2, jj3_dq_drho3, jj3_dq_dtt ! the intermediate J3 terms and corresponding derivatives, size [n_comp,n_comp,n_comp]
    real(dp), dimension(:, :, :, :), allocatable :: jj3_dq_dx ! the intermediate J2 composition derivative, size [n_comp,n_comp,n_comp,n_comp]
    real(dp), dimension(:, :), allocatable :: mol_rat2 ! the 2 molar ratios combination (x_i,x_j), size [n_comp,n_comp]
    real(dp), dimension(:, :, :), allocatable :: mol_rat2_dx ! the 2 molar ratios derivatives d/dx_i i,j=1,n_comp(x_i,x_j), size [n_comp,n_comp,n_comp]
    real(dp), dimension(:, :, :), allocatable :: mol_rat3 ! the 3 molar ratios combination (x_i,x_j,x_k), size [n_comp,n_comp,n_comp]
    real(dp), dimension(:, :, :, :), allocatable :: mol_rat3_dx ! the 3 molar ratios derivatives d/dx_i i,j,k=1,n_comp(x_i,x_j,x_k), size [n_comp,n_comp,n_comp,n_comp]
    real(dp), dimension(:), allocatable :: ts_2clj ! helper variable for conversion between TS and 2CLJ formulation = 4/m^2, size [n_comp]
    integer :: i, j, k, nn, l !,j_,k_! counters

    ! 	real(dp), dimension(n_comp, n_comp, n_comp) :: tmp_arr !TODO delete afterwards
    !------- PREPARING ------
    write (stddeb,*) "==> z_dq_sub" !DEBUG
    ! call for the required helper procedures
    if (init_flag(2) .eqv. .false.) then
        if (initialize_dzeta_fun() .eqv. .false.) then
            write (stderr, *) 'Error in helper function initialize_dzeta_fun'
            stop
        end if
    end if

    if (init_flag(6) .eqv. .false.) then
        if (initialize_polrules_fun() .eqv. .false.) then
            write (stderr, *) 'Error in helper function initialize_polrules_fun'
            stop
        end if
    end if

    ! ----- allocation section -----
    allocate(a2_dx(n_comp), a3_dx(n_comp))

    allocate(dzeta_dx_r(0:3, 1:n_comp))

    allocate(jj2_dq(n_comp, n_comp), jj2_dq_drho(n_comp, n_comp), jj2_dq_drho2(n_comp, n_comp), jj2_dq_drho3(n_comp, n_comp), jj2_dq_dtt(n_comp, n_comp))
    allocate(jj2_dq_dx(n_comp, n_comp, n_comp))

    allocate(jj3_dq(n_comp, n_comp, n_comp), jj3_dq_drho(n_comp, n_comp, n_comp), jj3_dq_drho2(n_comp, n_comp, n_comp), jj3_dq_drho3(n_comp, n_comp, n_comp), jj3_dq_dtt(n_comp, n_comp, n_comp))
    allocate(jj3_dq_dx(n_comp, n_comp, n_comp, n_comp))

    allocate(mol_rat2(n_comp, n_comp))
    allocate(mol_rat2_dx(n_comp, n_comp, n_comp))

    allocate(mol_rat3(n_comp, n_comp, n_comp))
    allocate(mol_rat3_dx(n_comp, n_comp, n_comp, n_comp))

    allocate(ts_2clj(n_comp))
    ! the actual dzeta (multiplied by density)
    dzeta_r = dzeta * rho_in
    dzeta_dx_r = dzeta_dx * rho_in
    dzeta_dtt_r = dzeta_dtt * rho_in

    ! DEBUG control computation
    write (stddeb,*) '-----  DQ  -----'  ! DEBUG just for brevity of debug
    call write_array_1d (dzeta,"dzeta      : ",stddeb)
    call write_array_1d (dzeta_r,"dzeta_r    : ",stddeb)

    call write_array_1d (dzeta_dtt,"dzeta_dtt  : ",stddeb)
    call write_array_1d (dzeta_dtt_r,"dzeta_dtt_r: ",stddeb)

    call write_array_2d (dzeta_dx,"dzeta_dx   : ",stddeb)
    call write_array_2d (dzeta_dx_r,"dzeta_dx_r : ",stddeb)
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug
    call write_array_2d (a__dq(1,:,:),"a__dq1  : ",stddeb)
    call write_array_2d (a__dq(2,:,:),"a__dq2  : ",stddeb)
    call write_array_2d (a__dq(3,:,:),"a__dq3  : ",stddeb)
    call write_array_2d (a__dq(4,:,:),"a__dq4  : ",stddeb)
    call write_array_2d (a__dq(5,:,:),"a__dq5  : ",stddeb)

    call write_array_2d (b__dq(1,:,:),"b__dq1  : ",stddeb)
    call write_array_2d (b__dq(2,:,:),"b__dq2  : ",stddeb)
    call write_array_2d (b__dq(3,:,:),"b__dq3  : ",stddeb)
    call write_array_2d (b__dq(4,:,:),"b__dq4  : ",stddeb)
    call write_array_2d (b__dq(5,:,:),"b__dq5  : ",stddeb)

    call write_array_3d (c__dq(1,:,:,:),"c__dq1 : ",stddeb)
    call write_array_3d (c__dq(2,:,:,:),"c__dq2 : ",stddeb)
    call write_array_3d (c__dq(3,:,:,:),"c__dq3 : ",stddeb)
    call write_array_3d (c__dq(4,:,:,:),"c__dq4 : ",stddeb)
    call write_array_3d (c__dq(5,:,:,:),"c__dq5 : ",stddeb)
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug
    ! DEBUG end of control computation

    ! initialization
    jj2_dq = 0.0_dp
    jj2_dq_drho = 0.0_dp
    jj2_dq_drho2 = 0.0_dp
    jj2_dq_drho3 = 0.0_dp
    jj2_dq_dtt = 0.0_dp
    jj2_dq_dx = 0.0_dp

    jj3_dq = 0.0_dp
    jj3_dq_drho = 0.0_dp
    jj3_dq_drho2 = 0.0_dp
    jj3_dq_drho3 = 0.0_dp
    jj3_dq_dtt = 0.0_dp
    jj3_dq_dx = 0.0_dp

    mol_rat2_dx = 0.0_dp
    mol_rat3_dx = 0.0_dp

    ! based upon paper Vrabec-Gross 2007 application of table 1
    ts_2clj(:) = 4.0_dp / param_m(:)**2

    ! J2_DQ and J3_DQ variables calculation - array multiplication is deprechated because of the initialization overhead and additional inner loops required for multiplication
    ! based upon paper Gross-Vrabec 2007 application of eq (16, 17)
    ! use of the symmetry of a__dq(i,j)=a__dq(j,i) & b__dq(i,j)=b__dq(j,i) & c__dq(i,j,k)=c__dq(i,k,j)=...=c__dq(k,j,i)
    ! computation is performed only for the nonsymetrycal elements and the symmetry is applied after the summation, which ensures that symmetrical elemets are coppied only once (but imply the use of next nested do loops)

    do nn = 1, 5
        do i = 1, n_comp ! --- diagonal terms i=j, i=j=k---
            ! exponent is preferred as integer because of subtraction
            jj2_dq(i, i) = jj2_dq(i, i) + (a__dq(nn, i, i) + b__dq(nn, i, i) * eps_ij(i, i) * ts_2clj(i) / tt) * dzeta_r(3)**(nn - 1)
            jj2_dq_drho(i, i) = jj2_dq_drho(i, i) + (a__dq(nn, i, i) + b__dq(nn, i, i) * eps_ij(i, i) * ts_2clj(i) / tt) * (nn - 1.0_dp) * dzeta_r(3)**(nn - 2)
            jj2_dq_drho2(i, i) = jj2_dq_drho2(i, i) + (a__dq(nn, i, i) + b__dq(nn, i, i) * eps_ij(i, i) * ts_2clj(i) / tt) * (nn - 1.0_dp) * (nn - 2.0_dp) * dzeta_r(3)**(nn - 3)

            jj3_dq(i, i, i) = jj3_dq(i, i, i) + c__dq(nn, i, i, i) * dzeta_r(3)**(nn - 1)
            jj3_dq_drho(i, i, i) = jj3_dq_drho(i, i, i) + c__dq(nn, i, i, i) * (nn - 1.0_dp) * dzeta_r(3)**(nn - 2)
            jj3_dq_drho2(i, i, i) = jj3_dq_drho2(i, i, i) + c__dq(nn, i, i, i) * (nn - 1.0_dp) * (nn - 2.0_dp) * dzeta_r(3)**(nn - 3)

            ! the second density derivative
            if (sec_der .eqv. .true.) then
                jj2_dq_drho3(i, i) = jj2_dq_drho3(i, i) + (a__dq(nn, i, i) + b__dq(nn, i, i) * eps_ij(i, i) / tt) * (nn - 1.0_dp) * (nn - 2.0_dp) * (nn - 3.0_dp) * dzeta_r(3)**(nn - 4)
                jj3_dq_drho3(i, i, i) = jj3_dq_drho3(i, i, i) + c__dq(nn, i, i, i) * (nn - 1.0_dp) * (nn - 2.0_dp) * (nn - 3.0_dp) * dzeta_r(3)**(nn - 4)
            end if
            ! the temperature derivative
            if (temperature_der .eqv. .true.) then
                jj2_dq_dtt(i, i) = jj2_dq_dtt(i, i) &
                        & + (b__dq(nn, i, i) * eps_ij(i, i) / tt &
                                + a__dq(nn, i, i)) * dzeta_r(3)**(nn - 2) * (nn - 1.0_dp) * dzeta_dtt_r(3) &
                        & - b__dq(nn, i, i) * eps_ij(i, i) / tt**2 * dzeta_r(3)**(nn - 1)
                jj3_dq_dtt(i, i, i) = jj3_dq_dtt(i, i, i) &
                        & + c__dq(nn, i, i, i) * dzeta_r(3)**(nn - 2) * (nn - 1.0_dp) * dzeta_dtt_r(3)
            end if
            do j = i + 1, n_comp ! --- planar terms i!=j, i!=j=k---

                jj2_dq(i, j) = jj2_dq(i, j) + (a__dq(nn, i, j) + b__dq(nn, i, j) * eps_ij(i, j) * ts_2clj(i) / tt) * dzeta_r(3)**(nn - 1)
                jj2_dq_drho(i, j) = jj2_dq_drho(i, j) + (a__dq(nn, i, j) + b__dq(nn, i, j) * eps_ij(i, j) * ts_2clj(i) / tt) * (nn - 1.0_dp) * dzeta_r(3)**(nn - 2)
                jj2_dq_drho2(i, j) = jj2_dq_drho2(i, j) + (a__dq(nn, i, j) + b__dq(nn, i, j) * eps_ij(i, j) * ts_2clj(i) / tt) * (nn - 1.0_dp) * (nn - 2.0_dp) * dzeta_r(3)**(nn - 3)

                jj3_dq(i, j, j) = jj3_dq(i, j, j) + c__dq(nn, i, j, j) * dzeta_r(3)**(nn - 1)
                jj3_dq_drho(i, j, j) = jj3_dq_drho(i, j, j) + c__dq(nn, i, j, j) * (nn - 1.0_dp) * dzeta_r(3)**(nn - 2)
                jj3_dq_drho2(i, j, j) = jj3_dq_drho2(i, j, j) + c__dq(nn, i, j, j) * (nn - 1.0_dp) * (nn - 2.0_dp) * dzeta_r(3)**(nn - 3)

                jj3_dq(i, i, j) = jj3_dq(i, i, j) + c__dq(nn, i, i, j) * dzeta_r(3)**(nn - 1)
                jj3_dq_drho(i, i, j) = jj3_dq_drho(i, i, j) + c__dq(nn, i, i, j) * (nn - 1.0_dp) * dzeta_r(3)**(nn - 2)
                jj3_dq_drho2(i, i, j) = jj3_dq_drho2(i, i, j) + c__dq(nn, i, i, j) * (nn - 1.0_dp) * (nn - 2.0_dp) * dzeta_r(3)**(nn - 3)

                ! the second density derivative
                if (sec_der .eqv. .true.) then
                    jj2_dq_drho3(i, j) = jj2_dq_drho3(i, j) &
                            & + (a__dq(nn, i, j) + b__dq(nn, i, j) * eps_ij(i, j) * ts_2clj(i) / tt) * (nn - 1.0_dp) * (nn - 2.0_dp) * (nn - 3.0_dp) * dzeta_r(3)**(nn - 4)

                    jj3_dq_drho3(i, j, j) = jj3_dq_drho3(i, j, j) &
                            & + c__dq(nn, i, j, j) * (nn - 1.0_dp) * (nn - 2.0_dp) * (nn - 3.0_dp) * dzeta_r(3)**(nn - 4)
                    jj3_dq_drho3(i, i, j) = jj3_dq_drho3(i, i, j) &
                            & + c__dq(nn, i, i, j) * (nn - 1.0_dp) * (nn - 2.0_dp) * (nn - 3.0_dp) * dzeta_r(3)**(nn - 4)
                end if

                ! the temperature derivative
                if (temperature_der .eqv. .true.) then
                    jj2_dq_dtt(i, j) = jj2_dq_dtt(i, j) &
                            & + (b__dq(nn, i, j) * eps_ij(i, j) * ts_2clj(i) / tt &
                                    + a__dq(nn, i, j)) * dzeta_r(3)**(nn - 2) * (nn - 1.0_dp) * dzeta_dtt_r(3) &
                            & - b__dq(nn, i, j) * eps_ij(i, j) * ts_2clj(i) / tt**2 * dzeta_r(3)**(nn - 1)
                    jj3_dq_dtt(i, j, j) = jj3_dq_dtt(i, j, j) &
                            & + c__dq(nn, i, j, j) * dzeta_r(3)**(nn - 2) * (nn - 1.0_dp) * dzeta_dtt_r(3)
                    jj3_dq_dtt(i, i, j) = jj3_dq_dtt(i, i, j) &
                            & + c__dq(nn, i, i, j) * dzeta_r(3)**(nn - 2) * (nn - 1.0_dp) * dzeta_dtt_r(3)
                end if
                do k = j + 1, n_comp
                    jj3_dq(i, j, k) = jj3_dq(i, j, k) + c__dq(nn, i, j, k) * dzeta_r(3)**(nn - 1)
                    jj3_dq_drho(i, j, k) = jj3_dq_drho(i, j, k) + c__dq(nn, i, j, k) * (nn - 1.0_dp) * dzeta_r(3)**(nn - 2)
                    jj3_dq_drho2(i, j, k) = jj3_dq_drho2(i, j, k) + c__dq(nn, i, j, k) * (nn - 1.0_dp) * (nn - 2.0_dp) * dzeta_r(3)**(nn - 3)
                    if (sec_der .eqv. .true.) then
                        jj3_dq_drho3(i, j, k) = jj3_dq_drho3(i, j, k) &
                                & + c__dq(nn, i, j, k) * (nn - 1.0_dp) * (nn - 2.0_dp) * (nn - 3.0_dp) * dzeta_r(3)**(nn - 4)
                    end if
                    if (temperature_der .eqv. .true.) then
                        jj3_dq_dtt(i, j, k) = jj3_dq_dtt(i, j, k) &
                                & + c__dq(nn, i, j, k) * dzeta_r(3)**(nn - 2) * (nn - 1.0_dp) * dzeta_dtt_r(3)
                    end if
                end do
            end do
        end do
    end do ! the sumation nn

    call write_array_2d (jj2_dq,"jj2_dq             :", stddeb) ! DEBUG control computation
    call write_array_3d (jj2_dq_dx,"jj2_dq_dx          :", stddeb) ! DEBUG control computation

    call write_array_3d (jj3_dq,"jj3_dq             :", stddeb) ! DEBUG control computation
    call write_array_3d (jj3_dq_drho,"jj3_dq_drho        :", stddeb) ! DEBUG control computation
    call write_array_3d (jj3_dq_dx(:,:,:,1),"jj3_dq_dx(:,:,:,1) :", stddeb) ! DEBUG control computation
    call write_array_3d (jj3_dq_dx(:,:,:,2),"jj3_dq_dx(:,:,:,2) :", stddeb) ! DEBUG control computation
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    !! apply the symmetry and multiplication by dzeta(3)**n for n-th derivative
    do i = 1, n_comp
        jj2_dq_drho(i, i) = jj2_dq_drho(i, i) * dzeta(3) ! multiplication problem
        jj2_dq_drho2(i, i) = jj2_dq_drho2(i, i) * dzeta(3)**2 ! multiplication problem

        jj3_dq_drho(i, i, i) = jj3_dq_drho(i, i, i) * dzeta(3) ! multiplication problem
        jj3_dq_drho2(i, i, i) = jj3_dq_drho2(i, i, i) * dzeta(3)**2 ! multiplication problem
        if (sec_der .eqv. .true.) then
            jj2_dq_drho3(i, i) = jj2_dq_drho3(i, i) * dzeta(3)**3 ! multiplication problem
            jj3_dq_drho3(i, i, i) = jj3_dq_drho3(i, i, i) * dzeta(3)**3 ! multiplication problem
        end if
        do j = i + 1, n_comp
            ! jj2 section
            jj2_dq(j, i) = jj2_dq(i, j)
            jj2_dq_drho(i, j) = jj2_dq_drho(i, j) * dzeta(3) ! multiplication problem
            jj2_dq_drho(j, i) = jj2_dq_drho(i, j)
            jj2_dq_drho2(i, j) = jj2_dq_drho2(i, j) * dzeta(3)**2 ! multiplication problem
            jj2_dq_drho2(j, i) = jj2_dq_drho2(i, j)
            ! jj3 section with iij terms
            jj3_dq(i, j, i) = jj3_dq(i, i, j)
            jj3_dq(j, i, i) = jj3_dq(i, i, j)
            jj3_dq_drho(i, i, j) = jj3_dq_drho(i, i, j) * dzeta(3) ! multiplication problem
            jj3_dq_drho(i, j, i) = jj3_dq_drho(i, i, j)
            jj3_dq_drho(j, i, i) = jj3_dq_drho(i, i, j)
            jj3_dq_drho2(i, i, j) = jj3_dq_drho2(i, i, j) * dzeta(3)**2 ! multiplication problem
            jj3_dq_drho2(i, j, i) = jj3_dq_drho2(i, i, j)
            jj3_dq_drho2(j, i, i) = jj3_dq_drho2(i, i, j)
            ! jj3 section with ijj terms
            jj3_dq(j, i, j) = jj3_dq(i, j, j)
            jj3_dq(j, j, i) = jj3_dq(i, j, j)
            jj3_dq_drho(i, j, j) = jj3_dq_drho(i, j, j) * dzeta(3) ! multiplication problem
            jj3_dq_drho(j, i, j) = jj3_dq_drho(i, j, j)
            jj3_dq_drho(j, j, i) = jj3_dq_drho(i, j, j)
            jj3_dq_drho2(i, j, j) = jj3_dq_drho2(i, j, j) * dzeta(3)**2 ! multiplication problem
            jj3_dq_drho2(j, i, j) = jj3_dq_drho2(i, j, j)
            jj3_dq_drho2(j, j, i) = jj3_dq_drho2(i, j, j)

            ! the second density derivative
            if (sec_der .eqv. .true.) then
                jj2_dq_drho3(i, j) = jj2_dq_drho3(i, j) * dzeta(3)**3 ! multiplication problem
                jj2_dq_drho3(j, i) = jj2_dq_drho3(i, j)

                jj3_dq_drho3(i, i, j) = jj3_dq_drho3(i, i, j) * dzeta(3)**3 ! multiplication problem
                jj3_dq_drho3(i, j, i) = jj3_dq_drho3(i, i, j)
                jj3_dq_drho3(j, i, i) = jj3_dq_drho3(i, i, j)

                jj3_dq_drho3(i, j, j) = jj3_dq_drho3(i, j, j) * dzeta(3)**3 ! multiplication problem
                jj3_dq_drho3(j, i, j) = jj3_dq_drho3(i, j, j)
                jj3_dq_drho3(j, j, i) = jj3_dq_drho3(i, j, j)
            end if
            ! the temperature derivative
            if (temperature_der .eqv. .true.) then
                jj2_dq_dtt(j, i) = jj2_dq_dtt(i, j)

                jj3_dq_dtt(i, j, i) = jj3_dq_dtt(i, i, j)
                jj3_dq_dtt(j, i, i) = jj3_dq_dtt(i, i, j)

                jj3_dq_dtt(j, i, j) = jj3_dq_dtt(i, j, j)
                jj3_dq_dtt(j, j, i) = jj3_dq_dtt(i, j, j)
            end if
            do k = j + 1, n_comp
                jj3_dq(i, k, j) = jj3_dq(i, j, k)
                jj3_dq(k, i, j) = jj3_dq(i, j, k)
                jj3_dq(k, j, i) = jj3_dq(i, j, k)
                jj3_dq(j, i, k) = jj3_dq(i, j, k)
                jj3_dq(j, k, i) = jj3_dq(i, j, k)
                jj3_dq_drho(i, j, k) = jj3_dq_drho(i, j, k) * dzeta(3) ! multiplication problem
                jj3_dq_drho(i, k, j) = jj3_dq_drho(i, j, k)
                jj3_dq_drho(k, i, j) = jj3_dq_drho(i, j, k)
                jj3_dq_drho(k, j, i) = jj3_dq_drho(i, j, k)
                jj3_dq_drho(j, i, k) = jj3_dq_drho(i, j, k)
                jj3_dq_drho(j, k, i) = jj3_dq_drho(i, j, k)
                jj3_dq_drho2(i, j, k) = jj3_dq_drho2(i, j, k) * dzeta(3)**2 ! multiplication problem
                jj3_dq_drho2(i, k, j) = jj3_dq_drho2(i, j, k)
                jj3_dq_drho2(k, i, j) = jj3_dq_drho2(i, j, k)
                jj3_dq_drho2(k, j, i) = jj3_dq_drho2(i, j, k)
                jj3_dq_drho2(j, i, k) = jj3_dq_drho2(i, j, k)
                jj3_dq_drho2(j, k, i) = jj3_dq_drho2(i, j, k)
                if (sec_der .eqv. .true.) then
                    jj3_dq_drho3(i, j, k) = jj3_dq_drho3(i, j, k) * dzeta(3)**3 ! multiplication problem
                    jj3_dq_drho3(i, k, j) = jj3_dq_drho3(i, j, k)
                    jj3_dq_drho3(k, i, j) = jj3_dq_drho3(i, j, k)
                    jj3_dq_drho3(k, j, i) = jj3_dq_drho3(i, j, k)
                    jj3_dq_drho3(j, i, k) = jj3_dq_drho3(i, j, k)
                    jj3_dq_drho3(j, k, i) = jj3_dq_drho3(i, j, k)
                end if
                if (temperature_der .eqv. .true.) then
                    jj3_dq_dtt(i, k, j) = jj3_dq_dtt(i, j, k)
                    jj3_dq_dtt(k, i, j) = jj3_dq_dtt(i, j, k)
                    jj3_dq_dtt(k, j, i) = jj3_dq_dtt(i, j, k)
                    jj3_dq_dtt(j, i, k) = jj3_dq_dtt(i, j, k)
                    jj3_dq_dtt(j, k, i) = jj3_dq_dtt(i, j, k)
                end if
            end do
        end do
    end do

    mol_rat2 = 0.0_dp
    mol_rat3 = 0.0_dp
    mol_rat2_dx = -1.0_dp
    mol_rat3_dx = -1.0_dp

    if (composition_der .eqv. .true.) then
        do i = 1, n_comp
            jj2_dq_dx(:, :, i) = jj2_dq_drho(:, :) / dzeta(3) * dzeta_dx_r(3, i)
            jj3_dq_dx(:, :, :, i) = jj3_dq_drho(:, :, :) / dzeta(3) * dzeta_dx_r(3, i)
            do j = 1, n_comp
                mol_rat2(i, j) = mol_rat(i) * mol_rat(j)
                do l = 1, n_comp
                    if (i==j) then ! ii terms for a2
                        if (i==l) then
                            mol_rat2_dx(i, j, l) = 2.0_dp * mol_rat(i)
                        else
                            mol_rat2_dx(i, j, l) = 0.0_dp ! obsolete
                        end if
                    else ! ij terms for a2
                        if (i==l) then
                            mol_rat2_dx(i, j, l) = mol_rat(j)
                        elseif (j==l) then
                            mol_rat2_dx(i, j, l) = mol_rat(i)
                        else
                            mol_rat2_dx(i, j, l) = 0.0_dp
                        end if
                    end if
                end do
                do k = 1, n_comp
                    mol_rat3(i, j, k) = mol_rat(i) * mol_rat(j) * mol_rat(k)
                    do l = 1, n_comp
                        if (i==j) then
                            if (j==k) then ! iii terms a3
                                if (i==l) then ! d/x_i
                                    mol_rat3_dx(i, j, k, l) = 3.0_dp * mol_rat(i)**2
                                else ! d/x_j
                                    mol_rat3_dx(i, j, k, l) = 0.0_dp ! obsolete
                                end if
                            else ! iik terms a3
                                if (i==l) then ! d/x_i
                                    mol_rat3_dx(i, j, k, l) = 2.0_dp * mol_rat(i) * mol_rat(k)
                                elseif (k==l) then ! d/x_k
                                    mol_rat3_dx(i, j, k, l) = mol_rat(i)**2
                                else ! d/x_j
                                    mol_rat3_dx(i, j, k, l) = 0.0_dp ! obsolete
                                end if
                            end if
                        else
                            if (j==k) then ! ijj terms a3
                                if (i==l) then
                                    mol_rat3_dx(i, j, k, l) = mol_rat(j)**2
                                elseif (j==l) then
                                    mol_rat3_dx(i, j, k, l) = 2.0_dp * mol_rat(i) * mol_rat(j)
                                else
                                    mol_rat3_dx(i, j, k, l) = 0.0_dp ! obsolete
                                end if
                            elseif (i==k) then ! iji terms a3
                                if (i==l) then
                                    mol_rat3_dx(i, j, k, l) = 2.0_dp * mol_rat(i) * mol_rat(j)
                                elseif (j==l) then
                                    mol_rat3_dx(i, j, k, l) = mol_rat(i)**2
                                else
                                    mol_rat3_dx(i, j, k, l) = 0.0_dp ! obsolete
                                end if
                            else ! ijk terms a3
                                if (i==l) then
                                    mol_rat3_dx(i, j, k, l) = mol_rat(j) * mol_rat(k)
                                elseif (j==l) then
                                    mol_rat3_dx(i, j, k, l) = mol_rat(i) * mol_rat(k)
                                elseif (k==l) then
                                    mol_rat3_dx(i, j, k, l) = mol_rat(i) * mol_rat(j)
                                else
                                    mol_rat3_dx(i, j, k, l) = 0.0_dp
                                end if
                            end if
                        end if
                    end do
                end do
            end do
        end do
    end if

    call write_array_2d(mol_rat2,"mol_rat2             :", stddeb) !DEBUG
    call write_array_3d(mol_rat3,"mol_rat3             :", stddeb) !DEBUG

    call write_array_3d(mol_rat2_dx,"mol_rat2_dx          :", stddeb) !DEBUG
    call write_array_3d(mol_rat3_dx(:,:,:,1),"mol_rat3_dx(:,:,:,1) :", stddeb) !DEBUG
    call write_array_3d(mol_rat3_dx(:,:,:,2),"mol_rat3_dx(:,:,:,2) :", stddeb) !DEBUG
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    a2 = 0.0_dp
    a2_drho = 0.0_dp
    a2_drho2 = 0.0_dp
    a2_drho3 = 0.0_dp
    a2_dtt = 0.0_dp
    a2_dx = 0.0_dp

    a3 = 0.0_dp
    a3_drho = 0.0_dp
    a3_drho2 = 0.0_dp
    a3_drho3 = 0.0_dp
    a3_dtt = 0.0_dp
    a3_dx = 0.0_dp

    ! write (stddeb,"('ees3s5_s5nq            : ',e)"), ees3s5_s5nq(1,1) ! DEBUG control computation
    ! write (stddeb,"('eees4s4s4_s2s2s2nnqnqq : ',e)"), eees4s4s4_s2s2s2nnqnqq(1,1,1) ! DEBUG control of computation
    call write_array_2d (ees3s5_s5nq,"ees3s5_s5nq            :",stddeb) ! DEBUG control computation
    call write_array_3d (eees4s4s4_s2s2s2nnqnqq,"eees4s4s4_s2s2s2nnqnqq :",stddeb) ! DEBUG control computation
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    ! based upon paper Gross-Vrabec 2007 application of eq (14, 15)
    ! NOTE alternative full iteration - symmetry cannot be used because the ees3s5_s5nq and eees4s4s4_s2s2s2nnqnqq are not symmetrycal
    do i = 1, n_comp
        do j = 1, n_comp
            a2 = a2 &
                    & + mol_rat(i) * mol_rat(j) * ees3s5_s5nq(i, j) * jj2_dq(i, j)
            a2_drho = a2_drho &
                    & + (mol_rat(i) * mol_rat(j) * ees3s5_s5nq(i, j) * jj2_dq(i, j) &
                            & + rho_in * (mol_rat(i) * mol_rat(j) * ees3s5_s5nq(i, j) * jj2_dq_drho(i, j)))
            a2_drho2 = a2_drho2 &
                    & + (2.0_dp * (mol_rat(i) * mol_rat(j) * ees3s5_s5nq(i, j) * jj2_dq_drho(i, j)) &
                            & + rho_in * (mol_rat(i) * mol_rat(j) * ees3s5_s5nq(i, j) * jj2_dq_drho2(i, j)))

            if (sec_der .eqv. .true.) then
                a2_drho3 = a2_drho3 &
                        & + (3.0_dp * (mol_rat(i) * mol_rat(j) * ees3s5_s5nq(i, j) * jj2_dq_drho2(i, j)) &
                                & + rho_in * (mol_rat(i) * mol_rat(j) * ees3s5_s5nq(i, j) * jj2_dq_drho3(i, j)))
            end if
            if (temperature_der .eqv. .true.) then
                a2_dtt = a2_dtt &
                        & + (2.0_dp / tt * (mol_rat(i) * mol_rat(j) * ees3s5_s5nq(i, j) * jj2_dq(i, j)) &
                                & - (mol_rat(i) * mol_rat(j) * ees3s5_s5nq(i, j) * jj2_dq_dtt(i, j)))
            end if

            do k = 1, n_comp
                a3 = a3 &
                        & + mol_rat(i) * mol_rat(j) * mol_rat(k) * eees4s4s4_s2s2s2nnqnqq(i, j, k) * jj3_dq(i, j, k)
                a3_drho = a3_drho &
                        & + (2.0_dp * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees4s4s4_s2s2s2nnqnqq(i, j, k) * jj3_dq(i, j, k)) &
                                & + rho_in * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees4s4s4_s2s2s2nnqnqq(i, j, k) * jj3_dq_drho(i, j, k)))
                a3_drho2 = a3_drho2 &
                        & + (2.0_dp * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees4s4s4_s2s2s2nnqnqq(i, j, k) * jj3_dq(i, j, k)) &
                                & + 4.0_dp * rho_in * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees4s4s4_s2s2s2nnqnqq(i, j, k) * jj3_dq_drho(i, j, k)) &
                                & + rho_in**2 * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees4s4s4_s2s2s2nnqnqq(i, j, k) * jj3_dq_drho2(i, j, k)))
                if (sec_der .eqv. .true.) then
                    a3_drho3 = a3_drho3 &
                            & + (6.0_dp * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees4s4s4_s2s2s2nnqnqq(i, j, k) * jj3_dq_drho(i, j, k)) &
                                    & + 6.0_dp * rho_in * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees4s4s4_s2s2s2nnqnqq(i, j, k) * jj3_dq_drho2(i, j, k)) &
                                    & + rho_in**2 * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees4s4s4_s2s2s2nnqnqq(i, j, k) * jj3_dq_drho3(i, j, k)))
                end if
                if (temperature_der .eqv. .true.) then
                    a3_dtt = a3_dtt &
                            & + (3.0_dp / tt * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees4s4s4_s2s2s2nnqnqq(i, j, k) * jj3_dq(i, j, k)) &
                                    & - (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees4s4s4_s2s2s2nnqnqq(i, j, k) * jj3_dq_dtt(i, j, k)))
                end if
            end do
        end do
    end do

    if (composition_der .eqv. .true.) then
        do i = 1, n_comp
            a2_dx(i) = sum(mol_rat2_dx(:, :, i) * ees3s5_s5nq(:, :) * jj2_dq(:, :) &
                    & + mol_rat2(:, :) * ees3s5_s5nq(:, :) * jj2_dq_dx(:, :, i))
            a3_dx(i) = sum(mol_rat3_dx(:, :, :, i) * eees4s4s4_s2s2s2nnqnqq(:, :, :) * jj3_dq(:, :, :) &
                    & + mol_rat3(:, :, :) * eees4s4s4_s2s2s2nnqnqq(:, :, :) * jj3_dq_dx(:, :, :, i))
        end do
    end if
    a2 = a2 * (-pi * rho_in / tt**2) * (3.0_dp / 4.0_dp)**2
    a2_drho = a2_drho * (-pi / tt**2) * (3.0_dp / 4.0_dp)**2
    a2_drho2 = a2_drho2 * (-pi / tt**2) * (3.0_dp / 4.0_dp)**2
    a2_drho3 = a2_drho3 * (-pi / tt**2) * (3.0_dp / 4.0_dp)**2
    a2_dtt = a2_dtt * (pi * rho_in / tt**2) * (3.0_dp / 4.0_dp)**2
    a2_dx = a2_dx * (-pi * rho_in / tt**2) * (3.0_dp / 4.0_dp)**2

    a3 = a3 * (pi**2 * rho_in**2 / tt**3) * (3.0_dp / 4.0_dp)**2
    a3_drho = a3_drho * (pi**2 * rho_in / tt**3) * (3.0_dp / 4.0_dp)**2
    a3_drho2 = a3_drho2 * (pi**2 / tt**3) * (3.0_dp / 4.0_dp)**2
    a3_drho3 = a3_drho3 * (pi**2 / tt**3) * (3.0_dp / 4.0_dp)**2
    a3_dtt = a3_dtt * (-pi**2 * rho_in**2 / tt**3) * (3.0_dp / 4.0_dp)**2
    a3_dx = a3_dx * (pi**2 * rho_in**2 / tt**3) * (3.0_dp / 4.0_dp)**2

    write(stddeb,"('a2      : ',e16.8)") a2 ! DEBUG control computation
    write(stddeb,"('a2_drho : ',e16.8)") a2_drho ! DEBUG control computation
    write(stddeb,"('a2_dtt  : ',e16.8)") a2_dtt ! DEBUG control computation
    call write_array_1d(a2_dx,"a2_dx   : ", stddeb) ! DEBUG control of computation

    write(stddeb,"('a3      : ',e16.8)") a3 ! DEBUG control computation
    write(stddeb,"('a3_drho : ',e16.8)") a3_drho ! DEBUG control computation
    write(stddeb,"('a3_dtt  : ',e16.8)") a3_dtt ! DEBUG control computation
    call write_array_1d(a3_dx,"a3_dx   : ", stddeb) ! DEBUG control of computation
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    ! dipol-quadrupol interaction contribution computation
    ! based upon paper Vrabec-Gross 2007 application of eq (13)

    a_dq_ = a2**2 / (a2 - a3) ! TODO decide if keep and make the a_dq obsolete
    ! following equalities are derived by hand, by maple and controlled by wolfram alpha
    z_dq_ = 2.0_dp * a2_drho * a2 / (a2 - a3) &
            & - a2**2 * (a2_drho - a3_drho) / (a2 - a3)**2
    z_dq_drho_ = 2.0_dp * a2_drho2 * a2 / (a2 - a3) &
            & + 2.0_dp * a2_drho**2 / (a2 - a3) &
            & - 4.0_dp * a2_drho * a2 * (a2_drho - a3_drho) / (a2 - a3)**2 &
            & - a2**2 * (a2_drho2 - a3_drho2) / (a2 - a3)**2 &
            & + 2.0_dp * a2**2 * (a2_drho - a3_drho)**2 / (a2 - a3)**3
    z_dq_drho2_ = 2.0_dp * a2_drho3 * a2 / (a2 - a3) &
            & + 6.0_dp * a2_drho2 * a2_drho / (a2 - a3) &
            & - 6.0_dp * a2_drho2 * a2 * (a2_drho - a3_drho) / (a2 - a3)**2 &
            & - 6.0_dp * a2_drho**2 * (a2_drho - a3_drho) / (a2 - a3)**2 &
            & - 6.0_dp * a2_drho * a2 * (a2_drho2 - a3_drho2) / (a2 - a3)**2 &
            & + 12.0_dp * a2_drho * a2 * (a2_drho - a3_drho)**2 / (a2 - a3)**3 &
            & - a2**2 * (a2_drho3 - a3_drho3) / (a2 - a3)**2 &
            & + 6.0_dp * a2**2 * (a2_drho - a3_drho) * (a2_drho2 - a3_drho2) / (a2 - a3)**3 &
            & - 6.0_dp * a2**2 * (a2_drho - a3_drho)**3 / (a2 - a3)**4
    if (temperature_der .eqv. .true.) then
        a_dq_dtt = (a2**2 * a2_dtt &
                & + a2**2 * a3_dtt &
                & - 2.0_dp * a2 * a3 * a2_dtt) / (a2 - a3)**2
    end if
    if (composition_der .eqv. .true.) then
        a_dq_dx(:) = (a2**2 * a2_dx(:) &
                & + a2**2 * a3_dx(:) &
                & - 2.0_dp * a2 * a3 * a2_dx(:)) / (a2 - a3)**2
    end if

    write (stddeb,*) "z_dq_      : ", z_dq_ ! DEBUG just for brevity of debug
    write (stddeb,*) "z_dq_drho_ : ", z_dq_drho_ ! DEBUG just for brevity of debug
    if (sec_der .eqv. .true.) then
        write (stddeb,*) "z_dq_drho2_: ", z_dq_drho2_ ! DEBUG just for brevity of debug
    end if

    write (stddeb,*) "a_dq_      : ", a_dq_ ! DEBUG just for brevity of debug
    if (temperature_der .eqv. .true.) then
        write (stddeb,*) "a_dq_dtt   : ", a_dq_dtt ! DEBUG just for brevity of debug
    end if
    if (composition_der .eqv. .true.) then
        write (stddeb,*) "a_dq_dx    : ", a_dq_dx ! DEBUG just for brevity of debug
    end if
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    a_dq = a_dq_
    ! NOTE a_dq_dtt direct modification in previous section
    ! NOTE a_dq_dx direct modification in previous section

    z_dq = rho_in * z_dq_
    z_dq_drho = z_dq_ + rho_in * z_dq_drho_
    z_dq_drho2 = 2.0_dp * z_dq_drho_ + rho_in * z_dq_drho2_

    init_flag_adq = .true.
    init_flag_zdq = .true.

    ! ----- deallocation section -----
    deallocate(a2_dx, a3_dx)

    deallocate(dzeta_dx_r)

    deallocate(jj2_dq, jj2_dq_drho, jj2_dq_drho2, jj2_dq_drho3, jj2_dq_dtt)
    deallocate(jj2_dq_dx)

    deallocate(jj3_dq, jj3_dq_drho, jj3_dq_drho2, jj3_dq_drho3, jj3_dq_dtt)
    deallocate(jj3_dq_dx)

    deallocate(mol_rat2)
    deallocate(mol_rat2_dx)

    deallocate(mol_rat3)
    deallocate(mol_rat3_dx)

    deallocate(ts_2clj)

    write (stddeb,*) "<== z_dq_sub" !DEBUG
end subroutine z_dq_sub
