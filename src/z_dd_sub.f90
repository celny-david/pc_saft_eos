subroutine z_dd_sub (rho_in)
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
            &    init_flag, init_flag_add, init_flag_zdd, &
            &    a_dd, a_dd_dtt, z_dd, z_dd_drho, z_dd_drho2, a_dd_dx
    use contrib_mod, only : a__dd, b__dd, c__dd, &
            &    dzeta, dzeta_dtt, dzeta_dx, &
            &    eps_ij, ees3s3_s3nnmm, eees3s3s3_sssnnnmmm, &
            &    initialize_dzeta_fun, initialize_polrules_fun
    ! 	use param_list_mod , only: param_m
    implicit none
    ! variable section IN/OUT
    real(dp), intent(in) :: rho_in ! the input density
    ! local variables outputed into control_mod
    real(dp) :: a_dd_, z_dd_, z_dd_drho_, z_dd_drho2_ ! the dispersive helmholtz contribution
    ! local variables
    real(dp) :: a2, a2_drho, a2_drho2, a2_drho3, a2_dtt ! the helper values for a_dd_ and derivatives computation
    real(dp) :: a3, a3_drho, a3_drho2, a3_drho3, a3_dtt ! the helper values for a_dd_ and derivatives computation
    real(dp), dimension(:), allocatable :: a2_dx, a3_dx ! helper values for a_dd_dx derivative, size [n_comp]
    real(dp), dimension(0:3) :: dzeta_r ! the dzeta vector dzeta_n for n=0,1,2,3 * rho !!
    real(dp), dimension(:, :), allocatable :: dzeta_dx_r ! the dzeta composition derivative vector dzeta_dx_r for n=0,1,2,3 * rho !!, size [0:3,1:n_comp]
    real(dp), dimension(0:3) :: dzeta_dtt_r ! the dzeta temperature derivative vector dzeta_dtt_r for n=0,1,2,3 * rho !!
    real(dp), dimension(:, :), allocatable :: jj2_dd, jj2_dd_drho, jj2_dd_drho2, jj2_dd_drho3, jj2_dd_dtt  ! the intermediate J2 terms and corresponding derivatives, size [n_comp,n_comp]
    real(dp), dimension(:, :, :), allocatable :: jj2_dd_dx ! the intermediate J2 composition derivative, size [n_comp,n_comp,n_comp]
    real(dp), dimension(:, :, :), allocatable :: jj3_dd, jj3_dd_drho, jj3_dd_drho2, jj3_dd_drho3, jj3_dd_dtt ! the intermediate J3 terms and corresponding derivatives, size [n_comp,n_comp,n_comp]
    real(dp), dimension(:, :, :, :), allocatable :: jj3_dd_dx ! the intermediate J2 composition derivative, size [n_comp,n_comp,n_comp,n_comp]
    real(dp), dimension(:, :), allocatable :: mol_rat2 ! the 2 molar ratios combination (x_i,x_j), size [n_comp,n_comp]
    real(dp), dimension(:, :, :), allocatable :: mol_rat2_dx ! the 2 molar ratios derivatives d/dx_i i,j=1,n_comp(x_i,x_j), size [n_comp,n_comp,n_comp]
    real(dp), dimension(:, :, :), allocatable :: mol_rat3 ! the 3 molar ratios combination (x_i,x_j,x_k), size [n_comp,n_comp,n_comp]
    real(dp), dimension(:, :, :, :), allocatable :: mol_rat3_dx ! the 3 molar ratios derivatives d/dx_i i,j,k=1,n_comp(x_i,x_j,x_k), size [n_comp,n_comp,n_comp,n_comp]
    integer :: i, j, k, nn, l !,j_,k_! counters
    ! 	real(dp), dimension(n_comp, n_comp, n_comp) :: tmp_arr !TODO delete afterwards

    !------- PREPARING ------
    write (stddeb,*) "==> z_dd_sub" !DEBUG
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

    allocate(jj2_dd(n_comp, n_comp), jj2_dd_drho(n_comp, n_comp), jj2_dd_drho2(n_comp, n_comp), jj2_dd_drho3(n_comp, n_comp), jj2_dd_dtt(n_comp, n_comp))
    allocate(jj2_dd_dx(n_comp, n_comp, n_comp))

    allocate(jj3_dd(n_comp, n_comp, n_comp), jj3_dd_drho(n_comp, n_comp, n_comp), jj3_dd_drho2(n_comp, n_comp, n_comp), jj3_dd_drho3(n_comp, n_comp, n_comp), jj3_dd_dtt(n_comp, n_comp, n_comp))
    allocate(jj3_dd_dx(n_comp, n_comp, n_comp, n_comp))

    allocate(mol_rat2(n_comp, n_comp))
    allocate(mol_rat2_dx(n_comp, n_comp, n_comp))

    allocate(mol_rat3(n_comp, n_comp, n_comp))
    allocate(mol_rat3_dx(n_comp, n_comp, n_comp, n_comp))

    ! the actual dzeta (multiplied by density)
    dzeta_r = dzeta * rho_in
    dzeta_dx_r = dzeta_dx * rho_in
    dzeta_dtt_r = dzeta_dtt * rho_in

    ! DEBUG control computation
    write (stddeb,*) '-----  DD  -----'  ! DEBUG just for brevity of debug
    call write_array_1d (dzeta,"dzeta      : ",stddeb)
    call write_array_1d (dzeta_r,"dzeta_r    : ",stddeb)

    call write_array_1d (dzeta_dtt,"dzeta_dtt  : ",stddeb)
    call write_array_1d (dzeta_dtt_r,"dzeta_dtt_r: ",stddeb)

    call write_array_2d (dzeta_dx,"dzeta_dx   : ",stddeb)
    call write_array_2d (dzeta_dx_r,"dzeta_dx_r : ",stddeb)
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug
    call write_array_3d (a__dd,"a__dd  : ",stddeb)
    call write_array_3d (b__dd,"b__dd  : ",stddeb)
    call write_array_3d (c__dd(1,:,:,:),"c__dd1 : ",stddeb)
    call write_array_3d (c__dd(2,:,:,:),"c__dd2 : ",stddeb)
    call write_array_3d (c__dd(3,:,:,:),"c__dd3 : ",stddeb)
    call write_array_3d (c__dd(4,:,:,:),"c__dd4 : ",stddeb)
    call write_array_3d (c__dd(5,:,:,:),"c__dd5 : ",stddeb)
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug
    ! DEBUG end of control computation

    ! initialization
    jj2_dd = 0.0_dp
    jj2_dd_drho = 0.0_dp
    jj2_dd_drho2 = 0.0_dp
    jj2_dd_drho3 = 0.0_dp
    jj2_dd_dtt = 0.0_dp
    jj2_dd_dx = 0.0_dp

    jj3_dd = 0.0_dp
    jj3_dd_drho = 0.0_dp
    jj3_dd_drho2 = 0.0_dp
    jj3_dd_drho3 = 0.0_dp
    jj3_dd_dtt = 0.0_dp
    jj3_dd_dx = 0.0_dp

    mol_rat2_dx = 0.0_dp
    mol_rat3_dx = 0.0_dp

    ! J2_DD and J3_DD variables calculation - array multiplication is deprechated because of the initialization overhead and additional inner loops required for multiplication
    ! based upon paper Gross-Vrabec 2006 application of eq (10, 11)
    ! use of the symmetry of a__dd(i,j)=a__dd(j,i) & b__dd(i,j)=b__dd(j,i) & c__dd(i,j,k)=c__dd(i,k,j)=...=c__dd(k,j,i)
    ! computation is performed only for the nonsymetrycal elements and the symmetry is applied after the summation, which ensures that symmetrical elemets are coppied only once (but imply the use of next nested do loops)

    do nn = 1, 5
        do i = 1, n_comp ! --- diagonal terms i=j, i=j=k---
            ! exponent is preferred as integer because of subtraction
            jj2_dd(i, i) = jj2_dd(i, i) + (a__dd(nn, i, i) + b__dd(nn, i, i) * eps_ij(i, i) / tt) * dzeta_r(3)**(nn - 1)
            jj2_dd_drho(i, i) = jj2_dd_drho(i, i) + (a__dd(nn, i, i) + b__dd(nn, i, i) * eps_ij(i, i) / tt) * (nn - 1.0_dp) * dzeta_r(3)**(nn - 2)
            jj2_dd_drho2(i, i) = jj2_dd_drho2(i, i) + (a__dd(nn, i, i) + b__dd(nn, i, i) * eps_ij(i, i) / tt) * (nn - 1.0_dp) * (nn - 2.0_dp) * dzeta_r(3)**(nn - 3)

            jj3_dd(i, i, i) = jj3_dd(i, i, i) + c__dd(nn, i, i, i) * dzeta_r(3)**(nn - 1)
            jj3_dd_drho(i, i, i) = jj3_dd_drho(i, i, i) + c__dd(nn, i, i, i) * (nn - 1.0_dp) * dzeta_r(3)**(nn - 2)
            jj3_dd_drho2(i, i, i) = jj3_dd_drho2(i, i, i) + c__dd(nn, i, i, i) * (nn - 1.0_dp) * (nn - 2.0_dp) * dzeta_r(3)**(nn - 3)

            ! the second density derivative
            if (sec_der .eqv. .true.) then
                jj2_dd_drho3(i, i) = jj2_dd_drho3(i, i) + (a__dd(nn, i, i) + b__dd(nn, i, i) * eps_ij(i, i) / tt) * (nn - 1.0_dp) * (nn - 2.0_dp) * (nn - 3.0_dp) * dzeta_r(3)**(nn - 4)
                jj3_dd_drho3(i, i, i) = jj3_dd_drho3(i, i, i) + c__dd(nn, i, i, i) * (nn - 1.0_dp) * (nn - 2.0_dp) * (nn - 3.0_dp) * dzeta_r(3)**(nn - 4)
            end if
            ! the temperature derivative
            if (temperature_der .eqv. .true.) then
                jj2_dd_dtt(i, i) = jj2_dd_dtt(i, i) &
                        & + (b__dd(nn, i, i) * eps_ij(i, i) / tt &
                                + a__dd(nn, i, i)) * dzeta_r(3)**(nn - 2) * (nn - 1.0_dp) * dzeta_dtt_r(3) &
                        & - b__dd(nn, i, i) * eps_ij(i, i) / tt**2 * dzeta_r(3)**(nn - 1)
                jj3_dd_dtt(i, i, i) = jj3_dd_dtt(i, i, i) &
                        & + c__dd(nn, i, i, i) * dzeta_r(3)**(nn - 2) * (nn - 1.0_dp) * dzeta_dtt_r(3)
            end if
            do j = i + 1, n_comp ! --- planar terms i!=j, i!=j=k---

                jj2_dd(i, j) = jj2_dd(i, j) + (a__dd(nn, i, j) + b__dd(nn, i, j) * eps_ij(i, j) / tt) * dzeta_r(3)**(nn - 1)
                jj2_dd_drho(i, j) = jj2_dd_drho(i, j) + (a__dd(nn, i, j) + b__dd(nn, i, j) * eps_ij(i, j) / tt) * (nn - 1.0_dp) * dzeta_r(3)**(nn - 2)
                jj2_dd_drho2(i, j) = jj2_dd_drho2(i, j) + (a__dd(nn, i, j) + b__dd(nn, i, j) * eps_ij(i, j) / tt) * (nn - 1.0_dp) * (nn - 2.0_dp) * dzeta_r(3)**(nn - 3)

                jj3_dd(i, j, j) = jj3_dd(i, j, j) + c__dd(nn, i, j, j) * dzeta_r(3)**(nn - 1)
                jj3_dd_drho(i, j, j) = jj3_dd_drho(i, j, j) + c__dd(nn, i, j, j) * (nn - 1.0_dp) * dzeta_r(3)**(nn - 2)
                jj3_dd_drho2(i, j, j) = jj3_dd_drho2(i, j, j) + c__dd(nn, i, j, j) * (nn - 1.0_dp) * (nn - 2.0_dp) * dzeta_r(3)**(nn - 3)

                jj3_dd(i, i, j) = jj3_dd(i, i, j) + c__dd(nn, i, i, j) * dzeta_r(3)**(nn - 1)
                jj3_dd_drho(i, i, j) = jj3_dd_drho(i, i, j) + c__dd(nn, i, i, j) * (nn - 1.0_dp) * dzeta_r(3)**(nn - 2)
                jj3_dd_drho2(i, i, j) = jj3_dd_drho2(i, i, j) + c__dd(nn, i, i, j) * (nn - 1.0_dp) * (nn - 2.0_dp) * dzeta_r(3)**(nn - 3)

                ! the second density derivative
                if (sec_der .eqv. .true.) then
                    jj2_dd_drho3(i, j) = jj2_dd_drho3(i, j) &
                            & + (a__dd(nn, i, j) + b__dd(nn, i, j) * eps_ij(i, j) / tt) * (nn - 1.0_dp) * (nn - 2.0_dp) * (nn - 3.0_dp) * dzeta_r(3)**(nn - 4)

                    jj3_dd_drho3(i, j, j) = jj3_dd_drho3(i, j, j) &
                            & + c__dd(nn, i, j, j) * (nn - 1.0_dp) * (nn - 2.0_dp) * (nn - 3.0_dp) * dzeta_r(3)**(nn - 4)
                    jj3_dd_drho3(i, i, j) = jj3_dd_drho3(i, i, j) &
                            & + c__dd(nn, i, i, j) * (nn - 1.0_dp) * (nn - 2.0_dp) * (nn - 3.0_dp) * dzeta_r(3)**(nn - 4)
                end if

                ! the temperature derivative
                if (temperature_der .eqv. .true.) then
                    jj2_dd_dtt(i, j) = jj2_dd_dtt(i, j) &
                            & + (b__dd(nn, i, j) * eps_ij(i, j) / tt &
                                    + a__dd(nn, i, j)) * dzeta_r(3)**(nn - 2) * (nn - 1.0_dp) * dzeta_dtt_r(3) &
                            & - b__dd(nn, i, j) * eps_ij(i, j) / tt**2 * dzeta_r(3)**(nn - 1)
                    jj3_dd_dtt(i, j, j) = jj3_dd_dtt(i, j, j) &
                            & + c__dd(nn, i, j, j) * dzeta_r(3)**(nn - 2) * (nn - 1.0_dp) * dzeta_dtt_r(3)
                    jj3_dd_dtt(i, i, j) = jj3_dd_dtt(i, i, j) &
                            & + c__dd(nn, i, i, j) * dzeta_r(3)**(nn - 2) * (nn - 1.0_dp) * dzeta_dtt_r(3)
                end if
                do k = j + 1, n_comp
                    jj3_dd(i, j, k) = jj3_dd(i, j, k) + c__dd(nn, i, j, k) * dzeta_r(3)**(nn - 1)
                    jj3_dd_drho(i, j, k) = jj3_dd_drho(i, j, k) + c__dd(nn, i, j, k) * (nn - 1.0_dp) * dzeta_r(3)**(nn - 2)
                    jj3_dd_drho2(i, j, k) = jj3_dd_drho2(i, j, k) + c__dd(nn, i, j, k) * (nn - 1.0_dp) * (nn - 2.0_dp) * dzeta_r(3)**(nn - 3)
                    if (sec_der .eqv. .true.) then
                        jj3_dd_drho3(i, j, k) = jj3_dd_drho3(i, j, k) &
                                & + c__dd(nn, i, j, k) * (nn - 1.0_dp) * (nn - 2.0_dp) * (nn - 3.0_dp) * dzeta_r(3)**(nn - 4)
                    end if
                    if (temperature_der .eqv. .true.) then
                        jj3_dd_dtt(i, j, k) = jj3_dd_dtt(i, j, k) &
                                & + c__dd(nn, i, j, k) * dzeta_r(3)**(nn - 2) * (nn - 1.0_dp) * dzeta_dtt_r(3)
                    end if
                end do
            end do
        end do
    end do ! the sumation nn


    call write_array_2d (jj2_dd,"jj2_dd             :", stddeb) ! DEBUG control computation
    call write_array_3d (jj2_dd_dx,"jj2_dd_dx          :", stddeb) ! DEBUG control computation

    call write_array_3d (jj3_dd,"jj3_dd             :", stddeb) ! DEBUG control computation
    call write_array_3d (jj3_dd_drho,"jj3_dd_drho        :", stddeb) ! DEBUG control computation
    call write_array_3d (jj3_dd_dx(:,:,:,1),"jj3_dd_dx(:,:,:,1) :", stddeb) ! DEBUG control computation
    call write_array_3d (jj3_dd_dx(:,:,:,2),"jj3_dd_dx(:,:,:,2) :", stddeb) ! DEBUG control computation
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    !! apply the symmetry and multiplication by dzeta(3)**n for n-th derivative
    do i = 1, n_comp
        jj2_dd_drho(i, i) = jj2_dd_drho(i, i) * dzeta(3) ! multiplication problem
        jj2_dd_drho2(i, i) = jj2_dd_drho2(i, i) * dzeta(3)**2 ! multiplication problem

        jj3_dd_drho(i, i, i) = jj3_dd_drho(i, i, i) * dzeta(3) ! multiplication problem
        jj3_dd_drho2(i, i, i) = jj3_dd_drho2(i, i, i) * dzeta(3)**2 ! multiplication problem
        if (sec_der .eqv. .true.) then
            jj2_dd_drho3(i, i) = jj2_dd_drho3(i, i) * dzeta(3)**3 ! multiplication problem
            jj3_dd_drho3(i, i, i) = jj3_dd_drho3(i, i, i) * dzeta(3)**3 ! multiplication problem
        end if
        do j = i + 1, n_comp
            ! jj2 section
            jj2_dd(j, i) = jj2_dd(i, j)
            jj2_dd_drho(i, j) = jj2_dd_drho(i, j) * dzeta(3) ! multiplication problem
            jj2_dd_drho(j, i) = jj2_dd_drho(i, j)
            jj2_dd_drho2(i, j) = jj2_dd_drho2(i, j) * dzeta(3)**2 ! multiplication problem
            jj2_dd_drho2(j, i) = jj2_dd_drho2(i, j)
            ! jj3 section with iij terms
            jj3_dd(i, j, i) = jj3_dd(i, i, j)
            jj3_dd(j, i, i) = jj3_dd(i, i, j)
            jj3_dd_drho(i, i, j) = jj3_dd_drho(i, i, j) * dzeta(3) ! multiplication problem
            jj3_dd_drho(i, j, i) = jj3_dd_drho(i, i, j)
            jj3_dd_drho(j, i, i) = jj3_dd_drho(i, i, j)
            jj3_dd_drho2(i, i, j) = jj3_dd_drho2(i, i, j) * dzeta(3)**2 ! multiplication problem
            jj3_dd_drho2(i, j, i) = jj3_dd_drho2(i, i, j)
            jj3_dd_drho2(j, i, i) = jj3_dd_drho2(i, i, j)
            ! jj3 section with ijj terms
            jj3_dd(j, i, j) = jj3_dd(i, j, j)
            jj3_dd(j, j, i) = jj3_dd(i, j, j)
            jj3_dd_drho(i, j, j) = jj3_dd_drho(i, j, j) * dzeta(3) ! multiplication problem
            jj3_dd_drho(j, i, j) = jj3_dd_drho(i, j, j)
            jj3_dd_drho(j, j, i) = jj3_dd_drho(i, j, j)
            jj3_dd_drho2(i, j, j) = jj3_dd_drho2(i, j, j) * dzeta(3)**2 ! multiplication problem
            jj3_dd_drho2(j, i, j) = jj3_dd_drho2(i, j, j)
            jj3_dd_drho2(j, j, i) = jj3_dd_drho2(i, j, j)

            ! the second density derivative
            if (sec_der .eqv. .true.) then
                jj2_dd_drho3(i, j) = jj2_dd_drho3(i, j) * dzeta(3)**3 ! multiplication problem
                jj2_dd_drho3(j, i) = jj2_dd_drho3(i, j)

                jj3_dd_drho3(i, i, j) = jj3_dd_drho3(i, i, j) * dzeta(3)**3 ! multiplication problem
                jj3_dd_drho3(i, j, i) = jj3_dd_drho3(i, i, j)
                jj3_dd_drho3(j, i, i) = jj3_dd_drho3(i, i, j)

                jj3_dd_drho3(i, j, j) = jj3_dd_drho3(i, j, j) * dzeta(3)**3 ! multiplication problem
                jj3_dd_drho3(j, i, j) = jj3_dd_drho3(i, j, j)
                jj3_dd_drho3(j, j, i) = jj3_dd_drho3(i, j, j)
            end if
            ! the temperature derivative
            if (temperature_der .eqv. .true.) then
                jj2_dd_dtt(j, i) = jj2_dd_dtt(i, j)

                jj3_dd_dtt(i, j, i) = jj3_dd_dtt(i, i, j)
                jj3_dd_dtt(j, i, i) = jj3_dd_dtt(i, i, j)

                jj3_dd_dtt(j, i, j) = jj3_dd_dtt(i, j, j)
                jj3_dd_dtt(j, j, i) = jj3_dd_dtt(i, j, j)
            end if
            do k = j + 1, n_comp
                jj3_dd(i, k, j) = jj3_dd(i, j, k)
                jj3_dd(k, i, j) = jj3_dd(i, j, k)
                jj3_dd(k, j, i) = jj3_dd(i, j, k)
                jj3_dd(j, i, k) = jj3_dd(i, j, k)
                jj3_dd(j, k, i) = jj3_dd(i, j, k)
                jj3_dd_drho(i, j, k) = jj3_dd_drho(i, j, k) * dzeta(3) ! multiplication problem
                jj3_dd_drho(i, k, j) = jj3_dd_drho(i, j, k)
                jj3_dd_drho(k, i, j) = jj3_dd_drho(i, j, k)
                jj3_dd_drho(k, j, i) = jj3_dd_drho(i, j, k)
                jj3_dd_drho(j, i, k) = jj3_dd_drho(i, j, k)
                jj3_dd_drho(j, k, i) = jj3_dd_drho(i, j, k)
                jj3_dd_drho2(i, j, k) = jj3_dd_drho2(i, j, k) * dzeta(3)**2 ! multiplication problem
                jj3_dd_drho2(i, k, j) = jj3_dd_drho2(i, j, k)
                jj3_dd_drho2(k, i, j) = jj3_dd_drho2(i, j, k)
                jj3_dd_drho2(k, j, i) = jj3_dd_drho2(i, j, k)
                jj3_dd_drho2(j, i, k) = jj3_dd_drho2(i, j, k)
                jj3_dd_drho2(j, k, i) = jj3_dd_drho2(i, j, k)
                if (sec_der .eqv. .true.) then
                    jj3_dd_drho3(i, j, k) = jj3_dd_drho3(i, j, k) * dzeta(3)**3 ! multiplication problem
                    jj3_dd_drho3(i, k, j) = jj3_dd_drho3(i, j, k)
                    jj3_dd_drho3(k, i, j) = jj3_dd_drho3(i, j, k)
                    jj3_dd_drho3(k, j, i) = jj3_dd_drho3(i, j, k)
                    jj3_dd_drho3(j, i, k) = jj3_dd_drho3(i, j, k)
                    jj3_dd_drho3(j, k, i) = jj3_dd_drho3(i, j, k)
                end if
                if (temperature_der .eqv. .true.) then
                    jj3_dd_dtt(i, k, j) = jj3_dd_dtt(i, j, k)
                    jj3_dd_dtt(k, i, j) = jj3_dd_dtt(i, j, k)
                    jj3_dd_dtt(k, j, i) = jj3_dd_dtt(i, j, k)
                    jj3_dd_dtt(j, i, k) = jj3_dd_dtt(i, j, k)
                    jj3_dd_dtt(j, k, i) = jj3_dd_dtt(i, j, k)
                end if
            end do
        end do
    end do

    mol_rat2 = 0.0_dp
    mol_rat3 = 0.0_dp
    mol_rat2_dx = -1.0_dp
    mol_rat3_dx = -1.0_dp
    ! composition defivative section which utilize similar shape of density derivation
    ! special treatment is required because o different composition-related derivatives
    ! NOTE the method was not succesfully included into symmetrical computation of previous jj2 jj3 integral due to above-mentioned problem
    if (composition_der .eqv. .true.) then
        do i = 1, n_comp
            jj2_dd_dx(:, :, i) = jj2_dd_drho(:, :) / dzeta(3) * dzeta_dx_r(3, i)
            jj3_dd_dx(:, :, :, i) = jj3_dd_drho(:, :, :) / dzeta(3) * dzeta_dx_r(3, i)
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

    ! write (stddeb,"('ees3s3_s3nnmm       : ',e)"), ees3s3_s3nnmm(1,1) ! DEBUG control computation
    ! write (stddeb,"('eees3s3s3_sssnnnmmm : ',e)"), eees3s3s3_sssnnnmmm(1,1,1) ! DEBUG control of computation
    call write_array_2d (ees3s3_s3nnmm,"ees3s3_s3nnmm       :",stddeb) ! DEBUG control computation
    call write_array_3d (eees3s3s3_sssnnnmmm,"eees3s3s3_sssnnnmmm :",stddeb) ! DEBUG control computation
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    ! a2, a3 terms computation with derivatives
    ! based upon paper Gross-Vrabec 2006 application of eq (8, 9)
    do i = 1, n_comp
        a2 = a2 &
                & + mol_rat(i)**2 * ees3s3_s3nnmm(i, i) * jj2_dd(i, i)
        a2_drho = a2_drho &
                & + mol_rat(i)**2 * ees3s3_s3nnmm(i, i) * jj2_dd(i, i) &
                & + rho_in * (mol_rat(i)**2 * ees3s3_s3nnmm(i, i) * jj2_dd_drho(i, i))
        a2_drho2 = a2_drho2 &
                & + 2.0_dp * (mol_rat(i)**2 * ees3s3_s3nnmm(i, i) * jj2_dd_drho(i, i)) &
                & + rho_in * (mol_rat(i)**2 * ees3s3_s3nnmm(i, i) * jj2_dd_drho2(i, i))

        a3 = a3 &
                & + mol_rat(i)**3 * eees3s3s3_sssnnnmmm(i, i, i) * jj3_dd(i, i, i)
        a3_drho = a3_drho &
                & + 2.0_dp * (mol_rat(i)**3 * eees3s3s3_sssnnnmmm(i, i, i) * jj3_dd(i, i, i)) &
                & + rho_in * (mol_rat(i)**3 * eees3s3s3_sssnnnmmm(i, i, i) * jj3_dd_drho(i, i, i))
        a3_drho2 = a3_drho2 &
                & + 2.0_dp * (mol_rat(i)**3 * eees3s3s3_sssnnnmmm(i, i, i) * jj3_dd(i, i, i)) &
                & + 4.0_dp * rho_in * (mol_rat(i)**3 * eees3s3s3_sssnnnmmm(i, i, i) * jj3_dd_drho(i, i, i)) &
                & + rho_in**2 * (mol_rat(i)**3 * eees3s3s3_sssnnnmmm(i, i, i) * jj3_dd_drho2(i, i, i))

        if (sec_der .eqv. .true.) then
            a2_drho3 = a2_drho3 &
                    & + 3.0_dp * (mol_rat(i)**2 * ees3s3_s3nnmm(i, i) * jj2_dd_drho2(i, i)) &
                    & + rho_in * (mol_rat(i)**2 * ees3s3_s3nnmm(i, i) * jj2_dd_drho3(i, i))
            a3_drho3 = a3_drho3 &
                    & + 6.0_dp * (mol_rat(i)**3 * eees3s3s3_sssnnnmmm(i, i, i) * jj3_dd_drho(i, i, i)) &
                    & + 6.0_dp * rho_in * (mol_rat(i)**3 * eees3s3s3_sssnnnmmm(i, i, i) * jj3_dd_drho2(i, i, i)) &
                    & + rho_in**2 * (mol_rat(i)**3 * eees3s3s3_sssnnnmmm(i, i, i) * jj3_dd_drho3(i, i, i))
        end if
        if (temperature_der .eqv. .true.) then
            a2_dtt = a2_dtt &
                    & + 2.0_dp / tt * (mol_rat(i)**2 * ees3s3_s3nnmm(i, i) * jj2_dd(i, i)) &
                    & - (mol_rat(i)**2 * ees3s3_s3nnmm(i, i) * jj2_dd_dtt(i, i))
            a3_dtt = a3_dtt &
                    & + 4.0_dp / tt * (mol_rat(i)**3 * eees3s3s3_sssnnnmmm(i, i, i) * jj3_dd(i, i, i)) &
                    & - 4.0_dp / 3.0_dp * (mol_rat(i)**3 * eees3s3s3_sssnnnmmm(i, i, i) * jj3_dd_dtt(i, i, i))
        end if
        do j = i + 1, n_comp
            a2 = a2 &
                    & + 2.0_dp * mol_rat(i) * mol_rat(j) * ees3s3_s3nnmm(i, j) * jj2_dd(i, j)
            a2_drho = a2_drho &
                    & + 2.0_dp * (mol_rat(i) * mol_rat(j) * ees3s3_s3nnmm(i, j) * jj2_dd(i, j) &
                            & + rho_in * (mol_rat(i) * mol_rat(j) * ees3s3_s3nnmm(i, j) * jj2_dd_drho(i, j)))
            a2_drho2 = a2_drho2 &
                    & + 2.0_dp * (2.0_dp * (mol_rat(i) * mol_rat(j) * ees3s3_s3nnmm(i, j) * jj2_dd_drho(i, j)) &
                            & + rho_in * (mol_rat(i) * mol_rat(j) * ees3s3_s3nnmm(i, j) * jj2_dd_drho2(i, j)))

            a3 = a3 &
                    & + 3.0_dp * mol_rat(i) * mol_rat(j)**2 * eees3s3s3_sssnnnmmm(i, j, j) * jj3_dd(i, j, j)
            a3_drho = a3_drho &
                    & + 3.0_dp * (2.0_dp * (mol_rat(i) * mol_rat(j)**2 * eees3s3s3_sssnnnmmm(i, j, j) * jj3_dd(i, j, j)) &
                            & + rho_in * (mol_rat(i) * mol_rat(j)**2 * eees3s3s3_sssnnnmmm(i, j, j) * jj3_dd_drho(i, j, j)))
            a3_drho2 = a3_drho2 &
                    & + 3.0_dp * (2.0_dp * (mol_rat(i) * mol_rat(j)**2 * eees3s3s3_sssnnnmmm(i, j, j) * jj3_dd(i, j, j)) &
                            & + 4.0_dp * rho_in * (mol_rat(i) * mol_rat(j)**2 * eees3s3s3_sssnnnmmm(i, j, j) * jj3_dd_drho(i, j, j)) &
                            & + rho_in**2 * (mol_rat(i) * mol_rat(j)**2 * eees3s3s3_sssnnnmmm(i, j, j) * jj3_dd_drho2(i, j, j)))

            a3 = a3 &
                    & + 3.0_dp * mol_rat(j) * mol_rat(i)**2 * eees3s3s3_sssnnnmmm(i, i, j) * jj3_dd(i, i, j)
            a3_drho = a3_drho &
                    & + 3.0_dp * (2.0_dp * (mol_rat(j) * mol_rat(i)**2 * eees3s3s3_sssnnnmmm(i, i, j) * jj3_dd(i, i, j)) &
                            & + rho_in * (mol_rat(j) * mol_rat(i)**2 * eees3s3s3_sssnnnmmm(i, i, j) * jj3_dd_drho(i, i, j)))
            a3_drho2 = a3_drho2 &
                    & + 3.0_dp * (2.0_dp * (mol_rat(j) * mol_rat(i)**2 * eees3s3s3_sssnnnmmm(i, i, j) * jj3_dd(i, i, j)) &
                            & + 4.0_dp * rho_in * (mol_rat(j) * mol_rat(i)**2 * eees3s3s3_sssnnnmmm(i, i, j) * jj3_dd_drho(i, i, j)) &
                            & + rho_in**2 * (mol_rat(j) * mol_rat(i)**2 * eees3s3s3_sssnnnmmm(i, i, j) * jj3_dd_drho2(i, i, j)))

            if (sec_der .eqv. .true.) then
                a2_drho3 = a2_drho3 &
                        & + 2.0_dp * (3.0_dp * (mol_rat(i) * mol_rat(j) * ees3s3_s3nnmm(i, j) * jj2_dd_drho2(i, j)) &
                                & + rho_in * (mol_rat(i) * mol_rat(j) * ees3s3_s3nnmm(i, j) * jj2_dd_drho3(i, j)))
                a3_drho3 = a3_drho3 &
                        & + 3.0_dp * (6.0_dp * (mol_rat(i) * mol_rat(j)**2 * eees3s3s3_sssnnnmmm(i, j, j) * jj3_dd_drho(i, j, j)) &
                                & + 6.0_dp * rho_in * (mol_rat(i) * mol_rat(j)**2 * eees3s3s3_sssnnnmmm(i, j, j) * jj3_dd_drho2(i, j, j)) &
                                & + rho_in**2 * (mol_rat(i) * mol_rat(j)**2 * eees3s3s3_sssnnnmmm(i, j, j) * jj3_dd_drho3(i, j, j)))

                a3_drho3 = a3_drho3 &
                        & + 3.0_dp * (6.0_dp * (mol_rat(j) * mol_rat(i)**2 * eees3s3s3_sssnnnmmm(i, i, j) * jj3_dd_drho(i, i, j)) &
                                & + 6.0_dp * rho_in * (mol_rat(j) * mol_rat(i)**2 * eees3s3s3_sssnnnmmm(i, i, j) * jj3_dd_drho2(i, i, j)) &
                                & + rho_in**2 * (mol_rat(j) * mol_rat(i)**2 * eees3s3s3_sssnnnmmm(i, i, j) * jj3_dd_drho3(i, i, j)))
            end if
            if (temperature_der .eqv. .true.) then
                a2_dtt = a2_dtt &
                        & + 2.0_dp * (2.0_dp / tt * (mol_rat(i) * mol_rat(j) * ees3s3_s3nnmm(i, j) * jj2_dd(i, j)) &
                                & - (mol_rat(i) * mol_rat(j) * ees3s3_s3nnmm(i, j) * jj2_dd_dtt(i, j)))

                a3_dtt = a3_dtt &
                        & + 3.0_dp * (4.0_dp / tt * (mol_rat(i) * mol_rat(j)**2 * eees3s3s3_sssnnnmmm(i, j, j) * jj3_dd(i, j, j)) &
                                & - 4.0_dp / 3.0_dp * (mol_rat(i) * mol_rat(j)**2 * eees3s3s3_sssnnnmmm(i, j, j) * jj3_dd_dtt(i, j, j)))

                a3_dtt = a3_dtt &
                        & + 3.0_dp * (4.0_dp / tt * (mol_rat(j) * mol_rat(i)**2 * eees3s3s3_sssnnnmmm(i, i, j) * jj3_dd(i, i, j)) &
                                & - 4.0_dp / 3.0_dp * (mol_rat(j) * mol_rat(i)**2 * eees3s3s3_sssnnnmmm(i, i, j) * jj3_dd_dtt(i, i, j)))
            end if
            do k = j + 1, n_comp
                a3 = a3 &
                        & + 6.0_dp * mol_rat(i) * mol_rat(j) * mol_rat(k) * eees3s3s3_sssnnnmmm(i, j, k) * jj3_dd(i, j, k)
                a3_drho = a3_drho &
                        & + 6.0_dp * (2.0_dp * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees3s3s3_sssnnnmmm(i, j, k) * jj3_dd(i, j, k)) &
                                & + rho_in * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees3s3s3_sssnnnmmm(i, j, k) * jj3_dd_drho(i, j, k)))
                a3_drho2 = a3_drho2 &
                        & + 6.0_dp * (2.0_dp * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees3s3s3_sssnnnmmm(i, j, k) * jj3_dd(i, j, k)) &
                                & + 4.0_dp * rho_in * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees3s3s3_sssnnnmmm(i, j, k) * jj3_dd_drho(i, j, k)) &
                                & + rho_in**2 * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees3s3s3_sssnnnmmm(i, j, k) * jj3_dd_drho2(i, j, k)))
                if (sec_der .eqv. .true.) then
                    a3_drho3 = a3_drho3 &
                            & + 6.0_dp * (6.0_dp * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees3s3s3_sssnnnmmm(i, j, k) * jj3_dd_drho(i, j, k)) &
                                    & + 6.0_dp * rho_in * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees3s3s3_sssnnnmmm(i, j, k) * jj3_dd_drho2(i, j, k)) &
                                    & + rho_in**2 * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees3s3s3_sssnnnmmm(i, j, k) * jj3_dd_drho3(i, j, k)))
                end if
                if (temperature_der .eqv. .true.) then
                    a3_dtt = a3_dtt &
                            & + 6.0_dp * (4.0_dp / tt * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees3s3s3_sssnnnmmm(i, j, k) * jj3_dd(i, j, k)) &
                                    & - 4.0_dp / 3.0_dp * (mol_rat(i) * mol_rat(j) * mol_rat(k) * eees3s3s3_sssnnnmmm(i, j, k) * jj3_dd_dtt(i, j, k)))
                end if
            end do
        end do
    end do
    ! 	write (*,"('rho_in: ',e)"), rho_in ! DEBUG control computation

    ! second section of special treatment of composition derivatives
    ! utilizing the vectorization to save do loop usage (and readability)
    ! NOTE this is also reason for tedious preparation in previous composition derivative section
    if (composition_der .eqv. .true.) then
        do i = 1, n_comp
            a2_dx(i) = sum(mol_rat2_dx(:, :, i) * ees3s3_s3nnmm(:, :) * jj2_dd(:, :) &
                    & + mol_rat2(:, :) * ees3s3_s3nnmm(:, :) * jj2_dd_dx(:, :, i))
            a3_dx(i) = sum(mol_rat3_dx(:, :, :, i) * eees3s3s3_sssnnnmmm(:, :, :) * jj3_dd(:, :, :) &
                    & + mol_rat3(:, :, :) * eees3s3s3_sssnnnmmm(:, :, :) * jj3_dd_dx(:, :, :, i))
        end do
    end if
    a2 = a2 * (-pi * rho_in / tt**2)
    a2_drho = a2_drho * (-pi / tt**2)
    a2_drho2 = a2_drho2 * (-pi / tt**2)
    a2_drho3 = a2_drho3 * (-pi / tt**2)
    a2_dtt = a2_dtt * (pi * rho_in / tt**2)
    a2_dx = a2_dx * (-pi * rho_in / tt**2)

    a3 = a3 * (-4.0_dp / 3.0_dp * pi**2 * rho_in**2 / tt**3)
    a3_drho = a3_drho * (-4.0_dp / 3.0_dp * pi**2 * rho_in / tt**3)
    a3_drho2 = a3_drho2 * (-4.0_dp / 3.0_dp * pi**2 / tt**3)
    a3_drho3 = a3_drho3 * (-4.0_dp / 3.0_dp * pi**2 / tt**3)
    a3_dtt = a3_dtt * (pi**2 * rho_in**2 / tt**3)
    a3_dx = a3_dx * (-4.0_dp / 3.0_dp * pi**2 * rho_in**2 / tt**3)

    write(stddeb,"('a2      : ',e16.8)") a2 ! DEBUG control computation
    write(stddeb,"('a2_drho : ',e16.8)") a2_drho ! DEBUG control computation
    write(stddeb,"('a2_dtt  : ',e16.8)") a2_dtt ! DEBUG control computation
    call write_array_1d(a2_dx,"a2_dx   : ", stddeb) ! DEBUG control of computation

    write(stddeb,"('a3      : ',e16.8)") a3 ! DEBUG control computation
    write(stddeb,"('a3_drho : ',e16.8)") a3_drho ! DEBUG control computation
    write(stddeb,"('a3_dtt  : ',e16.8)") a3_dtt ! DEBUG control computation
    call write_array_1d(a3_dx,"a3_dx   : ", stddeb) ! DEBUG control of computation
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    ! dipol-dipol interaction contribution computation
    ! based upon paper Gross-Vrabec 2006 application of eq (7)

    a_dd_ = a2**2 / (a2 - a3) ! TODO decide if keep and make the a_dd obsolete
    ! following equalities are derived by hand, by maple and controlled by wolfram alpha
    z_dd_ = 2.0_dp * a2_drho * a2 / (a2 - a3) &
            & - a2**2 * (a2_drho - a3_drho) / (a2 - a3)**2
    z_dd_drho_ = 2.0_dp * a2_drho2 * a2 / (a2 - a3) &
            & + 2.0_dp * a2_drho**2 / (a2 - a3) &
            & - 4.0_dp * a2_drho * a2 * (a2_drho - a3_drho) / (a2 - a3)**2 &
            & - a2**2 * (a2_drho2 - a3_drho2) / (a2 - a3)**2 &
            & + 2.0_dp * a2**2 * (a2_drho - a3_drho)**2 / (a2 - a3)**3
    if (sec_der .eqv. .true.) then
        z_dd_drho2_ = 2.0_dp * a2_drho3 * a2 / (a2 - a3) &
                & + 6.0_dp * a2_drho2 * a2_drho / (a2 - a3) &
                & - 6.0_dp * a2_drho2 * a2 * (a2_drho - a3_drho) / (a2 - a3)**2 &
                & - 6.0_dp * a2_drho**2 * (a2_drho - a3_drho) / (a2 - a3)**2 &
                & - 6.0_dp * a2_drho * a2 * (a2_drho2 - a3_drho2) / (a2 - a3)**2 &
                & + 12.0_dp * a2_drho * a2 * (a2_drho - a3_drho)**2 / (a2 - a3)**3 &
                & - a2**2 * (a2_drho3 - a3_drho3) / (a2 - a3)**2 &
                & + 6.0_dp * a2**2 * (a2_drho - a3_drho) * (a2_drho2 - a3_drho2) / (a2 - a3)**3 &
                & - 6.0_dp * a2**2 * (a2_drho - a3_drho)**3 / (a2 - a3)**4
    end if
    if (temperature_der .eqv. .true.) then
        a_dd_dtt = (a2**2 * a2_dtt &
                & + a2**2 * a3_dtt &
                & - 2.0_dp * a2 * a3 * a2_dtt) / (a2 - a3)**2
    end if
    if (composition_der .eqv. .true.) then
        a_dd_dx(:) = (a2**2 * a2_dx(:) &
                & + a2**2 * a3_dx(:) &
                & - 2.0_dp * a2 * a3 * a2_dx(:)) / (a2 - a3)**2
    end if

    write (stddeb,*) "z_dd_      : ", z_dd_ ! DEBUG just for brevity of debug
    write (stddeb,*) "z_dd_drho_ : ", z_dd_drho_ ! DEBUG just for brevity of debug
    if (sec_der .eqv. .true.) then
        write (stddeb,*) "z_dd_drho2_: ", z_dd_drho2_ ! DEBUG just for brevity of debug
    end if

    write (stddeb,*) "a_dd_      : ", a_dd_ ! DEBUG just for brevity of debug
    if (temperature_der .eqv. .true.) then
        write (stddeb,*) "a_dd_dtt   : ", a_dd_dtt ! DEBUG just for brevity of debug
    end if
    if (composition_der .eqv. .true.) then
        write (stddeb,*) "a_dd_dx    : ", a_dd_dx ! DEBUG just for brevity of debug
    end if
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug


    a_dd = a_dd_
    ! NOTE a_dd_dtt direct modification in previous section
    ! NOTE a_dd_dx direct modification in previous section

    z_dd = rho_in * z_dd_
    z_dd_drho = z_dd_ + rho_in * z_dd_drho_
    z_dd_drho2 = 2.0_dp * z_dd_drho_ + rho_in * z_dd_drho2_

    init_flag_add = .true.
    init_flag_zdd = .true.

    ! ----- deallocation section -----
    deallocate(a2_dx, a3_dx)

    deallocate(dzeta_dx_r)

    deallocate(jj2_dd, jj2_dd_drho, jj2_dd_drho2, jj2_dd_drho3, jj2_dd_dtt)
    deallocate(jj2_dd_dx)

    deallocate(jj3_dd, jj3_dd_drho, jj3_dd_drho2, jj3_dd_drho3, jj3_dd_dtt)
    deallocate(jj3_dd_dx)

    deallocate(mol_rat2)
    deallocate(mol_rat2_dx)

    deallocate(mol_rat3)
    deallocate(mol_rat3_dx)

    write (stddeb,*) "<== z_dd_sub" !DEBUG
end subroutine z_dd_sub
