subroutine rdf_sub (rho_in)
    ! compute the radial distribution function g_{i,j} and its derivatives for a single density
    ! - prepare the necesarry components such as dzeta function
    !	* calling only when necesarry
    ! - main loop throught the component number
    !	* loop is divided into two sections
    !		** first compute the diagonal terms
    !			*** these can be simplified
    !		** second compute the nondiagonal terms
    !			*** using the symmetry only half is computed
    !	* formulas were reviewed by maple program
    !
    !	== LAST MODIFICATIONS ==
    !	15.5.2017 - D.C - descriptions and comments
    !	20.6.2017 - D.C - init_flag_g dependency omitment -> moved to parent control module -> replaced by single flag

    use control_mod, only : dp, n_comp, sec_der, composition_der, temperature_der, &
            & stdout, stderr, stddeb, stdlog, &
            & write_array_1d, write_array_2d, write_array_3d, &
            & g_ij, g_ij_drho, g_ij_drho2, g_ij_drho3, & !g_ij_drho4, g_ij_drho5, & ! OLD presently unused
            & g_ij_dx, g_ij_dtt, &
            & init_flag, init_flag_gij
    use contrib_mod, only : d_seg, d_seg_dtt, dzeta, dzeta_dx, dzeta_dtt, &
            & initialize_dzeta_fun
    ! NOTE real(dp), dimension(0:3) :: dzeta ! the dzeta vector dzeta_n for n=0,1,2,3 without rho !!
    ! NOTE real(dp), dimension(n_comp) :: d_seg ! the segment diameter di (temperature dependant)
    implicit none
    ! variable section IN/OUT
    real(dp), intent(in) :: rho_in ! the input density
    ! local variables outputed into control_mod
    real(dp), dimension(:, :), allocatable :: g_ij_ ! the output argument
    real(dp), dimension(:, :), allocatable :: g_ij_drho_ ! the output argument first density derivative
    real(dp), dimension(:, :), allocatable :: g_ij_drho2_ ! the output argument second density derivative
    real(dp), dimension(:, :), allocatable :: g_ij_drho3_ ! the output argument third density derivative
    ! OLD presently unused
    ! 	real(dp), dimension(:,:), allocatable :: g_ij_drho5_ ! the output argument third density derivative
    ! 	real(dp), dimension(:,:), allocatable :: g_ij_drho4_ ! the output argument third density derivative

    ! compostion derivative has bigger dimension to accomodate derivatives with respect to components
    real(dp), dimension(:, :, :), allocatable :: g_ij_dx_ ! the output argument first composition derivative
    real(dp), dimension(:, :), allocatable :: g_ij_dtt_ ! the output argument first temperature derivative

    ! local variables
    real(dp), dimension(0:3) :: dzeta_r ! the dzeta vector dzeta_n for n=0,1,2,3 * rho !!
    real(dp), dimension(:, :), allocatable :: dzeta_dx_r ! the matrix ddzeta_n/dxi for (n=0,1,2,3 x n_comp) * rho !!
    real(dp), dimension(0:3) :: dzeta_dtt_r ! the matrix ddzeta_n/dT for (n=0,1,2,3) * rho !!
    real(dp) :: zms ! the VV name for 1 - dzeta_r(3)
    integer :: i, j, k ! counters
    real(dp) :: d ! the temporary variable that holds (di*dj)/(di+dj) !! is function of T !!

    !------- PREPARING ------
    write (stddeb,*) "===> rdf_sub" !DEBUG
    ! call for the required helper procedures
    if (init_flag(2) .eqv. .false.) then
        ! 		write (*,*), 'Call the initialize_dzeta_fun from rdf_sub' ! DEBUG control of computation
        if (initialize_dzeta_fun() .eqv. .false.) then
            write (stderr, *) 'Error in helper function initialize_dzeta_fun'
            stop
        end if
    end if

    ! ----- allocation section -----
    allocate(g_ij_(n_comp, n_comp))
    allocate(g_ij_drho_(n_comp, n_comp))
    allocate(g_ij_drho2_(n_comp, n_comp))
    allocate(g_ij_drho3_(n_comp, n_comp))
    !	allocate( g_ij_drho4_(n_comp,n_comp))
    !	allocate( g_ij_drho5_(n_comp,n_comp))

    allocate(g_ij_dx_(n_comp, n_comp, n_comp))
    allocate(g_ij_dtt_(n_comp, n_comp))

    allocate(dzeta_dx_r(0:3, 1:n_comp))

    ! the actual dzeta (multiplied by density)
    dzeta_r = dzeta * rho_in
    dzeta_dx_r = dzeta_dx * rho_in
    dzeta_dtt_r = dzeta_dtt * rho_in
    zms = 1.0_dp - dzeta_r(3)

    write (stddeb,*) '-----  RDF  -----'  ! DEBUG just for brevity of debug
    call write_array_1d (dzeta,"dzeta      : ",stddeb)
    call write_array_1d (dzeta_r,"dzeta_r    : ",stddeb)

    call write_array_1d (dzeta_dtt,"dzeta_dtt  : ",stddeb)
    call write_array_1d (dzeta_dtt_r,"dzeta_dtt_r: ",stddeb)

    call write_array_2d (dzeta_dx,"dzeta_dx   : ",stddeb)
    call write_array_2d (dzeta_dx_r,"dzeta_dx_r : ",stddeb)
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    ! the computation loop optimized for symmetry in di dj
    do i = 1, n_comp
        d = d_seg(i) / 2.0_dp
        ! 		write (*,"('d(',i2,',',i2,'): ',f12.8)"), i, i, d ! DEBUG test the computation
        ! ----- density derivative section - diagonal -----
        ! BUG g_ij is inconsistent with the paper Gross-sadowski-2001 PC-SAFT with eq (8) -error in paper
        ! NOTE this is Carnahan-Starling RDF for hard sphere mixtures
        g_ij_(i, i) = 1.0_dp / zms &
                & + (3.0_dp * d * dzeta_r(2)) / zms**2 &
                & +(2.0_dp * (dzeta_r(2) * d)**2) / zms**3 ! expressed by maple

        g_ij_drho_(i, i) = (dzeta(3) &
                & + 3.0_dp * d * dzeta(2)) / zms**2 &
                & +(6.0_dp * d * dzeta_r(2) * dzeta(3) &
                        & + 4.0_dp * dzeta(2) * dzeta_r(2) * d**2) / zms**3 &
                & + (6.0_dp * dzeta(3) * (dzeta_r(2) * d)**2) / zms**4 ! derived by maple

        g_ij_drho2_(i, i) = (2.0_dp * dzeta(3)**2 &
                & + 12.0_dp * d * dzeta(2) * dzeta(3) &
                & + 4.0_dp * (d * dzeta(2))**2) / zms**3  &
                & + (18.0_dp * d * dzeta_r(2) * dzeta(3)**2 &
                        & + 24.0_dp * dzeta_r(3) * (dzeta(2) * d)**2) / zms**4 &
                & + (24.0_dp * (dzeta_r(2) * dzeta(3) * d)**2) / zms**5 ! derived by maple

        if (sec_der .eqv. .true.) then
            g_ij_drho3_(i, i) = (6.0_dp * dzeta(3)**3 &
                    & + 54.0_dp * d * dzeta(2) * dzeta(3)**2 &
                    & + 36.0_dp * dzeta(3) * (dzeta(2) * d)**2) / zms**4 &
                    & +(72.0_dp * d * dzeta_r(2) * dzeta(3)**3 &
                            & + 144.0_dp * dzeta(2) * dzeta_r(2) * (dzeta(3) * d)**2) / zms**5 &
                    & + (120.0_dp * dzeta(3) * (d * dzeta_r(2) * dzeta(3))**2) / zms**6 ! derived by maple
            ! OLD presently unused
            ! 			g_ij_drho4_(i,i)=( 288.0_dp*(dzeta(2)*dzeta(3)*d)**2 &
            ! 						  &	  +288.0_dp*dzeta(2)*d*dzeta(3)**3 &
            ! 						  &	   +24.0_dp*dzeta(3)**4 )/zms**5 &
            ! 						  & +( 960.0_dp*dzeta_r(3)*(dzeta(2)*dzeta(3)*d)**2 &
            ! 						  &	  +360.0_dp*dzeta_r(2)*d*dzeta(3)**4 )/zms**6 &
            ! 						  & +( 720.0_dp*(dzeta_r(2)*d*dzeta(3)**2)**2 )/zms**7 ! derived by maple
            !
            ! 			g_ij_drho5_(i,i)=( 2400.0_dp*dzeta(3)*(dzeta(2)*dzeta(3)*d)**2 &
            ! 						  &	  +1800.0_dp*dzeta(2)*d*dzeta(3)**4 &
            ! 						  &	   +120.0_dp*dzeta(3)**5 )/zms**6 &
            ! 						  & +( 7200.0_dp*dzeta_r(2)*dzeta(2)*d**2 *dzeta(3)**4 &
            ! 						  &	  +2400.0_dp*dzeta(3)*(dzeta(2)*dzeta(3)*d)**2 &
            ! 						  &	  +2160.0_dp*dzeta_r(2)*d*dzeta(3)**5 )/zms**7 &
            ! 						  & +( 5040.0_dp*dzeta(3)*(dzeta_r(2)*d*dzeta(3)**2)**2 )/zms**8 ! derived by maple
        end if

        ! ----- compostition derivative section - diagonal -----
        if (composition_der .eqv. .true.) then
            do k = 1, n_comp
                g_ij_dx_(i, i, k) = (dzeta_dx_r(3, k) &
                        & + 3.0_dp * d * dzeta_dx_r(2, k)) / zms**2 &
                        & + (6.0_dp * d * dzeta_r(2) * dzeta_dx_r(3, k) &
                                & + 4.0_dp * d**2 * dzeta_r(2) * dzeta_dx_r(2, k)) / zms**3 &
                        & + (6.0_dp * dzeta_dx_r(3, k) * (d * dzeta_r(2))**2) / zms**4 ! BUG beware the difference between _r and _dx_r
            end do
        end if
        ! ----- temperature derivative section - diagonal -----
        if (temperature_der .eqv. .true.) then
            ! NOTE we use the d_seg that is more convenient here than d
            g_ij_dtt_(i, i) = (2.0_dp * dzeta_dtt_r(3) &
                    & + 3.0_dp * d_seg_dtt(i) * dzeta_r(2) &
                    & + 3.0_dp * d_seg(i) * dzeta_dtt_r(2)) / (2.0_dp * zms**2) &
                    & + (3.0_dp * d_seg(i) * dzeta_r(2) * dzeta_dtt_r(3) &
                            & + d_seg(i) * d_seg_dtt(i) * dzeta_r(2)**2 &
                            & + dzeta_r(2) * dzeta_dtt_r(2) * d_seg(i)**2) / zms**3 &
                    & + (3.0_dp * dzeta_dtt_r(3) * (d_seg(i) * dzeta_r(2))**2) / (2.0_dp * zms**4)
        end if
        ! ===== nondiagonal computation =====
        do j = i + 1, n_comp
            d = (d_seg(i) * d_seg(j)) / (d_seg(i) + d_seg(j)) ! BUG check correctness of d
            ! ----- density derivative section - diagonal -----
            ! 			write (*,"('d(',i2,',',i2,'): ',f12.8)"), i, j, d ! DEBUG test the computation
            g_ij_(i, j) = 1.0_dp / zms &
                    & + (3.0_dp * d * dzeta_r(2)) / zms**2 &
                    &  +(2.0_dp * (dzeta_r(2) * d)**2) / zms**3 ! expressed by maple
            g_ij_(j, i) = g_ij_(i, j) ! use of the symmetry in di dj

            g_ij_drho_(i, j) = (dzeta(3) &
                    & + 3.0_dp * d * dzeta(2)) / zms**2 &
                    &   +(6.0_dp * d * dzeta_r(2) * dzeta(3) &
                            & + 4.0_dp * dzeta(2) * dzeta_r(2) * d**2) / zms**3 &
                    & + (6.0_dp * dzeta(3) * (dzeta_r(2) * d)**2) / zms**4 ! derived by maple
            g_ij_drho_(j, i) = g_ij_drho_(i, j) ! use of the symmetry in di dj

            g_ij_drho2_(i, j) = (2.0_dp * dzeta(3)**2  &
                    & + 12.0_dp * d * dzeta(2) * dzeta(3) &
                    & + 4.0_dp * (d * dzeta(2))**2) / zms**3 &
                    & + (18.0_dp * d * dzeta_r(2) * dzeta(3)**2 &
                            & + 24.0_dp * dzeta_r(3) * (dzeta(2) * d)**2) / zms**4 &
                    & + (24.0_dp * (dzeta_r(2) * dzeta(3) * d)**2) / zms**5 ! derived by maple
            g_ij_drho2_(j, i) = g_ij_drho2_(i, j) ! use of the symmetry in di dj

            if (sec_der .eqv. .true.) then
                g_ij_drho3_(i, j) = (6.0_dp * dzeta(3)**3 &
                        & + 54.0_dp * d * dzeta(2) * dzeta(3)**2 &
                        & + 36.0_dp * dzeta(3) * (dzeta(2) * d)**2) / zms**4 &
                        & + (72.0_dp * d * dzeta_r(2) * dzeta(3)**3 &
                                & + 144.0_dp * dzeta(2) * dzeta_r(2) * (dzeta(3) * d)**2) / zms**5 &
                        & + (120.0_dp * dzeta(3) * (d * dzeta_r(2) * dzeta(3))**2) / zms**6 ! derived by maple
                g_ij_drho3_(j, i) = g_ij_drho3_(i, j) ! use of the symmetry in di dj
                ! OLD presently unused
                ! 				g_ij_drho4_(i,j)=( 288.0_dp*(dzeta(2)*dzeta(3)*d)**2 &
                ! 							 &	  +288.0_dp*dzeta(2)*d*dzeta(3)**3 &
                ! 							 &	   +24.0_dp*dzeta(3)**4 )/zms**5 &
                ! 							 &  +( 960.0_dp*dzeta_r(3)*(dzeta(2)*dzeta(3)*d)**2 &
                ! 							 &	  +360.0_dp*dzeta_r(2)*d*dzeta(3)**4 )/zms**6 &
                ! 							 &  +( 720.0_dp*(dzeta_r(2)*d*dzeta(3)**2)**2 )/zms**7 ! derived by maple
                ! 				g_ij_drho4_(j,i) = g_ij_drho4_(i,j) ! use of the symmetry in di dj
                !
                ! 				g_ij_drho5_(i,j)=( 2400.0_dp*dzeta(3)*(dzeta(2)*dzeta(3)*d)**2 &
                ! 							 &	  +1800.0_dp*dzeta(2)*d*dzeta(3)**4 &
                ! 							 &	   +120.0_dp*dzeta(3)**5 )/zms**6 &
                ! 							 &  +( 7200.0_dp*dzeta_r(2)*dzeta(2)*d**2 *dzeta(3)**4 &
                ! 							 &	  +2400.0_dp*dzeta(3)*(dzeta(2)*dzeta(3)*d)**2 &
                ! 							 &	  +2160.0_dp*dzeta_r(2)*d*dzeta(3)**5 )/zms**7 &
                ! 							 &  +( 5040.0_dp*dzeta(3)*(dzeta_r(2)*d*dzeta(3)**2)**2 )/zms**8 ! derived by maple
                ! 				g_ij_drho5_(j,i) = g_ij_drho5_(i,j) ! use of the symmetry in di dj
            end if
            ! OLD nondiagonal terms of g_ij composition and temperature derivatives are not required
            ! 			! ----- compostition derivative section - nondiagonal -----
            if (composition_der .eqv. .true.) then
                do k = 1, n_comp
                    g_ij_dx_(i, j, k) = (dzeta_dx_r(3, k) &
                            & + 3.0_dp * d * dzeta_dx_r(2, k)) / zms**2 &
                            & + (6.0_dp * d * dzeta_r(2) * dzeta_dx_r(3, k) &
                                    & + 4.0_dp * d**2 * dzeta_r(2) * dzeta_dx_r(2, k)) / zms**3 &
                            & + (6.0_dp * dzeta_dx_r(3, k) * (d * dzeta_r(2))**2) / zms**4 ! BUG beware the difference between _r and _r_dx
                    g_ij_dx_(j, i, k) = g_ij_dx_(i, j, k)
                end do
            end if
            !
            ! 			! ----- temperature derivative section - nondiagonal -----
            !			! NOTE correct this implementation for nondiagonal terms if ever required
            ! 			if (temperature_der==.true.) then
            ! 			    ! NOTE we use the d_seg that is more convenient here than d
            ! 				g_ij_dtt_(i,i) = ( 2.0_dp*dzeta_dtt_r(3) &
            ! 							&	 +3.0_dp*d_seg_dtt(i)*dzeta_r(3) &
            ! 							&	 +3.0_dp*d_seg(i)*dzeta_dtt_r(2) )/(2.0_dp*zms**2) &
            ! 							&  +( 3.0_dp*d_seg(i)*dzeta_r(2)*dzeta_dtt_r(3) &
            ! 							&	 + d_seg(i)*d_seg_dtt(i)*dzeta_r(2)**2 &
            ! 							&	 + dzeta_r(2)*dzeta_dtt_r(2)*d_seg(i)**2 )/zms**3 &
            ! 							&  +( 3.0_dp*dzeta_dtt_r(3)*(d_seg(i)*dzeta_r(2))**2 )/(2.0_dp*zms**4)
            !			end if
            !
        end do
    end do

    call write_array_2d (g_ij_,      "g_ij_      : ",stddeb)
    call write_array_2d (g_ij_drho_, "g_ij_drho  : ",stddeb)
    call write_array_2d (g_ij_drho2_,"g_ij_drho2 : ",stddeb)
    call write_array_2d (g_ij_drho3_,"g_ij_drho3 : ",stddeb)

    call write_array_2d (g_ij_dtt_,  "g_ij_dtt   : ",stddeb)

    call write_array_3d (g_ij_dx_,   "g_ij_dx    : ",stddeb)
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    ! output the results into control_mod for easy access
    g_ij = g_ij_
    g_ij_drho = g_ij_drho_
    g_ij_drho2 = g_ij_drho2_
    g_ij_drho3 = g_ij_drho3_
    ! 	g_ij_drho4 = g_ij_drho4_
    ! 	g_ij_drho5 = g_ij_drho5_
    if (composition_der .eqv. .true.) then
        g_ij_dx = g_ij_dx_
    end if
    if (temperature_der .eqv. .true.) then
        g_ij_dtt = g_ij_dtt_
    end if

    init_flag_gij = .true.

    ! ----- deallocation section -----
    deallocate(g_ij_)
    deallocate(g_ij_drho_)
    deallocate(g_ij_drho2_)
    deallocate(g_ij_drho3_)
    !	deallocate( g_ij_drho4_)
    !	deallocate( g_ij_drho5_)

    deallocate(g_ij_dx_)
    deallocate(g_ij_dtt_)

    deallocate(dzeta_dx_r)
    write (stddeb,*) "<=== rdf_sub" !DEBUG
end subroutine rdf_sub
