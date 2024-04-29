subroutine a_hc_sub (rho_in)
    ! subroutine for the reduced helmholtz energy contribution from hard chain HC term calculation
    ! please note that the HC term already incorporates the hard sphere HS term
    ! - preapre the necesarry components such as dzeta_fun, rdf_sub (because of g_ij) and a_hs_sub
    !	* calling only when necesarry
    ! - the actual initialization
    ! - computation of the associated sum needed for hard chain Helmholtz energy contribution
    !	* the loop also compute composition or temperature derivative if it is required
    !	* calling only when necesarry
    ! - parameter update in control module
    !
    !	== LAST MODIFICATIONS ==
    !	15.5.2017 - D.C - descriptions and comments
    !	20.6.2017 - D.C - init_flag_g dependency omitment(partial) -> moved to parent control module -> replaced by single flag
    !	04.7.2017 - D.C - moving the z_hs contribution to get method in control mod

    use control_mod, only : dp, n_comp, mol_rat, &
            &    stdout, stderr, stddeb, stdlog, &
            &    composition_der, temperature_der, &
            &    g_ij, g_ij_dx, g_ij_dtt, &
            &    a_hc, a_hc_dx, a_hc_dtt, &
            &    init_flag_gij, init_flag_ahc
    use param_list_mod, only : param_m
    implicit none
    ! variable section IN/OUT
    real(dp), intent(in) :: rho_in ! the input density
    ! local variables outputed into control_mod
    real(dp) :: a_hc_ ! the hard chain compressibility contribution
    real(dp), dimension(:), allocatable :: a_hc_dx_ ! the first composition derivative of hard chain compressibility contribution, size [n_comp]
    real(dp) :: a_hc_dtt_ ! the first temperature derivative of hard chain compressibility contribution
    ! local variables
    integer :: i, j ! counters

    !------- PREPARING ------
    write (stddeb,*) "==> a_hc_sub" !DEBUG
    ! call for the required level 1 subroutine rdf
    if (init_flag_gij .eqv. .false.) then
        call rdf_sub(rho_in)
    end if

    ! ----- allocation section -----
    allocate(a_hc_dx_(n_comp))

    ! initialization without hard sphere helmholtz energy - in oposite to (A.4) the operation is left for get function
    write (stddeb,*) '-----  HC  -----'  ! DEBUG just for brevity of debug

    a_hc_ = 0.0_dp
    a_hc_dtt_ = 0.0_dp
    a_hc_dx_ = 0.0_dp

    do i = 1, n_comp
        ! 	write (*,*), '--------- A-HC ---------'  ! DEBUG just for brevity of debug
        ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.4)
        a_hc_ = a_hc_ + mol_rat(i) * (1.0_dp - param_m(i)) * log(g_ij(i, i))

        ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.35)
        if (composition_der .eqv. .true.) then
            do j = 1, n_comp
                a_hc_dx_(i) = a_hc_dx_(i) + mol_rat(j) * (1.0_dp - param_m(j)) * g_ij_dx(j, j, i) / g_ij(j, j)
            end do
            ! correction for error in paper Gross-Sadowski 2001 application of PC-SAFT eq (A.35)
            ! the shape of derivative is identical with used fact that dzeta_dx(0)=0
            a_hc_dx_(i) = a_hc_dx_(i) + (1.0_dp - param_m(i)) * log(g_ij(i, i))
        end if

        ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.54)

        if (temperature_der .eqv. .true.) then
            a_hc_dtt_ = a_hc_dtt_ + mol_rat(i) * (1.0_dp - param_m(i)) * g_ij_dtt(i, i) / g_ij(i, i)
        end if
    end do
    ! output the results into control_mod for easy access

    write (stddeb,*) "a_hc_    : ",a_hc_ ! DC DEBUG
    if (composition_der .eqv. .true.) then
        write (stddeb,*) "a_hc_dx_ : ",a_hc_dx_ ! DC DEBUG
    end if
    if (temperature_der .eqv. .true.) then
        write (stddeb,*) "a_hc_dtt_: ",a_hc_dtt_ ! DC DEBUG
    end if
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    a_hc = a_hc_
    if (composition_der .eqv. .true.) then
        a_hc_dx = a_hc_dx_
    end if
    if (temperature_der .eqv. .true.) then
        a_hc_dtt = a_hc_dtt_
    end if

    init_flag_ahc = .true.

    ! ----- deallocation section -----
    deallocate(a_hc_dx_)

    write (stddeb,*) "<== a_hc_sub" !DEBUG
end subroutine a_hc_sub
