subroutine num_diff_calc_sub(rho_in, var_type, diff_step, num_diff)
    ! The evaluation of compresibility factor and helmholtz energy numerical difference
    ! var_type states what is varying: 1=density, 2=temperature, 3=composition
    ! employed difference shcema: CENTRAL
    ! purpose of this function is mainly for tests
    ! - the local variables definition
    ! - select case for operatin the differences
    !	* density difference has to set the new density via the set_density_sub()
    !	* others can use the initialize_***_sub for same operation
    !	* differences are computed directly without the need of additional function
    !
    !	== LAST MODIFICATIONS ==
    !	15.5.2017 - D.C - descriptions and comments
    !	04.7.2017 - D.C - utilising the get function to simplify numerical differences
    !	??.7.2017 - D.C - polarity implementation edits to output
    !	15.8.2017 - D.C - composition derivative output correction for 1,n-component systems

    use control_mod, only : n_comp, dp, n_a, k_b, m2angstrom, &
            &  tt, mol_rat, &
            &  z_hs, z_hc, z_disp, z_dd, z_qq, z_dq, z_asoc, &
            &  z_hs_drho, z_hc_drho, z_disp_drho, z_dd_drho, z_qq_drho, z_dq_drho, z_asoc_drho, &
            &  a_hs, a_hc, a_disp, a_dd, a_qq, a_dq, a_asoc, &
            &  a_hs_dtt, a_hc_dtt, a_disp_dtt, a_dd_dtt, a_qq_dtt, a_dq_dtt, a_asoc_dtt, &
            &  a_hs_dx, a_hc_dx, a_disp_dx, a_dd_dx, a_qq_dx, a_dq_dx, a_asoc_dx, &
            &  set_density_sub, set_temperature_sub, set_composition_sub, &
            &  get_zres_fun, get_zres_drho_fun, &
            &  get_ares_fun, get_ares_dtt_fun, get_ares_dx_fun
    !  	use contrib_mod ,only: m_mean
    implicit none
    ! variable section IN/OUT
    real(dp), intent(in) :: rho_in ! input density required always
    real(dp), intent(in) :: diff_step ! the difference step
    integer, intent(in) :: var_type    ! deremine the type of computed difference: density (=1), temperature (=2), composition (=3)
    real(dp), intent(out), dimension(14) :: num_diff ! output numerical difference
    ! local variables
    ! 	integer :: i ! counter
    real(dp), dimension(2, 14) :: diff_store ! storage of individual (z_hs, z_hc, z_disp, a_hs, a_hc, a_disp)values used for difference computation
    real(dp), dimension(3) :: z_res_store, a_res_store ! storage for residual properties
    real(dp), dimension(:), allocatable :: tmp_store ! temporary storage for work with the composition derivatives
    real(dp) :: tt_ ! the local copy of the global variable to hold the changes in global variable
    real(dp), dimension(:), allocatable :: mol_rat_, mol_rat_tmp ! the local copy of the global variable to hold the changes in global variable
    integer :: comp ! selected component index for numerical derivative of this component
    integer :: i ! counter

    ! ----- initialise comp section -----
    comp = 1
    ! ----- allocation section -----
    allocate(tmp_store(n_comp))
    allocate(mol_rat_(n_comp), mol_rat_tmp(n_comp))

    select case (var_type)
    case(1) ! density difference
        write (*, *) '=== Density difference section ==='
        ! set the new density and compute the required variables z_hs z_hc z_disp
        ! === values of x * (1-eps) ===
        call set_density_sub()
        z_res_store(1) = get_zres_fun(rho_in * (1.0_dp - diff_step))
        a_res_store(1) = get_ares_fun(rho_in * (1.0_dp - diff_step))
        diff_store(1, :) = (/z_hs, z_hc, z_disp, z_dd, z_qq, z_dq, z_asoc, a_hs, a_hc, a_disp, a_dd, a_qq, a_dq, a_asoc /)

        ! === values of x * (1+eps) ===
        call set_density_sub() ! computation reset of density dependent variables
        z_res_store(2) = get_zres_fun(rho_in * (1.0_dp + diff_step))
        a_res_store(2) = get_ares_fun(rho_in * (1.0_dp + diff_step))
        diff_store(2, :) = (/z_hs, z_hc, z_disp, z_dd, z_qq, z_dq, z_asoc, a_hs, a_hc, a_disp, a_dd, a_qq, a_dq, a_asoc /)

        ! === value of x ===
        call set_density_sub()
        z_res_store(3) = get_zres_drho_fun(rho_in)
        a_res_store(3) = get_ares_fun(rho_in)

        num_diff = (diff_store(2, :) - diff_store(1, :)) / (2.0_dp * diff_step * rho_in)

        ! === verification ===
        write (*, *) '=== plain output of values ==='
        ! 				write (*,"('z_hs:  ',e16.8,' | z_hc:  ',e16.8,' | z_disp:  ',e16.8,' | z_dd:  ',e16.8,' | z_qq:  ',e16.8,' | z_dq:  ',e16.8,' | z_asoc:  ',e16.8)") &
        ! 						&   z_hs,	 z_hc,	 z_disp,	z_dd,	z_qq,	z_dq,	z_asoc
        ! 				write (*,"('a_hs:',e16.8,' | a_hc:',e16.8,' | a_disp:',e16.8,' | a_dd:  ',e16.8,' | a_qq:  ',e16.8,' | a_dq:  ',e16.8,' | a_asoc:  ',e16.8)") &
        ! 						&   a_hs,	a_hc,	  a_disp,	a_dd, a_qq,	a_dq,	a_asoc
        write (*, "('z_hs_drho:',e16.8,' | z_hc_drho:',e16.8,' | z_disp_drho:',e16.8,' | z_dd_drho:  ',e16.8,' | z_qq_drho:  ',e16.8,' | z_dq_drho:  ',e16.8,' | z_asoc_drho:  ',e16.8)") &
                &   z_hs_drho, z_hc_drho, z_disp_drho, z_dd_drho, z_qq_drho, z_qq_drho, z_asoc_drho

        write (*, *) '=== analytical | numerical intermediate comparison ==='
        write (*, "('z_hs_drho  :',/e16.8,/e16.8)") z_hs_drho, num_diff(1)
        write (*, "('z_hc_drho  :',/e16.8,/e16.8)") z_hc_drho, num_diff(2)
        write (*, "('z_disp_drho:',/e16.8,/e16.8)") z_disp_drho, num_diff(3)
        write (*, "('z_dd_drho  :',/e16.8,/e16.8)") z_dd_drho, num_diff(4)
        write (*, "('z_qq_drho  :',/e16.8,/e16.8)") z_qq_drho, num_diff(5)
        write (*, "('z_dq_drho  :',/e16.8,/e16.8)") z_dq_drho, num_diff(6)
        write (*, "('z_asoc_drho:',/e16.8,/e16.8)") z_asoc_drho, num_diff(7)

        write (*, *) '=== analytical | numerical result comparison ==='
        write (*, "('z_res_drho  :',/e16.8,/e16.8)") z_res_store(3), (z_res_store(2) - z_res_store(1)) / (2.0_dp * diff_step * rho_in)

        ! 				write (*,"('=== analytical | numerical final result comparison ===)")
        ! 				write (*,"('p_drho  :',e16.8,' | ',e)"), z_res_store(3), ( z_res_store(2) - z_res_store(1) )/(2.0_dp*diff_step*rho_in)

    case(2) ! temperature difference
        tt_ = tt
        write (*, *) '=== Temperature difference section ==='
        ! set the new temeprature and compute the required variables z_hs z_hc z_disp z_dd
        ! === values of tt * (1-eps) ===
        call set_temperature_sub(tt_ * (1.0_dp - diff_step))
        z_res_store(1) = get_zres_fun(rho_in)
        a_res_store(1) = get_ares_fun(rho_in)
        diff_store(1, :) = (/z_hs, z_hc, z_disp, z_dd, z_qq, z_dq, z_asoc, a_hs, a_hc, a_disp, a_dd, a_qq, a_dq, a_asoc /)

        ! === values of tt * (1+eps) ===
        call set_temperature_sub(tt_ * (1.0_dp + diff_step))
        z_res_store(2) = get_zres_fun(rho_in)
        a_res_store(2) = get_ares_fun(rho_in)
        diff_store(2, :) = (/z_hs, z_hc, z_disp, z_dd, z_qq, z_dq, z_asoc, a_hs, a_hc, a_disp, a_dd, a_qq, a_dq, a_asoc /)

        ! === value of tt ===
        call set_temperature_sub(tt_)
        z_res_store(3) = get_zres_fun(rho_in)
        a_res_store(3) = get_ares_dtt_fun(rho_in)

        num_diff = (diff_store(2, :) - diff_store(1, :)) / (2.0_dp * diff_step * tt_)

        ! === verification ===
        write (*, *) '=== plain output of values ==='
        ! 				write (*,"('z_hs:  ',e16.8,' | z_hc:  ',e16.8,' | z_disp:  ',e16.8,' | z_dd:  ',e16.8,' | z_qq:  ',e16.8,' | z_dq:  ',e16.8,' | z_asoc:  ',e16.8)") &
        ! 						&   z_hs,	 z_hc,	 z_disp,	z_dd,	z_qq,	z_dq,	z_asoc
        ! 				write (*,"('a_hs:',e16.8,' | a_hc:',e16.8,' | a_disp:',e16.8,' | a_dd:  ',e16.8,' | a_qq:  ',e16.8,' | a_dq:  ',e16.8,' | a_asoc:  ',e16.8)") &
        ! 						&   a_hs,	a_hc,	  a_disp,	a_dd, a_qq,	a_dq,	a_asoc
        write (*, "('a_hs_dtt:',e16.8,' | a_hc_dtt:',e16.8,' | a_disp_dtt:',e16.8,' | a_dd_dtt:  ',e16.8,' | a_qq_dtt:  ',e16.8,' | a_dq_dtt:  ',e16.8,' | a_asoc_dtt:  ',e16.8)") &
                &   a_hs_dtt, a_hc_dtt, a_disp_dtt, a_dd_dtt, a_qq_dtt, a_dq_dtt, a_asoc_dtt

        write (*, *) '=== analytical | numerical intermediate comparison ==='
        write (*, "('a_hs_dtt  :',/e16.8,/e16.8)") a_hs_dtt, num_diff(8)
        write (*, "('a_hc_dtt  :',/e16.8,/e16.8)") a_hc_dtt, num_diff(9)
        write (*, "('a_disp_dtt:',/e16.8,/e16.8)") a_disp_dtt, num_diff(10)
        write (*, "('a_dd_dtt  :',/e16.8,/e16.8)") a_dd_dtt, num_diff(11)
        write (*, "('a_qq_dtt  :',/e16.8,/e16.8)") a_qq_dtt, num_diff(12)
        write (*, "('a_dq_dtt  :',/e16.8,/e16.8)") a_dq_dtt, num_diff(13)
        write (*, "('a_asoc_dtt:',/e16.8,/e16.8)") a_asoc_dtt, num_diff(14)

        write (*, *) '=== analytical | numerical result comparison ==='
        write (*, "('a_res_dtt  :',/e16.8,/e16.8)") a_res_store(3), (a_res_store(2) - a_res_store(1)) / (2.0_dp * diff_step * tt_)

    case(3)
        if (comp > n_comp) then
            write(*, *) "WARNING: Number of component derivative is higher than actual components: set it to 1."
            comp = 1
        end if

        mol_rat_ = mol_rat
        mol_rat_tmp = mol_rat
        write (*, "('=== Component (',i2,') composition difference section ===')") comp
        ! set the new composition and compute the required variables z_hs z_hc z_disp z_dd
        ! === values of mol_rat(comp) * (1-eps) ===
        mol_rat_tmp(comp) = mol_rat_(comp) * (1.0_dp - diff_step)
        ! 				mol_rat_tmp(modulo(comp,n_comp)+1) = 1.0_dp-mol_rat_(comp)*(1.0_dp-diff_step)
        call set_composition_sub(mol_rat_tmp)
        z_res_store(1) = get_zres_fun(rho_in)
        a_res_store(1) = get_ares_fun(rho_in)
        diff_store(1, :) = (/z_hs, z_hc, z_disp, z_dd, z_qq, z_dq, z_asoc, a_hs, a_hc, a_disp, a_dd, a_qq, a_dq, a_asoc /)

        ! === values of mol_rat(1) * (1+eps) ===
        mol_rat_tmp(comp) = mol_rat_(comp) * (1.0_dp + diff_step)
        ! 				mol_rat_tmp(modulo(comp,n_comp)+1) = 1.0_dp-mol_rat_(comp)*(1.0_dp+diff_step)
        call set_composition_sub(mol_rat_tmp)
        z_res_store(2) = get_zres_fun(rho_in)
        a_res_store(2) = get_ares_fun(rho_in)
        diff_store(2, :) = (/z_hs, z_hc, z_disp, z_dd, z_qq, z_dq, z_asoc, a_hs, a_hc, a_disp, a_dd, a_qq, a_dq, a_asoc /)

        ! === value of mol_rat ===
        call set_composition_sub(mol_rat_)
        z_res_store(3) = get_zres_fun(rho_in)
        call get_ares_dx_fun(rho_in, tmp_store)
        a_res_store(3) = tmp_store(comp)

        num_diff = (diff_store(2, :) - diff_store(1, :)) / (2.0_dp * diff_step * mol_rat_(comp))

        ! === verification ===
        write (*, *) '=== plain output of values ==='
        ! 				write (*,"('z_hs:  ',e16.8,' | z_hc:  ',e16.8,' | z_disp:  ',e16.8,' | z_dd:  ',e16.8,' | z_qq:  ',e16.8,' | z_dq:  ',e16.8,' | z_asoc:  ',e)"), &
        ! 						&   z_hs,	 z_hc,	 z_disp,	z_dd,	z_qq,	z_dq,	z_asoc
        ! 				write (*,"('a_hs:',e16.8,' | a_hc:',e16.8,' | a_disp:',e16.8,' | a_dd:  ',e16.8,' | a_qq:  ',e16.8,' | a_dq:  ',e16.8,' | a_asoc:  ',e)"), &
        ! 						&   a_hs,	a_hc,	  a_disp,	a_dd, a_qq,	a_dq,	a_asoc
        do i = 1, n_comp
            write (*, "('a_hs_dx(',i2,'):',e16.8,' | a_hc_dx(',i2,'):',e16.8,' | a_disp_dx(',i2,'):',e16.8,' | a_dd_dx(',i2,'):  ',e16.8,' | a_qq_dx(',i2,'):  ',e16.8,' | a_dq_dx(',i2,'):  ',e16.8,' | a_asoc_dx(',i2,'):  ',e16.8)") &
                    &   i, a_hs_dx(i), i, a_hc_dx(i), i, a_disp_dx(i), i, a_dd_dx(i), i, a_qq_dx(i), i, a_dq_dx(i), i, a_asoc_dx(i)
        end do

        write (*, *) '=== analytical | numerical intermediate comparison ==='
        write (*, "('a_hs_dx(',i2,')  :',/e16.8,/e16.8)") comp, a_hs_dx(comp), num_diff(8)
        write (*, "('a_hc_dx(',i2,')  :',/e16.8,/e16.8)") comp, a_hc_dx(comp), num_diff(9)
        write (*, "('a_disp_dx(',i2,'):',/e16.8,/e16.8)") comp, a_disp_dx(comp), num_diff(10)
        write (*, "('a_dd_dx(',i2,')  :',/e16.8,/e16.8)") comp, a_dd_dx(comp), num_diff(11)
        write (*, "('a_qq_dx(',i2,')  :',/e16.8,/e16.8)") comp, a_qq_dx(comp), num_diff(12)
        write (*, "('a_dq_dx(',i2,')  :',/e16.8,/e16.8)") comp, a_dq_dx(comp), num_diff(13)
        write (*, "('a_asoc_dx(',i2,'):',/e16.8,/e16.8)") comp, a_asoc_dx(comp), num_diff(14)

        write (*, *) '=== analytical | numerical result comparison ==='
        write (*, "('a_res_dx(',i2,')  :',/e16.8,/e16.8)") comp, a_res_store(3), (a_res_store(2) - a_res_store(1)) / (2.0_dp * diff_step * mol_rat_(comp))

    case default
        write (*, "('wrong var type: ',i2,' selected: nothing is done')") var_type
    end select

    ! ----- deallocation section -----
    deallocate(tmp_store)
    deallocate(mol_rat_, mol_rat_tmp)

end subroutine num_diff_calc_sub
