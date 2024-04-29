module bind_interface_mod
    ! binding of the importatnt interface function to C
    ! intended as only place exposing the iso_C_bound directives


    !	== LAST MODIFICATIONS ==
    !	28.01.2020 - D.C - implementation, descriptions and comments
    !   08.12.2021 - D.C - clone for win compilation

    use iso_c_binding
    use control_mod, only: rho, dp, n_comp, ch40_length, ch512_length, &
                         & n_a, m2angstrom, &
                         & set_density_sub, set_temperature_sub, set_composition_sub
    use interface_mod, only: &
                           & parse_input_universal, setup_flags_parse_sub, set_fluid!, &
!                           & allocate_all_sub, deallocate_all_sub
    implicit none
    real(dp) :: local_dens_packing ! the bind only relevant density to retain ssame form for setting input variables
contains
    ! ===== set methods =====
    subroutine bind_get_constants(r_gas_out, n_a_out, k_b_out, pi_out)
        ! return the constatnts the equation uses for calculation
        ! TODO consider allowing setting the constants
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"get_constants" :: bind_get_constants
        use control_mod, only: r_gas, n_a, k_b, pi

        implicit none
        real(C_DOUBLE), intent(out) :: r_gas_out ! output r_gas
        real(C_DOUBLE), intent(out) :: n_a_out ! output n_a
        real(C_DOUBLE), intent(out) :: k_b_out ! output k_b
        real(C_DOUBLE), intent(out) :: pi_out ! output pi

        r_gas_out = r_gas
        n_a_out = n_a
        k_b_out = k_b
        pi_out = pi

    end subroutine bind_get_constants

    subroutine bind_deallocate_memory()
        ! handle the file or command line input
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"deallocate_memory" :: bind_deallocate_memory
        use interface_mod, only: deallocate_all_sub
        implicit none
               
        call deallocate_all_sub()        

    end subroutine bind_deallocate_memory

    subroutine bind_parse_input(file_name_in, str_length)
        ! handle the file or command line input
        ! this is ISO C bound function
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"parse_input" :: bind_parse_input
        !DEC$ ATTRIBUTES VALUE :: str_length
       implicit none
        integer(c_int), intent(in), value :: str_length
        character (kind=c_char, len=1), dimension (*), intent (in) :: file_name_in
        character(len = ch512_length) :: file_name
        integer :: i, j

        if (str_length .gt. 0) then
            if( str_length .gt. ch512_length) then
                write (*, "('! WARNING input character lenght exceeds allowed -> may be unable to parse the path: ', i4)") ch512_length
            end if
        end if
        
        do i=1, ch512_length
            if ( file_name_in(i) == c_null_char ) then
                do j=i, ch512_length
                    file_name(j:j) = " "                    
                end do 
                exit
            else
                file_name(i:i) = file_name_in(i)
            end if
        end do
        ! print *, file_name
        call parse_input_universal(file_name, .true.)
        local_dens_packing = rho

    end subroutine bind_parse_input

   subroutine bind_set_n_comp_sub(n_comp_in)
       ! handle the file or command line input
       ! this is ISO C bound function
       !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"set_n_comp" :: bind_set_n_comp_sub
       !DEC$ ATTRIBUTES VALUE :: n_comp_in
       use param_list_mod, only: param_m
       use interface_mod, only: allocate_all_sub, deallocate_all_sub
       
       implicit none
       integer(c_int), intent(in), value ::n_comp_in

       if ((n_comp_in .le. 0) .or. (n_comp_in .gt. 99)) then
           write (*, "('! ERROR invalid number of components: ',i2,' is not in range <1,99>')")
           stop
       end if

       ! write (*,*) size(param_m) ! DEBUG print the original allocated size
       ! === de/allocate section ===
       if (n_comp .ne. n_comp_in) then
           call deallocate_all_sub()
           ! write (*,*) "DEBUG not same size "
       end if

       n_comp = n_comp_in

       ! write (*,*) size(param_m) ! DEBUG test deallocation

       call allocate_all_sub()

       ! write (*, "('! setting new number of components: |',i2,'|')") n_comp !DEBUG just for verification and user review of input

       ! write (*,*) size(param_m) ! DEBUG print the newly allocated size

   end subroutine bind_set_n_comp_sub

    subroutine bind_set_component_name_sub(fluid_name_in, str_length)
        ! handle the component setting
        ! this is ISO C bound function
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"set_component_name" :: bind_set_component_name_sub
        !DEC$ ATTRIBUTES VALUE :: str_length
        implicit none
        integer(c_int), intent(in), value :: str_length
!        character(kind=c_char), intent(in) :: fluid_name_in(str_length)
        character (kind=c_char, len=1), dimension (*), intent (in) :: fluid_name_in
        ! local variables
        character(len = :), allocatable :: fluid_name
        integer :: i, str_size

        str_size = str_length
        allocate(character(str_size) :: fluid_name)

        do i=1, str_size
            if ( fluid_name_in(i) == c_null_char ) then
                exit
            else
                fluid_name(i:i) = fluid_name_in(i)
            end if
        end do

        ! write (*, "('! parsing new fluid names: <',a,'>')") fluid_name

        call set_fluid(fluid_name)
        deallocate(fluid_name)

    end subroutine bind_set_component_name_sub

    subroutine bind_set_param_sub(sub_ind, par_array)
        ! procedure for direct setting of parameters: m | s | e/k | kapAB | eAB/k | u | n_u | q | n_q | M | bond
        ! operates with parameters for single substance
        ! accepts the substance index and parameter vector in prescribed order
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"set_component_param" :: bind_set_param_sub
        !DEC$ ATTRIBUTES VALUE :: sub_ind
        use param_list_mod, only: set_param_sub

        implicit none
        integer, parameter :: n_of_param = 11
        integer(c_int), intent(in), value :: sub_ind
        real(C_DOUBLE), dimension(n_of_param), intent(in) :: par_array
        real(dp), dimension(n_of_param) :: par_array_in
        integer :: i

        do i=1, n_of_param
            par_array_in(i) = par_array(i)
        end do  
        ! write(*,*) par_array_in

        call set_param_sub(sub_ind, par_array_in)
        ! write(*,*) "after call" ! DEBUG

    end subroutine bind_set_param_sub

    subroutine bind_set_density_sub(density_in)
        ! handle the file or command line input
        ! this is ISO C bound function
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"set_density" :: bind_set_density_sub
        !DEC$ ATTRIBUTES VALUE :: density_in
        implicit none
        real(C_DOUBLE), intent(in), value :: density_in ! input density

        if (density_in .le. 0_dp) then
            write (*, "('! ERROR invalid molar density(<=0): ',e16.8,' [mol/m^3]')") density_in
            stop
        end if
        ! write (*, "('! setting new molar density: ',e16.8,' [mol/m^3]')") density_in
        local_dens_packing = density_in * n_a / m2angstrom**3 ! conversion to the computation dimension
        
        call set_density_sub()

    end subroutine bind_set_density_sub

    subroutine bind_set_pacfrac_sub(pacfrac_in, density_out)
        ! handle the setting of packing fraction and return density
        ! this is ISO C bound function
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"set_pfrac" :: bind_set_pacfrac_sub
        !DEC$ ATTRIBUTES VALUE :: pacfrac_in
        use contrib_mod, only : dzeta, initialize_dzeta_fun, &
                            &   pi
        implicit none

        real(C_DOUBLE), intent(in), value :: pacfrac_in ! input density
        real(C_DOUBLE), intent(out) :: density_out ! input density

        if (pacfrac_in .le. 0.0_dp) then ! unfeasible low packing fraction
            write (*, "('! ERROR invalid packing fraction (<=0): ',e16.8,' [#]')") pacfrac_in
            density_out = -1.0
            return
        end if
        if (pacfrac_in .gt. pi/(3_dp*sqrt(2.0_dp)) ) then ! unfeasible high, more than closest packing fraction
            write (*, "('! ERROR packing fraction higher than closest packing (>pi/(3*sqrt(2))): ',e16.8,' [#]')") pacfrac_in
            density_out = -2.0
            return
        end if
    
        if (initialize_dzeta_fun().eqv..false.) then
            write (*, *) 'Error in helper function initialize_dzeta_fun'
            density_out = -3.0
            return
        end if
        local_dens_packing = pacfrac_in / dzeta(3);
        density_out = local_dens_packing / n_a * m2angstrom**3
        
        call set_density_sub()

    end subroutine bind_set_pacfrac_sub


    subroutine bind_set_temperature_sub(temperature_in)
        ! handle the file or command line input
        ! this is ISO C bound function
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"set_temperature" :: bind_set_temperature_sub
        !DEC$ ATTRIBUTES VALUE :: temperature_in
        implicit none
        real(C_DOUBLE), intent(in), value :: temperature_in ! input density

        if (temperature_in .le. 0_dp) then
            write (*, "('! ERROR invalid temperature(<=0): ',e16.8,' [K]')") temperature_in
            stop
        end if
        ! write (*, "('! setting new temperature: ',e16.8,' [K]')") temperature_in
        call set_temperature_sub(temperature_in)

    end subroutine bind_set_temperature_sub

    subroutine bind_set_molar_ratio_sub(mol_rat_in)
        ! handle the file or command line input
        ! this is ISO C bound function
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"set_molar_ratio" :: bind_set_molar_ratio_sub
        implicit none
        real(C_DOUBLE), dimension(n_comp), intent(in) :: mol_rat_in
        real(C_DOUBLE), dimension(n_comp) :: mol_rat_tmp
        integer :: i

        if (size(mol_rat_in) .ne. n_comp) then
            write (*, "('! ERROR invalid molar_ratio size: ',i4,' [mol/m^3]')") size(mol_rat_in)
            stop
        end if
        do i = 1,n_comp
            if (mol_rat_in(i) .le. 0_dp) then
                write (*, "('! ERROR invalid molar_ratio: ',e16.8,' [mol/m^3]')") mol_rat_in(i)
                stop
            end if
            mol_rat_tmp(i) = mol_rat_in(i)
        end do

        ! write (*, "('! setting new molar_ratio: ')")
        ! write (*,"(f16.8)") mol_rat_tmp
        call set_composition_sub(mol_rat_in)

    end subroutine bind_set_molar_ratio_sub

    subroutine bind_press_calc_sub(p_out, p_drho_out, p_drho2_out)
        ! handle the pressure calculation
        ! this is ISO C bound function
        ! required wrapping the core fortran functionality of calculating only the necessary
        ! set the flags so that all is calculated
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"press_calc" :: bind_press_calc_sub
        use control_mod, only: rho,tt

        implicit none
!        real(C_DOUBLE), intent(in), value :: rho_in ! input density
        real(C_DOUBLE), intent(out) :: p_out ! output pressure
        real(C_DOUBLE), intent(out) :: p_drho_out ! output pressure density derivative
        real(C_DOUBLE), intent(out) :: p_drho2_out ! output pressure second density derivative
        real(dp) :: p, p1, p2 ! output values of computational procedures
        real(dp) :: tmp ! DEBUG

        call setup_flags_parse_sub('p,pp,ppp') ! set flags for all pressure calculations

        ! error output initialization
        p_out = -1.0_dp
        p_drho_out = -1.0_dp
        p_drho2_out = -1.0_dp

        call press_calc_sub(local_dens_packing, p, p1, p2)

        ! DEBUG section
        ! write(*,*) "F> Temperature, packing density"
        ! write(*,*) tt,local_dens_packing
        ! write(*,*) "F> Pressure, dP/drho, dP^2/drho^2"
        ! write(*,*) p,p1,p2
        
        p_out = p
        p_drho_out = p1
        p_drho2_out = p2

    end subroutine bind_press_calc_sub

    subroutine bind_chem_pot_calc_sub(mu_out, fugacity_out)
        ! handle the pressure calculation
        ! this is ISO C bound function
        ! required wrapping the core fortran functionality of calculating only the necessary
        ! set the flags so that all is calculated
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"chem_pot_calc" :: bind_chem_pot_calc_sub
        use control_mod, only: rho,tt
        implicit none

        interface
            subroutine chem_pot_calc_sub(rho_in, chem_pot_out, fugacity_out)
                use control_mod
                real(dp), intent(in) :: rho_in ! input density
                real(dp), dimension(:), intent(out) :: chem_pot_out ! output chemical potential
                real(dp), dimension(:), intent(out) :: fugacity_out ! output fugacity coeffitient
            end subroutine chem_pot_calc_sub
        end interface

        real(C_DOUBLE), dimension(n_comp), intent(out) :: mu_out ! output pressure
        real(C_DOUBLE), dimension(n_comp), intent(out) :: fugacity_out ! output pressure density derivative

        real(dp), dimension(:),allocatable :: mu, fug ! output values of computational procedures

        call setup_flags_parse_sub('mu,fug') ! set flags for all pressure calculations

        ! write(*,*) "pack dens", local_dens_packing, "dens", rho, "temp", tt

        allocate(mu(n_comp), fug(n_comp))
        call chem_pot_calc_sub(local_dens_packing, mu, fug)

!        DEBUG section
!        write(*,*) "Chemical potential"
!        write(*,*) mu
!        write(*,*) "Fugacity"
!        write(*,*) fug

        mu_out = mu
        fugacity_out = fug

        deallocate(mu, fug)
    end subroutine bind_chem_pot_calc_sub

    subroutine bind_enthalpy_entropy_calc_sub(h_out, s_out, g_out)
        ! handle the pressure calculation
        ! this is ISO C bound function
        ! required wrapping the core fortran functionality of calculating only the necessary
        ! set the flags so that all is calculated
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"enthalpy_entropy_calc" :: bind_enthalpy_entropy_calc_sub
        implicit none
        !        real(C_DOUBLE), intent(in), value :: rho_in ! input density
        real(C_DOUBLE), intent(out) :: h_out ! output pressure
        real(C_DOUBLE), intent(out) :: s_out ! output pressure density derivative
        real(C_DOUBLE), intent(out) :: g_out ! output pressure second density derivative
        real(dp) :: h, s, g ! output values of computational procedures

        call setup_flags_parse_sub('h,s,g') ! set flags for all pressure calculations

        ! error output initialization
        h_out = -1.0_dp
        s_out = -1.0_dp
        g_out = -1.0_dp

        call enthalpy_entropy_calc_sub(local_dens_packing, h, s, g)

!        DEBUG section
!        write(*,*) "enthalphy, entropy, Gibbs energy"
!        write(*,*) h, s, g

        h_out = h
        s_out = s
        g_out = g

    end subroutine bind_enthalpy_entropy_calc_sub

    subroutine bind_helmholtz_contrib_calc_sub(a_res_out, a_hs_out, a_hc_out, a_disp_out, a_qq_out, a_dd_out, a_asoc_out)
        ! handle the pressure calculation
        ! this is ISO C bound function
        ! required wrapping the core fortran functionality of calculating only the necessary
        ! set the flags so that all is calculated
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"helmholtz_contrib_calc" :: bind_helmholtz_contrib_calc_sub
        use control_mod, only: a_hs, a_hc, a_disp, a_dd, a_qq, a_dq, a_asoc, get_ares_fun
        implicit none
        !        real(C_DOUBLE), intent(in), value :: rho_in ! input density
        real(C_DOUBLE), intent(out) :: a_res_out  ! output helmholtz res contribution
        real(C_DOUBLE), intent(out) :: a_hs_out   ! output helmholtz hs contribution
        real(C_DOUBLE), intent(out) :: a_hc_out   ! output helmholtz hc contribution
        real(C_DOUBLE), intent(out) :: a_disp_out ! output helmholtz disp contribution
        real(C_DOUBLE), intent(out) :: a_qq_out   ! output helmholtz qq contribution
        real(C_DOUBLE), intent(out) :: a_dd_out   ! output helmholtz dd contribution
        real(C_DOUBLE), intent(out) :: a_asoc_out ! output helmholtz asoc contribution

        call setup_flags_parse_sub('h') ! set flags for all pressure calculations

        ! error output initialization
        a_res_out  = -1.0_dp
        a_hs_out   = -1.0_dp
        a_hc_out   = -1.0_dp
        a_disp_out = -1.0_dp
        a_qq_out   = -1.0_dp
        a_dd_out   = -1.0_dp
        a_asoc_out = -1.0_dp

        a_res_out  = get_ares_fun(local_dens_packing)
        a_hs_out   = a_hs
        a_hc_out   = a_hc
        a_disp_out = a_disp
        a_qq_out   = a_qq
        a_dd_out   = a_dd
        a_asoc_out = a_asoc
!        DEBUG section
        ! write(*,*) "a_res : ",   a_res_out
        ! write(*,*) "a_hs  : ",    a_hs_out
        ! write(*,*) "a_hc  : ",    a_hc_out
        ! write(*,*) "a_hs+c: ",    a_hs_out + a_hc_out
        ! write(*,*) "a_disp: ",  a_disp_out
        ! write(*,*) "a_qq  : ",    a_qq_out
        ! write(*,*) "a_dd  : ",    a_dd_out
        ! write(*,*) "a_asoc: ",  a_asoc_out

    end subroutine bind_helmholtz_contrib_calc_sub

    subroutine bind_zz_contrib_calc_sub(z_res_out, z_hs_out, z_hc_out, z_disp_out, z_qq_out, z_dd_out, z_asoc_out)
        ! handle the pressure calculation
        ! this is ISO C bound function
        ! required wrapping the core fortran functionality of calculating only the necessary
        ! set the flags so that all is calculated
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"zz_contrib_calc" :: bind_zz_contrib_calc_sub
        use control_mod, only: z_hs, z_hc, z_disp, z_dd, z_qq, z_dq, z_asoc, get_zres_fun
        implicit none
        !        real(C_DOUBLE), intent(in), value :: rho_in ! input density
        real(C_DOUBLE), intent(out) :: z_res_out  ! output compresibility factor Z res contribution
        real(C_DOUBLE), intent(out) :: z_hs_out   ! output compresibility factor Z hs contribution
        real(C_DOUBLE), intent(out) :: z_hc_out   ! output compresibility factor Z hc contribution
        real(C_DOUBLE), intent(out) :: z_disp_out ! output compresibility factor Z disp contribution
        real(C_DOUBLE), intent(out) :: z_qq_out   ! output compresibility factor Z qq contribution
        real(C_DOUBLE), intent(out) :: z_dd_out   ! output compresibility factor Z dd contribution
        real(C_DOUBLE), intent(out) :: z_asoc_out ! output compresibility factor Z asoc contribution

        call setup_flags_parse_sub('h') ! set flags for all pressure calculations

        ! error output initialization
        z_res_out  = -1.0_dp
        z_hs_out   = -1.0_dp
        z_hc_out   = -1.0_dp
        z_disp_out = -1.0_dp
        z_qq_out   = -1.0_dp
        z_dd_out   = -1.0_dp
        z_asoc_out = -1.0_dp

        z_res_out  = get_zres_fun(local_dens_packing)
        z_hs_out   = z_hs
        z_hc_out   = z_hc
        z_disp_out = z_disp
        z_qq_out   = z_qq
        z_dd_out   = z_dd
        z_asoc_out = z_asoc
!        DEBUG section
        ! write(*,*) "z_res : ",   z_res_out
        ! write(*,*) "z_hs  : ",    z_hs_out
        ! write(*,*) "z_hc  : ",    z_hc_out
        ! write(*,*) "z_hs+c: ",    z_hs_out + z_hc_out
        ! write(*,*) "z_disp: ",  z_disp_out
        ! write(*,*) "z_qq  : ",    z_qq_out
        ! write(*,*) "z_dd  : ",    z_dd_out
        ! write(*,*) "z_asoc: ",  z_asoc_out

    end subroutine bind_zz_contrib_calc_sub

    subroutine bind_all_calc(p_out, p_drho_out, p_drho2_out, mu_out, fugacity_out, h_out, s_out, g_out)
        ! mimic the main function run and calculate all there is to be calculated
        ! this is ISO C bound function
        ! required wrapping the core fortran functionality of calculating only the necessary
        ! set the flags so that all is calculated
        !DEC$ ATTRIBUTES DLLEXPORT,DECORATE,ALIAS:"all_calc" :: bind_all_calc
        implicit none
        real(C_DOUBLE), intent(out) :: p_out ! output pressure
        real(C_DOUBLE), intent(out) :: p_drho_out ! output pressure density derivative
        real(C_DOUBLE), intent(out) :: p_drho2_out ! output pressure second density derivative
        real(C_DOUBLE), dimension(n_comp), intent(out) :: mu_out ! output pressure
        real(C_DOUBLE), dimension(n_comp), intent(out) :: fugacity_out ! output pressure density derivative
        real(C_DOUBLE), intent(out) :: h_out ! output pressure
        real(C_DOUBLE), intent(out) :: s_out ! output pressure density derivative
        real(C_DOUBLE), intent(out) :: g_out ! output pressure second density derivative

        call bind_press_calc_sub(p_out, p_drho_out, p_drho2_out)
        call bind_chem_pot_calc_sub(mu_out, fugacity_out)
        call bind_enthalpy_entropy_calc_sub(h_out, s_out, g_out)

    end subroutine bind_all_calc

end module bind_interface_mod
