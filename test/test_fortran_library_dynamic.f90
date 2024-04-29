
! minimal testing example of external library usage
program test_prog
    use iso_c_binding
    use bind_interface_mod
    use control_mod, only: stdout, stderr, stddeb, stdlog
    use param_list_mod, only: param_m
    use interface_mod, only: deallocate_all_sub
    implicit none

    integer, parameter :: name_length = 32
    real(dp) :: rho_in ! input density
    real(dp) :: p_out, p_drho_out, p_drho2_out ! pressure calculation results
    real(dp), dimension(:), allocatable :: mu, fug ! chem pot calculation results
    real(dp) :: h_out, s_out, g_out ! enthalpy entropy calculation results
    character (len=name_length) :: file_name ! name of the file that should be loaded
    character (kind=c_char, len=1), dimension (name_length) :: file_name_input ! name of the file that should be loaded
    integer :: i

    rho_in = 997.0
    file_name = "default.txt"
    do i=1, name_length - 1
        file_name_input(i) = file_name(i:i)
        ! print *, file_name_input(i) !DEBUG
    end do
    
    file_name_input(name_length) = c_null_char

    call bind_parse_input(file_name_input, len(file_name_input)) ! parse from default.txt file
!
!    call bind_press_calc_sub(p_out, p_drho_out, p_drho2_out)
!
!    allocate(mu(n_comp), fug(n_comp))
!    call bind_chem_pot_calc_sub(mu, fug)
!    deallocate(mu, fug)
!
!    call bind_enthalpy_entropy_calc_sub(h_out, s_out, g_out)

    allocate(mu(n_comp), fug(n_comp))
    call bind_all_calc(p_out, p_drho_out, p_drho2_out, &
                     & mu, fug, &
                     & h_out, s_out, g_out)

    write(*,"('!===== RESULTS =====')")
    write(*,"('!beware the properties are residual, i.e. no ideal gas contribution is added ')")
    write(*,*) "pressure          : ", p_out, "[Pa]"
    write(*,*) "dP/drho           : ", p_drho_out, "[Pa*m^3/mol]"
    write(*,*) "dP2/drho2         : ", p_drho2_out, "[Pa*m^6/mol^2]"
    write(*,*) "chemical potential: ", mu
    write(*,*) "fugacity coeffs   : ",fug
    write(*,*) "enthaply          :", h_out, "[J]"
    write(*,*) "entropy           :", s_out, "[J]"
    write(*,*) "gibbs energy      :", g_out, "[J]"

    deallocate(mu, fug)

    call deallocate_all_sub()
    
end program test_prog