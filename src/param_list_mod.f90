module param_list_mod
    ! module responsible for loading and maintaining the components parameters
    ! user is responsible for providing the resources in correct format
    ! - provide a simple acces to the parameters in parameter_list.txt
    !	* .txt file holds parameters important for pc-saft eos for various substances
    !	* file is stylyzed with "XML" like tags constructed for quick one-time read
    !		** file contain names first and then associated parameters (all floats) in same order
    !		** file use "!" for comments
    !
    !	== LAST MODIFICATIONS ==
    !	15.5.2017 - D.C - descriptions and comments

    use control_mod, only : dp, n_comp, ch120_length, ch40_length, n_of_file, default_res_path
    implicit none

    !private
    public param_m, param_sigma, param_eps_kb, param_kap_asoc, param_eps_asoc_kb, param_u, param_n_u, param_q, param_n_q, param_mm, param_bond, &
            & is_available, &
            & allocate_param_sub, deallocate_param_sub, initialize_param_sub, set_param_sub, get_param_sub!list the public methods

    ! variable section
    logical :: is_available ! determines the availability of parameters
    integer, dimension(:), allocatable :: pos ! positions of found components

    ! parameters form the table - for the clear denotion and filtering the param_ prefix is used
    character(len = ch120_length), dimension(:), allocatable :: param_name    ! names of components, size [n_comp]
    real(dp), dimension(:), allocatable :: param_m            ! segment number, size [n_comp]
    real(dp), dimension(:), allocatable :: param_sigma        ! sigma parameter from L-J potential, size [n_comp]
    real(dp), dimension(:), allocatable :: param_eps_kb        ! epsilon parameter from L-J potential divided by Boltzmann constant, size [n_comp]
    real(dp), dimension(:), allocatable :: param_kap_asoc    ! association volume parameter kappa AiBj, size [n_comp]
    real(dp), dimension(:), allocatable :: param_eps_asoc_kb! association volume energy epsilon AiBj divided by Boltzmann constant, size [n_comp]
    real(dp), dimension(:), allocatable :: param_u            ! dipole moment, size [n_comp]
    real(dp), dimension(:), allocatable :: param_n_u            ! fraction of dipolar segments in component, size [n_comp]
    real(dp), dimension(:), allocatable :: param_q            ! the quadrupole moment, size [n_comp]
    real(dp), dimension(:), allocatable :: param_n_q            ! fraction of quadrupolar segments in component, size [n_comp]
    real(dp), dimension(:), allocatable :: param_mm            ! molar mass, size [n_comp]
    real(dp), dimension(:), allocatable :: param_bond        ! type of bonding in associating fluids, size [n_comp]


contains
    ! ========================= ALLOCATION/ DEALLOCATION METHODS =============================
    ! ------------- ALLOCATION ---------------
    subroutine allocate_param_sub()
        ! allocate all related variables for controll module
        ! NOTE expect that n_comp is known nonzero
        if (n_comp == 0) then
            write (*, *) 'Error in function allocate_param n_comp==0'
            stop
        end if
        is_available = .false.
        allocate(pos(n_comp))

        allocate(param_name(n_comp))
        allocate(param_m(n_comp))
        allocate(param_sigma(n_comp))
        allocate(param_eps_kb(n_comp))
        allocate(param_kap_asoc(n_comp))
        allocate(param_eps_asoc_kb(n_comp))
        allocate(param_u(n_comp))
        allocate(param_n_u(n_comp))
        allocate(param_q(n_comp))
        allocate(param_n_q(n_comp))
        allocate(param_mm(n_comp))
        allocate(param_bond(n_comp))

    end subroutine allocate_param_sub
    ! ------------- DEALLOCATION -------------
    subroutine deallocate_param_sub()
        ! deallocate all related variables for controll module
        ! performs check if memory is allocated otherwise does nothing to unallocated memory

        if (allocated(pos) .eqv. .true.) then 
            deallocate(pos)
        end if

        if (allocated(param_name) .eqv. .true.) then 
            deallocate(param_name)
        end if
        if (allocated(param_m) .eqv. .true.) then 
            deallocate(param_m)
        end if
        if (allocated(param_sigma) .eqv. .true.) then 
            deallocate(param_sigma)
        end if
        if (allocated(param_eps_kb) .eqv. .true.) then 
            deallocate(param_eps_kb)
        end if
        if (allocated(param_kap_asoc) .eqv. .true.) then 
            deallocate(param_kap_asoc)
        end if
        if (allocated(param_eps_asoc_kb) .eqv. .true.) then 
            deallocate(param_eps_asoc_kb)
        end if
        if (allocated(param_u) .eqv. .true.) then 
            deallocate(param_u)
        end if
        if (allocated(param_n_u) .eqv. .true.) then 
            deallocate(param_n_u)
        end if
        if (allocated(param_q) .eqv. .true.) then 
            deallocate(param_q)
        end if
        if (allocated(param_n_q) .eqv. .true.) then 
            deallocate(param_n_q)
        end if
        if (allocated(param_mm) .eqv. .true.) then 
            deallocate(param_mm)
        end if
        if (allocated(param_bond) .eqv. .true.) then 
            deallocate(param_bond)
        end if
        
    end subroutine deallocate_param_sub
    ! =================================== SET UP METHODS ========================================
    !the subroutines and functions
    subroutine initialize_param_sub(in_names, success, is_write_to_console)
        ! search for the inputed names and initializes the parameters
        ! important notes and comments
        !use control_mod , only: n_comp ! NOTE see the hierarchy if it is not obsolete declaration
        implicit none

        ! variable section in/out
        character(len = ch40_length), dimension(:), intent(in) :: in_names ! the input names array
        integer, intent(out) :: success !  the status identificator
        logical, optional :: is_write_to_console

        logical :: write_to_console
        integer :: file_success !  the status identificator
        integer :: i, j, what  ! counters and indexers
        integer :: curline ! the current line counter
        integer :: shift !the amount of lines needed to jump to next component
        character(len = ch120_length) :: tmp_ch ! the temporary character

        if (present(is_write_to_console)) then
            write_to_console = is_write_to_console
        else
            write_to_console = .false.
        end if
        ! control section and initialization
        if (size(in_names, 1)/=n_comp) then
            write (*, *) 'Error the number of n_comp does not correspond with given lenght of inp_names'
            stop
        end if

        do i = 1, n_comp
            pos(i) = 0 ! TODO look at the intristic initialization if it is possible to use here
        end do

        ! file handling section
        open(20, file = trim(n_of_file), status = 'old', iostat = file_success)
        if (file_success /= 0) then
            if (write_to_console .eqv. .true. ) then
                write (*, "('INFO: file <',a,'> was not opened with iostat:',i4)") trim(n_of_file), file_success
                write (*, "('      trying another file: <',a,'>')") trim(default_res_path)
            end if
            close(20)
            ! second try on loading the file in res directory
            open(20, file = trim(default_res_path) // adjustl(n_of_file), status = 'old', iostat = file_success)

            if (file_success /= 0) then
                write (*, "('Error: even the other was not opened with iostat:', i4)") file_success

                success = 1
                stop
            end if
        end if

        ! search for the begining of the names section
!        write (*,*) file_success   					! DEBUG test the reading capability
        !write (*,*) '================' 				! DEBUG test the reading capability
        read(20, '(a)', iostat = success) tmp_ch
        ! 			write (*,*) success 							! DEBUG test the reading capability
        !				write (*,"('~',a,'~')") trim(tmp_ch) 	! DEBUG test the reading capability
        ! 			write (*,*), '================'
        do while ((index(tmp_ch, '<names>')==0).or.(tmp_ch(1:1)=='!'))
            read(20, '(a)', iostat = success) tmp_ch
            ! 				write (*,*) success 						! DEBUG test the reading capability
            !				write (*,"('~',a,'~')") trim(tmp_ch) 	! DEBUG test the reading capability
        end do

        ! search until the end of names section is reached
        curline = 0 ! used to number the names
        do while (index(tmp_ch, '</names>')==0)
            read(20, '(a)', iostat = success) tmp_ch
            ! write (*,*) success                     ! DEBUG test the reading capability
            ! write (*,"('~',a,'~')") trim(tmp_ch)    ! DEBUG test the reading capability
            curline = curline + 1
            if (tmp_ch(1:1)=='!') then ! to omit the commented lines
                cycle
            end if
            do i = 1, n_comp ! BUG substring is also accepted and search fails -SOLVED by change of txt format and ==
!                write (*,"('~',a,'~',a,'~')") trim(adjustl(tmp_ch)), trim(adjustl(in_names(i))) ! DEBUG test the reading capability
                if (trim(adjustl(tmp_ch))==trim(adjustl(in_names(i)))) then ! if found then note on the index
                    pos(i) = curline
                    ! DEBUG test if search is correct
                    ! write(*,"('found: ',a,' under line,coll number: ',i4,i4,' with: ',a)") trim(in_names(i)), curline, index( trim(adjustl(tmp_ch)),trim(adjustl(in_names(i))) ), trim(tmp_ch)
                end if
            end do
        end do
        ! 			write(*,"('pos=[',i4,',',i4,']')") pos(1), pos(2) ! DEBUG test the reading capability
        ! control if all found and then update the parameters
        if (minval(pos)==0) then
            write (*, *) 'Error some components were not found in the parameter list names:'
            do i = 1, n_comp
                if (pos(i) .eq. 0 ) then
                    write(*, *) in_names(i)
                end if
            end do
            stop
        else
            ! search for the begining of the parameters section
            read(20, '(a)', iostat = success) tmp_ch
            do while (index(tmp_ch, '<parameters>')==0)
                read(20, '(a)', iostat = success) tmp_ch
                ! 					write (*,"('~',a,'~')") trim(tmp_ch) 	! DEBUG test the reading capability
                if (tmp_ch(1:1)=='!') then ! to omit the commented lines
                    cycle
                end if
            end do
            ! do the search
            do i = 1, n_comp
                shift = minval(pos, mask = pos>0) ! find the minimal value
                what = minloc(pos, 1, mask = pos>0) ! locate the minimal value for correct parameter loading
                ! 					write(*,'(i4)') shift ! DEBUG test the reading capability
                do j = 1, shift - 1
                    read(20, *)
                end do
                ! 					read(20,'(a)',iostat=success) tmp_ch 	! DEBUG test the reading capability
                ! 					write (*,"('~',a,'~')") trim(tmp_ch) 	! DEBUG test the reading capability

                ! parse the line with format:
                !          m[-]           s[A]               e/k[K]              kapAB[-]            eAB/k[K]
                read(20, *) param_m(what), param_sigma(what), param_eps_kb(what), param_kap_asoc(what), param_eps_asoc_kb(what), &
                        !          u[Debye]       n_u[-]           q[DA]          n_q[-]           M[kg/mol]       bond
                        & param_u(what), param_n_u(what), param_q(what), param_n_q(what), param_mm(what), param_bond(what)

                ! correct the position vector
                pos = pos - shift
            end do
        end if
        ! initialize the interaction types
        call set_interaction_sub()

        ! 			! DEBUG -full printout of read parameters
        ! 				write (*,*) 'm | s | e/k | kapAB | eAB/k | u | n_u | q | n_q | M | bond'
        ! 			do i=1,n_comp
        ! 				write (*,"(f12.6,'|',f12.6,'|',f12.6,'|',f12.6,'|',f12.6,'|',f12.6,'|',f12.6,'|',f12.6,'|',f12.6,'|',f12.6,'|',f12.6)") &
        ! 				& param_m(i), param_sigma(i), param_eps_kb(i), param_kap_asoc(i), param_eps_asoc_kb(i), param_u(i), param_n_u(i), param_q(i), param_n_q(i), param_mm(i), param_bond(i)
        ! 			end do
        is_available = .true.
        close(20)
    end subroutine initialize_param_sub

    subroutine set_param_sub(subst_ind, par_in)
        ! procedure for direct setting of parameters: m | s | e/k | kapAB | eAB/k | u | n_u | q | n_q | M | bond
        ! operates with parameters for single substance
        ! accepts the substance index and parameter vector in prescribed order
        use control_mod, only : dp, &
                &   init_flag_gij, init_flag_zhs, init_flag_zhc, init_flag_zdisp, init_flag_zdd, init_flag_zqq, init_flag_zdq, init_flag_zasoc, &
                &   init_flag_ahs, init_flag_ahc, init_flag_adisp, init_flag_add, init_flag_aqq, init_flag_adq, init_flag_aasoc, &
                &   init_flag
        implicit none
        ! variable section
        integer, intent(in) :: subst_ind                ! substance number - where to save the parameters [1<=#<n_comp]
        real(dp), dimension(:), intent(in) :: par_in    ! parameter array, size [11]

        if (subst_ind .gt. n_comp) then 
            write (*, *) 'Error in function set_param_sub subst_ind > n_comp'
            stop
        end if 

        param_m(subst_ind) = par_in(1)
        param_sigma(subst_ind) = par_in(2)
        param_eps_kb(subst_ind) = par_in(3)
        param_kap_asoc(subst_ind) = par_in(4)
        param_eps_asoc_kb(subst_ind) = par_in(5)
        param_u(subst_ind) = par_in(6)
        param_n_u(subst_ind) = par_in(7)
        param_q(subst_ind) = par_in(8)
        param_n_q(subst_ind) = par_in(9)
        param_mm(subst_ind) = par_in(10)
        param_bond(subst_ind) = par_in(11)

        ! the dependedent variables are recomputed with the change of the param
        init_flag_gij = .false.
        init_flag_zhs = .false.
        init_flag_zhc = .false.
        init_flag_zdisp = .false.
        init_flag_zdd = .false.
        init_flag_zqq = .false.
        init_flag_zdq = .false.
        init_flag_zasoc = .false.

        init_flag_ahs = .false.
        init_flag_ahc = .false.
        init_flag_adisp = .false.
        init_flag_add = .false.
        init_flag_aqq = .false.
        init_flag_adq = .false.
        init_flag_aasoc = .false.
        ! this option ensures that all dependent function will be computed again
        init_flag(2) = .false.
        init_flag(3) = .false.
        init_flag(4) = .false.
        init_flag(6) = .false.

        ! can be made optional via parameter of dynamic -> use the knowledge of preivous+new values
        call set_interaction_sub()

    end subroutine set_param_sub

    subroutine get_param_sub(subst_ind, par_in)
        ! procedure for direct getting of parameters: m | s | e/k | kapAB | eAB/k | u | n_u | q | n_q | M | bond
        ! operates with parameters for single substance
        ! accepts the substance index and parameter vector that is filled in the prescribed order
        use control_mod, only : dp
        implicit none
        ! variable section
        integer, intent(in) :: subst_ind                ! substance number - where to save the parameters [1<=#<n_comp]
        real(dp), dimension(:), intent(inout) :: par_in    ! parameter array, size [11]

        par_in(1) = param_m(subst_ind)
        par_in(2) = param_sigma(subst_ind)
        par_in(3) = param_eps_kb(subst_ind)
        par_in(4) = param_kap_asoc(subst_ind)
        par_in(5) = param_eps_asoc_kb(subst_ind)
        par_in(6) = param_u(subst_ind)
        par_in(7) = param_n_u(subst_ind)
        par_in(8) = param_q(subst_ind)
        par_in(9) = param_n_q(subst_ind)
        par_in(10) = param_mm(subst_ind)
        par_in(11) = param_bond(subst_ind)

    end subroutine get_param_sub

    subroutine set_interaction_sub()
        ! sets up the available interaction types
        ! important notes and comments
        use control_mod, only : is_dd, is_qq, is_dq, is_asoc
        implicit none
        ! variable section

        ! count the available interaction count
        if (sum(param_u) > 0) then ! BUG assume the dipolar moment is non-negative
            is_dd = .true.
        else
            is_dd = .false.
        end if

        if (sum(param_q) >0) then ! BUG assume the quadrupolar moment is non-negative
            is_qq = .true.
        else
            is_qq = .false.
        end if

        if ((is_dd .eqv. .true.) .and. (is_qq .eqv. .true.)) then
            is_dq = .true.
            is_dq = .false. ! BUG DEBUG to disable dipolar-quadrupolar interaction
        else
            is_dq = .false.
        end if
        ! 		write (*,*) is_dd, is_qq ! DEBUG test the criteria result

        ! association type of interaction testing -> only when at least one is associating
        ! NOTE not optimal because it could be type not selfassociating
        if (sum(param_bond) > 0) then
            is_asoc = .true.
        else
            is_asoc = .false.
        end if
        ! 		write (*,*) is_asoc ! DEBUG test the criteria result
    end subroutine set_interaction_sub


end module param_list_mod
