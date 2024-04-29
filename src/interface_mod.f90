module interface_mod
    ! tough implementation of simple interface for pc-saft program.
    ! the interface is parsed from within this mod in rather lengthy and crude way open to extension
    ! mod also contains direct hardwired input
    !
    !	== LAST MODIFICATIONS ==
    !	23.5.2017 - D.C - implementation, descriptions and comments
    !	15.7.2017 - D.C - correction of stric float input of composition (valid for 1-component system)
    !	2?.8.2017 - D.C - modified interface with prepended number of components
    !	07.8.2018 - D.C - setup_flags_sub implementation

    use control_mod, only : dp, n_comp, ch40_length, ch120_length, ch512_length, &
            &    stdout, stderr, stddeb, stdlog, &
            &    rho, k_ij, tt, mol_rat, &
            &    r_gas, n_a, m2angstrom, &
            &    sec_der, composition_der, temperature_der, &
            &    set_k_ij, &
            &    default_res_path!, &

    use contrib_mod, only : dzeta, &
            &    initialize_dzeta_fun
    use param_list_mod, only : &
            &    initialize_param_sub
    
    implicit none
    !list the public methods
    public comp_opt, option, &
            & is_pfrac, is_dense, &
            & parse_input_universal, hardwired_input, &
            & setup_flags_parse_sub, set_fluid, &
            & allocate_all_sub, deallocate_all_sub, &
            & allocate_interface_sub, deallocate_interface_sub, &
            & my_getarg, &
            & names, parse_line ! this has to be visible for testing
    private
    ! parameter section (general)
!    character(len = 1), parameter :: lclose = "'", rclose = "'" ! the closure characters
    character(len = 1), parameter :: lclose = ":", rclose = ":" ! the closure characters temporary because of matlab messing quotes
    character(len = 1), parameter :: delimiter = "," ! the delimiter character
    ! parameter section (default values)

    integer, parameter :: n_comp_def = 2                                                            ! default number of components
    character(len = ch40_length), dimension(n_comp_def), parameter :: names_def = (/'methane', 'ethane '/)    ! default names for 2 components
    real(dp), parameter :: rho_def = 997.0_dp                                                        ! density [kg/m^3]
    real(dp), parameter :: tt_def = 300_dp                                                            ! temperature [°K]
    real(dp), dimension(n_comp_def), parameter :: mol_rat_def = (/0.57_dp, 0.43_dp/)                    ! molar ratios for 2 components
    real(dp), parameter :: k_ij_def = 0.0_dp                                                        ! binary interaction parameter
    !logical, dimension(9), parameter :: comp_opt_def = (/.true., .false., .false., .false., .false., .false., .false., .false., .false./)
    !integer, parameter :: option_def = 1															! option for calculation numerical derivative test
    ! VV notes
    logical, dimension(9), parameter :: comp_opt_def = (/.true., .true., .true., .true., .true., .true., .true., .true., .true./)
    integer, parameter :: option_def = 1                                                            ! option for calculation numerical derivative test
    ! variable section
    character(len = ch40_length), dimension(:), allocatable :: names    ! names input, size [n_comp]
    logical, dimension(9) :: comp_opt                ! compute options that holds what user requested to be computed
    ! contains: | p | dp/drho | dp/drho2 | mu | fug | A | H | S | G |
    ! can contain packing fraction input switch pfrac
    logical :: is_pfrac                                ! the switch controlling packing fraction recalculation to density
    logical :: is_dense                                ! the switch controlling if density is written to output
    integer :: option                                ! provide further runtime control over numerical derivatives if enabled
    !1= density, 2= temperature, 3= composition (default component set to first - change in num_diff_calc_sub)
    integer :: i                                    ! counter

contains
    ! ========================= ALLOCATION / DEALLOCATION / INITIALIZATION METHODS =============================
    ! ------------- ALLOCATION ---------------
    subroutine allocate_interface_sub()
        ! subroutine calls all other allocation subroutines of important modules
        ! NOTE expect that n_comp is known nonzero
        write (stddeb,*) "==> allocate_interface_sub" !DEBUG
        if (n_comp == 0) then
            write (stderr, *) 'Error in function allocate_control n_comp==0'
            stop
        end if
        allocate(names(n_comp))

        write (stddeb,*) "<== allocate_interface_sub" !DEBUG
    end subroutine allocate_interface_sub
    
    subroutine allocate_all_sub()
        ! subroutine calls all other allocation subroutines of important modules
        ! NOTE expect that n_comp is known nonzero
!        use interface_mod, only: allocate_interface_sub
        use param_list_mod, only: allocate_param_sub
        use control_mod, only: allocate_control_sub
        use contrib_mod, only: allocate_contrib_sub
                
        write (stddeb,*) "=> allocate_all_sub" !DEBUG
        call allocate_interface_sub()
        call allocate_param_sub()
        call allocate_control_sub()
        call allocate_contrib_sub()

        write (stddeb,*) "<= allocate_all_sub" !DEBUG
    end subroutine allocate_all_sub
    
    ! ------------- DEALLOCATION -------------
    subroutine deallocate_interface_sub()
        ! subroutine calls all other allocation subroutines of important modules
        ! performs check if memory is allocated otherwise does nothing to unallocated memory

        write (stddeb,*) "==> deallocate_interface_sub" !DEBUG
        if (allocated(names) .eqv. .true.) then 
            deallocate(names)
        end if

        write (stddeb,*) "<== deallocate_interface_sub" !DEBUG
    end subroutine deallocate_interface_sub
    
    subroutine deallocate_all_sub()
        ! subroutine calls all other allocation subroutines of important modules
!        use interface_mod, only: deallocate_interface_sub
        use param_list_mod, only: deallocate_param_sub
        use control_mod, only: deallocate_control_sub
        use contrib_mod, only: deallocate_contrib_sub
        
        write (stddeb,*) "=> deallocate_all_sub" !DEBUG
        call deallocate_interface_sub()
        call deallocate_param_sub()
        call deallocate_control_sub()
        call deallocate_contrib_sub()

        write (stddeb,*) "<= deallocate_all_sub" !DEBUG
    end subroutine deallocate_all_sub
    
    ! ------------- INITIALIZATION -------------
    subroutine setup_flags_hardwired_sub(comp_opt_in)
        ! provides unified subroutine for setting the execution flags
        ! in a sense overloaded to cover hardwired input
        ! to solve the problem with unsetted flags for hardwired input
        implicit none
        ! variable section
        logical, dimension(9), intent(in) :: comp_opt_in
        is_pfrac = .false. ! initialization of packing fraction switch

        if (comp_opt_in(3).eqv. .true.) then ! composition_dereviative
            sec_der = .true.
        end if

        if ((comp_opt_in(4).eqv. .true.) .or. (comp_opt_in(5).eqv. .true.)) then
            composition_der = .true.
        end if

        if ((comp_opt_in(7).eqv. .true.) .or. (comp_opt_in(8).eqv. .true.) .or. (comp_opt_in(9).eqv. .true.)) then
            temperature_der = .true.
        end if

    end subroutine setup_flags_hardwired_sub

    subroutine setup_flags_parse_sub(arg_input)
        ! provides unified subroutine for setting the execution flags
        ! in a sense overloaded to cover parsed input
        ! to solve the problem with unsetted flags for hardwired input
        implicit none
        ! variable section
        character(len = *), intent(in) :: arg_input
        character(len = 5) :: tmp_arg ! temporary argument helper for parsing the control options
        character(len = ch120_length) :: buff_arg ! buffer storage for the crunched arg_input
        integer :: item_count

        ! setting the initial values - effective reset
        comp_opt = .false. !set the whole vector

        is_dense = .false.
        is_pfrac = .false. ! initialization of packing fraction switch

        sec_der = .false.
        composition_der = .false.
        temperature_der = .false.

        option = 0

        ! pre-parsing ot the input
        buff_arg = arg_input ! put the original into buffer
        item_count = count(transfer(buff_arg, 'a', len(buff_arg)) == ",") + 1 ! count the number of ',' which equals number of component -1
        do i = 1, item_count ! parse the first item_count-1 control options
            if (i==item_count) then
                tmp_arg = trim(buff_arg) ! the last argument case
            else
                tmp_arg = buff_arg(1:index(buff_arg, ',', .false.) - 1) ! middle argument case
                buff_arg = buff_arg(index(buff_arg, ',', .false.) + 1:)
            end if

            select case (tmp_arg)
            case('rho')
                is_dense = .true.
            case('p')
                comp_opt(1) = .true.
            case('pp')
                comp_opt(2) = .true.
            case('ppp')
                sec_der = .true.
                comp_opt(3) = .true.
            case('mu')
                composition_der = .true.
                comp_opt(4) = .true.
            case('fug')
                composition_der = .true.
                comp_opt(5) = .true.
            case('a')
                comp_opt(6) = .true.
                ! NOTE because of pressure
                comp_opt(1) = .true.
!                 NOTE because of chemical potential
                composition_der = .true.
                comp_opt(4) = .true.
            case('h')
                temperature_der = .true.
                comp_opt(7) = .true.
            case('s')
                temperature_der = .true.
                comp_opt(8) = .true.
            case('g')
                temperature_der = .true.
                comp_opt(9) = .true.
            case('1')
                option = 1
            case('2')
                option = 2
            case('3')
                option = 3
            case('pfrac')
                is_pfrac = .true.
            case('')
                ! NOTE this case do nothing
                ! enables clearing of the parse flags
            case default
                write (stderr, "('Error in initialization with wrong control option: <',a,'>',/,'Check the header for help',/)") trim(tmp_arg)
                stop
            end select
        end do

    end subroutine setup_flags_parse_sub

    ! =================================== PARSER METHODS ========================================
    subroutine hardwired_input()
        ! supplies the hardwired variant of input
        ! used in the older versions of code
        implicit none
        ! variable section
        real(dp) :: k_ij_ ! the temporary binary interaction parameter
        integer :: is_initialized ! status values of param_list_mod initialization

        ! VV-change
        sec_der = .true.
        composition_der = .true.
        temperature_der = .true.

        n_comp = n_comp_def
        call allocate_all_sub()
        do i = 1, n_comp
            names(i) = names_def(i)
            mol_rat(i) = mol_rat_def(i)
        end do
        rho = rho_def * n_a / m2angstrom**3 ! conversion to the computation dimension
        tt = tt_def
        k_ij_ = k_ij_def
        comp_opt = comp_opt_def
        option = option_def
        call set_k_ij(k_ij_)

        ! output section
        write (*, "('!===== Default input parameters are: =====')") ! for input control

        write (*, "('! n_component: |',i2,'|')") n_comp !DEBUG just for verification and user review of input

        do i = 1, n_comp - 1 ! parse the first item_count-1 component names
            write (*, "('! names(',i2,'):   |',a,'| ')") i, trim(names(i)) !DEBUG just for verification and user review of input
        end do
        write (*, "('! names(',i2,'):   |',a,'|')") n_comp, trim(names(n_comp)) !DEBUG just for verification and user review of input
        call initialize_param_sub(names, is_initialized)
        if (is_initialized/=0) then
            write (stderr, *) 'Error in function initialize_param_sub'
            stop
        end if

        do i = 1, n_comp - 1 ! parse the first item_count-1 component names
            write (*, "('! mol_rat(',i2,'): |',f10.4,'| ')") i, mol_rat(i) !DEBUG just for verification and user review of input
        end do
        write (*, "('! mol_rat(',i2,'): |',f10.4,'|')") n_comp, mol_rat(n_comp) !DEBUG just for verification and user review of input

        write (*, "('! temperature: |',f10.4,'|')") tt !DEBUG just for verification and user review of input

        write (*, "('! input dens:  |',f10.4,'|')") rho_def !DEBUG just for verification ! Edited to display in same units as for parsed input

        write (*, "('! k_ij:        |',f10.4,'|')") k_ij_ !DEBUG just for verification

        write (*, "('! computation options: |p: ',l,' |pp: ',l,' |ppp: ',l,' |mu: ',l,' |fug: ',l,' |a: ',l,' |h: ',l,' |s: ',l,' |g: ',l,'|')") & !DEBUG just for verification and user review of input
                & comp_opt(1), comp_opt(2), comp_opt(3), comp_opt(4), comp_opt(5), comp_opt(6), comp_opt(7), comp_opt(8), comp_opt(9)
    end subroutine

    subroutine my_getarg(input_text, position, output_text)
        ! get appropriate argument at position from input_text to output_text
        implicit none
        character(len = ch512_length), intent(in) :: input_text
        character(len = 2 * ch40_length), intent(out) :: output_text
        integer, intent(inout) :: position
        ! variable section
        character, parameter :: quote = ":"
        integer :: start, stop, count, ii
        integer :: is_quoted

        ii = 1 ! skip the first space
        count = 0
        start = 1 ! start after the first space
        is_quoted = 0
        do while (ii < ch512_length)
            !            write (*,*) input_text(ii: ii) ! DEBUG
            if (input_text(ii:ii) == quote) then
                is_quoted = modulo(is_quoted + 1, 2) ! flip the true/false
            end if
            if ((input_text(ii:ii) == ' ') .and. (is_quoted == 0)) then
                count = count + 1
                if (count == position) then
                    stop = ii !does not matter as it will be trimmed
                    exit
                else
                    start = ii !omit the space character
                end if
            end if
            ii = ii + 1
        end do
        position = position+1
        output_text = input_text(start:stop)
        output_text = trim(adjustl(output_text))
    end subroutine my_getarg

    function get_n_substances_from_names(subst_text) result(n_subst)
        ! for given input subst_text perform counting of delimiters to
        ! determine the number of substances used later as n_comp
        !
        implicit none
        integer :: n_subst
        character(len = *) :: subst_text
        ! variable section
        integer :: found_index
        integer :: i, tmp_cnt

        n_subst = 0
        if(len(subst_text) == 0) then
            write (stderr, *) 'Error in subst_count: empty input char array'
            stop
        end if

        i = 1
        do
            found_index = index(subst_text(i:), delimiter)
            n_subst = n_subst + 1 ! the number of substances is one larger than found delimiters

            if(found_index == 0) then ! last delimiter found
                return
            end if
            i = i + found_index + 1 ! move after the search after the last found index
        end do

    end function get_n_substances_from_names

    function parse_line(line_input, is_write_to_console) result(error_flag)
        ! parse the given line input independent on the source
        ! parsing is sectioned into individual parsing blocks
        ! "number of component(s)"SPACE"component(s)"SPACE"molar ratio(s)"SPACE"temperature"SPACE"density"SPACE"k_ij"SPACE"computed properties flag(s)"
        implicit none
        integer :: error_flag
        character(len = ch512_length), intent(in) :: line_input
        logical, optional :: is_write_to_console
        logical :: write_to_console
        ! variable section
        integer :: i, item_count, tmp_subst_count
        character(len = 2 * ch40_length) :: arg_input
        real(dp) :: k_ij_ ! the temporary binary interaction parameter
        integer :: is_initialized ! status values of param_list_mod initialization
        integer :: arg_number ! the number of argument requested

        ! handle the optional argument
        if (present(is_write_to_console)) then
            write_to_console = is_write_to_console
        else
            write_to_console = .false.
        end if

        if (write_to_console .eqv. .true. ) then
            write (*, "('!===== Accepted input parameters are: =====')") ! for input control
        end if
        error_flag = 0
        arg_number = 1 ! initialization to request first argument

        ! === n_comp ===
        call my_getarg(line_input, arg_number, arg_input)
        arg_input = trim(arg_input)
        tmp_subst_count = get_n_substances_from_names(arg_input)
        if (tmp_subst_count .eq. 0) then
            write (stderr, *) 'Error in initialization: no substances given'
        end if

        n_comp = tmp_subst_count
        if (write_to_console .eqv. .true. ) then
            write (*, "('! n_component: |',i2,'|')") n_comp !DEBUG just for verification and user review of input
        end if
        ! === allocate section ===
        call allocate_all_sub()

        ! === names parse section ===
        ! NOTE arg_input is already loaded no need to reload
        ! 			write (*,"('|',a,'|')") arg_input !DEBUG just for verification of parsed input
        if (len_trim(arg_input)==0) then
            if (n_comp == 2) then
                names = ''
                do i = 1, n_comp
                    names(i) = names_def(i)
                end do
                write (*, *) 'Default value for names used.'
            else
                write (stderr, *) 'Error in initialization: wrong number of components'
                error_flag = 1                
            end if
        else
            ! strip the boundary brackets
!            write (*,"('|',a,'|')") arg_input !DEBUG just for verification of parsed input
            if ((index(arg_input, lclose, .false.) /= 0) .and. (index(arg_input, rclose, .true.) /= 0)) then
                arg_input = arg_input(index(arg_input, lclose, .false.) + 1:index(arg_input, rclose, .true.) - 1) ! OLD removed with change to "string" input format
!                write (*,"('|',a,'|')") arg_input !DEBUG just for verification of parsed input
            end if
            if (arg_input=='') then
                write (stderr, *) 'Error in initialization: wrong input form for names - check the header for help'
                error_flag = 1
            end if
            item_count = count(transfer(arg_input, 'a', len(arg_input)) == ",") + 1 ! count the number of ',' which equals number of component -1

            if (item_count/=n_comp) then
                write (stderr, *) 'Error in initialization: number of inputed components does not conform with expected number of components'
                error_flag = 1
            end if

            do i = 1, n_comp - 1 ! parse the first item_count-1 component names
                names(i) = arg_input(1:index(arg_input, ',', .false.) - 1) ! not including the ','
                arg_input = arg_input(index(arg_input, ',', .false.) + 1:)
                if (write_to_console .eqv. .true. ) then
                    write (*, "('! names(',i2,'):   |',a,'| ')") i, trim(names(i)) !DEBUG just for verification and user review of input
                end if
            end do
            names(n_comp) = trim(arg_input) ! the last component
            if (write_to_console .eqv. .true. ) then
                write (*, "('! names(',i2,'):   |',a,'|')") n_comp, trim(names(n_comp)) !DEBUG just for verification and user review of input
            end if
        end if
!        write (*, *) 'DEBUG: before initialize_param_sub' ! DEBUG
        ! == initialization of parameters ==
        call initialize_param_sub(names, is_initialized, write_to_console)

        if (is_initialized/=0) then
            write (stderr, *) 'Error in function initialize_param_sub'
            error_flag = 1
        end if
!         write (*, *) 'DEBUG: after initialize_param_sub' ! DEBUG
        ! === mol_rat parse section ===
        call my_getarg(line_input, arg_number, arg_input)
        ! 			write (*,"('|',a,'|')") arg_input !DEBUG just for verification of parsed input
        if (len_trim(arg_input)==0) then
            if (n_comp == 2) then
                mol_rat = 0.0_dp
                do i = 1, n_comp
                    mol_rat(i) = mol_rat_def(i)
                end do
                write (*, *) 'Default value for molrates used.'
            else
                write (stderr, *) 'Error in initialization: wrong number of components - default values are only for n_comp = 2'
                error_flag = 1
            end if
        else
            ! strip the boundary brackets
            ! 				arg_input = arg_input(index(arg_input,lclose,.false.)+1:index(arg_input,rclose,.true.)-1) ! OLD removed with change to "string" input format
            ! 		write (*,"('|',a,'|')") arg_input !DEBUG just for verification

            if (arg_input=='') then
                write (stderr, *) 'Error in initialization: wrong input form for mol ratios - check the header for help'
                error_flag = 1
            end if
            item_count = count(transfer(arg_input, 'a', len(arg_input)) == ",") + 1 ! count the number of ',' which equals number of component -1

            if (item_count/=n_comp) then
                write (stderr, *) 'Error in initialization: number of inputed mol ratios does not conform with expected number of components'
                error_flag = 1
            end if

            ! 				OLD this is obsolete with the list directed format read
            ! 				do i=1,n_comp-1 ! parse the first item_count-1 component names
            ! 					read(arg_input(1:index(arg_input,',',.false.)-1),*) mol_rat(i)
            ! 					arg_input = arg_input(index(arg_input,',',.false.)+1:)
            ! 					write (*,"('! mol_rat(',i2,'): |',f10.4,'| ')") i, mol_rat(i) !DEBUG just for verification and user review of input
            ! 				end do
            !
            ! 				read(arg_input,*) mol_rat(n_comp)  ! the last component

            read(arg_input, *) mol_rat  ! the last component
            do i = 1, n_comp - 1 ! parse the first item_count-1 component names
                if (write_to_console .eqv. .true. ) then
                    write (*, "('! mol_rat(',i2,'): |',f10.4,'| ')") i, mol_rat(i) !DEBUG just for verification and user review of input
                end if
            end do
            if (write_to_console .eqv. .true. ) then
                write (*, "('! mol_rat(',i2,'): |',f10.4,'|')") n_comp, mol_rat(n_comp) !DEBUG just for verification and user review of input
            end if
            if ((sum(mol_rat)>1.0_dp).or.(sum(mol_rat)<1.0_dp)) then ! additonal control
                write (stderr, *) 'Error in initialization: molar ratios does not sum to one'
                error_flag = 1
            end if
        end if

        ! === temperature ===
        call my_getarg(line_input, arg_number, arg_input)
        ! 			write (*,"('|',a,'|')") arg_input !DEBUG just for verification of parsed input
        if (len_trim(arg_input)==0) then
            tt = tt_def
            write (*, *) 'Default value for temperature used.'
        else
            ! strip the boundary brackets
            ! 				arg_input = arg_input(index(arg_input,lclose,.false.)+1:index(arg_input,rclose,.true.)-1) ! OLD removed with change to "string" input format
            ! 		write (*,"('|',a,'|')") arg_input !DEBUG just for verification

            if (arg_input=='') then
                write (stderr, *) 'Error in initialization: wrong input form for temperature - check the header for help'
                error_flag = 1
            end if

            read(arg_input, *) tt
            if (tt<=0.0_dp) then ! addittional control
                write (stderr, *) 'Error in initialization: temperature is not positive number'
                error_flag = 1
            end if
            if (write_to_console .eqv. .true. ) then
                write (*, "('! temperature: |',f10.4,'| [K]')") tt !DEBUG just for verification and user review of input
            end if
        end if


        ! === density ===
        call my_getarg(line_input, arg_number, arg_input)
        ! 			write (*,"('|',a,'|')") arg_input !DEBUG just for verification of parsed input
        if (len_trim(arg_input)==0) then
            rho = rho_def
            write (*, *) 'Default value for density used.'
        else
            ! strip the boundary brackets
            ! 				arg_input = arg_input(index(arg_input,lclose,.false.)+1:index(arg_input,rclose,.true.)-1) ! OLD removed with change to "string" input format
            ! 		write (*,"('|',a,'|')") arg_input !DEBUG just for verification

            if (arg_input=='') then
                write (stderr, *) 'Error in initialization: wrong input form for density - check the header for help'
                error_flag = 1
            end if

            read(arg_input, *) rho
            if (write_to_console .eqv. .true. ) then
                write (*, "('! input dens:  |',f10.4,'| [mol/m^3]')") rho !DEBUG just for verification
            end if
        end if


        ! === binary interaction parameter k_ij ===
        call my_getarg(line_input, arg_number, arg_input)
        ! 			write (*,"('|',a,'|')") arg_input !DEBUG just for verification of parsed input
        arg_input = trim(arg_input)
        if (arg_input=='') then
            k_ij_ = k_ij_def
        else
            ! strip the boundary brackets
            ! 				arg_input = arg_input(index(arg_input,lclose,.false.)+1:index(arg_input,rclose,.true.)-1) ! OLD removed with change to "string" input format
            ! 		write (*,"('|',a,'|')") arg_input !DEBUG just for verification

            if (arg_input=='') then
                write (stderr, *) 'Error in initialization: wrong input form for k_ij - check the header for help'
                error_flag = 1
            end if

            read(arg_input, *) k_ij_
            if (write_to_console .eqv. .true. ) then
                write (*, "('! k_ij:        |',f10.4,'|')") k_ij_ !DEBUG just for verification
            end if
        end if

        call set_k_ij(k_ij_) ! centralize inilialization of k_ij parameter

        ! === control options ===
        call my_getarg(line_input, arg_number, arg_input)
        ! 			write (*,"('|',a,'|')") arg_input !DEBUG just for verification of parsed input
        ! initialization of the helper variables
        option = 0
        !this contains: | p | dp/drho | dp/drho2 | mu | fug | H | S | G |
        comp_opt = (/.false., .false., .false., .false., .false., .false., .false., .false., .false./)

        if (len_trim(arg_input)==0) then
            comp_opt = comp_opt_def
            write (*, *) 'Default value for computation option used.'
        else
            ! strip the boundary brackets
            ! 				arg_input = arg_input(index(arg_input,lclose,.false.)+1:index(arg_input,rclose,.true.)-1) ! OLD removed with change to "string" input format
            ! 		write (*,*) arg_input !DEBUG just for verification

            if (arg_input=='') then
                write (stderr, *) 'Error in initialization: wrong input form for names - check the header for help'
                error_flag = 1                
            end if

            call setup_flags_parse_sub(arg_input)
            if (write_to_console .eqv. .true. ) then
                ! 			write (*,"('! options:  sd:|',l,' |',/,'          cd:|',l,' |',/,'          td:|',l,' |')") sec_der, composition_der, temperature_der !DEBUG just for verification and user review of input
                write (*, "('! computation options: |p: ',l,' |pp: ',l,' |ppp: ',l,' |mu: ',l,' |fug: ',l,' |a: ',l,' |h: ',l,' |s: ',l,' |g: ',l,'|')") & !DEBUG just for verification and user review of input
                        & comp_opt(1), comp_opt(2), comp_opt(3), comp_opt(4), comp_opt(5), comp_opt(6), comp_opt(7), comp_opt(8), comp_opt(9)
            end if
        end if

        if (is_pfrac .eqv. .true.) then ! correction of density if the input is in form of packing fraction
            if (initialize_dzeta_fun().eqv..false.) then
                write (stderr, *) 'Error in helper function initialize_dzeta_fun'
                error_flag = 1                
            end if
            rho = rho / dzeta(3);
        else
            ! NOTE: Input density - molar density (default settings)
            rho = rho * n_a / m2angstrom**3 ! conversion to the computation dimension
        end if
        ! 	stop ! TEST Just for testing the parser
    end function

    function parse_file(file_input, is_write_to_console) result(error_flag)
        ! handle the filename parsing
        ! parse the first uncommented line to as paremeters
        ! the comment line start with is ! or #
        implicit none
        integer :: error_flag !  the status identificator
        character(len = ch512_length), intent(in) :: file_input
        logical, optional :: is_write_to_console ! if the  output should be verbose - useful for batch execution
        logical :: write_to_console
        integer :: file_success !  the file operation identificator
        character(len = ch512_length) :: tmp_ch ! the temporary character
        logical :: file_existence

        ! handle the optional argument
        if (present(is_write_to_console)) then
            write_to_console = is_write_to_console
        else
            write_to_console = .false.
        end if

        open(21, file = trim(file_input), status = 'old', iostat = file_success)
        if (file_success /= 0) then
            if (write_to_console .eqv. .true.) then
                write (stdlog, "('INFO: file <',a,'> was not opened with iostat:',i4)") trim(file_input), file_success
                write (stdlog, "('      trying another file: <',a,'>')") trim(default_res_path)
            endif
            close(21)
            ! second try on loading the file in res directory
            open(21, file = trim(default_res_path) // adjustl(file_input), status = 'old', action='read', iostat = file_success)

            if (file_success /= 0) then
                ! BUG this section causes report error even though the file is read 
                ! if (write_to_console .eqv. .true.) then
                    ! write (stderr, "('Error: even the other was not opened with iostat:', i4)") file_success
                ! endif
                ! error_flag = file_success 
                error_flag = 0 ! BUG this bypass the error
            end if

        end if

        do ! read untill end of file
            read(21, '(a)', iostat = file_success) tmp_ch
            if (file_success /= 0) then
                exit ! exit the loop on EOF
            end if

            if (tmp_ch(1:10) == '          ') then ! to omit the empty lines
                cycle
            end if
            if ((tmp_ch(1:1) == '!') .or. (tmp_ch(1:1) == '#')) then ! to omit the commented lines
                cycle
            else
                error_flag = parse_line(tmp_ch,write_to_console)
                exit ! get the first line with calculation settings and exit
            end if
        end do

        close(21)

    end function parse_file

    subroutine parse_input_universal(input_filepath, is_write_to_console)
        ! handle the file or command line input
        implicit none
        character(len = ch512_length) :: arg_input
        ! variables
        character(len = ch512_length), optional :: input_filepath
        logical, optional :: is_write_to_console ! if the  output should be verbose - useful for batch execution
        logical :: write_to_console
        logical :: is_console_input
!        character(len = *), parameter :: name_of_program = "pcsaft_main_prog"
        character(len = *), parameter :: name_of_program = "pcsaft"
!        integer, parameter :: prog_name_length = 16
        integer, parameter :: prog_name_length = 6
        integer :: error_flag

        ! handle the optional argument
        if (present(input_filepath)) then
            is_console_input = .false.
        else
            is_console_input = .true.
        end if

        ! handle the optional argument
        if (present(is_write_to_console)) then
            write_to_console = is_write_to_console
        else
            write_to_console = .true.
        end if
        error_flag = 0

        if (is_console_input .eqv. .true.) then ! input from console
            call get_command(arg_input)
            ! NOTE do the search from back as there should be no name_of_program refrences
            arg_input = arg_input(index(arg_input, lclose, .false.):)
            ! write(*,*) arg_input !DEBUG
            if (command_argument_count()>=6) then ! the CLI input
                if (write_to_console .eqv. .true.) then
                    write (*, "('!===== Parsing parameters from CLI =====')") ! for input control
                endif
                error_flag = parse_line(arg_input, write_to_console)
            else ! the file input
                if (write_to_console .eqv. .true.) then
                    write (*, "('!===== Parsing parameters from FILE : ',(a),' =====')") trim(arg_input) ! for input control
                endif
                error_flag = parse_file(arg_input, write_to_console)
            end if
        else ! file input as argument
            if (write_to_console .eqv. .true.) then
                write (*, "('!===== Parsing parameters from FILE : ',(a),' =====')") trim(input_filepath) ! for input control
            endif
            error_flag = parse_file(input_filepath, write_to_console)
        end if
        if (error_flag /= 0) then
            write (stderr, *) "Error: general error in parse_input occured. ", error_flag
        else
            ! write (*, *) "error_flag", error_flag !DEBUG
        end if

    end subroutine parse_input_universal
    
    subroutine set_n_comp(n_comp_in)
        ! set the number of component and allocate them if required
        ! important notes and comments
        implicit none
        integer :: n_comp_in ! the number of components
        ! variable section        
        integer :: i, j ! counters
        
        if ((n_comp_in .lt. 1) .or. (n_comp_in .gt. 99)) then
            write (stderr, *) 'Error in set_n_comp: number of components is not in range <0,99>'
            stop
        end if
                
        if (n_comp .ne. n_comp_in) then
            if (allocated(mol_rat) .eqv. .true.) then
                call deallocate_all_sub()
            end if
            n_comp = n_comp_in
            call allocate_all_sub()
        end if
                
    end subroutine
        
    subroutine set_fluid(fluid_text, is_write_to_console)
        ! handle the file or command line input
        ! TODO this is duplicit code from parse_line function - factor the call of setting into the function
        implicit none
        character(len = *), intent(in) :: fluid_text
        logical, optional :: is_write_to_console! if the  output should be verbose - useful for batch execution
        logical :: write_to_console
        ! variable section
        integer :: i, item_count
        character(len = :), allocatable :: arg_input
        integer :: is_initialized ! status values of param_list_mod initialization
        integer :: n_comp_temp ! temporary value for number of components

        call set_k_ij(0.0_dp) ! initialization of k_ij

        ! handle the optional argument
        if (present(is_write_to_console)) then
          write_to_console = is_write_to_console
        else
          write_to_console = .false.
        end if

        arg_input = trim(fluid_text)
        ! 			write (*,"('|',a,'|')") arg_input !DEBUG just for verification of parsed input
        if (len_trim(arg_input) .eq. 0) then
            write (stderr, *) 'Error in initialization: wrong number of components'
            stop
        else
            n_comp_temp = get_n_substances_from_names(arg_input)            
            call set_n_comp(n_comp_temp)
            
            ! strip the boundary brackets
            ! write (*,"('|',a,'|')") arg_input !DEBUG just for verification of parsed input
            if ((index(arg_input, lclose, .false.) /= 0) .and. (index(arg_input, rclose, .true.) /= 0)) then
                arg_input = arg_input(index(arg_input, lclose, .false.) + 1:index(arg_input, rclose, .true.) - 1) ! OLD removed with change to "string" input format
                ! write (*,"('|',a,'|')") arg_input !DEBUG just for verification of parsed input
            end if
            if (arg_input=='') then
                write (stderr, *) 'Error in initialization: wrong input form for names - check the header for help'
                stop
            end if
!             NOTE this is obsolete as the counting of comas is used for number of component
!             item_count = count(transfer(arg_input, 'a', len(arg_input)) == ",") + 1 ! count the number of ',' which equals number of component -1
! 
!             if (item_count/=n_comp) then
!                 write (stderr, *) 'Error in initialization: number of inputed components does not conform with expected number of components'
!                 stop
!             end if

            do i = 1, n_comp - 1 ! parse the first item_count-1 component names
                names(i) = arg_input(1:index(arg_input, ',', .false.) - 1) ! not including the ','
                arg_input = arg_input(index(arg_input, ',', .false.) + 1:)
                if (write_to_console .eqv. .true.) then
                    write (*, "('! names(',i2,'):   |',a,'| ')") i, trim(names(i)) !DEBUG just for verification and user review of input
                endif
            end do
            names(n_comp) = trim(arg_input) ! the last component
            if (write_to_console .eqv. .true.) then
                write (*, "('! names(',i2,'):   |',a,'|')") n_comp, trim(names(n_comp)) !DEBUG just for verification and user review of input
            endif   
        end if
!        write (*, *) 'DEBUG: before initialize_param_sub' ! DEBUG
        ! == initialization of parameters ==
        call initialize_param_sub(names, is_initialized, write_to_console)
        if (is_initialized/=0) then
            write (stderr, *) 'Error in function initialize_param_sub'
            stop
        end if
!        write (*, *) 'DEBUG: after initialize_param_sub' ! DEBUG
    end subroutine

end module interface_mod
