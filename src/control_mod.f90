module control_mod
    ! the general global parameters and settings
    ! - define the fundamental technical parameters sp,dp,qp
    ! - set up the components, names character lengths
    ! - set up in run persistent variables: binary interaction parameter, temperature, molar ration
    ! - define universal constants: R, Na, Kbol, meter to angstrom conversion unit
    ! - set up the control falgs for intermediate variables in this fixed order:
    !	* | g_ij | z_hs | z_hc | z_disp | a_hs | a_hc | a_disp |
    !	* if deemed necesarry for brevity - this feature can be reimplemeted
    ! - set up the control flags for lower level operation in this order:
    !	* | initialize_a0b0_disp_fun | initialize_dzeta_fun | initialize_combrules_fun | compute_ab_disp_fun |
    !	* if deemed necesarry for brevity - this feature can be reimplemeted
    !	* offer the more complicated initializations to be omitted when not required
    ! - initilize the intermediate and derivative holding variables
    !	* because memory is not an issue even unused variables are initialized
    ! - initialization functions implementation
    !
    !	== LAST MODIFICATIONS ==
    !	15.5.2017 - D.C - descriptions and comments
    !	24.5.2017 - D.C - envelop the variable computation into get methods
    !					- get methods ensure all required is computed
    !					- unifies calls accross expanding contributions
    !	20.6.2017 - D.C - init_flag_g dependency omitment -> handling in this module
    !					- init_flag_g is no longer vetor (not practical & confusing)
    !	21.06.2017 - D.C - write_array subroutine implementation
    !	04.07.2017 - D.C - moving the z_hs contribution to get method in control mod
    !	??.07.2017 - D.C - polarity implemented with aditional checking and control in get methods
    !	2?.08.2017 - D.C - dynamic n_comp implementation via allocatable approach
    !	16.10.2017 - D.C - asociation interaface implemented
    !   23.04.2024 - D.C - add logging, debug messages features
    implicit none

    public !list the public methods

    ! precision settings

    integer, parameter :: sp = selected_real_kind(6, 37)
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: qp = selected_real_kind(33, 4931)

    ! output units definitions
    integer, parameter :: stdout = 6 !standard output
    integer, parameter :: stderr = 9 !standard error = 0
    integer, parameter :: stddeb = 10 !custom debug
    integer, parameter :: stdlog = 11 !custom log
    ! output units initialization control flag
    logical :: init_unit = .false.

    ! parameters section
    integer :: n_comp ! number of components
    ! TODO generalize n_comp into multi-component case with one of:
    !	* dynamic allocation - THIS OPTION IS CHOOSEN
    !	* static allocation with preprocessing (setting the correct n_comp)
    !	* static allocation with larger n_comp and implement limiter in cycles
    integer, parameter :: ch40_length = 40 ! number of characters for names of components
    integer, parameter :: ch120_length = 120 ! number of characters for the text parsing
    integer, parameter :: ch512_length = 512 ! number of characters for the file input parsing
    ! character(len = ch120_length), parameter :: n_of_file = "parameter_list.txt" ! the name of input text file
    ! character(len = ch120_length), parameter :: n_of_file = "parameter_list_trend4.txt" ! the name of input text file
    character(len = ch120_length), parameter :: n_of_file = "parameter_list.txt" ! the name of input text file
    character(len = ch120_length), parameter :: default_res_path = "../../res/" ! the path to resource directory from default execution location
    real(dp), parameter :: alpha_dq = 1.19374 ! fixed parameter for A3_DQ contribution

    ! variables persistent in one run
    real(dp) :: rho                                ! number density used in computation
    real(dp), dimension(:, :), allocatable :: k_ij ! the binary interaction parameter
    real(dp) :: tt ! temperature
    real(dp), dimension(:), allocatable :: mol_rat ! the imputed molar ratio

    logical :: sec_der  ! control variable that requests computation of the second density derivative
    logical :: composition_der ! control variable that request computation of the composition derivatives d/dxi
    logical :: temperature_der ! control variable that request computation of the temperature derivatives d/dT

    ! universal constant section
    ! 	real(dp), parameter :: r_gas= 8.3144598         ! gas constant [J/(K*mol)]
    ! 	real(dp), parameter :: n_a= 6.022140857e23      ! avogadro constant [1/mol]

    ! real(dp), parameter :: r_gas = 8.31441_dp         ! gas constant [J/(K*mol)] ! VV constants
    ! real(dp), parameter :: r_gas = 8.3144598_dp         ! gas constant [J/(K*mol)] ! TREND 4.1 constants
    real(dp), parameter :: r_gas = 8.31451_dp         ! gas constant [J/(K*mol)] ! TREND 4.1 constants for HFE calculation
    ! real(dp), parameter :: r_gas = 8.314462618_dp         ! gas constant [J/(K*mol)] ! CLAPEYRON
    ! real(dp), parameter :: n_a = 6.022045e23_dp       ! avogadro constant [1/mol] ! VV constants
    real(dp), parameter :: n_a = 6.022140857e23_dp      ! avogadro constant [1/mol] TREND 4.1 constant
    ! real(dp), parameter :: n_a = 6.02214076e23_dp      ! avogadro constant [1/mol] CLAPEYRON

    ! real(dp), parameter :: k_b = r_gas / n_a            ! boltzman constant [J/K]
    ! real(dp), parameter :: k_b = 1.38064852e-23_dp      ! boltzman constant [J/K] ! TREND 4.1 constant 
    real(dp), parameter :: k_b = 1.380656845862835e-23_dp      ! boltzman constant [J/K] ! adjust to agree with trend
    
    ! real(dp), parameter :: pi = 3.141592653589793_dp    ! The pie constant
    real(dp), parameter :: pi = 3.14159265359_dp    ! The pie constant  ! TREND 4.1 constant 
    real(dp), parameter :: m2angstrom = 1e10_dp         ! for the meter to angstrom conversion

    ! init_falg_* holds the information of availability of following variable *
    logical :: init_flag_gij = .false.
    logical :: init_flag_zhs = .false.
    logical :: init_flag_zhc = .false.
    logical :: init_flag_zdisp = .false.
    logical :: init_flag_zdd = .false.
    logical :: init_flag_zqq = .false.
    logical :: init_flag_zdq = .false.
    logical :: init_flag_zasoc = .false.

    logical :: init_flag_ahs = .false.
    logical :: init_flag_ahc = .false.
    logical :: init_flag_adisp = .false.
    logical :: init_flag_add = .false.
    logical :: init_flag_aqq = .false.
    logical :: init_flag_adq = .false.
    logical :: init_flag_aasoc = .false.

    ! the array that hold the initialization progress for the control of double initializing
    !| initialize_a0b0_disp_fun | initialize_dzeta_fun | initialize_combrules_fun | compute_ab_disp_fun |
    !  initialize_a0b0c0_polar_fun | initialize_polrules_fun |
    logical, dimension(6) :: init_flag = (/.false., .false., .false., .false., .false., .false./) ! TODO widen for other than dispersion

    ! interaction types indicators set up in param_list_mod initialization
    logical :: is_dd ! indicator if dipol-dipol contribution should be computed
    logical :: is_qq ! indicator if quadrupol-quadrupol contribution should be computed
    logical :: is_dq ! indicator if dipol-quadrupol contribution should be computed
    logical :: is_asoc ! indicator if association is taking place and should be computed

    ! tmp variables valid only as prevention of cyclic inclusion of modules
    real(dp) :: m_mean_tmp ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.5)

    ! runtime variables valid for single run
    ! BUG at first implemented with direct access -can cause problems
    ! ===== density derivative ====
    ! ----- RDF -----
    real(dp), dimension(:, :), allocatable :: g_ij ! the radial distribution function
    real(dp), dimension(:, :), allocatable :: g_ij_drho  ! the first density derivative of rdf
    real(dp), dimension(:, :), allocatable :: g_ij_drho2 ! the second density derivative of rdf
    real(dp), dimension(:, :), allocatable :: g_ij_drho3 ! the third density derivative of rdf
    ! 	real(dp), dimension(:,:), allocatable :: g_ij_drho4 ! the fourth density derivative of rdf
    ! 	real(dp), dimension(:,:), allocatable :: g_ij_drho5 ! the fifth density derivative of rdf
    ! ----- HS -----
    real(dp) :: z_hs       ! the hard sphere compressibility contribution
    real(dp) :: z_hs_drho  ! the first derivative of hard sphere contribution
    real(dp) :: z_hs_drho2 ! the second derivative of hard sphere contribution
    ! ----- HC -----
    real(dp) :: z_hc       ! the hard chain compressibility contribution
    real(dp) :: z_hc_drho  ! the first derivative of hard chain contribution
    real(dp) :: z_hc_drho2 ! the second derivative of hard chain contribution
    ! 	real(dp) :: z_hc_drho3 ! the third derivative of hard chain contribution
    ! 	real(dp) :: z_hc_drho4 ! the fourth derivative of hard chain contribution
    ! ----- DISP -----
    real(dp) :: z_disp       ! the dispersion compressibility contribution
    real(dp) :: z_disp_drho  ! the first derivative of dispersion contribution
    real(dp) :: z_disp_drho2 ! the second derivative of dispersion contribution
    ! ----- DD -----
    real(dp) :: z_dd       ! the dipol-dipol interaction compressibility contribution
    real(dp) :: z_dd_drho  ! the first derivative of dipol-dipol interaction compressibility contribution
    real(dp) :: z_dd_drho2 ! the second derivative of dipol-dipol interaction compressibility contribution
    ! ----- QQ -----
    real(dp) :: z_qq       ! the quadrupol-quadrupol interaction compressibility contribution
    real(dp) :: z_qq_drho  ! the first derivative of quadrupol-quadrupol interaction compressibility contribution
    real(dp) :: z_qq_drho2 ! the second derivative of quadrupol-quadrupol interaction compressibility contribution
    ! ----- DQ -----
    real(dp) :: z_dq       ! the dipol-quadrupol interaction compressibility contribution
    real(dp) :: z_dq_drho  ! the first derivative of dipol-quadrupol interaction compressibility contribution
    real(dp) :: z_dq_drho2 ! the second derivative of dipol-quadrupol interaction compressibility contribution
    ! ----- ASOC -----
    real(dp) :: z_asoc       ! the association interaction compressibility contribution
    real(dp) :: z_asoc_drho  ! the first derivative of association interaction compressibility contribution
    real(dp) :: z_asoc_drho2 ! the second derivative of association interaction compressibility contribution

    ! ===== plain value ====
    ! ----- HS -----
    real(dp) :: a_hs  ! the hard sphere helmholtz energy contribution
    ! ----- HC -----
    real(dp) :: a_hc  ! the hard chain helmholtz energy contribution
    ! ----- DISP -----
    real(dp) :: a_disp ! the dispersion helmholtz energy contribution
    ! ----- DD -----
    real(dp) :: a_dd   ! the dipol-dipol interaction helmholtz energy contribution
    ! ----- QQ -----
    real(dp) :: a_qq   ! the quadrupol-quadrupol interaction helmholtz energy contribution
    ! ----- DQ -----
    real(dp) :: a_dq   ! the dipol-quadrupol interaction helmholtz energy contribution
    ! ----- ASOC -----
    real(dp) :: a_asoc ! the association interaction helmholtz energy contribution

    ! ===== composition derivative ====
    ! ----- RDF -----
    real(dp), dimension(:, :, :), allocatable :: g_ij_dx  ! the first composition derivative of rdf
    ! ----- HS -----
    real(dp), dimension(:), allocatable :: a_hs_dx  ! the first composition derivative of hard sphere contribution vector in n_comp
    ! ----- HC -----
    real(dp), dimension(:), allocatable :: a_hc_dx  ! the first composition derivative of hard chain contribution vector in n_comp
    ! ----- DISP -----
    real(dp), dimension(:), allocatable :: a_disp_dx ! the first composition derivative of dispersion contribution vector in n_comp
    ! ----- DD -----
    real(dp), dimension(:), allocatable :: a_dd_dx   ! the first composition derivative of dipol-dipol contribution vector in n_comp
    ! ----- QQ -----
    real(dp), dimension(:), allocatable :: a_qq_dx   ! the first composition derivative of quadrupol-quadrupol contribution vector in n_comp
    ! ----- DQ -----
    real(dp), dimension(:), allocatable :: a_dq_dx   ! the first composition derivative of dipol-quadrupol contribution vector in n_comp
    ! ----- ASOC -----
    real(dp), dimension(:), allocatable :: a_asoc_dx ! the first composition derivative of association contribution vector in n_comp

    ! ===== temperature derivative ====
    ! ----- RDF -----
    real(dp), dimension(:, :), allocatable :: g_ij_dtt  ! the first temperature derivative of rdf
    ! ----- HS -----
    real(dp) :: a_hs_dtt  ! the first temperature derivative of hard sphere contribution
    ! ----- HC -----
    real(dp) :: a_hc_dtt  ! the first temperature derivative of hard chain contribution
    ! ----- DISP -----
    real(dp) :: a_disp_dtt ! the first temperature derivative of dispersion contribution
    ! ----- DD -----
    real(dp) :: a_dd_dtt   ! the first temperature derivative dipol-dipol interaction contribution
    ! ----- QQ -----
    real(dp) :: a_qq_dtt   ! the first temperature derivative quadrupol-quadrupol interaction contribution
    ! ----- DQ -----
    real(dp) :: a_dq_dtt   ! the first temperature derivative dipol-quadrupol interaction contribution
    ! ----- ASOC -----
    real(dp) :: a_asoc_dtt ! the first temperature derivative association interaction contribution

    ! ===== intermediate results ====
    ! accessed through the appropriate get methods
    ! NOTE not required with get functions - prohigit the direct unsafe acces
    ! 	real(dp) :: z_res ! the residual contribution of compresibilty factor - computed as sum of included contribution without ID
    ! 	real(dp) :: a_res ! the residual contribution of compresibilty factor - computed as sum of included contribution without ID
contains
    !the subroutines and functions
    ! ========================= ALLOCATION/ DEALLOCATION METHODS =============================
    ! ------------- ALLOCATION ---------------
    subroutine allocate_control_sub()
        ! allocate all related variables for controll module
        ! NOTE expect that n_comp is known nonzero
        if (n_comp == 0) then
            write (*, *) 'Error in function allocate_control n_comp==0'
            stop
        end if
        allocate(k_ij(n_comp, n_comp))
        allocate(mol_rat(n_comp))

        allocate(g_ij(n_comp, n_comp))
        allocate(g_ij_drho(n_comp, n_comp))
        allocate(g_ij_drho2(n_comp, n_comp))
        allocate(g_ij_drho3(n_comp, n_comp))
        ! 			allocate( g_ij_drho4(n_comp,n_comp))
        ! 			allocate( g_ij_drho5(n_comp,n_comp))

        allocate(g_ij_dx(n_comp, n_comp, n_comp))
        allocate(a_hs_dx(n_comp))
        allocate(a_hc_dx(n_comp))
        allocate(a_disp_dx(n_comp))
        allocate(a_dd_dx(n_comp))
        allocate(a_qq_dx(n_comp))
        allocate(a_dq_dx(n_comp))
        allocate(a_asoc_dx(n_comp))

        allocate(g_ij_dtt(n_comp, n_comp))

    end subroutine allocate_control_sub
        
    ! ------------- DEALLOCATION -------------
    subroutine deallocate_control_sub()
        ! deallocate all related variables for controll module
        ! performs check if memory is allocated otherwise does nothing to unallocated memory

        if (allocated(k_ij) .eqv. .true.) then 
            deallocate(k_ij)
        end if
        if (allocated(mol_rat) .eqv. .true.) then 
            deallocate(mol_rat)
        end if
        
        if (allocated(g_ij) .eqv. .true.) then 
            deallocate(g_ij)
        end if
        if (allocated(g_ij_drho) .eqv. .true.) then 
            deallocate(g_ij_drho)
        end if
        if (allocated(g_ij_drho2) .eqv. .true.) then 
            deallocate(g_ij_drho2)
        end if
        if (allocated(g_ij_drho3) .eqv. .true.) then 
            deallocate(g_ij_drho3)
        end if
        ! if (allocated(g_ij_drho4) .eqv. .true.) then 
        !     deallocate(g_ij_drho4)
        ! end if
        ! if (allocated(g_ij_drho5) .eqv. .true.) then 
        !     deallocate(g_ij_drho5)
        ! end if

        if (allocated(g_ij_dx) .eqv. .true.) then 
            deallocate(g_ij_dx)
        end if
        if (allocated(a_hs_dx) .eqv. .true.) then 
            deallocate(a_hs_dx)
        end if
        if (allocated(a_hc_dx) .eqv. .true.) then 
            deallocate(a_hc_dx)
        end if
        if (allocated(a_disp_dx) .eqv. .true.) then 
            deallocate(a_disp_dx)
        end if
        if (allocated(a_dd_dx) .eqv. .true.) then 
            deallocate(a_dd_dx)
        end if
        if (allocated(a_qq_dx) .eqv. .true.) then 
            deallocate(a_qq_dx)
        end if
        if (allocated(a_dq_dx) .eqv. .true.) then 
            deallocate(a_dq_dx)
        end if
        if (allocated(a_asoc_dx) .eqv. .true.) then 
            deallocate(a_asoc_dx)
        end if
        
        if (allocated(g_ij_dtt) .eqv. .true.) then 
            deallocate(g_ij_dtt)
        end if
              
        ! make sure unallocated values are not used
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
        
    end subroutine deallocate_control_sub
    
    
    ! =================================== GET METHODS ========================================
    ! ------------- compresibilty factor related -------------
    function get_zres_fun(rho_in) result (z_res)
        ! provide acces to the z_res in a safe way ensuring all is computed
        ! variable section
        real(dp), intent(in) :: rho_in ! input density number
        real(dp) :: z_res ! the residual contribution of compresibilty factor - coputed as sum of included contribution without ID

        ! control section of required z_res variables
        if (init_flag_zhs .eqv. .false.) then
            call z_hs_sub(rho_in) ! require gij beforehand
            ! NOTE already multiplied by m_mean
        end if
        if (init_flag_zhc .eqv. .false.) then
            call z_hc_sub(rho_in) ! require gij beforehand
        end if
        if (init_flag_zdisp .eqv. .false.) then
            call z_disp_sub(rho_in)
        end if

        if (is_dd .eqv. .true.) then
            if (init_flag_zdd .eqv. .false.) then
                call z_dd_sub(rho_in) ! this function also set the a_qq
            end if
        else
            z_dd = 0.0_dp
        end if
        if (is_qq .eqv. .true.) then
            if (init_flag_zqq .eqv. .false.) then
                call z_qq_sub(rho_in) ! this function also set the a_qq
            end if
        else
            z_qq = 0.0_dp
        end if
        if (is_dq .eqv. .true.) then
            if (init_flag_zdq .eqv. .false.) then
                call z_dq_sub(rho_in) ! this function also set the a_qq
            end if
        else
            z_dq = 0.0_dp
        end if
        if (is_asoc .eqv. .true.) then
            if (init_flag_zasoc .eqv. .false.) then
                call z_asoc_sub(rho_in) ! this function also set the a_asoc
            end if
        else
            z_asoc = 0.0_dp
        end if

        z_res = z_hs + z_hc + z_disp + z_dd + z_qq + z_dq + z_asoc
        ! write(*,"('z_res = ',e)"), z_res ! DEBUG

        ! write(*,*) "z_res : ",   z_res
        ! write(*,*) "z_hs  : ",    z_hs
        ! write(*,*) "z_hc  : ",    z_hc
        ! write(*,*) "z_hs+c: ",    z_hs + z_hc
        ! write(*,*) "z_disp: ",  z_disp
        ! write(*,*) "z_qq  : ",    z_qq
        ! write(*,*) "z_dd  : ",    z_dd
        ! write(*,*) "z_asoc: ",  z_asoc
        
    end function get_zres_fun

    function get_zres_drho_fun(rho_in) result (z_res_drho)
        ! provide acces to the z_res_drho in a safe way ensuring all is computed
        implicit none
        ! variable section
        real(dp), intent(in) :: rho_in ! input density number
        real(dp) :: z_res_drho ! the residual contribution of compresibilty factor - coputed as sum of included contribution without ID

        ! control section of required z_res variables
        if (init_flag_zhs .eqv. .false.) then
            call z_hs_sub(rho_in) ! require gij beforehand
            ! NOTE already multiplied by m_mean
        end if
        if (init_flag_zhc .eqv. .false.) then
            call z_hc_sub(rho_in) ! require gij beforehand
        end if
        if (init_flag_zdisp .eqv. .false.) then
            call z_disp_sub(rho_in)
        end if
        if (is_dd .eqv. .true.) then
            if (init_flag_zdd .eqv. .false.) then
                call z_dd_sub(rho_in) ! this function also set the a_qq
            end if
        else
            z_dd_drho = 0.0_dp
        end if
        if (is_qq .eqv. .true.) then
            if (init_flag_zqq .eqv. .false.) then
                call z_qq_sub(rho_in) ! this function also set the a_qq
            end if
        else
            z_qq_drho = 0.0_dp
        end if
        if (is_dq .eqv. .true.) then
            if (init_flag_zdq .eqv. .false.) then
                call z_dq_sub(rho_in) ! this function also set the a_qq
            end if
        else
            z_dq_drho = 0.0_dp
        end if
        if (is_asoc .eqv. .true.) then
            if (init_flag_zasoc .eqv. .false.) then
                call z_asoc_sub(rho_in) ! this function also set the a_asoc
            end if
        else
            z_asoc_drho = 0.0_dp
        end if

        z_res_drho = z_hs_drho + z_hc_drho + z_disp_drho + z_dd_drho + z_qq_drho + z_dq_drho + z_asoc_drho
        ! 			write(*,"('z_res_drho = ',e)") z_res_drho ! DEBUG
    end function get_zres_drho_fun

    function get_zres_drho2_fun(rho_in) result (z_res_drho2)
        ! provide acces to the z_res_drho2 in a safe way ensuring all is computed
        implicit none
        ! variable section
        real(dp), intent(in) :: rho_in ! input density number
        real(dp) :: z_res_drho2 ! the residual contribution of compresibilty factor - coputed as sum of included contribution without ID

        if (sec_der .eqv..false.) then
            write (*, *) 'Error in function get_zres_drho2_fun: cannot acces second derivative when it is not requested by setting sec_der = .true.'
            stop
        end if
        ! control section of required z_res variables
        if (init_flag_zhs .eqv. .false.) then
            call z_hs_sub(rho_in) ! require gij beforehand
            ! NOTE already multiplied by m_mean
        end if
        if (init_flag_zhs .eqv. .false.) then
            call z_hc_sub(rho_in) ! require gij beforehand
        end if
        if (init_flag_zdisp .eqv. .false.) then
            call z_disp_sub(rho_in)
        end if
        if (is_dd .eqv. .true.) then
            if (init_flag_zdd .eqv. .false.) then
                call z_dd_sub(rho_in) ! this function also set the a_qq
            end if
        else
            z_dd_drho2 = 0.0_dp
        end if
        if (is_qq .eqv. .true.) then
            if (init_flag_zqq .eqv. .false.) then
                call z_qq_sub(rho_in) ! this function also set the a_qq
            end if
        else
            z_qq_drho2 = 0.0_dp
        end if
        if (is_dq .eqv. .true.) then
            if (init_flag_zdq .eqv. .false.) then
                call z_dq_sub(rho_in) ! this function also set the a_qq
            end if
        else
            z_dq_drho2 = 0.0_dp
        end if
        if (is_asoc .eqv. .true.) then
            if (init_flag_zasoc .eqv. .false.) then
                call z_asoc_sub(rho_in) ! this function also set the a_asoc
            end if
        else
            z_asoc_drho2 = 0.0_dp
        end if

        z_res_drho2 = z_hs_drho2 + z_hc_drho2 + z_disp_drho2 + z_dd_drho2 + z_qq_drho2 + z_dq_drho2 + z_asoc_drho2
        ! 			write(*,"('z_res_drho2 = ',e)") z_res_drho2 ! DEBUG
    end function get_zres_drho2_fun

    ! ------------- reduced helmholtz energy related -------------
    function get_ares_fun(rho_in) result (a_res)
        ! provide acces to the a_res in a safe way ensuring all is computed
        ! NOTE TODO evaluate if a_res can be private
        implicit none
        ! variable section
        real(dp), intent(in) :: rho_in ! input density number
        real(dp) :: a_res ! the residual reduced helmholtz energy - coputed as sum of included contribution without ID

        ! control section of required a_res variables
        if (init_flag_ahs .eqv. .false.) then
            call a_hs_sub(rho_in) ! require gij beforehand
            ! NOTE not multiplied by m_mean - new version for numerical derivative control
        end if
        if (init_flag_ahc .eqv. .false.) then
            call a_hc_sub(rho_in) ! required gij
        end if
        if (init_flag_adisp .eqv. .false.) then
            call a_disp_sub(rho_in)
        end if
        if (is_dd .eqv. .true.) then
            if (init_flag_add .eqv. .false.) then
                call z_dd_sub(rho_in) ! this function also set the a_qq
            end if
        else
            a_dd = 0.0_dp
        end if
        if (is_qq .eqv. .true.) then
            if (init_flag_aqq .eqv. .false.) then
                call z_qq_sub(rho_in) ! this function also set the a_qq
            end if
        else
            a_qq = 0.0_dp
        end if
        if (is_dq .eqv. .true.) then
            if (init_flag_adq .eqv. .false.) then
                call z_dq_sub(rho_in) ! this function also set the a_qq
            end if
        else
            a_dq = 0.0_dp
        end if
        if (is_asoc .eqv. .true.) then
            if (init_flag_aasoc .eqv. .false.) then
                call z_asoc_sub(rho_in) ! this function also set the a_asoc
            end if
        else
            a_asoc = 0.0_dp
        end if

        a_res = a_hs + a_hc + a_disp + a_dd + a_qq + a_dq + a_asoc
        ! write(*,"('a_res = ',e)") a_res ! DEBUG

        ! write(*,*) "a_res : ",   a_res
        ! write(*,*) "a_hs  : ",    a_hs
        ! write(*,*) "a_hc  : ",    a_hc
        ! write(*,*) "a_hs+c: ",    a_hs + a_hc
        ! write(*,*) "a_disp: ",  a_disp
        ! write(*,*) "a_qq  : ",    a_qq
        ! write(*,*) "a_dd  : ",    a_dd
        ! write(*,*) "a_asoc: ",  a_asoc
        
    end function get_ares_fun

    subroutine get_ares_dx_fun(rho_in, a_res_dx)
        ! provide acces to the a_res_dx in a safe way ensuring all is computed
        ! NOTE TODO evaluate if a_res_dx can be private
        implicit none
        ! variable section
        real(dp), intent(in) :: rho_in ! input density number
        real(dp), dimension(:), intent(inout) :: a_res_dx ! the residual reduced helmholtz energy temperature derivative - coputed as sum of included contribution without ID
        integer :: i ! counter
        if (composition_der .eqv..false.) then
            write (*, *) 'Error in function get_ares_dx_fun: cannot acces composition derivative when it is not requested by setting composition_der = .true.'
            stop
        end if
        ! control section of required a_res variables
        if (init_flag_ahs .eqv. .false.) then
            call a_hs_sub(rho_in) ! require gij beforehand
            ! NOTE not multiplied by m_mean and added parts - new version for numerical derivative control
            ! derivative should be: dm_mean/dx_i + m_mean*dz_hs/dx_i
        end if
        if (init_flag_ahc .eqv. .false.) then
            call a_hc_sub(rho_in) ! require gij beforehand
        end if
        if (init_flag_adisp .eqv. .false.) then
            call a_disp_sub(rho_in)
        end if
        if (is_dd .eqv. .true.) then
            if (init_flag_add .eqv. .false.) then
                call z_dd_sub(rho_in) ! this function also set the a_qq
            end if
        else
            a_dd_dx = 0.0_dp
        end if
        if (is_qq .eqv. .true.) then
            if (init_flag_aqq .eqv. .false.) then
                call z_qq_sub(rho_in) ! this function also set the a_qq
            end if
        else
            a_qq_dx = 0.0_dp
        end if
        if (is_dq .eqv. .true.) then
            if (init_flag_adq .eqv. .false.) then
                call z_dq_sub(rho_in) ! this function also set the a_qq
            end if
        else
            a_dq_dx = 0.0_dp
        end if
        if (is_asoc .eqv. .true.) then
            if (init_flag_zasoc .eqv. .false.) then
                call z_asoc_sub(rho_in) ! this function also set the a_asoc
            end if
        else
            a_asoc_dx = 0.0_dp
        end if

        if (size(a_res_dx) == n_comp) then
            do i = 1, n_comp
                a_res_dx(i) = a_hs_dx(i) + a_hc_dx(i) + a_disp_dx(i) + a_dd_dx(i) + a_qq_dx(i) + a_dq_dx(i) + a_asoc_dx(i)
            end do
        else
            write (*, *) 'Error in function get_ares_dx_fun: the dimension of a_res_dx does not conform size(a_res_dx) /= n_comp '
            stop
        end if
        ! 			a_res_dx = a_hs_dx + a_hc_dx + a_disp_dx + a_dd_dx ! the array addition
        ! write(*,*) "a_res_dx : ",   a_res_dx
        ! write(*,*) "a_hs_dx : ",   a_hs_dx
        ! write(*,*) "a_hc_dx : ",   a_hc_dx
        ! write(*,*) "a_disp_dx : ",   a_disp_dx
        ! write(*,*) "a_dd_dx : ",   a_dd_dx
        ! write(*,*) "a_qq_dx : ",   a_qq_dx
        ! write(*,*) "a_dq_dx : ",   a_dq_dx
        ! write(*,*) "a_asoc_dx : ",   a_asoc_dx
    end subroutine get_ares_dx_fun

    function get_ares_dtt_fun(rho_in) result (a_res_dtt)
        ! provide acces to the a_res_dtt in a safe way ensuring all is computed
        ! NOTE TODO evaluate if a_res_dtt can be private
        implicit none
        ! variable section
        real(dp), intent(in) :: rho_in ! input density number
        real(dp) :: a_res_dtt ! the residual reduced helmholtz energy temperature derivative - coputed as sum of included contribution without ID

        if (temperature_der .eqv..false.) then
            write (*, *) 'Error in function get_ares_dtt_fun: cannot acces temperature derivative when it is not requested by setting temperature_der = .true.'
            stop
        end if
        ! control section of required a_res variables
        if (init_flag_ahs .eqv. .false.) then
            call a_hs_sub(rho_in) ! require gij beforehand
            ! NOTE not multiplied by m_mean - new version for numerical derivative control
        end if
        if (init_flag_ahc .eqv. .false.) then
            call a_hc_sub(rho_in) ! require rdf beforehand
        end if
        if (init_flag_adisp .eqv. .false.) then
            call a_disp_sub(rho_in)
        end if
        if (is_dd .eqv. .true.) then
            if (init_flag_add .eqv. .false.) then
                call z_dd_sub(rho_in) ! this function also set the a_qq
            end if
        else
            a_dd_dtt = 0.0_dp
        end if
        if (is_qq .eqv. .true.) then
            if (init_flag_aqq .eqv. .false.) then
                call z_qq_sub(rho_in) ! this function also set the a_qq
            end if
        else
            a_qq_dtt = 0.0_dp
        end if
        if (is_dq .eqv. .true.) then
            if (init_flag_adq .eqv. .false.) then
                call z_dq_sub(rho_in) ! this function also set the a_qq
            end if
        else
            a_dq_dtt = 0.0_dp
        end if
        if (is_asoc .eqv. .true.) then
            if (init_flag_aasoc .eqv. .false.) then
                call z_asoc_sub(rho_in) ! this function also set the a_asoc
            end if
        else
            a_asoc_dtt = 0.0_dp
        end if

        a_res_dtt = a_hs_dtt + a_hc_dtt + a_disp_dtt + a_dd_dtt + a_qq_dtt + a_dq_dtt + a_asoc_dtt
        ! 			write(*,"('a_res_dtt = ',e)") a_res_dtt ! DEBUG
    end function get_ares_dtt_fun
    ! ============================= INITIALIZATION SECTION ==================================
!     subroutine set_n_comp(n_comp_in)
!         ! set the number of component and allocate them if required
!         ! important notes and comments
!         implicit none
!         integer :: n_comp_in ! the number of components
!         ! variable section        
!         integer :: i, j ! counters
!         
!         if ((n_comp_in .lt. 1) .or. (n_comp_in .gt. 99)) then
!             write (*, *) 'Error in set_n_comp: number of components is not in range <0,99>'
!             stop
!         end if
!                 
!         if (n_comp .ne. n_comp_in) then
!             if (allocated(mol_rat) .eqv. .true.) then
!                 call deallocate_all_sub()
!             end if
!             n_comp = n_comp_in
!             call allocate_all_sub()
!         end if
!                 
!     end subroutine
!     
    subroutine set_k_ij(k_ij_in)
        ! initialize the k_ij matrix with the k_ij binary interaction parameter
        ! important notes and comments
        implicit none
        ! variable section
        real(dp) :: k_ij_in ! the binary interaction parameter
        integer :: i, j ! counters
        if (n_comp==1) then
            k_ij = 0.0_dp;
        else
            do i = 1, n_comp
                do j = 1, n_comp
                    if (i==j) then
                        k_ij(i, j) = 0.0_dp
                    else
                        k_ij(i, j) = k_ij_in
                    end if
                end do
            end do
            ! 				if (n_comp==2) then
            ! 					k_ij = reshape( (/0.0_dp,k_ij_in,k_ij_in,0.0_dp/), (/n_comp,n_comp/) ) ! init the 2*2 array
            ! 				else
            ! ! 					write (*,*) 'Error in function set k_ij: unimplemented for more then 2 components'
            ! ! 					stop
            ! 				end if
        end if

    end subroutine

    subroutine initialize_vars_sub(tt_in, mol_rat_in, k_ij_in, sec_der_in, composition_der_in, temperature_der_in)
        ! set up the dispersion parameters
        ! important notes and comments
        implicit none
        ! variable section
        real(dp) :: tt_in ! temperature K
        real(dp) :: k_ij_in ! the binary interaction parameter
        real(dp), dimension(n_comp) :: mol_rat_in ! the imputed molar ratio ! BUG not dynamic allocation
        logical, optional :: sec_der_in, composition_der_in, temperature_der_in ! second derivative computation flag default = true

        tt = tt_in
        mol_rat = mol_rat_in
        call set_k_ij(k_ij_in) ! initialization of k_ij

        if(present(sec_der_in)) then
            sec_der = sec_der_in
        else
            sec_der = .false.
        end if
        if(present(composition_der_in)) then
            composition_der = composition_der_in
        else
            composition_der = .false.
        end if
        if(present(temperature_der_in)) then
            temperature_der = temperature_der_in
        else
            temperature_der = .false.
        end if

        ! this option ensures that all global variables are recomputed
        init_flag_gij = .false.
        init_flag_zhs = .false.
        init_flag_zhc = .false.
        init_flag_zdisp = .false.
        init_flag_zdd = .false.
        init_flag_zqq = .false.
        init_flag_zdq = .false.

        init_flag_ahs = .false.
        init_flag_ahc = .false.
        init_flag_adisp = .false.
        init_flag_add = .false.
        init_flag_aqq = .false.
        init_flag_adq = .false.
        init_flag_aasoc = .false.
        ! this option ensures that all dependent function will be computed again
        ! init_flag(2) = .false.
        ! init_flag(3) = .false.
        ! init_flag(4) = .false.
        ! init_flag(6) = .false. ! TODO decide if necessary to initialize polar rules
    end subroutine initialize_vars_sub

    subroutine set_density_sub()
        ! reset the flag for density dependent variables
        ! important notes and comments
        ! usefull for resetting the flags
        implicit none
!        real(dp), intent(in) :: dens_in ! the imput density
!
!        rho = dens_in

        ! this option ensures that all global variables are recomputed
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

        ! 			write(*, "('Density dependent variable flags reset')") ! DEBUG control over function acces to flag
    end subroutine set_density_sub

    subroutine set_temperature_sub(tt_in)
        ! set up the temperature variable and ensures that all temperature dependent variables will be computed with new temperature
        ! important notes and comments
        implicit none
        ! variable section
        real(dp), intent(in) :: tt_in ! temperature K

        tt = tt_in
        ! 			write (*,"('Temperature set to: ',e)") tt ! DEBUG for the temperature changes observations
        ! this option ensures that all global variables are recomputed
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
        ! 			init_flag(1) = .false.! DEBUG further unnecesarry
        init_flag(2) = .false.
        init_flag(3) = .false.
        ! 			init_flag(4) = .false. ! DEBUG further unnecesarry
        init_flag(6) = .false. ! TODO decide if necessary to initialize polar rules for temperature der
        ! 			write(*, "('Temperature dependent variable flags reset')") ! DEBUG control over function acces to flag
    end subroutine set_temperature_sub

    subroutine set_composition_sub(mol_rat_in)
        ! set up the temperature variable and ensures that all temperature dependent variables will be computed with new temperature
        ! important notes and comments
        implicit none
        ! variable section
        real(dp), dimension(n_comp), intent(in) :: mol_rat_in ! the imputed molar ratio ! BUG not dynamic allocation
        
        if ((sum(mol_rat_in)>1.0_dp).or.(sum(mol_rat_in)<1.0_dp)) then ! additonal control
            write (*, *) 'Error in initialization: molar ratios does not sum to one'
            stop
        end if
        mol_rat = mol_rat_in

        ! 			write (*,"('Composition set to: ',e)") mol_rat ! DEBUG for the composition changes observations
        ! this option ensures that all global variables are recomputed
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
        ! 			init_flag(6) = .false. ! TODO decide if necessary to initialize polar rules for composition der
        ! 			write(*, "('Composition dependent variable flags reset')") ! DEBUG control over function acces to flag
    end subroutine set_composition_sub

    ! ============================= HELPER SECTION ==================================
    subroutine write_timestamp(write_unit)
        ! write the timestam into the write_unit
        implicit none
        integer, intent(in) :: write_unit
        ! variable section
        character(8)  :: date
        character(10) :: time
        character(5)  :: zone
        integer, dimension(8) :: values

        ! using keyword arguments
        call date_and_time(date,time,zone,values)
        call date_and_time(DATE=date,ZONE=zone)
        call date_and_time(TIME=time)
        ! call date_and_time(VALUES=values)

        write(write_unit,"('** TIMESTAMP: ',a,':',a,':',a)", advance="no") time(1:2), time(3:4), time(5:6)
        write(write_unit,"(' ',a,'-',a,'-',a)", advance="no") date(7:8), date(5:6), date(1:4)
        write(write_unit,"(' ',a,' **')") zone
        ! print '(8i5)', values

    end subroutine

    subroutine init_ouput_units()
        ! open and write the basic information into the output unit headers
        ! usually used for logging and debug information
        implicit none

        ! variable section
        integer :: file_success !  the status identificator

        if (init_unit .eqv. .false.) then
            ! stdout,stderr are handled by the system - now setup required
            open(stderr, file = "error.txt", status = 'replace', iostat = file_success)
            if (file_success /= 0) then
                write(stdout,*) "Unable to create file: error.txt"
            end if
            write(stderr,"('** ****************ERROR*************** **')")
            call write_timestamp(stderr)
            write(stderr,"('** ************************************ **')")

            ! stddeb custom file for debug logging
            open(stddeb, file = "debug.txt", status = 'old', position = "append", iostat = file_success, form = 'formatted', action='write')
            if (file_success /= 0) then
                write(stderr,*) "Unable to open file: debug.txt, make sure it exists"
            end if
            ! write(stddeb, *) "This is DEBUG logging file."
            write(stddeb,"('** ****************DEBUG*************** **')")
            call write_timestamp(stddeb)
            write(stddeb,"('** ************************************ **')")

            ! stdlog custom file for debug logging
            open(stdlog, file = "log.txt", status = 'old', position = "append", iostat = file_success, form = 'formatted', action='write')
            if (file_success /= 0) then
                write(stderr,*) "Unable to open file: log.txt, make sure it exists"
            end if
            ! write(stdlog, *) "This is LOG logging file."
            write(stdlog,"('** *****************LOG**************** **')")
            call write_timestamp(stdlog)
            write(stdlog,"('** ************************************ **')")

            init_unit = .true.
        end if

    end subroutine

    subroutine write_array_1d(array_in, header, write_unit_in)
        ! the debug/helper subroutine that writes the array to console with following header
        ! important notes and comments
        implicit none
        real(dp), dimension(:), intent(in) :: array_in
        character(*), intent(in) :: header
        integer, intent(in), optional :: write_unit_in
        ! variable section
        integer :: i
        integer :: write_unit !the actuall used write unit

        if (present(write_unit_in) .eqv. .true. ) then
            write_unit = write_unit_in
        else
            ! else use the default stdout
            write_unit = stdout
        end if

        write(write_unit, "('** ',a,' **')") header ! newline

        do i = 1, size(array_in, 1)
            write(write_unit, "(e16.8,', ')", advance = "no") array_in(i)
        end do
        ! 				write(write_unit"(e)",advance="yes") array_in(i+1)
        write(write_unit, *) ! newline
    end subroutine write_array_1d

    subroutine write_array_2d(array_in, header, write_unit_in)
        ! the debug/helper subroutine that writes the array to console with following header
        ! important notes and comments
        implicit none
        real(dp), dimension(:, :), intent(in) :: array_in
        character(*), intent(in) :: header
        integer, intent(in), optional :: write_unit_in
        ! variable section
        integer :: i, j
        integer :: write_unit !the actuall used write unit

        if (present(write_unit_in) .eqv. .true. ) then
            write_unit = write_unit_in
        else
            ! else use the default stdout
            write_unit = stdout
        end if

        write(write_unit, "('** ',a,' **')") header ! newline

        do i = 1, size(array_in, 1)
            do j = 1, size(array_in, 2)
                write(write_unit, "(e16.8,', ')", advance = "no") array_in(i, j)
            end do
            ! 					write(write_unit"(e)",advance="yes") array_in(i,j+1)
            write(write_unit, *) ! newline
        end do
        write(write_unit, *) ! newline
    end subroutine write_array_2d

    subroutine write_array_3d(array_in, header, write_unit_in)
        ! the debug/helper subroutine that writes the array to console with following header
        ! important notes and comments
        implicit none
        real(dp), dimension(:, :, :), intent(in) :: array_in
        character(*), intent(in) :: header
        integer, intent(in), optional :: write_unit_in
        ! variable section
        integer :: i, j, k
        integer :: write_unit !the actuall used write unit

        if (present(write_unit_in) .eqv. .true. ) then
            write_unit = write_unit_in
        else
            ! else use the default stdout
            write_unit = stdout
        end if

        write(write_unit, "('** ',a,' **')") header ! newline

        do i = 1, size(array_in, 1)
            write(write_unit, "('** i= ',i2,' **')") i
            do j = 1, size(array_in, 2)
                do k = 1, size(array_in, 3)
                    write(write_unit, "(e16.8,', ')", advance = "no")array_in(j, k, i)
                end do
                ! 							write(write_unit"(e)",advance="yes") array_in(i,j,k+1)
                write(write_unit, *) ! newline
            end do
            write(write_unit, *) ! newline
        end do
        write(write_unit, *) ! newline
    end subroutine write_array_3d

end module control_mod
