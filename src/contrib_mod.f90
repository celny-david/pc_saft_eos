module contrib_mod
    ! the preparation module for the low-level parameters used in dispersion and polar contributions
    ! this module contain the hardwired parameter tables for speed purposes (can be changed in modification)
    ! - declare public variables
    ! - define variables in groups: dispersion, dipolar, quadrupolar, dipolar-quadrupolar, helpers
    ! - initialization functions and computation functions implementation
    !
    !	== LAST MODIFICATIONS ==
    !	15.5.2017 - D.C - descriptions and comments
    !	--.7.2017 - D.C - implementation of polarity(dd,qq,dq) related function
    !					- ee,eee type variables calculation and polarity terms

    use control_mod, only : dp, n_comp, k_b, m2angstrom, pi, &
            &    stdout, stderr, stddeb, stdlog, &
            &    tt, mol_rat, k_ij, &
            &    sec_der, composition_der, temperature_der, &
            &    init_flag, &
            &    set_k_ij
    use param_list_mod, only : param_m, param_sigma, param_eps_kb, param_eps_asoc_kb, param_kap_asoc
    implicit none

    public a__disp, b__disp, a__disp_dx, b__disp_dx, & ! list of basic dispersion parameters
            & a__dd, b__dd, c__dd, a__qq, b__qq, c__qq, a__dq, b__dq, c__dq, & ! list of basic polar parameters
            & dzeta, dzeta_dx, dzeta_dtt, d_seg, d_seg_dtt, m_mean, mol_rat, &
            & sigma_ij, eps_ij, eps_asoc_ij, kapa_asoc_ij, &
            & m2es3, m2e2s3, m2es3_dx, m2e2s3_dx, &
            & ees3s3_s3nnmm, eees3s3s3_sssnnnmmm, ees5s5_s7nnqq, eees5s5s5_s3s3s3nnnqqq, & ! list of public intermediate parameters
            & allocate_contrib_sub, deallocate_contrib_sub, &
            & initialize_a0b0_disp_fun, initialize_dzeta_fun, check_dzeta_sub, initialize_combrules_fun, compute_ab_disp_fun, compute_qudratic_sub !list the public methods

    ! variable section
    real(dp), dimension(7, 3) :: a0_disp ! dispersion approximation parameters a0
    real(dp), dimension(7, 3) :: b0_disp ! dispersion approximation parameters b0
    real(dp), dimension(7) :: a__disp ! dispersion approximation parameters a
    real(dp), dimension(7) :: b__disp ! dispersion approximation parameters b
    real(dp), dimension(:, :), allocatable :: a__disp_dx ! dispersion approximation parameters a, size [7,n_comp]
    real(dp), dimension(:, :), allocatable :: b__disp_dx ! dispersion approximation parameters b, size [7,n_comp]

    real(dp), dimension(5, 3) :: a0_dd ! approximation parameters a0 for dipolar component
    real(dp), dimension(5, 3) :: b0_dd ! approximation parameters b0 for dipolar component
    real(dp), dimension(5, 3) :: c0_dd ! approximation parameters c0 for dipolar component
    real(dp), dimension(:, :, :), allocatable :: a__dd ! approximation parameters a for dipolar component, size [5,n_comp,n_comp]
    real(dp), dimension(:, :, :), allocatable :: b__dd ! approximation parameters b for dipolar component, size [5,n_comp,n_comp]
    real(dp), dimension(:, :, :, :), allocatable :: c__dd ! approximation parameters c for dipolar component, size [5,n_comp,n_comp,n_comp]

    real(dp), dimension(5, 3) :: a0_qq ! approximation parameters a0 for quadrupolar component
    real(dp), dimension(5, 3) :: b0_qq ! approximation parameters b0 for quadrupolar component
    real(dp), dimension(5, 3) :: c0_qq ! approximation parameters c0 for quadrupolar component
    real(dp), dimension(:, :, :), allocatable :: a__qq ! approximation parameters a for quadrupolar component, size [5,n_comp,n_comp]
    real(dp), dimension(:, :, :), allocatable :: b__qq ! approximation parameters b for quadrupolar component, size [5,n_comp,n_comp]
    real(dp), dimension(:, :, :, :), allocatable :: c__qq ! approximation parameters c for quadrupolar component, size [5,n_comp,n_comp,n_comp]

    real(dp), dimension(5, 3) :: a0_dq ! approximation parameters a0 for dipolar-quadrupolar component
    real(dp), dimension(5, 3) :: b0_dq ! approximation parameters b0 for dipolar-quadrupolar component
    real(dp), dimension(5, 2) :: c0_dq ! approximation parameters c0 for dipolar-quadrupolar component
    real(dp), dimension(:, :, :), allocatable :: b__dq ! approximation parameters b for dipolar-quadrupolar component, size [5,n_comp,n_comp]
    real(dp), dimension(:, :, :), allocatable :: a__dq ! approximation parameters a for dipolar-quadrupolar component, size [5,n_comp,n_comp]
    real(dp), dimension(:, :, :, :), allocatable :: c__dq ! approximation parameters c for dipolar-quadrupolar component, size [5,n_comp,n_comp,n_comp]

    real(dp), dimension(0:3) :: dzeta ! the dzeta vector dzeta_n for n=0,1,2,3 without rho !!
    real(dp), dimension(:, :), allocatable :: dzeta_dx ! the dzeta matrix with derivatives ddzeta_n/dxi for n=0,1,2,3 without rho !! , size [0:3,1:n_comp]
    real(dp), dimension(0:3) :: dzeta_dtt ! the dzeta vector with derivatives ddzeta_n/dtt for n=0,1,2,3 without rho !!
    real(dp), dimension(:), allocatable :: d_seg ! the segment diameter di (temperature dependant), size [n_comp]
    real(dp), dimension(:), allocatable :: d_seg_dtt ! the segment diameter di (temperature dependant), size [n_comp]
    real(dp) :: m_mean ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.5)

    real(dp), dimension(:, :), allocatable :: sigma_ij     ! the combining parameter sigma, size [n_comp,n_comp]
    real(dp), dimension(:, :), allocatable :: eps_ij       ! the combining parameter epsilon, size [n_comp,n_comp]
    real(dp), dimension(:, :), allocatable :: eps_asoc_ij  ! the combining associating parameter epsilon, size [n_comp,n_comp]
    real(dp), dimension(:, :), allocatable :: kapa_asoc_ij ! the combining associating parameter kappa, size [n_comp,n_comp]
    real(dp) :: m2es3, m2e2s3 ! the abbreviated terms (m**2*eps*sigma**3), (m**2*eps**2*sigma**3) for z_disp computation
    real(dp), dimension(:), allocatable :: m2es3_dx, m2e2s3_dx ! the abbreviated terms first composition derivative, size [n_comp]
    real(dp), dimension(:, :), allocatable :: ees3s3_s3nnmm, ees5s5_s7nnqq, ees3s5_s5nq ! the abreviated polar term used in a2 computation (both dipolar[mm] and quadrupolar[qq]), size [n_comp,n_comp]
    real(dp), dimension(:, :, :), allocatable :: eees3s3s3_sssnnnmmm, eees5s5s5_s3s3s3nnnqqq, eees4s4s4_s2s2s2nnqnqq ! the abreviated polar term used in the a3 composition (both dipolar[mmm] and quadrupolar[qqq]), size [n_comp,n_comp,n_comp]

    real(dp), dimension(:, :), allocatable :: m_ij ! the matrix of segment parameters, size [n_comp,n_comp]
    real(dp), dimension(:, :, :), allocatable :: m_ijk ! the matrix of segment parameters, size [n_comp,n_comp,n_comp]

contains
    !the subroutines and functions
    ! ========================= ALLOCATION/ DEALLOCATION METHODS =============================
    ! ------------- ALLOCATION ---------------
    subroutine allocate_contrib_sub()
        ! allocate all related variables for controll module
        ! NOTE expect that n_comp is known nonzero
        write (stddeb,*) "=> allocate_contrib_sub" !DEBUG
        if (n_comp == 0) then
            write (stderr, *) 'Error in function allocate_contrib n_comp==0'
            stop
        end if

        allocate(a__disp_dx(7, n_comp))
        allocate(b__disp_dx(7, n_comp))

        allocate(a__dd(5, n_comp, n_comp))
        allocate(b__dd(5, n_comp, n_comp))
        allocate(c__dd(5, n_comp, n_comp, n_comp))

        allocate(a__qq(5, n_comp, n_comp))
        allocate(b__qq(5, n_comp, n_comp))
        allocate(c__qq(5, n_comp, n_comp, n_comp))

        allocate(a__dq(5, n_comp, n_comp))
        allocate(b__dq(5, n_comp, n_comp))
        allocate(c__dq(5, n_comp, n_comp, n_comp))

        allocate(dzeta_dx(0:3, 1:n_comp))
        allocate(d_seg(n_comp))
        allocate(d_seg_dtt(n_comp))

        allocate(sigma_ij(n_comp, n_comp))
        allocate(eps_ij(n_comp, n_comp))
        allocate(eps_asoc_ij(n_comp, n_comp))
        allocate(kapa_asoc_ij(n_comp, n_comp))

        allocate(m2es3_dx(n_comp))
        allocate(m2e2s3_dx(n_comp))

        allocate(ees3s3_s3nnmm(n_comp, n_comp))
        allocate(ees5s5_s7nnqq(n_comp, n_comp))
        allocate(ees3s5_s5nq(n_comp, n_comp))

        allocate(eees3s3s3_sssnnnmmm(n_comp, n_comp, n_comp))
        allocate(eees5s5s5_s3s3s3nnnqqq(n_comp, n_comp, n_comp))
        allocate(eees4s4s4_s2s2s2nnqnqq(n_comp, n_comp, n_comp))

        allocate(m_ij(n_comp, n_comp))
        allocate(m_ijk(n_comp, n_comp, n_comp))

        write (stddeb,*) "<= allocate_contrib_sub" !DEBUG
    end subroutine allocate_contrib_sub

    ! ------------- DEALLOCATION ---------------
    subroutine deallocate_contrib_sub()
        ! deallocate all related variables for controll module
        ! performs check if memory is allocated otherwise does nothing to unallocated memory
        write (stddeb,*) "=> deallocate_contrib_sub" !DEBUG

        if (allocated(a__disp_dx) .eqv. .true.) then 
            deallocate(a__disp_dx)
        end if
        if (allocated(b__disp_dx) .eqv. .true.) then 
            deallocate(b__disp_dx)
        end if
        
        if (allocated(a__dd) .eqv. .true.) then 
            deallocate(a__dd)
        end if
        if (allocated(b__dd) .eqv. .true.) then 
            deallocate(b__dd)
        end if
        if (allocated(c__dd) .eqv. .true.) then 
            deallocate(c__dd)
        end if
        
        if (allocated(a__qq) .eqv. .true.) then 
            deallocate(a__qq)
        end if
        if (allocated(b__qq) .eqv. .true.) then 
            deallocate(b__qq)
        end if
        if (allocated(c__qq) .eqv. .true.) then 
            deallocate(c__qq)
        end if
        
        if (allocated(a__dq) .eqv. .true.) then 
            deallocate(a__dq)
        end if
        if (allocated(b__dq) .eqv. .true.) then 
            deallocate(b__dq)
        end if
        if (allocated(c__dq) .eqv. .true.) then 
            deallocate(c__dq)
        end if
        
        if (allocated(dzeta_dx) .eqv. .true.) then 
            deallocate(dzeta_dx)
        end if
        if (allocated(d_seg) .eqv. .true.) then 
            deallocate(d_seg)
        end if
        if (allocated(d_seg_dtt) .eqv. .true.) then 
            deallocate(d_seg_dtt)
        end if
        
        if (allocated(sigma_ij) .eqv. .true.) then 
            deallocate(sigma_ij)
        end if
        if (allocated(eps_ij) .eqv. .true.) then 
            deallocate(eps_ij)
        end if
        if (allocated(eps_asoc_ij) .eqv. .true.) then 
            deallocate(eps_asoc_ij)
        end if
        if (allocated(kapa_asoc_ij) .eqv. .true.) then 
            deallocate(kapa_asoc_ij)
        end if
        
        if (allocated(m2es3_dx) .eqv. .true.) then 
            deallocate(m2es3_dx)
        end if
        if (allocated(m2e2s3_dx) .eqv. .true.) then 
            deallocate(m2e2s3_dx)
        end if
        
        if (allocated(ees3s3_s3nnmm) .eqv. .true.) then 
            deallocate(ees3s3_s3nnmm)
        end if
        if (allocated(ees5s5_s7nnqq) .eqv. .true.) then 
            deallocate(ees5s5_s7nnqq)
        end if
        if (allocated(ees3s5_s5nq) .eqv. .true.) then 
            deallocate(ees3s5_s5nq)
        end if
        if (allocated(eees3s3s3_sssnnnmmm) .eqv. .true.) then 
            deallocate(eees3s3s3_sssnnnmmm)
        end if
        if (allocated(eees5s5s5_s3s3s3nnnqqq) .eqv. .true.) then 
            deallocate(eees5s5s5_s3s3s3nnnqqq)
        end if
        if (allocated(eees4s4s4_s2s2s2nnqnqq) .eqv. .true.) then 
            deallocate(eees4s4s4_s2s2s2nnqnqq)
        end if
        
        if (allocated(m_ij) .eqv. .true.) then 
            deallocate(m_ij)
        end if
        if (allocated(m_ijk) .eqv. .true.) then 
            deallocate(m_ijk)
        end if
      
        ! this option ensures that all dependent function will be computed again
        init_flag(1) = .false.
        init_flag(2) = .false.
        init_flag(3) = .false.
        init_flag(4) = .false.
        init_flag(5) = .false.
        init_flag(6) = .false.

        write (stddeb,*) "<= deallocate_contrib_sub" !DEBUG
    end subroutine deallocate_contrib_sub
    ! ============================= DISPERSIVE CONTRIBUTION ==================================
    logical function initialize_a0b0_disp_fun() result(success)
        ! set up the dispersion parameters
        ! DEPENDENCY: NONE
        implicit none
        ! variable section
        ! OLD logical :: success
        write (stddeb,*) "=> initialize_a0b0_disp_fun" !DEBUG
        if (init_flag(1) .eqv. .false.) then !to prevent unnecessary initialization
            ! ! initialization of the a0 dispersion approximation parameter
            ! a0_disp(1, :) = (/   0.91056314451539_dp, -0.30840169182720_dp, -0.09061483509767_dp/)
            ! a0_disp(2, :) = (/   0.63612814494991_dp, 0.18605311591713_dp, 0.45278428063920_dp/)
            ! a0_disp(3, :) = (/   2.68613478913903_dp, -2.50300472586548_dp, 0.59627007280101_dp/)
            ! a0_disp(4, :) = (/ -26.5473624914884_dp, 21.4197936296668_dp, -1.72418291311787_dp/)
            ! a0_disp(5, :) = (/  97.7592087835073_dp, -65.2558853303492_dp, -4.13021125311661_dp/)
            ! a0_disp(6, :) = (/-159.591540865600_dp, 83.3186804808856_dp, 13.7766318697211_dp/)
            ! a0_disp(7, :) = (/  91.2977740839123_dp, -33.7469229297323_dp, -8.67284703679646_dp/)
            ! ! initialization of the b0 dispersion approximation parameter
            ! b0_disp(1, :) = (/   0.72409469413165_dp, -0.57554980753450_dp, 0.09768831158356_dp/)
            ! b0_disp(2, :) = (/   2.23827918609380_dp, 0.69950955214436_dp, -0.25575749816100_dp/)
            ! b0_disp(3, :) = (/  -4.002584948463420_dp, 3.892567338953070_dp, -9.155856152973211_dp/)
            ! b0_disp(4, :) = (/ -21.003576814846479_dp, -17.215471647772119_dp, 20.642075974397240_dp/)
            ! b0_disp(5, :) = (/  26.855641362661501_dp, 192.6722644652495_dp, -38.804430052062848_dp/)
            ! b0_disp(6, :) = (/ 206.5513384066188_dp, -161.8264616487648_dp, 93.626774077014602_dp/)
            ! b0_disp(7, :) = (/-355.6023561220795_dp, -165.2076934555607_dp, -29.66690558514725_dp/)

            ! LITERATURE COEFITIENTS
            ! initialization of the a0 dispersion approximation parameter
            a0_disp(1, :) = (/0.9105631445_dp, -0.3084016918_dp, -0.0906148351_dp /)
            a0_disp(2, :) = (/0.6361281449_dp, 0.1860531159_dp, 0.4527842806_dp /)
            a0_disp(3, :) = (/2.6861347891_dp, -2.5030047259_dp, 0.5962700728_dp /)
            a0_disp(4, :) = (/-26.547362491_dp, 21.419793629_dp, -1.7241829131_dp /)
            a0_disp(5, :) = (/97.759208784_dp, -65.255885330_dp, -4.1302112531_dp /)
            a0_disp(6, :) = (/-159.59154087_dp, 83.318680481_dp, 13.776631870_dp /)
            a0_disp(7, :) = (/91.297774084_dp, -33.746922930_dp, -8.6728470368_dp /)
            ! initialization of the b0 dispersion approximation parameter
            b0_disp(1, :) = (/0.7240946941_dp, -0.5755498075_dp, 0.0976883116_dp /)
            b0_disp(2, :) = (/2.2382791861_dp, 0.6995095521_dp, -0.2557574982_dp /)
            b0_disp(3, :) = (/-4.0025849485_dp, 3.8925673390_dp, -9.1558561530_dp /)
            b0_disp(4, :) = (/-21.003576815_dp, -17.215471648_dp, 20.642075974_dp /)
            b0_disp(5, :) = (/26.855641363_dp, 192.67226447_dp, -38.804430052_dp /)
            b0_disp(6, :) = (/206.55133841_dp, -161.82646165_dp, 93.626774077_dp /)
            b0_disp(7, :) = (/-355.60235612_dp, -165.20769346_dp, -29.666905585_dp /)

            ! OLD initialization of the b0 dispersion approximation parameter according to VV code
            ! 				b0_disp(1,:) = (/  0.72409469413165_dp, -0.57554980753450_dp,  0.09768831158356_dp/)
            ! 				b0_disp(2,:) = (/  1.11913959304690_dp,  0.34975477607218_dp, -0.12787874908050_dp/) *2
            ! 				b0_disp(3,:) = (/ -1.33419498282114_dp,  1.29752244631769_dp, -3.05195205099107_dp/) *3
            ! 				b0_disp(4,:) = (/ -5.25089420371162_dp, -4.30386791194303_dp,  5.16051899359931_dp/) *4
            ! 				b0_disp(5,:) = (/  5.37112827253230_dp, 38.5344528930499_dp,  -7.76088601041257_dp/) *5
            ! 				b0_disp(6,:) = (/ 34.4252230677698_dp, -26.9710769414608_dp,  15.6044623461691_dp/) *6
            ! 				b0_disp(7,:) = (/-50.8003365888685_dp, -23.6010990650801_dp,  -4.23812936930675_dp/) *7

            init_flag(1) = .true.
            write (stddeb,*) "a0b0_disp is initialized " !DEBUG
        end if
        success = .true.
        write (stddeb,*) "<= initialize_a0b0_disp_fun" !DEBUG
    end function initialize_a0b0_disp_fun

    logical function initialize_dzeta_fun() result(success)
        ! set up the dzeta vector without the rho multiplication and m_mean value
        ! DEPENDENCY: Temperature, composition
        ! variable section
        integer :: j ! counters
        ! OLD logical :: success
        write (stddeb,*) "=> initialize_dzeta_fun" !DEBUG
        if (init_flag(2) .eqv. .false.) then !to prevent unnecessary initialization
            dzeta = (/0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/) ! initialization
            do j = 1, n_comp
                ! segment diameter di computation
                d_seg(j) = param_sigma(j) * (1.0_dp - 0.12_dp * exp(-3.0_dp * param_eps_kb(j) / tt))
                ! 					write (*,"('param_sigma(',i2,'): ',f12.8)") j, param_sigma(j) ! DEBUG test the computation
                ! 					write (*,"('param_eps(',i2,'): ',f12.8)") j, param_eps_kb(j) ! DEBUG test the computation
                ! 					write (*,"('d_seg(',i2,'): ',f12.8)") j, d_seg(j) ! DEBUG test the computation
                ! dzeta computation
                dzeta(0) = dzeta(0) + mol_rat(j) * param_m(j)
                dzeta(1) = dzeta(1) + mol_rat(j) * param_m(j) * d_seg(j)
                dzeta(2) = dzeta(2) + mol_rat(j) * param_m(j) * d_seg(j)**2
                dzeta(3) = dzeta(3) + mol_rat(j) * param_m(j) * d_seg(j)**3
            end do

            ! the derivatives for the d/dxi calculations
            do j = 1, n_comp
                ! 						dzeta_dx(:,j) = (/0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/) ! initialization
                dzeta_dx(0, j) = param_m(j)
                dzeta_dx(1, j) = param_m(j) * d_seg(j)
                dzeta_dx(2, j) = param_m(j) * d_seg(j)**2
                dzeta_dx(3, j) = param_m(j) * d_seg(j)**3
                ! 						write (*,"('param_m(',i2,'): ',f12.8)") j, param_sigma(j) ! DEBUG test the computation
            end do
            dzeta_dx = (pi / 6.0_dp) * dzeta_dx

            ! the derivatives for the d/dT calculations
            dzeta_dtt = (/0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp/) ! initialization
            do j = 1, n_comp
                ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.52)
                d_seg_dtt(j) = -param_sigma(j) * 0.12_dp * exp(-3.0_dp * param_eps_kb(j) / tt) &
                        & * (3.0_dp * param_eps_kb(j) / tt**2) ! the ddi/dT derivative
                ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.53)
                dzeta_dtt(1) = dzeta_dtt(1) + mol_rat(j) * param_m(j) * d_seg_dtt(j)
                dzeta_dtt(2) = dzeta_dtt(2) + mol_rat(j) * param_m(j) * 2.0_dp * d_seg(j) * d_seg_dtt(j)
                dzeta_dtt(3) = dzeta_dtt(3) + mol_rat(j) * param_m(j) * 3.0_dp * d_seg(j)**2 * d_seg_dtt(j)
            end do
            dzeta_dtt = (pi / 6.0_dp) * dzeta_dtt

            ! the mean number of segments computation and actual dzeta
            m_mean = dzeta(0)
            dzeta = (pi / 6.0_dp) * dzeta

            write (stddeb,*) "dzeta is initialized " !DEBUG
            init_flag(2) = .true.
        end if
        success = .true.
        write (stddeb,*) "<= initialize_dzeta_fun" !DEBUG
    end function initialize_dzeta_fun

    subroutine check_dzeta_sub(dzeta_r)
        ! perform the check of dzeta validity before computation proceeds
        ! intended for use in subroutines that may easily crash with incorrect dzeta_r
        ! method ensure that run does not lead to unexplained crash because of unphysical state
        implicit none
        ! variable section
        real(dp), intent(in), dimension(0:3) :: dzeta_r ! the dzeta vector dzeta_n for n=0,1,2,3 * rho !!

        write (stddeb,*) "=> check_dzeta_sub" !DEBUG
        ! look-ahead control of correct system setup
        if (1.0_dp - dzeta_r(3)<=0.0) then
            write (stderr, *) 'Error in 1.0-dzeta_r(3) nonpositive value lead to unphysical state, check the system setup.'
            stop
        end if
        if (dzeta_r(3)==0.0) then
            write (stderr, *) 'Warning dzeta_r(3) 0.0 value may lead to unphysical state.'
        end if

        write (stddeb,*) "<= check_dzeta_sub" !DEBUG
    end subroutine check_dzeta_sub

    logical function initialize_combrules_fun() result(success)
        ! set up the parameters sigma_ij eps_ij for pairs of unlike segments
        !  additionally sets up the abbreviated terms (m**2*eps*sigma**3), (m**2*eps**2*sigma**3)
        ! TODO decide if combrules will be private and called from dzeta function
        !      - if yes then change the implementation accordingly
        ! DEPENDENCY: temperature, composition
        implicit none
        ! variable section
        integer :: i, j ! counters
        real(dp) :: xxmm ! temporary variable for xi*xj*mi*mj
        success = .false.

        write (stddeb,*) "=> initialize_combrules_fun" !DEBUG
        if (init_flag(3) .eqv. .false.) then  !to prevent unnecessary initialization

            ! initialization
            m2es3 = 0.0_dp
            m2e2s3 = 0.0_dp

            ! DEBUG SECTION
            ! write(*,*) param_m
            ! write(*,*) param_sigma
            ! write(*,*) param_eps_kb
            ! write(*,*) param_kap_asoc
            ! write(*,*) param_eps_asoc_kb
            ! write(*,*) mol_rat
            ! write(*,*) k_ij
            
            if (abs(k_ij(1,1)) .gt. 100 ) then ! BUG TODO - resolve proper initialization
                call set_k_ij(0.0_dp)
            end if

            ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.12 - A.15)
            ! based upon paper Gross-Sadiwski 2002 application of perturbed-chain SAFT eos to associating systems eq (2,3)
            do i = 1, n_comp
                sigma_ij(i, i) = param_sigma(i) ! A.14 diagonal term
                eps_ij(i, i) = param_eps_kb(i) * (1.0_dp - k_ij(i, i)) ! A.15 diagonal term
                eps_asoc_ij(i, i) = param_eps_asoc_kb(i) ! eq (2) diagonal term
                kapa_asoc_ij(i, i) = param_kap_asoc(i) ! eq (3) diagonal term

                m2es3 = m2es3 + (eps_ij(i, i) / tt) * sigma_ij(i, i) * (sigma_ij(i, i) * mol_rat(i) * param_m(i))**2 ! A.12 diagonal term
                m2e2s3 = m2e2s3 + sigma_ij(i, i) * ((eps_ij(i, i) / tt) * sigma_ij(i, i) * mol_rat(i) * param_m(i))**2 ! A.13 diagonal term

                do j = i + 1, n_comp
                    sigma_ij(i, j) = (param_sigma(i) + param_sigma(j)) / 2.0_dp ! A.14 nondiagonal
                    sigma_ij(j, i) = sigma_ij(i, j) ! use of symmetry
                    eps_ij(i, j) = sqrt(param_eps_kb(i) * param_eps_kb(j)) * (1.0_dp - k_ij(i, j)) ! A.15 nondiagonal
                    eps_ij(j, i) = eps_ij(i, j) ! use of symmetry

                    eps_asoc_ij(i, j) = (param_eps_asoc_kb(i) + param_eps_asoc_kb(j)) / 2.0_dp ! eq (2) nondiagonal term
                    eps_asoc_ij(j, i) = eps_asoc_ij(i, j) ! use of symmetry
                    kapa_asoc_ij(i, j) = sqrt(param_kap_asoc(i) * param_kap_asoc(j) * param_sigma(i) * param_sigma(j)) &
                            & * param_sigma(i) * param_sigma(j) / sigma_ij(i, j)**3 ! eq (3) nondiagonal term
                    kapa_asoc_ij(j, i) = kapa_asoc_ij(i, j)! use of symmetry

                    xxmm = (param_m(i) * param_m(j) * mol_rat(i) * mol_rat(j)) ! temporary variable
                    m2es3 = m2es3 + 2.0_dp * xxmm * (eps_ij(i, j) / tt) * sigma_ij(i, j)**3 ! A.12 nodiagonal 2* is for symmetrical term
                    m2e2s3 = m2e2s3 + 2.0_dp * xxmm * (eps_ij(i, j) / tt)**2 * sigma_ij(i, j)**3 ! A.13 nodiagonal 2* is for symmetrical term
                end do
            end do
            ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.39, A.40)
            
            do i = 1, n_comp
                m2es3_dx(i) = 0.0_dp
                m2e2s3_dx(i) = 0.0_dp
                do j = 1, n_comp
                    m2es3_dx(i) = m2es3_dx(i) + mol_rat(j) * param_m(j) * (eps_ij(i, j) / tt) * sigma_ij(i, j)**3
                    m2e2s3_dx(i) = m2e2s3_dx(i) + mol_rat(j) * param_m(j) * (eps_ij(i, j) / tt)**2 * sigma_ij(i, j)**3
                end do
                m2es3_dx(i) = 2.0_dp * param_m(i) * m2es3_dx(i)
                m2e2s3_dx(i) = 2.0_dp * param_m(i) * m2e2s3_dx(i)
            end do        
            ! DEBUG SECTION START
            write(stddeb,"('sigma_ij     : ')", advance="no")
            write(stddeb,*) sigma_ij
            write(stddeb,"('eps_ij       : ')", advance="no")
            write(stddeb,*) eps_ij
            write(stddeb,"('eps_asoc_ij  : ')", advance="no")
            write(stddeb,*) eps_asoc_ij
            write(stddeb,"('kapa_asoc_ij : ')", advance="no")
            write(stddeb,*) kapa_asoc_ij
            write(stddeb,"('xxmm         : ')", advance="no")
            write(stddeb,*) xxmm
            write(stddeb,"('m2es3_dx     : ')", advance="no")
            write(stddeb,*) m2es3_dx
            write(stddeb,"('m2e2s3_dx    : ')", advance="no")
            write(stddeb,*) m2e2s3_dx
            ! DEBUG SECTION STOP

            write(stddeb,*) "combrules sigma_ij, eps_ij, eps_asoc_ij, kapa_asoc_ij, xxmm, m2es3_dx, m2e2s3_dx are initialized " !DEBUG
            init_flag(3) = .true.
        end if
        success = .true.
        write (stddeb,*) "<= initialize_combrules_fun" !DEBUG
    end function initialize_combrules_fun

    logical function compute_ab_disp_fun()    result(success)
        ! the function that computes the a,b dispersion parameters
        !  function also computes the derivatives of parameters and stores
        !  them in the module variables
        ! DEPENDENCY: none
        ! INDDIRECT DEPENDENCY: composition (m_mean)
        implicit none
        ! variable section
        integer :: i, j ! counters
        write (stddeb,*) "=> compute_ab_disp_fun" !DEBUG

        if (init_flag(4) .eqv. .false.) then  !to prevent unnecessary initialization

            ! compute the a,b parameters
            ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.18, A.19)
            do i = 1, 7
                ! A.18: a_disp = a01 + (m-1)/m *a02 + ((m-1)(m-2))/(m**2) *a03
                a__disp(i) = a0_disp(i, 1) + (1.0_dp - 1.0_dp / m_mean) * a0_disp(i, 2) &
                        & + (1.0_dp - 1.0_dp / m_mean) * (1.0_dp - 2.0_dp / m_mean) * a0_disp(i, 3)
                ! A.19: b_disp = b01 + (m-1)/m *b02 + ((m-1)(m-2))/(m**2) *b03
                b__disp(i) = b0_disp(i, 1) + (1.0_dp - 1.0_dp / m_mean) * b0_disp(i, 2) + &
                        & (1.0_dp - 1.0_dp / m_mean) * (1.0_dp - 2.0_dp / m_mean) * b0_disp(i, 3)

                ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.44, A.45)
                do j = 1, n_comp
                    !  compute composition derivative da/dxi = mi*a02/m**2 + 3*mi*a03/m**2 -4*mi*a03/m**3
                    a__disp_dx(i, j) = (param_m(j) * a0_disp(i, 2)) / (m_mean**2) + &
                            & (3.0_dp * param_m(j) * a0_disp(i, 3)) / (m_mean**2) - &
                            & (4.0_dp * param_m(j) * a0_disp(i, 3)) / (m_mean**3)
                    !  compute composition derivative db/dxi = mi*b02/m**2 + 3*mi*b03/m**2 -4*mi*b03/m**3
                    b__disp_dx(i, j) = (param_m(j) * b0_disp(i, 2)) / (m_mean**2) + &
                            & (3.0_dp * param_m(j) * b0_disp(i, 3)) / (m_mean**2) - &
                            & (4.0_dp * param_m(j) * b0_disp(i, 3)) / (m_mean**3)
                end do            
            end do

            write(stddeb,*) "combrules a,b dispersion parameters are initialized " !DEBUG
            init_flag(4) = .true.
        end if

        success = .true.
        write (stddeb,*) "<= compute_ab_disp_fun" !DEBUG
    end function compute_ab_disp_fun

    ! =============================  POLAR CONTRIBUTION  ===============================
    logical function initialize_a0b0c0_polar_fun() result(success)
        !! set up all polar parameters according to the tabelated values
        ! DEPENDENCY: NONE
        implicit none
        ! variable section
        ! OLD logical :: success

        write (stddeb,*) "=> initialize_a0b0c0_polar_fun" !DEBUG
        if (init_flag(5) .eqv. .false.) then  !to prevent unnecessary initialization
            !----------------- DIPOLAR-DIPOLAR -----------------------

            ! initialization of the a0 dd polar approximation parameter
            a0_dd(1, :) = (/  0.3043504_dp, 0.9534641_dp, -1.1610080_dp/)
            a0_dd(2, :) = (/ -0.1358588_dp, -1.8396383_dp, 4.5258607_dp/)
            a0_dd(3, :) = (/  1.4493329_dp, 2.0131180_dp, 0.9751222_dp/)
            a0_dd(4, :) = (/  0.3556977_dp, -7.3724958_dp, -12.2810380_dp/)
            a0_dd(5, :) = (/ -2.0653308_dp, 8.2374135_dp, 5.9397575_dp/)
            ! initialization of the b0 dd polar approximation parameter
            b0_dd(1, :) = (/  0.2187939_dp, -0.5873164_dp, 3.4869576_dp/)
            b0_dd(2, :) = (/ -1.1896431_dp, 1.2489132_dp, -14.9159740_dp/)
            b0_dd(3, :) = (/  1.1626889_dp, -0.5085280_dp, 15.3720220_dp/)
            b0_dd(4, :) = (/  0.0000000_dp, 0.0000000_dp, 0.0000000_dp/)
            b0_dd(5, :) = (/  0.0000000_dp, 0.0000000_dp, 0.0000000_dp/)
            ! initialization of the c0 dd polar approximation parameter
            c0_dd(1, :) = (/ -0.0646774_dp, -0.9520876_dp, -0.6260979_dp/)
            c0_dd(2, :) = (/  0.1975882_dp, 2.9924258_dp, 1.2924686_dp/)
            c0_dd(3, :) = (/ -0.8087562_dp, -2.3802636_dp, 1.6542783_dp/)
            c0_dd(4, :) = (/  0.6902849_dp, -0.2701261_dp, -3.4396744_dp/)
            c0_dd(5, :) = (/  0.0000000_dp, 0.0000000_dp, 0.0000000_dp/)

            !-------------- QUADRUPOLAR-QUADRUPOLAR--------------------

            ! initialization of the a0 qq polar approximation parameter
            a0_qq(1, :) = (/  1.2378308_dp, 1.2854109_dp, 1.7942954_dp/)
            a0_qq(2, :) = (/  2.4355031_dp, -11.4656150_dp, 0.7695103_dp/)
            a0_qq(3, :) = (/  1.6330905_dp, 22.0868930_dp, 7.2647923_dp/)
            a0_qq(4, :) = (/ -1.6118152_dp, 7.4691383_dp, 94.4866990_dp/)
            a0_qq(5, :) = (/  6.9771185_dp, -17.1977720_dp, -77.1484580_dp/)
            ! initialization of the b0 qq polar approximation parameter
            b0_qq(1, :) = (/  0.4542718_dp, -0.8137340_dp, 6.8682675_dp/)
            b0_qq(2, :) = (/ -4.5016264_dp, 10.0640300_dp, -5.1732238_dp/)
            b0_qq(3, :) = (/  3.5858868_dp, -10.8766310_dp, -17.2402070_dp/)
            b0_qq(4, :) = (/  0.0000000_dp, 0.0000000_dp, 0.0000000_dp/)
            b0_qq(5, :) = (/  0.0000000_dp, 0.0000000_dp, 0.0000000_dp/)
            ! initialization of the c0 qq polar approximation parameter
            c0_qq(1, :) = (/  -0.5000437_dp, 2.0002094_dp, 3.1358271_dp/)
            c0_qq(2, :) = (/   6.5318692_dp, -6.7838658_dp, 7.2475888_dp/)
            c0_qq(3, :) = (/ -16.0147800_dp, 20.3832460_dp, 3.0759478_dp/)
            c0_qq(4, :) = (/  14.4259700_dp, -10.8959840_dp, 0.0000000_dp/)
            c0_qq(5, :) = (/   0.0000000_dp, 0.0000000_dp, 0.0000000_dp/)

            !--------------- DIPOLAR-QUADRUPOLAR ----------------------

            ! initialization of the a0 dq polar approximation parameter
            a0_dq(1, :) = (/  0.6970950_dp, -0.6734593_dp, 0.6703408_dp/)
            a0_dq(2, :) = (/ -0.6335541_dp, -1.4258991_dp, -4.3384718_dp/)
            a0_dq(3, :) = (/  2.9455090_dp, 4.1944139_dp, 7.2341684_dp/)
            a0_dq(4, :) = (/ -1.4670273_dp, 1.0266216_dp, 0.0000000_dp/)
            a0_dq(5, :) = (/  0.0000000_dp, 0.0000000_dp, 0.0000000_dp/)
            ! initialization of the b0 dq polar approximation parameter
            b0_dq(1, :) = (/ -0.4840383_dp, 0.6765101_dp, -1.1675601_dp/)
            b0_dq(2, :) = (/  1.9704055_dp, -3.0138675_dp, 2.1348843_dp/)
            b0_dq(3, :) = (/ -2.1185727_dp, 0.4674266_dp, 0.0000000_dp/)
            b0_dq(4, :) = (/  0.0000000_dp, 0.0000000_dp, 0.0000000_dp/)
            b0_dq(5, :) = (/  0.0000000_dp, 0.0000000_dp, 0.0000000_dp/)
            ! initialization of the c0 dq polar approximation parameter
            c0_dq(1, :) = (/  7.8464310_dp, -20.7220200_dp/)
            c0_dq(2, :) = (/ 33.4270000_dp, -58.6390400_dp/)
            c0_dq(3, :) = (/  4.6891110_dp, -1.7648870_dp/)
            c0_dq(4, :) = (/  0.0000000_dp, 0.0000000_dp/)
            c0_dq(5, :) = (/  0.0000000_dp, 0.0000000_dp/)

            write(stddeb,*) "combrules a,b,c dipolar,quadrupolar and dq parameters are initialized " !DEBUG
            init_flag(5) = .true.
        end if

        success = .true.
        write (stddeb,*) "<= initialize_a0b0c0_polar_fun" !DEBUG
    end function initialize_a0b0c0_polar_fun

    logical function initialize_polrules_fun() result(success)
        ! the function that computes the set up the helper variables for polar calculation
        ! the helper variable names are abbreviation of terms condensed (e=epsilon/k_b, s=sigma, n=mu_dimles, q=qq_dimles)
        ! function also call the a__dd/qq/dq, ... computation method and set the abc dd/qq/dq parameters
        ! DEPENDENCY none
        use control_mod, only : is_dd, is_qq, alpha_dq, init_flag
        use param_list_mod, only : param_n_u, param_u, param_n_q, param_q
        implicit none
        ! variable section
        real(dp), dimension(:), allocatable :: mu_dimles, qq_dimles ! dimensionless dipolar moment an quadrupolar moment
        ! NOTE do not confuse with global sigma_ij
        real(dp) :: sigm_ij, sigm_jk, sigm_ki    ! the nondiagonal sigma parameters
        real(dp), parameter :: third = 1.0_dp / 3.0_dp ! for the 1/3 specification
        real(dp), dimension(:), allocatable :: ts_2clj ! helper variable for conversion between TS and 2CLJ formulation = 4/m**2
        integer :: i, j, k ! counters

        write (stddeb,*) "=> initialize_polrules_fun" !DEBUG
        ! ensure that a,b,c parameters are initialized but not more times than necesarry
        if (init_flag(5) .eqv. .false.) then
            ! 		write (*,*) 'Call the initialize_a0b0_disp_fun from z_dd_sub' ! DEBUG control of computation
            if (initialize_a0b0c0_polar_fun().eqv..false.) then
                write (stderr, *) 'Error in helper function initialize_a0b0c0_polar_fun'
                stop
            end if
        end if

        if ( (init_flag(6) .eqv. .false.) .and. ((is_dd .eqv. .true.) .or. (is_qq .eqv. .true.))) then  !to prevent unnecessary initialization in special cases
            ! + reduce initialization further when no dipolar or quadrupolar interaction occurs
            ! - due to the excesive complication dd/qq/dq are computed in one place which is wasting in some cases
            ! ----- allocation section ----
            allocate(mu_dimles(n_comp), qq_dimles(n_comp))
            allocate(ts_2clj(n_comp))

            ! based upon paper Gross-Vrabec 2006 application of text section after eq 9
            ! NOTE 1.0E19 is adopted from Gross 2005 eq. A14 - this information is missing in Gross-Vrabec 2006
            ! BUG possible wrong scaling because in paper is sigma**5 whereas here is sigma**3 -> applying same reasoning gives 1.0E39
            mu_dimles(:) = param_u(:)**2 / (param_m(:) * param_eps_kb(:) * param_sigma(:)**3) / (k_b * 1.0E19)
            ! based upon paper Gross 2005 eq. A14
            qq_dimles(:) = param_q(:)**2 / (param_m(:) * param_eps_kb(:) * param_sigma(:)**5) / (k_b * 1.0E19)
            ! BUG possible wrong scaling because in paper for dipol-quadrupol interaction is no scaling which contradicts the compatibility between dipolar and quadrupolar interaction

            ! 				write (*,"('param_m: ',e)")param_m ! DEBUG check that local changes remain local
            ! 				write (*,"('mu_dimles: ',e,e)")mu_dimles(1),mu_dimles(2) ! DEBUG printout of array element
            ! 				write (*,"('qq_dimles: ',e,e)")qq_dimles(1),qq_dimles(2) ! DEBUG printout of array element

            ! initialization
            ! based upon paper Vrabec-Gross 2007 application of table 1
            ts_2clj(:) = 4.0_dp / param_m(:)**2

            ! based upon paper Gross-Vrabec 2006 application of eq (8, 9)
            ees3s3_s3nnmm = 0.0_dp
            eees3s3s3_sssnnnmmm = 0.0_dp
            ! based upon paper Gross 2005 application of eq (9, 10)
            ees5s5_s7nnqq = 0.0_dp
            eees5s5s5_s3s3s3nnnqqq = 0.0_dp
            ! based upon paper Vrabec-Gross 2007 application of eq (14, 15)
            ees3s5_s5nq = 0.0_dp
            eees4s4s4_s2s2s2nnqnqq = 0.0_dp


            ! ees3s3_s3nnmm and eees3s3s3_sssnnnmmm computation and others corresponding to polar interaction
            ! evaluation of term depending on the input parameters
            ! computation use the symmetry and diagonal simplification in case dd/qq
            ! NOTE dq is not symmetrical the cycle holds but multiple elements are computed to account for that
            do i = 1, n_comp
                ! === DD interaction ===
                ees3s3_s3nnmm(i, i) = param_eps_kb(i)**2 * param_sigma(i)**3 * param_n_u(i)**2 * mu_dimles(i)**2
                eees3s3s3_sssnnnmmm(i, i, i) = param_eps_kb(i)**3 * param_sigma(i)**6 * param_n_u(i)**3 * mu_dimles(i)**3
                ! === QQ interaction ===
                ees5s5_s7nnqq(i, i) = param_eps_kb(i)**2 * param_sigma(i)**3 * param_n_q(i)**2 * qq_dimles(i)**2
                eees5s5s5_s3s3s3nnnqqq(i, i, i) = param_eps_kb(i)**3 * param_sigma(i)**6 * param_n_q(i)**3 * qq_dimles(i)**3
                ! === DQ interaction ===
                ees3s5_s5nq(i, i) = param_eps_kb(i)**2 * param_sigma(i)**3 * mu_dimles(i) * qq_dimles(i) ! NOTE ts_2clj cancel each out
                eees4s4s4_s2s2s2nnqnqq(i, i, i) = param_eps_kb(i)**3 * param_sigma(i)**6 &
                        & * (mu_dimles(i)**2 * qq_dimles(i) &
                                + alpha_dq * mu_dimles(i) * qq_dimles(i)**2) ! NOTE ts_2clj cancel each out

                ! fist set up the diagonal elements with all indexes same
                call calculate_abc_polar_sub(i, 0, 0, 1, param_m(i), param_m(i), 0.0_dp)

                do j = i + 1, n_comp
                    sigm_ij = (param_sigma(i) + param_sigma(j)) / 2.0_dp
                    ! === DD interaction ===
                    ees3s3_s3nnmm(i, j) = param_eps_kb(i) * param_eps_kb(j) &
                            & * param_sigma(i)**3 * param_sigma(j)**3 / sigm_ij**3 &
                            & * param_n_u(i) * param_n_u(j) &
                            & * mu_dimles(i) * mu_dimles(j)
                    ees3s3_s3nnmm(j, i) = ees3s3_s3nnmm(i, j)

                    eees3s3s3_sssnnnmmm(i, j, j) = param_eps_kb(i) * param_eps_kb(j)**2 &
                            & * param_sigma(i)**3 * param_sigma(j)**5 / sigm_ij**2 &
                            & * param_n_u(i) * param_n_u(j)**2 &
                            & * mu_dimles(i) * mu_dimles(j)**2
                    eees3s3s3_sssnnnmmm(j, i, j) = eees3s3s3_sssnnnmmm(i, j, j)
                    eees3s3s3_sssnnnmmm(j, j, i) = eees3s3s3_sssnnnmmm(i, j, j)

                    eees3s3s3_sssnnnmmm(i, i, j) = param_eps_kb(j) * param_eps_kb(i)**2 &
                            & * param_sigma(j)**3 * param_sigma(i)**5 / sigm_ij**2 &
                            & * param_n_u(j) * param_n_u(i)**2 &
                            & * mu_dimles(j) * mu_dimles(i)**2
                    eees3s3s3_sssnnnmmm(i, j, i) = eees3s3s3_sssnnnmmm(i, i, j)
                    eees3s3s3_sssnnnmmm(j, i, i) = eees3s3s3_sssnnnmmm(i, i, j)
                    ! === QQ interaction ===
                    ees5s5_s7nnqq(i, j) = param_eps_kb(i) * param_eps_kb(j) &
                            & * param_sigma(i)**5 * param_sigma(j)**5 / sigm_ij**7 &
                            & * param_n_q(i) * param_n_q(j) &
                            & * qq_dimles(i) * qq_dimles(j)
                    ees5s5_s7nnqq(j, i) = ees5s5_s7nnqq(i, j)

                    eees5s5s5_s3s3s3nnnqqq(i, j, j) = param_eps_kb(i) * param_eps_kb(j)**2 &
                            & * param_sigma(i)**5 * param_sigma(j)**7 / sigm_ij**6 &
                            & * param_n_q(i) * param_n_q(j)**2 &
                            & * qq_dimles(i) * qq_dimles(j)**2
                    eees5s5s5_s3s3s3nnnqqq(j, i, j) = eees5s5s5_s3s3s3nnnqqq(i, j, j)
                    eees5s5s5_s3s3s3nnnqqq(j, j, i) = eees5s5s5_s3s3s3nnnqqq(i, j, j)

                    eees5s5s5_s3s3s3nnnqqq(i, i, j) = param_eps_kb(j) * param_eps_kb(i)**2 &
                            & * param_sigma(j)**5 * param_sigma(i)**7 / sigm_ij**6 &
                            & * param_n_q(j) * param_n_q(i)**2 &
                            & * qq_dimles(j) * qq_dimles(i)**2
                    eees5s5s5_s3s3s3nnnqqq(i, j, i) = eees5s5s5_s3s3s3nnnqqq(i, i, j)
                    eees5s5s5_s3s3s3nnnqqq(j, i, i) = eees5s5s5_s3s3s3nnnqqq(i, i, j)
                    ! === DQ interaction ===
                    ! NOTE beware no symmetry can be used
                    ees3s5_s5nq(i, j) = param_eps_kb(i) * param_eps_kb(j) &
                            & * param_sigma(i)**3 * param_sigma(j)**5 / sigm_ij**5 &
                            & * mu_dimles(i) * qq_dimles(j) ! NOTE ts_2clj cancel each out
                    ees3s5_s5nq(j, i) = param_eps_kb(j) * param_eps_kb(i) &
                            & * param_sigma(j)**3 * param_sigma(i)**5 / sigm_ij**5 &
                            & * mu_dimles(j) * qq_dimles(i) ! NOTE ts_2clj cancel each out

                    eees4s4s4_s2s2s2nnqnqq(i, j, j) = param_eps_kb(i) * param_eps_kb(j)**2 &
                            & * param_sigma(i)**4 * param_sigma(j)**6 / sigm_ij**4 &
                            & * (mu_dimles(i) * mu_dimles(j) * qq_dimles(j) &
                                    & + alpha_dq * mu_dimles(i) * qq_dimles(j)**2) ! NOTE ts_2clj cancel each out
                    eees4s4s4_s2s2s2nnqnqq(j, i, j) = param_eps_kb(i) * param_eps_kb(j)**2 &
                            & * param_sigma(i)**4 * param_sigma(j)**6 / sigm_ij**4 &
                            & * (mu_dimles(j) * mu_dimles(i) * qq_dimles(j) &
                                    & + alpha_dq * mu_dimles(j) * qq_dimles(i) * qq_dimles(j)) ! NOTE ts_2clj cancel each out
                    eees4s4s4_s2s2s2nnqnqq(j, j, i) = param_eps_kb(i) * param_eps_kb(j)**2 &
                            & * param_sigma(i)**4 * param_sigma(j)**6 / sigm_ij**4 &
                            & * (mu_dimles(j)**2 * qq_dimles(i) &
                                    & + alpha_dq * mu_dimles(j) * qq_dimles(j) * qq_dimles(i)) ! NOTE ts_2clj cancel each out

                    eees4s4s4_s2s2s2nnqnqq(j, i, i) = param_eps_kb(j) * param_eps_kb(i)**2 &
                            & * param_sigma(j)**4 * param_sigma(i)**6 / sigm_ij**4 &
                            & * (mu_dimles(j) * mu_dimles(i) * qq_dimles(i) &
                                    & + alpha_dq * mu_dimles(j) * qq_dimles(i)**2) ! NOTE ts_2clj cancel each out
                    eees4s4s4_s2s2s2nnqnqq(i, j, i) = param_eps_kb(j) * param_eps_kb(i)**2 &
                            & * param_sigma(j)**4 * param_sigma(i)**6 / sigm_ij**4 &
                            & * (mu_dimles(i) * mu_dimles(j) * qq_dimles(i) &
                                    & + alpha_dq * mu_dimles(i) * qq_dimles(j) * qq_dimles(i)) ! NOTE ts_2clj cancel each out
                    eees4s4s4_s2s2s2nnqnqq(i, i, j) = param_eps_kb(j) * param_eps_kb(i)**2 &
                            & * param_sigma(j)**4 * param_sigma(i)**6 / sigm_ij**4 &
                            & * (mu_dimles(i)**2 * qq_dimles(j) &
                                    & + alpha_dq * mu_dimles(i) * qq_dimles(i) * qq_dimles(j)) ! NOTE ts_2clj cancel each out
                    ! second are planar therms with only two different indexes
                    call calculate_abc_polar_sub(i, j, 0, 2, sqrt(param_m(i) * param_m(j)), &
                            & (param_m(i) * param_m(j) * param_m(i))**(third), &
                            & (param_m(i) * param_m(j) * param_m(j))**(third))

                    do k = j + 1, n_comp
                        sigm_jk = (param_sigma(j) + param_sigma(k)) / 2.0_dp
                        sigm_ki = (param_sigma(k) + param_sigma(i)) / 2.0_dp
                        ! === DD interaction ===
                        eees3s3s3_sssnnnmmm(i, j, k) = param_eps_kb(i) * param_eps_kb(j) * param_eps_kb(k) &
                                & * param_sigma(i)**3 * param_sigma(j)**3 * param_sigma(k)**3 &
                                & / (sigm_ij * sigm_jk * sigm_ki) &
                                & * param_n_u(i) * param_n_u(j) * param_n_u(k) &
                                & * mu_dimles(i) * mu_dimles(j) * mu_dimles(k)
                        eees3s3s3_sssnnnmmm(i, k, j) = eees3s3s3_sssnnnmmm(i, j, k)
                        eees3s3s3_sssnnnmmm(j, i, k) = eees3s3s3_sssnnnmmm(i, j, k)
                        eees3s3s3_sssnnnmmm(j, k, i) = eees3s3s3_sssnnnmmm(i, j, k)
                        eees3s3s3_sssnnnmmm(k, i, j) = eees3s3s3_sssnnnmmm(i, j, k)
                        eees3s3s3_sssnnnmmm(k, j, i) = eees3s3s3_sssnnnmmm(i, j, k)
                        ! === QQ interaction ===
                        eees5s5s5_s3s3s3nnnqqq(i, j, k) = param_eps_kb(i) * param_eps_kb(j) * param_eps_kb(k) &
                                & * param_sigma(i)**5 * param_sigma(j)**5 * param_sigma(k)**5 &
                                & / (sigm_ij * sigm_jk * sigm_ki)**3 &
                                & * param_n_q(i) * param_n_q(j) * param_n_q(k) &
                                & * qq_dimles(i) * qq_dimles(j) * qq_dimles(k)
                        eees5s5s5_s3s3s3nnnqqq(i, k, j) = eees5s5s5_s3s3s3nnnqqq(i, j, k)
                        eees5s5s5_s3s3s3nnnqqq(j, i, k) = eees5s5s5_s3s3s3nnnqqq(i, j, k)
                        eees5s5s5_s3s3s3nnnqqq(j, k, i) = eees5s5s5_s3s3s3nnnqqq(i, j, k)
                        eees5s5s5_s3s3s3nnnqqq(k, i, j) = eees5s5s5_s3s3s3nnnqqq(i, j, k)
                        eees5s5s5_s3s3s3nnnqqq(k, j, i) = eees5s5s5_s3s3s3nnnqqq(i, j, k)
                        ! === DQ interaction ===
                        eees4s4s4_s2s2s2nnqnqq(i, j, k) = param_eps_kb(i) * param_eps_kb(j) * param_eps_kb(k) &
                                & * param_sigma(i)**4 * param_sigma(j)**4 * param_sigma(k)**4 / (sigm_ij * sigm_jk * sigm_ki)**2 &
                                & * (mu_dimles(i) * mu_dimles(j) * qq_dimles(k) &
                                        & + alpha_dq * mu_dimles(i) * qq_dimles(j) * qq_dimles(k)) ! NOTE ts_2clj cancel each out
                        eees4s4s4_s2s2s2nnqnqq(i, k, j) = param_eps_kb(i) * param_eps_kb(j) * param_eps_kb(k) &
                                & * param_sigma(i)**4 * param_sigma(j)**4 * param_sigma(k)**4 / (sigm_ij * sigm_jk * sigm_ki)**2 &
                                & * (mu_dimles(i) * mu_dimles(k) * qq_dimles(j) &
                                        & + alpha_dq * mu_dimles(i) * qq_dimles(k) * qq_dimles(j)) ! NOTE ts_2clj cancel each out
                        eees4s4s4_s2s2s2nnqnqq(j, i, k) = param_eps_kb(i) * param_eps_kb(j) * param_eps_kb(k) &
                                & * param_sigma(i)**4 * param_sigma(j)**4 * param_sigma(k)**4 / (sigm_ij * sigm_jk * sigm_ki)**2 &
                                & * (mu_dimles(j) * mu_dimles(i) * qq_dimles(k) &
                                        & + alpha_dq * mu_dimles(j) * qq_dimles(i) * qq_dimles(k)) ! NOTE ts_2clj cancel each out
                        eees4s4s4_s2s2s2nnqnqq(j, k, i) = param_eps_kb(i) * param_eps_kb(j) * param_eps_kb(k) &
                                & * param_sigma(i)**4 * param_sigma(j)**4 * param_sigma(k)**4 / (sigm_ij * sigm_jk * sigm_ki)**2 &
                                & * (mu_dimles(j) * mu_dimles(k) * qq_dimles(i) &
                                        & + alpha_dq * mu_dimles(j) * qq_dimles(k) * qq_dimles(i)) ! NOTE ts_2clj cancel each out
                        eees4s4s4_s2s2s2nnqnqq(k, i, j) = param_eps_kb(i) * param_eps_kb(j) * param_eps_kb(k) &
                                & * param_sigma(i)**4 * param_sigma(j)**4 * param_sigma(k)**4 / (sigm_ij * sigm_jk * sigm_ki)**2 &
                                & * (mu_dimles(k) * mu_dimles(i) * qq_dimles(j) &
                                        & + alpha_dq * mu_dimles(k) * qq_dimles(i) * qq_dimles(j)) ! NOTE ts_2clj cancel each out
                        eees4s4s4_s2s2s2nnqnqq(k, j, i) = param_eps_kb(i) * param_eps_kb(j) * param_eps_kb(k) &
                                & * param_sigma(i)**4 * param_sigma(j)**4 * param_sigma(k)**4 / (sigm_ij * sigm_jk * sigm_ki)**2 &
                                & * (mu_dimles(k) * mu_dimles(j) * qq_dimles(i) &
                                        & + alpha_dq * mu_dimles(k) * qq_dimles(j) * qq_dimles(i)) ! NOTE ts_2clj cancel each out

                        !last are calculated the remaining all three indexes different
                        call calculate_abc_polar_sub(i, j, k, 3, 0.0_dp, &
                                & (param_m(i) * param_m(j) * param_m(k))**(third), 0.0_dp)
                    end do
                end do
            end do
            ! 				eees4s4s4_s2s2s2nnqnqq = eees4s4s4_s2s2s2nnqnqq*1.0E19**3
            ! 				ees3s5_s5nq = ees3s5_s5nq*1.0E19**2
            write(stddeb,*) "polrules a,b,c ees3s3_s3nnmm, eees3s3s3_sssnnnmmm, ees5s5_s7nnqq, eees5s5s5_s3s3s3nnnqqq, ees3s5_s5nq, eees4s4s4_s2s2s2nnqnqq are initialized " !DEBUG
            init_flag(6) = .true.
            ! ----- deallocation section ----
            deallocate(mu_dimles, qq_dimles)
            deallocate(ts_2clj)
        end if

        if ( (init_flag(5) .eqv. .true.) .and. (init_flag(6) .eqv. .true.) ) then
            success = .true.
        else
            success = .false.
        end if

        write (stddeb,*) "<= initialize_polrules_fun" !DEBUG
    end function initialize_polrules_fun

    subroutine calculate_abc_polar_sub(i_in, j_in, k_in, type_in, val1_in, val2_in, val3_in)
        ! the function that computes the individual members of the a__dd,b__dd,c__dd; a__qq,b__qq,c__qq; a__dq,b__dq,c__dq matrixes
        ! NOTE TODO this procedure can be optimalized for simplier computation - use the array manipulation
        implicit none
        ! variable section
        integer, intent(in) :: i_in, j_in, k_in ! the input indexes
        integer, intent(in) :: type_in ! indicate type of execution for 1=iii&ii, 2=ijj&ij, 3=ijk
        real(dp), intent(in) :: val1_in, val2_in ! the value of the parameter 1=ij
        real(dp), intent(in) :: val3_in ! the value of the parameter  2,3=ijk
        real(dp) :: val1, val2, val3 ! the local counterparts to input variables -enable local change
        integer, parameter :: nn = 5 ! the sum over n from 1 to nn
        integer :: n ! counters

        write (stddeb,*) "=> calculate_abc_polar_sub" !DEBUG
        ! control that m_ij <2 m_ijk<2
        if (val1_in>2.0_dp) then
            val1 = 2.0_dp
        else
            val1 = val1_in
        end if
        if (val2_in>2.0_dp) then
            val2 = 2.0_dp
        else
            val2 = val2_in
        end if
        if (val3_in>2.0_dp) then
            val3 = 2.0_dp
        else
            val3 = val3_in
        end if

        ! the mij or mijk decision
        select case (type_in)
        case(1) ! diagonal terms i=j, i=j=k
            do n = 1, nn
                ! === DIPOLAR ======================================================
                a__dd(n, i_in, i_in) = a0_dd(n, 1) + (1.0_dp - 1.0_dp / val1) * a0_dd(n, 2) + & ! aii
                        & (1.0_dp - 1.0_dp / val1) * (1.0_dp - 2.0_dp / val1) * a0_dd(n, 3)
                b__dd(n, i_in, i_in) = b0_dd(n, 1) + (1.0_dp - 1.0_dp / val1) * b0_dd(n, 2) + & ! bii
                        & (1.0_dp - 1.0_dp / val1) * (1.0_dp - 2.0_dp / val1) * b0_dd(n, 3)
                !-------------------------------------------------------------------
                c__dd(n, i_in, i_in, i_in) = c0_dd(n, 1) + (1.0_dp - 1.0_dp / val2) * c0_dd(n, 2) + & ! ciii
                        & (1.0_dp - 1.0_dp / val2) * (1.0_dp - 2.0_dp / val2) * c0_dd(n, 3)
                ! === QUADRUPOLAR ==================================================
                a__qq(n, i_in, i_in) = a0_qq(n, 1) + (1.0_dp - 1.0_dp / val1_in) * a0_qq(n, 2) + & ! aii
                        & (1.0_dp - 1.0_dp / val1_in) * (1.0_dp - 2.0_dp / val1_in) * a0_qq(n, 3)
                b__qq(n, i_in, i_in) = b0_qq(n, 1) + (1.0_dp - 1.0_dp / val1_in) * b0_qq(n, 2) + & ! bii
                        & (1.0_dp - 1.0_dp / val1_in) * (1.0_dp - 2.0_dp / val1_in) * b0_qq(n, 3)
                !-------------------------------------------------------------------
                c__qq(n, i_in, i_in, i_in) = c0_qq(n, 1) + (1.0_dp - 1.0_dp / val2_in) * c0_qq(n, 2) + & ! ciii
                        & (1.0_dp - 1.0_dp / val2_in) * (1.0_dp - 2.0_dp / val2_in) * c0_qq(n, 3)
                ! === DIPOLAR-QUADRUPOLAR ==========================================
                a__dq(n, i_in, i_in) = a0_dq(n, 1) + (1.0_dp - 1.0_dp / val1_in) * a0_dq(n, 2) + & ! aii
                        & (1.0_dp - 1.0_dp / val1_in) * (1.0_dp - 2.0_dp / val1_in) * a0_dq(n, 3)
                b__dq(n, i_in, i_in) = b0_dq(n, 1) + (1.0_dp - 1.0_dp / val1_in) * b0_dq(n, 2) + & ! bii
                        & (1.0_dp - 1.0_dp / val1_in) * (1.0_dp - 2.0_dp / val1_in) * b0_dq(n, 3)
                !-------------------------------------------------------------------
                c__dq(n, i_in, i_in, i_in) = c0_dq(n, 1) + (1.0_dp - 1.0_dp / val2_in) * c0_dq(n, 2) ! ciii
            end do
        case(2) ! planar terms i~j, i~j i=k
            do n = 1, nn
                ! === DIPOLAR ======================================================
                a__dd(n, i_in, j_in) = a0_dd(n, 1) + (1.0_dp - 1.0_dp / val1) * a0_dd(n, 2) + & ! aij
                        & (1.0_dp - 1.0_dp / val1) * (1.0_dp - 2.0_dp / val1) * a0_dd(n, 3)
                a__dd(n, j_in, i_in) = a__dd(n, i_in, j_in) ! symmetry aji
                b__dd(n, i_in, j_in) = b0_dd(n, 1) + (1.0_dp - 1.0_dp / val1) * b0_dd(n, 2) + & ! bij
                        & (1.0_dp - 1.0_dp / val1) * (1.0_dp - 2.0_dp / val1) * b0_dd(n, 3)
                b__dd(n, j_in, i_in) = b__dd(n, i_in, j_in) ! symmetry bji
                !-------------------------------------------------------------------
                c__dd(n, i_in, j_in, i_in) = c0_dd(n, 1) + (1.0_dp - 1.0_dp / val2) * c0_dd(n, 2) + & ! ciji
                        & (1.0_dp - 1.0_dp / val2) * (1.0_dp - 2.0_dp / val2) * c0_dd(n, 3)
                c__dd(n, j_in, i_in, i_in) = c__dd(n, i_in, j_in, i_in) ! symmetry cjii
                c__dd(n, i_in, i_in, j_in) = c__dd(n, i_in, j_in, i_in) ! symmetry ciij
                !-------------------------------------------------------------------
                c__dd(n, i_in, j_in, j_in) = c0_dd(n, 1) + (1.0_dp - 1.0_dp / val3) * c0_dd(n, 2) + & ! cijj
                        & (1.0_dp - 1.0_dp / val3) * (1.0_dp - 2.0_dp / val3) * c0_dd(n, 3)
                c__dd(n, j_in, i_in, j_in) = c__dd(n, i_in, j_in, j_in) ! symmetry cjij
                c__dd(n, j_in, j_in, i_in) = c__dd(n, i_in, j_in, j_in) ! symmetry cjji

                ! === QUADRUPOLAR ==================================================
                a__qq(n, i_in, j_in) = a0_qq(n, 1) + (1.0_dp - 1.0_dp / val1_in) * a0_qq(n, 2) + & ! aij
                        & (1.0_dp - 1.0_dp / val1_in) * (1.0_dp - 2.0_dp / val1_in) * a0_qq(n, 3)
                a__qq(n, j_in, i_in) = a__qq(n, i_in, j_in) ! symmetry aji
                b__qq(n, i_in, j_in) = b0_qq(n, 1) + (1.0_dp - 1.0_dp / val1_in) * b0_qq(n, 2) + & ! bij
                        & (1.0_dp - 1.0_dp / val1_in) * (1.0_dp - 2.0_dp / val1_in) * b0_qq(n, 3)
                b__qq(n, j_in, i_in) = b__qq(n, i_in, j_in) ! symmetry bji
                !-------------------------------------------------------------------
                c__qq(n, i_in, j_in, i_in) = c0_qq(n, 1) + (1.0_dp - 1.0_dp / val2_in) * c0_qq(n, 2) + & ! ciji
                        & (1.0_dp - 1.0_dp / val2_in) * (1.0_dp - 2.0_dp / val2_in) * c0_qq(n, 3)
                c__qq(n, j_in, i_in, i_in) = c__qq(n, i_in, j_in, i_in) ! symmetry cjii
                c__qq(n, i_in, i_in, j_in) = c__qq(n, i_in, j_in, i_in) ! symmetry ciij
                !-------------------------------------------------------------------
                c__qq(n, i_in, j_in, j_in) = c0_qq(n, 1) + (1.0_dp - 1.0_dp / val3_in) * c0_qq(n, 2) + & ! cijj
                        & (1.0_dp - 1.0_dp / val3_in) * (1.0_dp - 2.0_dp / val3_in) * c0_qq(n, 3)
                c__qq(n, j_in, i_in, j_in) = c__qq(n, i_in, j_in, j_in) ! symmetry cjij
                c__qq(n, j_in, j_in, i_in) = c__qq(n, i_in, j_in, j_in) ! symmetry cjji

                ! === DIPOLAR-QUADRUPOLAR ==========================================
                a__dq(n, i_in, j_in) = a0_dq(n, 1) + (1.0_dp - 1.0_dp / val1_in) * a0_dq(n, 2) + & ! aij
                        & (1.0_dp - 1.0_dp / val1_in) * (1.0_dp - 2.0_dp / val1_in) * a0_dq(n, 3)
                a__dq(n, j_in, i_in) = a__dq(n, i_in, j_in) ! symmetry aji
                b__dq(n, i_in, j_in) = b0_dq(n, 1) + (1.0_dp - 1.0_dp / val1_in) * b0_dq(n, 2) + & ! bij
                        & (1.0_dp - 1.0_dp / val1_in) * (1.0_dp - 2.0_dp / val1_in) * b0_dq(n, 3)
                b__dq(n, j_in, i_in) = b__dq(n, i_in, j_in) ! symmetry bji
                !-------------------------------------------------------------------
                c__dq(n, i_in, j_in, i_in) = c0_dq(n, 1) + (1.0_dp - 1.0_dp / val2_in) * c0_dq(n, 2) ! ciji
                c__dq(n, j_in, i_in, i_in) = c__dq(n, i_in, j_in, i_in) ! symmetry cjii
                c__dq(n, i_in, i_in, j_in) = c__dq(n, i_in, j_in, i_in) ! symmetry ciij
                !-------------------------------------------------------------------
                c__dq(n, i_in, j_in, j_in) = c0_dq(n, 1) + (1.0_dp - 1.0_dp / val3_in) * c0_dq(n, 2) ! cijj
                c__dq(n, j_in, i_in, j_in) = c__dq(n, i_in, j_in, j_in) ! symmetry cjij
                c__dq(n, j_in, j_in, i_in) = c__dq(n, i_in, j_in, j_in) ! symmetry cjji
            end do
        case(3) ! nondiagonal terms i~j, i~j~k
            do n = 1, nn
                ! === DIPOLAR ======================================================
                c__dd(n, i_in, j_in, k_in) = c0_dd(n, 1) + (1.0_dp - 1.0_dp / val2) * c0_dd(n, 2) + & ! cijk
                        & (1.0_dp - 1.0_dp / val2) * (1.0_dp - 2.0_dp / val2) * c0_dd(n, 3)
                c__dd(n, i_in, k_in, j_in) = c__dd(n, i_in, j_in, k_in) ! symmetry cikj
                c__dd(n, j_in, i_in, k_in) = c__dd(n, i_in, j_in, k_in) ! symmetry cjik
                c__dd(n, j_in, k_in, i_in) = c__dd(n, i_in, j_in, k_in) ! symmetry cjki
                c__dd(n, k_in, i_in, j_in) = c__dd(n, i_in, j_in, k_in) ! symmetry ckij
                c__dd(n, k_in, j_in, i_in) = c__dd(n, i_in, j_in, k_in) ! symmetry ckji
                ! === QUADRUPOLAR ==================================================
                c__qq(n, i_in, j_in, k_in) = c0_qq(n, 1) + (1.0_dp - 1.0_dp / val2_in) * c0_qq(n, 2) + & ! cijk
                        & (1.0_dp - 1.0_dp / val2_in) * (1.0_dp - 2.0_dp / val2_in) * c0_qq(n, 3)
                c__qq(n, i_in, k_in, j_in) = c__qq(n, i_in, j_in, k_in) ! symmetry cikj
                c__qq(n, j_in, i_in, k_in) = c__qq(n, i_in, j_in, k_in) ! symmetry cjik
                c__qq(n, j_in, k_in, i_in) = c__qq(n, i_in, j_in, k_in) ! symmetry cjki
                c__qq(n, k_in, i_in, j_in) = c__qq(n, i_in, j_in, k_in) ! symmetry ckij
                c__qq(n, k_in, j_in, i_in) = c__qq(n, i_in, j_in, k_in) ! symmetry ckji
                ! === DIPOLAR-QUADRUPOLAR ==========================================
                c__dq(n, i_in, j_in, k_in) = c0_dq(n, 1) + (1.0_dp - 1.0_dp / val2_in) * c0_dq(n, 2)  ! cijk
                c__dq(n, i_in, k_in, j_in) = c__dq(n, i_in, j_in, k_in) ! symmetry cikj
                c__dq(n, j_in, i_in, k_in) = c__dq(n, i_in, j_in, k_in) ! symmetry cjik
                c__dq(n, j_in, k_in, i_in) = c__dq(n, i_in, j_in, k_in) ! symmetry cjki
                c__dq(n, k_in, i_in, j_in) = c__dq(n, i_in, j_in, k_in) ! symmetry ckij
                c__dq(n, k_in, j_in, i_in) = c__dq(n, i_in, j_in, k_in) ! symmetry ckji
            end do
        end select

        write (stddeb,*) "<= calculate_abc_polar_sub" !DEBUG
    end subroutine calculate_abc_polar_sub

    ! =============================  ASSOCIATION CONTRIBUTION  ===============================

    subroutine compute_qudratic_sub(a_in, b_in, c_in, xx_dd)
        ! the solver of quadratic equation for + case only the respective first derivative
        ! the input parameters are array 4 of coeffs and its (1,2,3) derivatives if necesarry
        ! the outputs are simple scalars xx_d?
        ! NOTE formulas were verified using wolfram mathematica
        ! TODO could be simplified for dominant case where b_in(2) = 0, c_in(0) = 0
        implicit none
        ! variable section
        real(dp), intent(in), dimension(4) :: a_in
        real(dp), intent(in), dimension(4) :: b_in
        real(dp), intent(in), dimension(4) :: c_in
        real(dp), intent(out) :: xx_dd
        ! local variable
        real(dp) :: dd, dd_dd

        write (stddeb,*) "=> compute_qudratic_sub" !DEBUG
        if (b_in(1)**2 - 4.0_dp * a_in(1) * c_in(1) <0) then
            write (stderr, *) 'Error in quadratic3 sub, sqrt negative number - check system setup.'
            stop
        end if

        dd = sqrt(b_in(1)**2 - 4.0_dp * a_in(1) * c_in(1))
        dd_dd = (-2.0_dp * c_in(1) * a_in(2) &
                & - 2.0_dp * a_in(1) * c_in(2) &
                & + b_in(1) * b_in(2)) / dd
        ! NOTE not perfectly optimal - some operation can be adjected to save compuation load
        xx_dd = (+b_in(1) * a_in(2) &
                & - dd * a_in(2) &
                & - a_in(1) * b_in(2) &
                & + a_in(1) * dd_dd) / (2.0_dp * a_in(1)**2)
        write (stddeb,*) "<= compute_qudratic_sub" !DEBUG
    end subroutine compute_qudratic_sub

    subroutine compute_qudratic3_sub(a_in, b_in, c_in, xx, xx_drho, xx_drho2, xx_drho3)
        ! the solver of quadratic equation for + case and respective first, second and third derivative
        ! the input parameters are array 4 of coeffs and its (1,2,3) derivatives
        ! the outputs are simple scalars xx, xx_drho, xx_drho2, xx_drho3
        ! NOTE formulas were verified using wolfram mathematica, partially tested against VV code
        ! TODO could be simplified for dominant case where b_in(2) = 0, c_in(0) = 0
        implicit none
        ! variable section
        real(dp), intent(in), dimension(4) :: a_in
        real(dp), intent(in), dimension(4) :: b_in
        real(dp), intent(in), dimension(4) :: c_in
        real(dp), intent(out) :: xx
        real(dp), intent(out) :: xx_drho
        real(dp), intent(out) :: xx_drho2
        real(dp), intent(out) :: xx_drho3
        ! local variable
        real(dp) :: dd, dd_drho, dd_drho2, dd_drho3

        write (stddeb,*) "=> compute_qudratic3_sub" !DEBUG
        if (b_in(1)**2 - 4.0_dp * a_in(1) * c_in(1) <0) then
            write (stderr, *) 'Error in quadratic3 sub, sqrt negative number - check system setup.'
            stop
        end if

        dd = sqrt(b_in(1)**2 - 4.0_dp * a_in(1) * c_in(1))
        dd_drho = (-2.0_dp * c_in(1) * a_in(2) &
                & - 2.0_dp * a_in(1) * c_in(2) &
                & + b_in(1) * b_in(2)) / dd
        dd_drho2 = +(-4.0_dp * c_in(1) * a_in(3) &
                & - 8.0_dp * a_in(2) * c_in(2) &
                & - 4.0_dp * a_in(1) * c_in(3) &
                & + 2.0_dp * b_in(1) * b_in(3) &
                & + 2.0_dp * b_in(2)**2) / (2.0_dp * dd) &
                & - (-4.0_dp * c_in(1) * a_in(2) &
                        & - 4.0_dp * a_in(1) * c_in(2) &
                        & + 2.0_dp * b_in(1) * b_in(2))**2 / (4.0_dp * dd**3)
        if (sec_der .eqv. .true.) then
            dd_drho3 = +(3.0_dp * (-4.0_dp * c_in(1) * a_in(2) &
                    & - 4.0_dp * a_in(1) * c_in(2) &
                    & + 2.0_dp * b_in(1) * b_in(2))**3) / (8.0_dp * dd**5) &
                    & - (3.0_dp * (-4.0_dp * c_in(1) * a_in(2) &
                            & - 4.0_dp * a_in(1) * c_in(2) &
                            & + 2.0_dp * b_in(1) * b_in(2)) &
                            & * (-4.0_dp * c_in(1) * a_in(3) &
                                    & - 8.0_dp * a_in(2) * c_in(2) &
                                    & - 4.0_dp * a_in(1) * c_in(3) &
                                    & + 2.0_dp * b_in(1) * b_in(3) &
                                    & + 2.0_dp * b_in(2)**2)) / (4.0_dp * dd**3) &
                    & + (-4.0_dp * a_in(4) * c_in(1) &
                            & - 12.0_dp * a_in(3) * c_in(2) &
                            & - 12.0_dp * a_in(2) * c_in(3) &
                            & - 4.0_dp * a_in(1) * c_in(4) &
                            & + 2.0_dp * b_in(1) * b_in(4) &
                            & + 6.0_dp * b_in(2) * b_in(3)) / (2.0_dp * dd)
        else
            dd_drho3 = 0.0_dp
        end if
        ! NOTE not perfectly optimal - some operation can be adjected to save compuation load
        xx = (-b_in(1) + dd) / (2.0_dp * a_in(1))
        xx_drho = (+b_in(1) * a_in(2) &
                & - dd * a_in(2) &
                & - a_in(1) * b_in(2) &
                & + a_in(1) * dd_drho) / (2.0_dp * a_in(1)**2)
        xx_drho2 = -(a_in(2) * (dd_drho - b_in(2))) / a_in(1)**2 &
                & + (dd - b_in(1)) * (+(2.0_dp * a_in(2)**2) / a_in(1)**3 &
                        & - (a_in(3)) / a_in(1)**2) / 2.0_dp &
                & + (dd_drho2 - b_in(3)) / (2.0_dp * a_in(1))
        if (sec_der .eqv. .true.) then
            xx_drho3 = (+3.0_dp * a_in(1)**2 * a_in(2) * (b_in(3) - dd_drho2) &
                    & + 3.0_dp * a_in(1) * (a_in(1) * a_in(3) - 2.0_dp * a_in(2)**2) * (b_in(2) - dd_drho) &
                    & + (b_in(1) - dd) * (a_in(1)**2 * a_in(4) + 6.0_dp * a_in(2)**3 - 6.0_dp * a_in(1) * a_in(2) * a_in(3)) &
                    & + a_in(1)**3 * (-b_in(4) + dd_drho3)) / (2.0_dp * a_in(1)**4)
        else
            xx_drho3 = 0.0_dp
        end if
        write (stddeb,*) "<= compute_qudratic3_sub" !DEBUG
    end subroutine compute_qudratic3_sub

end module contrib_mod
