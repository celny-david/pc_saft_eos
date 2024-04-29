subroutine z_disp_sub (rho_in)
    ! subroutine for the compresibility factor contribution from hard chain DISP term calculation
    ! - preapre the necesarry components such as initialize_a0b0_disp_fun, initialize_dzeta_fun, initialize_combrules_fun and compute_ab_disp_fun
    !	* calling only when necesarry
    ! - the actual initialization
    ! - computation of the dispersion compresibility contribution
    !	* helper terms such as I1, I2 are computed first
    !	* abbreviations C1, C2 are computed afterwards
    !	* also second derivative if required
    ! - parameter update in control module
    !
    !	== LAST MODIFICATIONS ==
    !	15.5.2017 - D.C - descriptions and comments
    !	20.6.2017 - D.C - init_flag_g dependency omitment -> moved to parent control module -> replaced by single flag

    use control_mod, only : dp, sec_der, pi, &
            &    stdout, stderr, stddeb, stdlog, &
            &    write_array_1d, &
            ! &   z_hc, z_hc_drho, z_hc_drho2, z_hc_drho3, z_hc_drho4, & ! OLD presently unused
            &   z_disp, z_disp_drho, z_disp_drho2, &
            &   init_flag, init_flag_zdisp
    use contrib_mod, only : a__disp, b__disp, dzeta, m_mean, m2es3, m2e2s3, &
            & initialize_a0b0_disp_fun, initialize_dzeta_fun, initialize_combrules_fun, compute_ab_disp_fun
    implicit none
    ! variable section IN/OUT
    real(dp), intent(in) :: rho_in ! the input density
    ! local variables outputed into control_mod
    real(dp) :: z_disp_       ! the dispersive compressibility contribution
    real(dp) :: z_disp_drho_  ! the first density derivative of dispersive compressibility contribution
    real(dp) :: z_disp_drho2_ ! the second density derivative of dispersive compressibility contribution
    ! local variables
    real(dp), dimension(0:3) :: dzeta_r ! the dzeta vector dzeta_n for n=0,1,2,3 * rho !!
    real(dp) :: ii1, ii2 ! the approximated integrals computed as power series
    real(dp) :: ii1_drho, ii2_drho ! the approximated integrals first density derivatives
    real(dp) :: ii1_drho2, ii2_drho2 ! the approximated integrals second density derivatives
    real(dp) :: ii1_drho3, ii2_drho3 ! the approximated integrals third density derivatives
    real(dp) :: cc1, cc2 ! the C1 C2 abbreiations used in Z_disp computation
    real(dp) :: cc1_drho, cc2_drho ! the C1 C2 first density derivative
    ! 	real(dp) :: cc1_drho2, cc2_drho2 ! the C1 C2 second density derivative ! OLD presently unused
    !	real(dp) :: tmp_div ! temporary abbreviated divisor (1 + z_hc + rho*z_hc_drho) ! OLD presently unused
    integer :: i ! counters

    !------- PREPARING ------
    write (stddeb,*) "==> z_disp_sub" !DEBUG
    ! call for the required helper procedures
    if (init_flag(1) .eqv. .false.) then
        if (initialize_a0b0_disp_fun() .eqv. .false.) then
            write (stderr, *) 'Error in helper function initialize_a0b0_disp_fun'
            stop
        end if
    end if
    if (init_flag(2) .eqv. .false.) then
        if (initialize_dzeta_fun() .eqv. .false.) then
            write (stderr, *) 'Error in helper function initialize_dzeta_fun'
            stop
        end if
    end if
    if (init_flag(3) .eqv. .false.) then
        if (initialize_combrules_fun() .eqv. .false.) then
            write (stderr, *) 'Error in helper function initialize_combrules_fun'
            stop
        end if
    end if
    if (init_flag(4) .eqv. .false.) then
        if (compute_ab_disp_fun() .eqv. .false.) then
            write (stderr, *) 'Error in helper function compute_ab_disp_fun'
            stop
        end if
    end if

    ! 	! OLD presently is not required due to the different implementation
    ! 	! call for the required level 1 subroutines z_hc
    ! 	if (z_hc == 0.0_dp) then
    ! 		call z_hc_sub(rho_in)
    ! 	end if

    ! the actual dzeta (multiplied by density)
    dzeta_r = dzeta * rho_in

    ! DEBUG control computation
    write (stddeb,*) '-----  QQ  -----'  ! DEBUG just for brevity of debug
    call write_array_1d (dzeta,"dzeta      : ",stddeb)
    call write_array_1d (dzeta_r,"dzeta_r    : ",stddeb)
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug
    ! DEBUG end of control computation

    ! initialization
    ii1 = 0.0_dp
    ii1_drho = 0.0_dp
    ii1_drho2 = 0.0_dp
    ii1_drho3 = 0.0_dp

    ii2 = 0.0_dp
    ii2_drho = 0.0_dp
    ii2_drho2 = 0.0_dp
    ii2_drho3 = 0.0_dp

    ! ----- PREPARING INTERMEDIATE VARIABLES -----
    ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.16, A.17)
    do i = 1, 7
        ! derived by maple with the property d/deta = rho/dzeta_r(3) *d/drho
        ii1 = ii1 + a__disp(i) * dzeta_r(3)**(i - 1.0_dp) ! A.16
        ii2 = ii2 + b__disp(i) * dzeta_r(3)**(i - 1.0_dp) ! A.17

        ii1_drho = ii1_drho + (i) * a__disp(i) * dzeta_r(3)**(i - 1.0_dp) ! corrected derivative d(eta*I1)/deta = 1/c_3 * d(rho*I1)/drho
        ii2_drho = ii2_drho + (i) * b__disp(i) * dzeta_r(3)**(i - 1.0_dp) ! corrected derivative d(eta*I2)/deta = 1/c_3 * d(rho*I1)/drho
        !  		write (*,"('adisp: ',f30.15,' | bdisp: ',f30.15,/)"), a__disp(i), b__disp(i) ! DEBUG testing of parameters
        ii1_drho2 = ii1_drho2 + (i) * (i - 1.0_dp) * a__disp(i) * dzeta_r(3)**(i - 2.0_dp) ! corrected derivative d2(eta*I1)/d2eta = 1/c_3**2 * d2(rho*I1)/d2rho
        ii2_drho2 = ii2_drho2 + (i) * (i - 1.0_dp) * b__disp(i) * dzeta_r(3)**(i - 2.0_dp) ! corrected derivative d2(eta*I2)/d2eta = 1/c_3**2 * d2(rho*I1)/d2rho

        if (sec_der .eqv. .true.) then
            ii1_drho3 = ii1_drho3 + (i) * (i - 1.0_dp) * (i - 2.0_dp) * a__disp(i) * dzeta_r(3)**(i - 3.0_dp) ! corrected derivative d2(eta*I1)/d2eta = 1/c_3**2 * d2(rho*I1)/d2rho
            ii2_drho3 = ii2_drho3 + (i) * (i - 1.0_dp) * (i - 2.0_dp) * b__disp(i) * dzeta_r(3)**(i - 3.0_dp) ! corrected derivative d2(eta*I2)/d2eta = 1/c_3**2 * d2(rho*I1)/d2rho
        end if
    end do

    write (stddeb,*) "ii1      : ",ii1 ! DC DEBUG
    write (stddeb,*) "ii2      : ",ii2 ! DC DEBUG
    write (stddeb,*) "ii1_drho : ",ii1_drho ! DC DEBUG
    write (stddeb,*) "ii2_drho : ",ii2_drho ! DC DEBUG
    write (stddeb,*) "ii1_drho2: ",ii1_drho ! DC DEBUG
    write (stddeb,*) "ii2_drho2: ",ii2_drho ! DC DEBUG
    if (sec_der .eqv. .true.) then
        write (stddeb,*) "ii1_drho3: ",ii1_drho3 ! DC DEBUG
        write (stddeb,*) "ii2_drho3: ",ii2_drho3 ! DC DEBUG
    end if
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    ! set up the temporary variable
    ! 	tmp_div = (1.0_dp + z_hc + rho_in*z_hc_drho) ! OLD the temporary invalid implementation in rho formalism ! TODO correct dzeta into rho
    ! 	write (*,"('tmp_div: ',f25.10)"), tmp_div

    ! ----- C1 ----- ! OLD until corrected
    ! BUG inconsistency with VV code in computation of the cc1, cc2
    ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.11) LHS
    ! BUG MAJOR numeric inconsistency - temporary implementation with in dzeta formalism
    ! 	cc1 = 1.0_dp/tmp_div ! OLD the temporary invalid implementation in rho formalism ! TODO correct dzeta into rho
    ! 	! derived differentiate and simplified by maple
    ! 	cc1_drho = (2*z_hc_drho + rho_in*z_hc_drho2)/tmp_div**2
    ! 	if (sec_der==.true.) then
    ! 		cc1_drho2= 2*(2*z_hc_drho + rho_in*z_hc_drho2)**2/tmp_div**3 &
    ! 				&  -(3*z_hc_drho2 + rho_in*z_hc_drho3)/tmp_div**2
    ! 	end if

    ! ----- C2 ----- ! OLD until corrected
    ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.31)
    ! derived and simplified by maple (computed via dC1/deta = rho/eta *dC1/drho)
    ! BUG MAJOR numeric inconsistency - temporary implementation with in dzeta formalism
    ! 	cc2 = cc1_drho/dzeta(3) ! OLD the temporary invalid implementation in rho formalism ! TODO correct dzeta into rho
    ! 	cc2_drho = cc1_drho2/dzeta(3)
    ! 	if (sec_der==.true.) then
    ! 		! TODO for implementation of the cc2_drho2 -> compute z_hc_drho3, z_hc_drho4
    ! 		! derived differentiate and simplified by maple
    ! 		cc2_drho2 = -6*(2*z_hc_drho + rho_in*z_hc_drho2)**3/tmp_div**4 &
    ! 				&   +6*(3*z_hc_drho2 + rho_in*z_hc_drho3)*(2*z_hc_drho + rho_in*z_hc_drho2)/tmp_div**3 &
    ! 				&   -(4*z_hc_drho3 + rho_in*z_hc_drho4)/tmp_div**2
    ! 		cc2_drho2 = cc2_drho2/dzeta(3)
    ! 	end if

    ! ----- C1 C2 VV implementation -----
    ! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.11) RHS
    ! used the VV implementation
    !	c1_con= 1/(1 + m_mean*(8*dense-2*dense^2)/zms^4 +
    !			+ (1-m_mean)*(20*dense-27*dense^2 + 12*dense^3-2*dense^4)/(zms*(2-dense))^2); ! VV code reference
    cc1 = 1.0_dp / (1.0_dp + m_mean * (8.0_dp * dzeta_r(3)&
            & - 2.0_dp * dzeta_r(3)**2) / (1.0_dp - dzeta_r(3))**4 &
            & + (1 - m_mean) * (20.0_dp * dzeta_r(3) &
                    & - 27.0_dp * dzeta_r(3)**2 &
                    & + 12.0_dp * dzeta_r(3)**3 &
                    & - 2.0_dp * dzeta_r(3)**4) / ((1.0_dp - dzeta_r(3)) * (2.0_dp - dzeta_r(3)))**2)

    !	c2_con= - c1_con*c1_con*(m_mean*(-4*dense^2+20*dense+8)/zms^5 +
    !			+(1-m_mean)*(2*dense^3+12*dense^2-48*dense+40)/(zms*(2-dense))^3); ! VV code reference
    cc2 = -(cc1**2) * (m_mean * (-4.0_dp * dzeta_r(3)**2 &
            & + 20.0_dp * dzeta_r(3) &
            & + 8.0_dp) / (1.0_dp - dzeta_r(3))**5 &
            & + (1.0_dp - m_mean) * (2.0_dp * dzeta_r(3)**3 &
                    & + 12.0_dp * dzeta_r(3)**2 &
                    & - 48.0_dp * dzeta_r(3) &
                    & + 40.0_dp) / ((1.0_dp - dzeta_r(3)) * (2.0_dp - dzeta_r(3)))**3)

    !	c3_con= 2*c2_con*c2_con/c1_con -
    !			- c1_con*c1_con*(m_mean*(-12*dense^2+72*dense+60)/zms^6 +
    !			+ (1-m_mean)*(-6*dense^4-48*dense^3+288*dense^2-480*dense+264)
    !				/(zms*(2-dense))^4); ! VV code reference
    cc1_drho = 2.0_dp * cc2**2 / cc1 &
            & - (cc1**2) * (m_mean * (-12.0_dp * dzeta_r(3)**2 &
                    & + 72.0_dp * dzeta_r(3) &
                    & + 60.0_dp) / (1.0_dp - dzeta_r(3))**6 &
                    & + (1.0_dp - m_mean) * (-6.0_dp * dzeta_r(3)**4 &
                            & - 48.0_dp * dzeta_r(3)**3 &
                            & + 288.0_dp * dzeta_r(3)**2 &
                            & - 480.0_dp * dzeta_r(3) &
                            & + 264.0_dp) / ((1.0_dp - dzeta_r(3)) * (2.0_dp - dzeta_r(3)))**4)

    ! 	c4_con = 4/c1_con*c2_con*c3_con
    !			- 2/c1_con^2*c2_con^3
    !			+ 2/c1_con*c2_con*(c3_con-2/c1_con*c2_con*c2_con)
    !			- c1_con*c1_con*( m_mean*(432+336*dense-48*dense^2)/zms^7
    !				+ (1-m_mean)*(2208-5280*dense+4800*dense^2-1920*dense^3+240*dense^4+24*dense^5)
    !				/(2-3*dense+dense^2)^5 ); ! VV code reference
    cc2_drho = 4.0_dp * (cc2 * cc1_drho) / cc1 &
            & - 2 * cc2**3 / cc1**2   &
            & + 2 * cc2 * (cc1_drho - &
                    &               2.0_dp * cc2**2 / cc1) / cc1 &
            & - cc1**2 * (m_mean * (-48.0_dp * dzeta_r(3)**2 &
                    & + 336.0_dp * dzeta_r(3) &
                    & + 432.0_dp) / (1.0_dp - dzeta_r(3))**7 &
                    & + (1 - m_mean) * (24.0_dp * dzeta_r(3)**5 &
                            & + 240.0_dp * dzeta_r(3)**4 &
                            & - 1920.0_dp * dzeta_r(3)**3 &
                            & + 4800.0_dp * dzeta_r(3)**2 &
                            & - 5280.0_dp * dzeta_r(3) &
                            & + 2208.0_dp) / (2.0_dp - 3.0_dp * dzeta_r(3) + dzeta_r(3)**2)**5)


    write (stddeb,*) "cc1      : ",cc1 ! DC DEBUG
    write (stddeb,*) "cc2      : ",cc2 ! DC DEBUG
    write (stddeb,*) "cc1_drho : ",cc1_drho ! DC DEBUG
    write (stddeb,*) "cc2_drho : ",cc2_drho ! DC DEBUG
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug

    ! 	! ----- COMPUTATION OF Z_disp ----- ! OLD
    ! 	! based upon paper Gross-Sadowski 2001 application of PC-SAFT eq (A.28)
    ! 	! used the property dii1/deta = rho/dzeta_r(3) *dii1/drho = 1/dzeta(3) *dii1/drho
    ! 	z_disp_ =  -2*pi*rho_in*ii1_drho*m2es3/dzeta(3) &
    ! 			& -pi*rho_in*m_mean*(cc1*ii2_drho/dzeta(3) + cc2*dzeta_r(3)*ii2)*m2e2s3
    !
    ! 	! derived differentiate and gradually simplified by maple
    ! 	! used the property d?/deta = rho/dzeta_r(3) *d?/drho = 1/dzeta(3) *d?/drho
    ! 	z_disp_drho_ = -2*pi*m2es3*(ii1_drho+rho_in*ii1_drho2)/dzeta(3)+ & ! terms with m2es3
    ! 				& -pi*m_mean*m2e2s3/dzeta(3)*( &               ! terms with m2e2s3
    ! 					&	cc1*ii2_drho & !constant terms in rho
    ! 					&	+rho_in *(2*cc2*ii2*dzeta(3)**3 + cc1_drho*ii2_drho +cc1*ii2_drho2) & !linear therms in rho_in
    ! 					&	+(rho_in*dzeta(3))**2 *(cc2_drho*ii2 + ii2_drho*cc2) ) ! quadratic terms in rho
    !
    ! 	if (sec_der==.true.) then
    ! 	! derived differentiate and gradually simplified by maple
    ! 	! used the property d?/deta = rho/dzeta_r(3) *d?/drho = 1/dzeta(3) *d?/drho
    ! 	z_disp_drho2_= -2*pi*m2es3*(2*ii1_drho2+rho_in*ii1_drho3)/dzeta(3)+ & ! terms with m2es3
    ! 				& -pi*m_mean*m2e2s3/dzeta(3)*( &               ! terms with m2e2s3
    ! 					&	2*cc2*ii2*dzeta(3) + 2*cc1_drho*ii2_drho + 2*cc1*ii2_drho2 & ! constant terms in rho
    ! 					&	+rho_in* ( 4*cc2_drho*ii2*dzeta(3)**2 +4*ii2_drho*cc2*dzeta(3)**2 &
    ! 						& +cc1_drho2*ii2_drho +ii2_drho3*cc1 +2*cc1_drho*ii2_drho2 ) & ! linear terms in rho
    ! 					&	+(rho_in*dzeta(3))**2 *(cc2_drho2*ii2 +ii2_drho2*cc2 +2*cc2_drho*ii2_drho) ) ! quadratic terms in rho
    !
    ! 	end if

    ! ----- COMPUTATION OF Z_disp VV implementation -----
    !	zdsp  = - 2*pi*rho*edI1dz*perP.order1
    !			- pi*rho*perP.order2*m_mean*(c2_con*I2*dense + c1_con*edI2dz); ! VV code reference
    z_disp_ = -2.0_dp * pi * rho_in * ii1_drho * m2es3 &
            & - pi * rho_in * m2e2s3 * m_mean * (cc1 * ii2_drho + cc2 * ii2 * dzeta_r(3))

    !	zdspdz= zdsp/dense - 2*pi*rho*edI1d2*perP.order1
    !			- pi*rho*perP.order2*m_mean*(c3_con*I2*dense + 2*c2_con*edI2dz + c1_con*edI2d2); ! VV code reference
    z_disp_drho_ = z_disp_ / dzeta_r(3) &
            & - 2.0_dp * pi * rho_in * ii1_drho2 * m2es3 &
            & - pi * rho_in * m2e2s3 * m_mean * (cc1_drho * ii2 * dzeta_r(3) + 2.0_dp * cc2 * ii2_drho + cc1 * ii2_drho2)

    !	zdspd2 = 2/dense*zdspdz - 2/dense^2*zdsp - 2*pi*rho*perP.order1*edI1d3
    !			- pi*rho*m_mean*perP.order2*(3*c3_con*edI2dz + 3*c2_con*edI2d2 + c1_con*edI2d3 + c4_con*I2*dense ); ! VV code reference
    z_disp_drho2_ = 2.0_dp * z_disp_drho_ / dzeta_r(3) &
            & - 2.0_dp * z_disp_ / dzeta_r(3)**2 &
            & - 2.0_dp * pi * rho_in * m2es3 * ii1_drho3 &
            & - pi * rho_in * m_mean * m2e2s3 * (3.0_dp * cc1_drho * ii2_drho + 3.0_dp * cc2 * ii2_drho2 + cc1 * ii2_drho3 + cc2_drho * ii2 * dzeta_r(3))

    write (stddeb,*) "z_disp_      : ",z_disp_ ! DC DEBUG
    write (stddeb,*) "z_disp_drho_ : ",z_disp_drho_ ! DC DEBUG
    if (sec_der .eqv. .true.) then
        write (stddeb,*) "z_disp_drho2_: ",z_disp_drho2_ ! DC DEBUG
    end if
    write (stddeb,*) '----- ++++ -----'  ! DEBUG just for brevity of debug


    z_disp = z_disp_
    z_disp_drho = z_disp_drho_ * dzeta(3) ! correction for derivation in dzeta terminology
    if (sec_der .eqv. .true.) then
        z_disp_drho2 = z_disp_drho2_ * dzeta(3)**2 ! correction for derivation in dzeta terminology
    end if

    init_flag_zdisp = .true.

    write (stddeb,*) "<== z_disp_sub" !DEBUG
end subroutine z_disp_sub
