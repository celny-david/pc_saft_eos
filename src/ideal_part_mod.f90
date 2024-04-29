module ideal_part_mod
    ! this module is part of individual piece of code responsible
    ! for calculation of ideal properties of substances
    ! mod is developed as independent on pc-saft eos with get methods for simple data withdrawal
    ! provided interface can be directly used in pc-saft eos (primary target)
    ! or further interfaced as independent program
    ! important notes and comments
    implicit none

    !	=== variable section ===
    private
    ! beware the sp, dp, qp, ch40_length parameters are mirrored from controll module
    integer, parameter :: sp = selected_real_kind(6, 37)
    integer, parameter :: dp = selected_real_kind(15, 307)
    integer, parameter :: qp = selected_real_kind(33, 4931)
    integer, parameter :: ch40_length = 40 ! number of characters for names of components
    ! beware n_comp parameter is taken from controll module
    ! TODO should be removed and replaced by allocable approach
    integer, parameter :: n_comp = 1 ! number of components

    !	=== variable section ===
    ! module specific controll variables
    real(dp) :: tt ! temperature
    real(dp), dimension(n_comp) :: mol_rat ! the imputed molar ratio
    character(len = ch40_length), dimension(n_comp) :: names    ! names input
    integer :: output_flag ! the mode of output

    ! temporary storage for the mixture computation purposes
    real(dp), dimension(n_comp) :: molar_mass_comp

    real(dp), dimension(n_comp) :: helmholtz_id_comp
    real(dp), dimension(n_comp) :: enthalpy_id_comp
    real(dp), dimension(n_comp) :: entropy_id_comp
    real(dp), dimension(n_comp) :: innerenergy_id_comp
    real(dp), dimension(n_comp) :: gibbs_id_comp
    real(dp), dimension(n_comp, n_comp) :: chem_pot_id_comp

    real(dp), dimension(n_comp) :: c_v_id_comp
    real(dp), dimension(n_comp) :: c_p_id_comp

    !list the public variables
    real(dp), public :: helmholtz_id
    real(dp), public :: enthalpy_id
    real(dp), public :: entropy_id
    real(dp), public :: innerenergy_id
    real(dp), public :: gibbs_id
    real(dp), public :: c_v_id
    real(dp), public :: c_p_id
    real(dp), dimension(n_comp), public :: chem_pot_id

    !list the public methods
    public :: set_output_flag_sub
    public :: set_properties_sub
    public :: calculate_ideal_sub

contains
    subroutine set_output_flag_sub(mode_in)
        ! set the mode of output
        ! intended as intermediate interface preferred in pc-saft
        implicit none
        ! variable section
        integer, intent(in) :: mode_in

        output_flag = mode_in

    end subroutine set_output_flag_sub

    subroutine set_properties_sub(tt_in, mol_rat_in, names_in)
        ! set the required properties into this module
        ! intended as intermediate interface common for both pc-saft and standalone module
        implicit none
        ! variable section
        real(dp), intent(in) :: tt_in ! temperature
        real(dp), dimension(n_comp), intent(in) :: mol_rat_in ! the imputed molar ratio
        character(len = ch40_length), dimension(n_comp), intent(in) :: names_in! names input

        ! initialisation
        tt = tt_in
        mol_rat = mol_rat_in
        names = names_in
        output_flag = 1 ! default value

    end subroutine set_properties_sub

    ! ================================ CALCULATE SECTION =====================================
    subroutine calculate_ideal_sub(rho_in)
        ! method for calculation of the ideal properties relaying the computation to individual methods specific for each substance
        ! NOTE TODO check if possible to create pointer to function to avoid the massive if statement
        implicit none
        ! variable section
        real(dp), intent(in) :: rho_in ! temperature
        ! local variables
        character(len = ch40_length) :: tmp_name ! temporary name for if statements
        real(dp) :: mm ! the total molar mass M
        real(dp) :: w ! the temporary molar fraction
        integer :: i ! counter

        ! initialization
        molar_mass_comp = 0.0_dp
        helmholtz_id = 0.0_dp
        enthalpy_id = 0.0_dp
        entropy_id = 0.0_dp
        innerenergy_id = 0.0_dp
        c_v_id = 0.0_dp
        c_p_id = 0.0_dp

        ! component ideal properties computation
        do i = 1, n_comp
            tmp_name = trim(adjustl(names(i)))
            if (tmp_name=='carbon dioxide' .or. tmp_name=='carbon dioxide NP') then
                call co2_ideal_sub(rho_in, i)
            elseif (tmp_name=='methane') then
                call ch4_ideal_sub(rho_in, i)
            elseif (tmp_name=='nitrogen') then
                call n2_ideal_sub(rho_in, i)
            else ! default action
                write (*, *) 'Error in subroutine calculate_ideal_sub: the name "', tmp_name, '" is not implemented.'
                stop
            end if

        end do

        ! mixture properties calculation

        do i = 1, n_comp
            helmholtz_id = helmholtz_id + helmholtz_id_comp(i) * mol_rat(i)
            enthalpy_id = enthalpy_id + enthalpy_id_comp(i) * mol_rat(i)
            entropy_id = entropy_id + entropy_id_comp(i) * mol_rat(i)
            innerenergy_id = innerenergy_id + innerenergy_id_comp(i) * mol_rat(i)
            gibbs_id = gibbs_id + gibbs_id_comp(i) * mol_rat(i)
            c_v_id = c_v_id + c_v_id_comp(i) * mol_rat(i)
            c_p_id = c_p_id + c_p_id_comp(i) * mol_rat(i)
        end do

    end subroutine calculate_ideal_sub
    ! ================================ INDIVIDUAL METHODS ====================================
    ! ---- ---- CO2 ---- ----
    subroutine co2_ideal_sub(rho_in, ii)
        ! Authors of paper: U. Setzmann and W. Wagner
        ! Paper title:      A New Equation of State and Tables of Thermodynamic Properties for
        !                   Methane Covering the Range from the Melting Line to 625 K at Pressures up to 1000 MPa
        ! Journal:          J. Phys. Chem. Ref. Data, Vol. 20, No. 6, 1991, pp. 1061-1155
        implicit none
        ! variable section
        real(dp), intent(in) :: rho_in ! density in [kg/m^3]
        integer, intent(in) :: ii ! counter to which submit the ideal part in case of mixtures
        ! local variables
        real(dp), parameter :: mm = 0.0440098_dp ! molar mass in [kg/mol]
        real(dp), parameter :: rr = 188.9241_dp ! R_gas in [J/kg*K]
        real(dp), parameter :: tt_crit = 304.1282_dp ! in [K]
        real(dp), parameter :: rho_crit = 467.6_dp ! in [kg*m^-3]
        real(dp), dimension(2, 8) :: eq_par ! the parameters of ideal equation for co2
        real(dp) :: delta, tau ! critically scaled values
        real(dp) :: phi_0, phi_0_dt, phi_0_dt2 ! phi0 intermediate parameters
        integer :: i ! counter

        molar_mass_comp(ii) = mm

        ! parameter initialisation
        eq_par(:, 1) = (/ 8.37304456_dp, 0.0_dp/)
        eq_par(:, 2) = (/-3.70454304_dp, 0.0_dp/)
        eq_par(:, 3) = (/ 2.50000000_dp, 0.0_dp/)
        eq_par(:, 4) = (/ 1.99427042_dp, 3.15163_dp/)
        eq_par(:, 5) = (/ 0.62105248_dp, 6.11190_dp/)
        eq_par(:, 6) = (/ 0.41195293_dp, 6.77708_dp/)
        eq_par(:, 7) = (/ 1.04028922_dp, 11.32384_dp/)
        eq_par(:, 8) = (/ 0.08327678_dp, 27.08792_dp/)

        delta = rho_in / rho_crit
        tau = tt_crit / tt

        ! based on Eq. (6.3), Table 28 from paper
        if (delta==0.0_dp) then
            phi_0 = 0.0_dp
            write(*, *) 'Warning: zero input density -related properties are inconsistent'
        else
            phi_0 = log(delta) + eq_par(1, 1) + eq_par(1, 2) * tau + eq_par(1, 3) * log(tau)
        end if
        phi_0_dt = eq_par(1, 2) + eq_par(1, 3) / tau
        phi_0_dt2 = -eq_par(1, 3) / tau**2
        do i = 4, 8
            phi_0 = phi_0 + eq_par(1, i) * log(1.0_dp - exp(-tau * eq_par(2, i)))
            phi_0_dt = phi_0_dt + eq_par(1, i) * eq_par(2, i) * (1.0_dp / (1.0_dp - exp(-tau * eq_par(2, i))) - 1.0_dp)
            phi_0_dt2 = phi_0_dt2 - eq_par(1, i) * eq_par(2, i)**2 * exp(-tau * eq_par(2, i)) / (1.0_dp - exp(-tau * eq_par(2, i)))**2
        end do

        ! based on Table 3 from paper

        ! modified for conversion to J/mol
        helmholtz_id_comp(ii) = phi_0 * rr * tt * mm
        enthalpy_id_comp(ii) = rr * tt * (1.0_dp + tau * phi_0_dt) * mm
        entropy_id_comp(ii) = rr * (tau * phi_0_dt - phi_0) * mm
        innerenergy_id_comp(ii) = rr * tt * tau * phi_0_dt * mm
        gibbs_id_comp(ii) = enthalpy_id_comp(ii) - tt * entropy_id_comp(ii)
        c_v_id_comp(ii) = -rr * tau**2 * phi_0_dt2 * mm
        c_p_id_comp(ii) = c_v_id_comp(ii) + rr * mm

        if (output_flag == 1) then
            write(*, *) '=== IDEAL MODULE ==='
            write(*, "('helmholtz(',i2,') = ',e16.8,' [J/mol]')") ii, helmholtz_id_comp(ii)
            write(*, "('entropy(',i2,') = ',e16.8,' [J/(mol*K)]')") ii, entropy_id_comp(ii)
            write(*, "('inner energy(',i2,') = ',e16.8,' [J/mol]')") ii, innerenergy_id_comp(ii)
            write(*, "('gibbs energy(',i2,') = ',e16.8,' [J/mol]')") ii, gibbs_id_comp(ii)
            write(*, "('c_v(',i2,') = ',e16.8,' [J/(mol*K)]')") ii, c_v_id_comp(ii)
            write(*, "('c_p(',i2,') = ',e16.8,' [J/(mol*K)]')") ii, c_p_id_comp(ii)
        end if

    end subroutine co2_ideal_sub

    ! ---- ---- CH4 ---- ----
    subroutine ch4_ideal_sub(rho_in, ii)
        ! Authors of paper: R. Span and W. Wagner
        ! Paper title:      A New Equation of State for Carbon Dioxide Covering the Fluid Region
        !                   from the Triple-Point Temperature to 1100 K at Pressures up to 800 MPa
        ! Journal:          J. Phys. Chern. Ref. Data, Vol. 25, No.6, 1996, pp. 1509 -- 1596
        implicit none
        ! variable section
        real(dp), intent(in) :: rho_in ! density in [kg/m^3]
        integer, intent(in) :: ii ! counter to which submit the ideal part in case of mixtures
        ! local variables
        real(dp), parameter :: mm = 0.0160428_dp ! molar mass in [kg/mol]
        real(dp), parameter :: rr = 518.2705_dp ! R_gas in [J/kg*K]
        real(dp), parameter :: tt_crit = 190.564_dp ! in [K]
        real(dp), parameter :: rho_crit = 162.66_dp ! in [kg*m^-3]
        real(dp), dimension(2, 8) :: eq_par ! the parameters of ideal equation for co2
        real(dp) :: delta, tau ! critically scaled values
        real(dp) :: phi_0, phi_0_dt, phi_0_dt2 ! phi0 intermediate parameters
        integer :: i ! counter

        molar_mass_comp(ii) = mm

        ! parameter initialisation
        eq_par(:, 1) = (/ 9.91243972_dp, 0.0_dp/)
        eq_par(:, 2) = (/-6.33270087_dp, 0.0_dp/)
        eq_par(:, 3) = (/ 3.0016_dp, 0.0_dp/)
        eq_par(:, 4) = (/ 0.008449_dp, 3.40043240_dp/)
        eq_par(:, 5) = (/ 4.6942_dp, 10.26951575_dp/)
        eq_par(:, 6) = (/ 3.4865_dp, 20.43932747_dp/)
        eq_par(:, 7) = (/ 1.6572_dp, 29.93744884_dp/)
        eq_par(:, 8) = (/ 1.4115_dp, 79.13351945_dp/)

        delta = rho_in / rho_crit
        tau = tt_crit / tt

        ! based on Eq. (5.2), Table 34 from paper
        if (delta==0.0_dp) then
            phi_0 = 0.0_dp
            write(*, *) 'Warning: zero input density -related properties are inconsistent'
        else
            phi_0 = log(delta) + eq_par(1, 1) + eq_par(1, 2) * tau + eq_par(1, 3) * log(tau)
        end if
        phi_0_dt = eq_par(1, 2) + eq_par(1, 3) / tau
        phi_0_dt2 = -eq_par(1, 3) / tau**2
        do i = 4, 8
            phi_0 = phi_0 + eq_par(1, i) * log(1.0_dp - exp(-tau * eq_par(2, i)))
            phi_0_dt = phi_0_dt + eq_par(1, i) * eq_par(2, i) * (1.0_dp / (1.0_dp - exp(-tau * eq_par(2, i))) - 1.0_dp)
            phi_0_dt2 = phi_0_dt2 - eq_par(1, i) * eq_par(2, i)**2 * exp(-tau * eq_par(2, i)) / (1.0_dp - exp(-tau * eq_par(2, i)))**2
        end do

        ! based on Table 3 from paper

        ! modified for conversion to J/mol
        helmholtz_id_comp(ii) = phi_0 * rr * tt * mm
        enthalpy_id_comp(ii) = rr * tt * (1.0_dp + tau * phi_0_dt) * mm
        entropy_id_comp(ii) = rr * (tau * phi_0_dt - phi_0) * mm
        innerenergy_id_comp(ii) = rr * tt * tau * phi_0_dt * mm
        gibbs_id_comp(ii) = enthalpy_id_comp(ii) - tt * entropy_id_comp(ii)
        c_v_id_comp(ii) = -rr * tau**2 * phi_0_dt2 * mm
        c_p_id_comp(ii) = c_v_id_comp(ii) + rr * mm

        if (output_flag == 1) then
            write(*, *) '=== IDEAL MODULE ==='
            write(*, "('helmholtz(',i2,') = ',e16.8,' [J/mol]')") ii, helmholtz_id_comp(ii)
            write(*, "('entropy(',i2,') = ',e16.8,' [J/(mol*K)]')") ii, entropy_id_comp(ii)
            write(*, "('inner energy(',i2,') = ',e16.8,' [J/mol]')") ii, innerenergy_id_comp(ii)
            write(*, "('gibbs energy(',i2,') = ',e16.8,' [J/mol]')") ii, gibbs_id_comp(ii)
            write(*, "('c_v(',i2,') = ',e16.8,' [J/(mol*K)]')") ii, c_v_id_comp(ii)
            write(*, "('c_p(',i2,') = ',e16.8,' [J/(mol*K)]')") ii, c_p_id_comp(ii)
        end if

    end subroutine ch4_ideal_sub

    ! ---- ---- N2 ---- ----
    subroutine n2_ideal_sub(rho_in, ii)
        ! Authors of paper: R. Span, E. W. Lemmon, R. T. Jcobsen, and W. Wagner
        ! Paper title:      A reference quality equation of state for nitrogen.
        ! Journal:          Int. J. of Thermophysics, vol. 19 No. 4 (1998), 1121-1132
        implicit none
        ! variable section
        real(dp), intent(in) :: rho_in ! density in [kg/m^3]
        integer, intent(in) :: ii ! counter to which submit the ideal part in case of mixtures
        ! local variables
        real(dp), parameter :: mm = 0.02801348_dp ! molar mass in [kg/mol]
        real(dp), parameter :: rr = 8.31451_dp / mm ! R_gas in [J/kg*K]
        real(dp), parameter :: tt_crit = 126.192_dp ! in [K]
        real(dp), parameter :: rho_crit = 11183.9_dp * mm ! in [kg*m^-3]
        real(dp), dimension(1, 7) :: eq_par ! the parameters of ideal equation for co2
        real(dp) :: c = 26.65788_dp ! constant parameters
        real(dp) :: delta, tau ! critically scaled values
        real(dp) :: alpha_0, alpha_1, alpha_0_dt, alpha_0_dt2 ! alpha intermediate parameters

        molar_mass_comp(ii) = mm

        ! parameter initialisation
        ! for format consistency in this form
        eq_par(:, 1) = (/ 2.5_dp        /)
        eq_par(:, 2) = (/-12.76953_dp   /)
        eq_par(:, 3) = (/-0.007841630_dp/)
        eq_par(:, 4) = (/-1.934819e-4_dp/)
        eq_par(:, 5) = (/-1.247742e-5_dp/)
        eq_par(:, 6) = (/ 6.678326e-8_dp/)
        eq_par(:, 7) = (/ 1.012941_dp   /)

        delta = rho_in / rho_crit
        tau = tt_crit / tt

        if (delta==0.0_dp) then
            alpha_0 = 0.0_dp
            write(*, *) 'Warning: zero input density -related properties are inconsistent'
        else
            alpha_1 = + eq_par(1, 1) * log(tt) &
                    & + eq_par(1, 2) &
                    & + eq_par(1, 3) * tt &
                    & + eq_par(1, 4) / tt &
                    & + eq_par(1, 5) / tt**2 &
                    & + eq_par(1, 6) / tt**3 &
                    & + eq_par(1, 7) * log(exp(c * tt) - 1.0_dp) - c * tt

            alpha_0 = log(delta) + alpha_1
        end if
        alpha_0_dt = + eq_par(1, 1) / tt &
                & + eq_par(1, 3) &
                & - eq_par(1, 4) / tt**2 &
                & - 2.0_dp * eq_par(1, 5) / tt**3 &
                & - 3.0_dp * eq_par(1, 6) / tt**4 &
                & + eq_par(1, 7) * c / (exp(c * tt) - 1.0_dp)
        alpha_0_dt2 = - eq_par(1, 1) / tt**2 &
                & + 2.0_dp * eq_par(1, 4) / tt**3 &
                & + 6.0_dp * eq_par(1, 5) / tt**4 &
                & + 12.0_dp * eq_par(1, 6) / tt**5 &
                & - eq_par(1, 7) * c * c * exp(c * tt) / (exp(c * tt) - 1.0_dp)**2

        ! based on Table 3 from paper

        ! modified for conversion to J/mol
        helmholtz_id_comp(ii) = alpha_0 * rr * tt * mm
        enthalpy_id_comp(ii) = rr * tt * (1.0_dp + tau * alpha_0_dt) * mm
        entropy_id_comp(ii) = rr * (tau * alpha_0_dt - alpha_0) * mm
        innerenergy_id_comp(ii) = rr * tt * tau * alpha_0_dt * mm
        gibbs_id_comp(ii) = enthalpy_id_comp(ii) - tt * entropy_id_comp(ii)
        c_v_id_comp(ii) = -rr * tau**2 * alpha_0_dt2 * mm
        c_p_id_comp(ii) = c_v_id_comp(ii) + rr * mm

        if (output_flag == 1) then
            write(*, *) '=== IDEAL MODULE ==='
            write(*, "('helmholtz(',i2,') = ',e16.8,' [J/mol]')") ii, helmholtz_id_comp(ii)
            write(*, "('entropy(',i2,') = ',e16.8,' [J/(mol*K)]')") ii, entropy_id_comp(ii)
            write(*, "('inner energy(',i2,') = ',e16.8,' [J/mol]')") ii, innerenergy_id_comp(ii)
            write(*, "('gibbs energy(',i2,') = ',e16.8,' [J/mol]')") ii, gibbs_id_comp(ii)
            write(*, "('c_v(',i2,') = ',e16.8,' [J/(mol*K)]')") ii, c_v_id_comp(ii)
            write(*, "('c_p(',i2,') = ',e16.8,' [J/(mol*K)]')") ii, c_p_id_comp(ii)
        end if

    end subroutine n2_ideal_sub
end module ideal_part_mod

! =========================================================================
! =========================================================================
! =========================================================================
! function [Consts,Props]=CO2WagnerSpan_Id1(T,rho)
! %R. Span and W. Wagner
! %A New Equation of State for Carbon Dioxide Covering the Fluid Region
! %from the Triple-Point Temperature to 1100 K at Pressures up to 800 MPa
! %J. Phys. Chern. Ref. Data, Vol. 25, No.6, 1996, pp. 1509 -- 1596
! 
! Props={}; Consts={};
! M=44.0098; %kg/kmol
! R=188.9241; %J/(kg*K)
! Tc=304.1282; %K
! rhoc=467.6; %kg/m^3
! Consts.M=M; Consts.R=R; Consts.Tc=Tc;  Consts.rhoc=rhoc;
! if isempty(T)
!     return
! end
! 
! tab27=[
! 1	8.37304456	NaN;
! 2	-3.70454304	NaN;
! 3	2.50000000	NaN;
! 4	1.99427042	3.15163;
! 5	0.62105248	6.11190;
! 6	0.41195293	6.77708;
! 7	1.04028922	11.32384;
! 8	0.08327678	27.08792];
! 
! a0=tab27(:,2);
! th0=tab27(:,3);
! 
! delta=rho/rhoc;
! tau=Tc./T;
! RT=R*T;
! 
! zt=zeros(size(T));
! zd=zeros(size(rho));
! zz=zt.*zd;
! 
! %Eq. (6.3), Table 28
! phi0=log(delta)+a0(1)+a0(2)*tau+a0(3)*log(tau);
! phi0t=a0(2)+a0(3)./tau;
! phi0tt=-a0(3)./tau.^2;
! for i=4:8
!     phi0=phi0+a0(i)*log(1-exp(-tau*th0(i)));
!     phi0t=phi0t+a0(i)*th0(i)*((1-exp(-tau*th0(i))).^(-1)-1);
!     phi0tt=phi0tt-a0(i)*th0(i)^2*exp(-tau*th0(i)).*(1-exp(-tau*th0(i))).^(-2);
! end
! 
! phi=phi0;
! 
! %table3
! f=phi.*RT;
! p=rho.*RT;
! s=R*(tau.*phi0t-phi0);
! u=RT.*tau.*phi0t;
! cv=-R*tau.^2.*phi0tt;
! h=RT.*(1+tau.*phi0t);
! cp=cv+R;
! w=1-1./tau.^2./phi0tt;
! w=(RT.*w).^0.5;
! kJT=0;
! prho=RT;
! pT=R*rhoc*delta;
! 
! Props.f=f;
! Props.p=p;
! Props.s=s;
! Props.u=u;
! Props.cv=cv;
! Props.h=h;
! Props.cp=cp;
! Props.w=w;
! Props.kJT=kJT;
! Props.prho=prho;
! Props.pT=pT;
! 
! end

! =========================================================================

! function [Consts,Props]=CH4id1(T,rho)
! %U. Setzmann and W. Wagner
! %A New Equation of State and Tables of Thermodynamic Properties for
! %Methane Covering the Range from the Melting Line to 625 K
! %at Pressures up to 1000 MPa
! %J. Phys. Chem. Ref. Data, Vol. 20, No. 6, 1991, pp. 1061-1155
! Props={}; Consts={};
! M=16.0428; %kg/kmol
! R=518.2705; %J/(kg*K)
! Tc=190.564; %K
! rhoc=162.66; %kg/m^3
! Consts.M=M; Consts.R=R; Consts.Tc=Tc;  Consts.rhoc=rhoc;
! Consts.pc=4.5922e6; %Pa
! Consts.Tt=90.6941;  %K
! Consts.pt=0.011696; %Pa
! Consts.Tb=111.668;  %K
! 
! Table33=[
!     1   9.91243972  NaN         
!     2   -6.33270087 NaN         
!     3   3.0016      NaN         
!     4   0.008449    3.40043240  
!     5   4.6942      10.26951575;
!     6   3.4865      20.43932747;
!     7   1.6572      29.93744884;
!     8   1.4115      79.13351945];
! a=Table33(:,2);
! th=Table33(:,3);
! 
! delta=rho/rhoc;
! tau=Tc./T;
! RT=R*T;
! 
! %Eq. (5.2), Table 34
! phi0=log(delta)+a(1)+a(2)*tau+a(3)*log(tau);
! phi0t=a(2)+a(3)./tau;
! phi0tt=-a(3)./tau.^2;
! for i=4:8
!     phi0=phi0+a(i)*log(1-exp(-tau*th(i)));
!     phi0t=phi0t+a(i)*th(i)*((1-exp(-tau*th(i))).^(-1)-1);
!     phi0tt=phi0tt-a(i)*th(i)^2*exp(-tau*th(i)).*(1-exp(-tau*th(i))).^(-2);
! end
! 
! %table28
! f=phi0.*RT;
! p=rho.*RT;
! s=R*(tau.*(phi0t)-phi0);
! u=RT.*tau.*(phi0t);
! cv=-R*tau.^2.*(phi0tt);
! h=RT.*(1+tau.*(phi0t));
! cp=cv+R;
! w=1-1./tau.^2./phi0tt;
! w=(RT.*w).^0.5;
! 
! Props.f=f;
! Props.p=p;
! Props.s=s;
! Props.u=u;
! Props.cv=cv;
! Props.h=h;
! Props.cp=cp;
! Props.w=w;
! 
! end

! =========================================================================

! function props = N2Id1(T,rho)
! %K, kg.m^-3
! % {R. Span, E. W. Lemmon, R. T. Jcobsen, and W. Wagner:
! %  A reference quality equation of state for nitrogen.
! %  Int. J. of Thermophysics, vol. 19 No. 4 (1998), 1121-1132}
! 
! Tt=63.151;
! Tc = 126.192;  %K
! pc = 3.3958e6; %Pa
! rhocM = 11.1839; %kmol / m^3
! R = 8314.51; %J / (kmol K)
! MM = 28.01348; %kg / (kmol K);
! r = R/MM; %J / (kg K)
! rhoc=rhocM*MM; %kg / m^3 
! 
! d=rho./rhoc;
! t=Tc./T;
! 
! M=[2.5;-12.76953;-0.007841630;-1.934819e-4;-1.247742e-5;6.678326e-8;1.012941];
! c = 26.65788;
! ee=exp(c*t);
! alph01 = M(1)*log(t)+M(2)+M(3)*t+M(4)*t.^(-1)+M(5)*t.^(-2)+M(6)*t.^(-3)+M(7)*log((ee-1)./ee);
! warning('off');
! alph0 =alph01+log(d);
! warning('on');
! alph0t =M(1)*t.^(-1)+M(3)-M(4)*t.^(-2)-2*M(5)*t.^(-3)-3*M(6)*t.^(-4)+M(7)*c./(ee-1); 
! alph0tt =-M(1)*t.^(-2)+2*M(4)*t.^(-3)+6*M(5)*t.^(-4)+12*M(6)*t.^(-5)-M(7)*c^2*ee./(ee-1).^2; 
! 
! rT=r*T;
! props.T=T;
! props.rho=rho; %for iteration
! props.f=rT.*(alph0);
! props.g=rT.*(1+alph0);
! props.p=rho.*rT;
! props.u=rT.*t.*alph0t;
! props.s=r*(t.*(alph0t)-alph0);
! props.sR=r*(t.*(alph0t)-alph01); %sR=s+log(rho/rhoc)
! props.h=rT.*(1+t.*(alph0t));
! props.cv=r*(-t.^2.*(alph0tt));
! props.cp=r*(-t.^2.*(alph0tt)+1);
! props.w=(rT.*(1-1./t.^2./alph0tt)).^0.5;
! props.p_rho=rT;
! props.p_T=r*rho;
! 
! props.r=r;
! props.M=MM;
! props.Tc=Tc;
! props.rhoc=rhoc;
! props.pc=pc;
! props.Tt=Tt;
! 
! end
