
#ifdef __cplusplus
extern "C" {
#endif /* __cplusplus */

void parse_input();
void set_n_comp(int n_comp_in);
void set_component_name(char * fluid_name_in, int str_length);
void set_density(double density_in );
void set_temperature(double temperature_in );
void set_molar_ratio(double* mol_rat_in );

void press_calc(double *p_out, double *p_drho_out, double *p_drho2_out);
void chem_pot_calc(double* mu_out, double*  fugacity_out);
void enthalpy_entropy_calc(double* h_out, double* s_out, double* g_out);
void all_calc(double *p_out, double *p_drho_out, double *p_drho2_out,
              double* mu_out, double*  fugacity_out,
              double* h_out, double* s_out, double* g_out);

#ifdef __cplusplus
}
#endif /* __cplusplus */
