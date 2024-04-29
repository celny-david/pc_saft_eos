digraph PC_SAFT {
  /* 
   Created by graphiz online available at:
   https://dreampuf.github.io/GraphvizOnline/
   */
/* prog is blac Mdiamond
   mods are blue squares
   subs are black box
   fun are black elipse*/
   
  /*MOD SECTION*/
  concentrate=true;
  //ranksep = 0.75;
  //splines=ortho;
  //splines=polyline;
  //unflatten=true;

    /* CONTROL MOD */  
    control_mod[shape=record,color=blue,label="{contol_mod.f90| {&#8594;| {<f0> allocate_control_sub |\
    {get \n methods | { <f11> get_zres_fun | <f12> get_zres_drho_fun | <f13> get_zres_drho2_fun | <f14> get_ares_fun | <f15> get_ares_dx_fun | <f16> get_ares_dtt_fun}}| \
    {init \n methods | { <f21> set_density_sub | <f22> initialize_vars_sub | {&#8627; |<f23> set_k_ij  }| <f24> set_temperature_sub | <f25> initialize_composition_sub}}|\
    {helper \n methods | { <f31> write_array1d | <f32> write_array2d | <f33> write_array3d }}}}}"];

    control_mod:f11 -> {z_hs_sub, z_hc_sub, z_disp_sub, z_dd_sub, z_qq_sub, z_dq_sub, z_assoc_sub};
    control_mod:f12 -> {z_hs_sub, z_hc_sub, z_disp_sub, z_dd_sub, z_qq_sub, z_dq_sub, z_assoc_sub};
    control_mod:f13 -> {z_hs_sub, z_hc_sub, z_disp_sub, z_dd_sub, z_qq_sub, z_dq_sub, z_assoc_sub};
    control_mod:f14 -> {a_hs_sub, a_hc_sub, a_disp_sub, z_dd_sub, z_qq_sub, z_dq_sub, z_assoc_sub};
    control_mod:f15 -> {a_hs_sub, a_hc_sub, a_disp_sub, z_dd_sub, z_qq_sub, z_dq_sub, z_assoc_sub};
    control_mod:f16 -> {a_hs_sub, a_hc_sub, a_disp_sub, z_dd_sub, z_qq_sub, z_dq_sub, z_assoc_sub};
    
    /* CONTRIB MOD */ 
    contrib_mod[shape=record,color=blue,label="{contrib_mod.f90| {&#8594;| {<f0> allocate_contrib_sub |\
    {disp \n methods | { <f11> initialize_a0b0_disp_fun | <f12> initialize_dzeta_fun | <f13> check_dzeta_sub | <f14> initialize_combrules_fun | <f15> compute_ab_disp_fun}}| \
    {polar \n methods | { <f21> initialize_polrules_fun | { &#8627; | {<f22> initialize_a0b0c0_polar_fun | <f23> calculate_abc_polar_sub}}}}|\
    {assoc \n methods | { <f31> compute_quadratic_sub | <f32> compute_quadratic3_sub}}}}}"];
 
    /* CONTRIB MOD */ 
    param_mod[shape=record,color=lightblue,label="{param_list_mod.f90|\
    { &#8594;| {<f0> allocate_param_sub | <f1> initialize_param_sub |{ &#8627; |<f2> set_interaction_sub }}}}"];
 
    /* CONTRIB MOD */ 
    interface_mod[shape=record,color=lightblue,label="{interface_mod.f90|\
    { &#8594;| {<f0> allocate_interface_sub | <f1> allocate_all_sub |{ parser \n methods |{<f21> hardwired_input |<f22> parsed_input }}}}}"];
 
    interface_mod:f1:e -> control_mod:f0:e;
    interface_mod:f1:e -> contrib_mod:f0:e;
    interface_mod:f1:e -> param_mod:f0:e;
    interface_mod:f1:e -> interface_mod:f0:e;
    
    interface_mod:f21:e -> interface_mod:f1:e;
    interface_mod:f22:e -> interface_mod:f1:e;
    interface_mod:f22:e -> param_mod:f1:e;
  
  /* HELPER METHODS SECTION */
  rdf_sub[shape=box];
      z_assoc_sub -> rdf_sub;
      a_hc_sub -> rdf_sub;
      z_hc_sub -> rdf_sub;
  /* INTERACTION METHODS SECTION */
  // HARD SPHERES CONTRIB
  subgraph cluster_hs {
    a_hs_sub[shape=box];
    z_hs_sub[shape=box];
        
    label="hard sphere";
    //style=bold;
    color=green;
  }
    a_hs_sub -> contrib_mod:f12:e;
    a_hs_sub -> contrib_mod:f13:e;
    z_hs_sub -> contrib_mod:f12:e;
    
  // HARD CHAIN CONTRIB
  subgraph cluster_hc {
    a_hc_sub[shape=box];
    z_hc_sub[shape=box];
    
    label="hard chain";
    //style=bold;
    color=green;
  }
  
  // HARD CHAIN CONTRIB
  subgraph cluster_disp {
    a_disp_sub[shape=box];
    z_disp_sub[shape=box];
    
    label="dispersion";
    //style=bold;
    color=green;
  }
    a_disp_sub -> contrib_mod:f11:e;
    a_disp_sub -> contrib_mod:f12:e;
    a_disp_sub -> contrib_mod:f13:e;
    a_disp_sub -> contrib_mod:f14:e;
    z_disp_sub -> contrib_mod:f11:e;
    z_disp_sub -> contrib_mod:f12:e;
    z_disp_sub -> contrib_mod:f13:e;
    z_disp_sub -> contrib_mod:f14:e;
  // HARD CHAIN CONTRIB
  subgraph cluster_polar {
    z_dd_sub[shape=box];
    z_qq_sub[shape=box];
    z_dq_sub[shape=box];
    
    label="polarity";
    //style=bold;
    color=green;
  }
    z_dd_sub -> contrib_mod:f12:e;
    z_dd_sub -> contrib_mod:f21:e;
    z_qq_sub -> contrib_mod:f12:e;
    z_qq_sub -> contrib_mod:f21:e;
    z_dq_sub -> contrib_mod:f12:e;
    z_dq_sub -> contrib_mod:f21:e;
    
  // HARD CHAIN CONTRIB
  subgraph cluster_asoc {
    z_assoc_sub[shape=box];
        
    label="association";
    //style=bold;
    color=green;
  }
    z_assoc_sub -> contrib_mod:f12:e;
    z_assoc_sub -> contrib_mod:f14:e;
    z_assoc_sub -> contrib_mod:f31:e;
    z_assoc_sub -> contrib_mod:f32:e;
  
  /* CALC METHODS SECTION */
  press_calc_sub[shape=box;color=red];
    press_calc_sub -> control_mod:f11:w;
    press_calc_sub -> control_mod:f12:w;
    press_calc_sub -> control_mod:f13:w;
  chem_pot_calc_sub[shape=box;color=red];
    chem_pot_calc_sub -> control_mod:f11:w;
    chem_pot_calc_sub -> control_mod:f14:w;
    chem_pot_calc_sub -> control_mod:f15:w;
  enthalpy_entropy_calc_sub[shape=box;color=red];
    enthalpy_entropy_calc_sub -> control_mod:f11:w;
    enthalpy_entropy_calc_sub -> control_mod:f14:w;
    enthalpy_entropy_calc_sub -> control_mod:f16:w;
  helmholtz_calc_sub[shape=box;color=red];
    helmholtz_calc_sub -> control_mod:f14:w;
  num_diff_calc_sub[shape=box;color=red];
    num_diff_calc_sub -> control_mod:f11:w;
    num_diff_calc_sub -> control_mod:f12:w;
    num_diff_calc_sub -> control_mod:f14:w;
    num_diff_calc_sub -> control_mod:f15:w;
    num_diff_calc_sub -> control_mod:f16:w;
    
  /* PROG SECTION */
  pcsaft_main_prog [shape=Mdiamond;style=filled;fillcolor=magenta;fontsize=20];
    pcsaft_main_prog -> interface_mod:f21:w;
    pcsaft_main_prog -> interface_mod:f22:w;
    pcsaft_main_prog -> press_calc_sub;
    pcsaft_main_prog -> chem_pot_calc_sub;
    pcsaft_main_prog -> enthalpy_entropy_calc_sub;
    pcsaft_main_prog -> helmholtz_calc_sub;
    pcsaft_main_prog -> num_diff_calc_sub;
  
  
  
  label = "PC_SAFT call structure"
}
