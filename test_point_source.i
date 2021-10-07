# species=0 is water
# species=1 is co2
# phase=0 is liquid, and since massfrac_ph0_sp0 = 1, this is all water
# phase=1 is gas, and since massfrac_ph1_sp0 = 0, this is all co2
#
[Mesh]
  type = FileMesh
  file ="mesh_coarse.msh"
[]

[GlobalParams]
  displacements = 'disp_x disp_y'
  PorousFlowDictator = dictator
  gravity = '0 0 0'
  biot_coefficient = 1.0
[]

[Variables]
  [pwater]
  []
  [sgas]
    initial_condition = 0.0
  []
  [disp_x]
  []
  [disp_y]
  []
[]

[ICs]
  [pwater]
    type = FunctionIC
    function = p_hydro
    variable = pwater
  []
[]

[AuxVariables]
  [massfrac_ph0_sp0]
    initial_condition = 1 # all H20 in phase=0
  []
  [massfrac_ph1_sp0]
    initial_condition = 0 # no H2O in phase=1
  []
  [pgas]
    family = MONOMIAL
    order = FIRST
  []
  [swater]
    family = MONOMIAL
    order = FIRST
  []
  [stress_xx]
    order = CONSTANT
    family = MONOMIAL
  []
  [stress_yy]
    order = CONSTANT
    family = MONOMIAL
  []
[]

[Kernels]
  [mass_water_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 0
    variable = pwater
  []
  [flux_water]
    type = PorousFlowAdvectiveFlux
    fluid_component = 0
    use_displaced_mesh = false
    variable = pwater
  []
  [mass_co2_dot]
    type = PorousFlowMassTimeDerivative
    fluid_component = 1
    variable = sgas
  []
  [flux_co2]
    type = PorousFlowAdvectiveFlux
    fluid_component = 1
    use_displaced_mesh = false
    variable = sgas
  []
  [grad_stress_x]
    type = StressDivergenceTensors
    variable = disp_x
    use_displaced_mesh = false
    component = 0
  []
  [grad_stress_y]
    type = StressDivergenceTensors
    variable = disp_y
    use_displaced_mesh = false
    component = 1
  []
  [poro_x]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_x
    use_displaced_mesh = false
    component = 0
  []
  [poro_y]
    type = PorousFlowEffectiveStressCoupling
    variable = disp_y
    use_displaced_mesh = false
    component = 1
  []
  [vol_exp_co2]
    type = PorousFlowMassVolumetricExpansion
    fluid_component = 1
    variable = sgas
  []
  [vol_exp_water]
    type = PorousFlowMassVolumetricExpansion
    fluid_component = 0
    variable = pwater
  []
[]

[AuxKernels]
  [pgas]
    type = PorousFlowPropertyAux
    property = pressure
    phase = 1
    variable = pgas
  []
  [swater]
    type = PorousFlowPropertyAux
    property = saturation
    phase = 0
    variable = swater
  []
  [stress_xx]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_xx
    index_i = 0
    index_j = 0
  []
  [stress_yy]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_yy
    index_i = 1
    index_j = 1
  []
[]


[UserObjects]
  [dictator]
    type = PorousFlowDictator
    porous_flow_vars = 'pwater sgas disp_x disp_y'
    number_fluid_phases = 2
    number_fluid_components = 2
  []
  [pc]
    type = PorousFlowCapillaryPressureConst
    pc = 0
  []
[]

[Modules]
  [FluidProperties]
    [water]
      type = SimpleFluidProperties
      bulk_modulus = 2.27e14
      density0 = 970.0
      viscosity = 0.3394e-3
      cv = 4149.0
      cp = 4149.0
      porepressure_coefficient = 0.0
      thermal_expansion = 0
    []
    [co2]
      type = SimpleFluidProperties
      bulk_modulus = 2.27e14
      density0 = 516.48
      viscosity = 0.0393e-3
      cv = 2920.5
      cp = 2920.5
      porepressure_coefficient = 0.0
      thermal_expansion = 0
    []
  []
[]

[Materials]
  [temperature]
    type = PorousFlowTemperature
    temperature = 60
  []
  [ppss]
    type = PorousFlow2PhasePS
    phase0_porepressure = pwater
    phase1_saturation = sgas
    capillary_pressure = pc
  []
  [massfrac]
    type = PorousFlowMassFraction
    mass_fraction_vars = 'massfrac_ph0_sp0 massfrac_ph1_sp0'
  []
  [water]
    type = PorousFlowSingleComponentFluid
    fp = water
    phase = 0
  []
  [gas]
    type = PorousFlowSingleComponentFluid
    fp = co2
    phase = 1
  []
  [porosity_reservoir]
    type =   PorousFlowPorosity
    fluid = true
    mechanical = true
#    thermal = true
    porosity_zero = 0.1
    reference_temperature = 330
    reference_porepressure = 10E6
#    thermal_expansion_coeff = 15E-6 # volumetric
    solid_bulk = 8E9 # unimportant since biot = 1
  []
  [permeability_reservoir]
#    type = PorousFlowPermeabilityKozenyCarman
    type = PorousFlowPermeabilityConst
    permeability = '2e-12 0 0  0 0 0  0 0 0'
  []
  [relperm_liquid]
    type = PorousFlowRelativePermeabilityCorey
    n = 4
    phase = 0
    s_res = 0.200
    sum_s_res = 0.405
  []
  [relperm_gas]
    type = PorousFlowRelativePermeabilityBC
    phase = 1
    s_res = 0.205
    sum_s_res = 0.405
    nw_phase = true
    lambda = 2
  []
  [elasticity_tensor]
    type = ComputeIsotropicElasticityTensor
    shear_modulus = 6.0E9
    poissons_ratio = 0.2
  []
  [strain]
    type = ComputeSmallStrain
    eigenstrain_names = 'ini_stress'
  []
  [ini_strain]
    type = ComputeEigenstrainFromInitialStress
    initial_stress = 'sigma_v_ini 0 0 0 sigma_h_ini 0 0 0 sigma_h_ini'
    eigenstrain_name = ini_stress
  []
#  [thermal_contribution]
#    type = ComputeThermalExpansionEigenstrain
#    temperature = temp
#    stress_free_temperature = 358
#    thermal_expansion_coeff = 5E-6
#    eigenstrain_name = thermal_contribution
#  []
  [stress]
    type = ComputeLinearElasticStress
#    type = ComputeFiniteStrainElasticStress
  []
  [eff_fluid_pressure]
    type = PorousFlowEffectiveFluidPressure
  []
  [vol_strain]
    type = PorousFlowVolumetricStrain
  []
[]

[Functions]
  [sigma_v_ini]
    type = ParsedFunction
    value = '0.1e6 - 2260*9.91 * y'        # y is -ve 
  []
  [sigma_h_ini]
    type = ParsedFunction
    value = '0.6*(0.1e6 - 2260*9.81*y)' # = 0.6*sigma_v_ini
  []
  [p_hydro]
    type = ParsedFunction
    value = '0.1e6 -1000*9.81*y' 
  []
[]

[BCs]
  # Bottom: No Flow and fixed uy  
  [fixed_y]
    type = DirichletBC
    variable = disp_y
    value = 0
    boundary = bottom
  []

  # Left & Right
  [outer_pressure_fixed]
    type = DirichletBC
    boundary = right
    value = 15.3e6      # TODO: Apply hydrostatic pressure
    variable = pwater
  []
  [outer_saturation_fixed]
    type = DirichletBC
    boundary = right
    value = 0.0
    variable = sgas
  []

  # Top:  No-flow and stress 
  [top]
    type = Pressure
    variable = disp_y
    boundary = 'top'
    component = 1
    factor = 0.1e6  # 0.1 MPa
  []
[]

[DiracKernels]
  [sink]
    type = PorousFlowSquarePulsePointSource
    point = '0 -1470 0'
    mass_flux = 0.05
    variable = sgas
  []
[]

[Preconditioning]
  active = 'smp'
  [smp]
    type = SMP
    full = true
    #petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -sub_pc_type -sub_pc_factor_shift_type -pc_asm_overlap -snes_atol -snes_rtol -snes_max_it'
    petsc_options_value = 'gmres      asm      lu           NONZERO                   2               1E2       1E-5        500'
  []
  [mumps]
    type = SMP
    full = true
    petsc_options = '-snes_converged_reason -ksp_diagonal_scale -ksp_diagonal_scale_fix -ksp_gmres_modifiedgramschmidt -snes_linesearch_monitor'
    petsc_options_iname = '-ksp_type -pc_type -pc_factor_mat_solver_package -pc_factor_shift_type -snes_rtol -snes_atol -snes_max_it'
    petsc_options_value = 'gmres      lu       mumps                         NONZERO               1E-5       1E2       50'
  []
[]

[Executioner]
  type = Transient
  solve_type = NEWTON
  end_time = 2
  #dtmax = 1e6
  [TimeStepper]
    type = IterationAdaptiveDT
    dt = 1e-2
    growth_factor = 1.1
  []
[]

[Outputs]
#  print_linear_residuals = false
#  sync_times = '3600 86400 2.592E6 1.5768E8'
#  perf_graph = true
  exodus = true
#  [csv]
#    type = CSV
#    sync_only = true
#  []
[]
