import subprocess
from random import randrange


params = {
    'fluid': 'TIP4P.fstprt', 
    'grow_file': 'TIP4P.grow',
    'framework': 'NU-1000.fstprt',
    'temperature': 298, 
    'min_particles': 0, 
    'max_particles': 400,  
    'min_sweeps': 20, 
    'beta_mu': -14.07,
    'trials_per': 1e5, 
    'movies_per': 1e7, 
    'hours_per_adjust': 0.2, 
    'hours_per_checkpoint': 1, 
    'seed': randrange(int(1e9)), 
    'equilibration': 1e6, 
    'procs_per_node': 48, 
    'script': __file__, 
    'min_window_size': 5,
    'restart': True
    }

params['alpha'] = 0.21912497166440476 
params['kmax'] = [8, 7, 6]

R = 1.3806488E-23 * 6.02214129E+23 # J/mol/K
params['beta'] = 1. / (params['temperature'] * R / 1e3) # mol/kJ
params['mu'] = params['beta_mu'] / params['beta']
params['hours_per_adjust'] = params['hours_per_adjust'] * params['procs_per_node']
params['hours_per_checkpoint'] = params['hours_per_checkpoint'] * params['procs_per_node']
params['mu_init'] = -7

# write TrialGrowFile for tip4p
with open(params['grow_file'], 'w') as f:
    f.write("""TrialGrowFile

particle_type 0 weight 2 transfer true site 0 num_steps 10 reference_index 0
bond true mobile_site 3 anchor_site 0 reference_index 0
branch true mobile_site 1 mobile_site2 2 anchor_site 0 anchor_site2 3 reference_index 0
""")

# write fst script
def input(params=params, simulation_file='simulation.in'):
    with open(simulation_file, "w") as myfile: myfile.write("""
# first, initialize multiple clones into windows
CollectionMatrixSplice hours_per {hours_per_adjust} ln_prob_file lnpi.txt bounds_file bounds.txt num_adjust_per_write 10
WindowExponential maximum {max_particles} minimum {min_particles} num {procs_per_node} overlap 0 alpha 1.1 min_size {min_window_size}
Checkpoint checkpoint_file checkpoint.fst num_hours {hours_per_checkpoint}

# begin description of each MC clone
RandomMT19937 seed {seed}
Configuration side_length0 39.3875 side_length1 34.11057559155958 side_length2 32.9658 xy -19.693749999999994 xz 0.0 yz 0.0 \
particle_type0 {fluid} particle_type1 {framework} group0 oxygen oxygen_site_type 0 group1 fluid \
fluid_particle_type 0 group2 framework framework_particle_type 1 add_particles_of_type1 1
Potential VisitModel Ewald alpha {alpha} kxmax {kmax[0]} kymax {kmax[1]} kzmax {kmax[2]}
Potential Model ModelTwoBodyFactory model0 LennardJones model1 ChargeScreened erfc_table_size 2e4 VisitModel VisitModelCutoffOuter energy_cutoff 1e11
RefPotential Model HardSphere group oxygen cutoff 3.1589
Potential Model ChargeScreenedIntra VisitModel VisitModelBond
Potential Model ChargeSelf
Potential VisitModel LongRangeCorrections group oxygen
ThermoParams beta {beta} chemical_potential {mu_init}
Metropolis
TrialTranslate particle_type 0 weight 0.5 tunable_param 0.2 tunable_target_acceptance 0.25
TrialParticlePivot particle_type 0 weight 0.5 particle_type 0 tunable_param 0.5 tunable_target_acceptance 0.25
Log trials_per_write {trials_per} output_file log_[sim_index].txt
Tune
CheckEnergy trials_per_update {trials_per} tolerance 1.0

# gcmc initialization and nvt equilibration
TrialAdd particle_type 0
Run until_num_particles [soft_macro_min] particle_type 0
RemoveTrial name TrialAdd
ThermoParams beta {beta} chemical_potential {mu}
Metropolis
Run num_trials {equilibration}
Movie trials_per_write {equilibration} output_file framework.xyz group framework
RemoveAnalyze name Movie
RemoveModify name Tune

# gcmc tm production
FlatHistogram Macrostate MacrostateNumParticles particle_type 0 width 1 max {max_particles} min {min_particles} \
soft_macro_max [soft_macro_max] soft_macro_min [soft_macro_min] \
Bias WLTM min_sweeps {min_sweeps} new_sweep 1 min_flatness 25 collect_flatness 20 min_collect_sweeps 20
TrialGrowFile grow_file {grow_file}
RemoveAnalyze name Log
Log trials_per_write {trials_per} output_file log_[sim_index].txt
Movie trials_per_write {movies_per} output_file fluid_[sim_index].xyz group fluid
Tune trials_per_write {trials_per} output_file tune_[sim_index].txt multistate true stop_after_iteration 1
Energy trials_per_write {trials_per} output_file energy_[sim_index].txt multistate true start_after_iteration 1
CriteriaUpdater trials_per_update 1e5
CriteriaWriter trials_per_write {trials_per} output_file crit_[sim_index].txt
""".format(**params))

# run the simulation
def run():
    restart = params.get('restart', False)
    if restart is True:
        syscode = subprocess.call('/net/people/plgrid/plgbamaz/software/feasst/build/bin/rst checkpoint.fst', shell=True, executable='/bin/bash')
    else:
        input(params)
        syscode = subprocess.call('/net/people/plgrid/plgbamaz/software/feasst/build/bin/fst < simulation.in > simulation.log', shell=True, executable='/bin/bash')

run()

