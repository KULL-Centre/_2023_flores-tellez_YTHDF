# DOCUMENTATION

Coarse-grained implicit-solvent simulations of biomolecules in the openMM framework.

Please cite the following references when using CALVADOS:

- G. Tesei, T. K. Schulze, R. Crehuet, K. Lindorff-Larsen. Accurate model of liquid-liquid phase behavior of intrinsically disordered proteins from optimization of single-chain properties. PNAS (2021), 118(44):e2111696118.
- G. Tesei, K. Lindorff-Larsen. Improved predictions of phase behaviour of intrinsically disordered proteins by tuning the interaction range. Open Research Europe (2022), 2(94).

# CONFIG parameters

The python script prepare.py (or prepare_minimal.py) should be used to set up the simulation. It prepares a config.yaml file alongside an auxiliary job.py file and a job submission file. Only few parameters are required (e.g. the name of the system, the box dimensions, temperature etc.) and most technical parameters can be omitted.

Description of the simulation parameters:

## job settings (required)

- `fconfig`: Name of the config yaml file.
- `ffasta`: Name of the file containing the sequence(s) in fasta format. The fasta file can contain any number of sequences, including sequences not used in the simulation. Only the sequences specified in the `proteins` subsection (or by `name` if `'proteins' : None`) are used.
- `batch_sys (options: 'PBS', 'SLURM')`: Type of submission system on cluster. Currently, only `'PBS'` or `'SLURM'` are supported and refer to specific job templates on the Danish Computerome (PBS) and ROBUST (SLURM) servers.
- `envname`: Name of conda environment to be loaded in the job submission script. This is often called `'calvados'`.
- `fbash`: Point to `.bashrc` or `.zshrc` to be loaded in the job submission script.

## general simulation settings (required)

- `box [nm]`: Simulation box dimensions.
- `name`: Name of the system. This does NOT have to be equal to the name of the protein used in the simulation. Only if no `'proteins'` dictionary is provided in the protein_settings section, then `name` is used to define the protein as well.
- `temp [K]`: Temperature.
- `ionic [M]`: Ionic strength.
- `pH`: pH value.

## runtime settings (optional)

This dictionary sets further runtime parameters. Set `runtime_settings = {}` to use default parameters.

- `'wfreq' (default: 100000)`: Writing frequency for .dcd trajectory file. The frequency is in units of 10 femtoseconds. `wfreq = 100000` thus corresponds to writing frames every 1 nanosecond.
- `'steps' (default: 100000000)`: Total number of steps written in the trajectory. The total number of recorded frames is `steps/wfreq`.
- `'runtime' (default: 0)` [hours]: Hours to run the simulation. If `runtime > 0`, the `steps` parameter is ignored.
- `'pfname' (default: 'CPU') (options: 'CPU', 'CUDA')`: Choose 'CPU' or 'CUDA' depending on your openMM configuration.
- `'restart' (default: 'checkpoint') (options: 'checkpoint', 'pdb', None)`: If and what checkpoint is read to start the simulation. At the end of the simulation, a checkpoint file `restart.chk` and a timestamped PDB file are written. The checkpoint file stores all simulation parameter and configurations, so use this for a continuation of the same system. If only the coordinates of the system should be preserved (e.g., when switching temperatures or going from slab equilibration to production run), then use the `pdb` option and provide a pdb file in `frestart` below.
- `'frestart' (default: 'restart.chk')`: Name of restart file corresponding to option provided in `restart`.

## protein settings (optional)

This dictionary sets further protein parameters. Set `protein_settings = {}` to use default parameters.

- `'proteins' (default: {name : 1})`: Dictionary of protein:number pairs to define the type and copy number of proteins in the system. Any combination of proteins can be used, e.g.: `'proteins': {'proteinA' : 10, 'proteinB' : 5}` to add 10 proteins of type proteinA and 5 proteins of type proteinB to the system.
- `'protein_topol' (default: 'single') (options: 'single', 'slab', 'grid', 'random')`: Protein placement in the simulation box. `single` places the single protein in the center (use this option only for a single protein in the system); `slab` arranges linear IDR strings in the z-center of the box, using a grid in xy direction (use only for IDRs); `grid` spaces proteins equally in the simulation box (use this e.g. for slab simulations of multidomain proteins); `random` places proteins randomly in the box, checking for clashes (harsh checks, not recommended).
- `'cutoff_lj' (default: 2.0) [nm]`: Cutoff for Ashbaugh-Hatch LJ-like potential. Refer to Tesei and Lindorff-Larsen, Open Research Europe (2022) for a discussion on this cutoff.
- `'calvados_version' (default: 2)`: CALVADOS version of original implementation (2021) or revised implementation (2022). CALVADOS 2 should be used for practically all cases.
- `'slab_eq' (default: False)`: Apply soft restraints towards the center of the box in z-direction. Useful when settings up multidomain proteins in a grid and equilibrating towards a single slab. The potential used is linear, `'k*abs(periodicdistance(x,y,z,x,y,z0))'`, where `z0` is half the box edge length.
- `'k_eq' (default: 0.01) [kJ/mol*nm]`: Force constant for slab equilibration potential (`slab_eq`).

## lipid settings (optional) (experimental!)

This dictionary sets lipid parameters. Work in progress! Set `lipid_settings = {}`.

## multidomain settings (optional)

This dictionary sets restraint parameters for folded or multidomain proteins. Set `restraint_settings = {}` to use default parameters (= no restraints). No fasta entries need to be provided for proteins in this list, as the sequence is instead read from the corresponding PDB file. The PDB file can be an atomistic structure or coarse-grained structure (if `'use_com' == True`)

- `'restraint_list' (default: [])`: List of proteins to be restrained. Protein names must be identical to the `proteins` defined in the protein settings dictionary (or to `name` if no `proteins` setting is provided).
- `'restraint_type' (default: 'harmonic') (options: 'harmonic', 'go')`: Use a harmonic potential or Go potential to keep folded domains folded. Depending on this choice, further parameters need to be provided. Harmonic potential restraints require manual domain definitions in `fdomains` below. Go potential restraints require an Alphafold 2 PAE matrix in the `fpdb` folder.
- `'k_restraint' (default: 700) [kJ/mol*nm^2]`: Force constant for restraints. This force constant has a different meaning depending on the type of restraint used (harmonic or Go). A value of `700` is ok for harmonic restraints.
- `'cutoff_restr' (default: 0.9) [nm]`: Cutoff below which to restrain bead pairs. Bonded beads (i->i+1) are never restrained, neither are beads not in the domain definitions if `'restraint_type' == 'harmonic'`.
- `'fpdb' (default: f'{cwd}/pdbs')`: Folder to store PDB files (and PAE matrices in .json format if `'restraint_type' == 'go'`). The PDB files in the folder NEED to have the same name as the protein name (+ '.pdb') (provided in `'protein'` of the protein settings, or `name` if no protein settings exist). If Go restraints are use, then the PAE files also need to have identical names (+ 'json').
- `'fdomains' (default: f'{cwd}/domains.yaml')`: Yaml file containing definitions of folded domains. An example file can be found in the `calvados/calvados/data` subfolder of the package.
- `'use_com' (default: True)`: Use center-of-mass coarse-graining. If `False`, then use coarse-graining to alpha-carbon atoms. Alpha-carbon coarse-graining can only be used if the PDB contains one CA entry for every residue.
- `'exclude_core' (default: False)`: Experimental! Exclude LJ interactions for residues in the core of folded domains. A weighted contact number (WCN) is computed to determine the core of folded domains.

## batch settings (optional)

This dictionary can be used to prepare multiple simulations at once. Several parameters from the above settings are overwritten by this dictionary if `'batch' == True`.

- `'batch' (default: False)`: Activate batch settings. If `True`, then overwrite the `'protein'` and `'protein_topol'` definitions in the protein_settings section.
- `'names' (default: [])`: List of protein names to be used for separate simulations.
- `'nproteins' (default: 1)`: Number of proteins per simulation.
- `'topology' (default: 'single')`: Protein topology per simulation.
