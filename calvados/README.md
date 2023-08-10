# CALVADOS

######## UNDER DEVELOPMENT !!!!! USE WITH CARE, NOT FOR PRODUCTION USE ##########

Coarse-grained implicit-solvent simulations of biomolecules in the openMM framework.

Please cite the following references when using the software:

- G. Tesei, T. K. Schulze, R. Crehuet, K. Lindorff-Larsen. Accurate model of liquid-liquid phase behavior of intrinsically disordered proteins from optimization of single-chain properties. PNAS (2021), 118(44):e2111696118.
- G. Tesei, K. Lindorff-Larsen. Improved predictions of phase behaviour of intrinsically disordered proteins by tuning the interaction range. Open Research Europe (2022), 2(94).

## Install

Check out the INSTALL.md documentation.

## Basic Usage (Single chain IDP simulation)

1. Copy prepare_minimal.py to your working folder and modify as needed.
2. Build config, run file, job submission file.

``` 
python prepare_minimal.py
```

This creates a config file (`config.yaml`), a run file (`run.py`) and a job submission file (`job.sh`) for use on a cluster. The files are created in a subfolder `{name}/{temp}/` with name and temp defined in the submit script.

3. The simulation must be run from the subfolder where `config.yaml` is located. It can be run either via

``` 
python job.py
```

or using the job submission script on a cluster. The config file `config.yaml` can be modified at any time before the start of the simulation.

4. The package has postprocessing and analysis scripts. For example:

```
import calvados as cal
import MDAnalysis as mda

u = mda.Universe('top.pdb','{your_file}.dcd')
uref = mda.Universe('top.pdb')

# RMSD to reference, RMSD to mean, RMSF
rref, rmean, rmsf = cal.analysis.calc_rmsd(u,uref,select={selection_string})

# scaling exponent
ij, dij, r0, v = cal.analysis.fit_scaling_exp(u,u.atoms)
```

or calculate centered concentration profile from slab simulations:

```
import calvados as cal
import MDAnalysis as mda

u = mda.Universe('top.pdb','{your_file}.dcd')

path = '.'
name = 'myprotein'
temp = 293. # K

hs, z = cal.postprocess.center_slab(path,name,temp)
# where hs = (n_frames,n_bins) and z = the centers of equally spaced 1 Angstrom bins
```

## Further Info

- prepare_minimal.py is the minimal starting point to build a config file, run file and job submission file
- prepare.py is a full version of the submission script, with examples for all simulation parameters
- EVERYTHING ELSE is part of the internal package and should only be modified when developing the package
