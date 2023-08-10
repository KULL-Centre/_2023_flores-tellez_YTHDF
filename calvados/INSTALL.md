# Installation Instructions

1. Make new conda environment for calvados
``` 
conda create -n calvados python=3.10
conda activate calvados
```
2. Install numba, mdtraj with conda and openmm (they have caused issues with pip install)
```
conda install numba
conda install -c conda-forge mdtraj
conda install -c conda-forge openmm cudatoolkit=11.2
```
3. Install CALVADOS and its dependencies using pip
``` 
pip install . 
```
4. Clean up faulty pip install of scipy:
```
conda install scipy
```
