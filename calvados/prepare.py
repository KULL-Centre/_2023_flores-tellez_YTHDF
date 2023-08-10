import subprocess
import time
import os

from calvados.cfg import write_config, write_runfile, write_job

from Bio import SeqIO

cwd = os.getcwd()

#########################################

fconfig = 'config.yaml'
ffasta = f'{cwd}/ect_idrs.fasta' # no fasta needed if pdb provided
batch_sys = 'PBS'
envname = 'calvados' # conda environment
fbash = '/home/people/sorbul/.bashrc'

# general simulation settings
box = [30., 30., 30.]
name = 'ECT1N'
temp = 293
ionic = 0.15
pH = 7.0

runtime_settings = {
    'wfreq' : 100000, # dcd writing frequency
    'steps' : 100000000, # number of simulation steps
    'runtime' : 0, # hours (overwrites steps if runtime is > 0)   
    'pfname' : 'CPU' # 'CUDA',
    # Start system from restart.chk (if it exists), from pdb, or from scratch:
    'restart' : 'checkpoint', # 'pdb', None 
    'frestart' : 'restart.chk', # restart file
}

# protein parameters (single simulation)
protein_settings = {
    'proteins' : {
        'ECT1N': 1,
        },
    'protein_topol' : 'single' # 'grid', 'random', 'slab'
    'cutoff_lj' : 2.0, # set the cutoff for the nonionic interactions
    'eps_factor' : 0.2, # eps_lj = eps_factor * 4.184
    'calvados_version' : 2,
    'slab_eq' : False, # slab equilibration flag (restrain towards center)
    'k_eq' : 0.01 # kJ/mol*nm for linear restraints towards box center in z
}

# lipid parameters
lipid_settings = {
    'lipids' : {
#        'DOPC': 100
        },
    'lipid_topol' : 'random', # 'random', 'sheet'
    'eps_lipids' : 0.91,
    'cutoff_wca' : 1.5,
    'cutoff_wcos' : 2.5
}

# multidomain settings
restraint_settings = {
    'restraint_list' : [],
    'restraint_type' : 'harmonic', # 'go'
    'k_restraint' : 700.,
    'cutoff_restr' : 0.9, # distance cutoff for restraints
    'fpdb' : f'{cwd}/pdbs', # pdb folder, must also contain pae if go restraints are used
    'fdomains' : f'{cwd}/domains.yaml', # domain defs for harmonic restraints
    'use_com' : True, # center of mass coarse-graining
    'exclude_core' : False
}

# batch runs of single chain simulations of different proteins
# This overwrites some protein parameters above if 'batch' == True
batch_settings = {
    'batch' : False,
    'names' : [], # overwrite protein_settings['proteins']
    'nproteins' : 1, # overwrite protein_settings['proteins']
    'topology' : 'single' # overwrite protein_settings['protein_topol']
}

#########################################

if batch_settings['batch']: # multiple single chain simulations
    if len(batch_settings['names']) == 0:
        records = SeqIO.to_dict(SeqIO.parse(ffasta, "fasta"))
        names = []
        for key in records.keys(): # read entire provided fasta
            names.append(key)
else: # single simulation
    names = [name]

for name in names:
    if batch_settings['batch']: # overwrite with batch settings
        protein_settings['proteins'] = {name : batch_settings['nproteins']}
        protein_settings['protein_topol'] = batch_settings['topology']

    fbase = f'{cwd}/{name}/{temp:d}'
    subprocess.run(f'mkdir -p {fbase}',shell=True)
    config = dict(
        fbase=fbase,
        name=name,temp=temp,ionic=ionic,
        box=box,pH=pH,
        ffasta=ffasta,
        runtime_settings=runtime_settings,
        protein_settings=protein_settings,
        lipid_settings=lipid_settings,
        restraint_settings=restraint_settings,
        )
    write_config(fbase,config,fconfig=fconfig)
    write_runfile(fbase)
    write_job(fbase,name,temp,fconfig,batch_sys=batch_sys,envname=envname,fbash=fbash)
