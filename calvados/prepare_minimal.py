import subprocess
import os

from calvados.cfg import write_config, write_runfile, write_job

cwd = os.getcwd()

#########################################

# job settings
fconfig = 'config.yaml'
ffasta = f'{cwd}/ect_idrs.fasta'
batch_sys = 'SLURM'
envname = 'calvados' # conda environment
fbash = '/home/sobuelow/.bashrc'

# general simulation settings
name = 'ECT1N'
box = [30., 30., 30.]
temp = 293
ionic = 0.15
pH = 7.0

runtime_settings = {
    'wfreq' : 100, # dcd writing frequency
    'steps' : 10000, # number of simulation steps
    'pfname' : 'CPU',
}

protein_settings = {}
lipid_settings = {}
restraint_settings = {}

#########################################

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
