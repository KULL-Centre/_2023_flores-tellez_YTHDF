import yaml

from jinja2 import Template

#### config and job functions

def write_config(fbase,config,fconfig='config.yaml'):
    with open(f'{fbase}/{fconfig}','w') as stream:
        yaml.dump(config,stream)

def write_runfile(fbase):
    stream = """import calvados as cv
from argparse import ArgumentParser

if __name__ == "__main__":
    parser = ArgumentParser()
    parser.add_argument('--config',nargs='?', default='config.yaml', const='config.yaml', type=str)
    args = parser.parse_args()

    fconfig = args.config

    cv.sim.run(fconfig=fconfig)
"""
    with open(f'{fbase}/run.py','w') as f:
        f.write(stream)

def write_job(fbase,name,temp,fconfig,batch_sys='PBS',fbash='/Users/sobuelow/.zshrc',envname='calvados'):
    """ write PBS or SLURM job """
    # PBS
    submission_PBS = Template("""#!/bin/sh
#PBS -W group_list=ku_10001 -A ku_10001
#PBS -N {{name}}_{{temp}}
### Only send mail when job is aborted or terminates abnormally
#PBS -m n
### Number of nodes
#PBS -l nodes=1:ppn=4:gpus=1
### Memory
#PBS -l mem=20gb
#PBS -l walltime=24:00:00
#PBS -o {{name}}_{{temp}}.out
#PBS -e {{name}}_{{temp}}.err

cd ${PBS_O_WORKDIR}

source {{fbash}}
conda activate {{envname}}

python run.py --config {{fconfig}}""")

    # SLURM
    submission_SLURM = Template("""#!/bin/bash
#SBATCH --job-name={{name}}_{{temp}}
#SBATCH --nodes=1
#SBATCH --partition=qgpu
#SBATCH --gres=gpu:1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=18
#SBATCH --ntasks=1
#SBATCH --mem-per-cpu=8888
#SBATCH -t 24:00:00
#SBATCH -o {{name}}_{{temp}}.out
#SBATCH -e {{name}}_{{temp}}.err

source {{fbash}}
module purge
module load cmake/3.9.4 gcc/6.5.0 openmpi/4.0.3 llvm/7.0.0 cuda/9.2.148
conda activate {{envname}}

echo $SLURM_CPUS_PER_TASK
echo $SLURM_JOB_NODELIST

python run.py --config {{fconfig}}""")

    if batch_sys == 'PBS':
        submission = submission_PBS
    elif batch_sys == 'SLURM':
        submission = submission_SLURM
    else:
        raise # space for more batch systems

    with open(f'{fbase}/job.sh', 'w') as submit:
        submit.write(submission.render(name=name,temp=f'{temp:d}',fconfig=fconfig,fbash=fbash,envname=envname))
