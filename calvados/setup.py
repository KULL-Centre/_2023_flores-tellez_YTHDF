from setuptools import setup

setup(
    name='calvados',
    version='0.2.0',    
    description='Coarse-grained implicit-solvent simulations of biomolecules',
    url='https://github.com/sobuelow/calvados',
    author='Soeren von Buelow',
    # authors = [
    #     {name='Giulio Tesei', email='giulio.tesei@bio.ku.dk'},
    #     {name='Sören von Bülow', email='soren.bulow@bio.ku.dk'},
    #     {name='Kresten Lindorff-Larsen', email='lindorff@bio.ku.dk'}
    # ]
    license='GNU GPL3',
    packages=['calvados'],
    install_requires=[
        'numpy',
        'pandas',
        'OpenMM',
        'mdtraj',
        'MDAnalysis',
        'scipy',
        'biopython',
        'Jinja2',
        'progressbar2',
        'matplotlib',
        'numba',
        'PyYAML',
        'statsmodels',
        'localcider'                            
    ],

    # include_package_data=True,
    package_data={'' : ['data/*.csv', 'data/fastabib.fasta', 'data/domains.yaml']},

    classifiers=[
        'Intended Audience :: Science/Research',
        'Operating System :: POSIX :: Linux',        
        'Programming Language :: Python :: 3>=3.7,<3.11',
    ],
)
