import numpy as np
import openmm
from openmm import app, unit
import pandas as pd
from datetime import datetime

import mdtraj as md
import MDAnalysis as mda
from MDAnalysis.analysis import distances

from scipy import constants

import Bio.SeqUtils
import progressbar

import json
import importlib

from jinja2 import Template
import os

import calvados.build as build
import calvados.interactions as interactions
from .sequence import *
from .analysis import self_distances, calc_wcn

class Sim:
    def __init__(self,config):
        """ 
        simulate openMM Calvados;
        parameters are provided by config dictionary """
        # parse config
        # self.fbase = config['fbase'] # where to run sim
        self.fbase = '.'
        self.pkg_base = importlib.resources.files('calvados')

        self.name = config['name']
        self.temp = config['temp']
        self.ionic = config['ionic']
        self.box = np.array(config['box'])
        
        self.pH = config.get('pH', 7.0)
        self.ffasta = config.get('ffasta', 'fastabib.fasta')

        self.verbose = config.get('verbose',False)

        # runtime settings
        crun = config.get('runtime_settings', {})
        self.steps = crun.get('steps', 100000000)
        self.wfreq = crun.get('wfreq', 100000)
        self.pfname = crun.get('pfname', 'CPU')
        self.runtime = crun.get('runtime', 0)
        self.restart = crun.get('restart', 'checkpoint')
        self.frestart = crun.get('frestart', 'restart.chk')
        # self.temp_increment = crun.get('temp_increment', 0)
        # self.temp_freq = crun.get('temp_freq', 0)

        # protein (and RNA) settings
        cp = config.get('protein_settings', {})
        self.proteins = cp.get('proteins', {self.name : 1})
        self.protein_topol = cp.get('protein_topol', 'single')
        self.cutoff_lj = cp.get('cutoff_lj', 2.0)
        self.calvados_version = cp.get('calvados_version', 2)
        self.eps_lj = cp.get('eps_factor', 0.2) * 4.184
        self.slab_eq = cp.get('slab_eq', False)
        self.k_eq = cp.get('k_eq', 0.)
        self.rna = cp.get('rna',{})
        self.bondlength = cp.get('bondlength','CA')
        self.kb = cp.get('kb',8033.0)

        self.nproteins = 0
        for key,val in self.proteins.items():
            self.nproteins += val
        print(f'Total number of proteins in the system: {self.nproteins}')

        self.nrna = 0
        for key,val in self.rna.items():
            self.nrna += val
        print(f'Total number of RNA in the system: {self.nrna}')

        # lipid settings
        cl = config.get('lipid_settings', {})
        self.lipids = cl.get('lipids', {})
        self.lipid_topol = cl.get('lipid_topol', 'sheet')
        self.eps_lipids = cl.get('eps_lipids', 0.91)
        self.cutoff_wca = cl.get('cutoff_wca', 1.5)
        self.cutoff_wcos = cl.get('cutoff_wcos', 2.5)

        # restraint settings
        crstr = config.get('restraint_settings', {})
        self.restraint_list = crstr.get('restraint_list', [])

        # init protein interactions
        self.load_parameters()
        self.eps_yu, self.k_yu = interactions.genParamsDH(self.temp,self.ionic)
        self.hb, self.ah, self.yu = interactions.init_protein_interactions(
            self.eps_lj,self.cutoff_lj,self.k_yu,calvados_version=self.calvados_version
            )

        # init lipid interactions
        if len(self.lipids.keys()) > 0:
            self.wbond, self.wbend, self.wca, self.wcos = interactions.init_lipid_interactions(
                self.temp,self.eps_lipids,self.residues,self.cutoff_wca,self.cutoff_wcos
                )

        if len(self.restraint_list) > 0:
            self.fpdb = crstr.get('fpdb', 'pdbs')
            self.use_com = crstr.get('use_com', True)
            self.restraint_type = crstr.get('restraint_type', 'harmonic')
            self.k_restraint = crstr.get('k_restraint', 700.)
            self.cutoff_restr = crstr.get('cutoff_restr', 0.9)
            self.exclude_core = crstr.get('exclude_core', False)
            self.sig_shift = crstr.get('sig_shift',0.8)
            self.colabfold = crstr.get('colabfold',0)
            self.exclude_lj = crstr.get('exclude_lj',False)
            self.cs = interactions.init_restraints(self.restraint_type)
            
            if self.restraint_type == 'harmonic':
                self.fdomains = crstr.get('fdomains', 'domains.yaml')
            elif self.restraint_type == 'go':
                self.scLJ = interactions.init_scaled_LJ(self.eps_lj,self.cutoff_lj)
                self.scYU = interactions.init_scaled_YU(self.k_yu)
        if self.slab_eq:
            self.rcent = interactions.init_eq_restraints(self.box,self.k_eq)

    def load_parameters(self):
        """ Load lambda, MW, charges for bead types """

        if self.calvados_version in [1,2]:
            nm = 'residues'
        elif self.calvados_version == 3:
            nm = 'residues_RHYFW'
        elif self.calvados_version == 4:
            nm = 'residues_RHYFW_075'
        
        f = f'{self.pkg_base}/data/{nm}.csv'
        self.residues = pd.read_csv(f).set_index('one')

        if self.calvados_version in [1,2]: 
            if self.calvados_version == 1:
                self.residues.lambdas = self.residues['CALVADOS1'] # select CALVADOS1 or CALVADOS2 stickiness parameters
            else:
                self.residues.lambdas = self.residues['CALVADOS2'] # select CALVADOS1 or CALVADOS2 stckiness parameters

    def build_system(self):
        """ 
        Set up system
        * component definitions
        * build particle coordinates
        * define interactions
        * set restraints
        """

        self.top = md.Topology()
        self.system = openmm.System()
        a, b, c = build.build_box(self.box[0],self.box[1],self.box[2])
        self.system.setDefaultPeriodicBoxVectors(a, b, c)

        self.nparticles = 0 # bead counter
        self.ctprot = 0 # protein molecule counter (for xy grid in slab sim)
        self.pos = []
        self.make_sys_dataframe() # pool component information into one df
        print(self.df_sys)
        if self.protein_topol == 'slab':
            self.xyzgrid = build.build_xygrid(self.nproteins+self.nrna,self.box)
        if self.protein_topol == 'grid':
            self.xyzgrid = build.build_xyzgrid(self.nproteins+self.nrna,self.box)
        for compname, comp in self.df_sys.iterrows():
            for idx in range(comp.nmol):
                # particle definitions
                self.add_mdtraj_topol(comp.fasta)
                self.add_particles_system(comp.fasta,protein=comp.protein)
                # add beads
                if comp.protein:
                    xs = self.place_protein(comp,d=0.38)
                elif comp.rna:
                    xs = self.place_protein(comp,d=0.5)
                elif comp.lipid:
                    xs = self.place_lipid(comp)
                else:
                    raise
                # protein interactions
                qs, qs_abs = get_qs(comp.fasta,flexhis=True,pH=self.pH,residues=self.residues)
                if comp.protein:
                    qs = patch_terminal_qs(qs) # N-term add 1; C-term subtract 1
                    self.add_protein_interactions(comp,xs,qs)
                elif comp.rna:
                    self.add_protein_interactions(comp,xs,qs)
                elif comp.lipid:
                    self.add_protein_interactions(comp,xs,qs)
                # lipid interactions
                if len(self.lipids.keys()) > 0:
                    self.add_lipid_interactions(comp)
                # core restraints
                if comp.name in self.restraint_list:
                    if self.exclude_core:
                        if self.restraint_type == 'harmonic':
                            fdomains = self.fdomains
                        else:
                            fdomains = None
                        input_pdb = f'{self.fpdb}/{comp.name}.pdb'
                        ncore = self.set_core_exclusions(comp,input_pdb,fdomains=fdomains)
                        print(f'Number of core exclusions: {ncore}')
                if self.slab_eq:
                    self.add_eq_restraints(comp)

        self.pdb_cg = f'{self.fbase}/top.pdb'
        a = md.Trajectory(self.pos, self.top, 0, self.box, [90,90,90])
        if self.restart != 'pdb': # only save new topology if no system pdb is given
            a.save_pdb(self.pdb_cg)
        # add forces to system
        for force in [self.hb, self.yu, self.ah]:
            self.system.addForce(force)
        if len(self.lipids.keys()) > 0:
            for force in [self.wbond, self.wbend, self.wca, self.wcos]:
                self.system.addForce(force)
        if len(self.restraint_list) > 0:
            self.system.addForce(self.cs)
            print(f'Number of restraints: {self.cs.getNumBonds()}')
            if (self.restraint_type == 'go') and (self.exclude_lj):
                self.system.addForce(self.scLJ)
                self.system.addForce(self.scYU)
                print(f'Number of scaled LJ pairs: {self.scLJ.getNumBonds()}')
                print(f'Number of scaled YU pairs: {self.scYU.getNumBonds()}')
        if self.slab_eq:
            self.system.addForce(self.rcent)

        print(f'{self.nparticles} particles in the system')
        print('---------- FORCES ----------')
        print(f'ah: {self.ah.getNumParticles()} particles, {self.ah.getNumExclusions()} exclusions')
        print(f'yu: {self.yu.getNumParticles()} particles, {self.yu.getNumExclusions()} exclusions')
        if len(self.lipids.keys()) > 0:
            print(f'wca: {self.wca.getNumParticles()} particles, {self.wca.getNumExclusions()} exclusions')
            print(f'wcos: {self.wcos.getNumParticles()} particles, {self.wcos.getNumExclusions()} exclusions')
        if self.slab_eq:
            print(f'Equilibration restraints (rcent) towards box center in z direction')
            print(f'rcent: {self.rcent.getNumParticles()} restraints')

    def set_core_exclusions(self,comp,input_pdb,fdomains=None,cutoff=25.):
        seq = comp.fasta
        N = len(seq)
        offset = self.nparticles - N # to get indices of current comp in context of system
        xprot = build.geometry_from_pdb(input_pdb,use_com=self.use_com) # read from pdb
        wcn = calc_wcn(comp,xprot,fdomains=fdomains)
        np.savetxt('wcn.txt', wcn)
        wcn_binary = np.where(wcn>cutoff,1,0)
        np.savetxt('wcn_binary.txt', wcn_binary)

        ncore = 0

        for idx, w in enumerate(wcn_binary):
            s = seq[idx]
            if w == 1:
                self.ah.setParticleParameters(idx+offset,
                [self.residues.loc[s].sigmas*unit.nanometer, 0.*unit.dimensionless])
                ncore += 1
        return ncore

    def place_protein(self,comp,ntries=10000,d=0.38):
        """
        place proteins based on topology:
        * random
        * slab
        """
        if comp.name in self.restraint_list: # contains folded domains
            input_pdb = f'{self.fpdb}/{comp.name}.pdb'
            xprot = build.geometry_from_pdb(input_pdb,use_com=self.use_com) # read from pdb
        else:
            if self.protein_topol in ['single','random','grid']:
                xprot = build.build_spiral(len(comp.fasta),arc=d) # Archimedes spiral
            elif self.protein_topol == 'slab':
                xprot = build.build_linear(len(comp.fasta),d=d) # linear chain in z direction
        if self.protein_topol == 'slab':
            x0 = self.xyzgrid[self.ctprot]
            self.ctprot += 1 # protein molecule counter
            x0[2] = self.box[2] / 2. # center in z
            xs = x0 + xprot
        elif self.protein_topol == 'grid':
            x0 = self.xyzgrid[self.ctprot]
            self.ctprot += 1 # protein molecule counter
            xs = x0 + xprot
        elif self.protein_topol == 'single':
            x0 = self.box * 0.5 # place in center of box
            xs = x0 + xprot
        else:
            ntry = 0
            while True: # random placement
                ntry += 1
                if ntry > ntries:
                    raise
                if self.protein_topol == 'random':
                    x0 = build.draw_starting_vec(self.box) # random point in box
                else:
                    raise
                xs = x0 + xprot # shift geometry to random point
                walls = build.check_walls(xs,self.box) # check if outside box
                if walls:
                    continue
                clash = build.check_clash(xs,self.pos,self.box) # check if clashes with existing pos
                if clash:
                    continue
                else:
                    break
        for x in xs:
            self.pos.append(x)
            self.nparticles += 1
        return xs # positions of the comp (to be used for restraints)

    def place_lipid(self,comp,d=1.3):
        """ Add lipids to pos list """
        while True:
            if self.lipid_topol == 'random':
                x0 = build.draw_starting_vec(self.box)
                vec = build.draw_vec(d)
                N = len(comp.fasta)
                xs = np.array([x0+n*vec for n in range(N)])
            walls = build.check_walls(xs,self.box)
            if walls:
                continue
            clash = build.check_clash(xs,self.pos,self.box)
            if clash:
                continue
            else:
                break
        for x in xs:
            self.pos.append(x)
            self.nparticles += 1
        return xs

    def add_lipid_interactions(self,comp):
        """ 
        Tesei, Vazdar, Lund, JCP 2018, 149
        Lipid interactions for one molecule of composition comp"""
        # set lipid interactions for a single molecules

        N = len(comp.fasta)
    
        begin = self.nparticles - N
        end = self.nparticles

        if comp.lipid:
            n = 1
        else:
            n = 0
        # nonbonded lipid interactions
        for a in comp.fasta:
            if a == 'Z': # tail
                m = 1
            else:
                m = 0
            self.wca.addParticle([self.residues.loc[a].sigmas*unit.nanometer,n])
            self.wcos.addParticle([self.residues.loc[a].sigmas*unit.nanometer,m,n])

        # bonded lipid interactions
        if comp.lipid:
            for i in range(begin,end-1):
                self.wbond.addBond(i,i+1, [0.5*self.residues.loc['Z'].sigmas*unit.nanometer])
                # print(f'adding bond between {i} and {i+1}')
                self.wca.addExclusion(i, i+1)
                self.wcos.addExclusion(i, i+1)
                self.ah.addExclusion(i, i+1)
                self.yu.addExclusion(i, i+1)

            self.wbend.addBond(begin,end-1, [0.5*self.residues.loc['Z'].sigmas*unit.nanometer])
            self.wca.addExclusion(begin,end-1)
            self.wcos.addExclusion(begin,end-1)
            # add exclusion for other forces
            self.ah.addExclusion(begin,end-1)
            self.yu.addExclusion(begin,end-1)

    def add_protein_interactions(self,comp,xs,qs):
        """ 
        Protein interactions for one molecule of composition comp
        """

        N = len(comp.fasta)
        offset = self.nparticles - N # to get indices of current comp in context of system

        # bond distance
        if comp.protein:
            d0 = 0.38
        elif comp.rna:
            d0 = 0.5

        # kscale = 1.

        # prep restraints
        if comp.name in self.restraint_list:
            restr = True
            restr_pairlist = []
            dmap = self_distances(xs) # in nm
            print(f'setting {self.restraint_type} restraints for {comp.name}')
            input_pdb = f'{self.fpdb}/{comp.name}.pdb'
            bfac = build.bfac_from_pdb(input_pdb,confidence=0.)
            if self.restraint_type == 'harmonic':
                ssdomains = build.get_ssdomains(comp.name,self.fdomains)
            elif self.restraint_type == 'go':
                input_pae = f'{self.fpdb}/{comp.name}.json'
                pae_inv = build.load_pae_inv(input_pae,colabfold=self.colabfold)
        else:
            restr = False

        # LJ and electrostatics
        for a,q in zip(comp.fasta,qs):
            self.yu.addParticle([q*self.eps_yu*unit.nanometer*unit.kilojoules_per_mole])
            if self.calvados_version in [3,4]:
                m = 1.0 if a in ['R','H','F','Y','W'] else 0.0
                self.ah.addParticle([self.residues.loc[a].sigmas*unit.nanometer, self.residues.loc[a].lambdas*unit.dimensionless, m*unit.dimensionless])
            else:
                # kscale = 0.99 * (1 - bfac[i]) + 0.01 # scale bonded term
                self.ah.addParticle([self.residues.loc[a].sigmas*unit.nanometer, self.residues.loc[a].lambdas*unit.dimensionless])

        bond_pairlist = []
        if comp.protein or comp.rna:
            # cutoff = 25.
            # sigmoid_shift = 6.
            # a = 0.5
            # if restr and (self.restraint_type == 'go'):
                # wcn = calc_wcn(comp,xs,ssonly=False)
                # print(wcn)
                
                # kscale = 0.9 * (1 - np.exp(a*(wcn-sigmoid_shift)) / (np.exp(a*(wcn-sigmoid_shift)) + 1)) + 0.1 # inverted shifted sigmoid
                # print(kscale)
            for i in range(0,N-1):
                ai = comp.fasta[i] # amino acid type of i
                si = self.residues.loc[ai].sigmas
                li = self.residues.loc[ai].lambdas
                qi = qs[i]
                for j in range(i+1,N):
                    aj = comp.fasta[j] # amino acid type of j
                    sj = self.residues.loc[aj].sigmas
                    lj = self.residues.loc[aj].lambdas
                    qj = qs[j]
                    ss = False
                    # restraints
                    width = 50 # shift width
                    if restr:
                        bfac_ij = (bfac[i] + bfac[j])/2 # mean pLDDT of i an j
                        sigmoid = np.exp(width*(bfac_ij-self.sig_shift)) / (np.exp(width*(bfac_ij-self.sig_shift)) + 1) # 0.7 --> 0, 0.8 --> 0.5, 0.9 --> 1.0
                        if dmap[i,j] < self.cutoff_restr: # nm
                            # print(i, j, bfac[i], bfac[j], bfac_ij, sigmoid)
                            if self.restraint_type == 'harmonic':
                                ss = check_ssdomain(ssdomains,i,j)
                                if ss: # both residues in structured domains
                                    self.cs.addBond(i+offset,j+offset, dmap[i,j]*unit.nanometer,
                                            self.k_restraint*unit.kilojoules_per_mole/(unit.nanometer**2))
                                    self.yu.addExclusion(i+offset, j+offset)
                                    self.ah.addExclusion(i+offset, j+offset)
                                    restr_pairlist.append([i+offset+1,j+offset+1, self.k_restraint, dmap[i,j]]) # 1-based
                                    restr = True
                            elif self.restraint_type == 'go':
                                if j > i+1: 
                                    if pae_inv[i,j]**2 * sigmoid > 0.01:
                                        k = self.k_restraint * pae_inv[i,j]**2 * sigmoid # scale by pae and pLDDT
                                        self.cs.addBond(i+offset,j+offset, [dmap[i,j]*unit.nanometer, k*unit.kilojoules_per_mole])
                                        if self.exclude_lj:
                                            self.yu.addExclusion(i+offset, j+offset)
                                            self.ah.addExclusion(i+offset, j+offset)
                                            s = 0.5*(si+sj)
                                            l = 0.5*(li+lj)
                                            n = 1. - pae_inv[i,j]**2 * sigmoid
                                            self.scLJ.addBond(i+offset,j+offset, [s*unit.nanometer, l*unit.dimensionless, n*unit.dimensionless])
                                            qij = qi * qj * self.eps_yu * self.eps_yu *unit.nanometer*unit.kilojoules_per_mole *unit.nanometer * unit.kilojoules_per_mole
                                            self.scYU.addBond(i+offset,j+offset, [qij])
                                        restr_pairlist.append([i+offset+1,j+offset+1, k, dmap[i,j]]) # 1-based
                                        restr = True
                    # bonds
                    if (j == i+1) and not ss: # deal with bonds, don't add bonds in ss for 'harmonic'
                        if restr:
                            if self.bondlength == 'CA':
                                d = d0
                            elif self.bondlength == 'COM':
                                dmap = self_distances(xs) # in nm
                                d = dmap[i,j]
                            elif (self.bondlength == 'flex') and (self.restraint_type == 'go'):
                                d = sigmoid * dmap[i,j] + (1. - sigmoid) * d0
                            else:
                                raise
                        else:
                            d = d0 # currently force to CA without restraints
                        bidx = self.hb.addBond(i+offset, j+offset, d*unit.nanometer, self.kb*unit.kilojoules_per_mole/(unit.nanometer**2))
                        bond_pairlist.append([i+offset+1,j+offset+1,bidx,d0,self.kb]) # 1-based
                        self.yu.addExclusion(i+offset, j+offset)
                        self.ah.addExclusion(i+offset, j+offset)
                        
                        if len(self.lipids.keys()) > 0:
                            self.wca.addExclusion(i+offset, j+offset)
                            self.wcos.addExclusion(i+offset, j+offset)
        # write lists
        if self.verbose:
            self.write_bonds(comp,bond_pairlist)
            if restr:
                self.write_restraints(comp,restr_pairlist)  

    def write_bonds(self,comp,bond_pairlist):
        # write bonds
        with open(f'{self.fbase}/bonds_{comp.name}_idx.txt','w') as f:
            f.write('i\tj\tb_idx\td[nm]\tk[kJ/mol/nm^2]\n')
            for b in bond_pairlist:
                f.write(f'{int(b[0])}\t{int(b[1])}\t{int(b[2])}\t{b[3]:.4f}\t{b[4]:.4f}\n')

    def write_restraints(self,comp,restr_pairlist):
        with open(f'{self.fbase}/restr_{comp.name}_idx.txt','w') as f:
            f.write('i j fc d[nm]\n')
            for r in restr_pairlist:
                f.write(f'{int(r[0])} {int(r[1])} {r[2]:.4f} {r[3]:.4f}\n')

    def add_eq_restraints(self,comp):
        N = len(comp.fasta)
        offset = self.nparticles - N # to get indices of current comp in context of system
        for i in range(0,N):
            self.rcent.addParticle(i+offset)

    def add_mdtraj_topol(self,fasta):
        """ add one molecule to mdtraj topology """
        chain = self.top.add_chain()
        for idx,resname in enumerate(fasta):
            res = self.top.add_residue(resname, chain, resSeq=idx+1)
            self.top.add_atom(resname, element=md.element.carbon, residue=res)
        for i in range(chain.n_atoms-1):
            self.top.add_bond(chain.atom(i),chain.atom(i+1))

    def add_particles_system(self,fasta,protein=True):
        """ add particles of one molecule to openMM system """
        if protein: # if component is protein
            self.system.addParticle((self.residues.loc[fasta[0]].MW+2)*unit.amu)
            for a in fasta[1:-1]:
                self.system.addParticle(self.residues.loc[a].MW*unit.amu)
            self.system.addParticle((self.residues.loc[fasta[-1]].MW+16)*unit.amu)
        else:
            for a in fasta:
                self.system.addParticle(self.residues.loc[a].MW*unit.amu)

    def make_sys_dataframe(self):
        self.df_sys = build.initDf() # initialize dataframe
        for p, nmol in self.proteins.items():
            self.addComponent(p,nmol,protein=True,rna=False,lipid=False)
        for p, nmol in self.rna.items():
            self.addComponent(p,nmol,protein=False,rna=True,lipid=False)
        for l, nmol in self.lipids.items():
            self.addComponent(l,nmol,protein=False,rna=False,lipid=True)

    def addComponent(self,cname,nmol,protein=True,rna=False,lipid=False):
        if cname in self.restraint_list:
            input_pdb = f'{self.fpdb}/{cname}.pdb'
            fasta = seq_from_pdb(input_pdb)
        else:
            records = read_fasta(self.ffasta)
            fasta = records[cname].seq
        self.df_sys.loc[cname] = dict(ionic=self.ionic,temp=self.temp,fasta=list(fasta),nmol=nmol,protein=protein,rna=rna,lipid=lipid)

    def simulate(self):
        if self.restart == 'pdb':
            pdb = app.pdbfile.PDBFile(self.frestart)
        else:
            pdb = app.pdbfile.PDBFile(self.pdb_cg)
        # use langevin integrator
        integrator = openmm.openmm.LangevinIntegrator(self.temp*unit.kelvin,0.01/unit.picosecond,0.01*unit.picosecond)
        print(integrator.getFriction(),integrator.getTemperature())

        # assemble simulation
        platform = openmm.Platform.getPlatformByName(self.pfname)
        simulation = app.simulation.Simulation(pdb.topology, self.system, integrator, platform=platform)

        fcheck_in = f'{self.fbase}/{self.frestart}'
        fcheck_out = f'{self.fbase}/restart.chk'

        if (os.path.isfile(fcheck_in)) and (self.restart == 'checkpoint'):
            print(f'Reading check point file {fcheck_in}')
            print(f'Appending trajectory to {self.fbase}/{self.name:s}.dcd')
            simulation.loadCheckpoint(fcheck_in)
            simulation.reporters.append(app.dcdreporter.DCDReporter(f'{self.fbase}/{self.name:s}.dcd',self.wfreq,append=True))
        else:
            if self.restart == 'pdb':
                print(f'Reading in system configuration {self.frestart}')
            elif self.restart == 'checkpoint':
                print(f'No checkpoint file {self.frestart} found: Starting from new system configuration')
            elif self.restart == None:
                print('Starting from new system configuration')
            else:
                raise

            if os.path.isfile(f'{self.fbase}/{self.name:s}.dcd'): # backup old dcd if not restarting from checkpoint
                now = datetime.now()
                dt_string = now.strftime("%Y%d%m_%Hh%Mm%Ss")
                print(f'Backing up existing {self.fbase}/{self.name:s}.dcd to {self.fbase}/backup_{self.name:s}_{dt_string}.dcd')
                os.system(f'mv {self.fbase}/{self.name:s}.dcd {self.fbase}/backup_{self.name:s}_{dt_string}.dcd')
            print(f'Writing trajectory to new file {self.fbase}/{self.name:s}.dcd')
            # print(pdb.positions)
            simulation.context.setPositions(pdb.positions)
            simulation.minimizeEnergy()
            state = simulation.context.getState(getPositions=True)
            pos2 = state.getPositions(asNumpy=True)
            simulation.reporters.append(app.dcdreporter.DCDReporter(f'{self.fbase}/{self.name:s}.dcd',self.wfreq))

        simulation.reporters.append(app.statedatareporter.StateDataReporter(f'{self.fbase}/{self.name}_{self.temp:d}.log',int(self.wfreq*10),
                step=True,speed=True,elapsedTime=True,separator='\t'))

        # run simulation
        print("STARTING SIMULATION")
        if self.runtime > 0: # in hours
            simulation.runForClockTime(self.runtime*unit.hour, checkpointFile=fcheck_out, checkpointInterval=30*unit.minute)
        else:
            # batch = self.bond_update_freq
            # nbatches = int(self.steps / batch)
            # temp_batch = 0
            # temps = [self.temp]
            # kscales = []
            nbatches = 10
            batch = int(self.steps / nbatches)
            for i in progressbar.progressbar(range(nbatches),min_poll_interval=1):
                simulation.step(batch)
                simulation.saveCheckpoint(fcheck_out)
                # state = simulation.context.getState(getPositions=True)
                # pos_tmp = np.array(state.getPositions(asNumpy=True))
                # update bonds
                # kscale = self.update_bonds(pos_tmp,i)
                # kscales.append(kscale)
                # self.hb.updateParametersInContext(simulation.context)

                # temp_batch += batch
                # if temp_batch >= self.temp_freq:
                # update temperature
                    # newT = self.temp + self.temp_increment * i
                    # temps.append(newT)
                    # self.update_temp(newT)
                    # self.yu.updateParametersInContext(simulation.context)
                    # simulation.context.reinitialize(preserveState=True)
                    # temp_batch = 0
            # temps = np.array(temps)
            # np.savetxt('temps.txt',temps.T)
            # kscales = np.array(kscales)
            # np.save(f'kscales.npy',kscales)
        simulation.saveCheckpoint(fcheck_out)
        
        now = datetime.now()
        dt_string = now.strftime("%Y%d%m_%Hh%Mm%Ss")

        state_final = simulation.context.getState(getPositions=True)
        rep = app.pdbreporter.PDBReporter(f'{self.name}_{dt_string}.pdb',0)
        rep.report(simulation,state_final)

    # def update_bonds(self,pos,step,d=0.38,sigmoid_shift=6.,a=0.5):
    #     offset = 0
    #     ct_bond = 0
    #     bond_pairlist = []
    #     for compname, comp in self.df_sys.iterrows():
    #         N = len(comp.fasta)
    #         for idx in range(comp.nmol):
    #             # print(offset,N)
    #             pos_tmp = pos[offset:offset+N]
    #             # print(pos_tmp.shape)
    #             wcn = calc_wcn(comp,pos_tmp,ssonly=False)
    #             kscale = 0.9 * (1 - np.exp(a*(wcn-sigmoid_shift)) / (np.exp(a*(wcn-sigmoid_shift)) + 1)) + 0.1 # inverted shifted sigmoid
    #             for i in range(0,N-1):
    #                 j = i+1
    #                 bond_pairlist.append([i+offset+1,j+offset+1,ct_bond,d,8033.0*kscale[i]*kscale[j], wcn[i], wcn[j]]) # 1-based
    #                 self.hb.setBondParameters(ct_bond, i+offset, j+offset, d*unit.nanometer, 8033.0*kscale[i]*kscale[j]*unit.kilojoules_per_mole/(unit.nanometer**2))
    #                 ct_bond += 1
    #             offset += N
    #     if self.verbose:
    #         with open(f'bonds_step_{step}.txt', 'w') as f:
    #             f.write('i\tj\tb_idx\td[nm]\tk\twcn[i]\t[wcn[j]\n')
    #             for b in bond_pairlist:
    #                 f.write(f'{int(b[0])}\t{int(b[1])}\t{int(b[2])}\t{b[3]:.4f}\t{b[4]:.4f}\t{b[5]:.4f}\t{b[6]:.4f}\n')
    #     return kscale

    # def update_temp(self,newT):
    #     self.integrator = openmm.openmm.LangevinIntegrator(newT*unit.kelvin,0.01/unit.picosecond,0.01*unit.picosecond)
    #     eps_yu_tmp, k_yu_tmp = interactions.genParamsDH(newT,self.ionic)
    #     self.yu.setGlobalParameterDefaultValue(0,k_yu_tmp/unit.nanometer)
    #     self.yu.setGlobalParameterDefaultValue(1,np.exp(-k_yu_tmp*4.0)/4.0/unit.nanometer)

    #     offset = 0
    #     yu_idx = 0
    #     for compname, comp in self.df_sys.iterrows():
    #         N = len(comp.fasta)
    #         qs, qs_abs = get_qs(comp.fasta,flexhis=True,pH=self.pH,residues=self.residues)
    #         if comp.protein:
    #             qs = patch_terminal_qs(qs) # N-term add 1; C-term subtract 1
    #         for idx in range(comp.nmol):
    #             for i, q in enumerate(qs):
    #                 self.yu.setParticleParameters(i+offset,[q*eps_yu_tmp*unit.nanometer*unit.kilojoules_per_mole])

def run(fconfig='config.yaml'):
    with open(f'{fconfig}','r') as stream:
        config = yaml.safe_load(stream)

    mysim = Sim(config)
    mysim.build_system()
    mysim.simulate()
    return mysim
