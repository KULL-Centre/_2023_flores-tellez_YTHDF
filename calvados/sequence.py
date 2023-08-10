import numpy as np
import Bio.SeqUtils
from scipy import optimize

import MDAnalysis as mda
from MDAnalysis.analysis import distances

import yaml
import random

import re

import MDAnalysis.analysis.rms as rms
import MDAnalysis.analysis.distances as distances
from MDAnalysis.analysis.align import AlignTraj
from MDAnalysis.analysis import align

from Bio import SeqIO, SeqUtils
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from localcider.sequenceParameters import SequenceParameters

import progressbar

### SEQUENCE INPUT / OUTPUT
def read_fasta(ffasta):
    records = SeqIO.to_dict(SeqIO.parse(ffasta, "fasta"))
    return records

def seq_from_pdb(pdb,selection='all',fmt='string'):
    """ Generate fasta from pdb entries """
    u = mda.Universe(pdb)
    ag = u.select_atoms(selection)
    if fmt == 'string':
        fastapdb = ""
    elif fmt == 'list':
        fastapdb = []
    else:
        raise
    res3 = ag.residues.resnames
    for res in res3:
        if len(res) == 3:
            res1 = SeqUtils.seq1(res)
        else:
            res1 = res
        if res1 == "":
            res1 = "X"
        fastapdb += res1
    return fastapdb

def write_fasta(seqs,fout,append=False):
    """ seqs: list of sequences """
    if append:
        records = []
        for record in SeqIO.parse(fout, "fasta"):
            records.append(record)
        for seq in seqs:
            records.append(seq)
        Bio.SeqIO.write(records, fout, "fasta")
    else:
        Bio.SeqIO.write(seqs, fout, "fasta")

def record_from_seq(seq,name):
    record = SeqRecord(
        Seq(seq),
        id=name,
        name='',
        description=''
    )
    return(record)

### SEQUENCE ANALYSIS

def get_qs(seq,flexhis=False,pH=7,calvados_version=2,residues=[]):
    """ charges and absolute charges vs. residues """
    qs, qs_abs = [], []
    # scaled charges
    if calvados_version == 4:
        qcoeff = 0.75
    else:
        qcoeff = 1.
    if len(residues) > 0: # residue charges provided
        for s in seq:
            q = residues.loc[s].q
            qs.append(q)
            qs_abs.append(abs(q))
    else: # guess residue charges
        # histidines
        if flexhis:
            qhis = qcoeff / ( 1 + 10**(pH-6) )
        else:
            qhis = 0.
        # loop through sequence
        for s in seq:
            if s in ['R','K']:
                qs.append(qcoeff)
            elif s in ['E','D']:
                qs.append(-1.*qcoeff)
            elif s == 'H':
                qs.append(qhis*qcoeff)
            else:
                qs.append(0.)
    qs = np.array(qs)
    qs_abs = np.abs(qs)
    return qs, qs_abs

def frac_charges(qs):
    N = len(qs)
    fpos = np.sum(np.where(qs>0, 1, 0)) / N
    fneg = np.sum(np.where(qs<0, 1, 0)) / N
    return fpos, fneg

def patch_terminal_qs(qs,calvados_version=2):
    qsnew = qs.copy()
    if calvados_version == 4:
        qcoeff = 0.75
    else:
        qcoeff = 1.0
    qsnew[0] += qcoeff
    qsnew[-1] -= qcoeff
    return qsnew

def seq_dipole(seq):
    """ 1D charge dipole along seq """
    # print(seq)
    qs, qs_abs = get_qs(seq)
    com = seq_com(qs_abs)
    dip = 0.
    # ct = 0
    for idx,q in enumerate(qs):
        dip += (com - idx) * q # positive if positive towards Nterm
    return com, dip

def seq_com(qs_abs):
    """ center of charges """
    com = 0.
    ct = 0.
    for idx,q_abs in enumerate(qs_abs):
        ct += int(q_abs)
        com += int(q_abs) * idx
    if ct > 0:
        com /= ct
    else:
        com = len(qs_abs) // 2
    return com

def calc_SCD(seq):
    """ Sequence charge decoration, eq. 14 in Sawle & Ghosh, JCP 2015 """
    qs, _ = get_qs(seq)
    N = len(seq)
    scd = 0.
    for idx in range(1,N):
        for jdx in range(0,idx):
            s = qs[idx] * qs[jdx] * (idx - jdx)**0.5
            scd += s
    scd /= N
    return scd

def calc_SHD(seq,residues,beta=-1.):
    """ Sequence hydropathy decoration, eq. 4 in Zheng et al., JPC Letters 2020"""
    N = len(seq)
    shd = 0.
    for idx in range(1,N):
        x = seq[idx]
        for jdx in range(0,idx):
            y = seq[jdx]
            s = (residues.lambdas[x] + residues.lambdas[y]) * (idx - jdx)**beta
            shd += s
    shd /= N
    return shd

def mean_lambda(seq,residues):
    """ Mean hydropathy """
    lambdas_sum = 0.
    for idx, x in enumerate(seq):
        lambdas_sum += residues.lambdas[x]
    lambdas_mean = lambdas_sum / len(seq)
    return  lambdas_mean

def calc_aromatics(seq):
    """ Fraction of aromatics """
    seq = str(seq)
    N = len(seq)    
    rY = len(re.findall('Y',seq)) / N
    rF = len(re.findall('F',seq)) / N
    rW = len(re.findall('W',seq)) / N
    return rY, rF, rW

def calc_mw(fasta):
    x = "".join(fasta)
    mw = Bio.SeqUtils.molecular_weight(x,seq_type='protein')
    return mw

def calc_kappa(seq):
    seq = "".join(seq)
    SeqOb = SequenceParameters(seq)
    k = SeqOb.get_kappa()
    return k

### SEQUENCE MANIPULATION
def shuffle_str(seq):
    l = list(seq)
    random.shuffle(l)
    seq_shuffled = "".join(l)
    return(seq_shuffled)

def construct_maxdipseq(seq):
    """ sequence permutation with largest dipole """
    seqpos, seqneg, seqneu = split_seq(seq)
    seqout = seqpos + seqneu + seqneg
    return seqout

def split_seq(seq):
    """ split sequence in positive, negative, neutral residues """
    seqpos = []
    seqneg = []
    seqneu = []
    for s in seq:
        if s in ['K','R']:
            seqpos.append(s)
        elif s in ['D','E']:
            seqneg.append(s)
        else:
            seqneu.append(s)
    seqpos = shuffle_str(seqpos)
    seqneg = shuffle_str(seqneg)
    seqneu = shuffle_str(seqneu)
    return seqpos, seqneg, seqneu

def mc_towards_kappa(seq,k_target=0.2,nsteps=100,dip_target=0.,a_kappa=1.,a_dip=1.,nswaps=1,
                    dipmax=0.):
    seq = "".join(seq)
    # clump charges
    k0 = calc_kappa(seq)
    print('-----------')
    if k_target >= 0.3:
        print('Calculating deltamax')
        seq = construct_deltamax(seq)
    else:
        print('Start from original seq')
    # print(seq)
    k0 = calc_kappa(seq)
    print('k before optimizing: ',f'{k0:.2f}')
    u0_k = k_energy(k0,k_target,a_kappa=a_kappa)
    u0_dip = dip_energy(seq,dip_target,a_dip=a_dip,dipmax=dipmax)
    u0 = u0_k + u0_dip
    # print(u0_k, u0_dip)
    # MC
    # print('-----------')
    print('RUNNING MC')
    # print(u0, seq)
    nacc = 0.
    nrej = 0.
    us = np.zeros((nsteps))
    ks = np.zeros((nsteps))
    dips = np.zeros((nsteps))
    umin = 10000.
    for idx in progressbar.progressbar(range(nsteps)):
        seqtemp = seq
        cswp = False
        for i in range(nswaps):
            seqtemp, charge_swap = trial_move(seqtemp)
            if charge_swap:
                cswp = True
        if cswp:
            k1 = calc_kappa(seqtemp)
            u1_k = k_energy(k1,k_target,a_kappa=a_kappa)
            u1_dip = dip_energy(seqtemp,dip_target,a_dip=a_dip,dipmax=dipmax)
            # print(u1_k, u1_dip)
            u1 = u1_k + u1_dip
            accept = metropolis(u0,u1,a=10*len(seq))
        else: # k should be identical for neutral swaps
            k1 = k0
            u1 = u0
            accept = True 
        if accept:
            nacc += 1
            seq = seqtemp
            k0 = k1
            u0 = u1
        else:
            nrej += 1
        ks[idx] = k0
        us[idx] = u0
        com, dip = seq_dipole(seq)
        dips[idx] = dip
        if u0 < umin:
            umin = u0
            seqmin = seq
    return seqmin, nacc/(nacc+nrej), us, ks, dips

def trial_move(seq):
    seq = list(seq)
    N = len(seq)
    while True:
        i = random.choice(range(N))
        j = random.choice(range(N))
        if i != j and seq[i] != seq[j]:
            break
    if (seq[i] in ['K','R','D','E']) or (seq[j] in ['K','R','D','E']):
        charge_swap = True
    else:
        charge_swap = False
    seq = swap_pos(seq,i,j)
    return seq, charge_swap

def k_energy(k,k_target,a_kappa=1.):
    # k = calc_kappa(seq)
    if k > k_target:
        u = a_kappa*(k-k_target)**2
    else:
        u = a_kappa*(k-k_target)**2
    return u

def dip_energy(seq,dip_target,a_dip=1.,dipmax=0.):
    com, dip = seq_dipole(seq)
    dip_red = dip / dipmax #dip_per_charge/len(seq)
    u = a_dip * (dip_red - dip_target)**2
    return u

def metropolis(u0,u1,a=100.):
    if u1 < u0:
        return True
    else:
        x = np.random.random()
        d = np.exp(a*(u0-u1))
        # print(x,d)
        if d > x:
            return True
        else:
            return False

def swap_pos(seq,i,j):
    seq = list(seq)
    itemp = seq[i]
    seq[i] = seq[j]
    seq[j] = itemp
    seq = "".join(seq)
    return seq

def construct_deltamax(seq):
    seqpos, seqneg, seqneu = split_seq(seq)
    # lr = np.random.randint(2)
    lr = 1 # force positive at N for now
    if len(seqneu) >= 18:
        neu3 = len(seqneu) // 3
        if lr == 1:
            seqout = seqneu[:neu3] + seqpos + seqneu[neu3:neu3*2] + seqneg + seqneu[neu3*2:]
        else:
            seqout = seqneu[:neu3] + seqneg + seqneu[neu3:neu3*2] + seqpos + seqneu[neu3*2:]
    else:
        if lr == 1:
            seqout = seqpos + seqneu + seqneg
        else:
            seqout = seqneg + seqneu + seqpos
    return seqout

class SeqFeatures:
    def __init__(self,seq,residues=[]):
        self.seq = seq
        self.N = len(seq)
        self.qs, self.qs_abs = get_qs(seq)
        self.charge = np.sum(self.qs)
        self.fpos, self.fneg = frac_charges(self.qs)
        self.ncpr = self.fpos-self.fneg
        self.fcr = self.fpos+self.fneg
        self.scd = calc_SCD(seq)
        self.rY, self.rF, self.rW = calc_aromatics(seq)
        self.kappa = calc_kappa(seq)
        self.com, self.dip = seq_dipole(seq)
        self.seqmax = construct_maxdipseq(seq)
        self.commax, self.dipmax = seq_dipole(self.seqmax)
        self.dipred = self.dip / self.dipmax
        if len(residues) > 0:
            self.lambdas_mean = mean_lambda(seq,residues)
            self.shd = calc_SHD(seq,residues,beta=-1.)