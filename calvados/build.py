import scipy as sp
import numpy as np
import pandas as pd
from openmm import unit

import MDAnalysis as mda
from MDAnalysis.analysis import distances

import json
import yaml

from .sequence import calc_mw

# BUILD DATAFRAME

def initDf():
    """ Initialize empty dataframe with default values """
    df = pd.DataFrame(columns=['ionic','temp','fasta','nmol','protein','rna','lipid'],dtype=object)
    return df

################ SYSTEM BUILDING FUNCTIONS ################

def build_box(Lx,Ly,Lz):
    # set box vectors
    a = unit.Quantity(np.zeros([3]), unit.nanometers)
    a[0] = Lx * unit.nanometers
    b = unit.Quantity(np.zeros([3]), unit.nanometers)
    b[1] = Ly * unit.nanometers
    c = unit.Quantity(np.zeros([3]), unit.nanometers)
    c[2] = Lz * unit.nanometers
    return a, b, c

def check_walls(x,box):
    """ x must be between 0 and [Lx,Ly,Lz] """
    if np.min(x) < 0: # clash with left box wall
        return True
    xbox = box
    d = xbox - x
    if np.min(d) < 0:
        return True # clash with right box wall
    return False # molecule in box

def check_clash(x,pos,box,cutoff=1.3):
    """ Check for clashes with other particles.
    Returns true if clash, false if no clash. """
    boxfull = np.append(box,[90,90,90])
    xothers = np.array(pos)
    if len(xothers) == 0:
        return False # no other particles
    d = distances.distance_array(x,xothers,boxfull)
    if np.amin(d) < cutoff:
        return True # clash with other particles
    else:
        return False # no clash

def draw_vec(l,ndim=3):
    """ 
    draw unbiased random vector with length l 
    """
    while True:
        vec = np.random.random(size=ndim) - 0.5
        if np.linalg.norm(vec) < 0.5:
            break
    vecnorm = vec / np.linalg.norm(vec)
    vecL = vecnorm * l
    return vecL

def draw_starting_vec(box):
    """
    draw random position within the simulation box
    """
    vec = np.random.random(size=3) * box
    return vec

def build_linear(n,d=0.38):
    """
    linear chain growing in z direction
    centered at [0,0,0]
    """
    x = np.array([[0,0,d*(z - 0.5*(n-1))] for z in range(n)])
    return x

def p2c(r, phi):
    """
    polar to cartesian
    """
    return (r * np.cos(phi), r * np.sin(phi))

def build_spiral(n, delta=0, arc=.38, separation=.7):
    """
    create points on an Archimedes' spiral
    with `arc` giving the length of arc between two points
    and `separation` giving the distance between consecutive 
    turnings
    """
    r = arc
    b = separation / (2 * np.pi)
    phi = float(r) / b
    coords = []
    for i in range(n):
        coords.append(list(p2c(r, phi))+[0])
        phi += float(arc) / r
        r = b * phi
    return np.array(coords)+delta

def build_xygrid(N,box,z=0.):
    """ Grid for xy slabs """
    if np.sqrt(N) % 1 > 0:
        b = 2
    else:
        b = 1
    nx = int(np.sqrt(N)) + b # nx spots in x dim
    ny = int(np.sqrt(N)) + b # ny spots in x dim

    dx = box[0] / nx
    dy = box[1] / ny

    xy = []
    x, y = 0., 0.
    ct = 0
    for n in range(N):
        ct += 1
        xy.append([x,y,z])
        if ct == ny:
            y = 0
            x += dx
            ct = 0
        else:
            y += dy
    xy = np.array(xy)
    return xy

def build_xyzgrid(N,box):
    """ 3D grid using HCP lattice """
    
    r = box / np.sum(box)
    a = np.cbrt(N / np.product(r))
    n = a * r
    nxyz = np.floor(n)
    while np.product(nxyz) < N:
        ndeviation = n / nxyz
        devmax = np.argmax(ndeviation)
        nxyz[devmax] += 1
    while np.product(nxyz) > N:
        nmax = np.argmax(nxyz)
        nxyz[nmax] -= 1
        if np.product(nxyz) < N:
            nxyz[nmax] += 1
            break
    
    xyz = []
    x, y, z = 0., 0., 0.

    ctx, cty, ctz = 0, 0, 0

    dx = box[0] / nxyz[0]
    dy = box[1] / nxyz[1]
    dz = box[2] / nxyz[2]

    zplane = 1
    xyplane = 1

    for n in range(N):
        if zplane > 0:
            xshift = 0
            yshift = 0
        else:
            xshift = dx/2
            yshift = dy/2

        if xyplane < 0:
            zshift = dz/2
        else:
            zshift = 0

        xyz.append([x+xshift,y+yshift,z+zshift])

        ctx += 1
        x += dx

        if (ctx % 2 == cty % 2):
            xyplane = 1
        else:
            xyplane = -1

        if ctx == nxyz[0]:
            ctx = 0
            x = 0

            cty += 1
            y += dy

            if cty == nxyz[1]:
                ctx = 0
                cty = 0
                x = 0.
                y = 0.
                
                ctz += 1
                z += dz

                zplane = -zplane

    xyz = np.array(xyz)
    return xyz

def slab_lipids(n_chains=100,d=0.85):
    npos = int(n_chains/2)
    print(f'npos: {npos}')
    n_per_dim = np.sqrt(npos)
    print(f'n_per_dim: {n_per_dim}')
    if n_per_dim % 1 == 0:
        margin = 0
    else:
        margin = 1
    n_per_dim = int(n_per_dim) + margin
    # n_per_dim = int(np.sqrt(npos))
    print(f'n_per_dim: {n_per_dim}')
    L = (n_per_dim)*d
    print('L: ',L)  
    xy = []
    x = 0.
    y = 0.
    ct = 0
    for n in range(npos):
        ct += 1
        xy.append([x,y])
        if ct == n_per_dim:
            y = 0
            x += d
            ct = 0
        else:
            y += d
    xy = np.array(xy)
    # print(xy)
    return L, xy

def build_vesicle(r,a_per_lipid=0.65,n_head=1,n_tail=2,h=4.,d=1.3,d_clash=0.5):
    x_lipids = []
    n_lipids = 0
    for m in [-1,1]:
        # m = 1: outer, m = -1: inner
        A = 4. * np.pi * (r+m*h/2.)**2
        n_shell = int(A / a_per_lipid)
        n_lipids += n_shell
        x_shell = place_on_sphere(r,h,d,d_clash,n_shell,m)
        for x in x_shell:
            x_lipids.append(x)
    x_lipids = np.array(x_lipids)
    return x_lipids, n_lipids

def place_on_sphere(r,h,d,d_clash,n,m):
    # m = 1: outer leaflet
    # m = -1: inner leaflet
    x = []
    for idx in range(n):
        # print(idx)
        while True:
            vec = draw_vec(r+m*h/2.)
            if len(x) > 0:
                clash = check_clash(vec,np.array(x),cutoff=d_clash)
            else:
                clash = False
            if not clash:
                vecnormal = vec / np.linalg.norm(vec)
                for jdx in [1,0,-1]:
                    x.append(vec+vecnormal*d*jdx*m)
                break
    return x

# FOLDED
def geometry_from_pdb(pdb,use_com=False):
    """ positions in nm"""
    u = mda.Universe(pdb)
    ag = u.atoms
    ag.translate(-ag.center_of_mass())
    if use_com:
        coms = []
        for res in u.residues:
            com = res.atoms.center_of_mass()
            coms.append(com)
        pos = np.array(coms) / 10.
    else:
        cas = u.select_atoms('name CA')
        pos = cas.positions / 10.
    return pos

def bfac_from_pdb(pdb,confidence=70.):
    """ get pLDDT encoded in pdb b-factor column """
    u = mda.Universe(pdb)
    bfac = np.zeros((len(u.residues)))
    for idx, res in enumerate(u.residues):
        bfac[idx] = np.mean(res.atoms.tempfactors) # average b-factor for residue
    bfac = np.where(bfac>confidence,bfac,0.) / 100. # high confidence filter
    return bfac

def load_pae_inv(input_pae,cutoff=0.1,colabfold=0):
    """ pae as json file (AF2 format) """
    pae = load_pae(input_pae,colabfold=colabfold)
    pae = np.where(pae < 1., 1, pae) # avoid division by zero (for i = j), min to 1
    pae_inv = 1/pae # inverse pae
    pae_inv = np.where(pae_inv > cutoff, pae_inv, 0)
    return pae_inv

def load_pae(input_pae,colabfold=0):
    """ pae as json file (AF2 format) """
    with open(input_pae) as f:
        pae = json.load(f)
        if colabfold == 0:
            pae = np.array(pae[0]['predicted_aligned_error'])
        elif colabfold == 1:
            pae = np.array(pae['predicted_aligned_error'])
        elif colabfold == 2:
            pae = np.array(pae['pae'])
    return pae

def get_ssdomains(name,fdomains):
    with open(f'{fdomains}','r') as f:
        stream = f.read()
        domainbib = yaml.safe_load(stream)

    domains = domainbib[name]
    print(f'Using domains {domains}')
    ssdomains = []

    for domain in domains:
        ssdomains.append(list(range(domain[0],domain[1])))
    return ssdomains

def check_ssdomain(ssdomains,i,j):
    ss = False
    for ssdom in ssdomains:
        if i in ssdom and j in ssdom:
            ss = True
    return ss

def conc_to_n(cinp,p,V,mw):
    """ g/L to n """
    n = p * 1e-24 * 1/mw * cinp * V * sp.constants.N_A
    return n

def n_to_conc(n,V,mw):
    """ n to g/L """
    c = n * mw * 1/sp.constants.N_A * 1e24 * 1/V
    return c

def calc_pair_n_in_box(cinp,pB,box,seqA,seqB):
    """ protein numbers in simulation box of two species """
    pA = 1.-pB
    V = box[0]*box[1]*box[2] # nm^3
    mwA = calc_mw(seqA) # Da = g/mol
    mwB = calc_mw(seqB) # Da = g/mol

    nA = round(conc_to_n(cinp,pA,V,mwA))
    nB = round(conc_to_n(cinp,pB,V,mwB))

    mA_rounded = nA*mwA
    mB_rounded = nB*mwB

    pB_rounded = mB_rounded / (mA_rounded + mB_rounded)
    c_rounded = n_to_conc(nA,V,mwA) + n_to_conc(nB,V,mwB)
    return nA, nB, pB_rounded, c_rounded

def calc_mixture_n_in_box(cinp,ps,box,seqs):
    """ 
    Protein numbers in simulation box of several species
    ps: list of mass proportions
    seqs: list of sequences
    """
    V = box[0]*box[1]*box[2] # nm^3

    ns = []
    mws = []

    for p, seq in zip(ps, seqs):
        mw = calc_mw(seq) # Da = g/mol
        mws.append(mw)
        n = round(conc_to_n(cinp,p,V,mw)) # number of proteins per type
        ns.append(n)

    cs_rounded = np.array([n_to_conc(n,V,mw) for n, mw in zip(ns,mws)]) # rounded g/L per type
    ctotal = np.sum(cs_rounded) # total g/L mass conc
    ps_rounded = cs_rounded / ctotal # mass fractions

    return ns, ps_rounded, ctotal