import numpy as np
import openmm
from openmm import unit

def genParamsDH(temp,ionic):
    """ Same for all components """
    kT = 8.3145*temp*1e-3
    # Calculate the prefactor for the Yukawa potential
    fepsw = lambda T : 5321/T+233.76-0.9297*T+0.1417*1e-2*T*T-0.8292*1e-6*T**3
    epsw = fepsw(temp)
    lB = 1.6021766**2/(4*np.pi*8.854188*epsw)*6.022*1000/kT
    eps_yu = np.sqrt(lB*kT)
    # Calculate the inverse of the Debye length
    k_yu = np.sqrt(8*np.pi*lB*ionic*6.022/10)
    return eps_yu, k_yu

def init_protein_interactions(eps_lj,cutoff_lj,k_yu,calvados_version=2):

    # harmonic bonds
    hb = openmm.openmm.HarmonicBondForce()
    hb.setUsesPeriodicBoundaryConditions(True)

    # intermolecular interactions
    energy_expression = 'select(step(r-2^(1/6)*s),4*eps*l*((s/r)^12-(s/r)^6-shift),4*eps*((s/r)^12-(s/r)^6-l*shift)+eps*(1-l))'
    if calvados_version in [1, 2]:
        ah = openmm.openmm.CustomNonbondedForce(energy_expression+'; s=0.5*(s1+s2); l=0.5*(l1+l2); shift=(0.5*(s1+s2)/rc)^12-(0.5*(s1+s2)/rc)^6')
    elif calvados_version == 3: # interactions scaled aromatics + R + H
        ah = openmm.openmm.CustomNonbondedForce(energy_expression+'; s=0.5*(s1+s2); l=sqrt(l1*l2)+m1*m2*0.8; shift=(0.5*(s1+s2)/rc)^12-(0.5*(s1+s2)/rc)^6')
    elif calvados_version == 4: # scaled charges
        ah = openmm.openmm.CustomNonbondedForce(energy_expression+'; s=0.5*(s1+s2); l=sqrt(l1*l2)+m1*m2*0.5; shift=(0.5*(s1+s2)/rc)^12-(0.5*(s1+s2)/rc)^6')
    else:
        raise

    ah.addGlobalParameter('eps',eps_lj*unit.kilojoules_per_mole)
    ah.addGlobalParameter('rc',float(cutoff_lj)*unit.nanometer)
    ah.addPerParticleParameter('s')
    ah.addPerParticleParameter('l')

    if calvados_version in [3,4]:
        ah.addPerParticleParameter('m')

    print('rc',cutoff_lj*unit.nanometer)

    yu = openmm.openmm.CustomNonbondedForce('q*(exp(-kappa*r)/r-shift); q=q1*q2')
    yu.addGlobalParameter('kappa',k_yu/unit.nanometer)
    yu.addGlobalParameter('shift',np.exp(-k_yu*4.0)/4.0/unit.nanometer)
    yu.addPerParticleParameter('q')
    
    ah.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
    yu.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)

    ah.setCutoffDistance(cutoff_lj*unit.nanometer)
    yu.setCutoffDistance(4.*unit.nanometer)
    
    ah.setForceGroup(0)
    yu.setForceGroup(1)

    return hb, ah, yu

def init_lipid_interactions(temp,eps_lipids,residues,cutoff_wca,cutoff_wcos):
    """ Define lipid interaction expressions """
    conv = constants.k * constants.N_A * temp / 1000. # kT in kJ/mol

    # wbond_expression = '-0.5*k*rinf^2*log(1-(r/rinf)^2)'
    # wbond = openmm.openmm.CustomBondForce(wbond_expression+'; k=10*eps/(4*R^2); rinf=3*R')
    # wbond.addGlobalParameter('eps',eps_lipids*conv*unit.kilojoules_per_mole)
    # wbond.addPerBondParameter('R')

    wbond_expression = '0.5*k*(r-req)^2'
    wbond = openmm.openmm.CustomBondForce(wbond_expression+'; k=10*eps/(4*R^2); req=8/2*R')
    wbond.addGlobalParameter('eps',eps_lipids*conv*unit.kilojoules_per_mole)
    wbond.addPerBondParameter('R')

    wbend_expression = '0.5*k*(r-req)^2'
    wbend = openmm.openmm.CustomBondForce(wbend_expression+'; k=10*eps/(4*R^2); req=8*R')
    wbend.addGlobalParameter('eps',eps_lipids*conv*unit.kilojoules_per_mole)
    wbend.addPerBondParameter('R')

    wca_expression = 'select(step(r-2^(1/6)*s),0,4*n*eps*((s/r)^12-(s/r)^6+1/4))'
    wca = openmm.openmm.CustomNonbondedForce(wca_expression+'; s=0.5*(s1+s2); n=n1*n2')
    wca.addGlobalParameter('eps',eps_lipids*conv*unit.kilojoules_per_mole)
    wca.addPerParticleParameter('s')
    wca.addPerParticleParameter('n') # check if lipid

    wcos_expression1 = 'select(step(r-2^(1/6)*s-wc),0,'
    wcos_expression2 = 'select(step(r-2^(1/6)*s),-m*n*eps*(cos(pi*(r-2^(1/6)*s)/(2*wc)))^2,-m*n*eps))'
    wcos = openmm.openmm.CustomNonbondedForce(wcos_expression1+wcos_expression2+'; s=0.5*(s1+s2); wc=1.78*btt; m=m1*m2; n=n1*n2')
    wcos.addGlobalParameter('eps',eps_lipids*conv*unit.kilojoules_per_mole)#eps_lipids*2.479*unit.kilojoules_per_mole)
    wcos.addGlobalParameter('btt',residues.loc['Z'].sigmas*unit.nanometer)
    wcos.addGlobalParameter('pi',np.pi)
    wcos.addPerParticleParameter('s')
    wcos.addPerParticleParameter('m') # check if tail
    wcos.addPerParticleParameter('n') # check if lipid

    wca.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
    wca.setCutoffDistance(cutoff_wca*unit.nanometer)
    wcos.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
    wcos.setCutoffDistance(cutoff_wcos*unit.nanometer)

    return wbond, wbend, wca, wcos

def init_restraints(restraint_type):
    """ Initialize restraints """
    if restraint_type == 'harmonic':
        cs = openmm.openmm.HarmonicBondForce()
    if restraint_type == 'go':
        go_expr = 'k*(5*(s/r)^12-6*(s/r)^10)'
        cs = openmm.openmm.CustomBondForce(go_expr+'; s=s; k=k')#; shift=(0.5*(s)/rc)^12-(0.5*(s)/rc)^6')
        cs.addPerBondParameter('s')
        cs.addPerBondParameter('k')
    cs.setUsesPeriodicBoundaryConditions(True)
    return cs

def init_scaled_LJ(eps_lj,cutoff_lj):
    """ Initialize restraints """
    energy_expression = 'select(step(r-2^(1/6)*s),n*4*eps*l*((s/r)^12-(s/r)^6-shift),n*4*eps*((s/r)^12-(s/r)^6-l*shift)+n*eps*(1-l))'
    scLJ = openmm.openmm.CustomBondForce(energy_expression+'; shift=(s/rc)^12-(s/rc)^6')
    scLJ.addGlobalParameter('eps',eps_lj*unit.kilojoules_per_mole)
    scLJ.addGlobalParameter('rc',float(cutoff_lj)*unit.nanometer)
    scLJ.addPerBondParameter('s')
    scLJ.addPerBondParameter('l')
    scLJ.addPerBondParameter('n')
    scLJ.setUsesPeriodicBoundaryConditions(True)
    return scLJ

def init_scaled_YU(k_yu):
    """ Initialize restraints """
    scYU = openmm.openmm.CustomBondForce('q*(exp(-kappa*r)/r-shift)')
    scYU.addGlobalParameter('kappa',k_yu/unit.nanometer)
    scYU.addGlobalParameter('shift',np.exp(-k_yu*4.0)/4.0/unit.nanometer)
    scYU.addPerBondParameter('q')
    scYU.setUsesPeriodicBoundaryConditions(True)
    return scYU

def init_eq_restraints(box,k):
    mindim = np.amin(box)
    rcent_expr = 'k*abs(periodicdistance(x,y,z,x,y,z0))'
    rcent = openmm.openmm.CustomExternalForce(rcent_expr)
    rcent.addGlobalParameter('k',k*unit.kilojoules_per_mole/unit.nanometer)
    rcent.addGlobalParameter('z0',box[2]/2.*unit.nanometer) # center of box in z
    # rcent.setNonbondedMethod(openmm.openmm.CustomNonbondedForce.CutoffPeriodic)
    # rcent.setCutoffDistance(mindim/2.*unit.nanometer)
    return rcent