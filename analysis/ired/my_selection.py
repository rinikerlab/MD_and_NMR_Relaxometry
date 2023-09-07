#!/usr/bin/env python3

from pyDR.Selection import select_tools as selt

import copy
import numpy as np

def select_caxis_iles(molsys, molsel):
    """
    This function selects the C_axis bonds of Isoleucines 
    used by Fabien Ferrage JACS 2018 to do the fitting
    
    Arguments
    ---------
        
    molsys: pyDIFRATE MolSys object
    molsel: pyDIFRATE MolSel object
    
    Both of these objects will be modified
 
    """
    sel0 = selt.sel_simple(molsys.uni.atoms)
    
    selec = sel0.select_atoms('resname ILE and name CG1')
    selec2 = sel0.select_atoms('resname ILE and name CD')
    
    l = len(selec)
    
    molsys.sel1 = []
    molsys.sel2 = []
    repr_sel = []
    
    for i in range(l):
        molsys.sel1.append(selec[i::l+1])
        molsys.sel2.append(selec2[i::l+1])
    
        repr_sel.append(molsys.sel1[i] + molsys.sel2[i])
            
    molsys._repr_sel = repr_sel
    molsys._label = [f'{int(atom.resids[0])}' for atom in molsys.sel1]
    
    molsel._sel1 = copy.copy(molsys.sel1)
    molsel._sel2 = copy.copy(molsys.sel2)
    
    return None

def select_5bb_bonds(molsys, molsel):
    """
    Selects 5 bonds from the backbone to do the analysis as in: 
    https://pubs.acs.org/doi/full/10.1021/ct500181v
    
    N-H, N-C_alpha, C_alpha - H_alpha, C_alpha - C_beta and C_alpha - C' 
    
    (except glycine which only has 4)
    
    N - H, and then all 4 bonds around the C_alpha (but not the C- side chain in glycine somehow)
    
    """
    
    # 1 : Make a basic selection object with the entire system

    system = selt.sel_simple(molsys.uni.atoms) 
    
    n_residues = len(system.residues)
    
    # To not have the problem of the sorting, we will do everything residue by residue
    
    selec = system.select_atoms('resid 1000')
    selec2 = system.select_atoms('resid 1000')
    
    for i in range(n_residues):
        selec  += system.select_atoms(f'resid {i+1} and name N and around 1.2 (name H or name HN)')
        selec2 += system.select_atoms(f'resid {i+1} and (name H or name HN) and around 1.2 name N ')
    
        # Next add the N - C_alpha 
        
        selec  += system.select_atoms(f'resid {i+1} and  name N and around 1.6 name CA')
        selec2 += system.select_atoms(f'resid {i+1} and name CA and around 1.6 name N ')
    
        # Next add the C_alpha to H_alpha 
    
        selec += system.select_atoms(f'resid {i+1} and name CA and around 1.2 (name HA )')
        selec2 += system.select_atoms(f'resid {i+1} and name HA and around 1.2 name CA')
    
        # Next add the C_alpha - C_beta
    
        selec += system.select_atoms(f'resid {i+1} and name CA and around 1.7 name CB')
        selec2 += system.select_atoms(f'resid {i+1} and name CB and around 1.7 name CA')
    
        # Finally add the C_alpha - C_carbonyl
    
        selec += system.select_atoms(f'resid {i+1} and name CA and around 1.7 name C')
        selec2 += system.select_atoms(f'resid {i+1} and name C and around 1.7 name CA')
    
    if len(selec) != len(selec2):
        print ('Warning - selections with different lengths !!!') 
    
    l = len(selec)
    
    repr_sel = [None] * l
    molsys._label = [None] * l
    
    for i, (s1, s2) in enumerate(zip(selec, selec2)):
        repr_sel[i] = s1 + s2
        molsys._label[i] = f'{s1.resname}_{s1.resid}_{s1.name}_{s2.name}'
  
    molsys.sel1 = selec
    molsys.sel2 = selec2
    
    molsel._sel1 = copy.copy(molsys.sel1)
    molsel._sel2 = copy.copy(molsys.sel2)
    
    return None

def find_anchored_5bbs(molsys, molsel):
    """
    Selects 5 bonds from the backbone to do the analysis as in: 
    https://pubs.acs.org/doi/full/10.1021/ct500181v
    
    N-H, N-C_alpha, C_alpha - H_alpha, C_alpha - C_beta and C_alpha - C' 
    
    (except glycine which only has 4)
    
    N - H, and then all 4 bonds around the C_alpha (but not the C- side chain in glycine somehow)
    
    """
    
    # 1 : Make a basic selection object with the entire system

    system = selt.sel_simple(molsys.uni.atoms) 
    
    # Start the lists with the N-H bonds (we will truncate two terminal residues in the end)
    
    selec = system.select_atoms('name N and around 1.2 (name H or name HN)')
    selec2 = system.select_atoms('(name H or name HN) and around 1.2 name N ')
    
    # Next add the N - C_alpha 
        
    selec += system.select_atoms('name N and around 1.6 name CA')
    selec2 += system.select_atoms('name CA and around 1.6 name N ')
    
    # Next add the C_alpha to H_alpha 
    
    selec += system.select_atoms('name CA and around 1.1 (name HA )')
    selec2 += system.select_atoms('name HA and around 1.1 name CA')
    
    # Next add the C_alpha - C_beta
    
    selec += system.select_atoms('name CA and around 1.8 name CB')
    selec2 += system.select_atoms('name CB and around 1.8 name CA')
    
    # Finally add the C_alpha - C_carbonyl
    
    selec += system.select_atoms('name CA and around 1.7 name C')
    selec2 += system.select_atoms('name C and around 1.7 name CA')
    
    if len(selec) != len(selec2):
        raise 
    
    # sort selections 1     
    sorted_idx = np.argsort(selec)

    l = len(selec)
    
    sel1 = [None] * l
    sel2 = [None] * l
    repr_sel = [None] * l
    labels = [None] * l 
    
    for k, i in enumerate(sorted_idx):
        sel1[k] = selec[i::l+1]
        sel2[k] = selec2[i::l+1]
        repr_sel[k] = sel1[k] + sel2[k]
    
    molsys._label = [None] * l
    
    for i, (a1, a2) in enumerate(zip(sel1, sel2)):
        if a1[0].resid != a2[0].resid:
            print ('problem in the sorting: ' + a1)
        labels[i] = f'{a1[0].resname}_{a1[0].resid}_{a1[0].name}-{a2[0].name}'
    


    return labels, (sel1, sel2), repr_sel



def find_anchored_methyls(molsys, molsel):
    """
    """
    
    # 1 : Make a basic selection object with the entire system

    system = selt.sel_simple(molsys.uni.atoms) 
    
    # Start the lists with the N-H bonds (we will truncate two terminal residues in the end)
    
    selec = system.select_atoms('name N and around 1.1 (name H or name HN)')
    selec2 = system.select_atoms('(name H or name HN) and around 1.1 name N ')
    
    # Next add the N - C_alpha 
        
    selec += system.select_atoms('name N and around 1.6 name CA')
    selec2 += system.select_atoms('name CA and around 1.6 name N ')
    
    # Next add the C_alpha to H_alpha 
    
    selec += system.select_atoms('name CA and around 1.2 (name HA)')
    selec2 += system.select_atoms('name HA and around 1.2 name CA')
    
    # Next add the C_alpha - C_beta
    
    selec += system.select_atoms('name CA and around 1.7 name CB')
    selec2 += system.select_atoms('name CB and around 1.7 name CA')
    
    # Finally add the C_alpha - C_carbonyl
    
    selec += system.select_atoms('name CA and around 1.7 name C')
    selec2 += system.select_atoms('name C and around 1.7 name CA')
    
    len_anchors = len(selec)
    
    # Now add the bonds for which we want to figure out the motion !!! 
    
    # Isoleucine gamma_2 and delta_1
    
    selec += system.select_atoms('resname ILE and name CG1')
    selec2 += system.select_atoms('resname ILE and name CD')
    
    selec += system.select_atoms('resname ILE and name CB')
    selec2 += system.select_atoms('resname ILE and name CG2')
   

    if len(selec) != len(selec2):
        raise 

    # sort selections 1     
    sorted_idx = np.argsort(selec)

    l = len(selec)
    
    sel1 = [None] * l
    sel2 = [None] * l
    repr_sel = [None] * l
    labels = [None] * l
    
    for k, i in enumerate(sorted_idx):
        sel1[k] = selec[i::l+1]
        sel2[k] = selec2[i::l+1]
        repr_sel[k] = sel1[k] + sel2[k]
        
        # Set label
        if sel1[k][0].resid != sel2[k][0].resid:
            print ('problem in the sorting: ' + sel1)
            return
        labels[i] = f'{sel1[k][0].resname}_{sel1[k][0].resid}_{sel1[k][0].name}-{sel2[k][0].name}'
    
   
    return labels, (sel1, sel2), repr_sel


def select_methyls_anchored(molsys, molsel):
    """
    Selects 5 bonds from the backbone to do the analysis as in: 
    https://pubs.acs.org/doi/10.1021/jacs.1c04687
    
    N-H, N-C_alpha, C_alpha - H_alpha, C_alpha - C_beta and C_alpha - C' 
    
    (except glycine which only has 4)
    
    N - H, and then all 4 bonds around the C_alpha (but not the C- side chain in glycine somehow)
    
    Then we have The Isoleucine, Valine, Leucine and Threonine terminal C - C bonds of the side chains.
    
    """
    
    # 1 : Make a basic selection object with the entire system

    system = selt.sel_simple(molsys.uni.atoms) 
    
    # Start the lists with the N-H bonds (we will truncate two terminal residues in the end)
    
    selec = system.select_atoms('name N and around 1.1 (name H or name HN)')
    selec2 = system.select_atoms('(name H or name HN) and around 1.1 name N ')
    
    # Next add the N - C_alpha 
        
    selec += system.select_atoms('name N and around 1.6 name CA')
    selec2 += system.select_atoms('name CA and around 1.6 name N ')
    
    # Next add the C_alpha to H_alpha 
    
    selec += system.select_atoms('name CA and around 1.2 (name HA)')
    selec2 += system.select_atoms('name HA and around 1.2 name CA')
    
    # Next add the C_alpha - C_beta
    
    selec += system.select_atoms('name CA and around 1.7 name CB')
    selec2 += system.select_atoms('name CB and around 1.7 name CA')
    
    # Finally add the C_alpha - C_carbonyl
    
    selec += system.select_atoms('name CA and around 1.7 name C')
    selec2 += system.select_atoms('name C and around 1.7 name CA')
    
    len_anchors = len(selec)
    
    # Now add the bonds for which we want to figure out the motion !!! 
    
    # Isoleucine gamma_2 and delta_1
    
    selec += system.select_atoms('resname ILE and name CG1')
    selec2 += system.select_atoms('resname ILE and name CD')
    
    selec += system.select_atoms('resname ILE and name CB')
    selec2 += system.select_atoms('resname ILE and name CG2')
   
    # Valine, gamma_1, gamma_2
    
    selec += system.select_atoms('resname VAL and name CB')
    selec2 += system.select_atoms('resname VAL and name CG1')
    
    selec += system.select_atoms('resname VAL and name CB')
    selec2 += system.select_atoms('resname VAL and name CG2')
    
    # Leucine, delta_1, delta_2
    selec += system.select_atoms('resname LEU and name CG')
    selec2 += system.select_atoms('resname LEU and name CD1')
    
    selec += system.select_atoms('resname LEU and name CG')
    selec2 += system.select_atoms('resname LEU and name CD2')
    
    # Threonine, gamma_2
    
    selec += system.select_atoms('resname THR and name CB')
    selec2 += system.select_atoms('resname THR and name CG2')
    
    if len(selec) != len(selec2):
        raise 
    
    # sort selections 1     
    sorted_idx = np.argsort(selec)

    l = len(selec)
    
    sel1 = [None] * l
    sel2 = [None] * l
    repr_sel = [None] * l
    labels = [None] * l
    
    for k, i in enumerate(sorted_idx):
        sel1[k] = selec[i::l+1]
        sel2[k] = selec2[i::l+1]
        repr_sel[k] = sel1[k] + sel2[k]
        
        # Set label
        if sel1[k][0].resid != sel2[k][0].resid:
            print ('problem in the sorting: ' + sel1)
            return
        labels[i] = f'{sel1[k][0].resname}_{sel1[k][0].resid}_{sel1[k][0].name}-{sel2[k][0].name}'
    
    molsys.sel1 = sel1
    molsys.sel2 = sel2
    
    molsel._sel1 = copy.copy(molsys.sel1)
    molsel._sel2 = copy.copy(molsys.sel2)
    
    molsys._label = labels
    
    return sorted_idx, labels, len_anchors

def selection_all_methyls(molsys, molsel):
    """
    This functions sets up the selection such that we have all the methyl group 
    that were measured by Rafael Bruschweiler. 
    
    """
    sel0 = selt.sel_simple(molsys.uni.atoms)
    
    # First Ile methyl
    selec = sel0.select_atoms('resname ILE and name CG1')
    selec2 = sel0.select_atoms('resname ILE and name CD')
    
    # Other Ile methyl 
    selec += sel0.select_atoms('resname ILE and name CG2')
    selec2 += sel0.select_atoms('resname ILE and name CB')
    
    # Leu methyl 
    selec += sel0.select_atoms('resname LEU and name CD1')
    selec2 += sel0.select_atoms('resname LEU and name CG')
    
    # Leu methyl 2
    selec += sel0.select_atoms('resname LEU and name CD2')
    selec2 += sel0.select_atoms('resname LEU and name CG')
    
    # Val methyl 
    selec += sel0.select_atoms('resname VAL and name CG1')
    selec2 += sel0.select_atoms('resname VAL and name CB')
    
    # Val methyl 2
    selec += sel0.select_atoms('resname VAL and name CG2')
    selec2 += sel0.select_atoms('resname VAL and name CB')
    
    # Threonine methyl
    selec += sel0.select_atoms('resname THR and name CG2')
    selec2 += sel0.select_atoms('resname THR and name CB') 
    
    # Alanine methyl
    selec += sel0.select_atoms('resname ALA and name CB')
    selec2 += sel0.select_atoms('resname ALA and name CA') 
    
    
    l = len(selec)
    
    if len(selec) != len(selec2):
        print ('problem')
        print (len(selec))
        print (len(selec2))
        return 
    molsys.sel1 = []
    molsys.sel2 = []
    repr_sel = []
    
    for i in range(l):
        #print (selec[i::l+1])
        #print (selec2[i::l+1])
        #print ('\n\n')
        
        molsys.sel1.append(selec[i::l+1])
        molsys.sel2.append(selec2[i::l+1])
    
        repr_sel.append(molsys.sel1[i] + molsys.sel2[i])
            
    molsys._repr_sel = repr_sel
    molsys._label = [f'{atom[0].resname}{int(atom[0].resid)}_{(atom[0].name)}' for atom in molsys.sel1]
    
    molsel._sel1 = copy.copy(molsys.sel1)
    molsel._sel2 = copy.copy(molsys.sel2)
