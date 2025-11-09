#!/usr/bin/env python
# coding: utf-8

# In[25]:


import os
import shutil
from pathlib import Path
import sys
from glob import glob
import numpy as np
import MDAnalysis


# In[2]:


import rdkit
from rdkit import Chem
from rdkit.Chem import AllChem


# In[3]:


from Bio import Align
aligner = Align.PairwiseAligner()
# https://www.biostars.org/p/424006/


# In[4]:


import logging

def _find_proc_path(procname):
    cmd = 'which %s'%procname
    RT = os.popen(cmd).read().strip()
    if not RT:
        logging.error('ERROR found %s'%procname)
        exit(1)
    return RT


# In[5]:


GROMACS_EXE = _find_proc_path('gmx')
#ACPYPE_EXE = _find_proc_path('acpype')
#ACPYPE_EXE = '/home/naf/software/miniconda3/envs/gbsa/bin/acpype'


# In[6]:


GROMACS_EXE


# In[7]:


capp_res = ['ACE','NME','NHE','NMA','NH2']
protein_body = ['ALA','GLY','SER','THR','LEU','ILE','VAL','ASN','GLN','ARG','HID','HIE','HIP','TRP',
                'PHE','TYR','GLU','ASP','LYS','ORN','DAB','LYN','PRO','HYP','CYS','CYM','CYX','MET','ASH','GLH']
water_body = ['HOH','SOL','HO4','WAT','TIP']


# In[8]:


am_mapping = {'ALA':'A','GLY':'G','SER':'S','THR':'T','LEU':'L','ILE':'I','VAL':'V',
              'ASN':'N','GLN':'Q','ARG':'R','HID':'H','HIE':'H','HIP':'H','TRP':'W',
              'PHE':'F','TYR':'Y','GLU':'E','ASP':'D','LYS':'K','LYN':'K','PRO':'P',
              'CYS':'C','CYM':'C','CYX':'C','MET':'M','ASH':'D','GLH':'E'}


# In[9]:


def SVD_align(Obj, Ref):
    c_o = Obj.mean(axis=0)
    c_r = Ref.mean(axis=0)
    coo_o = Obj - c_o
    coo_r = Ref - c_r
    C = np.dot(np.transpose(coo_o), coo_r)
    U, S, Vt = np.linalg.svd(C)
    flag = float(str(float(np.linalg.det(U) * np.linalg.det(Vt))))
    if flag < 0:
        S[-1] = -S[-1]
        U[:,-1] = -U[:,-1]
    R = np.dot(U, Vt)
    t = -np.dot(c_o, R) + c_r
    return R,t


# In[10]:


get_ipython().system('pwd')


# In[11]:


#os.chdir('/home/naf/Documents/ihuman/CXCR3/simulation/')


# In[12]:


# inuputs
mem_protein = './membrane/final-system-dry.pdb'
mem_dir = './membrane/'
#ori_protein = './docking_poses/Structure2.pdb'
#input_ligand_dir = './docking_poses/'


# In[13]:


input_protein_dir = './inputpdbs/'


# In[14]:


# tmp dirs
format_dir = './alignment/'
if not os.path.isdir(format_dir):
    os.mkdir(format_dir)


# In[15]:


# MD dirs
md_dir = './MD/'
if not os.path.isdir(md_dir):
    os.mkdir(md_dir)


# In[16]:


mem_protein = os.path.abspath(mem_protein)
mem_dir = os.path.abspath(mem_dir)
input_protein_dir = os.path.abspath(input_protein_dir)
format_dir = os.path.abspath(format_dir)
md_dir = os.path.abspath(md_dir)


# In[17]:


# process protein


# In[18]:


logfile = './log'
nt = 10


# In[19]:


paras = {
        'gmx':GROMACS_EXE,
        'logfile':os.path.abspath(logfile),
        'mem_protein':mem_protein,
        'format_dir':format_dir,
        'gmx_dir':mem_dir,
        'nt':nt,
        'amberff':'amber99sb-ildn',
        'lipidName':'POPC'
    }


# In[20]:


# for every protein in input_protein_dir, i need to align them to membrane protein, B chain
# but A chain should also moved


# In[21]:


# ---------- 1.  read membrane B-chain CA once ----------
mem_lines = open(mem_protein).readlines()
mem_CA, mem_seq = [], ''
for line in mem_lines:
    if line[0:4] == 'ATOM':
        if line[17:20] in am_mapping and line[12:16].strip() == 'CA' and line[21] == 'B':
            x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
            mem_CA.append([x, y, z])
            mem_seq += am_mapping[line[17:20]]
mem_CA = np.array(mem_CA)

# ---------- 2.  your original SVD ----------
def SVD_align(Obj, Ref):
    c_o = Obj.mean(axis=0)
    c_r = Ref.mean(axis=0)
    coo_o = Obj - c_o
    coo_r = Ref - c_r
    C = np.dot(np.transpose(coo_o), coo_r)
    U, S, Vt = np.linalg.svd(C)
    flag = float(str(float(np.linalg.det(U) * np.linalg.det(Vt))))
    if flag < 0:
        S[-1] = -S[-1]
        U[:, -1] = -U[:, -1]
    R = np.dot(U, Vt)
    t = -np.dot(c_o, R) + c_r
    return R, t

# ---------- 3.  process every PDB ----------
for pdb_path in glob(os.path.join(input_protein_dir, '*.pdb')):
    protein_lines = open(pdb_path).readlines()

    # 3a.  grab chain-B CA  (your old loop, chain-B only)
    pro_CA, pro_seq = [], ''
    for line in protein_lines:
        if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
            if line[17:20] in am_mapping and line[12:16].strip() == 'CA' and line[21] == 'B':
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                pro_CA.append([x, y, z])
                pro_seq += am_mapping[line[17:20]]
    pro_CA = np.array(pro_CA)

    # 3b.  align sequences  (your old code)
    alignments = aligner.align(mem_seq, pro_seq)
    seq1, seq2 = alignments[0][0], alignments[0][1]
    i = j = 0
    mem_CA_, pro_CA_ = [], []
    for k in range(len(seq1)):
        if seq1[k] == '-': j += 1; continue
        if seq2[k] == '-': i += 1; continue
        mem_CA_.append(mem_CA[i])
        pro_CA_.append(pro_CA[j])
        i += 1; j += 1
    ret_R, ret_t = SVD_align(np.array(pro_CA_), np.array(mem_CA_))

    # 3c.  apply same move to **every atom** in the file  (your old PDB write)
    base = os.path.basename(pdb_path).replace('.pdb', '.pdb')
    with open(pdb_path) as f, open(os.path.join(format_dir, base), 'w') as fw:
        for line in f:
            if line[0:4] == 'ATOM' or line[0:6] == 'HETATM':
                x = float(line[30:38]); y = float(line[38:46]); z = float(line[46:54])
                [r_x, r_y, r_z] = np.dot(np.array([x, y, z]), ret_R) + ret_t
                r_x = "%.3f" % r_x; r_y = "%.3f" % r_y; r_z = "%.3f" % r_z
                r_x = str(r_x) + ' ' * (8 - len(str(r_x)))
                r_y = str(r_y) + ' ' * (8 - len(str(r_y)))
                r_z = str(r_z) + ' ' * (8 - len(str(r_z)))
                fw.write(line[:30] + r_x + r_y + r_z + line[54:])
                continue
            fw.write(line)


# In[22]:


# now the alignment is ok
# I need to prepare MD dir for every aligned protein


# In[23]:


for protein_idx, pdb_path in enumerate(glob(os.path.join(format_dir, '*.pdb'))):
    print(protein_idx, pdb_path)
    mark = pdb_path.split('/')[-1][:-4]

    run_dir = os.path.join(md_dir, mark)
    if not os.path.isdir(run_dir):
        os.mkdir(run_dir)
    print(run_dir)

    # get ori chains
    f = open(pdb_path)
    chain_ = []
    chains_ = []
    lines = f.readlines()
    flag = 0
    lastresnum = 0
    for line in lines:
        if line[0:4] != 'ATOM':
            continue
        atomNam = line[12:16].strip()
        resnum = int(line[22:26])
        chainID = line[21:22]
        if flag == 0:
            lastresnum = resnum - 1
            lastchainID = chainID
            flag = 1
        if (resnum != lastresnum and resnum != lastresnum + 1) or chainID != lastchainID:
            chains_.append(chain_[:])
            chain_ = []
            chain_.append(line)
            lastresnum = resnum
            lastchainID = chainID
            continue
        chain_.append(line)
        lastresnum = resnum
        lastchainID = chainID
    chains_.append(chain_[:])
    f.close()

    for i in range(len(chains_)):
        fw = open(os.path.join(run_dir,'protein_'+str(i)+'.pdb'),'w')
        for line in chains_[i]:
            fw.write(line)
        fw.close()

    # generate top
    os.chdir(run_dir)
    for i in range(len(chains_)):
        paras['num'] = i
        logging.info('Generate top for each chain.')
        cmd = 'gmx pdb2gmx -f protein_{num}.pdb -o protein_{num}_ori.pdb -p protein_{num}.top -i protein_{num}_posre.itp -ff {amberff} -water tip3p -ignh>>{logfile} 2>&1'.format(**paras)
        RC = os.system(cmd)
        if RC!=0:
            logging.error('Failed generate top file, see the %s for details.'%os.path.abspath(logfile))
            exit(1)
    os.chdir('../../')

    # move top to itp
    for i in range(len(chains_)):
        input_file = os.path.join(run_dir, "protein_%s.top" % str(i))
        out_file = os.path.join(run_dir, "protein_%s.itp" % str(i))
        lines = open(input_file).readlines()
        fw = open(out_file, 'w')
        flag = False
        for line in lines:
            if "forcefield.itp" in line:
                continue
            if "Include water topology" in line:
                break
            if flag:
                line = "protein_%s     3\n" % str(i)
                flag = False
            if 'nrexcl' in line:
                flag = True
            fw.write(line)
        fw.close()

    for i in range(len(chains_)):
        input_file = os.path.join(run_dir, "protein_%s_posre.itp" % str(i))
        out_file = os.path.join(run_dir, "protein_%s_posre.itp" % str(i))
        lines = open(input_file).readlines()
        fw = open(out_file, 'w')
        flag = False
        for line in lines:
            fw.write(line)
        fw.close()

    # generate a new top
    with open(os.path.join(mem_dir,'topol_final_dry.top')) as fr:
        lines = fr.readlines()

    POPC_itp_flag = 0
    POPC_flag = 0
    flag = 0
    with open(os.path.join(run_dir, 'topol.top'), 'w') as fw:
        for line in lines:
            if line.startswith('#ifdef'):
                break
            if line.startswith('#include') and 'forcefield' in line:
                fw.write(line)
                for i in range(len(chains_)):
                    fw.write('#include "protein_%s.itp"\n' % str(i))
                continue
            if line.startswith('#include') and 'POPC' in line:
                POPC_itp_flag = 1
                fw.write(line)
                continue
            if line.startswith('#include') and POPC_itp_flag == 0:
                continue
            if line.startswith('POPC'):
                POPC_flag = 1
                fw.write(line)
                continue
            if line.startswith('TIP3P'):
                continue
            if line.startswith('NA'):
                continue
            if line.startswith('CL'):
                continue
            if line.startswith('PROTEIN'):
                continue
            if line.startswith('; Compound'):
                fw.write(line)
                flag = 1
                for i in range(len(chains_)):
                    fw.write("protein_%s    1\n" % str(i))
                continue
            if flag == 1 and POPC_flag == 0:
                continue
            fw.write(line)

    #proteinItpFiles = glob(mem_dir+'/protein_*.itp')
    #for itpfile in proteinItpFiles:
    #    filename = os.path.split(itpfile)[-1]
    #    if "disres" in filename or "ultra_posre" in filename or "ca_posre" in filename:
    #        continue
    #    with open(itpfile) as fr:
    #        lines = fr.readlines()
    #    POSRES = True
    #    with open(os.path.join(run_dir, filename), 'w') as fw:
    #        for line in lines:
    #            if line.startswith('#'):
    #                if POSRES:
    #                    fw.write('#ifdef POSRES\n#include "%s_posre.itp"\n#endif\n'%filename[:-4])
    #                    POSRES = False
    #                continue
    #            fw.write(line)   

    # single-file list
    extra_itps = ['CL.itp', 'NA.itp', 'POPC.itp', 'POPC_posre.itp', 'TIP3P.itp']

    for itp in extra_itps:
        src = os.path.join(mem_dir, itp)
        dst = os.path.join(run_dir, itp)
        if os.path.isfile(src):
            shutil.copy(src, dst)
        else:
            raise FileNotFoundError(f'{src} not found – check mem_dir')

    # force-field directory
    src_ff = os.path.join(mem_dir, 'amber99sb-ildn_slipids.ff')
    dst_ff = os.path.join(run_dir, 'amber99sb-ildn_slipids.ff')

    if os.path.isdir(src_ff):
        if os.path.exists(dst_ff):          # avoid “already exists” error
            shutil.rmtree(dst_ff)
        shutil.copytree(src_ff, dst_ff)
    else:
        raise FileNotFoundError(f'{src_ff} not found – check mem_dir')

    # generate a new pdb
    fw = open(os.path.join(run_dir, 'start_complex.pdb'), 'w')
    f = open(os.path.join(mem_dir, 'final-system-dry.pdb'))
    lines_ = f.readlines()
    fw.write(lines_[0])
    fw.write(lines_[1])
    fw.write(lines_[2])
    fw.write(lines_[3])
    for i in range(len(chains_)):
        input_file = os.path.join(run_dir, "protein_%s_ori.pdb" % str(i))
        lines = open(input_file).readlines()
        for line in lines:
            if line[0:4] == 'ATOM':
                fw.write(line)
        fw.write('TER\n')
    for line in lines_:
        if line[17:21] == paras['lipidName']:
            fw.write(line)
    f.close()
    fw.close()


# In[40]:


# prepare submit files


# In[26]:


src_path  = Path('mdps')
dest_root = Path('MD')

for dest_dir in dest_root.iterdir():
    if dest_dir.is_dir():
        for mdp in src_path.iterdir():
            if mdp.is_file():
                print(f"[COPY] {mdp.name}  ->  {dest_dir}/")
                shutil.copy2(mdp, dest_dir)
print("Finished.")


# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




