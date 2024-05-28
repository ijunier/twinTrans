#!/usr/local/bin/python3

import os, sys
from pathlib import Path
import subprocess
from joblib import Parallel, delayed
import argparse
import itertools

parser = argparse.ArgumentParser()
parser.add_argument('N', type=int, help='number of launched jobs')
parser.add_argument('-n', metavar='X', type=int, help='nice level', default=10)
args=parser.parse_args()

script='nice -%d '%args.n+os.path.dirname(os.path.realpath(__file__))+'/twin.py' 
parent_fo=str(Path(os.path.dirname(os.path.realpath(__file__))).parent) # absolute parent folder 

def run(cmd): # necessary for running Parallel()
    subprocess.run(cmd.split(' '))

LasGTs=[[0.7,1.4],[0.5,1.4]]
#ks=[0.01, 0.025, 0.05, 0.1, 0.25, 0.5, 1, 2.5, 5, 10]
ks=[0.025, 0.05, 0.1, 0.25, 0.5, 1, 2.5]
kbs,kos=ks[::-1],ks[::-1]
distances=[250, 500, 1000]
sos=[0, 0.025, 0.05]
    
cmds, fos = [], []
# cmds = list of commands using script
# fos: list of corresponding result folders (created on-the-fly)

for LasGT, kb, ko, so, d in itertools.product(LasGTs,kbs,kos,sos,distances):
    LasG,LasT=LasGT[0],LasGT[1]
    foo=parent_fo+'/results/2024_05/Degradation_Scenario/LasG%s_LasT%s/_kbkoso/'%(str(LasG),str(LasT))
    if not os.path.exists(foo): os.makedirs(foo) 
    fos.append(foo+'kb%s_ko%s_so%s_d%d'%(str(kb),str(ko),str(so),d)) 
    cmds.append('%s -LasG %s -LasT %s -kb %s -ko %s -so -%s -tss %d -Ld %d -f %s'%(script,str(LasG),str(LasT),str(kb),str(ko),str(so),d,d,fos[-1])) 

Parallel(n_jobs=args.N)(delayed(run)(cmd) for cmd in cmds)

