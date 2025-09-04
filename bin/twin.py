#!/usr/bin/python3

import numpy as np, os, sys

from pathlib import Path

fo = str(Path(os.path.dirname(os.path.realpath(__file__))).parent) + "/src"
if fo not in sys.path:
    sys.path.append(os.path.abspath(fo))

#
# adiciona automaticamente a pasta src/ (ao lado de bin/) ao PYTHONPATH,
# funcione onde for, seja local ou no cluster
#
# this_file = Path(__file__).resolve()
# project_root = this_file.parent.parent     # twinTrans/
# src_dir      = project_root / "src"
# if src_dir.is_dir() and str(src_dir) not in sys.path:
#     sys.path.insert(0, str(src_dir))

import param_var as pv, functions as f


# SETTING UP
args = pv.parsing_cmd()
# parsing command line

modelP = pv.ModelParam(args)
# modelling parameters
simuP = pv.SimuParam(args)
# simulation parameters

modelP._test()
# basic tests for the value of some parameters

cmd = " ".join(sys.argv)
f.output_variables(cmd, modelP, simuP)
# writing out cmd line, parameters and variables (and their value)


# RUNNING
if simuP.promoter_to_follow:
    f.generate_run_follow_promoter(modelP, simuP)
else:
    f.generate_run_multiple_transcrtipts(modelP, simuP)




# if simuP.promoter_to_follow:
#     traj = f.generate_run_follow_promoter(modelP, simuP)
# else:
#     traj = f.generate_run_multiple_transcrtipts(modelP, simuP)

# # agora que a simulação acabou, calcula rON e rOFF
# # rON  = np.mean(traj.vel_on)  if traj.vel_on  else float('nan')
# # rOFF = np.mean(traj.vel_off) if traj.vel_off else float('nan')

# # print(f"rON  = {rON:.2f} nt/s")
# # print(f"rOFF = {rOFF:.2f} nt/s")

# # comprimento do gene em nt
# L = modelP.gene.L  

# # listas que agora existem:
# tesc_list  = traj.termination_escape_times  
# tterm_list = traj.termination_times        

# # 1) caso ON (primeiro RNAP que terminou)
# rON = L / (tterm_list[0] - tesc_list[0])

# # 2) caso OFF (aquele RNAP que escapou antes de t_off e terminou depois)
# t_off = modelP.promoter._t_off
# rOFF = None
# for tesc, tterm in zip(tesc_list, tterm_list):
#     if tesc < t_off < tterm:
#         rOFF = L / (tterm - tesc)
#         break

# if rOFF is None:
#     raise RuntimeError("Nenhum RNAP atende tesc < t_off < tterm; sua simulação não passou do ponto de OFF.")

# print(f"rON  = {rON:.2f} nt/s")
# print(f"rOFF = {rOFF:.2f} nt/s")
