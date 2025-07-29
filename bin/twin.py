#!/usr/bin/python3

import numpy as np, os, sys

from pathlib import Path

# fo = str(Path(os.path.dirname(os.path.realpath(__file__))).parent) + "/src"
# if fo not in sys.path:
#     sys.path.append(os.path.abspath(fo))

#
# adiciona automaticamente a pasta src/ (ao lado de bin/) ao PYTHONPATH,
# funcione onde for, seja local ou no cluster
#
this_file = Path(__file__).resolve()
project_root = this_file.parent.parent     # twinTrans/
src_dir      = project_root / "src"
if src_dir.is_dir() and str(src_dir) not in sys.path:
    sys.path.insert(0, str(src_dir))

import param_var as pv, functions as f

from warnings import filterwarnings

filterwarnings("ignore")

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
