#!/usr/local/bin/python3
import numpy as np, os, sys

from pathlib import Path

fo = str(Path(os.path.dirname(os.path.realpath(__file__))).parent) + "/src"
if fo not in sys.path:
    sys.path.append(os.path.abspath(fo))
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
