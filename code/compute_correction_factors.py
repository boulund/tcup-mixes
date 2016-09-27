#!/usr/bin/env python3.5
# Compute correction factors for use with TCUP on mixed cultures
# Fredrik Boulund 2016

from sys import argv, exit, stderr
from collections import OrderedDict
import os
import argparse
import logging

import pandas as pd
import numpy as np
import scipy.stats
import matplotlib.pyplot as plt

os.environ["KMP_AFFINITY"] = "disabled" # Intel MKL is broken on LXSS
plt.style.use("ggplot")
logging.basicConfig(level=logging.DEBUG)
pd.options.display.width = None

columns = "Species  Sample  Peptides  Tot_dpeps  Correct_dpeps".split()
data  = [
    ("Staphylococcus aureus"      ,"QE_150611_128"  ,2392   ,552    ,549),
    ("Staphylococcus aureus"      ,"QE_150611_131"  ,2006   ,466    ,463),
    ("Staphylococcus aureus"      ,"QE_150611_137"  ,2479   ,571    ,565),
    ("Moraxella catarrhalis"      ,"QE_151124_11"   ,3531   ,2461   ,2617),
    ("Moraxella catarrhalis"      ,"QE_151124_08"   ,3469   ,2570   ,2539),
    ("Moraxella catarrhalis"      ,"QE_151124_05"   ,3119   ,2282   ,2253),
    ("Haemophilus influenzae"     ,"QE_151117_19"   ,2815   ,1137   ,1128),
    ("Haemophilus influenzae"     ,"QE_151117_16"   ,3283   ,1302   ,1289),
    ("Haemophilus influenzae"     ,"QE_151117_13"   ,3369   ,1239   ,1312),
    ("Streptococcus pneumoniae"   ,"QE_150508_38"   ,3283   ,356    ,329),
    ("Streptococcus pneumoniae"   ,"QE_150508_35"   ,3313   ,381    ,367),
    ("Streptococcus pneumoniae"   ,"QE_150508_32"   ,3491   ,417    ,401),
    ("Escherichia coli"   ,"QE_150611_140"   ,1996   ,262    ,253),
    ("Escherichia coli"   ,"QE_150611_143"   ,2825   ,316    ,294),
    ("Escherichia coli"   ,"QE_150611_146"   ,2911   ,331    ,311),
    ("Pseudomonas aeruginosa"   ,"QE_150611_152"   ,2756   ,1317    ,1303),
    #("Pseudomonas aeruginosa"   ,"QE_150611_155"   ,3268   ,1547    ,844),  # Weak sample; exclude!
    ("Pseudomonas aeruginosa"   ,"QE_150611_158"   ,1779   ,854    ,1527),
    ]


df = pd.DataFrame(data, columns=columns)

df.loc[:, "Prop_dpeps"] = df.Correct_dpeps / df.Peptides

print(df, file=stderr)

for species, row in df.groupby("Species").mean().iterrows():
    print(species, row.Prop_dpeps, sep="\t")

