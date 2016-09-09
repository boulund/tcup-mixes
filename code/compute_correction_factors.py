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

os.environ["KMP_AFFINITY"] = "disabled"
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
    #("Moraxella catarrhalis"      ,"QE_151124_05"   ,,,),
    #("Haemophilus influenzae"     ,"QE_151117_19"   ,,,),
    ("Haemophilus influenzae"     ,"QE_151117_16"   ,3283   ,1302   ,1289),
    ("Haemophilus influenzae"     ,"QE_151117_13"   ,3369   ,1239   ,1312),
    ("Streptococcus pneumoniae"   ,"QE_150508_38"   ,3283   ,356    ,329),
    ("Streptococcus pneumoniae"   ,"QE_150508_35"   ,3313   ,381    ,367),
    ("Streptococcus pneumoniae"   ,"QE_150508_32"   ,3491   ,417    ,401),
    ]


df = pd.DataFrame(data, columns=columns)

df.loc[:, "Prop_dpeps"] = df.Correct_dpeps / df.Peptides

print(df, file=stderr)

for species, row in df.groupby("Species").mean().iterrows():
    print(species, row.Prop_dpeps, sep="\t")

