import csv
import os
import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import pybedtools
import matplotlib.patheffects as path_effects
import scipy.stats
import statsmodels.api as sm
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
import statsmodels.stats.multitest as mt
import glob


df = pd.read_csv("GSM3756327_cast_pure.CAST.1.bedGraph",sep="\t")


plt.plot(df[df["#1_chr"]=="chr1"]["2_start"],
	df[df["#1_chr"]=="chr1"]["4_RT"])
plt.show()