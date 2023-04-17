# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 18:34:42 2023

@author: Marco
"""
import networkx as nx
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

file_name = "energy_file.csv"

df = pd.read_csv(file_name, comment='#', skipinitialspace=True)

plt.plot(df['time_step'],df['energy'])
plt.xlim([100000, None])
plt.ylabel('Energy')
plt.xlabel('time')
plt.ylim([None, 100])
