#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 29 11:56:02 2019

@author: alfaceor
"""

import IgorModel
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
pd_pgen_current = pd.read_csv("ReDoTRbeta/CurrentIgorPgenMouse_output/Pgen_counts.csv", delimiter=';', index_col=0)
pd_pgen_new     = pd.read_csv("ReDoTRbeta/NewTRbeta_output/Pgen_counts.csv", delimiter=';', index_col=0)
#pgen_data = np.loadtxt("ReDoTRbeta/NewTRbeta_output/Pgen_counts.csv", delimiter=';', skiprows=1)

pd_pgen_current = pd_pgen_current.sort_index().dropna()
pd_pgen_new     = pd_pgen_new.sort_index().dropna()
print(len(pd_pgen_current), len(pd_pgen_new))

#fig, ax = plt.subplots()
#
#ax.set_xlim(0,500)
#ax.set_yscale("log")
#ax.plot(pd_pgen_current.index.values, pd_pgen_current['Pgen_estimate'].values)
#ax.plot(pd_pgen_new.index.values, pd_pgen_new['Pgen_estimate'].values)

fig, ax = plt.subplots( figsize=(8,6))

log10_pgen_current = np.log10(pd_pgen_current['Pgen_estimate'].values)
counts, bins = np.histogram(log10_pgen_current, bins=np.linspace(-20,-1, 25))
binsX = 0.5*(bins[1:] + bins[:-1])
histo = counts/float(counts.sum())
#ax.set_yscale("log")
ax.set_xlabel("$\log_{10}(Pgen)$")
ax.set_ylabel("Density")
ax.plot(binsX, histo, marker='o', label="Current ")


log10_pgen_new = np.log10(pd_pgen_new['Pgen_estimate'].values)
counts, bins = np.histogram(log10_pgen_new, bins=np.linspace(-20,-1, 25))
binsX = 0.5*(bins[1:] + bins[:-1])
#binsX = bins[:-1]
histo = counts/float(counts.sum())
#ax.set_yscale("log")
ax.set_xlabel("$\log_{10}(Pgen)$")
ax.set_ylabel("Density")
ax.plot(binsX, histo, marker='s', label="New ")
ax.legend(bbox_to_anchor=(0.75, 0.72), loc="lower left", prop={'size':15})

fig.tight_layout()
fig.savefig("Pgen_compar02.pdf")










"""

fig, ax = plt.subplots(2,1, figsize=(8,6))

log10_pgen_current = np.log10(pd_pgen_current['Pgen_estimate'].values)
counts, bins = np.histogram(log10_pgen_current, bins=np.linspace(-20,-1, 25))
bins
#binsX = 0.5*(bins[:1] + bins[:-1])
binsX = bins[:-1]
counts.sum()
histo = counts #/float(counts.sum())
#ax.set_yscale("log")
ax[0].set_xlabel("$\log_{10}(Pgen)$")
ax[0].set_ylabel("Density")
ax[0].plot(binsX, histo, marker='o')

pgen_current = pd_pgen_current['Pgen_estimate'].values
counts, bins = np.histogram(pgen_current, bins=np.logspace(-20,-1, 25))
#binsX = 0.5*(bins[:1] + bins[:-1])
binsX = bins[:-1]
histo = counts #/float(counts.sum())
ax[1].set_xscale("log")
ax[1].set_xlabel("$Pgen$")
ax[1].set_ylabel("Density")
ax[1].plot(binsX, histo, marker='o')


log10_pgen_new = np.log10(pd_pgen_new['Pgen_estimate'].values)
counts, bins = np.histogram(log10_pgen_new, bins=np.linspace(-20,-1, 25))
#binsX = 0.5*(bins[:1] + bins[:-1])
binsX = bins[:-1]
histo = counts #/float(counts.sum())
#ax.set_yscale("log")
ax[0].set_xlabel("$\log_{10}(Pgen)$")
ax[0].set_ylabel("Density")
ax[0].plot(binsX, histo, marker='o')


pgen_new = pd_pgen_new['Pgen_estimate'].values
counts, bins = np.histogram(pgen_new, bins=np.logspace(-20,-1,25))
#binsX = 0.5*(bins[:1] + bins[:-1])
binsX = bins[:-1]
histo = counts #/float(counts.sum())
ax[1].set_xscale("log")
ax[1].set_xlabel("$Pgen$")
ax[1].set_ylabel("Density")
ax[1].plot(binsX, histo, marker='o')

fig.tight_layout()
"""


#np.linspace(-30,-1, 100)
#np.log(np.logspace(-30,-1,100))
#
#histo.max()

#
#pd_pgen.loc[0]
#pd_pgen.loc[99995]
#
#np.linspace()
#3.49*log10_pgen.std()*len(log10_pgen)**(-0.3333)
#print(log10_pgen.min(), log10_pgen.max() )
