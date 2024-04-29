#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
This is the interface for the PC-SAFT dynamic library call.
Expose the so/dll functionality for two supported systems {LINUX,WIN}.
"""
__author__ = 'David Celný'
__version__ = '0.0.1'

# Built-in/Generic Imports
import os
import platform
# […]

# Libs
from ctypes import *
from typing import Optional,List
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
# Own modules
import pcsaft_interface

# misc parameters
ppad = 20 # print padding

# main code
# NOTE make sure that libpcsaft.dll is visible by the pcsaft_interface.py
# NOTE if not set the path as pcsaft_interface.pcsaft_eos(lib_path_str=<path string>)
pcsaft = pcsaft_interface.pcsaft_eos()

# ===== ===== ===== ===== ===== =====
# ===== EXAMPLES & TEST SECTION =====
# ===== ===== ===== ===== ===== =====
# NOTE uncomment individual sections to try them out
#      sections are separated by the print("\n=== ... statements

# print("\n=== get eos calculation constants ===")
# constants = pcsaft.get_eos_constants()
# print(f"used constants : {constants}")

# print("\n=== check eos internal state ===")
# print(pcsaft)
# pcsaft.set_fluid(['methane'])
# print(pcsaft)
# pcsaft.set_molar_ratio([1.0])
# print(pcsaft)
# pcsaft.set_temperature(150)
# print(pcsaft)
# pcsaft.set_density(1000)
# print(pcsaft)


# print("\n=== 1 comp simple output, externall setting ===")
# d = 22334.0
# pcsaft.set_fluid(['methane'])
# pcsaft.set_molar_ratio([1.0])
# pcsaft.set_temperature(150)
# pcsaft.set_density(d)

# p_1, dp_dr_1, dp_dr2_1 = pcsaft.get_p()

# print(f"D  : {d:{ppad}} mol/m^3 \nP  : {p_1:{ppad}} Pa \ndP : {dp_dr_1:{ppad}} Pa*m^3/mol \nd2P: {dp_dr2_1:{ppad}} Pa*m^6/mol^2")

# print("\n=== 1 comp simple output, function setting ===")
# rho = 22334.0 # mol/m^3
# t = 150 # K
# c = [1.0] # list of values in rnage (0,1>, sum of the list==1
# pcsaft.set_fluid(['methane'])

# p_1, dp_dr_1, dp_dr2_1 = pcsaft.get_p(density=d, temperature=t, composition=c)

# print(f"D  : {d:{ppad}} mol/m^3 \nP  : {p_1:{ppad}} Pa \ndP : {dp_dr_1:{ppad}} Pa*m^3/mol \nd2P: {dp_dr2_1:{ppad}} Pa*m^6/mol^2")

# print("\n=== 1 comp by parameters input comparison case ===")
# # BEWARE parameters has to be supplied as float numbers i.e 0.0 instead of 0
# methane_par = [1.0, 3.7039, 150.03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0160428, 0.0]
# pcsaft.set_fluid([methane_par])
# pcsaft.set_molar_ratio([1.0])
# pcsaft.set_temperature(150)
# pcsaft.set_density(100)

# p1, pp1, ppp1 = pcsaft.get_p()
# p2, pp2, ppp2 = pcsaft.get_p(substance=["methane"])

# print(f"P   : {p1:{ppad}}| dP : {pp1:{ppad}}| d2P: {ppp1:{ppad}}")
# print(f"P   : {p2:{ppad}}| dP : {pp2:{ppad}}| d2P: {ppp2:{ppad}}")

# print(f"diff: {p1 - p2:{ppad}}|    : {pp1 -pp2:{ppad}}|    : {ppp1 -ppp2:{ppad}}")


# print("\n=== 2 comp case detailed print of basic properties ===")
# s = ['methane','ethane']
# c = [0.5,0.5]
# d = 22334.0
# t = 150

# print(f"subtance   : {s} ")
# print(f"molar ratio: {c} ")
# print(f"density    : {d:{ppad-8}} mol/m^3")
# print(f"temperature: {t:{ppad-8}} K")
# print("-"*(5+ppad))
# pcsaft.set_fluid(s)
# pcsaft.set_molar_ratio(c)
# pcsaft.set_temperature(t) 
# pcsaft.set_density(d)

# p, pp, ppp = pcsaft.get_p()
# mu, fug = pcsaft.get_mu_fug()
# h, s, g = pcsaft.get_h_s_g()

# print(f"P  : {p:{ppad}} Pa \ndP : {pp:{ppad}} Pa*m^3/mol \nd2P: {ppp:{ppad}} Pa*m^6/mol^2")
# print(f"mu : {mu} #")
# print(f"fug: {fug} #")
# print(f"H  : {h:{ppad}} J \nS  : {s:{ppad}} J \nG  : {g:{ppad}} J")

# print("\n=== 2 comp by parameters input ,possible comparison vith previous 2 comp ===")
# # BEWARE parameters has to be supplied as float numbers i.e 0.0 instead of 0
# methane_par = (1.0, 3.7039, 150.03, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0160428, 0.0)
# ethane_par = (1.6069, 3.5206, 191.42, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03007, 0.0)
# pcsaft.set_fluid([methane_par, ethane_par])
# pcsaft.set_molar_ratio([0.5,0.5])
# pcsaft.set_temperature(150)
# pcsaft.set_density(22334.0)

# p1, pp1, ppp1 = pcsaft.get_p()

# print(f"P   : {p1:{ppad}}| dP : {pp1:{ppad}}| d2P: {ppp1:{ppad}}")
# try:
# 	print(f"P   : {p:{ppad}}| dP : {pp:{ppad}}| d2P: {ppp:{ppad}}")
# 	print(f"diff: {p1 - p:{ppad}}|    : {pp1 -pp:{ppad}}|    : {ppp1 -ppp:{ppad}}")
# except:
# 	pass


# # === ======== ======= ===
# # === ADVANCED SECTION ===
# # === ======== ======= ===

# print("\n=== 1 comp comparison between ideal gas and pcsaft ===")
# # use this to visualize a pressure in the plot
# p_selected = 1e6
# rho_min = 1
# rho_max = 23000

# pcsaft.set_fluid(['methane'])
# pcsaft.set_molar_ratio([1.0])
# pcsaft.set_temperature(150)
# rho = np.linspace(rho_min,rho_max,1000)
# p_id = np.zeros(np.shape(rho))
# p_saft = np.zeros(np.shape(rho))

# for ind,d in enumerate(rho):
# 	p_id[ind] = 8.314462618* 150 * d
# 	p_saft[ind], _, _ = pcsaft.get_p(density=d)

# plt.figure()
# plt.plot(rho,p_id,label="ideal gas")
# plt.plot(rho,p_saft,label="pcsaft")
# plt.plot([rho_min,rho_max],[p_selected,p_selected],label="selected")
# plt.xlabel("Density [mol/m^3]")
# plt.ylabel("Pressure [Pa]")
# plt.legend()

# print("\n=== 1 comp check the density calculation by recalculating pressure ===")
# p_selected = 1e6
# pcsaft.set_fluid(['methane'])
# pcsaft.set_molar_ratio([1.0])
# pcsaft.set_temperature(150)
# pcsaft.set_pressure(p_selected)

# d_1 = pcsaft.get_d()
# print(d_1)

# for d in d_1:
# 	p, _, _ = pcsaft.get_p(density=d)
# 	print(f"D  : {d:{ppad}} mol/m^3 | P  : {p:{ppad}} Pa")

# print("\n=== 1 comp test of density search over pressure ===")
# count = 500
# tt = 150

# pressures = np.linspace(0,6e6,count)
# rho_v = [None]*count
# rho_l = [None]*count
# rho_m = [None]*count

# pcsaft.set_fluid(['methane'])
# pcsaft.set_molar_ratio([1.0])
# pcsaft.set_temperature(tt)

# for p_ind,p in enumerate(pressures):
# 	# NOTE use is silent to omit the warnings and failures to converge
# 	rho_v[p_ind], rho_m[p_ind], rho_l[p_ind] = pcsaft.get_d(pressure=p, is_silent=True)

# pressures *=1e-6 # fot MPa
# plt.figure()
# plt.plot(rho_v,pressures,"or",label="vap")
# plt.plot(rho_m,pressures,"og",label="metastable")
# plt.plot(rho_l,pressures,"ob",label="liq")
# plt.xlabel("Density [mol/m^3]")
# plt.ylabel("Pressure [MPa]")
# plt.legend()

# print("\n=== 1 calculate single saturation properties ===")
# tt = 90.0
# pcsaft.set_fluid(['methane'])
# pcsaft.set_molar_ratio([1.0])
# pcsaft.set_temperature(tt)

# p_sat, d_vap, d_liq = pcsaft.get_equil()
# print(f"t    : {tt:{ppad}} K ")
# print(f"p_sat: {p_sat:{ppad}} Pa")
# print(f"d_vap: {d_vap:{ppad}} mol/m^3")
# print(f"d_liq: {d_liq:{ppad}} mol/m^3")

# print("\n=== show Methane NIST prop ===")
# df = pd.read_csv("methane_t150_prop.csv",sep="\t")
# # print(df) #DEBUG
# df.plot(x="Density (mol/m3)", y="Pressure (MPa)", style="sg")


# print("\n=== show Methane NIST sat prop ===")
# df2 = pd.read_csv("methane_t150_sat.csv",sep="\t")
# # print(df2) # DEBUG
# df2.plot(y=["Density (v, mol/m3)","Density (l, mol/m3)"], x="Pressure (MPa)", style=["sr","sb"])

# print("\n=== compare calculation to Methane from NIST  ===")
# df2 = pd.read_csv("methane_t150_sat.csv",sep="\t")

# pcsaft.set_fluid(['methane'])
# pcsaft.set_molar_ratio([1.0])
# p_sat = []

# for tt in df2["Temperature (K)"]:
# 	p_tmp, _, _ = pcsaft.get_equil(temperature=tt, pressure_estimate=100, is_silent=True)
# 	if p_tmp is not None:
# 		p_sat.append(p_tmp*1e-6)
# 	else:
# 		p_sat.append(None)

# plt.figure()
# plt.plot(df2["Temperature (K)"], df2["Pressure (MPa)"],label="NIST")
# plt.plot(df2["Temperature (K)"], p_sat, label="saft")
# plt.xlabel("Temperature [K]")
# plt.ylabel("Pressure [MPa]")
# plt.legend()

# NOTE keep this uncomented so that each plot does not have to provide its own show
plt.show()