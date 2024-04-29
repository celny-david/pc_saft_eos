from ctypes import *
from ctypes.util import find_library
from subprocess import Popen, PIPE
import numpy as np
import os
import time

## this is check for existence of exported function
# out = Popen(
#     args="nm ./libpcsaft.so",
#     shell=True,
#     stdout=PIPE
# ).communicate()[0].decode("utf-8")
#
# attrs = [
#     i.split(" ")[-1].replace("\r", "")
#     for i in out.split("\n") if " T " in i
# ]
#
# functions = [i for i in attrs if hasattr(CDLL("./libpcsaft.so"), i)]
#
# print(functions)

# dllpath = os.path.realpath("../../cmake-build-debug/libpcsaft_lib.so") # or .dll or .so

t_init = time.time()
dllpath = os.path.realpath("libpcsaft.so") # or .dll or .so
libtest = cdll.LoadLibrary(dllpath)

# pause = input("press enter to continue:")

# print(libtest)
#

# print("\nresult = {0:d}".format(res))
n_comp = 2;
rho = c_double(997.0)
comp_names = c_char_p("methane,carbon dioxide".encode('utf-8'))
comp_names2 = c_char_p("ethane,carbon dioxide".encode('utf-8'))
temperature = c_double(350.0)
molar_ratio = np.asarray([0.1,0.9],dtype=c_double())

p = c_double()
pp = c_double()
ppp = c_double()

mu = np.empty((n_comp,1),dtype=c_double())
fug = np.empty((n_comp,1),dtype=c_double())

h = c_double()
s = c_double()
g = c_double()

parse_input = libtest.parse_input
parse_input.argtypes = []

set_component_name = libtest.set_component_name
set_component_name.argtypes = [c_char_p, c_int]

set_density = libtest.set_density
set_density.argtypes = [c_double]

set_temperature = libtest.set_temperature
set_temperature.argtypes = [c_double]

set_molar_ratio = libtest.set_molar_ratio
set_molar_ratio.argtypes = [POINTER(c_double)]

press_calc = libtest.press_calc
press_calc.argtypes = [POINTER(c_double),POINTER(c_double),POINTER(c_double)]

chem_pot_calc = libtest.chem_pot_calc
chem_pot_calc.argtypes = [POINTER(c_double),POINTER(c_double)]

enthalpy_entropy_calc = libtest.enthalpy_entropy_calc
enthalpy_entropy_calc.argtypes = [POINTER(c_double),POINTER(c_double),POINTER(c_double)]

all_calc = libtest.all_calc
all_calc.argtypes = [POINTER(c_double),POINTER(c_double),POINTER(c_double),
                     POINTER(c_double),POINTER(c_double),
                     POINTER(c_double),POINTER(c_double),POINTER(c_double)]

print("\n=== python call sequence === \n")
t_calc = time.time()
print("\n=== parse call === \n")
parse_input()

print("\n=== set call === \n")
print(comp_names.value)
print(len(comp_names.value))

set_component_name(comp_names,len(comp_names.value))
set_density(rho)
set_temperature(temperature)
set_molar_ratio(molar_ratio.ctypes.data_as(POINTER(c_double)))

print("\n=== function call === \n")
# print("press")
# press_calc(byref(p),byref(pp),byref(ppp))
#
# print("chem_pot")
# chem_pot_calc( mu.ctypes.data_as(POINTER(c_double)),
#                fug.ctypes.data_as(POINTER(c_double)))
#
# print("enthalpy_entropy")
# enthalpy_entropy_calc(byref(h),byref(s),byref(g))

all_calc(byref(p),byref(pp),byref(ppp),
         mu.ctypes.data_as(POINTER(c_double)),
         fug.ctypes.data_as(POINTER(c_double)),
         byref(h),byref(s),byref(g))

print("\n=== function call 2 === \n")
set_component_name(comp_names2,len(comp_names2.value))
all_calc(byref(p),byref(pp),byref(ppp),
         mu.ctypes.data_as(POINTER(c_double)),
         fug.ctypes.data_as(POINTER(c_double)),
         byref(h),byref(s),byref(g))

t_end = time.time()


print()
print('pressure = %8.4f [Pa]' % (p.value))
print('pressure derivative = %8.4f [Pa*m**3/mol]' % (pp.value))
print('second pressure derivative = %8.4f [Pa*m**6/mol**2]' % (ppp.value))
print()
print('chemical potential = ', mu)
print('fugacity = ', fug)
print()
print('helmholtz energy = %8.4f [?]' % (h.value))
print('entropy = %8.4f [?]' % (s.value))
print('gibbs energy = %8.4f [?]' % (g.value))
print()
print('Calculation took: %f(total) %f(calc)' % (t_end - t_init, t_end - t_calc) )